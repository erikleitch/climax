#include "gcp/datasets/PsfImageDataSet.h"

#include "gcp/fftutil/Generic2DAngularModel.h"

using namespace std;

using namespace gcp::datasets;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
PsfImageDataSet::PsfImageDataSet() 
{
  addParameter("frequency",   DataType::DOUBLE, "When simulating, the frequency at which to calculate the model");
  addParameter("noiserms",    DataType::DOUBLE, "Noise sigma to use (if dist=gauss)");
  addParameter("background",  DataType::DOUBLE, "Background count rate");

  addParameter("thetaMin",    DataType::DOUBLE, "If specified, data < thetaMin of the current model center will be excluded from the fit");
  addParameter("thetaMax",    DataType::DOUBLE, "If specified, data > thetaMax of the current model center will be excluded from the fit");

  addParameter("thetaMinErr", DataType::DOUBLE, "If specified, data outside of thetaMinErr will be used to estimate the error.  Note that this can conflict with 'region'");
  addParameter("errRegion",   DataType::STRING, "A region from which to estimate errors.  Note that this can conflict with 'region'");
  addParameter("errImage",    DataType::STRING, "The image from which errors will be estimated, if thetaMinErr or errRegion are set.  Set to 'pre' to estimate from the original image, 'post' to estimate from the final image");

  addParameter("psf",         DataType::STRING, "Type of psf: 'realistic', 'gauss' or 'none'.  Default is 'none'");
  addParameter("dist",        DataType::STRING, "Type of error distribution: 'gauss(ian)' or 'poiss(on)'.  Will default to 'poisson' for X-ray data, 'gauss' for all others");

  addParameter(imp_);

  setParameter("errImage", "post");

  psfType_ = PSF_NONE;
}

/**.......................................................................
 * Destructor.
 */
PsfImageDataSet::~PsfImageDataSet() {}

/**.......................................................................
 * Overload the base-class to initialize number of ants to 1
 */
void PsfImageDataSet::setName(std::string name)
{
  DataSet::setName(name);
  obs_.setParameter("nant", "1");
}

void PsfImageDataSet::loadData(bool simulate)
{
  if(simulate) {
    initializeForSim();
  } else {
    loadDataFromFile();
  }
}

/**.......................................................................
 * Initialize for simulating data
 */
void PsfImageDataSet::initializeForSim()
{
  if(!getParameter("size", false)->data_.hasValue_ ||
     !getParameter("npix", false)->data_.hasValue_)
    ThrowSimpleColorError("You must specify a size and npix for simulating data", "red");
  
  Distribution::Type type = getParameter("dist", false)->data_.hasValue() ? distType(getStringVal("dist")) : Distribution::DIST_POISS;
  
  if(type == Distribution::DIST_POISS && !getParameter("background", false)->data_.hasValue_)
    ThrowSimpleColorError("You must specify a background count rate if dist=poisson (i.e., " << name_ << ".background = xxx')", "red");

  if(type == Distribution::DIST_GAUSS && !getParameter("noiserms", false)->data_.hasValue_)
    ThrowSimpleColorError("You must specify a noise rms if dist=gauss (i.e., " << name_ << ".noiserms = xxx')", "red");

  Angle size;
  size.setDegrees(0.2);

  Image image(128, size);
  image.hasData_ = true;

  image = imp_.initializeImageParameters(image);

  initFromImage(image, image);
}

/**.......................................................................
 * Load data from a specified file
 */
void PsfImageDataSet::loadDataFromFile()
{
  //------------------------------------------------------------
  // Initialize our internal image from a FITS file
  //------------------------------------------------------------

  initializeFromFitsFile(getStringVal("file"));
}

/**.......................................................................
 * Initialize this data set from a fits file
 */
void PsfImageDataSet::initializeFromFitsFile(std::string fileName)
{
  //------------------------------------------------------------
  // Initialize our internal image from a FITS file.
  //------------------------------------------------------------

  Image preImage;
  Image image = imp_.initializeImage(fileName, preImage, &obs_);

  //------------------------------------------------------------
  // Call inheritor's initImage() method
  //------------------------------------------------------------

  initImage(fileName, image);

  //------------------------------------------------------------
  // And initialize our internal resources from the image
  //------------------------------------------------------------

  initFromImage(image, preImage);

  //------------------------------------------------------------
  // Finally, print what we know about the data
  //------------------------------------------------------------

  printFileStats();
}

/**.......................................................................
 * Initialize our internal data structures from the passed image
 */
void PsfImageDataSet::initFromImage(Image& image, Image& preImage)
{
  //------------------------------------------------------------
  // If this image has an absolute position, initialize our position
  // to match, for now (this can be overridden by user-specified
  // parameters)
  //------------------------------------------------------------

  if(!positionWasSpecified()) {
    setRa(image.getRa());
    setDec(image.getDec());
  }

  Antenna ant;
  Length diam;
  std::vector<Frequency> freqs(1);

  if(getParameter("psf", false)->data_.hasValue())
    psfType_ = PsfImageDataSet::psfType(getStringVal("psf"));

  if(psfType_ != PSF_NONE) {
    try {
      diam.setVal(obs_.getDoubleVal("ant.diameter"), obs_.getParameter("ant.diameter", true)->units_);
      ant.setDiameter(diam);
      freqs = obs_.getFrequencies();
    } catch(Exception& err) {
      ThrowSimpleColorError("Both an antenna diameter and a frequency/bandwidth must be specified to calculate the PSF.  Use, i.e.: " << std::endl << std::endl
			    << "  " << name_ << ".obs.ant.diameter = 10 m;" << std::endl
			    << "  " << name_ << ".obs.freqs = 90 GHz;"
			    << "  " << name_ << ".obs.bws = 1 GHz;"
			    , "red");
    }
  }

  try {
    freqs = obs_.getFrequencies();
  } catch(...) {
  }

  // For now, initialize to a single image

  data_.resize(1);
  data_[0].initialize(this, image, preImage, ant, freqs[0]);
}

//------------------------------------------------------------
// Overloaded Markov methods
//------------------------------------------------------------

void PsfImageDataSet::addModel(gcp::util::Model& model)
{
  //------------------------------------------------------------
  // Only proceed if this model applies to this data set
  //------------------------------------------------------------

  if(applies(model)) {
    for(unsigned i=0; i < data_.size(); i++)
      data_[i].addModel(model, dataSetType_);
  }
}

void PsfImageDataSet::clearModel()
{
  for(unsigned i=0; i < data_.size(); i++)
    data_[i].clearModel();
}

/**.......................................................................
 * Compute chisq over all images we are managing
 */
gcp::util::ChisqVariate PsfImageDataSet::computeChisq()
{
  //------------------------------------------------------------
  // Convolve models
  //------------------------------------------------------------

  convolveModels();

  //------------------------------------------------------------
  // Accumulate chi-squared over all frequencies
  //------------------------------------------------------------

  ChisqVariate chisq = accumulateChisq();

  //------------------------------------------------------------
  // Now clear models so that the next call to addModel()
  // re-initializes the model (rather than adds to it).
  //------------------------------------------------------------

  clearModel();

  //------------------------------------------------------------
  // Return the chi-squared
  //------------------------------------------------------------

  return chisq;
}

/**.......................................................................
 * Compute chisq over all images we are managing
 */
gcp::util::ChisqVariate PsfImageDataSet::computeChisq2()
{
  //------------------------------------------------------------
  // Convolve models
  //------------------------------------------------------------

  convolveModels();

  //------------------------------------------------------------
  // Accumulate chi-squared over all frequencies
  //------------------------------------------------------------

  ChisqVariate chisq = accumulateChisq2();

  //------------------------------------------------------------
  // Now clear models so that the next call to addModel()
  // re-initializes the model (rather than adds to it).
  //------------------------------------------------------------

  clearModel();

  //------------------------------------------------------------
  // Return the chi-squared
  //------------------------------------------------------------

  return chisq;
}

/**.......................................................................
 * Accumulate chisq over all images
 */
gcp::util::ChisqVariate PsfImageDataSet::accumulateChisq()
{
  ChisqVariate chisq;

  for(unsigned i=0; i < data_.size(); i++)
    chisq += data_[i].computeChisq();

  return chisq;
}

/**.......................................................................
 * Accumulate chisq over all images
 */
gcp::util::ChisqVariate PsfImageDataSet::accumulateChisq2()
{
  ChisqVariate chisq;

  for(unsigned i=0; i < data_.size(); i++)
    chisq += data_[i].computeChisq2();

  return chisq;
}

void PsfImageDataSet::convolveModels()
{
  for(unsigned i=0; i < data_.size(); i++)
    data_[i].convolveModel();
}

void PsfImageDataSet::display()
{
  setupForDisplay();

  zmin_ = data_[0].data_.min();
  zmax_ = data_[0].data_.max();

  if(getParameter("zmin", false)->data_.hasValue() && getParameter("zmax", false)->data_.hasValue()) {
    zmin_ = getDoubleVal("zmin");
    zmax_ = getDoubleVal("zmax");
  }

  PgUtil::setZmin(zmin_);
  PgUtil::setZmax(zmax_);

  for(unsigned i=0; i < data_.size(); i++)
    data_[i].display();
}

void PsfImageDataSet::displayResiduals() 
{
  setupForDisplay();

  PgUtil::setZmin(zmin_);
  PgUtil::setZmax(zmax_);

  for(unsigned i=0; i < data_.size(); i++)
    data_[i].displayResiduals();
}

void PsfImageDataSet::displayCompositeModel() 
{
  setupForDisplay();

  PgUtil::setZmin(zmin_);
  PgUtil::setZmax(zmax_);

  for(unsigned i=0; i < data_.size(); i++)
    data_[i].displayCompositeModel();
}

/**.......................................................................
 * Base-class method to initialize data that depend on absolute position
 */
void PsfImageDataSet::initializePositionDependentData()
{
  if(hasAbsolutePosition_) {
    for(unsigned i=0; i < data_.size(); i++) {
      data_[i].data_.setRaDec(ra_, dec_);
      data_[i].modelComponent_.setRaDec(ra_, dec_);
    }
  }
}

void PsfImageDataSet::simulateData(double sigma)
{
  for(unsigned i=0; i < data_.size(); i++) 
    data_[i].simulateData(sigma);
}

void PsfImageDataSet::writeCompositeModelToFile(std::string fileName, double sigma)
{
  data_[0].data_.writeToFitsFile(fileName);
}

//=======================================================================
// PdfImageData interface
//=======================================================================

PsfImageDataSet::PsfImageData::PsfImageData()
{
  hasNoiseRms_   = false;
  hasBackground_ = false;
  hasErrorVal_   = false;
  excMin_        = false;
  excMax_        = false;
  psfType_       = PSF_NONE;
  distType_      = Distribution::DIST_GAUSS;
}

PsfImageDataSet::PsfImageData::~PsfImageData()
{
}

/**.......................................................................
 * Initialize this object
 */
void PsfImageDataSet::PsfImageData::initialize(PsfImageDataSet* parent, Image& image, Image& preImage, Antenna& ant, Frequency& freq)
{
  parent_    = parent;
  antenna_   = ant;
  frequency_ = freq;
   
  //------------------------------------------------------------ 
  // See if a distribution or psf type were specified
  //------------------------------------------------------------ 

  if(parent_->getParameter("psf", false)->data_.hasValue())
    psfType_ = PsfImageDataSet::psfType(parent_->getStringVal("psf"));

  if(parent_->getParameter("dist", false)->data_.hasValue()) 
    distType_ = PsfImageDataSet::distType(parent_->getStringVal("dist"));

  //------------------------------------------------------------
  // Insert the image as our data set, padding if the image is not a
  // regular power of 2 grid
  //------------------------------------------------------------

  assignPaddedImageIfNecessary(image);

  //------------------------------------------------------------
  // Initialize the model components from the resulting image
  //------------------------------------------------------------

  initializeFromImage(compositeModel_, data_);
  initializeFromImage(modelComponent_, data_);
    
  //------------------------------------------------------------
  // If we are using a PSF, we also need to initialize the transform
  // arrays and primary beam.
  //
  // We calculate the primary beam for this
  // antenna type and frequency, and pre-compute the forward transform
  // of the beam.
  //------------------------------------------------------------
    
  if(psfType_ != PSF_NONE) {

    compositeModelDft_ = compositeModel_;
    
    initializeFromImage(primaryBeam_, data_);
    primaryBeamDft_ = primaryBeam_;

    std::cout << "\rComputing primary beam...                   ";
    fflush(stdout);
      
    if(psfType_ == PSF_GAUSS)
      primaryBeam_    = ant.getGaussianPrimaryBeam(primaryBeam_, freq);
    else
      primaryBeam_    = ant.getRealisticPrimaryBeam(primaryBeam_, freq);
      
    pbSum_ = primaryBeam_.sum();
      
    primaryBeamDft_.initialize(primaryBeam_);
    primaryBeamDft_.normalize(true);
    primaryBeamDft_.computeForwardTransform();
      
    std::cout << "\r                                             \r";

    if(parent_->debug_) {
      PgUtil::setWnad(true);
      primaryBeam_.display();
    }
  }

  //------------------------------------------------------------
  // If the parent has absolute position data, initialize it
  // in our images now
  //------------------------------------------------------------

  if(parent_->hasAbsolutePosition_) {
    data_.setRaDec(parent_->ra_, parent_->dec_);
    modelComponent_.setRaDec(parent_->ra_, parent_->dec_);
  }

  //------------------------------------------------------------
  // Initialize error estimates
  //------------------------------------------------------------

  initializeErrors(image, preImage);

  //------------------------------------------------------------
  // Also check if any exclusion radius was specified
  //------------------------------------------------------------

  if(parent_->getParameter("thetaMin", false)->data_.hasValue()) {
    thetaMin_.setVal(parent_->getDoubleVal("thetaMin"), parent_->getParameter("thetaMin", true)->units_);
    excMin_ = true;
  }

  if(parent_->getParameter("thetaMax", false)->data_.hasValue()) {
    thetaMax_.setVal(parent_->getDoubleVal("thetaMax"), parent_->getParameter("thetaMax", true)->units_);
    excMax_ = true;
  }

  //------------------------------------------------------------
  // Finally, assign the function for evaluating the likelihood
  //------------------------------------------------------------

  if(distType_ == Distribution::DIST_GAUSS) {

    if(hasErrorVal_) {
      lkFn_  = fixedErrorGauss;
      lkFn2_ = fixedErrorGauss;
    } else {
      lkFn_  = pixelErrorGauss;
      lkFn2_ = pixelErrorGauss;
    }

  } else {

    if(hasErrorVal_) {
      lkFn_  = fixedErrorPoiss;     // Use the real Poisson distribution for calculating likelihoods
      lkFn2_ = fixedErrorPoissTest; // Use the equivalent chi-squared distribution for quoting chi-squared
    } else {
      lkFn_  = pixelErrorPoiss;
      lkFn2_ = pixelErrorPoiss;
    }

  }
}

/**.......................................................................
 * Initialize errors now
 */
void PsfImageDataSet::PsfImageData::initializeErrors(Image& image, Image& preImage)
{
  //------------------------------------------------------------
  // Throw if we can't determine an error by any method
  //------------------------------------------------------------

  if(!parent_->getParameter("thetaMinErr", false)->data_.hasValue() && !parent_->getParameter("errRegion", false)->data_.hasValue()) {
    if(distType_ == Distribution::DIST_GAUSS && !parent_->getParameter("noiserms", false)->data_.hasValue())
      ThrowSimpleColorError("You must specify either a fixed noise rms via the 'noiserms' parameter, or use 'thetaMinErr' or 'errRegion' to estimate the error from the data", "red");
    if(distType_ == Distribution::DIST_GAUSS && !parent_->getParameter("background", false)->data_.hasValue())
      ThrowSimpleColorError("You must specify either a background count rate via the 'background' parameter, or use 'thetaMinErr' or 'errRegion' to estimate the error from the data", "red");

    if(distType_ == Distribution::DIST_POISS && !parent_->getParameter("background", false)->data_.hasValue())
      ThrowSimpleColorError("You must specify either a fixed background values via the 'background' parameter, or use 'thetaMinErr' or 'errRegion' to estimate the error from the data", "red");
  }

  //------------------------------------------------------------
  // Check if a minimum radius was specified for estimating the error.
  // If it was, estimate the error now.  The image we use is set by
  // the "errImage" parameter.
  //------------------------------------------------------------

  COUT("");
  if(parent_->getParameter("thetaMinErr", false)->data_.hasValue_) {

    if(parent_->getParameter("errRegion", false)->data_.hasValue_)
      ThrowSimpleColorError("You can use 'thetaMinErr' or 'errRegion' to estimate the error from the data, but not both", "red");

    thetaMinErr_.setVal(parent_->getDoubleVal("thetaMinErr"), parent_->getParameter("thetaMinErr", true)->units_);

    Image errImage;
    if(parent_->getStringVal("errImage") == "pre") {
      COUTCOLOR("Estimating noise from pixels > " << thetaMinErr_ << " from the center of the original image", "magenta");
      errImage = preImage;
    } else {
      COUTCOLOR("Estimating noise from pixels > " << thetaMinErr_ << " from the center of the extracted image", "magenta");
      errImage = data_;
    }

    Angle rmin, rmax;
    rmin.setVal(parent_->getDoubleVal("thetaMinErr"), parent_->getParameter("thetaMinErr", true)->units_);
    double xRad = errImage.xAxis().getAngularSize().radians();
    double yRad = errImage.yAxis().getAngularSize().radians();
    rmax.setRadians(2*sqrt(xRad*xRad + yRad*yRad));

    noiseRms_   = errImage.rms(rmin, rmax);
    background_ = errImage.mean(rmin, rmax);

    hasNoiseRms_   = true;
    hasBackground_ = true;
  }

  //------------------------------------------------------------
  // Check if a region was specified for estimating the error.
  //------------------------------------------------------------

  if(parent_->getParameter("errRegion", false)->data_.hasValue_) {

    if(parent_->getParameter("thetaMinErr", false)->data_.hasValue_)
      ThrowSimpleColorError("You can use 'thetaMinErr' or 'errRegion' to estimate the error from the data, but not both", "red");

    String region = parent_->getStringVal("errRegion");

    Image errImage;
    if(parent_->getStringVal("errImage") == "pre") {
      COUTCOLOR("Estimating noise from region " << region << " of the original image", "magenta");
      errImage = parent_->imp_.extractRegion(region, preImage);
    } else {
      COUTCOLOR("Estimating noise from region " << region << " of the extracted image", "magenta");
      errImage = parent_->imp_.extractRegion(region, image);
    }

    noiseRms_   = errImage.rms();
    background_ = errImage.mean();

    hasNoiseRms_   = true;
    hasBackground_ = true;
  }

  //------------------------------------------------------------
  // Print estimated values at this point
  //------------------------------------------------------------

  if(hasNoiseRms_) 
    COUTCOLOR("Estimating noise rms to be          " <<   noiseRms_, "magenta");

  if(hasBackground_)
    COUTCOLOR("Estimating background to be         " << background_, "magenta");
  
  if(distType_ == Distribution::DIST_POISS) {
    if(background_/(noiseRms_*noiseRms_) < 0.5)
      COUTCOLOR(std::endl << "Warning: this image has been processed in a way that is inconsistent with Poisson statistics "
		<< "(background counts are negative, or mean/var is significantly different than 1)."
		<< " Likelihood calculations will probably not give sensible results", "yellow");
  }

  //------------------------------------------------------------
  // Now check if a fixed error value was specified -- if any was, we
  // will override what we just estimated
  //------------------------------------------------------------

  // With Poisson statistics, we require an estimate of the background
  // count rate.  We can get this either from the image itself, or
  // from the 'background' keyword

  if(distType_ == Distribution::DIST_POISS && parent_->getParameter("background", false)->data_.hasValue()) {

    if(parent_->getParameter("thetaMinErr", false)->data_.hasValue() || parent_->getParameter("errRegion", false)->data_.hasValue())
      COUTCOLOR(std::endl << "Warning: a fixed value was specified for the 'background' parameter, and you have also specified 'thetaMinErr' or 'errRegion'.  The fixed value will override the estimate from the data.", "yellow");

    background_ = parent_->getDoubleVal("background");  
    hasBackground_ = true;
  }

  // With Gaussian statistics, we require an estimate of the
  // background count rate, and noise rms.  We can get this either
  // from the image itself, or from the 'background'and 'noiserms'
  // keywords.

  if(distType_ == Distribution::DIST_GAUSS) {

    if(parent_->getParameter("background", false)->data_.hasValue()) {
      if(parent_->getParameter("thetaMinErr", false)->data_.hasValue() || parent_->getParameter("errRegion", false)->data_.hasValue())
	COUTCOLOR(std::endl << "Warning: a fixed value was specified for the 'background' parameter, and you have also specified 'thetaMinErr' or 'errRegion'.  The fixed value will override the estimate from the data.", "yellow");
      
      background_ = parent_->getDoubleVal("background");  
      hasBackground_ = true;
    }
    
    if(parent_->getParameter("noiserms", false)->data_.hasValue()) {
      if(parent_->getParameter("thetaMinErr", false)->data_.hasValue() || parent_->getParameter("errRegion", false)->data_.hasValue())
	COUTCOLOR(std::endl << "Warning: a fixed value was specified for the 'noiserms' parameter, and you have also specified 'thetaMinErr' or 'errRegion'.  The fixed value will override the estimate from the data.", "yellow");
      
      noiseRms_   = parent_->getDoubleVal("noiserms");  
      hasNoiseRms_ = true;
    }
  }

  //------------------------------------------------------------
  // Finally, remove the background if we are using Gaussian
  // statistics.  Otherwise our likelihood calculations will be
  // meaningless
  //------------------------------------------------------------

  if(distType_ == Distribution::DIST_GAUSS) {
    data_ -= background_;
    hasErrorVal_ = hasBackground_ && hasNoiseRms_;
  } else {
    hasErrorVal_ = hasBackground_;
  }

  COUT("");
}

/**.......................................................................
 * Initialize our data image parameters from the passed image
 */
void PsfImageDataSet::PsfImageData::initializeFromImage(Image& imageToInit, Image& image)
{
  imageToInit = image;
  imageToInit.zero();
  imageToInit.hasData_ = true;
}

/**.......................................................................
 * Return an image that is a padded version of the original, if the
 * original image is not gridded to a power of 2 in each axis.
 */
void PsfImageDataSet::PsfImageData::assignPaddedImageIfNecessary(Image& image)
{
  data_ = image;

  iXStart_ = 0;
  iYStart_ = 0;
  nX_      = data_.xAxis().getNpix();
  nY_      = data_.xAxis().getNpix();

  //------------------------------------------------------------
  // If we are not using a psf, just return now, else compute the
  // nearest power of 2 grid needed to grid the data regularly
  //------------------------------------------------------------

  if(psfType_ == PSF_NONE)
    return;

  //------------------------------------------------------------
  // Figure out the size of this image
  //------------------------------------------------------------

  unsigned nx = image.xAxis().getNpix();
  unsigned ny = image.yAxis().getNpix();

  unsigned nxNear = Dft2d::nearestPowerOf2NotLessThan((double)nx);
  unsigned nyNear = Dft2d::nearestPowerOf2NotLessThan((double)ny);

  if(nx != nxNear || ny != nyNear) {

    COUTCOLOR("Note: padding " << nx << " x " << ny << " image to nearest power of 2 (" << nxNear << " x " << nyNear << ") since PSF is specified", "cyan");

    data_.xAxis().setNpix(nxNear);
    data_.yAxis().setNpix(nyNear);
    data_.resize();

    data_.zero();
    data_.assignDataFrom(image, nxNear/2-nx/2, nyNear/2-ny/2);

    iXStart_ = nxNear/2-nx/2;
    iYStart_ = nyNear/2-ny/2;
    nX_      = nx;
    nY_      = ny;
  }
}

/**.......................................................................
 * Add a model to this dataset
 */
void PsfImageDataSet::PsfImageData::addModel(gcp::util::Model& model, unsigned type)
{
  Generic2DAngularModel& model2D = dynamic_cast<Generic2DAngularModel&>(model);
    
  //------------------------------------------------------------
  // Now proceed with the calculation
  //------------------------------------------------------------
  
  model2D.fillImage(type, modelComponent_, &frequency_);
  
  if(!compositeModel_.hasData())
    compositeModel_.assignDataFrom(modelComponent_);
  else
    compositeModel_ += modelComponent_;
  
  //------------------------------------------------------------
  // And store the current value of the model offset
  //------------------------------------------------------------
  
  currentXoffRad_ = model2D.xOffset_.radians();
  currentYoffRad_ = model2D.yOffset_.radians();
}

/**.......................................................................
 * Generic method to compute chisq for a 2D data set (defaults to
 * image data set -- inheritors can override)
 *
 * We allow either a fixed error (hasErrorVal_ == true) to be
 * specified, or a per-pixel image (hasErrorVal_ == false).
 * 
 * Additionally, if an exclusion radius was specified, we exclude data
 * within this radius of the current model center.
 *
 * Finally, with image-plane models, we can have a situation where not
 * all pixels in the image are filled with valid data.  In this case,
 * we only iterate over valid pixels.
 */
gcp::util::ChisqVariate PsfImageDataSet::PsfImageData::computeChisq()
{
  ChisqVariate chisq;
  double lnLk = 0.0;
  unsigned nDof=0, ind;

  unsigned nx = data_.xAxis().getNpix();

  //------------------------------------------------------------
  // If we are excluding data, we must calculate the angular offsets
  // of the pixels and compare to the exclusion radius
  //------------------------------------------------------------

  if(excMin_ || excMax_) {

    unsigned nx     = data_.xAxis().getNpix();
    unsigned ny     = data_.yAxis().getNpix();

    int xSense      = data_.xAxis().getSense();
    int ySense      = data_.yAxis().getSense();

    double xResRad  = data_.xAxis().getAngularResolution().radians();
    double yResRad  = data_.yAxis().getAngularResolution().radians();

    unsigned ind;
    double r2;
    double rExcMin2 = thetaMin_.radians();
    double rExcMax2 = thetaMax_.radians();
    rExcMin2 *= rExcMin2;
    rExcMax2 *= rExcMax2;

    double xOffCenter, yOffCenter, dx, dy;

    for(unsigned iy=iYStart_; iy < nY_; iy++) {
      for(unsigned ix=iXStart_; ix < nX_; ix++) {

	ind = iy * nx + ix;

	if(!(data_.valid_[ind] && compositeModel_.valid_[ind]))
	  continue;

	//------------------------------------------------------------
	// Calculate the distance of this pixel from the current
	// model center
	//------------------------------------------------------------
	  
	xOffCenter = ((double)(ix) - double(nx)/2) * xResRad * xSense;
	yOffCenter = ((double)(iy) - double(ny)/2) * yResRad * ySense;

	dx         = xOffCenter - currentXoffRad_;
	dy         = yOffCenter - currentYoffRad_;

	r2 = dx*dx + dy*dy;

	if(!excMin_ || (r2 > rExcMin2)) {
	  if(!excMax_ || (r2 < rExcMax2)) {
	    lnLk += lkFn_(this, ind, nDof);
	  }
	}
      }
    }

    //------------------------------------------------------------
    // Else just include all data
    //------------------------------------------------------------
      
  } else {

    for(unsigned iy=iYStart_; iy < nY_; iy++) {
      for(unsigned ix=iXStart_; ix < nX_; ix++) {

	ind = iy * nx + ix;

	if(!(data_.valid_[ind] && compositeModel_.valid_[ind]))
	  continue;

	lnLk += lkFn_(this, ind, nDof);
      }
    }
  }

  chisq.setChisq(-2*lnLk, nDof);
  return chisq;
}

/**.......................................................................
 * Generic method to compute chisq for a 2D data set (defaults to
 * image data set -- inheritors can override)
 *
 * We allow either a fixed error (hasErrorVal_ == true) to be
 * specified, or a per-pixel image (hasErrorVal_ == false).
 * 
 * Additionally, if an exclusion radius was specified, we exclude data
 * within this radius of the current model center.
 *
 * Finally, with image-plane models, we can have a situation where not
 * all pixels in the image are filled with valid data.  In this case,
 * we only iterate over valid pixels.
 */
gcp::util::ChisqVariate PsfImageDataSet::PsfImageData::computeChisq2()
{
  ChisqVariate chisq;
  double lnLk = 0.0;
  unsigned nDof=0, ind;

  unsigned nx = data_.xAxis().getNpix();

  //------------------------------------------------------------
  // If we are excluding data, we must calculate the angular offsets
  // of the pixels and compare to the exclusion radius
  //------------------------------------------------------------

  if(excMin_ || excMax_) {

    unsigned nx     = data_.xAxis().getNpix();
    unsigned ny     = data_.yAxis().getNpix();

    int xSense      = data_.xAxis().getSense();
    int ySense      = data_.yAxis().getSense();

    double xResRad  = data_.xAxis().getAngularResolution().radians();
    double yResRad  = data_.yAxis().getAngularResolution().radians();

    unsigned ind;
    double r2;
    double rExcMin2 = thetaMin_.radians();
    double rExcMax2 = thetaMax_.radians();
    rExcMin2 *= rExcMin2;
    rExcMax2 *= rExcMax2;

    double xOffCenter, yOffCenter, dx, dy;

    for(unsigned iy=iYStart_; iy < nY_; iy++) {
      for(unsigned ix=iXStart_; ix < nX_; ix++) {

	ind = iy * nx + ix;

	if(!(data_.valid_[ind] && compositeModel_.valid_[ind]))
	  continue;

	//------------------------------------------------------------
	// Calculate the distance of this pixel from the current
	// model center
	//------------------------------------------------------------
	  
	xOffCenter = ((double)(ix) - double(nx)/2) * xResRad * xSense;
	yOffCenter = ((double)(iy) - double(ny)/2) * yResRad * ySense;

	dx         = xOffCenter - currentXoffRad_;
	dy         = yOffCenter - currentYoffRad_;

	r2 = dx*dx + dy*dy;

	if(!excMin_ || (r2 > rExcMin2)) {
	  if(!excMax_ || (r2 < rExcMax2)) {
	    lnLk += lkFn2_(this, ind, nDof);
	  }
	}
      }
    }

    //------------------------------------------------------------
    // Else just include all data
    //------------------------------------------------------------
      
  } else {

    for(unsigned iy=iYStart_; iy < nY_; iy++) {
      for(unsigned ix=iXStart_; ix < nX_; ix++) {

	ind = iy * nx + ix;

	if(!(data_.valid_[ind] && compositeModel_.valid_[ind]))
	  continue;

	lnLk += lkFn2_(this, ind, nDof);
      }
    }
  }

  chisq.setChisq(-2*lnLk, nDof);
  return chisq;
}

/**.......................................................................
 * Return ln Likelihood for gaussian data with fixed error
 */
PSF_LK_FN(PsfImageDataSet::PsfImageData::fixedErrorGauss)
{
  PsfImageData* psf = (PsfImageData*) obj;
  return psf->fixedErrorGauss(ind, nDof);
}

double PsfImageDataSet::PsfImageData::fixedErrorGauss(unsigned ind, unsigned& nDof)
{
  double val = (data_.data_[ind] - compositeModel_.data_[ind]) / noiseRms_;
  nDof += 1;
  return -val*val/2;
}

/**.......................................................................
 * Return ln Likelihood for gaussian data with per-pixel error
 */
PSF_LK_FN(PsfImageDataSet::PsfImageData::pixelErrorGauss)
{
  PsfImageData* psf = (PsfImageData*) obj;
  return psf->pixelErrorGauss(ind, nDof);
}

double PsfImageDataSet::PsfImageData::pixelErrorGauss(unsigned ind, unsigned& nDof)
{
  double val = (data_.data_[ind] - compositeModel_.data_[ind]) / error_.data_[ind];
  nDof += 1;
  return -val*val/2;
}

/**.......................................................................
 * Return ln Likelihood for poisson distribution with fixed background
 * count estimate
 */
PSF_LK_FN(PsfImageDataSet::PsfImageData::fixedErrorPoiss)
{
  PsfImageData* psf = (PsfImageData*) obj;
  return psf->fixedErrorPoiss(ind, nDof);
}

double PsfImageDataSet::PsfImageData::fixedErrorPoiss(unsigned ind, unsigned& nDof)
{
  double expectedCountRate = compositeModel_.data_[ind] + background_;
  nDof += 1;
  return Sampler::lnPoissPdf((unsigned)(data_.data_[ind]), expectedCountRate);
}

/**.......................................................................
 * Return ln Likelihood for Poisson distribution, using the gaussian
 * approximation
 */
PSF_LK_FN(PsfImageDataSet::PsfImageData::fixedErrorGaussApprox)
{
  PsfImageData* psf = (PsfImageData*) obj;
  return psf->fixedErrorGaussApprox(ind, nDof);
}

double PsfImageDataSet::PsfImageData::fixedErrorGaussApprox(unsigned ind, unsigned& nDof)
{
  double expectedCountRate = compositeModel_.data_[ind] + background_;
  nDof += 1;
  double val = expectedCountRate - data_.data_[ind];
  return -val*val/(2*expectedCountRate);
}

/**.......................................................................
 * Return ln Likelihood for poisson distribution with per-pixel background
 * count estimates
 */
PSF_LK_FN(PsfImageDataSet::PsfImageData::pixelErrorPoiss)
{
  PsfImageData* psf = (PsfImageData*) obj;
  return psf->pixelErrorPoiss(ind, nDof);
}

double PsfImageDataSet::PsfImageData::pixelErrorPoiss(unsigned ind, unsigned& nDof)
{
  double expectedCountRate = compositeModel_.data_[ind] + error_.data_[ind];
  nDof += 1;
  return Sampler::lnPoissPdf((unsigned)(data_.data_[ind]), expectedCountRate);
}

/**.......................................................................
 * Return ln Likelihood for poisson distribution with fixed background
 * counts, using the equivalent chi-squared estimate.
 *
 * This method notes that the Poisson likelihood for n given <n> and
 * chi-squared likelihood are equivalent if we identify:
 *
 *    chisq = 2*<n>  
 *
 * and
 *
 *    nDof = 2*(n + 1),
 *
 * whence lnLk = -chisq/2 = -<n>
 *
 */
PSF_LK_FN(PsfImageDataSet::PsfImageData::fixedErrorPoissTest)
{
  PsfImageData* psf = (PsfImageData*) obj;
  return psf->fixedErrorPoissTest(ind, nDof);
}

double PsfImageDataSet::PsfImageData::fixedErrorPoissTest(unsigned ind, unsigned& nDof)
{
  double expectedCountRate = compositeModel_.data_[ind] + background_;
  nDof += 2*(unsigned)(data_.data_[ind]) + 2;
  return -expectedCountRate;
}

void PsfImageDataSet::PsfImageData::clearModel()
{
  compositeModel_.hasData_ = false;
}

void PsfImageDataSet::setupForDisplay()
{
  //------------------------------------------------------------
  // Always adjust the window aspect ratio
  //------------------------------------------------------------

  PgUtil::setWnad(true);

  //------------------------------------------------------------
  // If a colormap was specified, set it now
  //------------------------------------------------------------

  if(getParameter("cmap", false)->data_.hasValue()) {
    PgUtil::setColormap(getStringVal("cmap"));
  } else {
    PgUtil::setColormap("grey");
  }
}

void PsfImageDataSet::PsfImageData::display()
{
  PgUtil::useHeader(true);
  PgUtil::setHeader(parent_->displayHeaderString("Data"), PgUtil::JUST_LEFT);

  data_.display();

  PgUtil::useHeader(false);

  if(parent_->getParameter("dataimage", false)->data_.hasValue()) {
    data_.writeToFitsFile(parent_->getStringVal("dataimage"));
  }
}

void PsfImageDataSet::PsfImageData::displayCompositeModel()
{
  PgUtil::useHeader(true);
  PgUtil::setHeader(parent_->displayHeaderString("Model"), PgUtil::JUST_LEFT);

  compositeModel_.hasData_ = true;
  compositeModel_.display();

  if(parent_->getParameter("modelimage", false)->data_.hasValue()) {
    compositeModel_.writeToFitsFile(parent_->getStringVal("modelimage"));
  }
}

void PsfImageDataSet::PsfImageData::displayResiduals()
{
  PgUtil::useHeader(true);
  PgUtil::setHeader(parent_->displayHeaderString("Residuals"), PgUtil::JUST_LEFT);

  Image res = data_ - compositeModel_;
  res.display();

  if(parent_->getParameter("resimage", false)->data_.hasValue()) {
    res.writeToFitsFile(parent_->getStringVal("resimage"));
  }
}

void PsfImageDataSet::PsfImageData::simulateData(double sigma)
{
#if 1
  displayCompositeModel();
#endif
  if(distType_ == Distribution::DIST_GAUSS) {
    COUT("Generating Gaussian data");
    for(unsigned i=0; i < data_.data_.size(); i++) {
      data_.data_[i] = compositeModel_.data_[i] + Sampler::generateGaussianSample(noiseRms_) + background_;
    }
  } else if(distType_ == Distribution::DIST_POISS) {
    COUT("Generating Poisson data with background = " << background_);
    for(unsigned i=0; i < data_.data_.size(); i++) {
      unsigned expectedCountRate = (unsigned) (compositeModel_.data_[i] + background_);
      if(expectedCountRate > 0)
	data_.data_[i] = Sampler::generatePoissonSample(expectedCountRate);
    }
  }
#if 1
#endif
}

void PsfImageDataSet::PsfImageData::convolveModel()
{
  if(psfType_ == PSF_NONE)
    return;

  compositeModelDft_.initialize(compositeModel_);
  compositeModelDft_.normalize(true);
  compositeModelDft_.computeForwardTransform();

  compositeModelDft_.complexMultiply(primaryBeamDft_, false);
  compositeModelDft_.shift();

  compositeModelDft_.computeInverseTransform();
  compositeModel_ = compositeModelDft_.getImage(false) / pbSum_;
}

PsfImageDataSet::PsfType PsfImageDataSet::psfType(std::string type)
{
  String typeStr(type);
  typeStr = typeStr.toLower();

  if(typeStr.contains("none"))
    return PSF_NONE;
  else if(typeStr.contains("real"))
    return PSF_REAL;
  if(typeStr.contains("gauss"))
    return PSF_GAUSS;

  ThrowSimpleColorError("Unknown type of psf: " << type, "red");

  return PSF_UNKNOWN;
}

Distribution::Type PsfImageDataSet::distType(std::string type)
{
  String typeStr(type);
  typeStr = typeStr.toLower();

  if(typeStr.contains("gauss"))
    return Distribution::DIST_GAUSS;
  else if(typeStr.contains("poiss"))
    return Distribution::DIST_POISS;

  ThrowSimpleColorError("Unsupported distribution type: " << type, "red");

  return Distribution::DIST_UNKNOWN;
}

std::string PsfImageDataSet::distToString(Distribution::Type type)
{
  switch(type) {
  case Distribution::DIST_GAUSS:
    return "Gaussian";
  case Distribution::DIST_POISS:
    return "Poisson";
  default:
    return "Unknown";
    break;
  }
}

void PsfImageDataSet::finalizeForDisplay()
{
  convolveModels();
}

void PsfImageDataSet::printFileStats()
{
  std::ostringstream os;
  os << "FITS File contains observations of: " << std::endl << std::endl;
    
  try {
    os << "  Object     = " << obs_.getSourceName() << std::endl;
  } catch(...) {
    os << "  Object     = UNKNOWN" << std::endl;
  }

  try {
    os << "  RA         = " << " " << obs_.getObsRa()      << std::endl;
  } catch(...) {
    os << "  RA         = UNKNOWN"  << std::endl;
  }

  try {
    os << "  DEC        = " << obs_.getObsDec()     << std::endl ;
  } catch(...) {
    os << "  DEC        = UNKNOWN"  << std::endl;
  }

  try {
    os << "  Equinox    = " << obs_.getObsEquinox() << std::endl << std::endl;
  } catch(...) {
    os << "  Equinox    = UNKNOWN" << std::endl << std::endl;
  }

  os << "with: " << std::endl << std::endl;

  try {
    os << "  Telescope  = " << obs_.getTelescopeName()  << std::endl;
  } catch(...) {
    os << "  Telescope  = UNKNOWN"  << std::endl;
  }

  try {
    os << "  Instrument = " << obs_.getInstrumentName() << std::endl << std::endl;;
  } catch(...) {
    os << "  Instrument = UNKNOWN" << std::endl << std::endl;;
  }
  
  COUTCOLOR(std::endl << os.str(), "cyan");
}

std::string PsfImageDataSet::displayHeaderString(std::string type)
{
  std::ostringstream os;
  String source = obs_.getSourceName();
  source.strip(' ');
  source.strip('\r');
  source.strip('\n');
  source.strip('"');
  source.strip('\0');

  os << type << " (" << source << ")";

  return os.str();
}
