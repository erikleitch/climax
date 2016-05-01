#include "gcp/datasets/ImageDataSet.h"

#include "gcp/fftutil/Generic2DAngularModel.h"

using namespace std;

using namespace gcp::datasets;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ImageDataSet::ImageDataSet() 
{
  hasErrorVal_         = false;
  excData_             = false;
  errScaleFactor_      = 1.0;

  addParameter("frequency", DataType::DOUBLE, "When simulating, the frequency at which to calculate the model");
  addParameter("error",     DataType::DOUBLE, "Noise sigma to use");
  addParameter("thetaExc",  DataType::DOUBLE, "If specified, data within thetaExc of the current model center will be excluded from the fit");

  addParameter(imageManager_);
}

/**.......................................................................
 * Destructor.
 */
ImageDataSet::~ImageDataSet() {}

void ImageDataSet::initializeErrors()
{
  if(getParameter("error", false)->data_.hasValue_) {
    hasErrorVal_ = true;
    errorVal_    = getDoubleVal("error") * errScaleFactor_;
  } else {
    ThrowSimpleError("No error value was defined for this data set");
  }
}

void ImageDataSet::loadData(bool simulate)
{
  if(simulate) {
    initializeForSim();
  } else {
    loadDataFromFile();
  }

  printFileStats();
}

void ImageDataSet::initializeForSim()
{
  if(!getParameter("npix", false)->data_.hasValue() || !getParameter("size", false)->data_.hasValue())
    ThrowColorError("You must specify both an image size ('" << name_ << ".size') and resolution ('" << name_ << ".npix') "
		    "when simulating data", "red");

  data_.setNpix(getUintVal("npix"));

  Angle size;
  size.setVal(getDoubleVal("size"), getParameter("size", true)->units_);

  data_.setAngularSize(size);

  //------------------------------------------------------------
  // Initialize model images to match the data size
  //------------------------------------------------------------

  modelComponent_ = data_;
  compositeModel_ = data_;

  //------------------------------------------------------------
  // Initialize the frequency at which we want to calculate the model
  //------------------------------------------------------------

  if(!getParameter("frequency", false)->data_.hasValue())
    ThrowColorError("You must specify the frequency ('" << name_ << ".frequency')"
		    "when simulating data", "red");

  frequency_.setVal(getDoubleVal("frequency"), getParameter("frequency", true)->units_);

  data_.setFrequency(frequency_);
  modelComponent_.setFrequency(frequency_);
  compositeModel_.setFrequency(frequency_);
}

/**.......................................................................
 * Load data from a specified file
 */
void ImageDataSet::loadDataFromFile()
{
  //------------------------------------------------------------
  // Initialize our internal image from a FITS file
  //------------------------------------------------------------

  initializeFromFitsFile(getStringVal("file"));

  //------------------------------------------------------------
  // Initialize errors for this data set
  //------------------------------------------------------------

  initializeErrors();

  //------------------------------------------------------------
  // Check if any absolute position was specified in the file.  If so,
  // set the position of this dataset to be the center position of the
  // data image.  But don't override a position that was set by the user
  //------------------------------------------------------------

  if(data_.hasAbsolutePosition_) {
    if(!getParameter("ra", false)->data_.hasValue())
      setRa(data_.getRa());

    if(!getParameter("dec", false)->data_.hasValue())
      setDec(data_.getDec());
  }

  //------------------------------------------------------------
  // Also check if any exclusion radius was specified
  //------------------------------------------------------------

  if(getParameter("thetaExc", false)->data_.hasValue()) {
    thetaExc_.setVal(getDoubleVal("thetaExc"), getParameter("thetaExc", true)->units_);
    excData_ = true;
  }
}

/**.......................................................................
 * Initialize this data set from a fits file
 */
void ImageDataSet::initializeFromFitsFile(std::string fileName)
{
  //------------------------------------------------------------
  // And initialize our internal image from a FITS file.
  //------------------------------------------------------------

  data_ = imageManager_.initializeImage(fileName, &obs_);


  //------------------------------------------------------------
  // Initialize model images to match the data size
  //------------------------------------------------------------

  modelComponent_ = data_;
  compositeModel_ = data_;
}

void ImageDataSet::clearModel()
{
  compositeModel_.hasData_ = false;
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
 */
gcp::util::ChisqVariate ImageDataSet::computeChisq()
{
  ChisqVariate chisq;
  double val;

  if(hasErrorVal_) {

    //------------------------------------------------------------
    // If we are excluding data, we must calculate the angular offsets
    // of the pixels and compare to the exclusion radius
    //------------------------------------------------------------

    if(excData_) {

      unsigned nx     = data_.xAxis().getNpix();
      unsigned ny     = data_.yAxis().getNpix();

      int xSense      = data_.xAxis().getSense();
      int ySense      = data_.yAxis().getSense();

      double xResRad  = data_.xAxis().getAngularResolution().radians();
      double yResRad  = data_.yAxis().getAngularResolution().radians();

      unsigned ind;
      double r2, rExc2 = thetaExc_.radians();
      rExc2 *= rExc2;
      double xOffCenter, yOffCenter, dx, dy;

      for(unsigned iy=0; iy < ny; iy++) {
	for(unsigned ix=0; ix < nx; ix++) {
	  ind = iy * nx + ix;

	  //------------------------------------------------------------
	  // Calculate the distance of this pixel from the current
	  // model center
	  //------------------------------------------------------------
	  
	  xOffCenter = ((double)(ix) - double(nx)/2 + 0.5) * xResRad * xSense;
	  yOffCenter = ((double)(iy) - double(ny)/2 + 0.5) * yResRad * ySense;

	  dx         = xOffCenter - currentXoffRad_;
	  dy         = yOffCenter - currentYoffRad_;

	  r2 = dx*dx + dy+dy;

	  if(r2 > rExc2) {
	    val = (data_.data_[ind] - compositeModel_.data_[ind]) / errorVal_;
	    chisq += val*val;
	  }
	}
      }

      //------------------------------------------------------------
      // Else just include all data
      //------------------------------------------------------------
      
    } else {

      for(unsigned iDat=0; iDat < data_.data_.size(); iDat++) {
	val = (data_.data_[iDat] - compositeModel_.data_[iDat]) / errorVal_;
	chisq += val*val;
      }
      
    }

    //------------------------------------------------------------
    // Else use an error 'image' (per-pixel errors)
    //------------------------------------------------------------

  } else {

    //------------------------------------------------------------
    // If we are excluding data, we must calculate the angular offsets
    // of the pixels and compare to the exclusion radius
    //------------------------------------------------------------

    if(excData_) {

      unsigned nx    = data_.xAxis().getNpix();
      unsigned ny    = data_.yAxis().getNpix();

      int     xSense = data_.xAxis().getSense();
      int     ySense = data_.yAxis().getSense();

      double xResRad = data_.xAxis().getAngularResolution().radians();
      double yResRad = data_.yAxis().getAngularResolution().radians();

      unsigned ind;
      double r2, rExc2 = thetaExc_.radians();
      rExc2 *= rExc2;
      double xOffCenter, yOffCenter, dx, dy;

      for(unsigned iy=0; iy < ny; iy++) {
	for(unsigned ix=0; ix < nx; ix++) {
	  ind = iy * nx + ix;

	  //------------------------------------------------------------
	  // Calculate the distance of this pixel from the current
	  // model center
	  //------------------------------------------------------------
	  
	  xOffCenter = ((double)(ix) - double(nx)/2 + 0.5) * xResRad * xSense;
	  yOffCenter = ((double)(iy) - double(ny)/2 + 0.5) * yResRad * ySense;

	  dx         = xOffCenter - currentXoffRad_;
	  dy         = yOffCenter - currentYoffRad_;

	  r2 = dx*dx + dy+dy;

	  if(r2 > rExc2) {
	    val = (data_.data_[ind] - compositeModel_.data_[ind]) / error_.data_[ind];
	    chisq += val*val;
	  }
	}
      }

      //------------------------------------------------------------
      // Else just include all data
      //------------------------------------------------------------
      
    } else {

      for(unsigned iDat=0; iDat < data_.data_.size(); iDat++) {
	val = (data_.data_[iDat] - compositeModel_.data_[iDat]) / error_.data_[iDat];
	chisq += val*val;
      }

    }

  }

  return chisq;
}

void ImageDataSet::display()
{
  data_.display();

  if(getParameter("dataimage", false)->data_.hasValue()) {
    data_.writeToFitsFile(getStringVal("dataimage"));
  }
}

void ImageDataSet::displayCompositeModel()
{
  compositeModel_.display();

  if(getParameter("modelimage", false)->data_.hasValue()) {
    compositeModel_.writeToFitsFile(getStringVal("modelimage"));
  }
}

void ImageDataSet::displayResiduals()
{
  PgUtil::setZmin(data_.min());
  PgUtil::setZmax(data_.max());

  Image res = data_ - compositeModel_;
  res.display();

  PgUtil::setZmin(0.0);
  PgUtil::setZmax(0.0);

  if(getParameter("resimage", false)->data_.hasValue()) {
    res.writeToFitsFile(getStringVal("resimage"));
  }
}

/**.......................................................................
 * Add a model to this dataset
 */
void ImageDataSet::addModel(gcp::util::Model& model)
{
  //------------------------------------------------------------
  // Only proceed if this model applies to this data set
  //------------------------------------------------------------

  if(applies(model)) {

    Generic2DAngularModel& model2D = dynamic_cast<Generic2DAngularModel&>(model);
    
    //------------------------------------------------------------
    // Now proceed with the calculation
    //------------------------------------------------------------
    
    COUT("MODEL (1) compon refpix = " << modelComponent_.raRefPix_ << " " << modelComponent_.decRefPix_);

    model2D.fillImage(dataSetType_, modelComponent_, &frequency_);
    
    if(!compositeModel_.hasData()) {
      if(model2D.modelType_ & ModelType::MODEL_ADDITIVE)
	compositeModel_.assignDataFrom(modelComponent_);
      else {
	ThrowSimpleColorError("Can't assign a multiplicative model without any additive models defined", "red");
      }
    } else {
      if(model2D.modelType_ & ModelType::MODEL_ADDITIVE)
	compositeModel_ += modelComponent_;
      else {
	compositeModel_ *= modelComponent_;
      }
    }

    // And store the current value of the model offset
    
    currentXoffRad_ = model2D.xOffset_.radians();
    currentYoffRad_ = model2D.yOffset_.radians();
  }

#if 0
  COUT("Added model");
  compositeModel_.display();
#endif
}

/**.......................................................................
 * Base-class method to initialize data that depend on absolute position
 */
void ImageDataSet::initializePositionDependentData()
{
  if(hasAbsolutePosition_) {
    modelComponent_.setRaDec(ra_, dec_);
    compositeModel_.setRaDec(ra_, dec_);
  }
}

void ImageDataSet::simulateData(double sigma)
{
  data_ = compositeModel_;

  double noiseRms = obs_.getFixedNoiseRms(data_.getUnits());

  for(unsigned i=0; i < data_.data_.size(); i++) {
    data_.data_[i] += Sampler::generateGaussianSample(noiseRms);
  }

}

void ImageDataSet::writeCompositeModelToFile(std::string file, double sigma)
{
  data_.writeToFitsFile(file);
}

/**.......................................................................
 * Print information about the file just read in
 */
void ImageDataSet::printFileStats()
{
  std::ostringstream os;
  FORMATCOLOR(os, "FITS File contains observations of: " << std::endl << std::endl
	      << "  Object     = " << obs_.getSourceName() << std::endl
	      << "  RA         = " << " " << obs_.getObsRa()      << std::endl
	      << "  Dec        = " << obs_.getObsDec()     << std::endl, "cyan");

  try {
    os << "  Equinox    = " << obs_.getObsEquinox() << std::endl << std::endl;
  } catch(...) {
    FORMATCOLOR(os, "  Equinox    = 2000 (assumed -- not specified in file)" << std::endl, "red");
  }

  FORMATCOLOR(os, "with: " << std::endl << std::endl
	      << "  Telescope  = " << obs_.getTelescopeName()  << std::endl
	      << "  Instrument = " << obs_.getInstrumentName() << std::endl, "cyan");

  COUT(os.str());
}
