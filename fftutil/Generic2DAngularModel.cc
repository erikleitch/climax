#include "gcp/fftutil/Generic2DAngularModel.h"
#include "gcp/fftutil/Image.h"
#include "gcp/fftutil/UvDataGridder.h"

#include "gcp/util/Astrometry.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::util;

//#define TIMER_TEST

#ifdef TIMER_TEST
#include "gcp/util/Timer.h"
Timer fit1, fit2, fit3, fit4;
double ft1=0.0, ft2=0.0, ft3=0.0, ft4=0.0;
#endif

/**.......................................................................
 * Constructor.
 */
Generic2DAngularModel::Generic2DAngularModel(ThreadPool* pool) :
  Model(pool)
{
  initialize();
}

/**.......................................................................
 * Destructor.
 */
Generic2DAngularModel::~Generic2DAngularModel() 
{
  // Delete any thread pool resources that were allocated

  for(unsigned iThread=0; iThread < execData_.size(); iThread++) {
    delete execData_[iThread];
  }
}

/**.......................................................................
 * Generic initialization method for all 2D models
 */
void Generic2DAngularModel::initialize()
{
  //------------------------------------------------------------
  // Set flags appropriate to this model type
  //------------------------------------------------------------

  dataSetType_ = 
    DataSetType::DATASET_2D |
    DataSetType::DATASET_GENERIC | 
    DataSetType::DATASET_RADIO | 
    DataSetType::DATASET_XRAY_IMAGE;

  execData_.resize(0);

  //------------------------------------------------------------
  // Add components common to all 2D angular models
  //------------------------------------------------------------

  addComponent(ra_);
  addComponent(dec_);
  addComponent(xOffset_);
  addComponent(yOffset_);
  addComponent(rotationAngle_);

  // Different normalizations for different dataset types

  addComponent(normalization_);      // Generic
  addComponent(radioNormalization_); // Radio
  addComponent(xrayNormalization_);  // Xray

  addComponent(spectralType_);
  addComponent(spectralIndex_);
  addComponent(normalizationFrequency_);
  addComponent(mu_);
  addComponent(mue_);
  addComponent(fGas_);

  //------------------------------------------------------------
  // Add handles for the various components, so that they can be
  // accessed by name.
  //------------------------------------------------------------

  addComponentName(ra_,                     "ra",                     "The center RA of this model: 'hh:mm:ss.ss'");
  addComponentName(dec_,                    "dec",                    "The center Dec of this model: 'dd.mm.ss.ss'");
  addComponentName(xOffset_,                "xoff",                   "The x-offset of the model center");
  addComponentName(yOffset_,                "yoff",                   "The y-offset of the model center");
  addComponentName(rotationAngle_,          "rotang",                 "A rotation angle for this model");
  addComponentName(normalization_,          "norm",                   "A normalization for the model");
  addComponentName(normalization_,          "normalization",          "A normalization for the model.  If fitting to a generic dataset, "
		                                                      "you will be required to specify this.");
  addComponentName(radioNormalization_,     "Sradio",                 "A radio normalization for the model.  If fitting to a radio dataset, "
		                                                      "you will be required to specify this.");
  addComponentName(xrayNormalization_,      "Sxray",                  "An xray normalization for the model.  If fitting to an xray dataset, "
		                                                      "you will be required to specify this.");
  addComponentName(spectralType_,           "spectralType",           "The spectral type of this model ('none', 'alpha', or 'sz')");
  addComponentName(spectralIndex_,          "spectralIndex",          "The spectral index, if spectralType = alpha");
  addComponentName(normalizationFrequency_, "normalizationFrequency", "The frequency at which the normalization applies");
  addComponentName(fGas_,                   "fgas",                   "The gas-mass fraction");
  addComponentName(mu_,                     "mu",                     "The mean molecular weight");
  addComponentName(mue_,                    "mue",                    "The mean molecular weight per free electron");

  ra_.hasValue_        = false;
  dec_.hasValue_       = false;
  hasAbsolutePosition_ = false;
  positionChecked_     = false;

  //------------------------------------------------------------
  // Allow the normalization to be specified without units (inheritors
  // can override this)
  //------------------------------------------------------------

  normalization_.allowUnitless(true);
  normalization_.addConversion("comptony", 1.0);
  normalization_.addConversion("muK",      1.0);
  normalization_.addConversion("Jy",       1.0);

  radioNormalization_.allowUnitless(true);
  radioNormalization_.addConversion("comptony", 1.0);
  radioNormalization_.addConversion("muK",      1.0);
  radioNormalization_.addConversion("Jy",       1.0);
  radioNormalization_.addConversion("Jy/bm",    1.0);
  radioNormalization_.addConversion("mJy",      1.0);
  radioNormalization_.addConversion("mJy/sr",   1.0);
  radioNormalization_.addConversion("MJy/sr",   1.0);
  radioNormalization_.addConversion("mK",       1.0);
  radioNormalization_.addConversion("K",        1.0);

  xrayNormalization_.allowUnitless(true);
  xrayNormalization_.addConversion("MJy/sr",    1.0);
  xrayNormalization_.addConversion("counts",    1.0);

  //------------------------------------------------------------
  // Allow the spectral index to be specified without units
  //------------------------------------------------------------

  spectralIndex_.allowUnitless(true);

  //------------------------------------------------------------
  // Frequency at which the source flux is specified is not allowed to
  // vary
  //------------------------------------------------------------

  normalizationFrequency_.allowToVary(false);

  //------------------------------------------------------------
  // Allow the gas-mass fraction to be specified without units, and
  // default it to something sensible
  //------------------------------------------------------------

  fGas_.allowUnitless(true);
  fGas_.setVal(0.175, "");

  //------------------------------------------------------------
  // Allow the mean molecular weight to be specified without units, and
  // default it to something sensible
  //------------------------------------------------------------

  mu_.allowUnitless(true);
  mu_.setVal(0.59, "");

  //------------------------------------------------------------
  // Allow the mean molecular weight per free electron to be specified
  // without units, and default it to something sensible
  //------------------------------------------------------------

  mue_.allowUnitless(true);
  mue_.setVal(1.14, "");

  //------------------------------------------------------------
  // Add parameters common to all 2D models
  //------------------------------------------------------------

  addParameter("innerRadius",         DataType::DOUBLE, "The inner physical radius (ie, Mpc) from which integrations will be performed");
  addParameter("outerRadius",         DataType::DOUBLE, "The outer physical radius (ie, Mpc) to which integrations will be performed");
  addParameter("electronTemperature", DataType::DOUBLE, "Electron temperature, where relevant");

  //------------------------------------------------------------
  // Initialize all base-class components as fixed.  
  //------------------------------------------------------------

  initializeComponentsToFixed();
}

/**.......................................................................
 * Set a right ascension for this model
 */
void Generic2DAngularModel::setRa(HourAngle ra)
{
  getVar("ra")->setVal(ra.hours(), "hours");
  hasAbsolutePosition_ = true;
}

/**.......................................................................
 * Set a declination for this model
 */
void Generic2DAngularModel::setDec(Declination dec)
{
  getVar("dec")->setVal(dec.degrees(), "degrees");
  hasAbsolutePosition_ = true;
}

/**.......................................................................
 * Set an X offset for this model
 */
void Generic2DAngularModel::setXOffset(Angle xOffset)
{
  xOffset_ = xOffset;
}

/**.......................................................................
 * Set a Y offset for this model
 */
void Generic2DAngularModel::setYOffset(Angle yOffset)
{
  yOffset_ = yOffset;
}

/**.......................................................................
 * Set a rotation angle for this model
 */
void Generic2DAngularModel::setRotationAngle(Angle rotationAngle)
{
  rotationAngle_ = rotationAngle;
}

/**.......................................................................
 * Envelope() just calls the appropriate method for this dataset
 */
double Generic2DAngularModel::envelope(unsigned type, double xRad, double yRad)
{
  if(type & DataSetType::DATASET_RADIO) {
    return radioEnvelope(xRad, yRad);
  } else if(type & DataSetType::DATASET_XRAY_IMAGE) {
    return xrayImageEnvelope(xRad, yRad);
  } else if(type & DataSetType::DATASET_GENERIC) {
    return genericEnvelope(xRad, yRad);
  } else {
    ThrowError("Unsupported dataset type: " << type);
    return 0.0;
  }
}

/**.......................................................................
 * Stub-out envelope in the base-class so that calling from undefined
 * inheritors will throw an error
 */
double Generic2DAngularModel::radioEnvelope(double xRad, double yRad)
{
  ThrowError("No radioEnvelope(double xRad, double yRad) method has been defined for this inherited class");
  return 0.0;
}

double Generic2DAngularModel::xrayImageEnvelope(double xRad, double yRad)
{
  ThrowError("No xrayImageEnvelope(double xRad, double yRad) method has been defined for this inherited class");
  return 0.0;
}

double Generic2DAngularModel::genericEnvelope(double xRad, double yRad)
{
  ThrowError("No genericEnvelope(double xRad, double yRad) method has been defined for this inherited class");
  return 0.0;
}

/**.......................................................................
 * Envelope() just calls the appropriate method for this dataset
 */
double Generic2DAngularModel::envelope(void* evalData, unsigned type, double xRad, double yRad)
{
  if(type & DataSetType::DATASET_RADIO) {
    return radioEnvelope(evalData, xRad, yRad);
  } else if(type & DataSetType::DATASET_XRAY_IMAGE) {
    return xrayImageEnvelope(evalData, xRad, yRad);
  } else if(type & DataSetType::DATASET_GENERIC) {
    return genericEnvelope(evalData, xRad, yRad);
  } else {
    ThrowError("Unsupported dataset type: " << type);
    return 0.0;
  }
}

/**.......................................................................
 * Stub-out envelope in the base-class so that calling from undefined
 * inheritors will throw an error
 */
double Generic2DAngularModel::radioEnvelope(void* evalData, double xRad, double yRad)
{
  ThrowError("No radioEnvelope(void* evalData) method has been defined for this inherited class");
  return 0.0;
}

double Generic2DAngularModel::xrayImageEnvelope(void* evalData, double xRad, double yRad)
{
  ThrowError("No xrayImageEnvelope(void* evalData) method has been defined for this inherited class");
  return 0.0;
}

double Generic2DAngularModel::genericEnvelope(void* evalData, double xRad, double yRad)
{
  ThrowError("No genericEnvelope(void* evalData) method has been defined for this inherited class");
  return 0.0;
}

//=======================================================================
// Methods to fill an image
//=======================================================================

/**.......................................................................
 * Fill an image with this model
 */
void Generic2DAngularModel::fillImage(unsigned type, Image& image, void* params)
{
  //------------------------------------------------------------
  // If no thread pool exists, just call for single thread
  //------------------------------------------------------------

  if(!pool_ || execData_.size() == 0) {

    fillImageSingleThread(type, image, params);
    
    //------------------------------------------------------------
    // Else split the image into segments
    //------------------------------------------------------------

  } else {

    unsigned nThreadExec = initializeExecData(type, image, params);

    synchronizer_.reset(nThreadExec);

    for(unsigned iThread=0; iThread < nThreadExec; iThread++) {

      ExecData* ed = execData_[iThread];

      //------------------------------------------------------------
      // Now execute it
      //------------------------------------------------------------

      synchronizer_.registerPending(iThread);
      pool_->execute(&execFillImageMultiThread, ed);
    }

    synchronizer_.wait();
  }

  image.setHasData(true);
  image.setUnits(units_);

#if 0
  Image image2 = image, delta;
  fillImageSingleThread(type, image2, params);
  
  delta = image - image2;
  PgUtil::setInteractive(true);
  delta.display();
#endif

}

/**.......................................................................
 * Initialize all data needed for multi-threaded computation
 */
unsigned Generic2DAngularModel::initializeExecData(unsigned type, Image& image, void* params)
{
  unsigned nx      = image.xAxis().getNpix();
  unsigned ny      = image.yAxis().getNpix();

  double dxRad     = image.xAxis().getAngularResolution().radians();
  double dyRad     = image.yAxis().getAngularResolution().radians();
		
  double xSense    = image.xAxis().getSense();
  double ySense    = image.yAxis().getSense();

  double cRotAng   = cos(rotationAngle_.radians());
  double sRotAng   = sin(rotationAngle_.radians());
	
  double prefactor = getEnvelopePrefactor(type, params);
	
  //------------------------------------------------------------ 
  // Get the separation (in the flat-sky approximation) between the
  // image center and this model's center (if both have absolute
  // positions specified, otherwise this function returns 0)
  //------------------------------------------------------------

  getAbsoluteSeparation(image);

  //------------------------------------------------------------
  // And compute a complete offset for this model component
  //------------------------------------------------------------
    
  double xOffRad  = xOffset_.radians() - xSep_.radians();
  double yOffRad  = yOffset_.radians() - ySep_.radians();

  unsigned nThreadTotal  = pool_->nThread();
  unsigned nThreadExec   = nThreadTotal > ny ? ny : nThreadTotal;
  unsigned nRowPerThread = ny / nThreadExec;
  unsigned iStart, iStop;
  
  for(unsigned iThread=0; iThread < nThreadExec; iThread++) {

    iStart = iThread * nRowPerThread;
    iStop  = iStart + nRowPerThread;

    if(iStop > ny)
      iStop = ny;

    ExecData* ed = execData_[iThread];

    ed->image_    = &image;
    ed->iSegment_ = iThread;
    ed->nSegment_ = nThreadExec;
    ed->iYStart_  = iStart;
    ed->iYStop_   = (iThread == nThreadExec-1 && iStop < ny) ? ny : iStop; // Last thread needs to finish

    ed->nx_       = nx;
    ed->ny_       = ny;

    ed->dxRad_    = dxRad;
    ed->dyRad_    = dyRad;

    ed->xSense_   = xSense;
    ed->ySense_   = ySense;
		    
    ed->cRotAng_  = cRotAng;
    ed->sRotAng_  = sRotAng;
		    
    ed->xOffRad_  = xOffRad;
    ed->yOffRad_  = yOffRad;

    ed->prefactor_= prefactor;
    ed->type_     = type;
    ed->params_   = params;

    initializeEvalData(ed->evalData_);
  }

  return nThreadExec;
}

/**.......................................................................
 * Single-threaded method to fill an image with this model
 */
void Generic2DAngularModel::fillImageSingleThread(unsigned type, Image& image, void* params)
{
#ifdef TIMER_TEST
  fit1.start();
#endif
  unsigned nx      = image.xAxis().getNpix();
  unsigned ny      = image.yAxis().getNpix();

  double raRefPix  = image.raRefPix_;
  double decRefPix = image.decRefPix_;

  double dxRad     = image.xAxis().getAngularResolution().radians();
  double dyRad     = image.yAxis().getAngularResolution().radians();

  int    xSense    = image.xAxis().getSense();
  int    ySense    = image.yAxis().getSense();

  double cRotAng   = cos(rotationAngle_.radians());
  double sRotAng   = sin(rotationAngle_.radians());

  double prefactor = getEnvelopePrefactor(type, params);

  //------------------------------------------------------------
  // Get the separation (in the flat-sky approximation) between the
  // image center and this model's center.
  //------------------------------------------------------------

  getAbsoluteSeparation(image);

  //------------------------------------------------------------
  // And compute a total offset for the center of this model
  // component
  //------------------------------------------------------------

  double xOffRad  = xOffset_.radians() - xSep_.radians();
  double yOffRad  = yOffset_.radians() - ySep_.radians();

  //  COUT("xOffset_ = " << xOffset_ << " yOffset_ = " << yOffset_ << " xOffRad = " << xOffRad << " yOffRad = " << yOffRad << " xSep = " << xSep_ << " ySep_ = " << ySep_);

  double x, y, xp, yp, xpp, ypp, arg;

  //------------------------------------------------------------
  // Now iterate over all pixels in this image
  //------------------------------------------------------------

#ifdef TIMER_TEST
  fit1.stop();
  ft1 += fit1.deltaInSeconds();
  fit2.start();
#endif

  for(int iy=0; iy < ny; iy++) {
    for(int ix=0; ix < nx; ix++) {

#ifdef TIMER_TEST
      fit3.start();
#endif

      // Get the coordinate of this pixel relative to the center of
      // the model image.

      x = (double)(ix) - raRefPix;
      y = (double)(iy) - decRefPix;

      x *= dxRad * xSense;
      y *= dyRad * ySense;

      // Compute its coordinate in the untranslated frame

      xp = x - xOffRad;
      yp = y - yOffRad;

      // Compute its coordinate in the unrotated frame

      xpp =   xp * cRotAng - yp * sRotAng;
      ypp =   xp * sRotAng + yp * cRotAng;

      //      COUT(&image << " rarefpix = " << raRefPix << " decrefpix = " << decRefPix << " nx = " << nx << " ny = " << ny << " ix = " << ix << " iy = " << iy << " x = " << x << " xpp = " << xpp << " y = " << y << " ypp = " << ypp);

#ifdef TIMER_TEST
      fit3.stop();
      ft3 += fit3.deltaInSeconds();
      fit4.start();
#endif

      // Now fill the image with the value of the model at this
      // location

      image.data_[iy * nx + ix] = prefactor * envelope(type, xpp, ypp);

#ifdef TIMER_TEST
  fit4.stop();
  ft4 += fit4.deltaInSeconds();
#endif
    }
  }

#ifdef TIMER_TEST
  fit2.stop();
  ft2 += fit2.deltaInSeconds();
#endif

#if 0
  image.hasData_ = true;
  PgUtil::setInteractive(true);
  image.display();
#endif
}

/**.......................................................................
 * Multi-threaded method used to fill part of an image with this model
 */
void Generic2DAngularModel::fillImageMultiThread(ExecData* ed)
{
  //------------------------------------------------------------
  // Now iterate over all pixels for this thread
  //------------------------------------------------------------

  for(ed->iy_ = ed->iYStart_; ed->iy_ < ed->iYStop_; ed->iy_++) {
    for(ed->ix_ = 0; ed->ix_ < ed->nx_; ed->ix_++) {

      // Get the coordinate of the center of this pixel

      ed->x_ = (double)(ed->ix_) - (double)(ed->nx_)/2;
      ed->y_ = (double)(ed->iy_) - (double)(ed->ny_)/2;

      ed->x_ *= ed->dxRad_ * ed->xSense_;
      ed->y_ *= ed->dyRad_ * ed->ySense_;

      // Compute its coordinate in the untranslated frame

      ed->xp_ = ed->x_ - ed->xOffRad_;
      ed->yp_ = ed->y_ - ed->yOffRad_;

      // Compute its coordinate in the unrotated frame

      ed->xpp_ = ed->xp_ * ed->cRotAng_ - ed->yp_ * ed->sRotAng_;
      ed->ypp_ = ed->xp_ * ed->sRotAng_ + ed->yp_ * ed->cRotAng_;

      ed->image_->data_[ed->iy_ * ed->nx_ + ed->ix_] = ed->prefactor_ * envelope(ed->evalData_, ed->type_, ed->xpp_, ed->ypp_);
    }
  }
}

/**.......................................................................
 * Method called by worker threads to execute their portion of the
 * calculation
 */
EXECUTE_FN(Generic2DAngularModel::execFillImageMultiThread)
{
  ExecData* ed = (ExecData*) args;
  Generic2DAngularModel* model = ed->model_;

  model->fillImageMultiThread(ed);
  model->synchronizer_.registerDone(ed->iSegment_, ed->nSegment_);
}

//=======================================================================
// Methods to fill an array
//=======================================================================

/**.......................................................................
 * Fill an array with this model
 */
void Generic2DAngularModel::fillArray(unsigned type, Angle& axisUnits, 
				      std::valarray<double>& x, std::valarray<double>& y, std::valarray<double>& d, 
				      void* params)
{
  //------------------------------------------------------------
  // If no thread pool exists, just call for single thread
  //------------------------------------------------------------

  if(!pool_ || execData_.size() == 0) {

    fillArraySingleThread(type, axisUnits, x, y, d, params);
    
    //------------------------------------------------------------
    // Else split the image into segments
    //------------------------------------------------------------

  } else {

    unsigned nThreadExec = initializeExecData(type, axisUnits, x, y, d, params);

    synchronizer_.reset(nThreadExec);

    for(unsigned iThread=0; iThread < nThreadExec; iThread++) {

      ExecData* ed = execData_[iThread];

      //------------------------------------------------------------
      // Now execute it
      //------------------------------------------------------------

      synchronizer_.registerPending(iThread);
      pool_->execute(&execFillArrayMultiThread, ed);
    }

    synchronizer_.wait();
  }
}

/**.......................................................................
 * Initialize all data needed for multi-threaded computation
 */
unsigned Generic2DAngularModel::initializeExecData(unsigned type, Angle& axisUnits, 
						   std::valarray<double>& x, std::valarray<double>& y, std::valarray<double>& d, 
						   void* params)
{
  double cRotAng   = cos(rotationAngle_.radians());
  double sRotAng   = sin(rotationAngle_.radians());
	
  double prefactor = getEnvelopePrefactor(type, params);
	
  //------------------------------------------------------------
  // And compute a complete offset for this model component
  //------------------------------------------------------------
    
  double xOffRad  = xOffset_.radians();
  double yOffRad  = yOffset_.radians();

  unsigned npt = x.size();
  unsigned nThreadTotal  = pool_->nThread();
  unsigned nThreadExec   = nThreadTotal > npt ? npt : nThreadTotal;
  unsigned nRowPerThread = npt / nThreadExec;
  unsigned iStart, iStop;
  
  axisUnits.setVal(1.0, axisUnits.units());
  double angleConv = axisUnits.radians();

  for(unsigned iThread=0; iThread < nThreadExec; iThread++) {

    iStart = iThread * nRowPerThread;
    iStop  = iStart + nRowPerThread;

    if(iStop > npt)
      iStop = npt;

    ExecData* ed = execData_[iThread];

    ed->xArr_      = &x;
    ed->yArr_      = &y;
    ed->dArr_      = &d;
    ed->angleConv_ = angleConv;

    ed->iSegment_  = iThread;
    ed->nSegment_  = nThreadExec;
    ed->iYStart_   = iStart;
    ed->iYStop_    = (iThread == nThreadExec-1 && iStop < npt) ? npt : iStop; // Last thread needs to finish

    ed->cRotAng_  = cRotAng;
    ed->sRotAng_  = sRotAng;
		    
    ed->xOffRad_  = xOffRad;
    ed->yOffRad_  = yOffRad;

    ed->prefactor_= prefactor;
    ed->type_     = type;
    ed->params_   = params;

    initializeEvalData(ed->evalData_);
  }

  return nThreadExec;
}

/**.......................................................................
 * Single-threaded method to fill an image with this model
 */
void Generic2DAngularModel::fillArraySingleThread(unsigned type, Angle& axisUnits, 
						  std::valarray<double>& x, std::valarray<double>& y, std::valarray<double>& d, 
						  void* params)
{
  axisUnits.setVal(1.0, axisUnits.units());
  double angleConv = axisUnits.radians();

  double prefactor = getEnvelopePrefactor(type, params);

  //------------------------------------------------------------
  // And compute a total offset for the center of this model
  // component
  //------------------------------------------------------------

  double xOffRad  = xOffset_.radians();
  double yOffRad  = yOffset_.radians();

  double cRotAng   = cos(rotationAngle_.radians());
  double sRotAng   = sin(rotationAngle_.radians());

  double xp, yp, xpp, ypp, arg;

  //------------------------------------------------------------
  // Now iterate over all elements in the  array
  //------------------------------------------------------------

  unsigned n = x.size();
  for(int i=0; i < n; i++) {

    // Compute its coordinate in the untranslated frame

    xp = x[i]*angleConv - xOffRad;
    yp = y[i]*angleConv - yOffRad;

    // Compute its coordinate in the unrotated frame

    xpp =   xp * cRotAng - yp * sRotAng;
    ypp =   xp * sRotAng + yp * cRotAng;

    // Now fill the array with the value of the model at this
    // location

    d[i] = prefactor * envelope(type, xpp, ypp);
  }
}

/**.......................................................................
 * Multi-threaded method used to fill part of an array with this model
 */
 void Generic2DAngularModel::fillArrayMultiThread(ExecData* ed)
{
  //------------------------------------------------------------
  // Now iterate over all pixels for this thread
  //------------------------------------------------------------

  for(ed->iy_ = ed->iYStart_; ed->iy_ < ed->iYStop_; ed->iy_++) {

    // Compute its coordinate in the untranslated frame

    ed->xp_ = (*ed->xArr_)[ed->iy_]*ed->angleConv_ - ed->xOffRad_;
    ed->yp_ = (*ed->yArr_)[ed->iy_]*ed->angleConv_ - ed->yOffRad_;
    
    // Compute its coordinate in the unrotated frame
    
    ed->xpp_ = ed->xp_ * ed->cRotAng_ - ed->yp_ * ed->sRotAng_;
    ed->ypp_ = ed->xp_ * ed->sRotAng_ + ed->yp_ * ed->cRotAng_;

    (*ed->dArr_)[ed->iy_] = ed->prefactor_ * envelope(ed->evalData_, ed->type_, ed->xpp_, ed->ypp_);
  }
}

/**.......................................................................
 * Method called by worker threads to execute their portion of the
 * calculation
 */
EXECUTE_FN(Generic2DAngularModel::execFillArrayMultiThread)
{
  ExecData* ed = (ExecData*) args;
  Generic2DAngularModel* model = ed->model_;

  model->fillArrayMultiThread(ed);
  model->synchronizer_.registerDone(ed->iSegment_, ed->nSegment_);
}

/**.......................................................................
 *  Base-class method for filling an image with this model
 */
void Generic2DAngularModel::fillUvData(unsigned type, gcp::util::UvDataGridder& gridder, void* params)
{
  ThrowColorError("No inherited fillUvData() method has been defined by this model", "red");
}

/**.......................................................................
 * Method to write fake 2D data from this model to a file
 */
void Generic2DAngularModel::generateFake2DData(std::string fileName, Image& image, double sigma, std::string units)
{
  fillImage(dataSetType_, image, 0);

  for(unsigned i=0; i < image.data_.size(); i++)
    image.data_[i] += Sampler::generateGaussianSample(sigma);

  image.writeToFitsFile(fileName);
}

/**.......................................................................
 * Method to set generic normalization with no units (throws if not
 * allowed for specific inherited models)
 */
void Generic2DAngularModel::setNormalization(double norm)
{
  if(normalization_.unitlessAllowed()) {
    normalization_.value() = norm;
  } else {
    ThrowColorError("You cannot specify the normalization without units: " << normalization_.getUnitsString(), "red");
  }
}

/**.......................................................................
 * Overload base-class setThreadPool method to define what happens
 * when running in a multi-threaded context
 */
void Generic2DAngularModel::setThreadPool(ThreadPool* pool)
{
  // Call the base-class method

  Model::setThreadPool(pool);

  // ExecData::ExecData() can throw if inheritor hasn't define an
  // allocateEvalData() method.  This just means that the inheritor
  // doesn't support multi-threaded evaluation.  Quietly catch this
  // and the size of execData_ will be used later to determine if
  // multi-threaded evalulation can be used.

  if(pool) {
    try {
      
      // Now initialize per-thread resources once
      
      for(unsigned iThread=0; iThread < pool_->nThread(); iThread++)
	execData_.push_back(new ExecData(this));
      
    } catch(...) {
      
    }
  }
}

/**.......................................................................
 * Base-class method to allocate data needed for multi-threaded running
 */
void* Generic2DAngularModel::allocateEvalData()
{
  ThrowError("This inheritor has not defined a method to allocate evaluation data");
  return 0;
}

/**.......................................................................
 * Base-class method to initialize data needed for multi-threaded running
 */
void Generic2DAngularModel::initializeEvalData(void* evalData)
{
  ThrowError("This inheritor has not defined a method to initialize evaluation data");
}

/**.......................................................................
 * Calculate absolute separation between the center of this image and
 * the center of this model.  If either is unspecified, assume they
 * are coincident.
 */
void Generic2DAngularModel::getAbsoluteSeparation(Image& image)
{
  if(!positionChecked_) {
    hasAbsolutePosition_ = getVar("ra")->hasValue_ || getVar("dec")->hasValue_;
    positionChecked_ = true;
  }

  if(image.hasAbsolutePosition_ && hasAbsolutePosition_) {
    gcp::util::Astrometry::flatSkyApproximationSeparations(xSep_, ySep_, image.getRa(), image.getDec(), ra_, dec_);
  } else {
    xSep_.setRadians(0.0);
    ySep_.setRadians(0.0);
  }
}

/**.......................................................................
 * Calculate absolute separation between the center of this DFT and
 * the center of this model.  If either is unspecified, assume they
 * are coincident.
 */
void Generic2DAngularModel::getAbsoluteSeparation(Dft2d& dft)
{
  if(!positionChecked_) {
    hasAbsolutePosition_ = getVar("ra")->hasValue_ || getVar("dec")->hasValue_;
    positionChecked_ = true;

    COUT("absolute position : " << hasAbsolutePosition_);
  }

  if(dft.hasAbsolutePosition_ && hasAbsolutePosition_) {
    gcp::util::Astrometry::flatSkyApproximationSeparations(xSep_, ySep_, dft.ra_, dft.dec_, ra_, dec_);
  } else {
    xSep_.setRadians(0.0);
    ySep_.setRadians(0.0);
  }
}

/**.......................................................................
 * Base-class version just uses the normalization as specified
 */
double Generic2DAngularModel::getEnvelopePrefactor(unsigned type, void* params)
{
  if(type & DataSetType::DATASET_RADIO) {

    radioNormalization_.throwIfNotSpecified();
    units_ = Unit::stringToUnits(radioNormalization_.units());
    return getRadioEnvelopePrefactor((Frequency*)params);

  } else if(type & DataSetType::DATASET_XRAY_IMAGE) {

    xrayNormalization_.throwIfNotSpecified();
    units_ = Unit::stringToUnits(xrayNormalization_.units());
    return getXrayImageEnvelopePrefactor((Frequency*)params);

  } else {

    normalization_.throwIfNotSpecified();
    units_ = Unit::stringToUnits(normalization_.units());
    return normalization_.value();

  }

}

/**.......................................................................
 * Return the prefactor for multiplying an Xray image model
 */
double Generic2DAngularModel::getXrayImageEnvelopePrefactor(Frequency* freq)
{
  std::string units = xrayNormalization_.units();

  if(units == "MJy/sr") {
    return xrayNormalization_.value();
  } else if(units == "counts") {
    return xrayNormalization_.value();
  } else {
    ThrowSimpleColorError("Unrecognized units: " << units << " for normalization of model " << name_, "red");
    return 0.0;
  }
}

/**.......................................................................
 * Return the prefactor for multiplying a radio model with spectral
 * shape
 */
double Generic2DAngularModel::getRadioEnvelopePrefactor(Frequency* freq)
{
  switch (spectralType_.type_) {
  case SpectralType::SPEC_ALPHA:
    return getSpectralIndexEnvelopePrefactor(freq);
    break;
  case SpectralType::SPEC_SZ:
    return getSzEnvelopePrefactor(freq);
    break;
  case SpectralType::SPEC_ITOH:
    return getItohEnvelopePrefactor(freq);
    break;
  default:
    return radioNormalization_.value();
    break;
  }
}

/**.......................................................................
 * Return the frequency-dependent prefactor for a spectral-index
 * spectrum
 */
double Generic2DAngularModel::getSpectralIndexEnvelopePrefactor(Frequency* freq)
{
  return radioNormalization_.value() * pow(freq->GHz()/normalizationFrequency_.GHz(), spectralIndex_.value());
}

/**.......................................................................
 * Return the frequency-dependent prefactor for an SZ spectrum.  We
 * don't change units here, but just return a prefactor to scale from
 * whatever units at the current frequency to the same units at the
 * normalization frequency.
 */
double Generic2DAngularModel::getSzEnvelopePrefactor(Frequency* freq)
{
  std::string units = radioNormalization_.units();

  //------------------------------------------------------------
  // If units are already in compton Y, there is no frequency
  // dependence
  //------------------------------------------------------------

  if(units == "comptony") {
    return radioNormalization_.value();

  //------------------------------------------------------------
  // If units are already in pressure, there is also no frequency
  // dependence (simple scaling from Y to pressure)
  //------------------------------------------------------------

  } else if(units == "keV/cm^3" || units == "erg/cm^3") {

    return pressureToComptonY() * radioNormalization_.value();

    //------------------------------------------------------------
    // Else return the ratio of the conversion from Y to T at the
    // current frequency, to the conversion at the normalization
    // frequency
    //------------------------------------------------------------

  } else if(units == "muK" || units == "mK" || units == "K") {

    double normConv;
    szCalculator_.comptonYToDeltaT(*freq, normTemperatureConv_);
    normConv = normTemperatureConv_.K();
    szCalculator_.comptonYToDeltaT(normalizationFrequency_, normTemperatureConv_);
    normConv /= normTemperatureConv_.K();
    return radioNormalization_.value() * normConv;

    //------------------------------------------------------------
    // Else return the ratio of the conversion from Y to I at the
    // current frequency, to the conversion at the normalization
    // frequency
    //------------------------------------------------------------

  } else if(units == "Jy") {

    double normConv;
    szCalculator_.comptonYToDeltaI(*freq, normTemperatureConv_, normIntensityConv_);
    normConv = normIntensityConv_.JyPerSr();
    szCalculator_.comptonYToDeltaI(normalizationFrequency_, normTemperatureConv_, normIntensityConv_);
    normConv /= normIntensityConv_.JyPerSr();
    return radioNormalization_.value() * normConv;

  } else {
    ThrowSimpleColorError("Unrecognized units: " << units << " for normalization of model " << name_, "red");
    return 0.0;
  }
}

/**.......................................................................
 * Return the frequency-dependent prefactor for the Itoh relativistic
 * SZ spectrum.  We don't change units here, but just return a
 * prefactor to scale from whatever units at the current frequency to
 * the same units at the normalization frequency.
 */
double Generic2DAngularModel::getItohEnvelopePrefactor(Frequency* freq)
{
  std::string units = radioNormalization_.units();

  //------------------------------------------------------------
  // If units are already in compton Y, there is no frequency
  // dependence
  //------------------------------------------------------------

  if(units == "comptony") {
    return radioNormalization_.value();

  //------------------------------------------------------------
  // If units are already in pressure, there is also no frequency
  // dependence (simple scaling from Y to pressure)
  //------------------------------------------------------------

  } else if(units == "keV/cm^3" || units == "erg/cm^3") {

    return pressureToComptonY() * radioNormalization_.value();

    //------------------------------------------------------------
    // Else return the ratio of the conversion from Y to T at the
    // current frequency, to the conversion at the normalization
    // frequency
    //------------------------------------------------------------

  } else if(units == "muK" || units == "mK" || units == "K") {

    double normConv;
    szCalculator_.comptonYToDeltaTItoh(electronTemperature_, *freq, normTemperatureConv_);
    normConv = normTemperatureConv_.K();
    szCalculator_.comptonYToDeltaTItoh(electronTemperature_, normalizationFrequency_, normTemperatureConv_);
    normConv /= normTemperatureConv_.K();
    return radioNormalization_.value() * normConv;

    //------------------------------------------------------------
    // Else return the ratio of the conversion from Y to I at the
    // current frequency, to the conversion at the normalization
    // frequency
    //------------------------------------------------------------

  } else if(units == "Jy") {

    double normConv;
    szCalculator_.comptonYToDeltaIItoh(electronTemperature_, *freq, normTemperatureConv_, normIntensityConv_);
    normConv = normIntensityConv_.JyPerSr();
    szCalculator_.comptonYToDeltaIItoh(electronTemperature_, normalizationFrequency_, normTemperatureConv_, normIntensityConv_);
    normConv /= normIntensityConv_.JyPerSr();
    return radioNormalization_.value() * normConv;

  } else {
    ThrowSimpleColorError("Unrecognized units: " << units << " for normalization of model " << name_, "red");
    return 0;
  }
}

/**.......................................................................
 * Check the setup of this model for logical errors
 */
void Generic2DAngularModel::checkSetup()
{
  Model::checkSetup();

  checkSpectralSetup();
  checkNormalizationSetup();
}

/**.......................................................................
 * By default, we check that at least one of the standard
 * normalization variates was specified
 */
void Generic2DAngularModel::checkNormalizationSetup()
{
  if(!loadedFromOutputFile_) {
    if(!(getVar("normalization")->wasSpecified_ ||
	 getVar("Sradio")->wasSpecified_ ||
	 getVar("Sxray")->wasSpecified_)) {
      
      ThrowColorError(std::endl << "At least one of the following normalizations must be specified: " 
		      << std::endl << std::endl
		      << "  " << name_ << ".normalization" << std::endl
		      << "  " << name_ << ".Sradio"        << std::endl
		      << "  " << name_ << ".Sxray"         << std::endl, "red");
    }
  }

  //------------------------------------------------------------
  // Also check units for unit'd quantities
  //------------------------------------------------------------

  if(getVar("Sradio")->wasSpecified_) {
    Frequency freq;
    freq.setGHz(30.0);

    try {
      Angle size(Angle::Degrees(), 1.0);
      SolidAngle beam;
      Image image(256, size);
      image.setUnits(Unit::stringToUnits(radioNormalization_.units()));
      double fac = image.nativeToJy(freq, beam);
    } catch(...) {
      ThrowColorError(std::endl << "You must specify one of the following units for " << name_ << ".Sradio: " 
		      << std::endl << std::endl
		      << "      muK   MicroKelvin"      << std::endl
		      << "       mK   MilliKelvin"      << std::endl
		      << "        K   Kelvin"           << std::endl
		      << " comptony   Compton Y"        << std::endl
		      << "    Jy/sr   Jansky/steradian" << std::endl
		      << "   MJy/sr   MegaJansky/steradian" << std::endl, "red");
    }
  }
}

/**.......................................................................
 * Check spectral setup for logical errors
 */
void Generic2DAngularModel::checkSpectralSetup()
{
  if(getVar("normalizationFrequency")->isVariable()) {
    ThrowColorError("'" << name_ << ".normalizationFrequency' cannot be variable", "red");
  }

  //------------------------------------------------------------
  // If no spectrum has been specified, see if we can determine what
  // it is intended to be
  //------------------------------------------------------------

  if(!getVar("spectralType")->hasValue_) {

    // If someone has specified a spectral index, then presumably they
    // intend a spectral index

    if(getVar("spectralIndex")->hasValue_) {
      spectralType_.type_ = SpectralType::SPEC_ALPHA;

      // if someone has specified a normalization in units of
      // comptonY, then presumably they intend an SZ spctrum

    } else if(getVar("normalization")->hasValue_ && getVar("normalization")->units() == "comptonY") {
      spectralType_.type_ = SpectralType::SPEC_SZ;
    }
  }

  //------------------------------------------------------------
  // Now that we have our best-guess at the type, check required
  // parameters for this type
  //------------------------------------------------------------

  switch (spectralType_.type_) {
  case SpectralType::SPEC_NONE: // No spectral type specified.  Make
		                // sure the user hasn't specified a
		                // normalization frequency in this
		                // case

    if(getVar("spectralIndex")->hasValue_) {
      ThrowColorError(std::endl << "You have specified a spectral index for model '" << name_ 
		      << "', but no spectral type.", "red");
    }

    if(getVar("normalizationFrequency")->hasValue_) {
      ThrowColorError(std::endl << "You have specified a normalization frequency for model '" << name_ 
		      << "', but no spectral type.", "red");
    }

    break;

  case SpectralType::SPEC_ALPHA: // Spectral index specified.  Make
          	                 // sure the user has specified a
          	                 // normalization frequency in this
          	                 // case and make sure it isn't
          	                 // variable -- that doesn't make
          	                 // sense.

    if(!getVar("spectralIndex")->hasValue_ && !getVar("spectralIndex")->isVariable()) {
      ThrowColorError(std::endl << "You have specified 'alpha' for the spectrum of model '" 
		      << name_ << "', but no spectral index", "red");
    }

    if(!getVar("normalizationFrequency")->hasValue_) {
      ThrowColorError(std::endl << "You have specified a spectral type for model '" << name_ 
		      << "', but no normalization frequency", "red");
    }

    break;

  case SpectralType::SPEC_SZ:    // SZ spectrum specified.  Make sure
                                 // the user has specified a
                                 // normalization frequency in this
                                 // case and make sure it isn't
                                 // variable -- that doesn't make
                                 // sense.

    if(getVar("spectralIndex")->hasValue_) {
      ThrowColorError(std::endl << "You have specified 'sz' for the spectrum of model '" 
		      << name_ << "', and also a spectral index", "red");
    }

    if(!getVar("normalizationFrequency")->hasValue_) {
      ThrowColorError(std::endl << "You have specified a spectral type for model '" << name_ 
		      << "', but no normalization frequency", "red");
    }

    break;

  case SpectralType::SPEC_ITOH:  // Itoh approximation to SZ spectrum
                                 // specified.  Make sure the user has
                                 // specified a normalization
                                 // frequency in this case and make
                                 // sure it isn't variable -- that
                                 // doesn't make sense.

    if(getVar("spectralIndex")->hasValue_) {
      ThrowColorError(std::endl << "You have specified 'sz' for the spectrum of model '" 
		      << name_ << "', and also a spectral index", "red");
    }

    if(!getVar("normalizationFrequency")->hasValue_) {
      ThrowColorError(std::endl << "You have specified a spectral type for model '" << name_ 
		      << "', but no normalization frequency", "red");
    }

    if(getVar("spectralIndex")->hasValue_) {
      ThrowColorError(std::endl << "You have specified 'sz' for the spectrum of model '" 
		      << name_ << "', and also a spectral index", "red");
    }

    if(!getParameter("electronTemperature", false)->data_.hasValue()) {
      ThrowColorError(std::endl << "You have specified the relativistic 'itoh' expansion for the spectrum of model '" 
		      << name_ << "', but no electron temperature (use '" << name_ << ".electronTemperature')", "red");
    }

    electronTemperature_.setVal(getDoubleVal("electronTemperature"), getParameter("electronTemperature", true)->units_);

    break;

  default:
    break;
  }

}

/**.......................................................................
 * Overloaded function just calls down to the base-class.  We define
 * this so that inheritors from us can overload this function too.
 */
void Generic2DAngularModel::sample()
{
  Model::sample();
}

/**.......................................................................
 * 
 */
bool Generic2DAngularModel::normalizationImpliesSz()
{
  // Was a normalization specified as comptonY?

  bool isComptonY = getVar("normalization")->hasValue_ && getVar("normalization")->units() == "comptonY";

  // Was a pressure specified?

  bool isPressure = getVar("pressure")->hasValue_;

  return isComptonY || isPressure;
}

void Generic2DAngularModel::debugPrint()
{
#ifdef TIMER_TEST
  COUTCOLOR("ft1 = " << ft1               << "s", "yellow");
  COUTCOLOR("ft2 = " << ft2               << "s", "yellow");
  COUTCOLOR("ft3 = " << ft3               << "s", "yellow");
  COUTCOLOR("ft4 = " << ft4               << "s", "yellow");
#endif
}

/**.......................................................................
 * Integrate this model from the center out to a radius specified in rad
 */
SolidAngle Generic2DAngularModel::solidAngleIntegral(unsigned type, Angle radius)
{
  SolidAngle angle;
  ThrowError("Inheritor has not defined an integrate() function");
  return angle;
}

void Generic2DAngularModel::solidAngleIntegral(unsigned type, Angle& radius, SolidAngle& sa)
{
  ThrowError("Inheritor has not defined an integrate() function");
}

/**.......................................................................
 * Calculate integrated Y out to the specified (physical) radius
 */
gcp::util::Area Generic2DAngularModel::integratedY(gcp::util::Cosmology& cosmo, gcp::util::Length& radius)
{
  //------------------------------------------------------------
  // Compute the angle on the sky that corresponds to the requested
  // radius
  //------------------------------------------------------------

  Length dA = cosmo.angularDiameterDistance();
  Angle theta(Angle::Radians(), radius/dA);

  //------------------------------------------------------------
  // Get the scale factor to go from native radio noarmalization units
  // to compton-Y
  //------------------------------------------------------------

  double scaleFactor = getComptonYScaleFactor(normalizationFrequency_);

  //------------------------------------------------------------
  // Now compute the integral out to the specified radius
  //------------------------------------------------------------

  Area ret;
  ret.setSquaredMpc(dA.Mpc() * dA.Mpc() * scaleFactor * radioNormalization_.value() * solidAngleIntegral(DataSetType::DATASET_RADIO, theta).sr());

  return ret;
}

/**.......................................................................
 * Calculate integrated Y out to the specified (physical) radius
 */
void Generic2DAngularModel::integratedY(double scaleFactor, Length& dA, Angle& theta, Area& ret)
{
  ret.setSquaredMpc(dA.Mpc() * dA.Mpc() * scaleFactor * radioNormalization_.value() * solidAngleIntegral(DataSetType::DATASET_RADIO, theta).sr());
}

/**.......................................................................
 * Return the (possibly) frequency-dependent scale factor to convert
 * from the current radio normalization units to Comptony-Y.
 */
double Generic2DAngularModel::getComptonYScaleFactor()
{
  return getComptonYScaleFactor(normalizationFrequency_);
}

/**.......................................................................
 * Get the scale factor to convert from whatever units the
 * normalization has been specified in, to Compton-y
 */
double Generic2DAngularModel::getComptonYScaleFactor(Frequency& freq)
{
  std::string units = radioNormalization_.units();

  //------------------------------------------------------------
  // If units are already in compton Y, there is no scale factor
  //------------------------------------------------------------

  if(units == "comptony") {
    return 1.0;

    //------------------------------------------------------------
    // Else return the ratio of the conversion from Y to T at the
    // current frequency, to the conversion at the normalization
    // frequency
    //------------------------------------------------------------

  } else if(units == "muK" || units == "mK" || units == "K") {

    szCalculator_.comptonYToDeltaT(freq, normTemperatureConv_);

    if(units == "muK")
      return 1.0/normTemperatureConv_.microK();
    else if(units == "mK")
      return 1.0/normTemperatureConv_.milliK();
    else
      return 1.0/normTemperatureConv_.K();

    //------------------------------------------------------------
    // Else return the ratio of the conversion from Y to I at the
    // current frequency, to the conversion at the normalization
    // frequency
    //------------------------------------------------------------

  } else if(units == "Jy" || units == "mJy") {

    szCalculator_.comptonYToDeltaI(freq, normTemperatureConv_, normIntensityConv_);

    if(units == "mJy")
      return 1.0/normIntensityConv_.mJyPerSr();
    else if(units == "Jy")
      return 1.0/normIntensityConv_.JyPerSr();

  }

  ThrowSimpleColorError("Unrecognized units: " << units << " for normalization of model " << name_, "red");
  return 0.0;
}

PgModel Generic2DAngularModel::pgModel()
{
  PgModel mod;

  mod.xMid_  = xOffset_.degrees();
  mod.yMid_  = yOffset_.degrees();
  mod.angle_ = rotationAngle_;
  mod.rot_   = rotationAngle_;
  mod.type_  = PgModel::TYPE_DELTA;

  return mod;
}

/**.......................................................................
 * Check for absolute position specification.  This is called externally
 * so that external position information from an initialization
 * script can override internal position information from a data file.
 */
void Generic2DAngularModel::checkPosition()
{
  if(!positionChecked_) {
    hasAbsolutePosition_ = getVar("ra")->hasValue_ || getVar("dec")->hasValue_;
    positionChecked_ = true;
  }
}

