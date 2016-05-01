#include "gcp/fftutil/Image.h"
#include "gcp/fftutil/UvDataGridder.h"

#include "gcp/models/PtSrcModel.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::models;

//#define TIMER_TEST

#ifdef TIMER_TEST
#include "gcp/util/Timer.h"
Timer fit1, fit2, fit3, fit4, fit5;
double ft1=0.0, ft2=0.0, ft3=0.0, ft4=0.0, ft5=0.0;
#endif

std::valarray<double> PtSrcModel::sinLookup_;
double PtSrcModel::lookupResRadians_;
unsigned PtSrcModel::lookupNpt_;

/**.......................................................................
 * Constructor.
 */
PtSrcModel::PtSrcModel()
{
  dataSetType_ |= DataSetType::DATASET_RADIO;
  dataSetType_ |= DataSetType::DATASET_XRAY_IMAGE;
  dataSetType_ |= DataSetType::DATASET_PTSRC;

  // And initialize all components to fixed values

  initializeComponentsToFixed();

  // Initialize the lookup table too

  initializeLookupTable();
}

/**.......................................................................
 * Destructor.
 */
PtSrcModel::~PtSrcModel() {}

/**.......................................................................
 * Set the flux of this source
 */
void PtSrcModel::setFlux(Flux& flux)
{
  setJy(flux.Jy());
}

void PtSrcModel::setJy(double fluxJy)
{
  radioNormalization_.setVal(fluxJy, "Jy");
  radioNormalization_.wasSpecified_ = true;
}

/**.......................................................................
 * Set the frequency of this source
 */
void PtSrcModel::setFrequency(Frequency& freq)
{
  setGHz(freq.GHz());
}

void PtSrcModel::setGHz(double ghz)
{
  normalizationFrequency_.setGHz(ghz);
  normalizationFrequency_.wasSpecified_ = true;
}

/**.......................................................................
 * Set the spectral index of this source
 */
void PtSrcModel::setSpectralIndex(double spectralIndex)
{
  spectralIndex_ = spectralIndex;
}

/**.......................................................................
 * Set an offset for this source
 */
void PtSrcModel::setOffset(Angle& xOff, Angle& yOff)
{
  setXOffset(xOff);
  setYOffset(yOff);
}

void PtSrcModel::fillImage(unsigned type, gcp::util::Image& image, void* params)
{
  Frequency* freq = (Frequency*)params;
  if(type & DataSetType::DATASET_RADIO) {
    fillSzImage(image, *freq);
  } else if(type & DataSetType::DATASET_XRAY_IMAGE) {
    fillXrayImage(image, *freq);
  } else {
    ThrowError("Unsupported image type for this model");
  }
}

/**.......................................................................
 * Construct an image suitable for use with an SZ dataset
 */
void PtSrcModel::fillSzImage(Image& image, Frequency& frequency)
{
  //------------------------------------------------------------
  // Zero the image
  //------------------------------------------------------------

  image.zero();

  //------------------------------------------------------------
  // Get the absolute separation between the image center and this
  // model's center
  //------------------------------------------------------------

  getAbsoluteSeparation(image);

  //------------------------------------------------------------
  // And compute a total offset for the center of this model
  // component
  //------------------------------------------------------------

  int    xSense    = image.xAxis().getSense();
  int    ySense    = image.yAxis().getSense();

  double xOffRad  = xOffset_.radians() - xSep_.radians();
  double yOffRad  = yOffset_.radians() - ySep_.radians();

  // Now calculate the pixel corresponding to the current source
  // location

  unsigned ix = (int)(xOffRad / image.xAxis().getAngularResolution().radians() * xSense) + image.raRefPix_;
  unsigned iy = (int)(yOffRad / image.yAxis().getAngularResolution().radians() * ySense) + image.decRefPix_;

  // Set the pixel corresponding to the current location of this
  // source to the scaled flux

  image.val(ix, iy) = getEnvelopePrefactor(DataSetType::DATASET_RADIO, &frequency);

#if 0
  COUT("xOff = " << xOffset_);
  COUT("Setting pixel " << ix << " " << iy << " to " << image.val(ix, iy) << " "
       << " xsize = " << image.xAxis().getAngularSize()
       << " ysize = " << image.yAxis().getAngularSize());
#endif

  // Set the units of this data

  image.setUnits(Unit::UNITS_JY);
  image.hasData_ = true;
}

/**.......................................................................
 * Construct an image suitable for use with an SZ dataset
 */
void PtSrcModel::fillXrayImage(Image& image, Frequency& frequency)
{
  //------------------------------------------------------------
  // Zero the image
  //------------------------------------------------------------

  image.zero();

  //------------------------------------------------------------
  // Get the absolute separation between the image center and this
  // model's center
  //------------------------------------------------------------

  getAbsoluteSeparation(image);

  //------------------------------------------------------------
  // And compute a total offset for the center of this model
  // component
  //------------------------------------------------------------

  int    xSense    = image.xAxis().getSense();
  int    ySense    = image.yAxis().getSense();

  double xOffRad  = xOffset_.radians() - xSep_.radians();
  double yOffRad  = yOffset_.radians() - ySep_.radians();

  // Now calculate the pixel corresponding to the current source
  // location

  unsigned ix = (int)(xOffRad / image.xAxis().getAngularResolution().radians() * xSense) + image.raRefPix_;
  unsigned iy = (int)(yOffRad / image.yAxis().getAngularResolution().radians() * ySense) + image.decRefPix_;

  // Set the pixel corresponding to the current location of this
  // source to the scaled flux

  image.val(ix, iy) = getEnvelopePrefactor(DataSetType::DATASET_XRAY_IMAGE, &frequency);

#if 0
  COUT("xOff = " << xOffset_);
  COUT("Setting pixel " << ix << " " << iy << " to " << image.val(ix, iy) << " "
       << " xsize = " << image.xAxis().getAngularSize()
       << " ysize = " << image.yAxis().getAngularSize());
#endif

  // Set the units of this data

  image.setUnits(Unit::UNITS_JY);
  image.hasData_ = true;
}

/**.......................................................................
 * Construct a DFT suitable for use with an SZ dataset
 */
void PtSrcModel::fillUvData(unsigned type, UvDataGridder& gridder, void* args)
{
#ifdef TIMER_TEST
  fit1.start();
#endif

  UvParams* params = (UvParams*)args;
  Image* beam      = params->beam_;
  Frequency* freq  = params->freq_;

  double beamVal;
  bool valid = false;

  //------------------------------------------------------------
  // Interpolate the beam to get the approximate point in the beam
  // envelope corresponding to the current point source offset
  //------------------------------------------------------------

  beam->interpolateData(xOffset_, yOffset_, beamVal, valid);

#if 0
  COUT("Inside fillUvData");
  PgUtil::setInteractive(true);
  beam->display();
#endif

  if(!valid) {
    ThrowColorError("Unable to interpolate the beam for frequency: " << (*freq)
		    << ". Image needs to be larger to accomodate offset: " << xOffset_ << ", " << yOffset_, "red");
  }

  //------------------------------------------------------------
  // Calculate the amplitude of the source, modulated by the beam
  //------------------------------------------------------------

  double ampJy = beamVal * getEnvelopePrefactor(type, freq);

  //------------------------------------------------------------
  // Get the current source offset in radians
  //------------------------------------------------------------

  getAbsoluteSeparation(gridder);

  double xRad  = xOffset_.radians() - xSep_.radians();
  double yRad  = yOffset_.radians() - ySep_.radians();

  unsigned dftInd;
  double u, v, arg;
  double sVal, cVal;

#ifdef TIMER_TEST
  fit1.stop();
  ft1 += fit1.deltaInSeconds();
  fit2.start();
#endif

  //------------------------------------------------------------
  // Now iterate (only) over populated indices in the dft, calculating
  // appropriate Fourier-space components
  //------------------------------------------------------------

  for(unsigned i=0; i < gridder.populatedIndices_.size(); i++) {

    dftInd = gridder.populatedIndices_[i];
    u      = gridder.populatedU_[i];
    v      = gridder.populatedV_[i];

    arg = 2*M_PI*(u * xRad + v * yRad);

#ifdef TIMER_TEST
  fit3.start();

  gridder.out_[dftInd][0] =  ampJy;
  gridder.out_[dftInd][1] = -ampJy;

  fit3.stop();
  ft3 += fit3.deltaInSeconds();

  fit4.start();

  PtSrcModel::sincos(arg, sVal, cVal);
  gridder.out_[dftInd][0] =  ampJy * cVal;
  gridder.out_[dftInd][1] = -ampJy * sVal;

  fit4.stop();
  ft4 += fit4.deltaInSeconds();
  fit5.start();
#endif

  gridder.out_[dftInd][0] =  ampJy * ::cos(arg);
  gridder.out_[dftInd][1] = -ampJy * ::sin(arg);

#ifdef TIMER_TEST
  fit5.stop();
  ft5 += fit5.deltaInSeconds();
#endif

  }

  // Mark this object as containing valid data

  gridder.setHasData(true);

  if(type & DataSetType::DATASET_RADIO)
    units_ = Unit::stringToUnits(radioNormalization_.units());
  else if(type & DataSetType::DATASET_XRAY_IMAGE)
    units_ = Unit::stringToUnits(xrayNormalization_.units());
  else
    units_ = Unit::stringToUnits(normalization_.units());

  gridder.setUnits(units_);

#ifdef TIMER_TEST
  fit2.stop();
  ft2 += fit2.deltaInSeconds();
#endif
}

void PtSrcModel::debugPrint()
{
#ifdef TIMER_TEST
  COUTCOLOR("ft1 = " << ft1               << "s", "yellow");
  COUTCOLOR("ft2 = " << ft2               << "s", "yellow");
  COUTCOLOR("ft3 = " << ft3               << "s", "yellow");
  COUTCOLOR("ft4 = " << ft4               << "s", "yellow");
  COUTCOLOR("ft5 = " << ft5               << "s", "yellow");
#endif
}


/**.......................................................................
 * Initialize the sin() lookup table if it hasn't already been
 */
void PtSrcModel::initializeLookupTable()
{
  // Only initialize the lookup table if it hasn't already been initialized

  if(sinLookup_.size() == 0)
    initializeLookupTable(1e-5);
}

/**.......................................................................
 * Initialize the sin() lookup table to support the required precision
 */
void PtSrcModel::initializeLookupTable(double precision)
{
  // Sin function changes most about zero,when the function value is
  // also zero.  We want to know what resolution we require to achieve
  // the requested precision:

  double dxRad = asin(precision);
  unsigned nPt = (unsigned)ceil(2*M_PI / dxRad);

  // And the final resolution at which we will precompute the lookup
  // table is 2*M_PI divided by nPt

  dxRad = 2*M_PI/nPt;

  sinLookup_.resize(nPt);
  for(unsigned i=0; i < nPt; i++) {
    sinLookup_[i] = ::sin(dxRad * i);
  }

  lookupResRadians_ = dxRad;
  lookupNpt_        = nPt;
}

/**.......................................................................
 * Class method for getting sin() from our internal lookup table
 */
double PtSrcModel::sin(double argRad)
{
  int sign = 1;

  // Convert negative arguments to their positive counterparts

  if(argRad < 0.0) {
    argRad = -argRad;
    sign = -1;
  }

  // Get the nearest index in the lookup table

  unsigned ind = (argRad / lookupResRadians_);

  // And truncate for arg > 2*M_PI

  ind %= lookupNpt_;

  return sign * sinLookup_[ind];
}

/**.......................................................................
 * Class method for getting cos() from our internal lookup table
 */
double PtSrcModel::cos(double argRad)
{
  // Convert negative arguments to their positive counterparts

  if(argRad < 0.0)
    argRad = -argRad;

  // Get the nearest index in the sin lookup table.  Add M_PI/2 to get
  // the cos value

  unsigned ind = ((argRad + M_PI/2) / lookupResRadians_);

  // And truncate for arg > 2*M_PI

  ind %= lookupNpt_;

  return sinLookup_[ind];
}

/**.......................................................................
 * Compute both sin and cos in one function call
 */
void PtSrcModel::sincos(double argRad, double& sVal, double& cVal)
{
  int sign = 1;

  // Convert negative arguments to their positive counterparts

  if(argRad < 0.0) {
    argRad = -argRad;
    sign = -1;
  }

  // Get the nearest indices in the sin lookup table.  Add M_PI/2 to get
  // the cos value

  unsigned sinInd = ((argRad) / lookupResRadians_);
  unsigned cosInd = ((argRad + M_PI/2) / lookupResRadians_);

  // Truncate for arg > 2*M_PI

  sinInd %= lookupNpt_;
  cosInd %= lookupNpt_;
  
  sVal = sign * sinLookup_[sinInd];
  cVal = sinLookup_[cosInd];
}
