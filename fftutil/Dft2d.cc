#include "gcp/fftutil/Dft2d.h"
#include "gcp/pgutil/PgUtil.h"

#include <vector>

using namespace std;
using namespace gcp::util;

const double Dft2d::convSigInPixels_  = 0.594525;
const int    Dft2d::convMaskInPixels_ = 2;
Mutex Dft2d::planGuard_;

/**.......................................................................
 * Constructors
 */
Dft2d::Dft2d(bool optimize)
{
  initialize();
  optimize_  = optimize;
}

Dft2d::Dft2d(const Dft2d& dft)
{
  *this = dft;
}

void Dft2d::operator=(const Dft2d& dft)
{
  *this = (Dft2d&)dft;
}

void Dft2d::operator=(Dft2d& dft)
{
  initialize();

  xAxis_ = dft.xAxis_;
  yAxis_ = dft.yAxis_;
}

/**.......................................................................
 * Initialize this object to defaults
 */
void Dft2d::initialize() 
{
  debugPrint_ = false;

  in_          = 0;
  out_         = 0;

  optimize_      = true;
  normalize_     = false;
  isTransformed_ = false;

  nx_          = 0;
  ny_          = 0;
  nIn_         = 0;
  nOut_        = 0;

  nxZeroPad_   = 0;
  nyZeroPad_   = 0;
  nInZeroPad_  = 0;
  nOutZeroPad_ = 0;

  nxZeroPadPrev_ = 0;
  nyZeroPadPrev_ = 0;

  fwdPlanTmp_  = 0;
  invPlanTmp_  = 0;
  precomputedPlan_ = false;

  fwdPlanComputed_  = false;
  invPlanComputed_  = false;

  axes_ = ImageAxis::AXIS_NONE;

  xAxis_.setAxisType(Axis::AXIS_X);
  yAxis_.setAxisType(Axis::AXIS_Y);

  xAxis_.parent_ = this;
  yAxis_.parent_ = this;
}

/**.......................................................................
 * Constructor with size arguments only
 */
Dft2d::Dft2d(int nx, int ny, bool optimize) 
{
  initialize();
  optimize_  = optimize;

  resize(nx, ny);
}

/**.......................................................................
 * Constructor with initialization of input data from an Image
 */
Dft2d::Dft2d(Image& image, bool optimize) 
{
  initialize();
  optimize_  = optimize;

  initialize(image);  
}

/**.......................................................................
 * Destructor.
 */
Dft2d::~Dft2d() 
{
  if(in_) {
    fftw_free(in_);
    in_ = 0;
  }

  if(out_) {
    fftw_free(out_);
    out_ = 0;
  }

  if(fwdPlanTmp_) {
    fftw_destroy_plan(forwardPlan_);
  }

  if(invPlanTmp_) {
    fftw_destroy_plan(inversePlan_);
  }
}

/**.......................................................................
 * True to normalize the transform
 */
void Dft2d::normalize(bool norm)
{
  normalize_ = norm;
}

void Dft2d::zeropad(bool pad, unsigned fac)
{
  if(pad) {
    xAxis_.setZeropadFactor(fac);
    yAxis_.setZeropadFactor(fac);
  } else {
    xAxis_.setZeropadFactor(1);
    yAxis_.setZeropadFactor(1);
  }
}

/**.......................................................................
 * Compute a plan for this fft
 */
void Dft2d::computePlan(unsigned flag) 
{
  planGuard_.lock();

  if(fwdPlanTmp_ == 0 && invPlanTmp_ == 0) {

    //    COUT(pthread_self() << " Computing plans for nxZeroPad_ = " << nxZeroPad_);
    //    COUT(pthread_self() << " Computing plans for nyZeroPad_ = " << nyZeroPad_);
    
    if(fwdPlanComputed_) {
      fftw_destroy_plan(forwardPlan_);
    }

    if(invPlanComputed_) {
      fftw_destroy_plan(inversePlan_);
    }

    forwardPlan_ = fftw_plan_dft_r2c_2d(nxZeroPad_, nyZeroPad_, in_,  out_, flag);
    inversePlan_ = fftw_plan_dft_c2r_2d(nxZeroPad_, nyZeroPad_, out_, in_,  flag);

    fwdPlanComputed_ = true;
    invPlanComputed_ = true;
    
    //    COUT(pthread_self() << " Computed " << forwardPlan_);
  }

  planGuard_.unlock();
}

/**.......................................................................
 * Compute a plan for this fft
 */
void Dft2d::computePlan(fftw_plan& fwdPlan, fftw_plan& invPlan) 
{
  planGuard_.lock();

  unsigned flag = optimize_ ? FFTW_MEASURE : FFTW_ESTIMATE;

  if(fwdPlanTmp_ == 0 && invPlanTmp_ == 0) {

    COUT(pthread_self() << " Pre-computing plans for nxZeroPad_ = " << nxZeroPad_);
    COUT(pthread_self() << " Pre-computing plans for nyZeroPad_ = " << nyZeroPad_);
    
    fwdPlan = fftw_plan_dft_r2c_2d(nxZeroPad_, nyZeroPad_, in_,  out_, flag);
    invPlan = fftw_plan_dft_c2r_2d(nxZeroPad_, nyZeroPad_, out_, in_,  flag);
    
    COUT(pthread_self() << " Computing plans for nyZeroPad_ = " << nyZeroPad_ << " ...done");
  }

  planGuard_.unlock();
}

/**.......................................................................
 * Initialize the input data to this transform from an image
 */
void Dft2d::initialize(Image& image)
{
  xAxis_ = image.xAxis();
  yAxis_ = image.yAxis();

  setInput(image);
}

/**.......................................................................
 * Initialize the input data to this transform from an image.  This
 * method assumes the image and dft already have the same size
 */
void Dft2d::setInput(Image& image)
{
  if(image.hasData_) {

    if(nx_*ny_ != image.data_.size()) {
      ThrowError("image = " << &image << " Cannot set input from image of a different size: nx_ = " << nx_ << " ny_ = " << ny_ << " size = " << image.data_.size());
    }

    unsigned imInd, dftInd;
    for(unsigned ix=0; ix < nx_; ix++) {
      for(unsigned iy=0; iy < ny_; iy++) {
	dftInd = (ix + xOffset_) * nyZeroPad_ + (iy + yOffset_);
	imInd  = iy * nx_ + ix;
	in_[dftInd] = image.data_[imInd];
      }
    }
   
    hasData_ = true;
  }
}

/**.......................................................................
 * Resize for FFT of a different length
 */
void Dft2d::resize(unsigned nx, unsigned ny)
{
  xAxis_.setNpix(nx);
  yAxis_.setNpix(ny);
}

/**.......................................................................
 * Resize for FFT of a different length
 */
void Dft2d::resize()
{
  nx_   = xAxis_.logicalAxis_.getNpix();
  ny_   = yAxis_.logicalAxis_.getNpix();
  nIn_  = nx_ * ny_;
  nOut_ = nx_ * (ny_/2+1);

  nxZeroPad_   = xAxis_.zeropadFactor_ * nx_;
  nyZeroPad_   = yAxis_.zeropadFactor_ * ny_;
  nInZeroPad_  = nxZeroPad_ * nyZeroPad_;
  nOutZeroPad_ = nxZeroPad_ * (nyZeroPad_/2+1);

  //  COUT("Current " << nxZeroPad_ << " and " << nyZeroPad_);
  //  COUT("Previos " << nxZeroPadPrev_ << " and " << nyZeroPadPrev_);

  if(nxZeroPad_ == nxZeroPadPrev_ && nyZeroPad_ == nyZeroPadPrev_)
    return;

  // Store these values for later comparison if we resize

  //  COUT("Storing " << nxZeroPad_ << " and " << nyZeroPad_);

  nxZeroPadPrev_ = nxZeroPad_;
  nyZeroPadPrev_ = nyZeroPad_;

  xOffset_ = (xAxis_.zeropadFactor_ == 1) ? 0 : (xAxis_.zeropadFactor_/2 - 1) * nx_ + nx_/2;
  yOffset_ = (yAxis_.zeropadFactor_ == 1) ? 0 : (yAxis_.zeropadFactor_/2 - 1) * ny_ + ny_/2;

  if(in_) {
    fftw_free(in_);
    in_ = 0;
  }

  if(out_) {
    fftw_free(out_);
    out_ = 0;
  }

  // Allocate arrays

  if((in_ = (double*)fftw_malloc(nInZeroPad_ * sizeof(double)))==0)
    ThrowError("Couldn't allocate input data array");

  if((out_ = (fftw_complex*)fftw_malloc(nOutZeroPad_ * sizeof(fftw_complex)))==0)
    ThrowError("Couldn't allocate output data array");

  for(unsigned i=0; i < nInZeroPad_; i++)
    in_[i] = 0.0;

  for(unsigned i=0; i < nOutZeroPad_; i++) {
    out_[i][0] = 0.0;
    out_[i][1] = 0.0;
  }

  // Now compute a plan for this machine architecture
      
  if(!precomputedPlan_) 
    computePlan(optimize_ ? FFTW_MEASURE : FFTW_ESTIMATE);
}

/**.......................................................................
 * Remove the mean
 */
void Dft2d::removeMean()
{
  double mean=0.0;
  unsigned dftInd;

  // Only compute the mean over the non-zeropadded part of the array

  for(unsigned i=0, iy=0; iy < ny_; iy++) {
    for(unsigned ix=0; ix < nx_; ix++, i++) {
      dftInd = (ix + xOffset_) * nyZeroPad_ + (iy + yOffset_);
      mean += (in_[dftInd] - mean)/(i + 1);
    }
  }

  for(unsigned i=0, iy=0; iy < ny_; iy++) {
    for(unsigned ix=0; ix < nx_; ix++, i++) {
      dftInd = (ix + xOffset_) * nyZeroPad_ + (iy + yOffset_);
      in_[dftInd] -= mean;
    }
  }
}

/**.......................................................................
 * Compute the forward transform
 */
void Dft2d::computeForwardTransform(fftw_plan* plan)
{
  if(plan) {
    fftw_execute(*plan);
  } else {
    if(precomputedPlan_) {
      fftw_execute_dft_r2c(forwardPlan_, in_, out_);
    } else {
      fftw_execute(forwardPlan_);
    }
  }
}

/**.......................................................................
 * Compute the reverse transform
 */
void Dft2d::computeInverseTransform(fftw_plan* plan)
{
  if(plan) {
    fftw_execute(*plan);
  } else {

    if(inversePlan_ == 0)
      ThrowError("Gridder contains no data");

    if(precomputedPlan_) {
      fftw_execute_dft_c2r(inversePlan_, out_, in_);
    } else {
      fftw_execute(inversePlan_);
    }
  }

  normalizeTransform();
}

/**.......................................................................
 * Normalize, if requested
 */
void Dft2d::normalizeTransform()
{
  if(normalize_) {
    for(unsigned i=0; i < nInZeroPad_; i++) {
      in_[i] /= nInZeroPad_;
    }
  }
}

fftw_complex* Dft2d::getTransformDataPtr()
{
  return out_;
}

double* Dft2d::getImageDataPtr()
{
  return in_;
}

Image Dft2d::getImage(bool includeZeroPaddedRegion)
{
  Image image;
  unsigned imInd, dftInd;
  float min, max;
  bool first = true;

  if(!includeZeroPaddedRegion) {
    
    image.resize(xAxis_.logicalAxis_.getNpix(), yAxis_.logicalAxis_.getNpix());

    try {
      image.xAxis().setAngularSize(xAxis_.logicalAxis_.getAngularSize());
      image.yAxis().setAngularSize(yAxis_.logicalAxis_.getAngularSize());
    } catch(...) {
    }
    
    for(unsigned ix=0; ix < nx_; ix++) {
      for(unsigned iy=0; iy < ny_; iy++) {
	dftInd = (ix + xOffset_) * nyZeroPad_ + (iy + yOffset_);
	imInd  = iy * nx_ + ix;
	image.data_[imInd] = in_[dftInd];

	if(first) {
	  min = max = image.data_[imInd];
	  first = false;
	} else {
	  min = min < image.data_[imInd] ? min : image.data_[imInd];
	  max = max > image.data_[imInd] ? max : image.data_[imInd];
	}

      }
    }

  } else {

    image.resize(xAxis_.getNpix(), yAxis_.getNpix());

    try {
      image.xAxis().setAngularSize(xAxis_.getAngularSize());
      image.yAxis().setAngularSize(yAxis_.getAngularSize());
    } catch(...) {
    }

    for(unsigned ix=0; ix < nxZeroPad_; ix++) {
      for(unsigned iy=0; iy < nyZeroPad_; iy++) {
	dftInd = ix * nyZeroPad_ + iy;
	imInd  = iy * nxZeroPad_ + ix;
	image.data_[imInd] = in_[dftInd];

	if(first) {
	  min = max = image.data_[imInd];
	  first = false;
	} else {
	  min = min < image.data_[imInd] ? min : image.data_[imInd];
	  max = max > image.data_[imInd] ? max : image.data_[imInd];
	}

      }
    }
  }

  image.dataMin_ = min;
  image.dataMax_ = max;

  image.initializeRefPixFftConvention();

  if(hasAbsolutePosition_)
    image.setRaDecFft(ra_, dec_);

  image.hasData_  = true;

  image.units_    = units_;
  image.hasUnits_ = hasUnits_;

  return image;
}

/**.......................................................................
 * Perform a complex multiplication of this dft with another.
 *
 * Note:
 *
 *  F(f conv  g) =  F(f) x F(g)
 *
 * ie, the transform of the convolution is given by the complex
 * multiplication of the transforms, and
 *
 *  F(f xcorr g) = F*(f) x F(g)
 *
 * ie, the cross-correlation is given by the complex multiplication of
 * the conjugate transform of f with the tansform of g.
 */
void Dft2d::complexMultiply(Dft2d& dft, bool conjugate)
{
  checkConsistency(dft);

  double re1, re2;
  double im1, im2;
  int fac = conjugate ? - 1: 1;

  for(unsigned i=0; i < nOutZeroPad_; i++) {

    re1 =       out_[i][0];
    im1 = fac * out_[i][1];

    re2 = dft.out_[i][0];
    im2 = dft.out_[i][1];

    // Re = re1 * re2 - im1 * im2
    // Im = re1 * im2 + re2 * im1

    out_[i][0] = re1 * re2 - im1 * im2; 
    out_[i][1] = re1 * im2 + re2 * im1; 
  }
}

/**.......................................................................
 * Check consistency between two transforms
 */
void Dft2d::checkConsistency(Dft2d& dft)
{
  if(nOutZeroPad_ != dft.nOutZeroPad_) {
    ThrowError("Transforms are not of the same size");
  }
}

/**.......................................................................
 * Shift the data in this transform so that the inverse-transformed
 * image will be centered.
 */
void Dft2d::shift()
{
  int yfac, xfac, fac;

  for(unsigned iOut=0, ix=0; ix < nxZeroPad_; ix++) {
    xfac = ix%2==0 ? 1 : -1;
    for(unsigned iy=0; iy <= nyZeroPad_/2; iy++, iOut++) {
      yfac = iy%2==0 ? 1 : -1;
      fac = xfac * yfac;
      out_[iOut][0] *= fac;
      out_[iOut][1] *= fac;
    }
  }
}

/**.......................................................................
 * Shift this DFT by the specified angular offset:
 *
 *       f(x,y) = F(g(u,v))
 *
 *       f(x + dx, y + dy) = F(g(u, v) * exp(-2*pi*(u*dx + v*dy))
 *
 * So we transform:
 *
 *       g(u, v) --> g(u, v) * exp(-2*pi*(u*dx + v*dy))
 *
 * or 
 *
 *       re --> re * cos() - im * sin()
 *       im --> re * sin() + im * cos()
 *
 */
void Dft2d::shiftBy(Angle& xoff, Angle& yoff)
{
  double u, v, val, cosFac, sinFac, re, im;
  double xRad = xoff.radians();
  double yRad = yoff.radians();

  for(unsigned iOut=0, ix=0; ix < nxZeroPad_; ix++) {

    for(unsigned iy=0; iy <= nyZeroPad_/2; iy++, iOut++) {

      getUVData(ix, iy, DATA_REAL, u, v, val);

      cosFac = cos(-2*M_PI*(u * xRad + v * yRad));
      sinFac = sin(-2*M_PI*(u * xRad + v * yRad));

      re = out_[iOut][0];
      im = out_[iOut][1];

      out_[iOut][0] =  re*cosFac - im*sinFac;
      out_[iOut][1] =  re*sinFac + im*cosFac;
    }
  }
}

void Dft2d::plotInput(bool includeZeroPad)
{
  unsigned xOff = includeZeroPad ? 0 : xOffset_;
  unsigned yOff = includeZeroPad ? 0 : yOffset_;
  unsigned nx   = includeZeroPad ? nxZeroPad_ : nx_;
  unsigned ny   = includeZeroPad ? nyZeroPad_ : ny_;

  std::vector<double> arr(nx * ny);

  unsigned dftInd, imInd;
  for(unsigned ix=0; ix < nx; ix++) {
    for(unsigned iy=0; iy < ny; iy++) {
      dftInd = (ix + xOff) * nyZeroPad_ + (iy + yOff);
      imInd  = iy * nx + ix;
      arr[imInd] = in_[dftInd];
    }
  }

  PgUtil::greyScale(arr, nx, ny);
}

/**.......................................................................
 * Extract the U-row of the transform corresponding to this V index
 */
void Dft2d::plotRealU(unsigned iV)
{
  std::vector<double> xarr;
  std::vector<double> yarr;

  double spatialFrequencyResolution;
  std::string label;

  try {
    spatialFrequencyResolution = xAxis_.getSpatialFrequencyResolution();
    label = "U (inverse radians)";
  } catch(...) {
    spatialFrequencyResolution = 1.0;
    label = "U index";
  }

  getUData(iV, DATA_REAL, xarr, yarr);

  PgUtil::linePlot(xarr.size(), &xarr[0], &yarr[0], 0, (char*)label.c_str(), "Real", "");
}

/**.......................................................................
 * Extract the U-row of the transform corresponding to this V index
 */
void Dft2d::plotAbsU(unsigned iV)
{
  std::vector<double> xarr;
  std::vector<double> yarr;

  double spatialFrequencyResolution;
  std::string label;

  try {
    spatialFrequencyResolution = xAxis_.getSpatialFrequencyResolution();
    label = "U (inverse radians)";
  } catch(...) {
    spatialFrequencyResolution = 1.0;
    label = "U index";
  }

  getUData(iV, DATA_ABS, xarr, yarr);

  PgUtil::linePlot(xarr.size(), &xarr[0], &yarr[0], 0, (char*)label.c_str(), "Abs", "");
}

/**.......................................................................
 * Extract the U-row of the transform corresponding to this V index
 */
void Dft2d::plotImagU(unsigned iV)
{
  std::vector<double> xarr;
  std::vector<double> yarr;

  double spatialFrequencyResolution;
  std::string label;

  try {
    spatialFrequencyResolution = xAxis_.getSpatialFrequencyResolution();
    label = "U (inverse radians)";
  } catch(...) {
    spatialFrequencyResolution = 1.0;
    label = "U index";
  }

  getUData(iV, DATA_IMAG, xarr, yarr);

  PgUtil::linePlot(xarr.size(), &xarr[0], &yarr[0], 0, (char*)label.c_str(), "Imaginary", "");
}

/**.......................................................................
 * Extract the U-row of the transform corresponding to this V index
 *
 * Note the standard ordering of the DFT, where indices 0 -> n/2
 * correspond to positive frequencies, and indices n/2+1 -> n-1
 * correspond to negative frequencies.
 *
 * I.e., consider the transform of an n = 8 element array.  The
 * frequency resolution is given by 1/n, where the ordering is as
 * follows:
 *
 *                ------------------------------------------------
 *   index:      | 0  |  1  |  2  |  3  |   4   |  5  |  6  |  7  |
 *                ------------------------------------------------
 *               |    |   1 |   2 |   3 |    4  |   3 |   2 |   1 |
 *   freq:       |  0 | + - | + - | + - | +- -  | - - | - - | - - |
 *               |    |   n |   n |   n |    n  |   n |   n |   n |
 *                ------------------------------------------------
 *
 *  Note that the elements corresponding to the extremal frequencies
 *  (+- 4/n) are identical.
 *
 * In the zero-padded version of this array:
 *
 *                -------------------------------------------------------------------------------------------------
 *   index:      | 0  |  1  |  2  |  3  |  4  |  5  |  6  |  7  |   8   |  9  |  10 |  11 |  12  |  13 |  14 |  15 | 
 *                -------------------------------------------------------------------------------------------------
 *               |    |   1 |   2 |   3 |   4 |   5 |   6 |   7 |    8  |   7 |   6 |   5 |   4  |   3 |   2 |   1 |
 *   freq:       |  0 | + - | + - | + - | + - | + - | + - | + - | +- -  | - - | - - | - - | - -  | - - | - - | - - |
 *               |    |  2n |  2n |  2n |  2n |  2n |  2n |  2n |   2n  |  2n |  2n |  2n |  2n  |  2n |  2n |  2n |
 *                -------------------------------------------------------------------------------------------------
 *
 * And our positive frequencies correspond to 0 -> n/2, while negative frequencies correspond to 2n-1 -> 2n-1 - n/2,
 * or, if n == nx_ and 2n = nxZeroPad_, then nxZeroPad_-1 -> 
 */
void Dft2d::getUData(unsigned iV, DataType type, std::vector<double>& xarr, std::vector<double>& yarr)
{
  xarr.resize(nxZeroPad_+1);
  yarr.resize(nxZeroPad_+1);

  double spatialFrequencyResolution;

  try {
    spatialFrequencyResolution = xAxis_.getSpatialFrequencyResolution();
  } catch(...) {
    spatialFrequencyResolution = 1.0;
  }

  unsigned dftInd, imInd;

  // First the negative frequencies -- these are stored in the top
  // half of the array

  unsigned iU=0;
  for(unsigned iX=nxZeroPad_/2; iX < nxZeroPad_; iU++, iX++) {
    dftInd = iX * (nyZeroPad_/2+1) + iV;

    if(type == DATA_REAL) {
      yarr[iU] = out_[dftInd][0];
    } else if(type == DATA_IMAG) {
      yarr[iU] = out_[dftInd][1];
    } else {
      yarr[iU] = sqrt(out_[dftInd][0] * out_[dftInd][0] + out_[dftInd][1] * out_[dftInd][1]);
    }

    xarr[iU] = - ((double)(nxZeroPad_ - iX) * spatialFrequencyResolution);
  }

  // Now the positive frequencies -- these are stored in the first
  // half of the array

  for(unsigned iX=0; iX <= nxZeroPad_/2; iU++, iX++) {
    dftInd = iX * (nyZeroPad_/2+1) + iV;

    if(type == DATA_REAL) {
      yarr[iU] = out_[dftInd][0];
    } else if(type == DATA_IMAG) {
      yarr[iU] = out_[dftInd][1];
    } else {
      yarr[iU] = sqrt(out_[dftInd][0] * out_[dftInd][0] + out_[dftInd][1] * out_[dftInd][1]);
    }

    xarr[iU] = iX * spatialFrequencyResolution;
  }
}

/**.......................................................................
 * Extract the V-row of the transform corresponding to this U index.
 *
 * The standard ordering of the DFT has indices 0 -> n/2 corresponding
 * to positive frequencies.
 *
 * I.e., the v-axis of the transform of an n = 8 element array is
 * given by:
 *
 *                ------------------------------
 *   index:      | 0  |  1  |  2  |  3  |   4   |
 *                -------------------------------
 *               |    |   1 |   2 |   3 |    4  |
 *   freq:       |  0 | + - | + - | + - | +- -  |
 *               |    |   n |   n |   n |    n  |
 *                ------------------------------
 *
 *  Note that the elements corresponding to the extremal frequencies
 *  (+- 4/n) are identical.
 *
 * In the zero-padded version of this array:
 *
 *                ------------------------------------------------------
 *   index:      | 0  |  1  |  2  |  3  |  4  |  5  |  6  |  7  |   8   |
 *                ------------------------------------------------------
 *               |    |   1 |   2 |   3 |   4 |   5 |   6 |   7 |    8  |
 *   freq:       |  0 | + - | + - | + - | + - | + - | + - | + - | +- -  |
 *               |    |  2n |  2n |  2n |  2n |  2n |  2n |  2n |   2n  |
 *                -------------------------------------------------------
 *
 * Unlike the u-axis (see getUData() above), the negative frequencies
 * are completely determined by Hermitian symmetry F(-u, -v) = F*(u, v), 
 * whence:
 *
 *  F(u, -v) = F*(-u, v)
 *
 * Given the ordering of the data discussed in notes to getUData(),
 * the negative u frequency corresponding to index iU corresponds to
 * index n - iU, and we have schematically:
 *
 *  F(u, -v) = F*(N-iU, v)
 *  F(u,  v) = F(   iU, v)
 *
 */
void Dft2d::getVData(unsigned iU, DataType type, std::vector<double>& xarr, std::vector<double>& yarr)
{
  xarr.resize(nyZeroPad_+1);
  yarr.resize(nyZeroPad_+1);

  double spatialFrequencyResolution;

  try {
    spatialFrequencyResolution = yAxis_.getSpatialFrequencyResolution();
  } catch(...) {
    spatialFrequencyResolution = 1.0;
  }

  unsigned dftInd, imInd;

  // First the negative part of the transform

  unsigned iV = 0;

  unsigned iUNeg = (iU > 0) ? nxZeroPad_ - iU : iU;
  for(unsigned iY=0; iY <= nyZeroPad_/2; iY++, iV++) {
    dftInd = iUNeg * (nyZeroPad_/2+1) + (nyZeroPad_/2 - iY);

    if(type == DATA_REAL) {
      yarr[iV] =  out_[dftInd][0];
    } else if(type == DATA_IMAG) {
      yarr[iV] = -out_[dftInd][1];
    } else if(type == DATA_ABS) {
      yarr[iV] = sqrt(out_[dftInd][0] * out_[dftInd][0] + out_[dftInd][1] * out_[dftInd][1]);
    } else if(type == DATA_ABS_LOG) {
      yarr[iV] = log(sqrt(out_[dftInd][0] * out_[dftInd][0] + out_[dftInd][1] * out_[dftInd][1]));
    }

    xarr[iV] = - ((double)(nyZeroPad_/2 - iY) * spatialFrequencyResolution);
  }

  for(unsigned iY=0; iY <= nyZeroPad_/2; iY++, iV++) {
    dftInd = iU * (nyZeroPad_/2+1) + iY;

    if(type == DATA_REAL) {
      yarr[iV] =  out_[dftInd][0];
    } else if(type == DATA_IMAG) {
      yarr[iV] = -out_[dftInd][1];
    } else if(type == DATA_ABS) {
      yarr[iV] = sqrt(out_[dftInd][0] * out_[dftInd][0] + out_[dftInd][1] * out_[dftInd][1]);
    } else if(type == DATA_ABS_LOG) {
      yarr[iV] = log(sqrt(out_[dftInd][0] * out_[dftInd][0] + out_[dftInd][1] * out_[dftInd][1]));
    }

    xarr[iV] = iY * spatialFrequencyResolution;
  }
}

void Dft2d::getUVData(double u, double v, DataType type, double& uNearest, double& vNearest, double& val)
{
  unsigned iU, iV;
  bool conj;
  uvIndex(u, v, iU, iV, conj);

  getUVData(iU, iV, type, uNearest, vNearest, val);

  if(conj && type == DATA_IMAG)
    val = -val;
}

void Dft2d::getUVData(unsigned dftInd, DataType type, double& u, double& v, double& val)
{
  unsigned iU = dftInd / (nyZeroPad_/2+1);
  unsigned iV = dftInd - iU * (nyZeroPad_/2+1);

  return getUVData(iU, iV, type, u, v, val);
}

void Dft2d::getUVData(unsigned iU, unsigned iV, DataType type, double& u, double& v, double& val)
{
  double xSpatialFrequencyResolution;
  double ySpatialFrequencyResolution;
  
  if(xAxis_.hasAngularSize())
    xSpatialFrequencyResolution = xAxis_.getSpatialFrequencyResolution();
  else
    xSpatialFrequencyResolution = 1.0/xAxis_.getNpix();

  if(yAxis_.hasAngularSize())
    ySpatialFrequencyResolution = yAxis_.getSpatialFrequencyResolution();
  else
    ySpatialFrequencyResolution = 1.0/yAxis_.getNpix();

  unsigned dftInd = iU * (nyZeroPad_/2+1) + iV;

  if(type == DATA_REAL) {
    val = out_[dftInd][0];
  } else if(type == DATA_IMAG) {
    val = out_[dftInd][1];
  } else if(type == DATA_ABS) {
    val = sqrt(out_[dftInd][0] * out_[dftInd][0] + out_[dftInd][1] * out_[dftInd][1]);
  } else if(type == DATA_ABS_LOG) {
    val = log(sqrt(out_[dftInd][0] * out_[dftInd][0] + out_[dftInd][1] * out_[dftInd][1]));
  }

  if(iU > nxZeroPad_/2) {
    u = ((int) iU - (int) nxZeroPad_) * xSpatialFrequencyResolution;
  } else {
    u = iU * xSpatialFrequencyResolution;
  }

  v = iV * ySpatialFrequencyResolution;
}
  
void Dft2d::plotData2D(DataType type)
{
  std::vector<double> arr = getUVData(type);

  double xSpatialFrequencyResolution, xmin, xmax;
  double ySpatialFrequencyResolution, ymin, ymax;
  std::string xlabel, ylabel, title;

  try {
    xSpatialFrequencyResolution = xAxis_.getSpatialFrequencyResolution();
    xlabel = "U (inverse radians)";
  } catch(...) {
    xSpatialFrequencyResolution = 1.0/xAxis_.getNpix();
    xlabel = "U (inverse pixels)";
  }

  xmin = 0.0;
  xmax = xSpatialFrequencyResolution * nxZeroPad_;

  try {
    ySpatialFrequencyResolution = yAxis_.getSpatialFrequencyResolution();
    ylabel = "V (inverse radians)";
  } catch(...) {
    ySpatialFrequencyResolution = 1.0/yAxis_.getNpix();
    ylabel = "V (inverse pixels)";
  }

  ymin = 0.0;
  ymax = ySpatialFrequencyResolution * nyZeroPad_/2+1;

  switch(type) {
  case DATA_ABS_LOG:
    title = "log(Absolute)";
    break;
  case DATA_ABS:
    title = "Absolute";
    break;
  case DATA_REAL:
    title = "Real";
    break;
  case DATA_IMAG:
    title = "Imaginary";
    break;
  default:
    title = "";
    break;
  }

  PgUtil::setTitle(true);
  PgUtil::greyScale(arr.size(), &arr[0], nxZeroPad_, nyZeroPad_/2+1, xmin, xmax, ymin, ymax, 0, 0, 0, (char*)xlabel.c_str(), (char*)ylabel.c_str(), (char*)title.c_str());
}

/**.......................................................................
 * Plot the real component of the V-row of the transform corresponding
 * to this U index
 */
void Dft2d::plotRealV(unsigned iU)
  {
  std::vector<double> xarr;
  std::vector<double> yarr;

  double spatialFrequencyResolution;
  std::string label;

  try {
    spatialFrequencyResolution = yAxis_.getSpatialFrequencyResolution();
    label = "V (inverse radians)";
  } catch(...) {
    spatialFrequencyResolution = 1.0;
    label = "V index";
  }

  getVData(iU, DATA_REAL, xarr, yarr);

  PgUtil::linePlot(xarr.size(), &xarr[0], &yarr[0], 0, (char*)label.c_str(), "Real", "");
}

/**.......................................................................
 * Plot the abs component of the V-row of the transform corresponding
 * to this U index
 */
void Dft2d::plotAbsV(unsigned iU)
{
  std::vector<double> xarr;
  std::vector<double> yarr;

  double spatialFrequencyResolution;
  std::string label;

  try {
    spatialFrequencyResolution = yAxis_.getSpatialFrequencyResolution();
    label = "V (inverse radians)";
  } catch(...) {
    spatialFrequencyResolution = 1.0;
    label = "V index";
  }

  getVData(iU, DATA_ABS, xarr, yarr);

  PgUtil::linePlot(xarr.size(), &xarr[0], &yarr[0], 0, (char*)label.c_str(), "Real", "");
}

/**.......................................................................
 * Plot the real component of the V-row of the transform corresponding
 * to this U index
 */
void Dft2d::plotImagV(unsigned iU)
{
  std::vector<double> xarr;
  std::vector<double> yarr;

  double spatialFrequencyResolution;
  std::string label;

  try {
    spatialFrequencyResolution = yAxis_.getSpatialFrequencyResolution();
    label = "V (inverse radians)";
  } catch(...) {
    spatialFrequencyResolution = 1.0;
    label = "V index";
  }

  getVData(iU, DATA_IMAG, xarr, yarr);

  PgUtil::linePlot(xarr.size(), &xarr[0], &yarr[0], 0, (char*)label.c_str(), "Imaginary", "");
}

std::vector<double> Dft2d::getUVData(DataType type)
{
  std::vector<double> arr(nOutZeroPad_);
  double u, v;
  
  unsigned dftInd, imInd;
  for(unsigned ix=0; ix < nxZeroPad_; ix++) {
    for(unsigned iy=0; iy < nyZeroPad_/2+1; iy++) {
      dftInd = ix * (nyZeroPad_/2+1) + iy;
      imInd  = iy * nxZeroPad_   + ix;
      getUVData(dftInd, type, u, v, arr[imInd]);
    }
  }

  return arr;
}

/**.......................................................................
 * Interpolate the data in this container onto the requested (u,v)
 * point.
 *
 * Valid will record if the data could be interpolated
 */
void Dft2d::interpolateReImData(double u, double v, double& re, double& im, bool& valid)
{
  // Get the spatial frequency resolution of the axes

  int nu = xAxis_.getNpix();
  int nv = yAxis_.getNpix();

  double du;
  double dv;

  try {
    du = xAxis_.getSpatialFrequencyResolution();
  } catch(...) {
    du = 1.0/nu;
  }

  try {
    dv = yAxis_.getSpatialFrequencyResolution();
  } catch(...) {
    dv = 1.0/nv;
  }

  // Now initialize return values to zero

  re = 0.0;
  im = 0.0;

  valid = false;

  double s2    = 2*convSigInPixels_*convSigInPixels_;
  double wtSum = 0.0;

  // Get the UV index closest to the current point

  bool conjugateResult=false;
  unsigned uInd, vInd;

  uvIndex(u, v, uInd, vInd, conjugateResult, valid);

  if(!valid)
    return;

  // If the closest index is the conjugate point to the requested uv
  // coordinate, we must do all calculations relative to the conjugate
  // point

  if(conjugateResult) {
    u = -u;
    v = -v;
  }

  //------------------------------------------------------------
  // Get the UV coordinate of the nearest point
  //------------------------------------------------------------

  double uVal0, vVal0, uVal, vVal;
  uvCoord(uInd, vInd, uVal0, vVal0);

  // If we are closer than convMaskInPixels_ from the maximum spatial
  // frequency, treat the point as invalid, since we cannot do a
  // symmetric convolution

  long int uDelt = uInd-nu/2;
  long int vDelt = vInd-nv/2;
  if((abs(uDelt) < convMaskInPixels_) || (abs(vDelt) < convMaskInPixels_)) {
    valid = false;
    return;
  }

  // If the maximum spatial frequency was requested, check whether we
  // are after the negative or positive spatial frequency

  if(uInd == nu/2) {
    if((u < 0 && uVal0 > 0) || (u > 0 && uVal0 < 0)) {
      uVal0 = -uVal0;
    }
  }

  if(vInd == nv/2) {
    if((v < 0 && vVal0 > 0) || (v > 0 && vVal0 < 0)) {
      vVal0 = -vVal0;
    }
  }

  //------------------------------------------------------------
  // Now iterate over a mask of pixels centered on the nearest pixel
  //------------------------------------------------------------

  bool doConj=false;
  double reVal, imVal;

  // If we are on the edge of the frequency range sampled by the
  // array, only iterate in on direction.  If uInd == nu/2 or vInd ==
  // nv/2, getUVData() below has no way of distinguishing which
  // direction is valid to interpolate

  int iUStop = (uInd == nu/2) ? 0 : convMaskInPixels_;
  int iVStop = (vInd == nv/2) ? 0 : convMaskInPixels_;

  for(int iU = -convMaskInPixels_; iU <= iUStop; iU++) {
    for(int iV = -convMaskInPixels_; iV <= iVStop; iV++) {

      // Calculate the coordinate of the center of the current pixel

      uVal = uVal0 + iU * du;
      vVal = vVal0 + iV * dv;

      // Get the UV data of the index that is (iU, iV) away from the
      // nearest point.  

      try {
	getUVData(uInd, iU, vInd, iV, reVal, imVal);
      } catch(...) {
	continue;
      }

      // Calculate the delta of the requested (u,v) point from the
      // center of this pixel, in fractional pixels

      double duPix = (u - uVal)/du;
      double dvPix = (v - vVal)/dv;

      // And the weight corresponding to it

      double arg = (duPix * duPix + dvPix * dvPix)/s2;

      double wt = exp(-arg);

      // Accumulate running means for re and im
	
      re += ((reVal - re) * wt) / (wtSum + wt);
      im += ((imVal - im) * wt) / (wtSum + wt);

      wtSum += wt;

      // If any actual value went into the sum, mark it as valid

      valid = true;
    }
  }

  // Finally, if we were calculating sums for the conjugate point,
  // conjugate the result

  if(conjugateResult)
    im = -im;
}

/**.......................................................................
 * Return the complex value of the point at uInd + iU, vInd + iV
 */
void Dft2d::getUVData(int uInd, int iU, int vInd, int iV, double& re, double& im)
{
  int nu = xAxis_.getNpix();
  int nv = yAxis_.getNpix();

  int uIndex, vIndex;
  bool doConj=false;

  // If this is a positive u-frequency, and the offset would
  // position us beyond the highest frequency sampled by the transform

  if(uInd < nu/2 && (uInd + iU) > nu/2) {
    ThrowError("Requesting a U frequency > the largest frequency in the transform");
  }

  // If this is a negative u-frequency, and the offset would position us
  // beyond the highest negative frequency sampled by the transform

  if(uInd > nu/2 && (uInd + iU) < nu/2) {
    ThrowError("Requesting a U frequency < the largest negative frequency in the transform");
  }

  // If the v-offset would position us beyond the highest frequency
  // sampled by the transform

  if(abs(vInd+iV) > nv/2) {
    ThrowError("Requesting a V frequency > the largest frequency in the transform");
  }

  // Now account for wrapping of negative frequencies about the end of
  // the array

  uIndex = uInd+iU;

  if(uIndex < 0) {
    uIndex += nu;
  }

  if(uIndex > nu-1) {
    uIndex -= nu;
  }

  // And account for negative V points: F(u, -v) = F*(-u, v)

  vIndex = vInd+iV;
  if(vIndex < 0) {
    vIndex = abs(vIndex);

    if(uIndex != 0 && uIndex != nu/2) {
      uIndex = (nu - uIndex);
    }

    doConj = true;
  }

  // Now get the data

  unsigned dftInd = uIndex * (nv/2+1) + vIndex;

  re = out_[dftInd][0];
  im = out_[dftInd][1];

  if(doConj) 
    im = -im;
}

void Dft2d::plotReal()
{
  plotData2D(DATA_REAL);
}

void Dft2d::plotImag()
{
  plotData2D(DATA_IMAG);
}

void Dft2d::plotAbs()
{
  plotData2D(DATA_ABS);
}

void Dft2d::plotAbsLog()
{
  plotData2D(DATA_ABS_LOG);
}

/**........................................................................
 * Return the Fourier-space correlation length corresponding to the
 * requested correlation coefficient, for the specified diameter and
 * frequency.
 */
double Dft2d::correlationLength(Length& diameter, Frequency& freq, double correlationCoeff)
{
  return correlationLength(diameter, diameter, freq, correlationCoeff);
}

/**........................................................................
 * Return the Fourier-space correlation length corresponding to the
 * requested correlation coefficient, for the specified diameters and
 * frequency.
 */
double Dft2d::correlationLength(Length& diameter1, Length& diameter2, Frequency& freq, double correlationCoeff)
{
  // The fourier-space correlation is equivalent to the square of the
  // primary beam in the image plane.  This function has a sigma that
  // is sqrt(2) smaller than the primary beam in the case of the same
  // diameters, and in general is given by the product of the two
  // gaussians.

  double gs1 = Image::gaussianSigma(diameter1, freq).radians();
  double gs2 = Image::gaussianSigma(diameter2, freq).radians();
  double gs = gs1 * gs2 / sqrt(gs1 * gs1 + gs2 * gs2);

  // Since this is just another gaussian, the transform is also a
  // gaussian, with analytically calculable width:

  double fourierSigma = inverseGaussianWidth(gs);

  // Now the correlationLength is just the UV radius that corresponds
  // to the requested correlation coefficient:

  double length = sqrt(-log(correlationCoeff) * 2 * fourierSigma * fourierSigma);

  return length;
}

double Dft2d::inverseGaussianWidth(double gaussianWidth)
{
  return 1.0/(2 * M_PI * gaussianWidth);
}

/**.......................................................................
 * Return the x-axis descriptor of this image
 */
Dft2d::Axis& Dft2d::xAxis()
{
  return xAxis_;
}

/**.......................................................................
 * Return the y-axis descriptor of this image
 */
Dft2d::Axis& Dft2d::yAxis()
{
  return yAxis_;
}

/**.......................................................................
 * Fill the transform array with a uniform plane
 */
void Dft2d::createUniformDft(unsigned nx, unsigned ny, double value)
{
  resize(nx, ny);

  unsigned dftInd = 0;
  for(unsigned iV=0; iV < nyZeroPad_/2+1; iV++) {
    for(unsigned iU=0; iU < nxZeroPad_; iU++) {
      dftInd = iU * (nyZeroPad_/2+1) + iV;

      out_[dftInd][0] = isZeroPad(iU, iV) ? 0.0 : value;
      out_[dftInd][1] = 0.0;
    }
  }

  shift();
}

/**.......................................................................
 * Get the UV radius corresponding to the passsed indices
 */
void Dft2d::uvCoord(unsigned iU, unsigned iV, double& u, double& v)
{
  double val;
  getUVData(iU, iV, DATA_UV, u, v, val);
}

/**.......................................................................
 * Get the UV radius corresponding to the passsed indices
 */
double Dft2d::uvRadius(unsigned iU, unsigned iV)
{
  double u, v;
  uvCoord(iU, iV, u, v);
  return sqrt(u*u + v*v);
}

/**.......................................................................
 * Get the UV radius corresponding to the passsed indices
 */
void Dft2d::uvCoord(unsigned dftIndex, double& u, double& v)
{
  double val;
  getUVData(dftIndex, DATA_UV, u, v, val);
}

/**.......................................................................
 * Get the UV radius corresponding to the passsed indices
 */
double Dft2d::uvRadius(unsigned dftIndex)
{
  double u, v;
  uvCoord(dftIndex, u, v);
  return sqrt(u*u + v*v);
}

/**.......................................................................
 * Return the indices in our dft array corresponding to the specified uv point
 */
void Dft2d::uvIndex(double u, double v, unsigned& iU, unsigned& iV, bool& conj)
{
  int iu, iv;
  double uRes, vRes;
  conj = false;

  try {
    uRes = xAxis_.getSpatialFrequencyResolution();
  } catch(...) {
    uRes = 1.0/xAxis_.getNpix();
  }

  try {
    vRes = yAxis_.getSpatialFrequencyResolution();
  } catch(...) {
    vRes = 1.0/yAxis_.getNpix();
  }

  // Note that our array contains positive and negative U frequencies,
  // but only positive V frequencies

  if(v < 0) {
    u = -u;
    v = -v;
    conj = true;
  }

  // Convert to effective u and v indices

  iu = (int)floor(u/uRes + 0.5);
  iv = (int)floor(v/vRes + 0.5);

  if(abs(iu) > nxZeroPad_/2 || iv > nyZeroPad_/2) {
    ThrowColorError("\nUV point u = " << u << ", v = " << v << " lies outside the frequency span covered by this dft: "
		    << "(umax = " << (xAxis_.getSpatialFrequencyResolution() * nxZeroPad_/2) << ", "
		    << " vmax = " << (yAxis_.getSpatialFrequencyResolution() * nyZeroPad_/2) << ")", "red");
  }

  if(iu < 0)
    iu = nxZeroPad_ + iu;

  iU = iu;
  iV = iv;
}

void Dft2d::uvIndex(double u, double v, unsigned& iU, unsigned& iV, bool& conj, bool& valid)
{
  int iu, iv;
  double uRes, vRes;
  conj = false;

  try {
    uRes = xAxis_.getSpatialFrequencyResolution();
  } catch(...) {
    uRes = 1.0/xAxis_.getNpix();
  }

  try {
    vRes = yAxis_.getSpatialFrequencyResolution();
  } catch(...) {
    vRes = 1.0/yAxis_.getNpix();
  }

  // Note that our array contains positive and negative U frequencies,
  // but only positive V frequencies

  if(v < 0) {
    u = -u;
    v = -v;
    conj = true;
  }

  // Convert to effective u and v indices

  iu = (int)floor(u/uRes + 0.5);
  iv = (int)floor(v/vRes + 0.5);

  if(abs(iu) > nxZeroPad_/2 || iv > nyZeroPad_/2) {
    valid = false;
    return;
  }

  if(iu < 0)
    iu = nxZeroPad_ + iu;

  iU = iu;
  iV = iv;
  valid = true;
}

/**.......................................................................
 * Return the indices in our dft array corresponding to the specified uv point
 */
void Dft2d::dftIndex(double u, double v, unsigned& dftInd, bool& conj)
{
  unsigned iU, iV;

  uvIndex(u, v, iU, iV, conj);

  dftInd = iU * (nyZeroPad_/2+1) + iV;
}

/**.......................................................................
 * Return true if this point lies in the zeropadded part of the array
 */
bool Dft2d::isZeroPad(unsigned iU, unsigned iV)
{
  if(zeropad_)
    return (iU > nx_/2 && iU < nxZeroPad_ - 1 - nx_/2) || iV > ny_/2;
 else
   return false;
}

/**.......................................................................
 * Fill the transform array with a uniform disk
 */
void Dft2d::createUniformDiskDft(double radius)
{
  unsigned dftInd = 0;
  for(unsigned iV=0; iV < nyZeroPad_/2+1; iV++) {
    for(unsigned iU=0; iU < nxZeroPad_; iU++) {

      dftInd = iU * (nyZeroPad_/2+1) + iV;

      double r = uvRadius(iU, iV);

      if(isZeroPad(iU, iV)) {
	out_[dftInd][0] = 0.0;
	out_[dftInd][1] = 0.0;
      } else {
	out_[dftInd][0] = (r > radius) ? 0.0 : 1.0;
	out_[dftInd][1] = 0.0;
      }
    }
  }

  shift();
}

/**.......................................................................
 * Fill the transform array with a uniform disk
 */
void Dft2d::createUniformDiskDft(unsigned nx, unsigned ny, double radius)
{
  resize(nx, ny);
  createUniformDiskDft(radius);
}

/**.......................................................................
 * Fill the transform array with a uniform disk
 */
void Dft2d::createBlockedApertureUniformDiskDft(double innerRadius, double outerRadius)
{
  unsigned dftInd = 0;
  for(unsigned iV=0; iV < nyZeroPad_/2+1; iV++) {
    for(unsigned iU=0; iU < nxZeroPad_; iU++) {

      dftInd = iU * (nyZeroPad_/2+1) + iV;

      double r = uvRadius(iU, iV);

      if(isZeroPad(iU, iV)) {
	out_[dftInd][0] = 0.0;
	out_[dftInd][1] = 0.0;
      } else if(r < innerRadius) {
	out_[dftInd][0] = 0.0;
	out_[dftInd][1] = 0.0;
      } else {
	out_[dftInd][0] = (r > outerRadius) ? 0.0 : 1.0;
	out_[dftInd][1] = 0.0;
      }
    }
  }

  shift();
}

/**.......................................................................
 * Fill the transform array with a uniform disk
 */
void Dft2d::createBlockedApertureUniformDiskDft(unsigned nx, unsigned ny, double innerRadius, double outerRadius)
{
  resize(nx, ny);
  createBlockedApertureUniformDiskDft(innerRadius, outerRadius);
}

/**.......................................................................
 * Fill the transform array with a tapered disk, a la Rohlfs and
 * Wilson, 6.34
 */
void Dft2d::createTaperedApertureDft(Wavelength& lambda, Length& diameter, double K, unsigned p)
{
  // Calculate the UV radius corresponding to the effective edge of
  // the dish

  double rhoMax = (diameter / lambda) / 2;

  unsigned dftInd = 0;
  for(unsigned iV=0; iV < nyZeroPad_/2+1; iV++) {
    for(unsigned iU=0; iU < nxZeroPad_; iU++) {
      
      dftInd = iU * (nyZeroPad_/2+1) + iV;
      
      double r = uvRadius(iU, iV) / rhoMax;
      
      if(isZeroPad(iU, iV)) {
	out_[dftInd][0] = 0.0;
	out_[dftInd][1] = 0.0;
      } else {
	out_[dftInd][0] = (r > 1.0) ? 0.0 : (K + pow( (1.0 - (r*r)), (double)p));
	out_[dftInd][1] = 0.0;
      }
    }
  }
  
  shift();
}

/**.......................................................................
 * Fill the transform array with a tapered disk, a la Rohlfs and
 * Wilson, 6.34
 */
void Dft2d::createTaperedApertureDft(unsigned nx, unsigned ny, Wavelength& lambda, Length& diameter, double K, unsigned p)
{
  resize(nx, ny);
  createTaperedApertureDft(lambda, diameter, K, p);
}

void Dft2d::updateAxis(unsigned axis)
{
  unsigned axesMask = (unsigned)axes_;
  unsigned axisMask = (unsigned)axis;
  axes_ = (axesMask | axisMask);

  checkAxes();
}

void Dft2d::checkAxes()
{
  if(axes_ == ImageAxis::AXIS_BOTH)
    resize();
}

//-----------------------------------------------------------------------
// Dft2d::Axis keeps track of the actual size of the transform arrays,
// including zeropadding
//-----------------------------------------------------------------------

Dft2d::Axis::Axis()
{
  zeropadFactor_ = 1;
}

/**.......................................................................
 * Change the zeropadding factor for this axis
 */
void Dft2d::Axis::setZeropadFactor(unsigned zeropadFactor)
{
  zeropadFactor_ = zeropadFactor;

  try {
    Axis::setNpix(logicalAxis_.getNpix() * zeropadFactor_);
    Axis::setAngularSize(logicalAxis_.getAngularSize() * zeropadFactor_);
  } catch(...) {
  }

  parent_->checkAxes();
}

/**.......................................................................
 * Set the number of (logical) pixels in this axis
 */
void Dft2d::Axis::setNpix(unsigned n)
{
  logicalAxis_.setNpix(n);
  ImageAxis::setNpix(n * zeropadFactor_);
  parent_->updateAxis(type_);
}

/**.......................................................................
 * Set the (logical) size of this axis
 */
void Dft2d::Axis::setAngularSize(Angle size)
{
  logicalAxis_.setAngularSize(size);
  ImageAxis::setAngularSize(size * zeropadFactor_);
}

/**.......................................................................
 * Set the (logical) spatial frequency resolution  of this axis
 */
void Dft2d::Axis::setSpatialFrequencyResolution(double inverseRadians)
{
  logicalAxis_.setSpatialFrequencyResolution(inverseRadians);
  ImageAxis::setSpatialFrequencyResolution(inverseRadians / zeropadFactor_);
}

/**.......................................................................
 * Get the actual spatial frequency resolution of this axis
 */
double Dft2d::Axis::getSpatialFrequencyResolution()
{
  return logicalAxis_.getSpatialFrequencyResolution() / zeropadFactor_;
}

void Dft2d::Axis::operator=(ImageAxis& axis)
{
  //  COUT(pthread_self() << " Inside Dft2d::ImageAxis oprator= 0");

  if(axis.hasAngularSize()) {
    setAngularSize(axis.getAngularSize());
  }

  if(axis.hasNpix()) {
    setNpix(axis.getNpix());
  }

  //  COUT(pthread_self() << " Inside Dft2d::ImageAxis oprator= 1");
}

void Dft2d::Axis::operator=(const Dft2d::Axis& axis)
{
  *this = (Dft2d::Axis&)axis;
}

void Dft2d::Axis::operator=(Dft2d::Axis& axis)
{
  //  COUT("Inside Dft2d::Axis oprator= 0");
  zeropadFactor_ = axis.zeropadFactor_;

  if(axis.hasAngularSize())
    setAngularSize(axis.getAngularSize());

  if(axis.hasNpix())
    setNpix(axis.getNpix());

  //  COUT("Inside Dft2d::Axis oprator= 1");
}

/**.......................................................................
 * Fill the transform array with a bessel function, tapered to a
 * specific point in the main lobe of the bessel function
 */
void Dft2d::createJ0BesselFnDft(unsigned nx, unsigned ny, Wavelength& lambda, Length& diameter, double fractionalPower)
{
  resize(nx, ny);
  createJ0BesselFnDft(lambda, diameter, fractionalPower);
}

void Dft2d::createJ0BesselFnDft(Wavelength& lambda, Length& diameter, double fractionalPower)
{
  // Binary search to find the radius corresponding to the specified
  // fractional power point

  double rhoEdge = j0BesselFnRadius(fractionalPower);
  double rhoMax  = (diameter / lambda) / 2;

  unsigned dftInd = 0;
  for(unsigned iV=0; iV < nyZeroPad_/2+1; iV++) {
    for(unsigned iU=0; iU < nxZeroPad_; iU++) {
      
      dftInd = iU * (nyZeroPad_/2+1) + iV;
      
      double r = uvRadius(iU, iV) / rhoMax;
      
      if(isZeroPad(iU, iV)) {
	out_[dftInd][0] = 0.0;
	out_[dftInd][1] = 0.0;
      } else {
	out_[dftInd][0] = (r > 1.0) ? 0.0 : j0BesselFn(r * rhoEdge);
	out_[dftInd][1] = 0.0;
      }
    }
  }
  
  shift();
}

/**.......................................................................
 * Fill the transform array with a bessel function, tapered to a
 * specific point in the main lobe of the bessel function
 */
void Dft2d::createBlockedApertureJ0BesselFnDft(unsigned nx, unsigned ny, Wavelength& lambda, 
					       Length& innerDiameter, Length& outerDiameter, 
					       double fractionalPower)
{
  resize(nx, ny);
  createBlockedApertureJ0BesselFnDft(lambda, innerDiameter, outerDiameter, fractionalPower);
}

void Dft2d::createBlockedApertureJ0BesselFnDft(Wavelength& lambda, 
					       Length& innerDiameter, Length& outerDiameter,
					       double fractionalPower)
{
  // Binary search to find the radius corresponding to the specified
  // fractional power point

  double rhoEdge = j0BesselFnRadius(fractionalPower);
  double rhoMin  = (innerDiameter / lambda) / 2;
  double rhoMax  = (outerDiameter / lambda) / 2;

  unsigned dftInd = 0;
  for(unsigned iV=0; iV < nyZeroPad_/2+1; iV++) {
    for(unsigned iU=0; iU < nxZeroPad_; iU++) {
      
      dftInd = iU * (nyZeroPad_/2+1) + iV;
      
      double rho = uvRadius(iU, iV);

      if(isZeroPad(iU, iV)) {

	out_[dftInd][0] = 0.0;
	out_[dftInd][1] = 0.0;

      } else {

	if(rho <= rhoMin) {
	  out_[dftInd][0] = 0.0;
	} else {
	  out_[dftInd][0] = (rho > rhoMax) ? 0.0 : j0BesselFn(rho * (rhoEdge/rhoMax));
	}

	out_[dftInd][1] = 0.0;
      }
    }
  }
  
  shift();
}

/**.......................................................................
 * Compute the zero-th order Bessel function (NR routine)
 */
double Dft2d::j0BesselFn(double x)
{
  double ax,z;
  double xx,y,ans,ans1,ans2;
  
  if ((ax=fabs(x)) < 8.0) {
    y=x*x;
    ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
					    +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
    ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
					  +y*(59272.64853+y*(267.8532712+y*1.0))));
    ans=ans1/ans2;
  } else {
    z=8.0/ax;
    y=z*z;
    xx=ax-0.785398164;
    ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
				    +y*(-0.2073370639e-5+y*0.2093887211e-6)));
    ans2 = -0.1562499995e-1+y*(0.1430488765e-3
			       +y*(-0.6911147651e-5+y*(0.7621095161e-6
						       -y*0.934935152e-7)));
    ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  }
  return ans;
}

/**.......................................................................
 * Compute the point at which the square of the J0 Bessel function is
 * equal to the requested value
 */
double Dft2d::j0BesselFnRadius(double val, double eps)
{
  double firstNull=2.4048255;     // first zero of Bessel fcn J0

  double low  = 0.0;
  double high = firstNull;
  double mid;
  double j0b2;

  do {
    mid = (low + high) / 2;
    j0b2 = j0BesselFn(mid);
    j0b2 *= j0b2;

    if(j0b2 > val) {
      low = mid;
    } else if(j0b2 <= val) {
      high = mid;
    }

  } while(fabs(j0b2 - val) > eps);

  return mid;
}

void Dft2d::setNpix(unsigned npix)
{
  xAxis().setNpix(npix);
  yAxis().setNpix(npix);
}

void Dft2d::setAngularSize(Angle size)
{
  xAxis().setAngularSize(size);
  yAxis().setAngularSize(size);
}

void Dft2d::setSpatialFrequencyResolution(double inverseRadians)
{
  xAxis().setSpatialFrequencyResolution(inverseRadians);
  yAxis().setSpatialFrequencyResolution(inverseRadians);
}

unsigned Dft2d::nearestPowerOf2NotLessThan(double val)
{
  return (unsigned)pow(2.0, ceil(log(val) / log(2.0)));
}

unsigned Dft2d::getTransformLength()
{
  return nOutZeroPad_;
}

void Dft2d::setPlan(Dft2d* dft)
{
  if(dft) {
    fwdPlanTmp_ = &dft->forwardPlan_;
    invPlanTmp_ = &dft->inversePlan_;
  } else {
    fwdPlanTmp_ = 0;
    invPlanTmp_ = 0;
  }
}

/**.......................................................................
 * Apply a gaussian taper to the high-frequency side of this Dft
 */
void Dft2d::lowPass(double uvrPeak, double uvrSigma)
{
  unsigned dftInd = 0;
  for(unsigned iV=0; iV < nyZeroPad_/2+1; iV++) {
    for(unsigned iU=0; iU < nxZeroPad_; iU++) {

      dftInd = iU * (nyZeroPad_/2+1) + iV;

      double r = uvRadius(iU, iV);

      if(r > uvrPeak) {
	double dr    = (r - uvrPeak)/uvrSigma;
	double taper = exp(-0.5*(dr*dr));
	out_[dftInd][0] = out_[dftInd][0] * taper;
	out_[dftInd][1] = out_[dftInd][1] * taper;
      }

    }
  }
}

/**.......................................................................
 * Apply a gaussian taper to the low-frequency side of this Dft
 */
void Dft2d::highPass(double uvrPeak, double uvrSigma)
{
  unsigned dftInd = 0;
  for(unsigned iV=0; iV < nyZeroPad_/2+1; iV++) {
    for(unsigned iU=0; iU < nxZeroPad_; iU++) {

      dftInd = iU * (nyZeroPad_/2+1) + iV;

      double r = uvRadius(iU, iV);

      if(r < uvrPeak) {
	double dr    = (r - uvrPeak)/uvrSigma;
	double taper = exp(-0.5*(dr*dr));
	out_[dftInd][0] = out_[dftInd][0] * taper;
	out_[dftInd][1] = out_[dftInd][1] * taper;
      }

    }
  }
}

/**.......................................................................
 * Apply a gaussian notch
 */
void Dft2d::notch(double uvrCenter, double uvrSigma)
{
  unsigned dftInd = 0;
  for(unsigned iV=0; iV < nyZeroPad_/2+1; iV++) {
    for(unsigned iU=0; iU < nxZeroPad_; iU++) {

      dftInd = iU * (nyZeroPad_/2+1) + iV;

      double r = uvRadius(iU, iV);

      double dr    = (r - uvrCenter)/uvrSigma;
      double taper = 1.0 - exp(-0.5*(dr*dr));
      out_[dftInd][0] = out_[dftInd][0] * taper;
      out_[dftInd][1] = out_[dftInd][1] * taper;

    }
  }
}

ImageAxis& Dft2d::xImageAxis()
{
  return xAxis_;
}

ImageAxis& Dft2d::yImageAxis()
{
  return yAxis_;
}

void Dft2d::operator*=(double mult)
{
  if(!hasData_) {
    ThrowError("Dft contains no data");
  }

  for(unsigned i=0; i < nOutZeroPad_; i++) {
    out_[i][0] *= mult;
    out_[i][1] *= mult;
  }
}

void Dft2d::setPlan(fftw_plan forwardPlan, fftw_plan inversePlan)
{
  precomputedPlan_ = true;
  forwardPlan_ = forwardPlan;
  inversePlan_ = inversePlan;
}

void Dft2d::initializeForVis(Image& image)
{
  setUnits(Unit::UNITS_JYBEAM);
  zeropad(true);

  xAxis().setNpix(image.xAxis().getNpix());
  yAxis().setNpix(image.yAxis().getNpix());
  xAxis().setAngularSize(image.xAxis().getAngularSize());
  yAxis().setAngularSize(image.yAxis().getAngularSize());
}

void Dft2d::initializeForVis(Angle& xSize, Angle& ySize, unsigned nXPix, unsigned nYPix)
{
  setUnits(Unit::UNITS_JYBEAM);
  zeropad(true);

  xAxis().setNpix(nXPix);
  yAxis().setNpix(nYPix);
  xAxis().setAngularSize(xSize);
  yAxis().setAngularSize(ySize);
}

void Dft2d::zero()
{
  unsigned nInd = nOutZeroPad_;
  unsigned dftInd;

  for(unsigned i=0; i < nInd; i++) {
    out_[i][0] = 0.0;
    out_[i][1] = 0.0;
  }
}
