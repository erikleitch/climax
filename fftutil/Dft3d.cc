#include "gcp/fftutil/Dft3d.h"
#include "gcp/pgutil/PgUtil.h"

#include <vector>

using namespace std;
using namespace gcp::util;

/**.......................................................................
 * Constructors
 */
Dft3d::Dft3d(bool optimize)
{
  initialize();
  optimize_  = optimize;
}

void Dft3d::initialize() {
  in_          = 0;
  out_         = 0;

  optimize_    = true;
  normalize_   = false;

  nx_          = 0;
  ny_          = 0;
  nz_          = 0;
  nIn_         = 0;
  nOut_        = 0;
}

Dft3d::Dft3d(int nx, int ny, int nz, bool optimize) 
{
  initialize();
  optimize_  = optimize;
  normalize_ = false;
  resize(nx, ny, nz);
}

Dft3d::Dft3d(Image& image, bool optimize) 
{
  initialize();
  optimize_  = optimize;
  normalize_ = false;
  initialize(image);  
}

/**.......................................................................
 * Destructor.
 */
Dft3d::~Dft3d() {}

/**.......................................................................
 * Compute a plan for this fft
 */
void Dft3d::computePlan(unsigned flag) 
{
  forwardPlan_ = fftw_plan_dft_r2c_3d(nx_, ny_, nz_, in_,  out_, flag);
  inversePlan_ = fftw_plan_dft_c2r_3d(nx_, ny_, nz_, out_, in_,  flag);
}

/**.......................................................................
 * Initialize the input data to this transform from an image
 */
void Dft3d::initialize(Image& image)
{
  resize(16, image.xAxis().getNpix(), image.yAxis().getNpix());

  unsigned imInd, dftInd;
  for(unsigned ix=0; ix < nx_; ix++) {
    for(unsigned iy=0; iy < ny_; iy++) {
      for(unsigned iz=0; iz < nz_; iz++) {
      dftInd = (ix * ny_ + iy) * nz_ + iz;
      imInd  = iz * ny_ + iy;
      in_[dftInd] = image.data_[imInd];
      }
    }
  }
}

/**.......................................................................
 * Resize for FFT of a different length
 */
void Dft3d::resize(unsigned nx, unsigned ny, unsigned nz)
  {
  nx_   = nx;
  ny_   = ny;
  nz_   = nz;
  nIn_  = nx_ * ny_ * nz_;
  nOut_ = nx_ * ny_ * (nz_/2+1);

  if(in_) {
    fftw_free(in_);
    in_ = 0;
  }

  if(out_) {
    fftw_free(out_);
    out_ = 0;
  }

  // Allocate arrays

  COUT("About to allocate array: nIn = " << nIn_);

  if((in_ = (double*)fftw_malloc(nIn_ * sizeof(double)))==0)
    ThrowError("Couldn't allocate input data array");

  COUT("About to allocate array: nOut = " << nOut_);

  if((out_ = (fftw_complex*)fftw_malloc(nOut_ * sizeof(fftw_complex)))==0)
    ThrowError("Couldn't allocate output data array");

  for(unsigned i=0; i < nIn_; i++)
    in_[i] = 0.0;

  for(unsigned i=0; i < nOut_; i++) {
    out_[i][0] = 0.0;
    out_[i][1] = 0.0;
  }

  // Now compute a plan for this machine architecture
      
  computePlan(optimize_ ? FFTW_MEASURE : FFTW_ESTIMATE);
}

/**.......................................................................
 * Compute the forward transform
 */
void Dft3d::computeForwardTransform()
{
  //  removeMean();
  fftw_execute(forwardPlan_);
}

/**.......................................................................
 * Compute the reverse transform
 */
void Dft3d::computeInverseTransform()
{
  fftw_execute(inversePlan_);
}


