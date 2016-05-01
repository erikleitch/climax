#include "gcp/fftutil/UvDataGridder.h"
#include "gcp/pgutil/PgUtil.h"

#include <vector>

using namespace std;
using namespace gcp::util;

#define EMLTEST

UvDataGridder::UvDataGridder(bool optimize) : Dft2d(optimize)
{
  nPt_.resize(0);
  wtSum_.resize(0);
  wt2Sum_.resize(0);

  store_                     = 0;
  errorInMean_               = 0;
  estimateErrInMeanFromData_ = false;
  errorInMeanIsValid_        = false;
  debugPrint_                = false;
}

UvDataGridder::UvDataGridder(const UvDataGridder& gridder)
{
  *this = (UvDataGridder&) gridder;
}

UvDataGridder::UvDataGridder(UvDataGridder& gridder)
{
  *this = gridder;
}

void UvDataGridder::operator=(const UvDataGridder& gridder)
{
  *this = (UvDataGridder&) gridder;
}

void UvDataGridder::operator=(UvDataGridder& gridder)
{
  Dft2d::operator=(gridder);

  if(errorInMean_) {
    fftw_free(errorInMean_);
    errorInMean_ = 0;
  }

  if(store_) {
    fftw_free(store_);
    store_ = 0;
  }
}

void UvDataGridder::sizeToMatch(UvDataGridder& gridder)
{
  COUT("Insize sizeToMatch");
  operator=(gridder);
  resize();
  COUT("leaving Insize sizeToMatch");
}

void UvDataGridder::duplicate(UvDataGridder& gridder)
{
  COUT("Inside duplicate: this = " << this << " gridd =  " << &gridder << " uu = " << gridder.uuSum_);

  operator=(gridder);
  resize();
  assignDataFrom(gridder);


  wtSumTotal_ = gridder.wtSumTotal_;
  uuSum_      = gridder.uuSum_;
  vvSum_      = gridder.vvSum_;
  uvSum_      = gridder.uvSum_;

  wtSum_      = gridder.wtSum_;
  wt2Sum_     = gridder.wt2Sum_;

  populatedIndices_ = gridder.populatedIndices_;
  populatedU_       = gridder.populatedU_;
  populatedV_       = gridder.populatedV_;

  errorInMeanIsValid_ = gridder.errorInMeanIsValid_;
}

void UvDataGridder::operator+=(const UvDataGridder& gridder)
{
  *this += (UvDataGridder&) gridder;
}

void UvDataGridder::operator+=(UvDataGridder& gridder)
{
  if(nOut_ != gridder.nOut_ || populatedIndices_.size() != gridder.populatedIndices_.size()) {
    ThrowColorError("Cannot assign data from dfts of different sizes", "red");
  }

  unsigned nInd = gridder.populatedIndices_.size();
  unsigned dftInd;

  for(unsigned i=0; i < nInd; i++) {
    dftInd = gridder.populatedIndices_[i];
    out_[dftInd][0] += gridder.out_[dftInd][0];
    out_[dftInd][1] += gridder.out_[dftInd][1];
  }
}

void UvDataGridder::assignDataFrom(const UvDataGridder& gridder)
{
  assignDataFrom((UvDataGridder&) gridder);
}

void UvDataGridder::assignDataFrom(UvDataGridder& gridder)
{
  if(gridder.hasData_) {
    if(nOut_ != gridder.nOut_ || populatedIndices_.size() != gridder.populatedIndices_.size()) {
      ThrowColorError("Cannot assign data from dfts of different sizes", "red");
    }
    
    unsigned nInd = gridder.populatedIndices_.size();
    unsigned dftInd;
    
    for(unsigned i=0; i < nInd; i++) {
      dftInd = gridder.populatedIndices_[i];

      out_[dftInd][0] = gridder.out_[dftInd][0];
      out_[dftInd][1] = gridder.out_[dftInd][1];

      errorInMean_[dftInd][0] = gridder.errorInMean_[dftInd][0];
      errorInMean_[dftInd][1] = gridder.errorInMean_[dftInd][1];

      store_[dftInd][0] = gridder.store_[dftInd][0];
      store_[dftInd][1] = gridder.store_[dftInd][1];
    }

    hasData_ = true;
  }
}

void UvDataGridder::setPopulatedIndicesToValue(double val)
{
  unsigned nInd = populatedIndices_.size();
  unsigned dftInd;
  
  for(unsigned i=0; i < nInd; i++) {
    dftInd = populatedIndices_[i];
    out_[dftInd][0] = val;
    out_[dftInd][1] = 0.0;
  }
  
  hasData_ = true;
}

UvDataGridder::~UvDataGridder() 
{
  if(errorInMean_) {
    fftw_free(errorInMean_);
    errorInMean_ = 0;
  }

  if(store_) {
    fftw_free(store_);
    store_ = 0;
  }
}

void UvDataGridder::resize()
{
  Dft2d::resize();

  nPt_.resize(nOutZeroPad_);
  wtSum_.resize(nOutZeroPad_);
  wt2Sum_.resize(nOutZeroPad_);

  if(errorInMean_) {
    fftw_free(errorInMean_);
    errorInMean_ = 0;
  }

  if((errorInMean_ = (fftw_complex*)fftw_malloc(nOutZeroPad_ * sizeof(fftw_complex)))==0) {
    ThrowError("Couldn't allocate output data array");
  }

  if(store_) {
    fftw_free(store_);
    store_ = 0;
  }

  if((store_ = (fftw_complex*)fftw_malloc(nOutZeroPad_ * sizeof(fftw_complex)))==0) {
    ThrowError("Couldn't allocate output data array");
  }
}

/**.......................................................................
 * Reinitialize for first-moment calculation
 */
void UvDataGridder::initializeForFirstMoments()
{
  populatedIndices_.resize(0);

  for(unsigned i=0; i < nInZeroPad_; i++) {
    in_[i] = 0.0;
    in_[i] = 0.0;
  }

  for(unsigned i=0; i < nOutZeroPad_; i++) {
    out_[i][0] = 0.0;
    out_[i][1] = 0.0;
    wtSum_[i]  = 0.0;
    nPt_[i]    = 0;
  }

  wtSumTotal_ = 0.0;

  uuSum_ = 0.0;
  vvSum_ = 0.0;
  uvSum_ = 0.0;
}

/**.......................................................................
 * Reinitialize for second-moment calculation
 */
void UvDataGridder::initializeForSecondMoments()
{
  for(unsigned i=0; i < nOutZeroPad_; i++) {
    errorInMean_[i][0] = 0.0;
    errorInMean_[i][1] = 0.0;
    wtSum_[i]          = 0.0;
    wt2Sum_[i]         = 0.0;
    nPt_[i]            = 0;
  }

  wtSumTotal_ = 0.0;
}

/**.......................................................................
 * Add a visibility to the running mean of weighted visibilities.
 * 
 * Here, by weight I do mean weight, in the sense that the mean is given by:
 *
 *         sum(wt * val)
 *   mn = ---------------
 *            sum(wt)
 *
 * So do NOT pass a variance instead.  The variance is generally the
 * inverse of what I mean by the weight.
 */
void UvDataGridder::accumulateFirstMoments(double u, double v, double re, double im, double wt)
{
  unsigned dftInd;
  bool conj=false;

  //------------------------------------------------------------
  // Get the index in our dft array corresponding to this uv point
  //------------------------------------------------------------

  dftIndex(u, v, dftInd, conj);

  //------------------------------------------------------------
  // Conjugate the data if required, to convert from (u, v) to (-u, -v)
  //------------------------------------------------------------

  if(conj)
    im = -im;

  //------------------------------------------------------------
  // Recompute the running mean
  //------------------------------------------------------------

  double remean = out_[dftInd][0];
  double immean = out_[dftInd][1];

  out_[dftInd][0] += (re - remean) * wt / (wtSum_[dftInd] + wt);
  out_[dftInd][1] += (im - immean) * wt / (wtSum_[dftInd] + wt);

  if(debugPrint_) {
    COUTCOLOR("dftInd = " << dftInd << " wtSum = " << wtSum_[dftInd] << " wt = " << wt, "cyan");
  }

  //------------------------------------------------------------
  // Increment weight sums needed for running averages and errors
  //------------------------------------------------------------

  nPt_[dftInd]++;
  wtSum_[dftInd]  += wt;
  wt2Sum_[dftInd] += wt*wt;
  wtSumTotal_     += wt;

  //------------------------------------------------------------ 
  // Also compute moments we will use to estimate the synthesized beam
  //------------------------------------------------------------

  uuSum_ += wt * (u*u - uuSum_) / wtSumTotal_;
  vvSum_ += wt * (v*v - vvSum_) / wtSumTotal_;
  uvSum_ += wt * (u*v - uvSum_) / wtSumTotal_;
}

void UvDataGridder::accumulateFirstMomentsTest(double u, double v, double re, double im, double wt)
{
  unsigned dftInd;
  bool conj=false;

  //------------------------------------------------------------
  // Get the index in our dft array corresponding to this uv point
  //------------------------------------------------------------

  dftIndex(u, v, dftInd, conj);

  //------------------------------------------------------------
  // Conjugate the data if required, to convert from (u, v) to (-u, -v)
  //------------------------------------------------------------

  if(conj)
    im = -im;

  //------------------------------------------------------------
  // Recompute the running mean
  //------------------------------------------------------------

  double remean = out_[dftInd][0];
  double immean = out_[dftInd][1];

  out_[dftInd][0] = (remean * wtSum_[dftInd] + re * wt) / (wtSum_[dftInd] + wt);
  out_[dftInd][1] = (immean * wtSum_[dftInd] + im * wt) / (wtSum_[dftInd] + wt);

  //------------------------------------------------------------
  // Increment weight sums needed for running averages and errors
  //------------------------------------------------------------

  nPt_[dftInd]++;
  wtSum_[dftInd]  += wt;
  wt2Sum_[dftInd] += wt*wt;
  wtSumTotal_     += wt;

  if(debugPrint_) {
    COUTCOLOR("dftInd = " << dftInd << " wtSum = " << wtSum_[dftInd] << " total = " << wtSumTotal_, "cyan");
  }

  //------------------------------------------------------------ 
  // Also compute moments we will use to estimate the synthesized beam
  //------------------------------------------------------------

  uuSum_ += wt * (u*u - uuSum_) / wtSumTotal_;
  vvSum_ += wt * (v*v - vvSum_) / wtSumTotal_;
  uvSum_ += wt * (u*v - uvSum_) / wtSumTotal_;
}

/**.......................................................................
 * Add a visibility to the running mean of weighted visibilities.
 * 
 * Here, by weight I mean weight, in the sense that the mean is given by:
 *
 *         sum(wt * val)
 *   mn = ---------------
 *            sum(wt)
 *
 * So do NOT pass a variance instead, for example.
 */
void UvDataGridder::accumulateFirstMoments(unsigned dftInd, double re, double im, double wt)
{
  // Recompute the running mean

  double remean = out_[dftInd][0];
  double immean = out_[dftInd][1];

  out_[dftInd][0] += (re - remean) * wt / (wtSum_[dftInd] + wt);
  out_[dftInd][1] += (im - immean) * wt / (wtSum_[dftInd] + wt);

  // Increment weight sums needed for running averages and errors

  nPt_[dftInd]++;
  wtSum_[dftInd]  += wt;
  wt2Sum_[dftInd] += wt*wt;
  wtSumTotal_     += wt;

  // If this is the first point in this cell, add to the populated index array

  if(nPt_[dftInd] == 1)
    populatedIndices_.push_back(dftInd);

  // Also compute moments we will use to estimate the synthesized beam

  double u, v, val;
  getUVData(dftInd, DATA_UV, u, v, val);

  uuSum_ += wt * (u*u - uuSum_) / wtSumTotal_;
  vvSum_ += wt * (v*v - vvSum_) / wtSumTotal_;
  uvSum_ += wt * (u*v - uvSum_) / wtSumTotal_;
}

/**.......................................................................
 * Add a visibility to the running mean of weighted visibilities.
 * 
 * Here, by weight I mean weight, in the sense that the mean is given by:
 *
 *         sum(wt * val)
 *   mn = ---------------
 *            sum(wt)
 *
 * So do NOT pass a variance instead, for example, which is generally
 * the inverse of what I mean by the weight.
 */
void UvDataGridder::accumulateFirstMomentsWithInterpolation(double u, double v, double re, double im, double dataWt)
{
  //------------------------------------------------------------
  // Get the spatial frequency resolution of the axes
  //------------------------------------------------------------

  int nu = xAxis_.getNpix();
  int nv = yAxis_.getNpix();

  double du;
  double dv;
#ifdef EMLTEST
  static bool first = true;
#endif

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

  double s2    = 2*convSigInPixels_*convSigInPixels_;
  double wtSum = 0.0;

  //------------------------------------------------------------
  // Get the UV index closest to the current point
  //------------------------------------------------------------

  bool conjugateResult=false;
  unsigned uInd, vInd;

  uvIndex(u, v, uInd, vInd, conjugateResult);

#ifdef EMLTEST
  if(first) {
    int dftInd = (uInd * (nv/2+1) + vInd);
    COUT("(2) u = " << u << " v = " << v << " re = " << re << " im = " << im << " wt = " << dataWt 
	 << " uInd = " << uInd << " vInd = " << vInd << " dftInd = " <<  dftInd << " conj = " << conjugateResult
	 << " nu = " << nu << " nv = " << nv << " du = " << du << " dv = " << dv);
  }

  bool print = false;
  static int dftIndFid = -1;
  if(fabs(u + 912.253) < 0.1) {
    print = true;
    fprintf(stdout, "Found min u: uu = %f vv = %f upix = %d, vpix = %d\n", u, v, uInd, vInd);
  }
  if(fabs(u - 1262.086) < 0.1) {
    fprintf(stdout, "Found min u: uu = %f vv = %f upix = %d, vpix = %d\n", u, v, uInd, vInd);
  }
  if(fabs(v + 742.09) < 0.1) {
    fprintf(stdout, "Found min u: uu = %f vv = %f upix = %d, vpix = %d\n", u, v, uInd, vInd);
  }
  if(fabs(v - 1330.95) < 0.1) {
    fprintf(stdout, "Found min u: uu = %f vv = %f upix = %d, vpix = %d\n", u, v, uInd, vInd);
  }

#endif

  //------------------------------------------------------------
  // Get the UV coordinate of this nearest point
  //------------------------------------------------------------

  double uVal0, vVal0;
  uvCoord(uInd, vInd, uVal0, vVal0);

  // If the closest point was in the conjugate half of the array,
  // conjugate to make sure we do the separation calculation
  // correctly

  if(conjugateResult) {
    uVal0 = -uVal0;
    vVal0 = -vVal0;
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

#ifdef EMLTEST
  if(first) {
    COUT("uVal of closest point is: " << uVal0 << " v = " << vVal0);
  }
#endif

  //------------------------------------------------------------
  // Now iterate over a mask of pixels centered on the nearest pixel.
  //------------------------------------------------------------

  double reVal, imVal;
  double uVal, vVal;
  unsigned dftInd;
  for(int iU = -convMaskInPixels_; iU <= convMaskInPixels_; iU++) {
    for(int iV = -convMaskInPixels_; iV <= convMaskInPixels_; iV++) {

      //------------------------------------------------------------
      // Get the uv coordinate of this point
      //------------------------------------------------------------
  
      uVal = uVal0 + iU * du;
      vVal = vVal0 + iV * dv;

      //------------------------------------------------------------
      // Get the dft index corresponding to this point. 
      //------------------------------------------------------------

      try {
	dftIndex(uVal, vVal, dftInd, conjugateResult);
	reVal = re;
	imVal = conjugateResult ? -im : im;
      } catch(...) {

	//------------------------------------------------------------
	// If this point lies outside of our dft grid, just ignore it
	//------------------------------------------------------------

	continue;
      }

      //------------------------------------------------------------
      // Calculate the delta of the requested (u,v) point from the
      // center of this pixel, in fractional pixels
      //------------------------------------------------------------

      double duPix = (u - uVal)/du;
      double dvPix = (v - vVal)/dv;

#ifdef EMLTEST

      // Save the fiducial dft index

      if(print && dftIndFid < 0 && iU == 0 && iV == 0)
	dftIndFid = dftInd;

      if(first) {
	COUT("u = " << u << " v = " << v << " uVal = " << uVal << " vVal = " << vVal << " du = " << du << " dv = " << dv);
      }      
#endif

      //------------------------------------------------------------
      // And the weight corresponding to it
      //------------------------------------------------------------

      double arg = (duPix * duPix + dvPix * dvPix)/s2;

      //------------------------------------------------------------
      // The combined weight is the convolution weight multiplied by
      // the data weight
      //------------------------------------------------------------

      double wt = exp(-arg) * dataWt;

#if 0
#ifdef EMLTEST
      if(print || dftInd == dftIndFid) {
	fprintf(stdout, "Adding u = %f v = %f re = %f im=%f vis->wt = %f wttotal = %f to pixel iu = %d iv = %d dftInd = %d dftIndFid = %d wtSum = %f\n",
		u, v, reVal, imVal, dataWt, wt, iU, iV, dftInd, dftIndFid, wtSum_[dftInd]);
      }
#endif
#endif

      //------------------------------------------------------------
      // Finally, update the running mean for this pixel
      //------------------------------------------------------------

      double remean = out_[dftInd][0];
      double immean = out_[dftInd][1];

#if 1
      out_[dftInd][0] += ((reVal - remean) * wt / (wtSum_[dftInd] + wt));
      out_[dftInd][1] += ((imVal - immean) * wt / (wtSum_[dftInd] + wt));
#else
      out_[dftInd][0] += reVal * wt;
      out_[dftInd][1] += imVal * wt;
#endif

#if 0
#ifdef EMLTEST
      if(print || dftInd == dftIndFid) {
	COUT("reVal = " << reVal << " imVal = " << imVal << " wt = " << wt << " remean = " << remean << " immean = " << immean << " Re = " << out_[dftInd][0] << " Im = " << out_[dftInd][1] << " wtsum = " << wtSum_[dftInd] << " dftInd = " << dftInd);
      }
#endif
#endif

      //------------------------------------------------------------
      // Increment weight sums needed for running averages and errors
      //------------------------------------------------------------
      
      nPt_[dftInd]++;
      wtSum_[dftInd] += wt;
      wtSumTotal_    += wt;
    }
  }

#ifdef EMLTEST
  {
    if(first) {
      first = false;
      unsigned dftInd;
      dftIndex(u, v, dftInd, conjugateResult);
      COUT("Value of closest point is now: re = " << out_[dftInd][0] << " im = " << out_[dftInd][1] << " re = " << re << " im = " << im);    }
  }
#endif
}

/**.......................................................................
 * Add a visibility to the weighted running mean of the second moment
 * 
 * Here, by weight I mean weight, in the sense that the mean is given by:
 *
 *         sum(wt * (val - mean)^2)
 *   mn = ---------------
 *            sum(wt)
 *
 * So do NOT pass a variance instead, for example.
 */
void UvDataGridder::accumulateSecondMoments(double u, double v, double re, double im, double wt)
{
  unsigned dftInd;
  bool conj=false;

  // Get the index in our dft array corresponding to this uv point

  dftIndex(u, v, dftInd, conj);

  // Conjugate the data if required, to convert from (u, v) to (-u, -v)

  if(conj)
    im = -im;

  // Get the previously calculated mean for this index:

  double remean = out_[dftInd][0];
  double immean = out_[dftInd][1];

  // Construct the current value of (val- mean)^2

  double re2 = (re - remean) * (re - remean);
  double im2 = (im - immean) * (im - immean);

  // And the running average of the second moment so far
  
  double re2mean =  errorInMean_[dftInd][0];
  double im2mean =  errorInMean_[dftInd][1];

  errorInMean_[dftInd][0] += (re2 - re2mean) * wt / (wtSum_[dftInd] + wt);
  errorInMean_[dftInd][1] += (im2 - im2mean) * wt / (wtSum_[dftInd] + wt);

  // Increment weight sums needed for running averages and errors

  nPt_[dftInd]++;
  wtSum_[dftInd]  += wt;
  wt2Sum_[dftInd] += wt*wt;
  wtSumTotal_     += wt;
}

/**.......................................................................
 * Add a visibility to the weighted running mean of the second moment
 * 
 * Here, by weight I mean weight, in the sense that the mean is given by:
 *
 *         sum(wt * (val - mean)^2)
 *   mn = ---------------
 *            sum(wt)
 *
 * So do NOT pass a variance instead, for example.
 */
void UvDataGridder::accumulateSecondMomentsWithInterpolation(double u, double v, double re, double im, double dataWt)
{
  //------------------------------------------------------------
  // Get the spatial frequency resolution of the axes
  //------------------------------------------------------------

  int nu = xAxis_.getNpix();
  int nv = yAxis_.getNpix();

  double du;
  double dv;
#ifdef EMLTEST
  static bool first = true;
#endif

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

  double s2    = 2*convSigInPixels_*convSigInPixels_;
  double wtSum = 0.0;

  //------------------------------------------------------------
  // Get the UV index closest to the current point
  //------------------------------------------------------------

  bool conjugateResult=false;
  unsigned uInd, vInd;

  uvIndex(u, v, uInd, vInd, conjugateResult);

#ifdef EMLTEST
  if(first) {
    int dftInd = (uInd * (nv/2+1) + vInd);
    COUT("(2) u = " << u << " v = " << v << " re = " << re << " im = " << im << " wt = " << dataWt 
	 << " uInd = " << uInd << " vInd = " << vInd << " dftInd = " <<  dftInd << " conj = " << conjugateResult);
  }
#endif

  //------------------------------------------------------------
  // Get the UV coordinate of this nearest point
  //------------------------------------------------------------

  double uVal0, vVal0;
  uvCoord(uInd, vInd, uVal0, vVal0);

  // If the closest point was in the conjugate half of the array,
  // conjugate to make sure we do the separation calculation
  // correectly

  if(conjugateResult) {
    uVal0 = -uVal0;
    vVal0 = -vVal0;
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
  // Now iterate over a mask of pixels centered on the nearest pixel.
  //------------------------------------------------------------

  double reVal, imVal;
  double uVal, vVal;
  unsigned dftInd;

  for(int iU = -convMaskInPixels_; iU <= convMaskInPixels_; iU++) {
    for(int iV = -convMaskInPixels_; iV <= convMaskInPixels_; iV++) {

      //------------------------------------------------------------
      // Get the uv coordinate of this point
      //------------------------------------------------------------
  
      uVal = uVal0 + iU * du;
      vVal = vVal0 + iV * dv;

      //------------------------------------------------------------
      // Get the dft index corresponding to this point. 
      //------------------------------------------------------------

      try {
	dftIndex(uVal, vVal, dftInd, conjugateResult);
	reVal = re;
	imVal = conjugateResult ? -im : im;
      } catch(...) {

	//------------------------------------------------------------
	// If this point lies outside of our dft grid, just ignore it
	//------------------------------------------------------------

	continue;
      }

      //------------------------------------------------------------
      // Calculate the delta of the requested (u,v) point from the
      // center of this pixel, in fractional pixels
      //------------------------------------------------------------

      double duPix = (u - uVal)/du;
      double dvPix = (v - vVal)/dv;

#ifdef EMLTEST
      if(first) {
	COUT("u = " << u << " v = " << v << " uVal = " << uVal << " vVal = " << vVal << " du = " << du << " dv = " << dv);
      }      
#endif

      //------------------------------------------------------------
      // And the weight corresponding to it
      //------------------------------------------------------------

      double arg = (duPix * duPix + dvPix * dvPix)/s2;

      //------------------------------------------------------------
      // The combined weight is the convolution weight multiplied by
      // the data weight
      //------------------------------------------------------------

      double wt = exp(-arg) * dataWt;

#ifdef EMLTEST
      if(first) {
	COUT("arg = " << arg << " dataWt = " << dataWt << " wt = " << wt);
      }      
#endif

      //------------------------------------------------------------
      // Finally, update the running mean for this pixel
      //------------------------------------------------------------

      double remean = out_[dftInd][0];
      double immean = out_[dftInd][1];

      // Construct the current value of (val- mean)^2

      double re2 = (reVal - remean) * (reVal - remean);
      double im2 = (imVal - immean) * (imVal - immean);

      // And the running average of the second moment so far
      
      double re2mean =  errorInMean_[dftInd][0];
      double im2mean =  errorInMean_[dftInd][1];
      
      errorInMean_[dftInd][0] += (re2 - re2mean) * wt / (wtSum_[dftInd] + wt);
      errorInMean_[dftInd][1] += (im2 - im2mean) * wt / (wtSum_[dftInd] + wt);
      
      // Increment weight sums needed for running averages and errors
      
      nPt_[dftInd]++;
      wtSum_[dftInd]  += wt;
      wt2Sum_[dftInd] += wt*wt;
      wtSumTotal_     += wt;
    }
  }
}

/**.......................................................................
 * Convert from moment sums to error in the mean.  Should only be
 * called after all data have been gridded into this object.
 */
void UvDataGridder::calculateErrorInMean()
{
  double prefac;
  double refm, refm2, resm;
  double imfm, imfm2, imsm;

  static bool count=0;

  for(unsigned i=0; i < nOutZeroPad_; i++) {
    
    if(nPt_[i] > 0) {

      //------------------------------------------------------------
      // Mark this index as containing data
      //------------------------------------------------------------

      populatedIndices_.push_back(i);

      //------------------------------------------------------------
      // Store the UV coordinate of this point for fast access later
      //------------------------------------------------------------

      double u, v, val;
      getUVData(i, DATA_UV, u, v, val);

      populatedU_.push_back(u);
      populatedV_.push_back(v);

      //------------------------------------------------------------
      // If estimating errors from the variance of the data, get the
      // error in the mean from the appropriate rescaling by weight
      // sums
      //------------------------------------------------------------

      if(estimateErrInMeanFromData_) {

	//------------------------------------------------------------
	// Calculate the prefactor to the variance (weighted equivalent to 1/(N-1)
	//------------------------------------------------------------

	prefac = wt2Sum_[i] / (wtSum_[i] * wtSum_[i] - wt2Sum_[i]);
	
	//------------------------------------------------------------
	// Now overwrite the second moments stored in the errorInMean_
	// array with the true error in the mean.  If only one point
	// made it into this bin, set the error to the sqrt of the
	// inverse weight (assumed to be an estimate of the variance).
	//------------------------------------------------------------

	if(nPt_[i] == 1) {
	  double sigma = sqrt(1.0/wtSum_[i]);
	  errorInMean_[i][0] = sigma;
	  errorInMean_[i][1] = sigma;
	} else {
	  errorInMean_[i][0] = sqrt(prefac * errorInMean_[i][0]);
	  errorInMean_[i][1] = sqrt(prefac * errorInMean_[i][1]);
	}

	//------------------------------------------------------------
	// Else take the error in the mean directly from the weights.
	//------------------------------------------------------------

      } else {
	double sigma = sqrt(1.0/(wtSum_[i]));
	errorInMean_[i][0] = sigma;
	errorInMean_[i][1] = sigma;
      }

    }

  }

  errorInMeanIsValid_ = true;
}

void UvDataGridder::estimateErrorInMeanFromData(bool estimate)
{
  estimateErrInMeanFromData_ = estimate;
}

/**.......................................................................
 * Return an array of populated indices
 */
std::vector<double> UvDataGridder::getPopulatedIndices()
{
  std::vector<double> arr(nOutZeroPad_);
  double u, v;
  
  unsigned dftInd, imInd;
  for(unsigned ix=0; ix < nxZeroPad_; ix++) {
    for(unsigned iy=0; iy < nyZeroPad_/2+1; iy++) {
      dftInd = ix * (nyZeroPad_/2+1) + iy;
      imInd  = iy * nxZeroPad_   + ix;
      arr[imInd] = nPt_[dftInd] > 0 ? 1.0 : 0.0;
    }
  }

  return arr;
}

void UvDataGridder::zero()
{
  unsigned nInd = populatedIndices_.size();
  unsigned dftInd;

  for(unsigned i=0; i < nInd; i++) {
    dftInd = populatedIndices_[i];
    out_[dftInd][0] = 0.0;
    out_[dftInd][1] = 0.0;
  }
}

void UvDataGridder::assignPopulatedIndicesFrom(UvDataGridder& gridder)
{
  populatedIndices_ = gridder.populatedIndices_;
  populatedU_       = gridder.populatedU_;
  populatedV_       = gridder.populatedV_;
}

/**.......................................................................
 * Return the estimated rms of a map made from inverting the UV data
 * in this container
 */
double UvDataGridder::estimateMapPixelRms()
{
  if(!errorInMeanIsValid_) {
    ThrowColorError("No errors have been calculated yet -- use calculateErrorInMean() first", "red");
  }

  double wtsum = 0.0;
  double sigma, wt;
  unsigned ind;

  for(unsigned i=0; i < populatedIndices_.size(); i++) {
    ind = populatedIndices_[i];
    sigma = errorInMean_[ind][0];
    wt = 1.0/(sigma * sigma);
    wtsum += wt;
  }

  return sqrt(1.0/wtsum);
}

void UvDataGridder::initializePopulatedIndicesToAll()
{
  populatedIndices_.resize(0);
  populatedU_.resize(0);
  populatedV_.resize(0);

  for(unsigned i=0; i < nOutZeroPad_; i++) {
    double u, v, val;
    getUVData(i, DATA_UV, u, v, val);
    populatedIndices_.push_back(i);
    populatedU_.push_back(u);
    populatedV_.push_back(v);
  }
}

/**.......................................................................
 * Make a plot of the occupied UV cells in this gridder
 */
void UvDataGridder::plotOccupiedUV()
{
  std::vector<double> x, y;

  x.resize(populatedIndices_.size());
  y.resize(populatedIndices_.size());

  for(unsigned i=0; i < populatedIndices_.size(); i++) {
    unsigned dftInd = populatedIndices_[i];
    double u,v,val;
    uvCoord(dftInd, u, v);
    x[i] = u;
    y[i] = v;
  }

  PgUtil::linePlot(x.size(), &x[0], &y[0], 0, "U", "V", "", false);
}

/**.......................................................................
 * Overload this from the base class. While the Dft2d version just
 * does a straight inverse transform, this overloaded method reflects
 * the dual-purpose nature of the UvDataGridder object.
 *
 * Its main use is to store a gridded array of weighted means in the
 * normal out_ array, which we store along with their weights in the
 * wtSum_ array.  However, if we want to transform this data to make
 * an image, we need to re-weight the data accordingly.
 */
void UvDataGridder::computeInverseTransform(fftw_plan* invPlan)
{
  unsigned nInd = populatedIndices_.size();
  unsigned dftInd;

  //------------------------------------------------------------
  // Renormalize with weights appropriate for making maps
  //------------------------------------------------------------

  bool first = true;
  for(unsigned i=0; i < nInd; i++) {
    dftInd = populatedIndices_[i];

    if(debugPrint_) {
      if(first) {
	COUT("dftInd = " << dftInd << " re = " << out_[dftInd][0] << " im = " << out_[dftInd][1] 
	     << " wt = " << wtSum_[dftInd] << " total = " << wtSumTotal_);
	first = false;
      }
    }

    out_[dftInd][0] *= (wtSum_[dftInd] / (2*wtSumTotal_));
    out_[dftInd][1] *= (wtSum_[dftInd] / (2*wtSumTotal_));
  }

  Dft2d::computeInverseTransform(invPlan);
}

void UvDataGridder::renormalize()
{
  unsigned nInd = populatedIndices_.size();
  unsigned dftInd;
  bool first = true;

  for(unsigned i=0; i < nInd; i++) {
    dftInd = populatedIndices_[i];
    out_[dftInd][0] *= (wtSum_[dftInd] / (2*wtSumTotal_));
    out_[dftInd][1] *= (wtSum_[dftInd] / (2*wtSumTotal_));
  }
}

void UvDataGridder::computeInverseTransformNoRenorm(fftw_plan* invPlan)
{
  Dft2d::computeInverseTransform(invPlan);
}

void UvDataGridder::operator*=(double mult)
{
  if(!hasData_) {
    ThrowError("Dft contains no data");
  }

  for(unsigned i=0; i < populatedIndices_.size(); i++) {
    unsigned dftInd = populatedIndices_[i];
    out_[dftInd][0] *= mult;
    out_[dftInd][1] *= mult;
  }
}

void UvDataGridder::copyToStore()
{
  for(unsigned i=0; i < nOutZeroPad_; i++) {
    store_[i][0] = out_[i][0];
    store_[i][1] = out_[i][1];
  }
}

void UvDataGridder::copyFromStore()
{
  for(unsigned i=0; i < nOutZeroPad_; i++) {
    out_[i][0] = store_[i][0];
    out_[i][1] = store_[i][1];
  }
}


void UvDataGridder::debugPrint(bool print)
{
  debugPrint_ = print;
}

void UvDataGridder::getEstimatedSynthesizedBeam(Angle& fwhmMaj, Angle& fwhmMin, Angle& posAngle)
{
  const float fudge=0.7f; // Empirical fudge factor of TJP's algorithm

  double uu = uuSum_;
  double vv = vvSum_;
  double uv = uvSum_;

  float ftmp = sqrt((uu-vv)*(uu-vv) + 4.0 * uv * uv);
  
  // First the position angle of the equivalent elliptical gaussian
  // distribution.

  posAngle.setRadians(-0.5*atan2(2.0*uv, uu - vv));
  posAngle.setRadians(2*M_PI - posAngle.radians());
  
  double radVal = posAngle.radians();

  Angle mod(Angle::Radians(), posAngle.radians(), true);
  posAngle.setRadians(mod.radians());

  // Then the equivalent elliptical beam widths in radians.

  fwhmMin.setRadians(fudge/(sqrt(2.0*(uu+vv) + 2.0*ftmp)));
  fwhmMaj.setRadians(fudge/(sqrt(2.0*(uu+vv) - 2.0*ftmp)));

  return;
}

/**.......................................................................
 * Co-add two sets of gridded data.  Note that this is NOT the same as
 * operator+=(), which simply adds two grids together, used to
 * accumulate Fourier-space models, for example.  This constructs the
 * statistically correct average of two visibility datasets, keeping
 * track of the running average and the total weight sum.
 */
void UvDataGridder::mergeData(UvDataGridder& gridder)
{
  if(estimateErrInMeanFromData_)
    ThrowSimpleColorError("Unable to estimate the errors from the data when combining gridders", "red");

  unsigned nInd = gridder.populatedIndices_.size();
  unsigned dftInd;

  double u, v, re, im, wt, val;

  //------------------------------------------------------------
  // Iterate over all populated indices in the passed gridder,
  // accumulating its data into ours
  //------------------------------------------------------------

  bool first = true;
  for(unsigned i=0; i < nInd; i++) {
    dftInd = gridder.populatedIndices_[i];
    gridder.getUVData(dftInd, DATA_UV, u, v, val);

    re = gridder.out_[dftInd][0];
    im = gridder.out_[dftInd][1];
    wt = gridder.wtSum_[dftInd]; 

    if(first)
      debugPrint(true);
    else
      debugPrint(false);

#if 1
    accumulateFirstMoments(u, v, re, im, wt);
#else
    accumulateFirstMomentsTest(u, v, re, im, wt);
#endif

    if(first) {
      COUT("Accumulating with wt = " << wt);
      first = false;
    }
  }

  //------------------------------------------------------------
  // Now recalculate our errors
  //------------------------------------------------------------

  calculateErrorInMean();
}

void UvDataGridder::plotWt()
{
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

  PgUtil::setTitle(true);
  PgUtil::greyScale(wtSum_.size(), &wtSum_[0], nxZeroPad_, nyZeroPad_/2+1, xmin, xmax, ymin, ymax, 0, 0, 0, (char*)xlabel.c_str(), (char*)ylabel.c_str(), (char*)title.c_str());
}
