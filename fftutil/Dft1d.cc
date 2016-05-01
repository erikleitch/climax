#include "gcp/fftutil/Dft1d.h"

#include "gcp/util/Exception.h"

#include<iostream>
#include<cmath>

using namespace std;

using namespace gcp::util;

#define EPS 1e-12

/**.......................................................................
 * Constructor.
 */
Dft1d::Dft1d(int n, bool optimize, Apodization apod)
{
  in_           = 0;
  out_          = 0;
  nAv_          = 0;
  haveRes_      = false;
  mean_         = 0.0;
  optimize_     = optimize;
  outAv_        = 0;
  doAv_         = false;
  subtractMean_ = true;

  resize(n);

  setAverage(false);
  setVectorAverage(true);

  assertApodizationType(apod);
}

/**.......................................................................
 * Const Copy Constructor.
 */
Dft1d::Dft1d(const Dft1d& objToBeCopied)
{
  *this = (Dft1d&)objToBeCopied;
};

/**.......................................................................
 * Copy Constructor.
 */
Dft1d::Dft1d(Dft1d& objToBeCopied)
{
  *this = objToBeCopied;
};

/**.......................................................................
 * Const Assignment Operator.
 */
void Dft1d::operator=(const Dft1d& objToBeAssigned)
{
  *this = (Dft1d&)objToBeAssigned;
};

/**.......................................................................
 * Assignment Operator.
 */
void Dft1d::operator=(Dft1d& objToBeAssigned)
{
  std::cout << "Calling default assignment operator for class: Dft1d" << std::endl;
};

/**.......................................................................
 * Output Operator.
 */
std::ostream& gcp::util::operator<<(std::ostream& os, Dft1d& obj)
{
  os << "Default output operator for class: Dft1d" << std::endl;
  return os;
};

/**.......................................................................
 * Destructor.
 */
Dft1d::~Dft1d() 
{
  if(in_ != 0) {
    fftw_free(in_);
    in_ = 0;
  }

  if(out_ != 0) {
    fftw_free(out_);
    out_ = 0;
  }
}

/**.......................................................................
 * Compute a plan for this fft
 */
void Dft1d::computePlan(unsigned flag) 
{
  plan_ = fftw_plan_dft_r2c_1d(n_, in_, out_, flag);
}

/**.......................................................................
 * Resize for FFTs of a different length
 */
void Dft1d::resize(unsigned n)
{
  n_ = n;
      
  if(in_) {
    fftw_free(in_);
    in_ = 0;
  }

  if(out_) {
    fftw_free(out_);
    out_ = 0;
  }

  if(outAv_) {
    fftw_free(outAv_);
    outAv_ = 0;
  }

  // Allocate arrays

  if((in_ = (double*)fftw_malloc(n_ * sizeof(double)))==0)
    ThrowError("Couldn't allocate input data array");

  if((out_ = (fftw_complex*)fftw_malloc(n_ * sizeof(fftw_complex)))==0)
    ThrowError("Couldn't allocate output data array");

  if((outAv_ = (fftw_complex*)fftw_malloc(n_ * sizeof(fftw_complex)))==0)
    ThrowError("Couldn't allocate averaged output data array");

  // And resize the buffer into which we average power spectra

  absBuf_.resize(transformSize());
  
  // Resize the buffers we will use for storing samples

  sampBuf_.resize(n);
  apodBuf_.resize(n);

  // Now compute a plan for this machine architecture
      
  computePlan(optimize_ ? FFTW_MEASURE : FFTW_ESTIMATE);

  // And reset pertinent members

  transformIsReady_    = false;
  nSinceLastTransform_ = 0;
  nAv_                 = 0;
}

/**.......................................................................
 * Push a sample onto the input array
 */
void Dft1d::pushSample(double sample) 
{
  // Push the sample and the corresponding apodization coefficient
  // into the buffers

  sampBuf_[nSinceLastTransform_] = sample;
  apodBuf_[nSinceLastTransform_] =  apodizationCoefficient(nSinceLastTransform_);

  // Keep a running mean of samples

  mean_ += (sample - mean_)/(nSinceLastTransform_ + 1);

  // Check if we should compute a transform

  ++nSinceLastTransform_;
  transformIsReady_ = false;

  if(nSinceLastTransform_ == n_) {

    copySamples();
    computeTransform();

    nSinceLastTransform_ = 0;
    mean_                = 0.0;

    // If a request for a new apodization type is pending, assert it
    // now

    if(apodTypePending_) {
      assertApodizationType(pendingApodType_);
    }

  }
}

double Dft1d::apodizationCoefficient(unsigned iSamp)
{
  return apodFn_(iSamp, n_);
}

/**.......................................................................
 * Fill the input data array
 */
void Dft1d::fillInputArray(double* data) 
{
  double* from = data;
  double* to   = &in_[0];

  for(unsigned i=0; i < n_; i++)
    *(to+i) = *(from+i);
}

/**.......................................................................
 * Compute the transform
 */
void Dft1d::removeMean()
{
  double mean=0.0;

  for(unsigned i=0; i < n_; i++) {
    mean += (in_[i] - mean)/(i + 1);
  }

  for(unsigned i=0; i < n_; i++) {
    in_[i] -= mean;
  }
}

/**.......................................................................
 * Compute the transform
 */
void Dft1d::computeTransform()
{
  removeMean();
  fftw_execute(plan_);

  if(doAv_) {
    addToAverage();
  }

  transformIsReady_ = true;
}

/**.......................................................................
 * Add a transform to the averaege
 */
void Dft1d::addToAverage()
{
  double re, im, absVal;

  for(unsigned i=0; i < transformSize(); i++) {

    re = out_[i][0];
    im = out_[i][1];

    // If vector-averaging, take the running mean of the re and im

    if(vecAv_) { 
      outAv_[i][0] += (re - outAv_[i][0])/(nAv_ + 1);
      outAv_[i][1] += (im - outAv_[i][1])/(nAv_ + 1);

      // Else just average the variance

    } else {
      absVal = re*re + im*im;
      absBuf_[i] += (absVal - absBuf_[i])/(nAv_ + 1);
    }

  }

  ++nAv_;
}

/**.......................................................................
 * Return a pointer to the input array
 */
double* Dft1d::getInputData()
{
  return in_;
}

/**.......................................................................
 * Return a pointer to the transformed data
 */
fftw_complex* Dft1d::getTransform()
{
  return out_;
}

/**.......................................................................
 * Return the size of the transformed array
 */
unsigned Dft1d::inputSize()
{
  return n_;
}

/**.......................................................................
 * Return the size of the transformed array
 */
unsigned Dft1d::transformSize()
{
  return n_/2 + 1;
}

/**.......................................................................
 * Return the size of the transformed array
 */
unsigned Dft1d::transformSize(unsigned n)
{
  return n/2 + 1;
}

/**.......................................................................
 * Return the absolute value of the output array
 */
std::vector<double> Dft1d::abs()
{
  std::vector<double> output;

  unsigned nOut = transformSize();

  output.resize(nOut);

  double* outputPtr = &output[0];

  if(!doAv_) {

    for(unsigned i=0; i < nOut; i++)
      outputPtr[i] = sqrt(out_[i][0] * out_[i][0] + 
			  out_[i][1] * out_[i][1]);

  } else {

    if(vecAv_) {
      for(unsigned i=0; i < nOut; i++) 
	outputPtr[i] = sqrt(outAv_[i][0] * outAv_[i][0] + 
			    outAv_[i][1] * outAv_[i][1]);
    } else {

      for(unsigned i=0; i < nOut; i++)
	outputPtr[i] = sqrt(absBuf_[i]); 

    }
  }

  return output;
}

/**.......................................................................
 * Return the powerd
 */
std::vector<double> Dft1d::abs2()
{
  std::vector<double> output;

  unsigned nOut = transformSize();

  output.resize(nOut);

  double* outputPtr = &output[0];

  if(!doAv_) {

    for(unsigned i=0; i < nOut; i++)
      outputPtr[i] = out_[i][0] * out_[i][0] + 
	out_[i][1] * out_[i][1];

  } else {

    if(vecAv_) {
      for(unsigned i=0; i < nOut; i++) 
	outputPtr[i] = outAv_[i][0] * outAv_[i][0] + 
	  outAv_[i][1] * outAv_[i][1];
    } else {

      for(unsigned i=0; i < nOut; i++)
	outputPtr[i] = absBuf_[i]; 

    }
  }

  return output;
}

/**.......................................................................
 * Update the min max
 */
void Dft1d::updateMinMax(unsigned i, double& val, bool doLog, 
		       bool& first, double& min, double& max)
{
  // If the log of the data is being returned, check that val > 0.0.
  // If it's not, set it to some tiny number and don't include this
  // point in the min/max calculation.  Also, if doLog = true, we
  // ignore this value if it is the zero-frequency component, since
  // this will be very close to zero.

  if(!doLog || (i > 0 && val > 0.0)) {
    if(first) {
      min = max = val;
      first = false;
    } else {
      min = val < min ? val : min;
      max = val > max ? val : max;
    }
  }

  if(doLog) {
    if(val > 0.0) {
      val = log10(val);
    } else {
      val = EPS;
      return;
    }
  }

}

/**.......................................................................
 * Return the absolute value of the output array
 */
void Dft1d::abs(float* outputPtr, bool& first, 
	      double& min, double& max, bool doLog)
{
  unsigned nOut = transformSize();

  double val;

  // If not averaging, just return sqrt(transform) in the output array

  if(!doAv_) {

    for(unsigned i=0; i < nOut; i++) {
      val = sqrt(out_[i][0] * out_[i][0] + out_[i][1] * out_[i][1]);
      updateMinMax(i, val, doLog, first, min, max);
      outputPtr[i] = val;
    }

    // Else if averaging

  } else {

    // If vector-averaging

    if(vecAv_) {

      for(unsigned i=0; i < nOut; i++) {
	val = sqrt(outAv_[i][0] * outAv_[i][0] + outAv_[i][1] * outAv_[i][1]);
	updateMinMax(i, val, doLog, first, min, max);
	outputPtr[i] = val;
      }

    // If rms-averaging

    } else {

      for(unsigned i=0; i < nOut; i++) {
	val = sqrt(absBuf_[i]); 
	updateMinMax(i, val, doLog, first, min, max);
	outputPtr[i] = val;
      }

    }
  }
}

/**.......................................................................
 * Return the power spectrum of the input array
 */
std::vector<double> Dft1d::powerSpectrum(double* data)
{
  fillInputArray(data);
  computeTransform();
  return abs();
}

/**.......................................................................
 * Return true if the input array is full
 */
bool Dft1d::transformIsReady()
{
  return transformIsReady_;
}

/**.......................................................................
 * If true, we will store on-the-fly averages of the transforms
 */
void Dft1d::setAverage(bool doAv)
{
  nAv_   = 0;
  doAv_  = doAv;

  for(unsigned i=0; i < n_; i++) {
    outAv_[i][0] = 0.0;
    outAv_[i][1] = 0.0;
  }

  for(unsigned i=0; i < transformSize(); i++) 
    absBuf_[i] = 0.0;
}

/**.......................................................................
 * If true, we will store on-the-fly averages of the transforms
 */
void Dft1d::setVectorAverage(bool vectorAverage)
{
  // If we are currently averaging, and this call changes the
  // integration type, rezero the average counter

  if(doAv_ && vecAv_ != vectorAverage)
    nAv_ = 0;

  // And set the averaging type

  vecAv_ = vectorAverage;
}

/**.......................................................................
 * Set the time resolution of the x-axis.  This will be used to
 * calculate frequency span and resolution of the transform.
 */
void Dft1d::setTimeRes(TimeVal tVal)
{
  tRes_ = tVal;
  haveRes_ = true;
}

/**.......................................................................
 * Get the frequency resolution of the transform
 */
Frequency Dft1d::getFrequencyResolution()
{
  Frequency freq;

  if(haveRes_) {
    freq.setHz(1.0/(n_ * tRes_.getTimeInSeconds()));
  } else {
    ThrowError("No time resolution has been specified");
  }

  return freq;
}

/**.......................................................................
 * Get the frequency resolution of the transform
 */
Frequency Dft1d::getFrequencyResolution(unsigned n, TimeVal timeResolution)
{
  Frequency freq;

  freq.setHz(1.0/(n * timeResolution.getTimeInSeconds()));

  return freq;
}

/**.......................................................................
 * User-level method to request an apodization type.  This will be
 * asserted at the start of the next transform
 */
void Dft1d::setApodizationType(Apodization type)
{
  pendingApodType_ = type;
  apodTypePending_ = true;
}

void Dft1d::assertApodizationType(Apodization type)
{
  switch (type) {
  case APOD_RECTANGLE:
    apodFn_ = apodRectangle;
    break;
  case APOD_TRIANGLE:
    apodFn_ = apodTriangle;
    break;
  case APOD_HAMMING:
    apodFn_ = apodHamming;
    break;
  case APOD_HANN:
    apodFn_ = apodHann;
    break;
  case APOD_COS:
    apodFn_ = apodCos;
    break;
  case APOD_SINC:
    apodFn_ = apodSinc;
    break;
  default:
    apodFn_ = apodRectangle;
    break;
  }

  currentApodType_ = type;
  pendingApodType_ = type;
  apodTypePending_ = false;
}

//-----------------------------------------------------------------------
// Apodization functions
//-----------------------------------------------------------------------

APOD_FN(Dft1d::apodRectangle)
{
  return 1.0;
}

APOD_FN(Dft1d::apodHamming)
{
  return 0.54 - 0.46 * cos((2*M_PI*iSamp)/(n - 1));
}

APOD_FN(Dft1d::apodHann)
{
  return 0.5 * (1.0 - cos((2*M_PI*iSamp)/(n - 1)));
}

APOD_FN(Dft1d::apodCos)
{
  return sin((M_PI * iSamp)/(n - 1));
}

APOD_FN(Dft1d::apodTriangle)
{
  // For now just default to triangle apodization

  double iVal    = (double)iSamp;
  double halfVal = (double)(n)/2;

  if(iVal < halfVal) {
    return 2.0/(n - 1) * iVal;
  } else {
    return 2.0/(n - 1) * (n - 1 - iSamp);
  }
}

APOD_FN(Dft1d::apodSinc)
{
  if(2*iSamp == n-1) {
    return 1.0;
    } else {
    double darg = (2.0 * iSamp) / (n - 1) - 1.0;
    return sin(darg) / darg;
  }
}

/**.......................................................................
 * Move samples from the buffer to the FFTW input array
 */
void Dft1d::copySamples()
{
  if(subtractMean_) {
    sampBuf_ = (sampBuf_ - mean_) * apodBuf_;
  } else {
    sampBuf_ = sampBuf_ * apodBuf_;
  }
  
  fillInputArray(&sampBuf_[0]);
}
