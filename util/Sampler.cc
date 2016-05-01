#include "gcp/util/Exception.h"
#include "gcp/util/Sampler.h"

#include <cmath>

using namespace std;
using namespace gcp::util;

static int gcf(double *gammcf, double a, double x, double *gln);
static int gser(double *gamser, double a, double x, double *gln);

#define EPSFRAC 1e-12

const double Sampler::posInf_ =  1.0/0.0;
const double Sampler::negInf_ = -1.0/0.0;
const double Sampler::sqrt2pi_= sqrt(2*M_PI);

/**.......................................................................
 * Constructor.
 */
Sampler::Sampler() 
{
  haveFn_ = false;
  isFn_   = false;
}

/**.......................................................................
 * Destructor.
 */
Sampler::~Sampler() {}

void Sampler::setYX(SAMPLER_FN(fn))
{
}

/**.......................................................................
 * Set the function we will use to generate samples
 */
void Sampler::setdYdX(SAMPLER_FN(*fn), double xmin, double xmax, double dx, double* args)
{
  unsigned nSamp = (unsigned)((xmax - xmin)/dx) + 1;

  std::vector<double> x, y;

  x.resize(nSamp);
  y.resize(nSamp);

  for(unsigned i = 0; i < nSamp; i++) {
    x[i] = xmin + i*dx;
    y[i] = fn(x[i], args);
  }

  setdYdX(x, y);
}

/**.......................................................................
 * Set the function we will use to generate samples
 */
void Sampler::setdYdX(std::vector<double>& x, std::vector<double>& y)
{
  setdYdX(x.size(), &x[0], &y[0]);
}

/**.......................................................................
 * Set the function we will use to generate samples
 */
void Sampler::setdYdX(unsigned n, double* x, double* y)
{
  haveFn_ = true;
  isFn_   = false;

  nPt_ = n+1;

  // Resize the differential and integral versions

  x_.resize(nPt_);
  y_.resize(nPt_);

  xInt_.resize(nPt_);
  yInt_.resize(nPt_);

  // Integrate the function now

  double norm = 0.0;
  double dx = x[1]-x[0];

  for(unsigned i=0; i < n; i++) {
    x_[i] = x[i];
    y_[i] = y[i];
    norm += y[i];
  }

  xInt_[0] = x[0]-dx/2;
  for(unsigned i=1; i <= n; i++) {
    xInt_[i] = x[i-1]+dx/2;
  }

  yInt_[0] = 0.0;
  for(unsigned i=1; i <= n; i++) {
    yInt_[i] = yInt_[i-1] + y[i-1]/norm;
  }

#if 0
  for(unsigned i=1; i < n; i++) {
    COUT("xint[i] = " << xInt_[i] << " yint[i] = " << yInt_[i]);
  }
#endif

}

double Sampler::generateSample()
{
  if(nPt_ == 0)
    ThrowError("No sampling function has been specified");

  double sample;

  if(!isFn_) {
    sample = binSearchForSample();
  } else {
    ThrowError("No functional form is currently allowed");
  }

  return sample;
}

/**.......................................................................
 * Generate samples according to the sampling function we were passed
 */
std::vector<double> Sampler::generateSamples(unsigned nSamp)
{
  if(nPt_ == 0)
    ThrowError("No sampling function has been specified");

  std::vector<double> samples(nSamp);

  if(!isFn_) {
    for(unsigned i=0; i < nSamp; i++) {
      samples[i] = binSearchForSample();
    }
  } else {
    
  }

  return samples;
}

/**.......................................................................
 * Binary search for the closest sample in the y-array.
 */
double Sampler::binSearchForSample()
{
  double y=(double)(rand())/RAND_MAX;

  double yHi, yLo, yMid;
  unsigned lo=0, hi=nPt_-1, mid;

  while(hi - lo > 1) {

    mid = (hi+lo)/2;

    yHi  = yInt_[hi];
    yLo  = yInt_[lo];
    yMid = yInt_[mid];

    if(y > yMid) {
      lo = mid;
    } else if(y < yMid) {
      hi = mid;
    } else {
      lo = hi = mid;
    }
  }

  if(lo == hi)
    return yInt_[hi];
  else
    return xInt_[lo] + (xInt_[hi] - xInt_[lo])/(yInt_[hi] - yInt_[lo]) * 
      (y-yInt_[lo]);
}

void Sampler::seed(unsigned int s)
{
  srand(s);
}

void Sampler::seedRandom()
{
  int r = rand();
  srand((unsigned int)r);
}

/**.......................................................................
 * Generate poisson samples
 */
unsigned Sampler::
generatePoissonSample(double mean)
{
  std::vector<unsigned> samps = generatePoissonSamples(mean, 1);
  return samps[0];
}

/**.......................................................................
 * Generate poisson samples
 */
std::vector<unsigned> Sampler::
generatePoissonSamples(double mean, unsigned ndev)
{
  int go=1,waserr=0,i;
  double y,ysmall,ymid,ylarge;
  
  std::vector<unsigned> outvals(ndev);

  // Generate the deviates here.
  // Test by doing a binary search.
  
  for(i=0;i < ndev && !waserr;i++) {
    
    double nlarge=10*mean,nsmall=1;
    double nmid=(nlarge+nsmall)/2,nrange=(nlarge-nsmall)/2;
    double nclose;
    
    y = (double)rand()/RAND_MAX;
    
    // First find values of nsmall and nlarge which bracket the value.

    do {
      waserr |= gammp(nsmall, mean, &ysmall);
      waserr |= gammp(nlarge, mean, &ylarge);
      
      ysmall = 1.0-ysmall;
      ylarge = 1.0-ylarge;
      
      if(!waserr) {
	
	// If y is above the range, increase the range until it is
	// bracketed.  Note that y can't be below the range, since
	// gammp(0,mean) = 0.0

	if(y > ylarge) {
	  nsmall = nlarge;
	  nlarge += nrange;
	} else {
	  go = 0;
	}
	
	nrange = (nlarge-nsmall)/2;
	nmid = (nlarge+nsmall)/2;
	
      }
    } while(go && !waserr);
    
    // Now we have bracketing values.  DO a binary search to narrow
    // down the range.

    while(nrange > 0.1 && !waserr) {
      
      waserr |= gammp(nsmall, mean, &ysmall);
      waserr |= gammp(nmid,   mean, &ymid);
      waserr |= gammp(nlarge, mean, &ylarge);
      
      ysmall = 1.0-ysmall;
      ymid = 1.0-ymid;
      ylarge = 1.0-ylarge;
      
      if(!waserr) {
	if(y > ymid)
	  nsmall = nmid;
	else if(y < ymid)
	  nlarge = nmid;
	else
	  nsmall = nlarge = nmid;
	
	nrange = (nlarge-nsmall);
	nmid = (nlarge+nsmall)/2;
	nclose = fabs(y-ysmall) < fabs(y-ylarge) ? nsmall : nlarge;
      }
    }
    
    // When we exit the above loop, nmid should be the value we are
    // looking for.

    if(!waserr)
      outvals[i] = (unsigned)(nclose);
  }

  if(waserr) {
    ThrowError("Error in generatePoissonSamples");
  }

  return outvals;
}


/*.......................................................................
 * Incomplete gamma function.
 * Cribbed from NR.
 */
int Sampler::gammp(double a, double x, double *val)
{
  double gln;
  int waserr=0;
  
  if (x < 0.0 || a <= 0.0) {
    fprintf(stderr,"gammp: Invalid arguments.\n");
    return 1;
  }

  if (x < (a+1.0)) {
    waserr = gser(val,a,x,&gln);
  } else {
    waserr = gcf(val,a,x,&gln);
    *val = 1.0-*val;
  }

  return waserr;
}

/*.......................................................................
 * Evaluates the incomplete gamma function by its series representation.
 * Cribbed from NR.
 */
static int gser(double *gamser, double a, double x, double *gln)
{
#define GSER_ITMAX 1000
#define GSER_EPS 3.0e-7
  int n;
  double sum,del,ap;

  *gln=Sampler::lnGamma(a);
  if (x <= 0.0) {
    if (x < 0.0) {
      fprintf(stderr,"gser: x less than 0.\n");
      return 1;
    }
    *gamser=0.0;
    return 0;
  }
  else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=GSER_ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*GSER_EPS) {
        *gamser=sum*exp(-x+a*log(x)-(*gln));
        return 0;
      }
    }
    fprintf(stderr,"gser: a too large, GSER_ITMAX too small.\n");
    return 1;
  }
}


/*.......................................................................
 * Evaluate the incomplete gamma function by continued fraction
 * representation.
 *
 * Cribbed from NR.
 */
static int gcf(double *gammcf, double a, double x, double *gln)
{
#define GCF_FPMIN 1.0e-30
#define GCF_ITMAX 100
#define GCF_EPS 3.0e-7
  int i;
  double an,b,c,d,del,h;

  *gln=Sampler::lnGamma(a);
  b=x+1.0-a;
  c=1.0/GCF_FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=GCF_ITMAX;i++) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < GCF_FPMIN) d=GCF_FPMIN;
    c=b+an/c;
    if (fabs(c) < GCF_FPMIN) c=GCF_FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < GCF_EPS) break;
  }
  if (i > GCF_ITMAX) {
    //    fprintf(stderr,"gcf: a too large, GCF_ITMAX too small\n.");
    return 1;
  }
  *gammcf=exp(-x+a*log(x)-(*gln))*h;
  return 0;
}

/**.......................................................................
 * Return the Chebyshev fit to the complementary Error function:
 *
 *
 *                            2      / +inf    2     
 * erfc(x)  = 1 - erf(x) = --------  |   exp(-t ) dt 
 *                         sqrt(pi)  / x             
 *                               
 */
double Sampler::erfcc(double x)
{
  //------------------------------------------------------------
  // Trap extremal values which don't require computation, so we can
  // validly call this function with infinities
  //------------------------------------------------------------

  if(x == negInf_)
    return 1.0;

  if(x == posInf_)
    return 0.0;

  double t,z,ans;

  z=fabs(x);
  t=1.0/(1.0+0.5*z);
  ans=t*exp(-z*z-1.26551223 
	    + t*( 1.00002368 + t*( 0.37409196 + t*(0.09678418 + 
              t*(-0.18628806 + t*( 0.27886807 + t*(-1.13520398+
	      t*( 1.48851587 + t*(-0.82215223 + t*0.17087277)))))))));
  return x >= 0.0 ? ans : 2.0-ans;
}

/*.......................................................................
 * This function returns a random number, from a gaussian distribution
 * of standard deviation, num.
 */
double Sampler::generateGaussianSample(double sigma)
{
  double rrad;
  double aval, bval, rand_num;
  
  double sample;

  // Acquire two uniform random numbers between -1 and 1.
    
  do {
    aval=(double)(rand())/RAND_MAX*2-1;
    bval=(double)(rand())/RAND_MAX*2-1;
    
    // The Box-Muller transformation to convert uniform to gaussian
    // deviates requires that the two deviates be converted to a
    // radius squared. The value of the radius must be less than one
    // and not equal to zero.
    
    rrad = aval*aval + bval*bval;
    
  } while(rrad >= 1 || rrad == 0.0);
  
  // Apply the Box-Muller transformation to turn the random square
  // radius found above into two Gaussian deviates.
  
  sample = sigma * aval * sqrt(-2.0 * log((double) rrad)/rrad);

  return sample;
}

/*.......................................................................
 * This function returns a random number, from a truncated gaussian
 * distribution of the specified standard deviation and mean.
 */
double Sampler::generateTruncatedGaussianSample(double mean, double sigma, double xmin, double xmax)
{
  double val;
  do {
    val = generateGaussianSample(sigma) + mean;
  } while(val < xmin || val > xmax);

  return val;
}

/**.......................................................................
 * Return the normalization factor for a truncated gaussian variate.
 *
 * The normalizing factor for a gaussian is:
 *
 *   sqrt(2 pi) * sig
 *
 * For a truncated Gaussian distribution, therefore, the normalizing
 * factor must be multiplied by the integral from xmin to xmax. 
 */
double Sampler::truncatedGaussIntegral(double mean, double sigma, double xmin, double xmax, double& norm, bool& normIsCalculated)
{
  //------------------------------------------------------------
  // Don't recompute the integral unless we have to
  //------------------------------------------------------------

  if(!normIsCalculated) {
    normIsCalculated = true;
    norm = gaussCdf(xmax, mean, sigma) - gaussCdf(xmin, mean, sigma);
  }
  
  return norm;
}

/*.......................................................................
 * This function returns a random number, from a gaussian distribution
 * of standard deviation, num.
 */
double
Sampler::generateGaussianSampleNonStatic(double sigma)
{
  double rrad;
  double aval, bval, rand_num;
  
  double sample;

  // Acquire two uniform random numbers between -1 and 1.
    
  do {
    aval=(double)(rand())/RAND_MAX*2-1;
    bval=(double)(rand())/RAND_MAX*2-1;
    
    // The Box-Muller transformation to convert uniform to gaussian
    // deviates requires that the two deviates be converted to a
    // radius squared. The value of the radius must be less than one
    // and not equal to zero.
    
    rrad = aval*aval + bval*bval;
    
  } while(rrad >= 1 || rrad == 0.0);
  
  // Apply the Box-Muller transformation to turn the random square
  // radius found above into two Gaussian deviates.
  
  sample = sigma * aval * sqrt(-2.0 * log((double) rrad)/rrad);

  return sample;
}

/**.......................................................................
 * This function returns a random number, from a gaussian distribution
 * of standard deviation, num.
 */
std::vector<double> Sampler::
generateGaussianSamples(double sigma, unsigned nSamp)
{
  std::vector<double> samples(nSamp);

  for(unsigned i=0; i < nSamp; i++)
    samples[i] = generateGaussianSample(sigma);

  return samples;
}

std::vector<double> 
Sampler::generateTruncatedGaussianSamples(double mean, double sigma, double xmin, double xmax, unsigned nSamp)
{
  std::vector<double> samples(nSamp);

  for(unsigned i=0; i < nSamp; i++)
    samples[i] = generateTruncatedGaussianSample(mean, sigma, xmin, xmax);
  
  return samples;
}

/**.......................................................................
 * Return a uniform sample
 */
double Sampler::generateUniformSample(double xmin, double xmax)
{
  return xmin + (xmax - xmin)* rand()/RAND_MAX;
}

/**.......................................................................
 * Return uniform samples
 */
std::vector<double> Sampler::
generateUniformSamples(double xmin, double xmax, unsigned nSamp)
{
  std::vector<double> samples(nSamp);

  for(unsigned i=0; i < nSamp; i++)
    samples[i] = generateUniformSample(xmin, xmax);

  return samples;
}

/*.......................................................................
 * Series representation for the gamma function.
 */
double Sampler::lnGamma(double xx)
{
  double x,y,tmp,ser;

  static double cof[6]={76.18009172947146,-86.50532032941677,
                          24.01409824083091,-1.231739572450155,
                          0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

/**.......................................................................
 * Return the log of the factorial of n.  For large n, we use 
 *
 *   n! = Gamma(n+1)
 *
 * with a series approximation for the gamma function
 */
double Sampler::lnFactrl(unsigned n)
{
  static int ntop=4;
  int j;
  
  if(n > 32) {
    return lnGamma(n+1.0);
  }

  double sum = 0.0;
  for(unsigned i=0; i < n;i++) {
    sum += log((double)(n-i));
  }

  return sum;
}

/**.......................................................................
 * Return the natural log of the pdf of a user-specified variate
 */
double Sampler::lnUserPdf(double x)
{
  double updf = userPdf(x);

  // Note -- I'm deliberately not checking for log(0.0), since the
  // compiler will return -inf in this case and that's what I want

  return log(updf);
}

/**.......................................................................
 * Return the pdf of a user-specified variate
 */
double Sampler::userPdf(double x)
{
  //------------------------------------------------------------
  // If we have more than one point, find bracketing values and
  // interpolate
  //------------------------------------------------------------

  if(nPt_ > 1) {
    double dx = x_[1] - x_[0];
    double xmin = x_[0] - dx/2;
    double xmax = x_[nPt_-1] + dx/2;

    if(x < xmin || x >= xmax)
      return 0.0;

    unsigned lo = floor((x-xmin)/dx);
    unsigned hi;
    
    if(lo==nPt_-1) {
      hi = nPt_-1;
      lo = hi-1;
    } else {
      hi = lo+1;
    }
    
    return y_[lo] + (y_[hi] - y_[lo])/(x_[hi] - x_[lo]) * (x - x_[lo]);

    //------------------------------------------------------------
    // Else see if we are within EPSFRAC of the one specified point
    //------------------------------------------------------------

  } else if(nPt_ == 1) {

    return fabs((x - x_[0])/x_[0]) < EPSFRAC ? y_[0] : 0.0;

    //------------------------------------------------------------
    // Else assume zero
    //------------------------------------------------------------

  } else {
    return 0.0;
  }
}

/**.......................................................................
 * Return the pdf of the chi-square variate (x), for nDof (k) degrees of
 * freedom.
 *
 *                           1
 * Note:  f(x|k) =  -------------------- x^(k/2-1) e^(-x/2)
 *                  2^(k/2) * Gamma(k/2)
 *
 * So that:
 *
 *   ln(f) = (k/2-1) * ln(x) - x/2 - (k/2) * ln(2) - ln(Gamma(k/2))
 */
double Sampler::lnChisqPdf(double chisq, unsigned nDof)
{
  //  COUT("Inside lnChisqpdf with chisq = " << chisq << " nDof = " << nDof);

  double k2  = ((double)nDof)/2;
  double lk2 = log(k2);
  double lg  = lnGamma(k2);
  double lx  = log(chisq);
  static double l2 = log(2.0);
  
  double lf = (k2 - 1) * lx - chisq/2 - k2 * l2 - lg;

  return lf;
}

/**.......................................................................
 * Return the pdf of the chi-square variate, for nDof degrees of
 * freedom.
 *
 *                           1
 * Note:  f(x|k) =  -------------------- x^(k/2-1) e^(-x/2)
 *                  2^(k/2) * Gamma(k/2)
 *
 * So that:
 *
 *   ln(f) = (k/2-1) * ln(x) - x/2 - (k/2) * ln(2) - ln(Gamma(k/2))
 */
double Sampler::chisqPdf(double chisq, unsigned nDof)
{
  return exp(lnChisqPdf(chisq, nDof));
}

/**.......................................................................
 * Return the natural log of the pdf of a gaussian variate
 */
double Sampler::lnTruncatedGaussPdf(double x, double mean, double sigma, double xmin, double xmax, double& norm, bool& normIsCalculated)
{
  if(x >= xmin && x <= xmax) {
    double var    = (x - mean)/sigma;
    double arg    = (var * var)/2;
    double prefac = sqrt2pi_ * truncatedGaussIntegral(mean, sigma, xmin, xmax, norm, normIsCalculated);
    normIsCalculated = true;
    
    return -arg - log(prefac);
  } else {
    return log(0.0);
  }
}

/**.......................................................................
 * Return the pdf of a gaussian variate
 */
double Sampler::truncatedGaussPdf(double x, double mean, double sigma, double xmin, double xmax, double& norm, bool& normIsCalculated)
{
  if(x >= xmin && x <= xmax) {
    double var    = (x - mean)/sigma;
    double arg    = (var * var)/2;
    double prefac = sqrt2pi_ * truncatedGaussIntegral(mean, sigma, xmin, xmax, norm, normIsCalculated);
    
    return exp(-arg)/prefac;
  } else {
    return 0.0;
  }
}

/**.......................................................................
 * Return the natural log of the pdf of a gaussian variate
 */
double Sampler::lnGaussPdf(double x, double mean, double sigma)
{
  double var    = (x - mean)/sigma;
  double arg    = (var * var)/2;
  double prefac = sqrt2pi_ * sigma;

  return -arg - log(prefac);
}

/**.......................................................................
 * Return the pdf of a gaussian variate
 */
double Sampler::gaussPdf(double x, double mean, double sigma)
{
  double var    = (x - mean)/sigma;
  double arg    = (var * var)/2;
  double prefac = sqrt2pi_ * sigma;

  return exp(-arg)/prefac;
}

/**.......................................................................
 * Return the Poisson pdf for k events given a mean of lambda:
 *
 *          e^-l * l^k
 * p(k,l) = ----------
 *               k!
 *
 * So ln(p) = -l + k*ln(l) - ln(k!)
 *
 */ 
double Sampler::lnPoissPdf(unsigned k, double lambda)
{
  double lnp = -lambda + k * log(lambda) - lnFactrl(k);
  return lnp;
}

double Sampler::poissPdf(unsigned k, double lambda)
{
  return exp(lnPoissPdf(k, lambda));
}

double Sampler::lnUniformPdf(double x, double xMin, double xMax)
{
  if(x < xMin || x >= xMax) {
    return log(0.0); // Note that I'm deilberately returning -inf
		     // here.  This will produce valid comparisons
		     // (ie, -inf < some number) which is what I
		     // intend.
  } else {
    //    return -log(xMax - xMin);
    return 0.0;
  }
}

double Sampler::uniformPdf(double x, double xMin, double xMax)
{
  if(x < xMin || x >= xMax) {
    return 0.0;
  } else {
    //    return 1.0/(xMax - xMin);
    return 1.0;
  }
}

/**.......................................................................
 * Return the cdf of a user-specified variate
 */
double Sampler::userCdf(double x)
{
  //------------------------------------------------------------
  // If we have more than one point, find bracketing values and
  // interpolate
  //------------------------------------------------------------

  if(nPt_ > 1) {

    double dx = x_[1] - x_[0];
    double xmin = x_[0] - dx/2;
    double xmax = x_[nPt_-1] + dx/2;
    
    // For the CDF, return 0.0 for x < xmin, 1.0 for all values >=
    // xmax

    if(x < xmin)
      return 0.0;

    if(x >= xmax)
      return 1.0;
    
    unsigned lo = floor((x-xmin)/dx);
    unsigned hi;

    if(lo==nPt_-1) {
	hi = nPt_-1;
	lo = hi-1;
      } else {
	hi = lo+1;
      }

    return yInt_[lo] + (yInt_[hi] - yInt_[lo])/(x_[hi] - x_[lo]) * (x - x_[lo]);

    //------------------------------------------------------------
    // Else see if we are within EPSFRAC of the one specified point
    //------------------------------------------------------------

  } else if(nPt_ == 1) {

    return fabs((x - x_[0])/x_[0]) < EPSFRAC ? yInt_[0] : 0.0;

    //------------------------------------------------------------
    // Else assume zero
    //------------------------------------------------------------

  } else {
    return 0.0;
  }
}

/**.......................................................................
 * Return the probability to exceed a given value of chisq
 */
double Sampler::chisqCdf(double chisq, unsigned nDof)
{
  double a = (double)(nDof)/2;
  double x = chisq/2;
  
  double val;

  //------------------------------------------------------------
  // Try to compute the exact expression.  If we can't, return the
  // gaussian approximation instead
  //------------------------------------------------------------

  if(gammp(a, x, &val)) {
    double k = (double) nDof;
    return gaussCdf(chisq, k, sqrt(2*k));
  }

  return val;
}

/**.......................................................................
 * Return the cdf of a gaussian variate.  We have:
 * 
 *             2      / x       2
 * erf(x) = --------  |   exp(-t ) dt
 *          sqrt(pi)  / 0
 *
 * and
 *
 *                            2      / +inf      2     
 * erfc(x)  = 1 - erf(x) = --------  |     exp(-t ) dt 
 *                         sqrt(pi)  / x             
 *                               
 *                               
 * The gaussian cdf is given by: 
 *                               
 *
 *                                                  2       
 *                  1          / x        /   (v - m)   \     
 * cdf(x) =  ----------------  |     exp  | - ---------  | dv 
 *           sig * sqrt(2*pi)  / -inf     \  2 sig*sig  /     
 *                                                         
 *                                                         
 *                                                        
 *             1      / x           2                  1      / +inf        2                       
 *        = --------  |     exp  (-t ) dt  =  1.0 - --------  |     exp  (-t ) dt                                             
 *          sqrt(pi)  / -inf                        sqrt(pi)  / x                                                                
 *                       
 *                               
 *                          v - m
 * with substitution t = -----------
 *                       sqrt(2) sig
 *
 *
 *                   erfc(x) 
 * Thus cdf(x) = 1 - -------
 *                      2
 */                                                           
double Sampler::gaussCdf(double x, double mean, double sigma)
{
  double arg = (x - mean)/(sqrt(2.0)*sigma);
  return (2.0 - erfcc(arg))/2;
}

/**.......................................................................
 * Return the cdf of a truncated gaussian variate.
 *
 * Because this function is not required for fast evaluation, we can
 * disregard any normalizing factor by computing the integral from
 * xmin to xmax each time (highLim - lowLim) below
 */
double Sampler::truncatedGaussCdf(double x, double mean, double sigma, double xmin, double xmax, double& norm, bool& normIsCalculated)
{
  //------------------------------------------------------------
  // If the requested value is below our cutoff, the cumulative
  // probability is zero
  //------------------------------------------------------------

  if(x < xmin)
    return 0.0;

  //------------------------------------------------------------
  // If the requested value is above our cutoff, the cumulative
  // probability is 1.0
  //------------------------------------------------------------
  
  if(x > xmax)
    return 1.0;

  double lowLim  = gaussCdf(xmin, mean, sigma);
  double val     = gaussCdf(x, mean, sigma);
  double normVal = truncatedGaussIntegral(mean, sigma, xmin, xmax, norm, normIsCalculated);

  return (val - lowLim) / normVal;
}

double Sampler::poissCdf(unsigned k, double lambda)
{
  double val;
  gammp((double)k, lambda, &val);
  return 1.0 - val;
}

double Sampler::uniformCdf(double x, double xMin, double xMax)
{
  if(x <= xMin)
    return 0.0;
  else if(x > xMax)
    return 1.0;
  else {
    return (x - xMin)/(xMax - xMin);
  }
}

/**.......................................................................
 * Return the probability to exceed a given value of a user-specified variate
 */
double Sampler::userPte(double x)
{
  return 1.0 - userCdf(x);
}

/**.......................................................................
 * Return the probability to exceed a given value of chisq
 */
double Sampler::chisqPte(double chisq, unsigned nDof)
{
  return 1.0 - chisqCdf(chisq, nDof);
}

/**.......................................................................
 * Return the pte for a gaussian variate
 */
double Sampler::gaussPte(double x, double mean, double sigma)
{
  return 1.0 - gaussCdf(x, mean, sigma);
}

/**.......................................................................
 * Return the pte for a truncated gaussian variate
 */
double Sampler::truncatedGaussPte(double x, double mean, double sigma, double xmin, double xmax, double& norm, bool& normIsCalculated)
{
  return 1.0 - truncatedGaussCdf(x, mean, sigma, xmin, xmax, norm, normIsCalculated);
}

/**.......................................................................
 * Return the Poisson pdf for k events given a mean of lambda:
 *
 *          e^-l * l^k
 * p(k,l) = ----------
 *               k!
 *
 * So ln(p) = -l + k*ln(l) - ln(k!)
 *
 */ 
double Sampler::poissPte(unsigned k, double lambda)
{
  return 1.0 - poissCdf(k, lambda);
}

double Sampler::uniformPte(double x, double xMin, double xMax)
{
  return 1.0 - uniformCdf(x, xMin, xMax);
}

void Sampler::histogram(std::vector<double>& vals, unsigned nbin, std::vector<float>& x, std::vector<float>& y)
{
  double xmin, xmax;

  for(unsigned i=0; i < vals.size(); i++) {

    float val = vals[i];

    if(i==0) {
      xmin = xmax = val;
    }

    xmin = (xmin < val) ? xmin : val;
    xmax = (xmax > val) ? xmax : val;
  }

  float dx    = (xmax - xmin) / (nbin - 1);
  float hxmin = xmin - dx/2;

  x.resize(nbin);
  y.resize(nbin);

  for(unsigned i=0; i < x.size(); i++) {
    x[i] = hxmin + dx/2 + i*dx;
    y[i] = 0.0;
  }

  float ymax = 0.0;
  for(unsigned i=0; i < vals.size(); i++) {

    float val = vals[i];
    unsigned ind = (unsigned)((val - hxmin) / dx);

    //    COUT("hxmin = " << hxmin << " dx = " << dx << " ind = " << ind);

    y[ind] += 1.0;

    ymax = (ymax > y[ind]) ? ymax : y[ind];
  }

  for(unsigned i=0; i < x.size(); i++) {
    y[i] /= ymax;
  }

  COUT("Leaving...");
}

SAMPLER_FN(Sampler::gaussian)
{
  double mean  = args[0];
  double sigma = args[1];

  double arg = (x - mean) * (x - mean) / (2*sigma*sigma);

  return exp(-arg);
}

/**.......................................................................
 * Generate a sample from a multivariate normal distribution
 */
Vector<double> 
Sampler::generateMultiVariateGaussianSample(Vector<double>& val, 
					    Vector<double>& mean, 
					    Vector<double>& sigma,
					    Matrix<double>& correl)
{
  Vector<double> sample = val;

  // First calculate the covariance matrix from the correlation matrix:
  //
  // Given:
  //
  //        |  1   r01  r02 |
  // Corr = | r10   1   r12 |
  //        | r20  r21   1  |
  //
  //       | sig0^2    0      0   | |  1   r01  r02 | | sig0^2    0      0   |
  // Cov = |    0   sig1^2    0   | | r10   1   r12 | |    0   sig1^2    0   |  
  //       |    0      0   sig2^2 | | r20  r21   1  | |    0      0   sig2^2 |

  Matrix<double> s(sigma.size(), sigma.size());

  for(unsigned iSig=0; iSig < sigma.size(); iSig++)
    s[iSig][iSig] = sigma[iSig];

  Matrix<double> cov = s * correl * s;
  Matrix<double> invCov = cov.inverse();

  //  COUT("det(cov) =  " << cov.det());

  unsigned nParam = sigma.size();

  for(unsigned iParam=0; iParam < nParam; iParam++) {
    multiVariateSampleIterator(iParam, sample, mean, invCov);
  }

  return sample;
}

/**.......................................................................
 * Generate a sample from a multivariate normal distribution, passing
 * in a pre-calculated inverse covariance matrix
 */
Vector<double> 
Sampler::generateMultiVariateGaussianSample(Vector<double>& val, 
					    Vector<double>& mean, 
					    Matrix<double>& invCov)
{
  Vector<double> sample = val;
  unsigned nParam = sample.size();

  for(unsigned iParam=0; iParam < nParam; iParam++) 
    multiVariateSampleIterator(iParam, sample, mean, invCov);

  return sample;
}

/**.......................................................................
 * Replace val[k] with a sample from its conditional multivariate
 * normal distribution
 */
void 
Sampler::multiVariateSampleIterator(unsigned k,
				    Vector<double>& val, 
				    Vector<double>& mean, 
				    Matrix<double>& invCov)
{
  Vector<double> xm = val - mean;

  double M   = xm * invCov * xm;
  double Ckk = invCov[k][k];
  double muk = mean[k];
  unsigned n = mean.size();

  //------------------------------------------------------------
  // Now calculate A_k:
  //------------------------------------------------------------

  double Ak = 0.0;

  //  COUT("xm = " << std::endl << xm);
  //  COUT("val = " << std::endl << val);

  for(unsigned i=0; i < n; i++) {
    if(i != k)
      Ak += xm[i] * invCov[i][k];
  }

  //  COUT("Ak = " << Ak << " muk = " << muk);

  Ak = 2 * (Ak/Ckk - muk);

  //------------------------------------------------------------
  // Now generate a sample from a distribution with sigma =
  // 1.0/sqrt(Ckk) and mean -Ak/2
  //------------------------------------------------------------

  //  COUT("About to generate sample with Ak = " << Ak << " Ckk = " << Ckk);

  val[k] = -Ak/2 + generateGaussianSample(1.0/sqrt(Ckk));

  //  COUT("val is now = " << std::endl << val);
}

/**.......................................................................
 *
 * The pdf of a multivariate gaussian is given by: 
 *
 *   _            1               1  _   _          _   _
 * p(x) = ----------------- exp - - (x - m)^T C^-1 (x - m)
 *        (2pi)^k/2 |C|^1/2       2
 */
double Sampler::lnMultiVariateGaussPdf(Vector<double>& x, 
				       Vector<double>& mean, 
				       Vector<double>& sigma, 
				       Matrix<double>& correl)
{
  // First calculate the covariance matrix from the correlation matrix:
  //
  // Given:
  //
  //        |  1   r01  r02 |
  // Corr = | r10   1   r12 |
  //        | r20  r21   1  |
  //
  //       | sig0^2    0      0   | |  1   r01  r02 | | sig0^2    0      0   |
  // Cov = |    0   sig1^2    0   | | r10   1   r12 | |    0   sig1^2    0   |  
  //       |    0      0   sig2^2 | | r20  r21   1  | |    0      0   sig2^2 |

  Matrix<double> s(sigma.size(), sigma.size());

  for(unsigned iSig=0; iSig < sigma.size(); iSig++)
    s[iSig][iSig] = sigma[iSig];

  Matrix<double> cov = s * correl * s;
  Matrix<double> invCov = cov.inverse();
  double detC = cov.determinant();

  COUT("invCov = " << std::endl << invCov);
  COUT("detC = " << detC);

  return lnMultiVariateGaussPdf(x, mean, invCov, detC);
}

/**.......................................................................
 *
 * The pdf of a multivariate gaussian is given by: 
 *
 *   _            1               1  _   _          _   _
 * p(x) = ----------------- exp - - (x - m)^T C^-1 (x - m)
 *        (2pi)^k/2 |C|^1/2       2
 */
double Sampler::lnMultiVariateGaussPdf(Vector<double>& x, 
				       Vector<double>& mean, 
				       Matrix<double>& invCov,
				       double detC)
{
  Vector<double> xm = x - mean;

  double lnExp  = -(xm * invCov * xm)/2;
  double lnNorm = -(double)(x.size())/2 * log(2*M_PI) - 0.5 * log(detC);

  return lnNorm + lnExp;
}

double Sampler::multiVariateGaussPdf(Vector<double>& x, 
				     Vector<double>& mean, 
				     Matrix<double>& invCov,
				     double detC)
{
  return exp(lnMultiVariateGaussPdf(x, mean, invCov, detC));
}

double Sampler::multiVariateGaussPdf(Vector<double>& x, 
				     Vector<double>& mean, 
				     Vector<double>& sigma, 
				     Matrix<double>& correl)
{
  return exp(lnMultiVariateGaussPdf(x, mean, sigma, correl));
}
