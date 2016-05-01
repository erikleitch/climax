#include "gcp/util/Distribution.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Variate.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Distribution::Distribution() 
{
  type_         = DIST_UNKNOWN;
  requiredMask_ = PARAM_UNKNOWN;
  paramMask_    = PARAM_NONE;
  pdfFn_        = 0x0;
  cdfFn_        = 0x0;
  sampleFn_     = 0x0;

  setType(DIST_UNSPEC);
}

/**.......................................................................
 * Destructor.
 */
Distribution::~Distribution() {}

void Distribution::setType(Type type)
{
  //  COUT("this = " << this << " Inside setType with type = " << type << " unspec = " << DIST_UNSPEC);

  switch(type) {
  case DIST_UNSPEC:
    requiredMask_ = PARAM_NONE;
    pdfFn_    = Distribution::unspecPdf;
    cdfFn_    = Distribution::unspecCdf;
    sampleFn_ = Distribution::unspecSample;
    break;
  case DIST_USER:
    requiredMask_ = PARAM_USER_DYDX;
    pdfFn_    = Distribution::userPdf;
    cdfFn_    = Distribution::userCdf;
    sampleFn_ = Distribution::userSample;
    break;
  case DIST_GAUSS:
    requiredMask_  = PARAM_GAUSS_MEAN;
    requiredMask_ |= PARAM_GAUSS_SIGMA;
    pdfFn_    = Distribution::gaussPdf;
    cdfFn_    = Distribution::gaussCdf;
    sampleFn_ = Distribution::gaussSample;
    break;
  case DIST_TRUNC_GAUSS:
    requiredMask_  = PARAM_GAUSS_MEAN;
    requiredMask_ |= PARAM_GAUSS_SIGMA;
    pdfFn_    = Distribution::truncGaussPdf;
    cdfFn_    = Distribution::truncGaussCdf;
    sampleFn_ = Distribution::truncGaussSample;
    break;
  case DIST_POISS:
    requiredMask_  = PARAM_POISS_MEAN;
    pdfFn_    = Distribution::poissPdf;
    cdfFn_    = Distribution::poissCdf;
    sampleFn_ = Distribution::poissSample;
    break;
  case DIST_CHISQ:
    requiredMask_ = PARAM_CHISQ_NDOF;
    pdfFn_    = Distribution::chisqPdf;
    cdfFn_    = Distribution::chisqCdf;
    sampleFn_ = Distribution::chisqSample;
    break;
  case DIST_JOINT_GAUSS:
    requiredMask_ = PARAM_CHISQ_NDOF;
    pdfFn_    = Distribution::jointGaussPdf;
    cdfFn_    = Distribution::jointGaussCdf;
    sampleFn_ = Distribution::jointGaussSample;
    break;
  case DIST_UNIFORM:
    requiredMask_  = PARAM_UNIFORM_XMIN;
    requiredMask_ |= PARAM_UNIFORM_XMAX;
    pdfFn_    = Distribution::uniformPdf;
    cdfFn_    = Distribution::uniformCdf;
    sampleFn_ = Distribution::uniformSample;
    break;
  default:
    ThrowError("Unrecognized distribution type: " << type);
    break;
  }

  type_      = type;
  paramMask_ = PARAM_NONE;
}

//------------------------------------------------------------
// Parameters for uniform distributions
//------------------------------------------------------------

void Distribution::setUniformXMin(Variate& xMin)
{
  setUniformXMin(xMin.val_);
}

void Distribution::setUniformXMax(Variate& xMax)
{
  setUniformXMax(xMax.val_);
}

void Distribution::setUniformXMin(double xMin)
{
  checkType(DIST_UNIFORM);
  uniform_.xMin_ = xMin;
  paramMask_ |= PARAM_UNIFORM_XMIN;

  invalidateTruncatedGaussianParameters();
}

void Distribution::setUniformXMax(double xMax)
{
  checkType(DIST_UNIFORM);
  uniform_.xMax_ = xMax;
  paramMask_ |= PARAM_UNIFORM_XMAX;

  invalidateTruncatedGaussianParameters();
}

double Distribution::getUniformXMin()
{
  checkType(DIST_UNIFORM);
  return uniform_.xMin_;
}

double Distribution::getUniformXMax()
{
  checkType(DIST_UNIFORM);
  return uniform_.xMax_;
}

//------------------------------------------------------------
// Parameters for poisson distributions
//------------------------------------------------------------

void Distribution::setPoissMean(double mean)
{
  checkType(DIST_POISS);
  poiss_.mean_ = mean;
  paramMask_ |= PARAM_POISS_MEAN;
}

double Distribution::getPoissMean()
{
  checkType(DIST_POISS);
  return poiss_.mean_;
}

//------------------------------------------------------------
// Parameters for gaussian distributions
//------------------------------------------------------------

void Distribution::setGaussMean(Variate& var)
{
  setGaussMean(var.val_);
}

void Distribution::setGaussSigma(Variate& var)
{
  setGaussSigma(var.val_);
}

double Distribution::getGaussMean()
{
  checkType(DIST_GAUSS);
  return  gauss_.mean_;
}

double Distribution::getGaussSigma()
{
  checkType(DIST_GAUSS);
  return  gauss_.sigma_;
}

void Distribution::setGaussMean(double mean)
{
  checkType(DIST_GAUSS);
  gauss_.mean_ = mean;
  paramMask_ |= PARAM_GAUSS_MEAN;

  invalidateTruncatedGaussianParameters();
}

void Distribution::setGaussSigma(double sigma)
{
  checkType(DIST_GAUSS);
  gauss_.sigma_ = sigma;
  paramMask_ |= PARAM_GAUSS_SIGMA;

  invalidateTruncatedGaussianParameters();
}

//------------------------------------------------------------
// Parameters for chisq distributions
//------------------------------------------------------------

void Distribution::setChisqNdof(unsigned nDof)
{
  checkTypes(DIST_CHISQ, DIST_JOINT_GAUSS);
  chisq_.nDof_ = nDof;
  paramMask_ |= PARAM_CHISQ_NDOF;
}

unsigned Distribution::chisqNdof()
{
  checkTypes(DIST_CHISQ, DIST_JOINT_GAUSS);
  return chisq_.nDof_;
}

//------------------------------------------------------------
// Parameters for user-specified distributions
//------------------------------------------------------------

/**.......................................................................
 * Methods for specifying the sampling function, as two arrays
 * specifying dy/dx
 */
void Distribution::setdYdX(unsigned n, double* x, double* y)
{
  checkType(DIST_USER);
  user_.sampler_.setdYdX(n, x, y);
  paramMask_ |= PARAM_USER_DYDX;
}

void Distribution::setdYdX(std::vector<double>& x, std::vector<double>& y)
{
  checkType(DIST_USER);
  user_.sampler_.setdYdX(x, y);
  paramMask_ |= PARAM_USER_DYDX;
}


/**.......................................................................
 * Method for specifying the sampling dy/dx as a function
 */
void Distribution::setdYdX(SAMPLER_FN(fn), double xMin, double xMax, double dx, double* args)
{
  checkType(DIST_USER);
  user_.sampler_.setdYdX(fn, xMin, xMax, dx, args);
  paramMask_ |= PARAM_USER_DYDX;
}

/**.......................................................................
 * Method for specifying the integral of dy/dx (Y(x' > x))
 * directly
 */
void Distribution::setYX(SAMPLER_FN(fn))
{
  checkType(DIST_USER);
  user_.sampler_.setYX(fn);
  paramMask_ |= PARAM_USER_DYDX;
}

void Distribution::checkType(Type type)
{
  unsigned tThis = (unsigned) type_;
  unsigned tThat = (unsigned) type;

  if(!(tThis & tThat)) {
    ThrowError("You are attempting to specify parameters for distribution type " 
	       << type << " but this distribution is of type " << type_);
  }
}

void Distribution::checkTypes(Type type1, Type type2)
{
  if(!(type1 == type_ || type2 == type_)) {
    ThrowError("You are attempting to specify parameters for distribution type " 
	       << type1 << " or " << type2 << " but this distribution is of type " << type_);
  }
}

Distribution::Type Distribution::getType()
{
  return type_;
}

void Distribution::checkParams()
{
  if((paramMask_ & requiredMask_) != requiredMask_)
    ThrowError("Required parameters have not been specified for this distribution");
}

/**.......................................................................
 * Return the value of the pdf at the specified value of x
 */
Probability Distribution::pdf(double x)
{
  if(pdfFn_ == 0) {
    ThrowError("No pdf has been specified");
  }

  checkParams();

  return pdfFn_(x, this);
}

/**.......................................................................
 * Return the value of the cdf at the specified value of x
 */
Probability Distribution::cdf(double x)
{
  if(cdfFn_ == 0) {
    ThrowError("No cdf has been specified");
  }

  checkParams();

  return cdfFn_(x, this);
}

Probability Distribution::pte(double x)
{
  Probability prob;
  prob.setValue(1.0 - cdf(x).value());
  return prob;
}

double Distribution::sample()
{
  if(sampleFn_ == 0) {
    ThrowError("No sample function has been specified");
  }

  checkParams();

  return sampleFn_(this);
}

//-----------------------------------------------------------------------
// Unspecified distribution functions
//-----------------------------------------------------------------------

EVAL_FN(Distribution::unspecPdf)
{
  Probability prob;
  prob.setLnValue(0.0);
  return prob;
}

EVAL_FN(Distribution::unspecCdf)
{
  Probability prob;
  ThrowError("CDF for an unspecified distribution function makes no sense");
  return prob;
}

SAMP_FN(Distribution::unspecSample)
{
  ThrowError("Sample FN for an unspecified distribution function makes no sense");
  return 0.0;
}

//-----------------------------------------------------------------------
// User-specified functions
//-----------------------------------------------------------------------

EVAL_FN(Distribution::userPdf)
{
  Probability prob;
  prob.setLnValue(dist->user_.sampler_.lnUserPdf(x));
  return prob;
}

EVAL_FN(Distribution::userCdf)
{
  Probability prob;
  prob.setValue(dist->user_.sampler_.userCdf(x));
  return prob;
}

SAMP_FN(Distribution::userSample)
{
  return dist->user_.sampler_.generateSample();
}

//-----------------------------------------------------------------------
// Gaussian functions
//-----------------------------------------------------------------------

EVAL_FN(Distribution::gaussPdf)
{
  Probability prob;
  prob.setLnValue(Sampler::lnGaussPdf(x, dist->gauss_.mean_, dist->gauss_.sigma_));
  return prob;
}

EVAL_FN(Distribution::gaussCdf)
{
  Probability prob;
  prob.setValue(Sampler::gaussCdf(x, dist->gauss_.mean_, dist->gauss_.sigma_));
  return prob;
}

SAMP_FN(Distribution::gaussSample)
{
  COUT("Generating gaussian sample from " << dist->gauss_.sigma_ << " mean = " << dist->gauss_.mean_);
  return Sampler::generateGaussianSample(dist->gauss_.sigma_) + dist->gauss_.mean_;
}

//-----------------------------------------------------------------------
// Truncated Gaussian functions
//-----------------------------------------------------------------------

EVAL_FN(Distribution::truncGaussPdf)
{
  Probability prob;
  prob.setLnValue(Sampler::lnTruncatedGaussPdf(x, dist->gauss_.mean_, dist->gauss_.sigma_, dist->uniform_.xMin_, dist->uniform_.xMax_, dist->gauss_.norm_, dist->gauss_.normIsCalculated_));
  return prob;
}

EVAL_FN(Distribution::truncGaussCdf)
{
  Probability prob;
  prob.setValue(Sampler::truncatedGaussCdf(x, dist->gauss_.mean_, dist->gauss_.sigma_, dist->uniform_.xMin_, dist->uniform_.xMax_, dist->gauss_.norm_, dist->gauss_.normIsCalculated_));
  return prob;
}

SAMP_FN(Distribution::truncGaussSample)
{
  return Sampler::generateTruncatedGaussianSample(dist->gauss_.mean_, dist->gauss_.sigma_, dist->uniform_.xMin_, dist->uniform_.xMax_);
}

//-----------------------------------------------------------------------
// Poisson functions
//-----------------------------------------------------------------------

EVAL_FN(Distribution::poissPdf)
{
  Probability prob;
  prob.setLnValue(Sampler::lnPoissPdf((unsigned)x, dist->poiss_.mean_));
  return prob;
}

EVAL_FN(Distribution::poissCdf)
{
  Probability prob;
  prob.setValue(Sampler::poissCdf((unsigned)x, dist->poiss_.mean_));
  return prob;
}

SAMP_FN(Distribution::poissSample)
{
  return Sampler::generatePoissonSample(dist->poiss_.mean_);
}

//-----------------------------------------------------------------------
// Chisq functions
//-----------------------------------------------------------------------

EVAL_FN(Distribution::chisqPdf)
{
  Probability prob;

  // For chisq distributions, the passed x-value is actually the
  // reduced chisq, so multiply by nDof before passing to Sampler

  prob.setLnValue(Sampler::lnChisqPdf(x * dist->chisq_.nDof_, dist->chisq_.nDof_));

  return prob;
}

EVAL_FN(Distribution::chisqCdf)
{
  Probability prob;

  // For chisq distributions, the passed x-value is actually the
  // reduced chisq, so multiply by nDof before passing to Sampler

  prob.setValue(Sampler::chisqCdf(x * dist->chisq_.nDof_, dist->chisq_.nDof_));
  return prob;
}

SAMP_FN(Distribution::chisqSample)
{
  ThrowError("chisqSample not defined");
  return 0.0;
}

//-----------------------------------------------------------------------
// Uniform functions
//-----------------------------------------------------------------------

EVAL_FN(Distribution::uniformPdf)
{
  Probability prob;
  prob.setLnValue(Sampler::lnUniformPdf(x, dist->uniform_.xMin_, dist->uniform_.xMax_));
  return prob;
}

EVAL_FN(Distribution::uniformCdf)
{
  Probability prob;
  prob.setValue(Sampler::uniformCdf(x, dist->uniform_.xMin_, dist->uniform_.xMax_));
  return prob;
}

SAMP_FN(Distribution::uniformSample)
{
  COUT("Generating uniform sample from " << dist->uniform_.xMin_ << " to " <<  dist->uniform_.xMax_);
  return Sampler::generateUniformSample(dist->uniform_.xMin_, dist->uniform_.xMax_);
}

//-----------------------------------------------------------------------
// Joint Gaussian functions
//-----------------------------------------------------------------------

EVAL_FN(Distribution::jointGaussPdf)
{
  Probability prob;

  // For chisq distributions, the passed x-value is actually the
  // reduced chisq, so multiply by nDof before setting

  prob.setLnValue(- 0.5 * x * dist->chisq_.nDof_);
  return prob;
}

EVAL_FN(Distribution::jointGaussCdf)
{
  Probability prob;
  ThrowError("jointGaussCdf not defined");
  return prob;
}

SAMP_FN(Distribution::jointGaussSample)
{
  ThrowError("jointGaussSample not defined");
  return 0.0;
}

ostream& 
gcp::util::operator<<(ostream& os, const Distribution& dist)
{
  operator<<(os, (Distribution&) dist);
}

ostream& 
gcp::util::operator<<(ostream& os, Distribution& dist)
{
  COUT("Distribution: type = " << dist.type_);

  if(dist.type_ == Distribution::DIST_UNIFORM) {
    COUT("min = " << dist.getUniformXMin() << " max = " << dist.getUniformXMax());
  }
}

void Distribution::invalidateTruncatedGaussianParameters()
{
  gauss_.norm_ = 0.0;
  gauss_.normIsCalculated_ = false;
}
