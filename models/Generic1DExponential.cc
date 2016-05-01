#include "gcp/models/Generic1DExponential.h"

#include "gcp/fftutil/DataSetType.h"

using namespace std;

using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Generic1DExponential::Generic1DExponential() 
{
  addComponent(peak_);
  addComponent(xPeak_);
  addComponent(lambdaRise_);
  addComponent(lambdaDecay_);

  addComponentName(peak_,        "peak",        "The amplitude of the exponential peak");
  addComponentName(xPeak_,       "xpeak",       "The x-value of the peak");
  addComponentName(lambdaRise_,  "lambdarise",  "The rise-time of the exponential envelope");
  addComponentName(lambdaDecay_, "lambdadecay", "The decay time of the exponential envelope");

  initializeComponentsToFixed();
}

/**.......................................................................
 * Destructor.
 */
Generic1DExponential::~Generic1DExponential() {}

void Generic1DExponential::setPeak(double peak)
{
  peak_ = peak;
}

void Generic1DExponential::setXPeak(double xPeak)
{
  xPeak_ = xPeak;
}

void Generic1DExponential::setLambdaRise(double lambda)
{
  lambdaRise_ = lambda;
}

void Generic1DExponential::setLambdaDecay(double lambda)
{
  lambdaDecay_ = lambda;
}

double Generic1DExponential::eval(double x)
{
  double dx = x - xPeak_.value();
  double lr = lambdaRise_.value();
  double ld = lambdaDecay_.value();

  if(dx < 0.0)
    return peak_.value() * exp(dx/lr);
  else
    return peak_.value() * exp(-dx/ld);
}

double Generic1DExponential::eval(void* evalData, double x)
{
  Generic1DExponentialEvalData* stack = (Generic1DExponentialEvalData*)evalData;

  stack->dx_ = x - stack->xPeak_;

  if(stack->dx_ < 0.0)
    return stack->peak_ * exp(stack->dx_/stack->lr_);
  else
    return stack->peak_ * exp(-stack->dx_/stack->ld_);
}

void* Generic1DExponential::allocateEvalData()
{
  return new Generic1DExponentialEvalData();
}

void Generic1DExponential::initializeEvalData(void* evalData)
{
  Generic1DExponentialEvalData* stack = (Generic1DExponentialEvalData*)evalData;
  stack->lr_    = lambdaRise_.value();
  stack->ld_    = lambdaDecay_.value();
  stack->peak_  = peak_.value();
  stack->xPeak_ = xPeak_.value();
}

void Generic1DExponential::checkSetup()
{
  checkVar("peak");
  checkVar("xpeak");
  checkVar("lambdarise");
  checkVar("lambdadecay");
}
