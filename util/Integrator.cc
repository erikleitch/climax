#include "gcp/util/Exception.h"
#include "gcp/util/Integrator.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Integrator::Integrator() 
{
  initializeGslMembers();
}

/**.......................................................................
 * Destructor.
 */
Integrator::~Integrator() 
{
  if(gslWork_) {
    gsl_integration_workspace_free(gslWork_);
    gslWork_ = 0;
  }
}

/**.......................................................................
 * Initialize resources needed for integration with the GSL library
 */
void Integrator::initializeGslMembers()
{
  gslLimit_ = 100;

  gslEpsAbs_ = 0;
  gslEpsRel_ = 1e-7;

  gslWork_ = 0;
  gslWork_ = gsl_integration_workspace_alloc(gslLimit_);

  if(!gslWork_) {
    ThrowError("Unable to allocate GSL workspace");
  }

  gslFn_.function   = 0;
  gslFn_.params     = 0;
}

void Integrator::installGslFn(INT_FN(*gslIntFn), void* params)
{
  gslFn_.function = gslIntFn;
  gslFn_.params   = params;
}

double Integrator::integrateFromLowlimToHighlim(INT_FN(*gslIntFn), void* params, double lowLim, double highLim)
{
  installGslFn(gslIntFn, params);
  return integrateFromLowlimToHighlim(lowLim, highLim);
}

double Integrator::integrateFromLowlimToHighlim(double lowLim, double highLim)
{
  double result, abserr;
  gsl_integration_qag(&gslFn_, lowLim, highLim, gslEpsAbs_, gslEpsRel_, gslLimit_, 6, gslWork_, &result, &abserr);
  return result;
}

double Integrator::integrateFromNegInftyToPosInfty(INT_FN(*gslIntFn), void* params)
{
  installGslFn(gslIntFn, params);
  return integrateFromNegInftyToPosInfty();
}

double Integrator::integrateFromNegInftyToPosInfty()
{
  double result, abserr;
  gsl_integration_qagi(&gslFn_, gslEpsAbs_, gslEpsRel_, gslLimit_, gslWork_, &result, &abserr);
  return result;
}

double Integrator::integrateFromLowlimToInfty(INT_FN(*gslIntFn), void* params, double lowLim)
{
  installGslFn(gslIntFn, params);
  return integrateFromLowlimToInfty(lowLim);
}

double Integrator::integrateFromLowlimToInfty(double lowLim)
{
  double result, abserr;
  gsl_integration_qagiu(&gslFn_, lowLim, gslEpsAbs_, gslEpsRel_, gslLimit_, gslWork_, &result, &abserr);
  return result;
}

double Integrator::integrateFromNegInftyToHighlim(INT_FN(*gslIntFn), void* params, double highLim)
{
  installGslFn(gslIntFn, params);
  return integrateFromNegInftyToHighlim(highLim);
}

double Integrator::integrateFromNegInftyToHighlim(double highLim)
{
  double result, abserr;
  gsl_integration_qagil(&gslFn_, highLim, gslEpsAbs_, gslEpsRel_, gslLimit_, gslWork_, &result, &abserr);
  return result;
}
