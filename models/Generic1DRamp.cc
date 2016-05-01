#include "gcp/models/Generic1DRamp.h"

using namespace std;

using namespace gcp::models;

/**.......................................................................
 * Constructor.
 */
Generic1DRamp::Generic1DRamp() 
{
  addComponent(off_);
  addComponent(slope_);

  addComponentName(off_,        "off",   "An offset (y-intercept)");
  addComponentName(slope_,      "slope", "The slope of the line");

  initializeComponentsToFixed();
}

void Generic1DRamp::setOffset(double off)
{
  off_ = off;
}

void Generic1DRamp::setSlope(double slope)
{
  slope_ = slope;
}

double Generic1DRamp::eval(double x)
{
  return off_.value() + slope_.value() * x;
}

/**.......................................................................
 * Destructor.
 */
Generic1DRamp::~Generic1DRamp() {}
