#include "gcp/models/GenericLine.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
GenericLine::GenericLine() 
{
  addComponent(width_);
  addComponent(length_);

  addComponentName(width_,  "width",  "The width of the line");
  addComponentName(length_, "length", "The length of the line");

  initializeComponentsToFixed();
}

/**.......................................................................
 * Destructor.
 */
GenericLine::~GenericLine() {}

/**.......................................................................
 * Set the length of the line
 */
void GenericLine::setLength(gcp::util::Angle len)
{
  length_.setVal(len.arcmin(), "arcmin");
  length_.wasSpecified_ = true;
}

/**.......................................................................
 * Set the width of the line
 */
void GenericLine::setWidth(gcp::util::Angle len)
{
  width_.setVal(len.arcmin(), "arcmin");
  width_.wasSpecified_ = true;
}

/**.......................................................................
 * Evaluate the dimensionless shape of this profile
 */
double GenericLine::envelope(unsigned type, double xRad, double yRad)
{
  COUT("Evaluating env with xRad = " << width_.radians() << " yRad = " << length_.radians());
  COUT("Evaluating env with xRad = " << xRad << " yRad = " << yRad);
  if(fabs(xRad) < width_.radians()/2 && fabs(yRad) < length_.radians()/2)
    return 1.0;
  else
    return 0.0;
}

/**.......................................................................
 * Evaluate the dimensionless shape of this profile
 */
double GenericLine::envelope(void* evalData, unsigned type, double xRad, double yRad)
{
  GenericLineEvalData* stack = (GenericLineEvalData*)evalData;

  if(fabs(xRad) < stack->widthRad_/2 && fabs(yRad) < stack->lengthRad_/2)
    return 1.0;
  else
    return 0.0;
}

void* GenericLine::allocateEvalData()
{
  return new GenericLineEvalData();
}

void GenericLine::initializeEvalData(void* evalData)
{
  GenericLineEvalData* stack = (GenericLineEvalData*)evalData;

  stack->widthRad_ = width_.radians();
  stack->lengthRad_ = length_.radians();

  return;
}

/**.......................................................................
 * Check this model's setup for sense
 */
void GenericLine::checkSetup()
{
  Generic2DAngularModel::checkSetup();
  checkVar("width");
  checkVar("length");
}

PgModel GenericLine::pgModel()
{
  PgModel mod;
  return mod;
}
