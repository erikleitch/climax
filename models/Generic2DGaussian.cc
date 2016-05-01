#include "gcp/models/Generic2DGaussian.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Generic2DGaussian::Generic2DGaussian() 
{
  addComponent(majSigma_);
  addComponent(axialRatio_);

  addComponentName(axialRatio_, "axialRatio", "The axial ratio (min/maj), for non-circularly symmetric models");
  addComponentName(majSigma_,   "majSigma",   "Width of the major axis");

  axialRatio_.allowUnitless(true);
  axialRatio_.setVal(1.0, "");

  initializeComponentsToFixed();
}

/**.......................................................................
 * Destructor.
 */
Generic2DGaussian::~Generic2DGaussian() {}

/**.......................................................................
 * Set the major axis sigma
 */
void Generic2DGaussian::setMajSigma(gcp::util::Angle majSigma)
{
  majSigma_.setVal(majSigma.arcmin(), "arcmin");
  majSigma_.wasSpecified_ = true;
}

/**.......................................................................
 * Specify the major axis as a FWHM
 */
void Generic2DGaussian::setMajFwhm(gcp::util::Angle majFwhm)
{
  setMajSigma(majFwhm/sqrt(8*log(2.0)));
}

/**.......................................................................
 * Set the axial ratio
 */
void Generic2DGaussian::setAxialRatio(double ratio)
{
  axialRatio_.setVal(ratio, "");
  axialRatio_.wasSpecified_ = true;
}

/**.......................................................................
 * Set both sigmas
 */
void Generic2DGaussian::setSigma(gcp::util::Angle sigma)
{
  setMajSigma(sigma);
  setAxialRatio(1.0);
}

/**.......................................................................
 * Set both FWHMs
 */
void Generic2DGaussian::setFwhm(gcp::util::Angle fwhm)
{
  setSigma(fwhm/sqrt(8*log(2.0)));
}

/**.......................................................................
 * Evaluate the dimensionless shape of this profile
 */
double Generic2DGaussian::envelope(unsigned type, double xRad, double yRad)
{
  double majSigRad = majSigma_.radians();
  double minSigRad = majSigRad * axialRatio_.val_;
  double arg = 0.5 * ((xRad * xRad)/(majSigRad * majSigRad) + (yRad * yRad)/(minSigRad * minSigRad));

  return exp(-arg);
}

/**.......................................................................
 * Evaluate the dimensionless shape of this profile
 */
double Generic2DGaussian::envelope(void* evalData, unsigned type, double xRad, double yRad)
{
  Generic2DGaussianEvalData* stack = (Generic2DGaussianEvalData*)evalData;
  stack->arg_ = 0.5 * ((xRad * xRad)/(stack->majSigRad_ * stack->majSigRad_) + (yRad * yRad)/(stack->minSigRad_ * stack->minSigRad_));
  return exp(-stack->arg_);
}

void* Generic2DGaussian::allocateEvalData()
{
  return new Generic2DGaussianEvalData();
}

void Generic2DGaussian::initializeEvalData(void* evalData)
{
  Generic2DGaussianEvalData* stack = (Generic2DGaussianEvalData*)evalData;

  stack->majSigRad_ = majSigma_.radians();
  stack->minSigRad_ = stack->majSigRad_ * axialRatio_.val_;

  return;
}

/**.......................................................................
 * Check this model's setup for sense
 */
void Generic2DGaussian::checkSetup()
{
  Generic2DAngularModel::checkSetup();
  checkVar("axialRatio");
  checkVar("majSigma");
}

PgModel Generic2DGaussian::pgModel()
{
  PgModel mod;

  mod.xMid_  = xOffset_.degrees();
  mod.yMid_  = yOffset_.degrees();
  mod.angle_ = rotationAngle_;
  mod.rot_   = rotationAngle_;
  mod.type_  = PgModel::TYPE_GAUSS;

  double rad1 = majSigma_.degrees();
  mod.xRad1_  = mod.xMid_ + rad1; 
  mod.yRad1_  = mod.yMid_;

  double rad2 = axialRatio_.val_ * majSigma_.degrees();
  mod.yRad2_  = mod.yMid_ + rad2;
  mod.xRad2_  = mod.xMid_;

  mod.rad1_ = rad1;
  mod.rad2_ = rad2;

  return mod;
}
