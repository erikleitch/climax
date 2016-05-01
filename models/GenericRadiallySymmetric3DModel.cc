#include "gcp/models/GenericRadiallySymmetric3DModel.h"
#include "gcp/models/CosmologyModel.h"

#include "gcp/util/Constants.h"
#include "gcp/util/Exception.h"

#include "gcp/pgutil/PgUtil.h"

#include "gsl/gsl_errno.h"

using namespace std;

using namespace gcp::models;
using namespace gcp::util;

#ifdef TIMER_TEST
#include "gcp/util/Timer.h"
Timer qit0, qit1, qit2, qit3;
double qt0=0.0, qt1=0.0, qt2=0.0, qt3=0.0;
#endif

/**.......................................................................
 * Constructor.
 */
GenericRadiallySymmetric3DModel::GenericRadiallySymmetric3DModel() 
{
  initialize();
}

/**.......................................................................
 * Constructor.
 */
void GenericRadiallySymmetric3DModel::initialize()
{
  //------------------------------------------------------------
  // Default to evaluating line integrals to x = 0. Inheritors with
  // non-integrable kernels can redefine this
  //------------------------------------------------------------

  xLim_ = 0.0;

  setGslLimit(100);
  initializeGslMembers();

  //------------------------------------------------------------
  // Default these to sensible values -- inheritors are free to
  // redefine
  //------------------------------------------------------------

  newSample_[DataSetType::DATASET_RADIO] = true;
  newSample_[DataSetType::DATASET_XRAY_IMAGE] = true;
  needsRecomputing_[DataSetType::DATASET_RADIO] = true;
  needsRecomputing_[DataSetType::DATASET_XRAY_IMAGE] = true;

  //------------------------------------------------------------
  // Add prerequisites for derived variate computation
  //------------------------------------------------------------

  scaleFactorY_.allowUnitless(true);
  addPrerequisite(scaleFactorY_,  "");
  addComponentName(scaleFactorY_, "scaleFactorY");
  scaleFactorY_.deriveWith(GenericRadiallySymmetric3DModel::deriveScaleFactorY, this);
  scaleFactorY_.dependsOn(normalizationFrequency_);

  addPrerequisite(innerRadius_,  "Mpc");
  addComponentName(innerRadius_, "radInner");
  innerRadius_.deriveWith(GenericRadiallySymmetric3DModel::deriveInnerRadius, this);

  addPrerequisite(outerRadius_,  "Mpc");
  addComponentName(outerRadius_, "radOuter");
  outerRadius_.deriveWith(GenericRadiallySymmetric3DModel::deriveOuterRadius, this);

  restMassThermalEnergyRatio_.allowUnitless(true);
  addPrerequisite(restMassThermalEnergyRatio_,  "");
  addComponentName(restMassThermalEnergyRatio_,  "restMassThermalEnergyRatio");
  restMassThermalEnergyRatio_.deriveWith(GenericRadiallySymmetric3DModel::deriveRestMassThermalEnergyRatio, this);

  //============================================================
  // Define derived variates
  //============================================================

  addDerivedComponent(scaleFactorArea_, "Mpc^2", false);
  addComponentName(scaleFactorArea_,    "scaleFactorArea");
  scaleFactorArea_.deriveWith(GenericRadiallySymmetric3DModel::deriveScaleFactorArea, this);

  addDerivedComponent(scaleFactorMass_, "Msolar", false);
  addComponentName(scaleFactorMass_,    "scaleFactorMass");
  scaleFactorMass_.deriveWith(GenericRadiallySymmetric3DModel::deriveScaleFactorMass, this);
  scaleFactorMass_.dependsOn(restMassThermalEnergyRatio_);
  scaleFactorMass_.dependsOn(scaleFactorArea_);
  scaleFactorMass_.dependsOn(mue_);

  addDerivedComponent(scaleFactorPressure_, "keV/cm^3", false);
  addComponentName(scaleFactorPressure_,    "scaleFactorPressure");
  scaleFactorPressure_.deriveWith(GenericRadiallySymmetric3DModel::deriveScaleFactorPressure, this);

  addDerivedComponent(thetaInner_, "arcmin", false);
  addComponentName(thetaInner_,    "thetaInner");
  thetaInner_.deriveWith(GenericRadiallySymmetric3DModel::deriveThetaInner, this);
  thetaInner_.dependsOn(innerRadius_);

  addDerivedComponent(thetaOuter_, "arcmin", false);
  addComponentName(thetaOuter_,    "thetaOuter");
  thetaOuter_.deriveWith(GenericRadiallySymmetric3DModel::deriveThetaOuter, this);
  thetaOuter_.dependsOn(outerRadius_);

  //-----------------------------------------------------------
  // Line integral isn't really a derived variate, but it does depend
  // on thetaCore_ being set to a sensible value (ie, non-negative).
  // Since thetaCore_ can be derived, we add the line integral
  // normalization as a derived variate so that its dependency on
  // thetaCore will be taken account of correctly
  //------------------------------------------------------------

  lineIntegralNormalization_[DataSetType::DATASET_RADIO].allowUnitless(true);
  addDerivedComponent(lineIntegralNormalization_[DataSetType::DATASET_RADIO], "", false);
  addComponentName(lineIntegralNormalization_[DataSetType::DATASET_RADIO],    "lin");
  lineIntegralNormalization_[DataSetType::DATASET_RADIO].deriveWith(GenericRadiallySymmetric3DModel::deriveRadioLineIntegralNormalization, this);
  lineIntegralNormalization_[DataSetType::DATASET_RADIO].dependsOn(thetaCore_);

  //------------------------------------------------------------
  // Likewise volumeIntegralFactor isn't a true derived variate, but
  // we want the ability to insert it into the dependency tree
  //------------------------------------------------------------

  volumeIntegralFactor_.allowUnitless(true);
  addDerivedComponent(volumeIntegralFactor_, "", false);
  addComponentName(volumeIntegralFactor_,    "vif");
  volumeIntegralFactor_.deriveWith(GenericRadiallySymmetric3DModel::deriveVolumeIntegralFactor, this);
  volumeIntegralFactor_.dependsOn(thetaCore_);
  volumeIntegralFactor_.dependsOn(thetaInner_);
  volumeIntegralFactor_.dependsOn(thetaOuter_);

  //------------------------------------------------------------
  // Add Pe0 and set default units to erg/cm^3
  //------------------------------------------------------------

  addDerivedComponent(pe0_,  "erg/cm^3");
  addComponentName(pe0_,     "Pe0", "A pressure normalization for this model (derived)");
  pe0_.deriveWith(GenericRadiallySymmetric3DModel::derivePe0, this);
  pe0_.dependsOn(lineIntegralNormalization_[DataSetType::DATASET_RADIO]);
  pe0_.dependsOn(radioNormalization_);
  pe0_.dependsOn(thetaCore_);
  pe0_.dependsOn(scaleFactorY_);
  pe0_.dependsOn(scaleFactorPressure_);

  //------------------------------------------------------------
  // Add Ysph and set default units to Mpc^2
  //------------------------------------------------------------

  addDerivedComponent(ySph_, "Mpc^2");
  addComponentName(ySph_,    "Ysph", "Compton Y, spherically integrated between the specified innerRadius and outerRadius (derived)");
  ySph_.deriveWith(GenericRadiallySymmetric3DModel::deriveYsph, this);
  ySph_.dependsOn(lineIntegralNormalization_[DataSetType::DATASET_RADIO]);
  ySph_.dependsOn(volumeIntegralFactor_);
  ySph_.dependsOn(radioNormalization_);
  ySph_.dependsOn(thetaCore_);
  ySph_.dependsOn(scaleFactorY_);
  ySph_.dependsOn(scaleFactorArea_);

  //------------------------------------------------------------
  // Add Mgas and default units to Msolar
  //------------------------------------------------------------

  addDerivedComponent(mGas_, "Msolar");
  addComponentName(mGas_,    "Mgas", "Gas mass, spherically integrated between the specified innerRadius and outerRadius (derived)");
  mGas_.deriveWith(GenericRadiallySymmetric3DModel::deriveMgas, this);
  mGas_.dependsOn(lineIntegralNormalization_[DataSetType::DATASET_RADIO]);
  mGas_.dependsOn(volumeIntegralFactor_);
  mGas_.dependsOn(radioNormalization_);
  mGas_.dependsOn(thetaCore_);
  mGas_.dependsOn(scaleFactorY_);
  mGas_.dependsOn(scaleFactorMass_);

  //------------------------------------------------------------
  // Add Mtot and default units to Msolar
  //------------------------------------------------------------

  addDerivedComponent(mTot_, "Msolar");
  addComponentName(mTot_,    "Mtot", "Total mass, spherically integrated between the specified innerRadius and outerRadius (derived)");
  mTot_.deriveWith(GenericRadiallySymmetric3DModel::deriveMtot, this);
  mTot_.dependsOn(lineIntegralNormalization_[DataSetType::DATASET_RADIO]);
  mTot_.dependsOn(volumeIntegralFactor_);
  mTot_.dependsOn(radioNormalization_);
  mTot_.dependsOn(thetaCore_);
  mTot_.dependsOn(fGas_);
  mTot_.dependsOn(scaleFactorY_);
  mTot_.dependsOn(scaleFactorMass_);

  addParameter("thetaCoreMin",         DataType::DOUBLE, "The minimum thetaCore for which line integrals will be precomputed");
  addParameter("thetaCoreMax",         DataType::DOUBLE, "The maximum thetaCore for which line integrals will be precomputed");

  //------------------------------------------------------------
  // And initialize all components to fixed
  //------------------------------------------------------------

  initializeComponentsToFixed();
}

/**.......................................................................
 * Initialize limits needed for integration with the GSL library
 */
void GenericRadiallySymmetric3DModel::setGslLimit(size_t limit)
{
  gslLimit_ = limit;
}

/**.......................................................................
 * Initialize resources needed for integration with the GSL library
 */
void GenericRadiallySymmetric3DModel::initializeGslMembers()
{
  gslWork_ = 0;
  gslWork_ = gsl_integration_workspace_alloc(gslLimit_);

  if(!gslWork_) {
    ThrowError("Unable to allocate GSL workspace");
  }

  gslLineIntegralFn_.function   = &GenericRadiallySymmetric3DModel::evaluateLineIntegralKernel2;
  gslLineIntegralFn_.params     = (void*)this;

  gslVolumeIntegralFn_.function = &GenericRadiallySymmetric3DModel::evaluateVolumeIntegralKernel;
  gslVolumeIntegralFn_.params   = (void*)this;

  gsl_set_error_handler(&intErrHandler);
}

/**.......................................................................
 * Destructor.
 */
GenericRadiallySymmetric3DModel::~GenericRadiallySymmetric3DModel() 
{
  if(gslWork_) {
    gsl_integration_workspace_free(gslWork_);
    gslWork_ = 0;
  }
}

/**.......................................................................
 * Single-parameter radial model to evaluate.
 */
double GenericRadiallySymmetric3DModel::radialModel(unsigned type, double x, void* params)
{
  ThrowError("Inheritor has not defined an inherited radialModel() method");
  return 0.0;
}

/**.......................................................................
 * Evaluate the single-parameter radial model at the requested
 * xSky-coordinate.  If P(x) is defined for x, where x = r/r_c, or
 * equivalently x = theta/theta_c, then xSky = theta_sky/theta_c, ie,
 * it is the dimensionless coordinate given by the ratio of the angle
 * on the sky to theta_c.
 *                                             / +infty
 * We define the line-integral of P(x) to be: |    P(x(l)) dl
 *                                            /  -infty
 *
 * If R is the cylindrical radius (physical length) corresponding to
 * an angle on the sky, then each x value along the line L is given
 * by:
 *
 *  x(l) = r/r_c = sqrt(l^2 + R^2)/r_c 
 *               = sqrt((D_A*theta_l)^2 + (D_A * theta_sky)^2)/(D_A * theta_c)
 *               = sqrt(theta_l^2 + theta_sky^2)/theta_c
 *               = sqrt(xl^2 + xSky^2)
 *    
 * And the line integral becomes:
 *
 *            / +infty
 * I(xSky) = |     P(sqrt(xl^2 + xSky^2)) dl
 *           /  -infty
 *
 * We can then make the transformation: theta_l = l/D_A, 
 * or dl = D_A * dtheta_l to obtain:
 * 
 *                / +infty
 * I(xSky) = D_A |  P(sqrt(xl^2 + xSky^2)) dtheta_l
 *               /  -infty
 * 
 * evaluateLineIntegralKernel(xl) is called to evaluate the kernel of
 * this integral at a fixed value of xl, so we calculate the x that
 * corresponds to this, per the above, and return P(x).
 *
 * In practice, this class normalizes the integrated radial model to 1
 * at xSky = 0.  
 *
 * For pressure models, normalization of the model in units of Compton
 * y means:
 *
 *           sigma_T		   
 *  y0 =  ------------- * D_A * P0 
 *          m_e * c^2              
 *
 * For density models, we have:
 *
 *        k_B * sigma_T		   
 *  y0 =  ------------- * D_A * n_e0 * T_e 
 *	    m_e * c^2              
 *
 */
double GenericRadiallySymmetric3DModel::evaluateLineIntegralKernel(double xl, void* params)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)params;

  //------------------------------------------------------------
  // Convert from line-of-sight xl to 3D x
  //------------------------------------------------------------

  double x = sqrt(xl*xl + model->xSky_ * model->xSky_);

  //------------------------------------------------------------
  // And evaluate the radial model at this x
  //------------------------------------------------------------

  return model->radialModel(model->currentIntegrationType_, x, model->params_);
}

double GenericRadiallySymmetric3DModel::evaluateLineIntegralKernel2(double x, void* params)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)params;

  //------------------------------------------------------------
  // Convert from line-of-sight xl to 3D x
  //------------------------------------------------------------

  double denom = sqrt(x * x - model->xSky_ * model->xSky_);

  //------------------------------------------------------------
  // And evaluate the radial model at this x
  //------------------------------------------------------------

  //  COUT("Evaluating kernel at x = " << x << " xSky_ = " << model->xSky_);

  return model->radialModel(model->currentIntegrationType_, x, model->params_) * x / denom;
}

/**.......................................................................
 * Evaluate the line integral of the 3D radial model at the specified
 * angular radius (xSky).  Integral is evaluated from -infinity to
 * +infinity using GSL adaptive integration.
 */
double GenericRadiallySymmetric3DModel::lineIntegral(unsigned type, double xSky, void* params)
{
  double result, abserr;

  params_ = params;
  xSky_   = xSky;
  currentIntegrationType_ = type;

  //------------------------------------------------------------
  // Just integrate the model from -infty to +infty: i.e., we assume
  // the center is effectively infinitely away from us.  In practice,
  // we integrate from 0 to +infty, and double the result
  //------------------------------------------------------------

  gsl_integration_qagiu(&gslLineIntegralFn_, xSky_, 0, 1e-7, gslLimit_, gslWork_, &result, &abserr);

  return 2*result;
}

/**.......................................................................
 * 
 *                                              / x(R)
 * We define the volume-integral of P(x) to be: |   P(r) 4pi r^2 dr
 *                                              / 0
 *
 *
 * Here we evaluate the kernel of this integral (P(r) 4pi r^2) 
 */
double GenericRadiallySymmetric3DModel::evaluateVolumeIntegralKernel(double x, void* params)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)params;

  //------------------------------------------------------------
  // And evaluate the radial model at this x
  //------------------------------------------------------------

  return 4 * M_PI * model->radialModel(model->currentIntegrationType_, x, model->params_) * x * x;
}

/**.......................................................................
 * Evaluate the volume integral of the 3D radial model at the
 * specified angular radius.  Integral is evaluated from 0 to xSky
 * using GSL adaptive integration.
 */
double GenericRadiallySymmetric3DModel::volumeIntegral(unsigned type, double xSky, void* params)
{
  double result, abserr;

  params_ = params;
  currentIntegrationType_ = type;

  //------------------------------------------------------------
  // Just integrate the model from 0 to xSky = theta/theta_c
  //------------------------------------------------------------

  gsl_integration_qag(&gslVolumeIntegralFn_, 1e-9, xSky, 0, 1e-7, gslLimit_, 6, gslWork_, &result, &abserr);

  return result;
}

/**.......................................................................
 * Return the value of the envelope at the specified angular point on
 * the sky.
 */
double GenericRadiallySymmetric3DModel::radioEnvelope(double xRad, double yRad)
{
  //------------------------------------------------------------
  // Get the dimensionless x-coordinate on the sky at which we want
  // the line integral
  //------------------------------------------------------------

#ifdef TIMER_TEST
  qit0.start();
#endif

  double xSky = sqrt(xRad * xRad + yRad * yRad) / thetaCore_.radians();

#ifdef TIMER_TEST
  qit0.stop();
  qt0 += qit0.deltaInSeconds();
#endif

  //------------------------------------------------------------
  // And interpolate
  //------------------------------------------------------------

  return interpolateLineIntegral(DataSetType::DATASET_RADIO, xSky);
}

/**.......................................................................
 * Return the value of the envelope at the specified angular point on
 * the sky.
 */
double GenericRadiallySymmetric3DModel::xrayImageEnvelope(double xRad, double yRad)
{
  //------------------------------------------------------------
  // Get the dimensionless x-coordinate on the sky at which we want
  // the line integral
  //------------------------------------------------------------

  double xSky = sqrt(xRad * xRad + yRad * yRad) / thetaCore_.radians();

  //------------------------------------------------------------
  // And interpolate
  //------------------------------------------------------------

  double retval = interpolateLineIntegral(DataSetType::DATASET_XRAY_IMAGE, xSky);

  return retval;
}

/**.......................................................................
 * Return the value of the envelope at the specified angular point on
 * the sky.
 */
double GenericRadiallySymmetric3DModel::genericEnvelope(double xRad, double yRad)
{
  //------------------------------------------------------------
  // Get the dimensionless x-coordinate on the sky at which we want
  // the line integral
  //------------------------------------------------------------

  double xSky = sqrt(xRad * xRad + yRad * yRad) / thetaCore_.radians();

  //------------------------------------------------------------
  // And interpolate
  //------------------------------------------------------------

  return interpolateLineIntegral(DataSetType::DATASET_GENERIC, xSky);
}

/**.......................................................................
 * Overload the base-class method to set a flag whenever a new sample
 * has just been generated.
 */
void GenericRadiallySymmetric3DModel::sample()
{
  Generic2DAngularModel::sample();

  for(std::map<unsigned, bool>::iterator iter=newSample_.begin(); iter != newSample_.end(); iter++)
    iter->second = true;
}

/**.......................................................................
 * Overload the base-class method to calculate the line integral of
 * the radial model if this is the first time we've been called since
 * a new sample has been generated.
 */
void GenericRadiallySymmetric3DModel::
fillImage(unsigned type, Image& image, void* params)
{
  //------------------------------------------------------------
  // Select the specific type of this dataset (ie, we only want to
  // know if it is radio or xray for example, not whether it is a 1D
  // or 2D dataset, which we already know at this point)
  //------------------------------------------------------------

  type = type & (DataSetType::DATASET_GENERIC | DataSetType::DATASET_RADIO | DataSetType::DATASET_XRAY_IMAGE);

  //------------------------------------------------------------
  // Recalculate interpolation values if needed prior to filling the
  // image
  //------------------------------------------------------------

  if(newSample_[type] && needsRecomputing_[type]) {
    calculateInterpolationValues(type, image);
    newSample_[type] = false;
  }

  //------------------------------------------------------------
  // Now call down to the base-class to fill the image
  //------------------------------------------------------------

  Generic2DAngularModel::fillImage(type, image, params);
}

/**.......................................................................
 * Calculate the line integral of our model on a grid fine enough to
 * interpolate anywhere in the passed image.
 */
void GenericRadiallySymmetric3DModel::calculateInterpolationValues(unsigned type, Image& image)
{
  //------------------------------------------------------------
  // If shape parameters are fixed, we can calculate the line integral
  // of the model on a grid fine enough for all values of the core
  // radius once, and we don't have to calculate them again.
  //------------------------------------------------------------

  if(shapeParametersAreFixed()) {

    calculateInterpolationValuesForAllScaleRadii(type, image);
    needsRecomputing_[type] = false;

    //------------------------------------------------------------
    // Otherwise we will in general have to recalculate the
    // interpolation values each time a new sample of the shape
    // parameters is drawn, in which case we only compute the line
    // integral on a grid appropriate for the current value of the
    // core radius
    //------------------------------------------------------------

  } else {
    calculateInterpolationValuesForCurrentScaleRadius(type, image);
    needsRecomputing_[type] = true;
  }
}

/**.......................................................................
 * If the shape parameters can change from sample to sample, we
 * calculate interpolation values for the current value of the scale
 * radius only, since that requires fewer evaluations of the line
 * integral.
 */
void GenericRadiallySymmetric3DModel::
calculateInterpolationValuesForCurrentScaleRadius(unsigned type, Image& image)
{
  //------------------------------------------------------------
  // Make the interpolation array large enough to cover the corners of
  // the image
  //------------------------------------------------------------
  
  nInterp_ = sqrt(2.0) * image.xAxis().getNpix();

  //------------------------------------------------------------
  // Compute at least three points, so that we can perform a quadratic
  // interpolation
  //------------------------------------------------------------

  if(nInterp_ < 3)
    nInterp_ = 3;
  
  //------------------------------------------------------------
  // Get the delta-r as a dimensionless ratio to the scale radius
  //------------------------------------------------------------

  deltaInterp_ = image.xAxis().getAngularResolution().radians() / thetaCore_.radians();
  
  calculateInterpolationValues(type);
}

/**.......................................................................
 * If the shape parameters cannot change from sample to sample, we
 * calculate interpolation values suitable for interpolating all
 * possible values of the scale radius, even though it is inefficient,
 * since we can do this once for all time and never again.
 */
void GenericRadiallySymmetric3DModel::
calculateInterpolationValuesForAllScaleRadii(unsigned type, Image& image)
{
  //------------------------------------------------------------
  // Make the interpolation array large enough to extend out to the
  // corners of the image, at the lowest resolution of the model that
  // will ever be required.  This xMax will be sampled at the highest
  // resolution that will ever be required.  
  // 
  // Note that this calculation assumes a uniform prior has been
  // specified for thetaCore_.  If not, then we will have to come up
  // with a heuristic that accounts for other types of priors. (+- N
  // sigma for Gaussian, etc).
  //------------------------------------------------------------

  double thetaCoreMinRad;
  double thetaCoreMaxRad;

  if(thetaCore_.isVariable()) {

    //------------------------------------------------------------
    // If a uniform prior has been assigned, then we know a priori the
    // range of values we have to consider
    //------------------------------------------------------------

    if(thetaCore_.prior().getType() == Distribution::DIST_UNIFORM) {

      thetaCoreMinRad = thetaCore_.getVal(thetaCore_.prior().getUniformXMin(), "radians");
      thetaCoreMaxRad = thetaCore_.getVal(thetaCore_.prior().getUniformXMax(), "radians");

      //------------------------------------------------------------
      // If a gaussian prior has been assigned, we consider the range
      // of values likely to be encountered in a chain of this length
      //------------------------------------------------------------

    } else if(thetaCore_.prior().getType() == Distribution::DIST_GAUSS) {

      unsigned nSigma = getMaxNSigma();

      double thetaMin = 
	thetaCore_.getVal(thetaCore_.prior().getGaussMean() - nSigma * thetaCore_.prior().getGaussSigma(), "radians");

      double thetaMax = 
	thetaCore_.getVal(thetaCore_.prior().getGaussMean() + nSigma * thetaCore_.prior().getGaussSigma(), "radians");

      //------------------------------------------------------------
      // Limit the minimum radius to sigma/10 if the range extends
      // below zero
      //------------------------------------------------------------

      thetaCoreMinRad = thetaMin > 0.0 ? thetaMin : thetaCore_.prior().getGaussSigma()/10;
      thetaCoreMaxRad = thetaMax;

      //------------------------------------------------------------
      // If a truncated gaussian prior has been assigned, we consider
      // the truncation range, if finite, else the range of values
      // likely to be encountered in a chain of this length
      //------------------------------------------------------------

    } else if(thetaCore_.prior().getType() == Distribution::DIST_TRUNC_GAUSS) {

      unsigned nSigma = getMaxNSigma();

      double thetaMin = isfinite(thetaCore_.prior().getUniformXMin()) ? thetaCore_.prior().getUniformXMin() :
	thetaCore_.getVal(thetaCore_.prior().getGaussMean() - nSigma * thetaCore_.prior().getGaussSigma(), "radians");
      
      double thetaMax = isfinite(thetaCore_.prior().getUniformXMax()) ? thetaCore_.prior().getUniformXMax() :
	thetaCore_.getVal(thetaCore_.prior().getGaussMean() + nSigma * thetaCore_.prior().getGaussSigma(), "radians");
      
      //------------------------------------------------------------
      // Limit the minimum radius to sigma/10 if the range extends
      // below zero
      //------------------------------------------------------------
      
      thetaCoreMinRad = thetaMin > 0.0 ? thetaMin : thetaCore_.prior().getGaussSigma()/10;
      thetaCoreMaxRad = thetaMax;
      
    } else {
      ThrowError("No algorithm has been specified for how to compute limits for distribution of type " 
		 << thetaCore_.prior().getType());
    }

    //------------------------------------------------------------
    // If thetaCoreMin or thetaCoreMax were specified, overwrite the
    // value we just calculated
    //------------------------------------------------------------

    if(getParameter("thetaCoreMin", false)->data_.hasValue()) {
      Angle tMin;
      tMin.setVal(getDoubleVal("thetaCoreMin"), getParameter("thetaCoreMin", true)->units_);
      thetaCoreMinRad = tMin.radians();
    }

    if(getParameter("thetaCoreMax", false)->data_.hasValue()) {
      Angle tMax;
      tMax.setVal(getDoubleVal("thetaCoreMax"), getParameter("thetaCoreMax", true)->units_);
      thetaCoreMaxRad = tMax.radians();
    }

  } else {
    thetaCoreMinRad = thetaCore_.getVal(thetaCore_.value(), "radians");
    thetaCoreMaxRad = thetaCore_.getVal(thetaCore_.value(), "radians");
  }

  //------------------------------------------------------------
  // Get the length of the diagonal of the image, in dimensionless 
  // x = r/r_c, at the smallest value (lowest resolution) that the 
  // prior can assume:
  //------------------------------------------------------------

  double xMax = (sqrt(2.0) * image.xAxis().getAngularSize().radians()) / thetaCoreMinRad;

  //------------------------------------------------------------
  // Set deltaInterp_ to the per-pixel increment in dimensionless 
  // x = r/r_c, at the largest value (highest resolution) that the 
  // prior can assume:
  //------------------------------------------------------------

  deltaInterp_ = image.xAxis().getAngularResolution().radians() / thetaCoreMaxRad;

  //------------------------------------------------------------
  // Compute at least three points, so that we can perform a quadratic
  // interpolation
  //------------------------------------------------------------

  nInterp_ = (unsigned)ceil(xMax/deltaInterp_);

  if(nInterp_ < 3)
    nInterp_ = 3;

  calculateInterpolationValues(type);
}

/**.......................................................................
 * Use the current nInterp_ and deltaInterp_ to compute interpolation
 * values for the passed type
 */
void GenericRadiallySymmetric3DModel::
calculateInterpolationValues(unsigned type)
{
  if(nInterp_ != interpolatedLineIntegral_[type].size()) {
    interpolatedLineIntegral_[type].resize(nInterp_);
    arealIntegral_[type].resize(nInterp_);
  }
  
  double x;

  lineIntegralNormalization_[type].val_ = lineIntegral(type, xLim_);

  //  COUT("nInterp_ = " << nInterp_ << " deltaInterp_ = " << deltaInterp_ << " norm = " << lineIntegralNormalization_[type].val_ << " AV = " << lineIntegral(type, 0.03));

  for(unsigned i=0; i < nInterp_; i++) {
    x = i * deltaInterp_;

    //------------------------------------------------------------
    // Get the line integral at the specified point
    //------------------------------------------------------------

    if(x <= xLim_) {
      interpolatedLineIntegral_[type][i] = 1.0;
    } else {
      interpolatedLineIntegral_[type][i] = lineIntegral(type, x) / lineIntegralNormalization_[type].val_;
    }

    //------------------------------------------------------------
    // Since the profiles for this class by definition don't depend on
    // angle (but only radius), we can calculate the areal integration
    // of this function at the same time by performing a simple line
    // integral of the function we've just evaluated.
    //
    // Given f(r), we want the integral:
    //
    //                / x
    //          1     |
    // F(x) = ------  |  f(r) 2pi*r dr
    //        pi*x^2  |
    //                / 0
    //
    //------------------------------------------------------------

    if(i == 0) {

      //------------------------------------------------------------
      //  Integral at zero radius is zero
      //------------------------------------------------------------

      arealIntegral_[type][i] = 0.0;

    } else {

      //------------------------------------------------------------
      // Take the value of the function at the midpoint
      //------------------------------------------------------------

      unsigned i2 = i;
      unsigned i1 = i-1;
      x -= deltaInterp_/2;
      
      //------------------------------------------------------------
      // For the current value of i, we take f(r) * 2*pi*r * dr, with r a
      // the midpoint of the last two samples
      //------------------------------------------------------------

      double val = (interpolatedLineIntegral_[type][i2] + interpolatedLineIntegral_[type][i1])/2 * 2 * M_PI * x * deltaInterp_;
      arealIntegral_[type][i] = val + arealIntegral_[type][i-1];
    }
  }

  //------------------------------------------------------------
  // Now we normalize the integral by the area at each radius
  //------------------------------------------------------------

  for(unsigned i=0; i < nInterp_; i++) {
    x = i * deltaInterp_;
    arealIntegral_[type][i] = arealIntegral_[type][i]/(M_PI*x*x);
  }
}

/**.......................................................................
 * Use the current nInterp_ and deltaInterp_ to compute interpolation
 * values for the passed type
 */
void GenericRadiallySymmetric3DModel::
calculateAllInterpolationValues(unsigned type)
{
  if(nInterp_ != interpolatedLineIntegral_[type].size()) {
    interpolatedLineIntegral_[type].resize(nInterp_);
    arealIntegral_[type].resize(nInterp_);
    volumeIntegral_[type].resize(nInterp_);
  }
  
  double x;
  lineIntegralNormalization_[type].val_ = lineIntegral(type, xLim_);

  for(unsigned i=0; i < nInterp_; i++) {
    x = i * deltaInterp_;

    //------------------------------------------------------------
    // Get the line and volume integrals at the specified point
    //------------------------------------------------------------

    if(x <= xLim_) {
      interpolatedLineIntegral_[type][i] = 1.0;
    } else {
      interpolatedLineIntegral_[type][i] = lineIntegral(type, x) / lineIntegralNormalization_[type].val_;
    }

    volumeIntegral_[type][i]           = volumeIntegral(type, x);

    //------------------------------------------------------------
    // Since the profiles for this class by definition don't depend on
    // angle (but only radius), we can calculate the areal integration
    // of this function at the same time by performing a simple line
    // integral of the function we've just evaluated.  
    //
    // Given f(r), we want the integral:
    //
    //          1     /x
    // F(x) = ------  |  f(r) 2pi*r dr
    //        pi*x^2  /0
    //
    //------------------------------------------------------------

    if(i == 0) {

      //------------------------------------------------------------
      //  Integral at zero radius is zero
      //------------------------------------------------------------

      arealIntegral_[type][i] = 0.0;

    } else {

      //------------------------------------------------------------
      // Take the value of the function at the midpoint
      //------------------------------------------------------------

      unsigned i2 = i;
      unsigned i1 = i-1;
      x -= deltaInterp_/2;
      
      //------------------------------------------------------------
      // For the current value of i, we take f(r) * 2*pi*r * dr, with r a
      // the midpoint of the last two samples
      //------------------------------------------------------------

      double val = (interpolatedLineIntegral_[type][i2] + interpolatedLineIntegral_[type][i1])/2 * 2 * M_PI * x * deltaInterp_;
      arealIntegral_[type][i] = val + arealIntegral_[type][i-1];
    }
  }

  //------------------------------------------------------------
  // Now we normalize the integral by the area at each radius
  //------------------------------------------------------------

  for(unsigned i=0; i < nInterp_; i++) {
    x = i * deltaInterp_;
    arealIntegral_[type][i] = arealIntegral_[type][i]/(M_PI*x*x);
  }
}

/**.......................................................................
 * Use the current nInterp_ and deltaInterp_ to compute interpolation
 * values for the passed type
 */
void GenericRadiallySymmetric3DModel::
calculateVolumeInterpolationValues(unsigned type)
{
  if(nInterp_ != volumeIntegral_[type].size()) {
    volumeIntegral_[type].resize(nInterp_);
  }
  
  double x;

  for(unsigned i=0; i < nInterp_; i++) {
    x = i * deltaInterp_;

    //------------------------------------------------------------
    // Get the line and volume integrals at the specified point
    //------------------------------------------------------------

    volumeIntegral_[type][i]           = volumeIntegral(type, x);
  }
}

/**.......................................................................
 * Perform a quadratic interpolation of the value of the line integral at
 * the requested x-value
 */
double GenericRadiallySymmetric3DModel::interpolateLineIntegral(unsigned type, double xSky)
{
#ifdef TIMER_TEST
  qit1.start();
#endif

  //  unsigned iNear = (unsigned)floor(xSky / deltaInterp_);
  unsigned iNear = (unsigned)(xSky / deltaInterp_);
  unsigned i1, i2, i3;

  //------------------------------------------------------------
  // Protect against values off the ends of our array
  //------------------------------------------------------------

#ifdef TIMER_TEST
  qit1.stop();
  qt1 += qit1.deltaInSeconds();
  qit2.start();
#endif

  if(iNear >= nInterp_-1) {
    i1 = nInterp_ - 3;
    i2 = nInterp_ - 2;
    i3 = nInterp_ - 1;
  } else if(iNear <= 0) {
    i1 = 0;
    i2 = 1;
    i3 = 2;
  } else {
    i1 = iNear-1;
    i2 = iNear;
    i3 = iNear+1;
  }

#ifdef TIMER_TEST
  qit2.stop();
  qt2 += qit2.deltaInSeconds();
  qit3.start();
#endif

  double ret = quadInterp_.directEval(i1 * deltaInterp_, interpolatedLineIntegral_[type][i1],
				      i2 * deltaInterp_, interpolatedLineIntegral_[type][i2],
				      i3 * deltaInterp_, interpolatedLineIntegral_[type][i3], 
				      xSky);

#ifdef TIMER_TEST
  qit3.stop();
  qt3 += qit3.deltaInSeconds();
#endif

  return ret;
}

/**.......................................................................
 * Perform a quadratic interpolation of the areal integration of the
 * profile function at the requested x-value
 */
double GenericRadiallySymmetric3DModel::interpolateArealIntegral(unsigned type, double xSky)
{
  unsigned iNear = (unsigned)(xSky / deltaInterp_);
  unsigned i1, i2, i3;

  //------------------------------------------------------------
  // Protect against values off the ends of our array
  //------------------------------------------------------------

  if(iNear >= nInterp_-1) {
    i1 = nInterp_ - 3;
    i2 = nInterp_ - 2;
    i3 = nInterp_ - 1;
  } else if(iNear <= 0) {
    i1 = 0;
    i2 = 1;
    i3 = 2;
  } else {
    i1 = iNear-1;
    i2 = iNear;
    i3 = iNear+1;
  }

  if(i1 == 0)
    ThrowError("Request for integral at radius " << xSky << " which is smaller than can be evaluated");

  double ret = quadInterp_.directEval(i1 * deltaInterp_, arealIntegral_[type][i1],
				      i2 * deltaInterp_, arealIntegral_[type][i2],
				      i3 * deltaInterp_, arealIntegral_[type][i3], 
				      xSky);

  return ret;
}

/**.......................................................................
 * Perform a quadratic interpolation of the volume integration of the
 * profile function at the requested x-value
 */
double GenericRadiallySymmetric3DModel::interpolateVolumeIntegral(unsigned type, double xSky)
{
  unsigned iNear = (unsigned)(xSky / deltaInterp_);
  unsigned i1, i2, i3;

  //------------------------------------------------------------
  // Protect against values off the ends of our array
  //------------------------------------------------------------

  if(iNear >= nInterp_-1) {
    i1 = nInterp_ - 3;
    i2 = nInterp_ - 2;
    i3 = nInterp_ - 1;
  } else if(iNear <= 0) {
    i1 = 0;
    i2 = 1;
    i3 = 2;
  } else {
    i1 = iNear-1;
    i2 = iNear;
    i3 = iNear+1;
  }

  if(i1 == 0)
    ThrowError("Request for integral at radius " << xSky << " which is smaller than can be evaluated");

  double ret = quadInterp_.directEval(i1 * deltaInterp_, volumeIntegral_[type][i1],
				      i2 * deltaInterp_, volumeIntegral_[type][i2],
				      i3 * deltaInterp_, volumeIntegral_[type][i3], 
				      xSky);

  return ret;
}

/**.......................................................................
 * Return true if the shape parameters of the radial model are constant
 */
bool GenericRadiallySymmetric3DModel::shapeParametersAreFixed()
{
  ThrowColorError("Inheritor has not defined shapeParametersAreFixed()", "red");
  return false;
}

void GenericRadiallySymmetric3DModel::debugPrint() 
{
#ifdef TIMER_TEST
  COUTCOLOR("3D qt0 = " << qt0               << "s", "yellow");
  COUTCOLOR("3D qt1 = " << qt1               << "s", "yellow");
  COUTCOLOR("3D qt2 = " << qt2               << "s", "yellow");
  COUTCOLOR("3D qt3 = " << qt3               << "s", "yellow");
#endif

  Generic2DAngularModel::debugPrint();
}

/**.......................................................................
 * Return the 3D integral of the shape function, as a dimensionless number
 */
double GenericRadiallySymmetric3DModel::volumeIntegral(unsigned type, Angle& radius)
{
  //------------------------------------------------------------
  // Get the dimensionless ratio corresponding to this radius
  //------------------------------------------------------------

  double x = radius / thetaCore_;

  //------------------------------------------------------------
  // Integral is stored as a ratio to the solid angle pi*x^2 at the
  // dimensionless radius x.  Convert this to an actual solid angle
  // at the current radius
  //------------------------------------------------------------


  return interpolateVolumeIntegral(type, x);
}

/**.......................................................................
 * Return the 2D integral of the shape function, as a solid angle
 */
void GenericRadiallySymmetric3DModel::solidAngleIntegral(unsigned type, Angle& radius, SolidAngle& sa)
{
  //------------------------------------------------------------
  // Get the dimensionless ratio corresponding to this radius
  //------------------------------------------------------------

  double x = radius / thetaCore_;

  //------------------------------------------------------------
  // Integral is stored as a ratio to the solid angle pi*x^2 at the
  // dimensionless radius x.  Convert this to an actual solid angle
  // at the current radius
  //------------------------------------------------------------

  double radRad = radius.radians();
  sa.setSr(M_PI * radRad * radRad * interpolateArealIntegral(type, x));
}

/**.......................................................................
 * Return the 2D integral of the shape function, as a solid angle
 */
SolidAngle GenericRadiallySymmetric3DModel::solidAngleIntegral(unsigned type, Angle radius)
{
  //------------------------------------------------------------
  // Get the dimensionless ratio corresponding to this radius
  //------------------------------------------------------------

  double x = radius / thetaCore_;

  //------------------------------------------------------------
  // Integral is stored as a ratio to the solid angle pi*x^2 at the
  // dimensionless radius x.  Convert this to an actual solid angle
  // at the current radius
  //------------------------------------------------------------

  SolidAngle sa;
  double radRad = radius.radians();
  sa.setSr(M_PI * radRad * radRad * interpolateArealIntegral(type, x));

  return sa;
}

/**.......................................................................
 * Method to perform internal initialization prior to calculating
 * derived variates.  We split this into a cosmology-independent and
 * cosmology-dependent initialization for efficiency; i.e., if the
 * cosmology can vary, we only have to redo the cosmology dependent
 * part.
 */
void GenericRadiallySymmetric3DModel::initializeDerivedVariates(Model* caller)
{
}

VAR_DERIVE_FN(GenericRadiallySymmetric3DModel::deriveRadioLineIntegralNormalization)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)args;
  model->computeLineIntegralNormalization(DataSetType::DATASET_RADIO);
}

VAR_DERIVE_FN(GenericRadiallySymmetric3DModel::deriveXrayLineIntegralNormalization)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)args;
  model->computeLineIntegralNormalization(DataSetType::DATASET_XRAY_IMAGE);
}

void GenericRadiallySymmetric3DModel::computeLineIntegralNormalization(unsigned type)
{
  lineIntegralNormalization_[type].val_ = lineIntegral(type, xLim_);
}

//=======================================================================
// Methods to calculate derived variates
//=======================================================================

//-----------------------------------------------------------------------
// Volume integral factor needed for quantities integrated out to
// outerRadius_
//-----------------------------------------------------------------------

VAR_DERIVE_FN(GenericRadiallySymmetric3DModel::deriveVolumeIntegralFactor)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)args;
  model->fillVolumeIntegralFactor();
}

void GenericRadiallySymmetric3DModel::fillVolumeIntegralFactor()
{
  volumeIntegralFactor_.val_ = computeVolumeIntegralFactor(thetaInner_, thetaOuter_);
}

//-----------------------------------------------------------------------
// Ysph integrated out to outerRadius_
//-----------------------------------------------------------------------

VAR_DERIVE_FN(GenericRadiallySymmetric3DModel::deriveYsph)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)args;
  model->fillYsph();
}

void GenericRadiallySymmetric3DModel::fillYsph()
{
  integratedYsph(ySph_, volumeIntegralFactor_.val_);
}

//-----------------------------------------------------------------------
// Mtot integrated out to outerRadius_
//-----------------------------------------------------------------------

VAR_DERIVE_FN(GenericRadiallySymmetric3DModel::deriveMtot)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)args;
  model->fillMtot();
}

void GenericRadiallySymmetric3DModel::fillMtot()
{
  integratedTotalMass(mTot_, volumeIntegralFactor_.val_);
}

//-----------------------------------------------------------------------
// Mgas integrated out to outerRadius_
//-----------------------------------------------------------------------

VAR_DERIVE_FN(GenericRadiallySymmetric3DModel::deriveMgas)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)args;
  model->fillMgas();
}

void GenericRadiallySymmetric3DModel::fillMgas()
{
  integratedGasMass(mGas_, volumeIntegralFactor_.val_);
}

//-----------------------------------------------------------------------
// Pressure implied by the current normalization
//-----------------------------------------------------------------------

VAR_DERIVE_FN(GenericRadiallySymmetric3DModel::derivePe0)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)args;
  model->fillPe0();
}

/**.......................................................................
 * Calculate the pressure implied by the current normalization
 */
void GenericRadiallySymmetric3DModel::fillPe0()
{
  pressure(pe0_);
}

//-----------------------------------------------------------------------
// Mass scale factor
//-----------------------------------------------------------------------

VAR_DERIVE_FN(GenericRadiallySymmetric3DModel::deriveScaleFactorMass)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)args;
  model->fillScaleFactorMass();
}

/**.......................................................................
 * Calculate the pressure implied by the current normalization
 */
void GenericRadiallySymmetric3DModel::fillScaleFactorMass()
{
  scaleFactorMass_ = Constants::protonMass_ * (mue_.val_ * (scaleFactorArea_ / Constants::sigmaT_) * restMassThermalEnergyRatio_.val_);
}

//-----------------------------------------------------------------------
// Area scale factor
//-----------------------------------------------------------------------

VAR_DERIVE_FN(GenericRadiallySymmetric3DModel::deriveScaleFactorArea)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)args;
  model->fillScaleFactorArea();
}

/**.......................................................................
 * Calculate the pressure implied by the current normalization
 */
void GenericRadiallySymmetric3DModel::fillScaleFactorArea()
{
  scaleFactorArea_ = getCosmology()->dA_ * getCosmology()->dA_;
}

//-----------------------------------------------------------------------
// Scale factor to convert from current normalization to Compton Y
//-----------------------------------------------------------------------

VAR_DERIVE_FN(GenericRadiallySymmetric3DModel::deriveScaleFactorY)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)args;
  model->fillScaleFactorY();
}

/**.......................................................................
 * Calculate the pressure implied by the current normalization
 */
void GenericRadiallySymmetric3DModel::fillScaleFactorY()
{
  scaleFactorY_.val_ = getComptonYScaleFactor(normalizationFrequency_);
}

//-----------------------------------------------------------------------
// Scale factor for pressure calculations
//-----------------------------------------------------------------------

VAR_DERIVE_FN(GenericRadiallySymmetric3DModel::deriveScaleFactorPressure)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)args;
  model->fillScaleFactorPressure();
}

/**.......................................................................
 * Return a scale factor to convert from pressure to Compton-Y.  Call
 * this PN
 *
 * Given a pressure P, we convert to y by:
 *
 * y = (P / PN) * thetaCore
 *
 */
void GenericRadiallySymmetric3DModel::fillScaleFactorPressure()
{
  Energy restMassEnergy;
  restMassEnergy = Constants::electronMass_;

  Volume vol = Constants::sigmaT_ * getCosmology()->dA_;

  scaleFactorPressure_ = restMassEnergy / vol;
}

//-----------------------------------------------------------------------
// Inner radius
//-----------------------------------------------------------------------

VAR_DERIVE_FN(GenericRadiallySymmetric3DModel::deriveInnerRadius)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)args;
  model->fillInnerRadius();
}

/**.......................................................................
 * Get the inner radius
 */
void GenericRadiallySymmetric3DModel::fillInnerRadius()
{
  innerRadius_.setVal(getDoubleVal("innerRadius"), getParameter("innerRadius", true)->units_);
}

VAR_DERIVE_FN(GenericRadiallySymmetric3DModel::deriveThetaInner)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)args;
  model->fillThetaInner();
}

/**.......................................................................
 * Get the inner radius
 */
void GenericRadiallySymmetric3DModel::fillThetaInner()
{
  //  COUT("Inside fillThetaInner: innerRadius_ = " << innerRadius_.Mpc() << " " << getCosmology()->dA_);
  thetaInner_.setRadians(innerRadius_ / getCosmology()->dA_);
}

//-----------------------------------------------------------------------
// Outer radius
//-----------------------------------------------------------------------

VAR_DERIVE_FN(GenericRadiallySymmetric3DModel::deriveOuterRadius)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)args;
  model->fillOuterRadius();
}

/**.......................................................................
 * Get the outer radius
 */
void GenericRadiallySymmetric3DModel::fillOuterRadius()
{
  outerRadius_.setVal(getDoubleVal("outerRadius"), getParameter("outerRadius", true)->units_);
}

VAR_DERIVE_FN(GenericRadiallySymmetric3DModel::deriveThetaOuter)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)args;
  model->fillThetaOuter();
}

/**.......................................................................
 * Get the inner radius
 */
void GenericRadiallySymmetric3DModel::fillThetaOuter()
{
  thetaOuter_.setRadians(outerRadius_ / getCosmology()->dA_);
}

//-----------------------------------------------------------------------
// Rest mass / thermal energy ratio
//-----------------------------------------------------------------------

VAR_DERIVE_FN(GenericRadiallySymmetric3DModel::deriveRestMassThermalEnergyRatio)
{
  GenericRadiallySymmetric3DModel* model = (GenericRadiallySymmetric3DModel*)args;
  model->fillRestMassThermalEnergyRatio();
}

/**.......................................................................
 * Get the rest mass / thermal rnergy ratio
 */
void GenericRadiallySymmetric3DModel::fillRestMassThermalEnergyRatio()
{
  //------------------------------------------------------------
  // Calculate the scale factor to get from Compton-Y to mass
  //------------------------------------------------------------
  
  Temperature temp;
  temp.setVal(getDoubleVal("electronTemperature"), getParameter("electronTemperature", true)->units_);
  
  Energy thermalEnergy;
  thermalEnergy = temp;
  Energy restMassEnergy;
  restMassEnergy = Constants::electronMass_;
  
  restMassThermalEnergyRatio_.val_ = restMassEnergy / thermalEnergy;
}

double GenericRadiallySymmetric3DModel::computeVolumeIntegralFactor(Angle& thetaInner, Angle& thetaOuter)
{
  double xInner = thetaInner / thetaCore_;
  double xOuter = thetaOuter / thetaCore_;

  //  COUT("thetaInner = " << thetaInner << " thetaOuter = " << thetaOuter << " thetaCore = " << thetaCore_);
  //  COUT("xInner     = " << xInner << " xOuter = " << xOuter);
  //  COUT("xInner     = " << thetaInner.radians() << " " << thetaCore_.radians());


#if 0
  double vOuter = volumeIntegral(DataSetType::DATASET_RADIO, xOuter, 0);
  double vInner = volumeIntegral(DataSetType::DATASET_RADIO, xInner, 0);

  COUT("vInner = " << vInner << " vOuter = " << vOuter << " delta = " << vOuter - vInner);
  return vOuter - vInner;
#else
  return volumeIntegral(DataSetType::DATASET_RADIO, xOuter, 0) - volumeIntegral(DataSetType::DATASET_RADIO, xInner, 0);
#endif
}

/**.......................................................................
 * Calculate spherically integrated Y
 */
void GenericRadiallySymmetric3DModel::integratedYsph(Area& yInt, double& volumeIntegralFactor)
{
  double thetaCoreRad  = thetaCore_.radians();
  double thetaCoreRad2 = thetaCoreRad * thetaCoreRad;

  //  COUT("Inside integrated Ysph: " << " sfy = " << scaleFactorY_.val_ << " rn - " << radioNormalization_.val_
  //       << " scalefactorarea = " << scaleFactorArea_.squaredMpc() << " vif = " << volumeIntegralFactor
  //       << " lin = " << lineIntegralNormalization_[DataSetType::DATASET_RADIO].val_ << " thetaCore = " << thetaCore_.radians());


  double scaleFactor = (scaleFactorY_.val_ * radioNormalization_.val_ * thetaCoreRad2 / lineIntegralNormalization_[DataSetType::DATASET_RADIO].val_) * volumeIntegralFactor;
  yInt.setSquaredMpc(scaleFactorArea_.squaredMpc() * scaleFactor);
}

/**.......................................................................
 * Calculate spherically integrated mass out to the specified theta
 */
void GenericRadiallySymmetric3DModel::integratedGasMass(Mass& mass, Angle& theta)
{
  //  COUT("thetainner = " << thetaInner_.arcmin() << " theta = " << theta.arcmin());
  double volumeIntegralFactor = computeVolumeIntegralFactor(thetaInner_, theta);
  integratedGasMass(mass, volumeIntegralFactor);
}

/**.......................................................................
 * Calculate spherically integrated mass, with a precomputed volume
 * integral factor
 */
void GenericRadiallySymmetric3DModel::integratedGasMass(Mass& mass, double& volumeIntegralFactor)
{
  //  COUT("Inside integrated mass: " << " sfy = " << scaleFactorY_.val_ << " rn - " << radioNormalization_.val_
  //       << " scalefactormass = " << scaleFactorMass_.solarMass() << " vif = " << volumeIntegralFactor
  //       << " lin = " << lineIntegralNormalization_[DataSetType::DATASET_RADIO].val_);

  double thetaCoreRad  = thetaCore_.radians();
  double thetaCoreRad2 = thetaCoreRad * thetaCoreRad;

  double scaleFactor = (scaleFactorY_.val_ * radioNormalization_.val_ * thetaCoreRad2 / lineIntegralNormalization_[DataSetType::DATASET_RADIO].val_) * volumeIntegralFactor;
  mass.setSolarMass(scaleFactorMass_.solarMass() * scaleFactor);
}

/**.......................................................................
 * Calculate spherically integrated mass out to the specified theta
 */
void GenericRadiallySymmetric3DModel::integratedTotalMass(Mass& mass, Angle& theta)
{
  integratedGasMass(mass, theta);
  mass /= fGas_.val_;
}

/**.......................................................................
 * Calculate spherically integrated mass, with a precomputed volume
 * integral factor
 */
void GenericRadiallySymmetric3DModel::integratedTotalMass(Mass& mass, double& volumeIntegralFactor)
{
  integratedGasMass(mass, volumeIntegralFactor);
  mass /= fGas_.val_;
}

/**.......................................................................
 * Calculate the pressure implied by the current normalization
 */
void GenericRadiallySymmetric3DModel::pressure(Pressure& pressure)
{
  //  COUT("Inside pressure: " << " sfy = " << scaleFactorY_.val_ << " rn - " << radioNormalization_.val_
  //       << " scalefactorpress = " << scaleFactorPressure_.keVPerCm3() 
  //       << " lin = " << lineIntegralNormalization_[DataSetType::DATASET_RADIO].val_);

  double thetaCoreRad = thetaCore_.radians();
  double scaleFactor = (scaleFactorY_.val_ * radioNormalization_.val_) / (thetaCoreRad * lineIntegralNormalization_[DataSetType::DATASET_RADIO].val_);
  pressure.setKeVPerCm3(scaleFactorPressure_.keVPerCm3() * scaleFactor);
}

void GenericRadiallySymmetric3DModel::checkSetup()
{
  Generic2DAngularModel::checkSetup();

  if(!loadedFromOutputFile_)
    checkVar("thetaCore");
}

unsigned GenericRadiallySymmetric3DModel::getMaxNSigma()
{
  unsigned maxNSigma = 100;

  for(unsigned i=1; i < maxNSigma-1; i++) {
    double cdf = Sampler::gaussCdf((double)(i), 0.0, 1.0);

    double nOccur = (1.0 - cdf) * nTotal_;

    if(nOccur < 1)
      return i+1;
  }

  return maxNSigma;
}

GSL_HANDLER_FN(GenericRadiallySymmetric3DModel::intErrHandler)
{
  ThrowSimpleColorError("Function is not integrable", "red");
}
