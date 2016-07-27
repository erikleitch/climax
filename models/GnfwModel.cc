#include "gcp/models/GnfwModel.h"

#include "gcp/fftutil/DataSetType.h"

#include "gcp/util/Constants.h"
#include "gcp/util/Exception.h"

using namespace std;

using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
GnfwModel::GnfwModel() 
{
  initialize();
}

/**.......................................................................
 * Destructor.
 */
GnfwModel::~GnfwModel() {}

/**.......................................................................
 * Initialize components of this model
 */
void GnfwModel::initialize()
{
  // The three parameters of the shape function are dimensionless
  // exponents

  alpha_.allowUnitless(true);
  beta_.allowUnitless(true);
  gamma_.allowUnitless(true);
  c_.allowUnitless(true);
  p0_.allowUnitless(true);
  fac_.allowUnitless(true);

  //------------------------------------------------------------
  // This is not really a parameter of the GNFW model, but I include
  // it here so that the same code can generate the beta-model
  //------------------------------------------------------------

  fac_.setVal(1.0, "");

  addComponent(alpha_);
  addComponent(beta_);
  addComponent(gamma_);
  addComponent(c_);
  addComponent(p0_);

  addComponentName(alpha_, "alpha", "The alpha exponent of the generalized NFW model");
  addComponentName(beta_,  "beta",  "The beta exponent of the generalized NFW model");
  addComponentName(gamma_, "gamma", "The gamma exponent of the generalized NFW model");
  addComponentName(c_,     "c",     "The concentration parameter of the generalized NFW model");
  addComponentName(p0_,    "p0",    "The (dimensionless) normalization of the generalized NFW model");

  //------------------------------------------------------------
  // Add a component to allow rescaling of the Arnaud P500 / M500
  // relation
  //------------------------------------------------------------

  addComponent(rescale_);
  addComponentName(rescale_, "rescale", "A rescale factor for the Arnaud pressure normalization");

  rescale_.allowUnitless(true);
  rescale_.setVal(1.0, "");

  //------------------------------------------------------------
  // Add malleable components
  //------------------------------------------------------------

  //------------------------------------------------------------
  // The M500 implied by the GNFW model pressure normalization.  Note
  // that in the GNFW model, both the M500 and the radio normalization
  // are malleable -- that is, either can be specified as the primary
  // variate, and the other can be derived from it.  If M500 is
  // specified, we also derive R500 from it and therefore thetaCore,
  // so thetaCore is also declared as malleable
  //------------------------------------------------------------

  addMalleableComponent(m500_, "Msolar");
  addComponentName(m500_,  "m500",  "Total mass at R500.  Can be a primary variate, or derived from the pressure normalization");
  addComponentName(m500_,  "M500",  "Total mass at R500.  Can be a primary variate, or derived from the pressure normalization");
  m500_.deriveWith(GnfwModel::deriveM500, this);

  m500_.dependsOn(lineIntegralNormalization_[DataSetType::DATASET_RADIO]);
  m500_.dependsOn(thetaCore_);
  m500_.dependsOn(scaleFactorY_);
  m500_.dependsOn(scaleFactorPressure_);

  //------------------------------------------------------------
  // If m500 is specified, then thetaCore will be derived from rhoCrit
  // via m500 and r500, so it can now depend on these variates too
  //------------------------------------------------------------

  initializeMalleableComponent(thetaCore_, "arcmin");
  thetaCore_.deriveWith(GnfwModel::deriveThetaCore, this);
  thetaCore_.dependsOn(r500_);

  //------------------------------------------------------------
  // Likewise the radio normalization has dependencies if m500 is
  // specified
  //------------------------------------------------------------

  initializeMalleableComponent(radioNormalization_, "comptony");
  radioNormalization_.deriveWith(GnfwModel::deriveRadioNormalization, this);
  radioNormalization_.dependsOn(lineIntegralNormalization_[DataSetType::DATASET_RADIO]);
  radioNormalization_.dependsOn(thetaCore_);
  radioNormalization_.dependsOn(m500_);
  radioNormalization_.dependsOn(scaleFactorY_);
  radioNormalization_.dependsOn(scaleFactorPressure_);

  //------------------------------------------------------------
  // Add derived variates specific to the GNFW model to our map of
  // variates
  //------------------------------------------------------------

  volumeIntegralFactorR500_.allowUnitless(true);
  addDerivedComponent(volumeIntegralFactorR500_,    "", false);
  addComponentName(volumeIntegralFactorR500_,  "vif500", "");
  volumeIntegralFactorR500_.deriveWith(GnfwModel::deriveVolumeIntegralFactorR500, this);
  volumeIntegralFactorR500_.dependsOn(r500_);
  volumeIntegralFactorR500_.dependsOn(thetaCore_);

  //------------------------------------------------------------
  // The R500 implied by the GNFW model fits.  If derived from M500,
  // it depends only on M500.  Else it depends on thetaCore
  //------------------------------------------------------------

  addDerivedComponent(r500_,    "Mpc");
  addComponentName(r500_,  "R500", "The self-similar R500 implied by the normalization of the Arnaud model (derived)");
  r500_.deriveWith(GnfwModel::deriveR500, this);
  r500_.dependsOn(m500_);
  r500_.dependsOn(thetaCore_);

  //------------------------------------------------------------
  // The T500 implied by the GNFW model fits, when derived from M500.
  //------------------------------------------------------------

  addDerivedComponent(t500_,    "keV");
  addComponentName(t500_,  "T500", "The self-similar T500 implied by the normalization of the Arnaud model (derived)");
  t500_.deriveWith(GnfwModel::deriveT500, this);
  t500_.dependsOn(m500_);
  t500_.dependsOn(r500_);
  t500_.dependsOn(mu_);

  //------------------------------------------------------------
  // Add Ysph500 and set default units to Mpc^2
  //------------------------------------------------------------

  addDerivedComponent(ySph500_, "Mpc^2");
  addComponentName(ySph500_,    "Ysph500", "Compton Y, spherically integrated out to the self-similar R500 implied by the GNFW model (derived)");
  ySph500_.deriveWith(GnfwModel::deriveYsph500, this);
  ySph500_.dependsOn(lineIntegralNormalization_[DataSetType::DATASET_RADIO]);
  ySph500_.dependsOn(volumeIntegralFactorR500_);
  ySph500_.dependsOn(radioNormalization_);
  ySph500_.dependsOn(thetaCore_);
  ySph500_.dependsOn(scaleFactorY_);
  ySph500_.dependsOn(scaleFactorArea_);

  //------------------------------------------------------------
  // Mgas integrated out to this R500 requires whatever integrated mas
  // requires plus R500
  //------------------------------------------------------------

  addDerivedComponent(mGas500_, "Msolar");
  addComponentName(mGas500_, "Mgas500", "Gas mass, spherically integrated out to the self-similar R500 implied by the GNFW model (derived)");
  mGas500_.deriveWith(GnfwModel::deriveMgas500, this);
  mGas500_.dependsOn(r500_);
  mGas500_.dependsOn(lineIntegralNormalization_[DataSetType::DATASET_RADIO]);
  mGas500_.dependsOn(volumeIntegralFactorR500_);
  mGas500_.dependsOn(radioNormalization_);
  mGas500_.dependsOn(thetaCore_);
  mGas500_.dependsOn(scaleFactorMass_);


  addDerivedComponent(scaleFactorMassT500_, "Msolar", false);
  addComponentName(scaleFactorMassT500_,    "scaleFactorMassT500");
  scaleFactorMassT500_.deriveWith(GnfwModel::deriveScaleFactorMassT500, this);
  scaleFactorMassT500_.dependsOn(scaleFactorArea_);
  scaleFactorMassT500_.dependsOn(mue_);
  scaleFactorMassT500_.dependsOn(t500_);

  //------------------------------------------------------------
  // Mtot integrated out to this R500 requires whatever integrated
  // mass requires, plus R500
  //------------------------------------------------------------

  addDerivedComponent(mTot500_, "Msolar");
  addComponentName(mTot500_, "Mtot500", "Total mass, spherically integrated out to the self-similar R500 implied by the GNFW model (derived)");
  mTot500_.deriveWith(GnfwModel::deriveMtot500, this);
  mTot500_.dependsOn(r500_);
  mTot500_.dependsOn(lineIntegralNormalization_[DataSetType::DATASET_RADIO]);
  mTot500_.dependsOn(volumeIntegralFactorR500_);
  mTot500_.dependsOn(radioNormalization_);
  mTot500_.dependsOn(thetaCore_);
  mTot500_.dependsOn(fGas_);
  mTot500_.dependsOn(mue_);
  mTot500_.dependsOn(scaleFactorMassT500_);

  //------------------------------------------------------------
  // And initialize all components to fixed
  //------------------------------------------------------------

  initializeComponentsToFixed();
}

/**.......................................................................
 * Check setup of this model
 */
void GnfwModel::checkSetup()
{
  //------------------------------------------------------------
  // If the run file specifically called out an m500, then we are
  // deriving the model normalization from the Arnaud mass-pressure
  // relationship.  In this case, check that the radioNormalization_
  // and thetaCore are not also specified, because they are now both
  // derived parameters
  //------------------------------------------------------------

  if(m500_.wasSpecified_ && !m500_.isDerived_) {

    if(radioNormalization_.wasSpecified_ && !radioNormalization_.isDerived_) {
      ThrowSimpleColorError(std::endl << "You can't specify both " << name_ << ".M500 and " << name_ 
                            << ".Sradio, but must specify one or the other", "red");
    }

    if(thetaCore_.wasSpecified_ && !thetaCore_.isDerived_) {
      ThrowSimpleColorError(std::endl << "If " << name_ << ".M500 is specified, you cannot also specify " 
                            << name_ << ".thetaCore, since it will be derived from " << name_  << ".M500",
		      "red");
    }

    //------------------------------------------------------------
    // Remove the dependency of R500 on thetaCore_: if M500 is specified,
    // R500 is determined by M500 alone
    //------------------------------------------------------------

    r500_.doesntDependOn(thetaCore_);

    //------------------------------------------------------------
    // Remove the dependency of M500 on the scale factor -- only
    // needed if m500 is derived from the radio normalization
    //------------------------------------------------------------

    m500_.doesntDependOn(scaleFactorY_);
    m500_.doesntDependOn(scaleFactorPressure_);

    //------------------------------------------------------------
    // And setup the radioNormalization_ variate so that subsequent
    // checks will succeed
    //
    // If no-one has requested to display this variate, don't display
    // it by default
    //------------------------------------------------------------

    if(radioNormalization_.isDefaultDisplay_)
      radioNormalization_.setDisplay(false);

    if(!loadedFromOutputFile_)
      radioNormalization_.wasSpecified_ = true;

    initializeDerivedComponent(radioNormalization_);
    radioNormalization_.isVariable() = m500_.isVariable();

    if(thetaCore_.isDefaultDisplay_)
      thetaCore_.setDisplay(false);

    if(!loadedFromOutputFile_)
      thetaCore_.wasSpecified_ = true;

    initializeDerivedComponent(thetaCore_);
    thetaCore_.isVariable() = m500_.isVariable();

    m500_.isDerived_ = false;

    //------------------------------------------------------------
    // Else if the radioNormalization was explicitly called out, then
    // we are fitting the normalization directly
    //------------------------------------------------------------

  } else if(radioNormalization_.wasSpecified_ && !radioNormalization_.isDerived_) {

    if(m500_.wasSpecified_ && !m500_.isDerived_) {
      ThrowSimpleColorError(std::endl << "You can't specify both " << name_ << ".M500 and " << name_ 
                            << ".Sradio, but must specify one or the other", "red");
    }

    radioNormalization_.isDerived_ = false;
    m500_.isDerived_ = true;
    m500_.isVariable() = radioNormalization_.isVariable();

    //------------------------------------------------------------
    // Remove the dependency of the radio normalization on the scale factor -- only
    // needed if the radio normalization is derived from m500
    //------------------------------------------------------------

    thetaCore_.doesntDependOn(r500_);

    radioNormalization_.doesntDependOn(m500_);
    radioNormalization_.doesntDependOn(scaleFactorY_);
    radioNormalization_.doesntDependOn(scaleFactorPressure_);

  } else if(xrayNormalization_.wasSpecified_) {

    //------------------------------------------------------------
    // Remove the dependency of the radio normalization on the scale factor -- only
    // needed if the radio normalization is derived from m500
    //------------------------------------------------------------

    thetaCore_.doesntDependOn(r500_);

    radioNormalization_.doesntDependOn(m500_);
    radioNormalization_.doesntDependOn(scaleFactorY_);
    radioNormalization_.doesntDependOn(scaleFactorPressure_);
  }

  //------------------------------------------------------------
  // Now check that exponents have been specified for this model
  //------------------------------------------------------------

  checkVar("alpha");
  checkVar("beta");
  checkVar("gamma");

  //------------------------------------------------------------
  // Finally, perform any additional base-class checks
  //------------------------------------------------------------

  GenericRadiallySymmetric3DModel::checkSetup();
}

/**.......................................................................
 * Return the radial Generalized NFW (GNFW) model.  Here x is the
 * dimensionless r/r_c parameter.
 *
 * The empirical shape function that Nagai et al come up with is:
 * 
 *         -g  /       a \ (g - b)/a
 * p(x) = x   |  1 + x    |
 *             \         /
 *
 * (So that with g = 0, a = 2 and b = 3*beta, for example, we would
 *  recover the beta model)
 */
double GnfwModel::radialModel(unsigned type, double x, void* params)
{
  switch (type) {
  case DataSetType::DATASET_RADIO:
    return radialRadioModel(x, params);
    break;
  case DataSetType::DATASET_XRAY_IMAGE:
    return radialXrayModel(x, params);
    break;
  case DataSetType::DATASET_GENERIC:
    return radialRadioModel(x, params);
    break;
  default:
    ThrowColorError("Unsupported dataset type: " << type, "red");
    return 0.0;
    break;
  }
}

double GnfwModel::radialRadioModel(double x, void* params)
{
  double bga = fac_.val_*(beta_.val_ - gamma_.val_)/alpha_.val_;
  double xg  = pow(x, gamma_.val_);
  double xa  = pow(x, alpha_.val_);

  return 1.0/(xg * pow(1.0 + xa, bga));
}

double GnfwModel::radialXrayModel(double x, void* params)
{
  double bga = 2*fac_.val_*(beta_.val_ - gamma_.val_)/alpha_.val_;
  double xg  = pow(x, 2*gamma_.val_);
  double xa  = pow(x, alpha_.val_);

  return 1.0/(xg * pow(1.0 + xa, bga));
}

/**.......................................................................
 * Return true if all parameters affecting the shape (note: not the
 * scale) of the model are fixed.  However, if thetaCore has no prior
 * specified, this means that we can't precompute its possible range,
 * so return false so that the interpolation grid is recomputed on
 * each trial of thetaCore
 */
bool GnfwModel::shapeParametersAreFixed()
{
  return !(alpha_.isVariable() || beta_.isVariable() || gamma_.isVariable()) &&
    thetaCore_.prior().getType() != Distribution::DIST_UNSPEC;
}

/**.......................................................................
 * The relation between my x and Gnfw x is:
 *
 * x_climax = c500 * x_gnfw,
 *
 * where x_gnfw = r/R500
 * 
 * so to get the x_gnfw that corresponds to x_climax, we need to
 * divide by c500
 */
double GnfwModel::gnfwX(double x)
{
  return x * c_.val_;
}

//=======================================================================
// Methods for computing derived variates
//=======================================================================

VAR_DERIVE_FN(GnfwModel::deriveM500)
{
  GnfwModel* gnfwModel = (GnfwModel*) args;
  gnfwModel->fillM500();
}

/**.......................................................................
 * Calculate the M500 implied by the current normalization
 */
void GnfwModel::fillM500()
{
  double thetaCoreRad  = thetaCore_.radians();
  double scaleFactor = (scaleFactorY_.val_ * radioNormalization_.val_) / (thetaCoreRad * lineIntegralNormalization_[DataSetType::DATASET_RADIO].val_);
  pressure_.setKeVPerCm3(scaleFactorPressure_.keVPerCm3() * scaleFactor);

  //------------------------------------------------------------
  // Pressure is now GNFW P500 * P0.  Convert to P500
  //------------------------------------------------------------

  pressure_ /= getGnfwNormalization();

  //------------------------------------------------------------
  //  This is now directly comparable to Equation (3) of Nagai07
  //------------------------------------------------------------

  double E = getCosmology()->cosmology_.E();
  double h100 = getCosmology()->H0_.h(100.0);

  double rat = (pressure_.ergPerCm3()) / (rescale_.val_ * 1.45e-11 * pow(E, 8.0/3) * h100 * h100);
  double rat1 = pow(rat, 3.0/2);

  m500_.setSolarMass(rat1 * 1e15 / h100);
}

VAR_DERIVE_FN(GnfwModel::deriveRadioNormalization)
{
  GnfwModel* gnfwModel = (GnfwModel*) args;
  gnfwModel->fillRadioNormalization();
}

/**.......................................................................
 * Method to calculate radioNormalization_ when this variate is
 * derived from m500_
 */
void GnfwModel::fillRadioNormalization()
{
  //------------------------------------------------------------
  // Set the radio normalization according to the pressure
  //------------------------------------------------------------
  
  double Ez83  = pow(getCosmology()->cosmology_.E(), 8.0/3);
  double h100  = getCosmology()->H0_.h(100.0);

  //------------------------------------------------------------
  // Calculate P = P500 * P0 according to the Nagai normalization
  // 
  // Note that fabs() allows m500 to be negative
  //------------------------------------------------------------
  
  double massPrefactor = pow(fabs(m500_.solarMass()) / (1e15 / h100), 2.0/3);
  pressure_.setErgPerCm3(rescale_.val_ * 1.45e-11 * massPrefactor * Ez83 * getGnfwNormalization() * h100 * h100);
  
  //------------------------------------------------------------
  // Now get the radioNormalization in units of Compton Y
  //------------------------------------------------------------
  
  double val = (pressure_ / scaleFactorPressure_) * thetaCore_.radians() * lineIntegralNormalization_[DataSetType::DATASET_RADIO].val_;

  //------------------------------------------------------------
  // Fill the normalization in whatever units it was specified in (our
  // value is calculated in ComptonY, so divide by the scale factor to
  // convert)
  //------------------------------------------------------------

  radioNormalization_.val_ = val / scaleFactorY_.val_;
}

VAR_DERIVE_FN(GnfwModel::deriveT500)
{
  GnfwModel* gnfwModel = (GnfwModel*) args;
  gnfwModel->fillT500();
}

/**.......................................................................
 * Calculate the T500 implied by the current M500
 */
void GnfwModel::fillT500()
{
  // Note that fabs() allows m500_ to go negative

  t500_.setK(mu_.val_ * Constants::protonMass_.g() * Constants::gravitationalConstantCgs_ * fabs(m500_.g()) / (2 * r500_.cm() * Constants::kBoltzCgs_));
}

Temperature GnfwModel::getTemperature(Mass mDelta, double delta)
{
  Length rDelta = getRDeltaFromRhoCrit(mDelta, delta);

  Temperature tDelta;
  tDelta.setK(mu_.val_ * Constants::protonMass_.g() * Constants::gravitationalConstantCgs_ * mDelta.g() / (2 * rDelta.cm() * Constants::kBoltzCgs_));
  return tDelta;
}

VAR_DERIVE_FN(GnfwModel::deriveR500)
{
  GnfwModel* gnfwModel = (GnfwModel*) args;
  gnfwModel->fillR500();
}

/**.......................................................................
 * Calculate the R500 implied by the current thetaCore
 */
void GnfwModel::fillR500()
{
  if(m500_.wasSpecified_ && !m500_.isDerived_) {
    fillR500FromRhoCrit();
  } else {
    fillR500FromThetaCore();
  }
}

void GnfwModel::fillR500FromThetaCore()
{
  double thetaRad = thetaCore_.radians();
  double c500     = c_.val_;
  double r500Mpc  = getCosmology()->dA_.Mpc() * thetaRad * c500;

  r500_.setMpc(getCosmology()->dA_.Mpc() * thetaRad * c500);
}

void GnfwModel::fillR500FromRhoCrit()
{
  Density rho = getCosmology()->cosmology_.criticalDensity();
  Volume vol = (m500_ / (500 * 4*M_PI/3)) / rho;

  // Note that fabs() of the volume allows for m500_ to be negative

  r500_.setCentimeters(pow(fabs(vol.cubicCentimeters()), 1.0/3));
}

Length GnfwModel::getRDeltaFromRhoCrit(Mass mDelta, double delta)
{
  Density rho = getCosmology()->cosmology_.criticalDensity();
  Volume vol = (mDelta / (delta * 4*M_PI/3)) / rho;

  Length rDelta;
  rDelta.setCentimeters(pow(vol.cubicCentimeters(), 1.0/3));
  return rDelta;
}

/**.......................................................................
 * If M500 was specified, this implies an R500, which specifies
 * thetaCore, via thetaCore = rCore / dA, and rCore = R500/c500.  This
 * functions does a binary search on the mass integral to find the
 * R500 that corresponds to the specified mass
 */
void GnfwModel::interpolateForThetaCore()
{
  Angle thetaLow, thetaHigh, thetaMid;
  Mass massTest, massLow, massHigh, massMid;
  double eps = 0.01;

  thetaLow  = thetaInner_;
  thetaHigh = 10*thetaInner_;
  thetaMid  = (thetaHigh + thetaLow)/2.0;
  Angle thetaInt;

  bool stop = false;

  COUT("Theta inner = " << thetaInner_.arcmin());
  COUT("Theta core  = " << thetaCore_.arcmin());

  // Low mass is by definition 0

  massLow.setSolarMass(0.0);

  thetaCore_.setRadians(thetaHigh.radians());
  thetaInt.setRadians(thetaHigh.radians() * c_.val_);
  integratedTotalMassT500(massHigh, thetaInt);

  do {

    thetaCore_.setRadians(thetaMid.radians());
    thetaInt.setRadians(thetaMid.radians() * c_.val_);
    integratedTotalMassT500(massMid,  thetaInt);

    COUT("thetaLow  = " << thetaLow.arcmin()  << " Mass Low  = " << massLow.solarMass());
    COUT("thetaMid  = " << thetaMid.arcmin()  << " Mass Mid  = " << massMid.solarMass());
    COUT("thetaHigh = " << thetaHigh.arcmin() << " Mass High = " << massHigh.solarMass());
    COUT("Mass M500 = " << m500_.solarMass());

    if(m500_ > massLow && m500_ <= massMid) {


      if(fabs((massLow.solarMass() - m500_.solarMass())/m500_.solarMass()) < eps) {
	thetaCore_ = thetaLow;
	massTest = massLow;
	stop = true;
      }

      if(fabs((massMid.solarMass() - m500_.solarMass())/m500_.solarMass()) < eps) {
	thetaCore_ = thetaMid;
	massTest = massMid;
	stop = true;
      }

      massHigh  = massMid;
      thetaHigh = thetaMid;
      thetaMid  = (thetaLow + thetaHigh)/2.0;

    } else if(m500_ > massMid && m500_ <= massHigh) {


      if(fabs((massHigh.solarMass() - m500_.solarMass())/m500_.solarMass()) < eps) {
	thetaCore_ = thetaHigh;
	massTest = massHigh;
	stop = true;
      }

      massLow   = massMid;
      thetaLow  = thetaMid;
      thetaMid  = (thetaLow + thetaHigh)/2.0;

    } else if(m500_ > massHigh) {

      thetaLow  = thetaHigh;
      massLow   = massHigh;
      thetaHigh = 2*thetaHigh;

      thetaCore_.setRadians(thetaHigh.radians());
      thetaInt.setRadians(thetaHigh.radians() * c_.val_);
      integratedTotalMassT500(massHigh,  thetaInt);

      thetaMid  = (thetaLow + thetaHigh)/2.0;

    }

  } while (!stop);

#if 1
  COUT("(1) Theta core is now: " << thetaCore_.arcmin() << " arcmin ");

  fillR500();

  COUT("(1) R500 is now: " << r500_.Mpc() << " Mpc");
#endif
}

/**.......................................................................
 * If M500 was specified, this implies an R500, which specifies
 * thetaCore, via thetaCore = rCore / dA, and rCore = R500/c500.  This
 * functions does a binary search on the mass integral to find the
 * R500 that corresponds to the specified mass
 */
void GnfwModel::interpolateForRescale()
{
  double scaleLow, scaleHigh, scaleMid, rescale;
  Mass massTest, massLow, massHigh, massMid;
  double eps = 0.01;

  scaleLow  = 0.1;
  scaleHigh = 10;
  scaleMid  = (scaleHigh + scaleLow)/2.0;

  bool stop = false;

  rescale_.val_ = scaleLow;
  fillRadioNormalization();
  integratedTotalMassT500(massLow,  volumeIntegralFactorR500_.val_);

  rescale_.val_ = scaleHigh;
  fillRadioNormalization();
  integratedTotalMassT500(massHigh, volumeIntegralFactorR500_.val_);

  do {

    rescale_.val_ = scaleMid;
    fillRadioNormalization();
    integratedTotalMassT500(massMid, volumeIntegralFactorR500_.val_);

    if(m500_ > massLow && m500_ <= massMid) {

      if(fabs((massLow.solarMass() - m500_.solarMass())/m500_.solarMass()) < eps) {
	rescale = scaleLow;
	massTest = massLow;
	stop = true;
      }

      if(fabs((massMid.solarMass() - m500_.solarMass())/m500_.solarMass()) < eps) {
	rescale = scaleMid;
	massTest = massMid;
	stop = true;
      }

      massHigh  = massMid;
      scaleHigh = scaleMid;
      scaleMid  = (scaleLow + scaleHigh)/2.0;

    } else if(m500_ > massMid && m500_ <= massHigh) {

      if(fabs((massHigh.solarMass() - m500_.solarMass())/m500_.solarMass()) < eps) {
	rescale  = scaleHigh;
	massTest = massHigh;
	stop = true;
      }

      massLow  = massMid;
      scaleLow = scaleMid;
      scaleMid = (scaleLow + scaleHigh)/2.0;

    } else if(m500_ > massHigh) {

      scaleLow  = scaleHigh;
      massLow   = massHigh;
      scaleHigh = 2*scaleHigh;

      rescale_.val_ = scaleHigh;
      fillRadioNormalization();
      integratedTotalMassT500(massHigh, volumeIntegralFactorR500_.val_);

      scaleMid  = (scaleLow + scaleHigh)/2.0;
    }

  } while (!stop);

  rescale_.val_ = rescale;
}

double GnfwModel::getGnfwNormalization()
{
  return p0_.val_;
}

//=======================================================================
// Methods to compute derived variates
//=======================================================================

//-----------------------------------------------------------------------
// Volume integral factor needed for R500 calculations
//-----------------------------------------------------------------------

VAR_DERIVE_FN(GnfwModel::deriveVolumeIntegralFactorR500)
{
  GnfwModel* gnfwModel = (GnfwModel*) args;
  gnfwModel->fillVolumeIntegralFactorR500();
}

/**.......................................................................
 * Calculate the MGas500 implied by the current normalization
 */
void GnfwModel::fillVolumeIntegralFactorR500() 
{
  thetaR500_.setRadians(r500_ / getCosmology()->dA_ * c_.val_);
  volumeIntegralFactorR500_.val_ = computeVolumeIntegralFactor(thetaInner_, thetaR500_);
}

//-----------------------------------------------------------------------
// Ysph integrated out to R500
//-----------------------------------------------------------------------

VAR_DERIVE_FN(GnfwModel::deriveYsph500)
{
  GnfwModel* model = (GnfwModel*)args;
  model->fillYsph500();
}

void GnfwModel::fillYsph500()
{
  integratedYsph(ySph500_, volumeIntegralFactorR500_.val_);
}

//-----------------------------------------------------------------------
// Mgas integrated to R500
//-----------------------------------------------------------------------

VAR_DERIVE_FN(GnfwModel::deriveMgas500)
{
  GnfwModel* gnfwModel = (GnfwModel*) args;
  gnfwModel->fillMgas500();
}

void GnfwModel::fillMgas500() 
{
  integratedGasMass(mGas500_, volumeIntegralFactorR500_.val_);
}

//-----------------------------------------------------------------------
// Mtot integrated to R500
//-----------------------------------------------------------------------

VAR_DERIVE_FN(GnfwModel::deriveMtot500)
{
  GnfwModel* gnfwModel = (GnfwModel*) args;
  gnfwModel->fillMtot500();
}

/**.......................................................................
 * Calculate the MTot500 implied by the current normalization
 */
void GnfwModel::fillMtot500() 
{
  integratedTotalMassT500(mTot500_, volumeIntegralFactorR500_.val_);
}

/**.......................................................................
 * Calculate spherically integrated mass out to the specified theta,
 * assuming characteristic temperature of an isothermal sphere derived
 * from the specified m500_
 */
void GnfwModel::integratedGasMassT500(Mass& mass, Angle& theta)
{
  integratedGasMassT500(mass, volumeIntegralFactorR500_.val_);
}

void GnfwModel::integratedTotalMassT500(Mass& mass, Angle& theta)
{
  integratedGasMassT500(mass, theta);
  mass /= fGas_.val_;
}

/**.......................................................................
 * Calculate spherically integrated mass, with a precomputed volume
 * integral factor
 */
void GnfwModel::integratedGasMassT500(Mass& mass, double& volumeIntegralFactor)
{
  //  COUT("Inside integrated mass: " << " sfy = " << scaleFactorY_.val_ << " rn - " << radioNormalization_.val_
  //       << " scalefactormass = " << scaleFactorMassT500_.solarMass() << " vif = " << volumeIntegralFactor
  //       << " lin = " << lineIntegralNormalization_[DataSetType::DATASET_RADIO].val_);

  double thetaCoreRad  = thetaCore_.radians();
  double thetaCoreRad2 = thetaCoreRad * thetaCoreRad;

  double scaleFactor = (scaleFactorY_.val_ * radioNormalization_.val_ * thetaCoreRad2 / lineIntegralNormalization_[DataSetType::DATASET_RADIO].val_) * volumeIntegralFactor;

  mass.setSolarMass(scaleFactorMassT500_.solarMass() * scaleFactor);
}

void GnfwModel::integratedTotalMassT500(Mass& mass, double& volumeIntegralFactor)
{
  integratedGasMassT500(mass, volumeIntegralFactor);
  mass /= fGas_.val_;
}

VAR_DERIVE_FN(GnfwModel::deriveScaleFactorMassT500)
{
  GnfwModel* model = (GnfwModel*)args;
  model->fillScaleFactorMassT500();
}

/**.......................................................................
 * Return the scale factor to convert to mass
 */
void GnfwModel::fillScaleFactorMassT500()
{
  Energy thermalEnergy;
  thermalEnergy = t500_;

  Energy restMassEnergy;
  restMassEnergy = Constants::electronMass_;

  scaleFactorMassT500_ = Constants::protonMass_ * (mue_.val_ * (scaleFactorArea_ / Constants::sigmaT_) * (restMassEnergy / thermalEnergy));
}

VAR_DERIVE_FN(GnfwModel::deriveThetaCore)
{
  GnfwModel* gnfwModel = (GnfwModel*) args;
  gnfwModel->fillThetaCore();
}

void GnfwModel::fillThetaCore()
{
  calculateThetaCoreFromRhoCrit();
}

void GnfwModel::calculateThetaCoreFromRhoCrit()
{
  thetaCore_.setRadians((r500_ / getCosmology()->dA_) / c_.val_);
}
