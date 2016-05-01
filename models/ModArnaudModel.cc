#include "gcp/models/ModArnaudModel.h"
#include "gcp/models/fastonebigheader.h"

#include "gcp/util/Variate.h"

using namespace std;

using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ModArnaudModel::ModArnaudModel() 
{
  //------------------------------------------------------------
  // Default these to the values Arnaud et al find are a good fit to
  // REXCESS clusters
  //------------------------------------------------------------

  p0_.setVal(1.0, "");
  c_.setVal(1.177, "");
  alpha_.setVal(1.0510, "");
  beta_.setVal(5.4905,  "");
  gamma_.setVal(0.3081, "");

  rfrac_.allowUnitless(true);
  rfrac_.setVal(1.0, "");

  addComponentName(rfrac_, "rfrac", "Value at which the model diverges to the Planck model (units of R500)");

  // Make the models piecewise continuous at rfrac:

  radioRat_ = 1.0;
  xrayRat_   = 1.0;

  radioRat_ = radialRadioModel(rfrac_.val_, 0) / radialRadioModel(rfrac_.val_ + 1e-12, 0);
  xrayRat_  = radialXrayModel(rfrac_.val_, 0) / radialXrayModel(rfrac_.val_ + 1e-12, 0);
}

/**.......................................................................
 * Destructor.
 */
ModArnaudModel::~ModArnaudModel() {}

/**.......................................................................
 * Method to calculate m500 when this variate is derived from
 * radioNormalization_
 */
void ModArnaudModel::fillM500()
{
  //------------------------------------------------------------
  // Calculate the pressure from the radio normalization
  //------------------------------------------------------------

  double val = scaleFactorY_.val_ * radioNormalization_.val_;
  double prat = val / (thetaCore_.radians() * lineIntegralNormalization_[DataSetType::DATASET_RADIO].val_);

  pressure_.setKeVPerCm3(scaleFactorPressure_.keVPerCm3() * prat);

  //------------------------------------------------------------
  // Pressure is now GNFW P500 * P0.  Convert to P500
  //------------------------------------------------------------

  pressure_ /= getGnfwNormalization();

  //------------------------------------------------------------
  //  This should now be directly comparable to Equation (3) of Nagai07
  //------------------------------------------------------------

  double h = getCosmology()->cosmology_.dimensionlessHubbleConstant();
  double h70 = getCosmology()->H0_.h(70.0);

  double rat = (pressure_.keVPerCm3() / (rescale_.value() * 1.65e-3 * pow(h, 8.0/3) * h70 * h70));
  double lnrat = log(rat) / (2.0/3 + 0.12);

  m500_.setSolarMass(exp(lnrat) * 3e14 / h70);
}

/**.......................................................................
 * Method to calculate radioNormalization_ when this variate is
 * derived from m500_
 */
void ModArnaudModel::fillRadioNormalization()
{
  //------------------------------------------------------------
  // Set the radio normalization according to the pressure
  //------------------------------------------------------------
  
  double hz    = getCosmology()->H_ / getCosmology()->H0_;
  double h70   = getCosmology()->H0_.h(70.0);
  double h7032 = sqrt(h70) * h70;
  double hz83  = pow(hz, 8.0/3);

  //------------------------------------------------------------
  // Calculate P = P500 * P0 according to the Arnaud normalization
  //------------------------------------------------------------

  pressure_.setKeVPerCm3(rescale_.value() * 1.65e-3 * hz83 * massPrefactorFirstOrder(h70) * h70 * h70 * getGnfwNormalization());
  
  //------------------------------------------------------------
  // Now get the radioNormalization in units of Compton Y
  //------------------------------------------------------------

  double prat = (pressure_ / scaleFactorPressure_);
  double val = prat * (thetaCore_.radians() * lineIntegralNormalization_[DataSetType::DATASET_RADIO].val_);

  //------------------------------------------------------------
  // Fill the normalization in whatever units it was specified in (our
  // value is calculated in ComptonY, so divide by the scale factor to
  // convert)
  //------------------------------------------------------------

  radioNormalization_.val_ = val / scaleFactorY_.val_;
}

/**.......................................................................
 * Return the mass prefactor defined in Arnaud et al, Equation 13
 *
 * Note that fabs() allows m500_ to go negative
 */
double ModArnaudModel::massPrefactorFirstOrder(double h70)
{
  double ex  = 2.0/3 + 0.12;
  double rat = (fabs(m500_.solarMass()) / 3e14) * h70;

  return pow(rat, ex);
}

double ModArnaudModel::getGnfwNormalization()
{
  double h70 = getCosmology()->H0_.h(70.0);
  return 8.403 / pow(h70, 3.0/2);
}

PgModel ModArnaudModel::pgModel()
{
  PgModel mod;

  mod.xMid_  = xOffset_.degrees();
  mod.yMid_  = yOffset_.degrees();
  mod.angle_.setDegrees(90);
  mod.rot_.setDegrees(90);
  mod.type_  = PgModel::TYPE_ARNAUD;

  double rad = thetaCore_.degrees();
  mod.xRad1_ = mod.xMid_ + rad; 
  mod.yRad1_ = mod.yMid_;

  mod.yRad2_ = mod.yMid_ + rad;
  mod.xRad2_ = mod.xMid_;

  mod.rad1_ = rad;
  mod.rad2_ = rad;

  return mod;
}

double ModArnaudModel::radialRadioModel(double x, void* params)
{
  double rat = 1.0;

  if(x > rfrac_.val_) {
    p0_.setVal(6.41, "");
    c_.setVal(1.81, "");
    alpha_.setVal(1.33, "");
    beta_.setVal(4.13,  "");
    gamma_.setVal(0.31, "");

    // Make the model piecewise continuous at rfrac:
    
    rat = radioRat_;
    
  } else {
    p0_.setVal(1.0, "");
    c_.setVal(1.177, "");
    alpha_.setVal(1.0510, "");
    beta_.setVal(5.4905,  "");
    gamma_.setVal(0.3081, "");
  }

  double bga = fac_.val_*(beta_.val_ - gamma_.val_)/alpha_.val_;
  double xg  = pow(x, gamma_.val_);
  double xa  = pow(x, alpha_.val_);

  return rat/(xg * pow(1.0 + xa, bga));
}

double ModArnaudModel::radialXrayModel(double x, void* params)
{
  double rat = 1.0;

  if(x > rfrac_.val_) {
    p0_.setVal(6.41, "");
    c_.setVal(1.81, "");
    alpha_.setVal(1.33, "");
    beta_.setVal(4.13,  "");
    gamma_.setVal(0.31, "");

    rat = xrayRat_;

  } else {
    p0_.setVal(1.0, "");
    c_.setVal(1.177, "");
    alpha_.setVal(1.0510, "");
    beta_.setVal(5.4905,  "");
    gamma_.setVal(0.3081, "");
  }

  double bga = 2*fac_.val_*(beta_.val_ - gamma_.val_)/alpha_.val_;
  double xg  = pow(x, 2*gamma_.val_);
  double xa  = pow(x, alpha_.val_);

  return rat/(xg * pow(1.0 + xa, bga));
}
