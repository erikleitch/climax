#include "gcp/models/ClusterModel.h"
#include "gcp/models/CosmologyModel.h"

#include "gcp/util/Constants.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ClusterModel::ClusterModel() 
{
  initialize();
}

/**.......................................................................
 * Destructor.
 */
ClusterModel::~ClusterModel() {}

/**.......................................................................
 * Initialize this model
 */
void ClusterModel::initialize() 
{
  addComponent(thetaCore_);
  addComponentName(thetaCore_,  "thetaCore",  "An angular core radius");

  //------------------------------------------------------------
  // For cluster models, we allow the normalization to be specified in
  // pressure as well
  //------------------------------------------------------------

  radioNormalization_.addConversion("keV/cm^3",        1.0);
  radioNormalization_.addConversion("erg/cm^3",        1.0);

  pressureToComptonYScaleFactorNeedsInitializing_ = true;

  //------------------------------------------------------------
  // Initialize all components as fixed.  
  //------------------------------------------------------------

  initializeComponentsToFixed();
}

/**.......................................................................
 * Initialize the conversion from Pressure to ComptonY
 */
double ClusterModel::initializePressureToComptonYScaleFactor()
{
  Energy restMassEnergy;
  restMassEnergy = Constants::electronMass_;

  //------------------------------------------------------------
  // Calculate the factor by which to multiply the pressure to get
  // Compton-y
  //------------------------------------------------------------

  Length dA;
  dA.setGpc(1.0);

  Volume vol = dA * Constants::sigmaT_;

  Pressure press;
  press.setKeVPerCm3(restMassEnergy.keV() / vol.cubicCentimeters());

  std::string units = radioNormalization_.units();

  if(units == "keV/cm^3") {
    pressureToComptonYScaleFactor_ = 1.0/press.keVPerCm3();
  } else if(units == "erg/cm^3") {
    pressureToComptonYScaleFactor_ = 1.0/press.ergPerCm3();
  }

  COUT("Pressure scale factor is now: " <<     pressureToComptonYScaleFactor_);

  pressureToComptonYScaleFactorNeedsInitializing_ = false;
}

/**.......................................................................
 * Return the factor to convert from pressure in whatever units it
 * was specified, to Compton-Y
 */
double ClusterModel::pressureToComptonY()
{
  //------------------------------------------------------------
  // Initialize first, if we haven't done so already
  //------------------------------------------------------------

  if(pressureToComptonYScaleFactorNeedsInitializing_)
    initializePressureToComptonYScaleFactor();

  //------------------------------------------------------------
  // Final units of the image will be Compton-Y
  //------------------------------------------------------------

  units_ = Unit::UNITS_Y;

  //------------------------------------------------------------
  // Now return the current value of the scaling
  //------------------------------------------------------------

  return pressureToComptonYScaleFactor_ * getCosmology()->dA_.Gpc() * thetaCore_.radians();
}
