#include "gcp/util/Constants.h"

#include <cmath>

using namespace std;

using namespace gcp::util;

//------------------------------------------------------------
// Speed of light, cm/s
//------------------------------------------------------------

Speed Constants::lightSpeed_ = 
  Speed(Speed::CentimetersPerSec(), 2.99792458e10);

//------------------------------------------------------------
// Electron mass (g)
//------------------------------------------------------------

Mass Constants::electronMass_ = 
  Mass(Mass::Gram(), 9.1095e-28);

//------------------------------------------------------------
// Proton mass (g)
//------------------------------------------------------------

Mass Constants::protonMass_ = 
  Mass(Mass::Gram(), 1.672621778e-24);

//------------------------------------------------------------
// Solar mass (g)
//------------------------------------------------------------

Mass Constants::solarMass_ = 
  Mass(Mass::Gram(), 1.98855e33);

//------------------------------------------------------------
// Planck's constant, in cgs units (cm^2 g / s = erg s)
//------------------------------------------------------------

const double Constants::hPlanckCgs_          = 6.626176e-27;

//------------------------------------------------------------
// Boltzman constant:
//
// CGS units are (cm^2 g / s^2 / K) = erg / K)
// SI units are  (m^2 kg / s^2 / K)
//------------------------------------------------------------

const double Constants::kBoltzCgs_           = 1.3807e-16;
const double Constants::kBoltzSi_            = 1.3807e-23;

//------------------------------------------------------------
// Gravitational constant:
//
// CGS units are  cm^3 /  (g s^2)
// SI  units are   m^3 / (kg s^2)
//------------------------------------------------------------

const double Constants::gravitationalConstantCgs_ = 6.67384e-8;
const double Constants::gravitationalConstantSi_  = 6.67384e-11;

const double Constants::JyPerCgs_            = 1e23;

const double Constants::electronMassCgs_     = 9.1095e-28;
const double Constants::protonMassCgs_       = 1.6726e-24;
const double Constants::utSecPerSiderealSec_ = 1.002737909350795;

Area        Constants::sigmaT_               = Area(Area::SquaredCentimeters(), 6.6524e-25);
Temperature Constants::Tcmb_                 = Temperature(Temperature::Kelvin(), 2.726);
Length      Constants::au_                   = Length(Length::Meters(), 1.49597870e11);
Length      Constants::defaultEarthRadius_   = Length(Length::Meters(), 6378137.000);

const double Constants::milesPerDegree_ = 2*M_PI*Constants::defaultEarthRadius_.kilometers() / (360) * Length::milesPerKm_;

/**.......................................................................
 * Constructor.
 */
Constants::Constants() {}

/**.......................................................................
 * Destructor.
 */
Constants::~Constants() {}
