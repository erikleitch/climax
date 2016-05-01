#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/Constants.h"
#include "gcp/util/Cosmology.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Debug.h"
#include "gcp/util/Energy.h"
#include "gcp/util/Mass.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "m500",     "5e15",   "d", "m500 (in solar masses)"},
  { "Te",       "10",     "d", "Te (keV)"},
  { "R",        "1.34",    "d", "R (Mpc)"},
  { "y",        "3.7e-4", "d", "y"},
  { "z",        "0.168", "d", "z"},
  { "p",        "2e-10",  "d", "P (erg/cm^3)"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

Mass nagaiM500(Pressure& p, double z);
Mass arnaudM500(Pressure& p, double z);
Pressure nagaiP500(Mass& m, double z);
Pressure arnaudP500(Mass& m500, double z);
double arnaudPrefactor();

int Program::main()
{
  Pressure p;
  p.setErgPerCm3(Program::getDoubleParameter("p"));
  double z = Program::getDoubleParameter("z");

  Mass m500 = nagaiM500(p, z);
  COUT("Nagai mass = " << m500.solarMass());

  m500 = arnaudM500(p, z);
  COUT("Arnaud mass = " << m500.solarMass());

  m500.setSolarMass(Program::getDoubleParameter("m500"));

  p = nagaiP500(m500, z);
  COUT("Nagai P500 = " << p.ergPerCm3());
  p = arnaudP500(m500, z);
  COUT("Arnaud P500 = " << p.ergPerCm3());

  COUT("Arnaud prefactor should be: " << arnaudPrefactor());
  return 0;

  HubbleConstant H0;
  H0.setH100(0.7);

  Cosmology cosmo;
  cosmo.setRedshift(0.168);
  cosmo.setH0(H0);
  cosmo.setOmegaM(0.3);
  cosmo.setOmegaL(0.7);

  Energy restMassEnergy;
  restMassEnergy = Constants::electronMass_;

  Volume vol = Constants::sigmaT_ * cosmo.angularDiameterDistance();

  Angle theta;
  theta.setArcMinutes(3.8);

  Pressure norm = restMassEnergy / vol;

  m500.setSolarMass(Program::getDoubleParameter("m500"));

  double hz = cosmo.dimensionlessHubbleConstant();

  double val  = 1.65e-3 * pow(hz, 8.0/3) * pow(m500.solarMass() / 3e14, 2.0/3+0.12) * 8.4;

  double val1 = 1.65e-03 * pow(hz, 8.0/3) * pow(m500.solarMass() / 3e14, 2.0/3);
  double val2 = 1.45e-11 * pow(cosmo.E(), 8.0/3) * pow(m500.solarMass() / (1e15 * cosmo.dimensionlessHubbleConstant()), 2.0/3);
  Pressure aP500, nP500;

  aP500.setKeVPerCm3(val1);
  nP500.setErgPerCm3(val2);

  Pressure arnaudPressure;
  arnaudPressure = aP500 * 8.4;

  Pressure nagaiPressure;
  nagaiPressure  = nP500 * 3.3;

  COUT("Arnaud P500 = " << aP500.keVPerCm3());
  COUT("Nagai  P500 = " << nP500.keVPerCm3());

  COUT("hz = " << hz << " DA = " << cosmo.angularDiameterDistance().Gpc() << " Pressure is now: " << arnaudPressure.keVPerCm3() << " keV/cm^3");

  double y = arnaudPressure / norm * theta.radians();

  COUT("y = " << y);

  theta.setArcMinutes(3.6);
  double p0 = (-2100e-6/2.73)/(-2) * norm.keVPerCm3() / theta.radians();
  COUT("y of 3.7e-4 implies P0 = " << p0);

  double apf = 8.403 * 1.65e-3 * pow(1.089, 8.0/3);

  nagaiPressure.setErgPerCm3(3.3   * 1.45e-11 * pow(cosmo.E(), 8.0/3));

  COUT("Arnaud prefactors = " << apf);
  COUT("Nagai  prefactors = " << nagaiPressure.keVPerCm3());

  double exa = log(p0/apf) / (2.0/3 + 0.12);
  double msolara = exp(exa) * 3e14;

  double exn = log(p0/nagaiPressure.keVPerCm3()) / (2.0/3 + 0.12);
  double msolarn = exp(exn) * 1e15 / cosmo.dimensionlessHubbleConstant();

  COUT("\nMy norm (" << p0 << " keVPerCm3) implies M500 = " << msolara << " for Arnaud normalization");
  COUT("\nMy norm (" << p0 << " keVPerCm3) implies M500 = " << msolarn << " for Nagai normalization");

  //------------------------------------------------------------
  // Test Tony's chain
  //------------------------------------------------------------

  Pressure tonyNorm;
  tonyNorm.setErgPerCm3(2e-10);

  p0 = tonyNorm.keVPerCm3();

  exn = log(p0/nagaiPressure.keVPerCm3()) / (2.0/3 + 0.12);
  msolarn = exp(exn) * 1e15 / cosmo.dimensionlessHubbleConstant();

  COUT("\nTony's norm (" << p0 << " keVPerCm3) implies M500 = " << msolarn << " for the Nagai07 model");

  //------------------------------------------------------------
  // Now do another calculation
  //------------------------------------------------------------

  Temperature temp;
  temp.setKeV(Program::getDoubleParameter("Te"));

  Length R;
  R.setMegaParsec(Program::getDoubleParameter("R"));

  Angle thetaCore;
  thetaCore.setArcSec(20);

  Energy thermalEnergy;
  thermalEnergy = temp;

  Length dA = cosmo.angularDiameterDistance();
  Area a = dA * dA;

  double rat = a / Constants::sigmaT_;

  COUT(std::endl << "restmass/thermal = " << restMassEnergy/thermalEnergy);
  COUT("dA^2/sigmaT      = " << rat);
  COUT("m_p = " << Constants::protonMass_.solarMass());

  double trad = thetaCore.radians();
  double t2 = trad * trad;

  double x = (R/dA) / trad;
  double y0 = Program::getDoubleParameter("y");

  COUT(std::endl << "Mass = " 
       <<  1.18 * Constants::protonMass_.solarMass() * y0 * (restMassEnergy/thermalEnergy) * rat * t2 * 4 * (x - atan(x)));

  //------------------------------------------------------------
  // Now do another calculation
  //------------------------------------------------------------

  //  temp.setKeV(1.0);

  thetaCore.setArcSec(1);

  thermalEnergy = temp;

  dA.setGigaParsec(1.0);
  a = dA * dA;

  rat = a / Constants::sigmaT_;

  COUT(std::endl << "restmass/thermal = " << restMassEnergy/thermalEnergy);
  COUT("dA^2/sigmaT      = " << rat);
  COUT("m_p = " << Constants::protonMass_.solarMass());

  trad = thetaCore.radians();
  t2 = trad * trad;

  COUT(std::endl << "Mass scaling = " <<  1.18 * Constants::protonMass_.solarMass() * (restMassEnergy/thermalEnergy) * rat * t2 * 4);

  COUT(std::endl << "Y to Mass scaling = " <<  1.18 * Constants::protonMass_.solarMass() * (restMassEnergy/thermalEnergy) / Constants::sigmaT_.squaredMpc() << " Msolar/Mpc^2");

  return 0;
}

Mass arnaudM500(Pressure& p, double z)
{
  HubbleConstant H0;
  H0.setH100(0.7);

  Cosmology cosmo;
  cosmo.setRedshift(z);
  cosmo.setH0(H0);
  cosmo.setOmegaM(0.3);
  cosmo.setOmegaL(0.7);

  double hz = cosmo.dimensionlessHubbleConstant();

  double prat = p.keVPerCm3() / (5 * 1.65e-3 * pow(hz, 8.0/3) * 8.403);
  double lnprat = log(prat) / (2.0/3 + 0.12);

  Mass m500;

  m500.setSolarMass(exp(lnprat) * 3e14);

  return m500;
}

Mass nagaiM500(Pressure& p, double z)
{
  HubbleConstant H0;
  H0.setH100(0.7);

  Cosmology cosmo;
  cosmo.setRedshift(z);
  cosmo.setH0(H0);
  cosmo.setOmegaM(0.3);
  cosmo.setOmegaL(0.7);
  COUT("DA = " << cosmo.angularDiameterDistance().Mpc());

  double Ez = cosmo.E();
  double h = 0.7;

  double prat = p.ergPerCm3() / (5 * 1.45e-11 * pow(Ez, 8.0/3) * 3.3);
  double lnprat = log(prat) / (2.0/3);

  Mass m500;

  COUT("h = " << h);
  m500.setSolarMass(exp(lnprat) * 1e15 * h);

  

  return m500;
}

Pressure nagaiP500(Mass& m500, double z)
{
  HubbleConstant H0;
  H0.setH100(0.7);

  Cosmology cosmo;
  cosmo.setRedshift(z);
  cosmo.setH0(H0);
  cosmo.setOmegaM(0.3);
  cosmo.setOmegaL(0.7);
  COUT("DA = " << cosmo.angularDiameterDistance().Mpc());

  double Ez = cosmo.E();
  double h = 0.7;

  double mrat = m500.solarMass() / 1e15 * h;

  Pressure p500;
  p500.setErgPerCm3(1.45e-11 * pow(mrat, 2.0/3) * pow(Ez, 8.0/3));

  return p500;
}

Pressure arnaudP500(Mass& m500, double z)
{
  HubbleConstant H0;
  H0.setH100(0.7);

  Cosmology cosmo;
  cosmo.setRedshift(z);
  cosmo.setH0(H0);
  cosmo.setOmegaM(0.3);
  cosmo.setOmegaL(0.7);
  COUT("DA = " << cosmo.angularDiameterDistance().Mpc());

  double hz = cosmo.h();
  double h70 = 1.0;

  double mrat = m500.solarMass() / 3e14 * h70;

  Pressure p500;
  p500.setKeVPerCm3(1.65e-3 * pow(mrat, 2.0/3) * pow(hz, 8.0/3) * h70 * h70);

  return p500;
}

Pressure arnaudCharP500(Mass& m500, double z)
{
  HubbleConstant H0;
  H0.setH100(0.7);

  Cosmology cosmo;
  cosmo.setRedshift(z);
  cosmo.setH0(H0);
  cosmo.setOmegaM(0.3);
  cosmo.setOmegaL(0.7);
  COUT("DA = " << cosmo.angularDiameterDistance().Mpc());

  double hz = cosmo.h();
  double h = 0.7;

  double mrat = m500.solarMass() / 3e14;

  Pressure p500;
  p500.setKeVPerCm3(1.65e-3 * pow(mrat, 2.0/3) * pow(hz, 8.0/3));

  return p500;
}

double arnaudPrefactor()
{
  HubbleConstant H0;
  H0.setKmPerSecPerMpc(70.0);

  Mass msolar;
  msolar.setSolarMass(1.0);

  double G = Constants::gravitationalConstantCgs_;
  
  Pressure press1, press2;

  press1.setErgPerCm3(4197.54 * pow(H0.inverseSeconds(), 8.0/3) * pow(3e14, 2.0/3) * pow(msolar.g(), 2.0/3));

  press2.setErgPerCm3(4197.54 * pow(H0.inverseSeconds(), 8.0/3) * pow(1e15, 2.0/3) * pow(msolar.g(), 2.0/3));

  COUT("press2 = " << press2.ergPerCm3());

  return press1.keVPerCm3();
}
