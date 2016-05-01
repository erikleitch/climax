#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Constants.h"
#include "gcp/util/Cosmology.h"
#include "gcp/util/HubbleConstant.h"
#include "gcp/util/Speed.h"

#include "gcp/pgutil/PgUtil.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "rad",        "0", "d", "Radians"},
  { "z",          "0.2", "d", "redshift"},
  { "kpc",       "70", "d", "Distance in kpc"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  HubbleConstant hc;
  hc.setKmPerSecPerMpc(100);
  Length d = Constants::lightSpeed_ / hc;

  Length size(Length::KiloParsec(), Program::getDoubleParameter("kpc"));
  double z = Program::getDoubleParameter("z");

  COUT("Distance is: " << d.gigaParsec() << " Gpc");

  Cosmology cosmo;
  cosmo.setH0(hc);
  
  COUT("Hubble constant at z = " << z << " 0: "    << cosmo.H(z));
  
  // Now calculate some distance measures

  hc.setKmPerSecPerMpc(72.0);
  cosmo.setH0(hc);
  cosmo.setOmegaL(0.7);
  cosmo.setOmegaM(0.3);

  COUT("AD(" << z << ") = " << cosmo.angularDiameterDistance(z).megaParsec() << " Mpc");

  Angle asize;
  asize.setRadians(size.Kpc() / cosmo.angularDiameterDistance(z).Kpc());
  COUT("size(" << size.Kpc() << " kpc at z = " << z << ": " << asize.arcmin() << " arcmin");

  COUT("AD(" << z << ") = " << cosmo.angularDiameterDistance(z).lightYear()  << " light-years");
  COUT("LD(" << z << ") = " << cosmo.luminosityDistance(z).lightYear()  << " light-years");
  COUT("LTD(" << z << ") = " << cosmo.lightTravelDistance(z).lightYear()  << " light-years");
  COUT("CD(" << z << ") = " << cosmo.comovingDistance(z).lightYear()  << " light-years");

  // Let's calculate a typical pressure, just for kicks:

  double y = 1.5e-5;
  double c = Constants::lightSpeed_.centimetersPerSec(); 
  double mc2 = Constants::electronMassCgs_ * c * c;
  //  double da = 3.7e27; // Typical DA in cm
  double da = 3.7e24; // Typical DA in cm

  double presscgs = (y * mc2) / (Constants::sigmaT_.squaredCentimeters() * da);
  Pressure pressure;
  pressure.setPascal(presscgs * 1e-3 / 1e-2); // To convert (g -> kg / cm -> m)
  Temperature temp;
  temp.setK(5e7);
  
  COUT("Pressure       = " << pressure);
  NumberDensity nd = pressure/temp;
  COUT("Number density = " << nd.inverseCubicCentimeters() << " cm-3");
  COUT("Number density = " << nd.inverseCubicMeters() << " m-3");

  // Setup the calculation from the other side

  nd.setInverseCubicCentimeters(1e-3);

  pressure = nd * temp;
  COUT("Pressure for nd = " << nd << " would be " << pressure);

  nd = pressure / temp;
  COUT("ND = " << nd);

  pressure = nd * temp;
  COUT("Pressure for nd = " << nd << " would be " << pressure);


  unsigned nz = 1000;
  double zmax = 10000, zmin = 0.001;
  double dlz = (log10(zmax) - log10(zmin))/nz;

  std::vector<double> zs(nz);
  std::vector<double> distance(nz);

  PgUtil::setXmin(1e-4);
  PgUtil::setXmax(1e5);
  PgUtil::setYmin(0.001);
  PgUtil::setYmax(10000);
  PgUtil::setUsedefs(true);

  PgUtil::setTraceColor(5);

  for(unsigned i=0; i < nz; i++) {
    zs[i] = pow(10, log10(zmin) + i*dlz);
    distance[i] = cosmo.angularDiameterDistance(zs[i]).lightYear()/1e9;
  }

  PgUtil::setLogPlot(true);
  PgUtil::linePlot(zs, distance);
  PgUtil::setOverplot(true);

  PgUtil::setTraceColor(2);

  for(unsigned i=0; i < nz; i++) {
    zs[i] = pow(10, log10(zmin) + i*dlz);
    distance[i] = cosmo.luminosityDistance(zs[i]).lightYear()/1e9;
  }

  PgUtil::linePlot(zs,distance);

  PgUtil::setTraceColor(3);

  for(unsigned i=0; i < nz; i++) {
    zs[i] = pow(10, log10(zmin) + i*dlz);
    distance[i] = cosmo.comovingDistance(zs[i]).lightYear()/1e9;
  }

  PgUtil::linePlot(zs,distance);

  PgUtil::setTraceColor(4);

  for(unsigned i=0; i < nz; i++) {
    zs[i] = pow(10, log10(zmin) + i*dlz);
    distance[i] = cosmo.lightTravelDistance(zs[i]).lightYear()/1e9;
  }

  PgUtil::linePlot(zs,distance);

  
  return 0;
}
