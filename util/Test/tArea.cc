#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Area.h"
#include "gcp/util/Constants.h"
#include "gcp/util/Energy.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/Volume.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Debug.h"
#include "gcp/util/NumberDensity.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "rad",        "0", "d", "Radians"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  Volume vol;
  vol.setCubicCentimeters(1000);

  Area area;
  area.setSquaredCentimeters(100);

  Length length;
  length.setCentimeters(10);

  COUT(vol);
  COUT(vol/area);
  COUT(vol/length);
  COUT(1.0/vol);

  COUT(area);
  COUT(area/length);

  COUT(area*length);
  COUT(length*area);

  COUT(length*length);
  COUT(length*length*length);

  NumberDensity nd;
  nd.setInverseCubicCentimeters(13.0/1000);

  COUT(nd * vol);

  Energy energy;
  energy = Constants::electronMass_;
  COUT("ergs = " << energy.ergs());
  COUT(energy.eV() << " eV");

  Wavelength wave;
  wave.setAngstroms(10);
  COUT("Wave = " << wave);
  Frequency freq;
  freq = wave;
  COUT("Freq = " << freq);
  energy = wave;
  COUT("Energy = " << energy.eV() << " eV");
  
  area.setUnits("Mpc^2");
  area.setSquaredMpc(1.0);
  COUT("Area is now: mpc = " << area.squaredMpc() << " val = " << area.value() << " unitval = " << area.getUnitVal());

  return 0;
}
