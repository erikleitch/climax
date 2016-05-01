#include <iostream>
#include <sstream>
#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Debug.h"
#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/Astrometry.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

void printDate();

KeyTabEntry Program::keywords[] = {
  { "ra1",        "00:00:00", "s", "Hours"},
  { "ra2",        "00:00:00", "s", "Hours"},
  { "dec1",       "45:00:00", "s", "Degress"},
  { "dec2",       "45:00:00", "s", "Degress"},
  {END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

void testCoordConversions();
void testAngularSep(HourAngle& ra1, Declination& dec1, HourAngle& ra2, Declination& dec2);

int Program::main()
{
  //  testCoordConversions();

  HourAngle ra1, ra2;
  Declination dec1, dec2;

  ra1.setHours(Program::getStringParameter("ra1"));
  ra2.setHours(Program::getStringParameter("ra2"));
  dec1.setDegrees(Program::getStringParameter("dec1"));
  dec2.setDegrees(Program::getStringParameter("dec2"));

  testAngularSep(ra1, dec1, ra2, dec2);

  return 0;
}

void testAngularSep(HourAngle& ra1, Declination& dec1, HourAngle& ra2, Declination& dec2)
{
  Astrometry astro;

  Angle sep = astro.angularSeparation(ra1, dec1, ra2, dec2);
  COUT("sep = " << sep);

  return;
}

void testCoordConversions()
{
  Astrometry astro;
  TimeVal mjd;

  mjd.setMjd(53627.0);

  HourAngle mRa,  aRa;
  Declination  mDec, aDec;

  mRa.setHours(20.3043);
  mDec.setDegrees(50.546);

  COUT("Mean RA:     " << mRa << " Mean Dec:     " << mDec);

  astro.j2000ToApparentPlace(mRa, mDec, mjd, aRa, aDec);

  COUT("Apparent RA: " << mRa << " Apparent Dec: " << mDec);

  astro.apparentToJ2000Place(aRa, aDec, mjd, mRa, mDec);

  COUT("Mean RA:     " << mRa << " Mean Dec:     " << mDec);

  // 14.2569 Hours
  // 36.1782 Degrees
  // MJD 53611.08126758102

  aRa.setHours(14.2569);
  aDec.setDegrees(36.1782);
  mjd.setMjd(53611.08126758102);

  COUT("Apparent RA: " << aRa << " Apparent Dec: " << aDec);

  astro.apparentToJ2000Place(aRa, aDec, mjd, mRa, mDec);

  COUT("Mean RA:     " << mRa << " Mean Dec:     " << mDec);
  
}

