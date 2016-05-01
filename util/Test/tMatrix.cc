#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Matrix.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "rad",        "0", "d", "Radians"},
  { "ixpix",        "0", "d", "Radians"},
  { "iypix",        "0", "d", "Radians"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  Matrix<double> m(3,3);

  m[0][0] = 1.0; m[0][1] = 0.1; m[0][2] = 0.1;
  m[1][0] = 0.1; m[1][1] = 1.0; m[1][2] = 0.1;
  m[2][0] = 0.1; m[2][1] = 0.1; m[2][2] = 1.0;

  COUT("det = " << m.determinant());

  Matrix<double> wcs(2,2);

  wcs[0][0] = -1.0744180904515E-6; wcs[0][1] = -9.2477152608260E-5;
  wcs[1][0] = -9.2472923034340E-5; wcs[1][1] = 1.11903810112625E-6;

  Vector<double> vec(2);
  vec[0] = -348;
  vec[1] = 3020;

  //  COUT(wcs * vec);

  double ixpix = Program::getDoubleParameter("ixpix");
  double iypix = Program::getDoubleParameter("iypix");

  vec[0] = (ixpix + 348);
  vec[1] = (iypix - 3020);

  COUT(wcs * vec);

  return 0;
}
