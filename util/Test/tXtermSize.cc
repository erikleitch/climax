#include <iostream>
#include <iomanip>
#include <termios.h>
#include <sys/ioctl.h>
#include <stdio.h>
#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Debug.h"

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
  XtermManip xtm;

  COUT("Background is " << xtm.getBg());
  COUT("Foreground is " << xtm.getFg());

  return 0;
}
