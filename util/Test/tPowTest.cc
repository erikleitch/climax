#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Debug.h"
#include "gcp/util/Timer.h"

#include "gcp/models/fastonebigheader.h"

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
  double x, p, y;
  Timer timer;
  double delta = 0.0;

  unsigned iMax = 5e9;
  timer.start();

  for(unsigned i=0; i < iMax; i++) {
    x = ((double)i)/ iMax;
    p = 0.8;

    //    y = fastpow(x, p);
    y = pow(x, p);
      //    COUT("x = " << x << " y = " << y);
    //    y = exp(p*log(x));
  }
  timer.stop();

  COUT("delta = " << timer.deltaInSeconds());

  return 0;
}
