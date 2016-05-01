#include <iostream>
#include <iomanip>
#include <unistd.h>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Mass.h"

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
  std::string path;
  path.resize(PATH_MAX);

  getcwd(&path[0], PATH_MAX);

  COUT("Cwd = " << path);
  return 0;
}
