#include <iostream>
#include <iomanip>

#include <cmath>

#include <errno.h>

#include "gcp/program/Program.h"

#include "gcp/util/String.h"
#include "gcp/util/Exception.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "debuglevel", "0", "i", "Debug level"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  String str("this is another test Mpc^2");

  COUT("str = '" << str << "'");
  str.replace("^", "\\\\^");
  COUT("str = '" << str << "'");

  return 0;
}
