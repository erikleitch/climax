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
  { "val", "0", "i", "Value to print"},
  { "n",   "0", "i", "nchar"},
  { "str", " ", "s", "String to print"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  std::string str("~/projects/carma/madcows/reduc/reduc_19sep/W0010+31/WB/src_W0010+31_20140907_050452_D_WB.uvf");
  std::string::iterator start = str.begin();
  std::string::iterator stop  = start + 1;

  str.replace(start, stop, "eml");
  COUT("str = " << str);

  String Str("~/projects/carma/madcows/reduc/reduc_19sep/W0010+31/WB/src_W0010+31_20140907_050452_D_WB.uvf");
  Str.expandTilde();

  COUT("Str = " << Str);

  return 0;
}
