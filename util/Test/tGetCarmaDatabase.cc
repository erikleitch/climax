#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/CurlUtils.h"
#include "gcp/util/Exception.h"
#include "gcp/util/String.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "src",        "", "s", "Source to fetch"},
  { "print",        "f", "b", "True to print record"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

void getIt(std::string prefix)
{
  std::ostringstream os;

  os.str("");
  os << "curl http://www.ovro.caltech.edu/~eml/carma/" << prefix << ".txt > " << prefix << ".txt";
  system(os.str().c_str());
}
int Program::main()
{
  getIt("allprojs");
  getIt("2007");
  getIt("2008");
  getIt("2009");
  getIt("2010");
  getIt("2011");
  getIt("2012");
  getIt("2013");
  getIt("2014");

  return 0;
}

