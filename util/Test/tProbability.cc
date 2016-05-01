#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Exception.h"
#include "gcp/util/Probability.h"

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
  Probability p1;
  p1.setValue(0.0);

  Probability p2;
  p2.setValue(0.5);

  Probability p3;
  p3.setValue(0.5);

  COUT("p1 " << p1);
  COUT("p2 " << p2);
  COUT("p3 " << p3);

  COUT("p1*p2 " << p1*p2);
  COUT("p2*p3 " << p2*p3);

  COUT("p1/p2 " << p1/p2);
  COUT("p2/p3 " << p2/p3);
  COUT("p2/p1 " << p2/p1);

  return 0;
}
