#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Timer.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "n",        "100000", "i", "Number of iterations"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

void doSomething(unsigned i, double& result)
{
  result += pow(i, 2.0);
}

double doSomething(unsigned i)
{
  return pow(i, 2.0);
}

void fn1(unsigned n)
{
  double result=0.0;

  for(unsigned i=0; i < n; i++) {
    result += pow(i, 2.0);
  }

  COUT("Result = " << result);
}

void fn2(unsigned n)
{
  double result=0.0;

  for(unsigned i=0; i < n; i++) {
    doSomething(i, result);
  }

  COUT("Result = " << result);
}

void fn3(unsigned n)
{
  double result=0.0;

  for(unsigned i=0; i < n; i++) {
    result += doSomething(i);
  }

  COUT("Result = " << result);
}

int Program::main()
{
  Timer timer;

  timer.start();
  fn1(Program::getIntegerParameter("n"));
  timer.stop();
  COUT(timer.deltaInSeconds());

  timer.start();
  fn2(Program::getIntegerParameter("n"));
  timer.stop();
  COUT(timer.deltaInSeconds());

  timer.start();
  fn3(Program::getIntegerParameter("n"));
  timer.stop();
  COUT(timer.deltaInSeconds());

  return 0;
}
