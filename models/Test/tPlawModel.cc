#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Timer.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/models/PowerlawProfile.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "d",        "0", "d", "distance at which the 3D model is centered"},
  { "r",        "0", "d", "cylindrical radius at which to evaluate the line integral"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

void plotModel(double fac, double r1);

int Program::main()
{
  PgUtil::open("/xs");

  plotModel(0.0, 1.0);

  PgUtil::setOverplot(true);
  PgUtil::setWin(false);

  plotModel(0.0, 0.2);
  plotModel(1.0, 0.2);
  plotModel(0.2, 0.2);


  return 0;
}

void plotModel(double fac, double r1)
{
  PowerlawProfile   model;

  model.getVar("thetaCore")->setVal(2.0, "\'");
  model.getVar("thetaCore")->wasSpecified_ = true;

  model.getVar("Sradio")->setVal(1.0, "mK");
  model.getVar("Sradio")->setUnits("mK");
  model.getVar("Sradio")->wasSpecified_ = true;

  model.setParameter("n", "1", "");

  model.getVar("a0")->setVal(2.0, "");
  model.getVar("a0")->wasSpecified_ = true;

  model.getVar("fac")->setVal(fac, "");
  model.getVar("fac")->wasSpecified_ = true;

  model.getVar("a1")->setVal(1.0, "");
  model.getVar("a1")->wasSpecified_ = true;

  model.getVar("r1")->setVal(r1, "");
  model.getVar("r1")->wasSpecified_ = true;

  model.checkSetup();

  unsigned n = 100;
  double xmin = 0.1;
  double xmax = 4.0;
  double dx = (xmax - xmin)/(n-1);

  std::vector<double> x(n);
  std::vector<double> y(n);

  for(unsigned i=0; i < n; i++) {
    x[i] = dx * i;
    y[i] = model.radialRadioModelEml(x[i], 0);
  }

  PgUtil::linePlot(x,y, "", "", "", true);
}
