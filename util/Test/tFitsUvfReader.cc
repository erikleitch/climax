#include <stdio.h>
#include <iostream>
#include <vector>
#include <iostream>

#include "gcp/util/Directives.h"

#include "gcp/util/FitsUvfReader.h"

#include "gcp/program/Program.h"

#include "gcp/util/Coordinates.h"
#include "gcp/util/Vector.h"

#include "fitsio.h"
#include "cpgplot.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "file",      "test.fits", "s", "The input file"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

unsigned plot(std::vector<double>& x, std::vector<double>& y);

int Program::main()
{
  gcp::util::FitsUvfReader reader(Program::getStringParameter("file"));
  COUT("File contains " << reader.nGroup() << " groups");

  FitsUvfReader::Vis vis;

  COUT("Here main 0");

  std::vector<double> u;
  u.resize(reader.nGroup());
  std::vector<double> v;
  v.resize(reader.nGroup());

  //  for(unsigned i=1; i < reader.nGroup(); i++) {

  COUT("Here main 1");

  for(unsigned i=1; i < reader.nGroup(); i++) {
    reader.readData(i, vis);
    u[i] = vis.u_;
    v[i] = vis.v_;
  }

  COUT("OBject is: " << reader.object());

  plot(u, v);

  return 0;
}

unsigned plot(std::vector<double>& x, std::vector<double>& y)
{
  std::vector<float> xvec;
  std::vector<float> yvec;

  xvec.resize(x.size());
  yvec.resize(y.size());

  for(unsigned i=0; i < x.size(); i++) {
    xvec[i] = x[i];
    yvec[i] = y[i];
  }

  float xmin, xmax;
  float ymin, ymax, ymean;
  bool first = true;

  for(unsigned i=0; i < yvec.size(); i++) {

    if(first) {

      xmin = xvec[i];
      xmax = xvec[i];
      ymin = yvec[i];
      ymax = yvec[i];

      ymean = 0.0;

      first = false;
    }

    xmin = xvec[i] < xmin ? xvec[i] : xmin;
    ymin = yvec[i] < ymin ? yvec[i] : ymin;

    xmax = xvec[i] > xmax ? xvec[i] : xmax;
    ymax = yvec[i] > ymax ? yvec[i] : ymax;

    ymean += (yvec[i] - ymean)/(i+1);
    
  }

  first = true;
  for(unsigned i=0; i < yvec.size(); i++) {
    //    yvec[i] = log(yvec[i]);

    if(first) {
      ymin = yvec[i];
      ymax = yvec[i];
      first = false;
    }

    ymin = yvec[i] < ymin ? yvec[i] : ymin;
    ymax = yvec[i] > ymax ? yvec[i] : ymax;
  } 

  float xrange = xmax - xmin;
  float yrange = ymax - ymin;

  // pgplot stuff

  if(cpgbeg(0,"?",1,1)!=1)
    return 1;

  cpgswin(xmin-0.1*xrange, xmax+0.1*xrange, ymin-0.1*yrange, ymax+0.1*yrange);
  cpgwnad(xmin-0.1*xrange, xmax+0.1*xrange, ymin-0.1*yrange, ymax+0.1*yrange);

  cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
  cpgpt(yvec.size(), &xvec[0], &yvec[0], 1);

  return 0;
}
