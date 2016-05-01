#include <iostream>
#include <vector>
#include <cmath>

#include "cpgplot.h"

#include "gcp/program/Program.h"

#include "gcp/fftutil/Dft1d.h"

#include "gcp/util/Exception.h"
#include "gcp/util/Stats.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

void Program::initializeUsage() {};

KeyTabEntry Program::keywords[] = {
  { "amp",     "1",                         "d", "amplitude of the random noise component"},
  { "floor",   "0.0",                       "d", "amplitude of the offset"},
  { "n",       "100",                       "i", "n"},
  { "manual",  "f",                         "b", "manual or auto?"},
  { "av",      "t",                         "b", "average?"},
  { "apod",    "f",                         "b", "apodize?"},
  { "nav",     "10",                        "i", "number of transforms to average"},
  { "period1", "20",                        "i", "period of the first component"},
  { "period2", "40",                        "i", "period of the second component"},
  { "ratio",   "1",                         "d", "ratio of the amplitudes of the two components"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

unsigned plot(std::vector<double>& yvec, double dx=1.0);

double apodizationCoefficient(unsigned n, unsigned iSamp)
{
  // For now just default to triangle apodization

  double iVal    = (double)iSamp;
  double halfVal = (double)(n)/2;

  if(iVal < halfVal) {
    return 2.0/(n - 1) * iVal;
  } else {
    return 2.0/(n - 1) * (n - 1 - iSamp);
  }

}

int Program::main()
{
  double amp       = Program::getDoubleParameter("amp");
  double floor     = Program::getDoubleParameter("floor");
  double ratio     = Program::getDoubleParameter("ratio");
  unsigned n       = Program::getIntegerParameter("n");
  unsigned nav     = Program::getIntegerParameter("nav");
  bool manual      = Program::getBooleanParameter("manual");
  bool av          = Program::getBooleanParameter("av");
  bool apod        = Program::getBooleanParameter("apod");
  unsigned period1 = Program::getIntegerParameter("period1");
  unsigned period2 = Program::getIntegerParameter("period2");

  TimeVal tVal;
  tVal.setMilliSeconds(10);

  Dft1d dft(n, true);

  dft.setTimeRes(tVal);

  dft.setAverage(av);

  if(apod)
    dft.setApodizationType(Dft1d::APOD_TRIANGLE);

  nav = av ? nav : 1;

  std::vector<double> input(n);

  for(unsigned iav=0; iav < nav; iav++) {
    for(unsigned i=0; i < n; i++) {
      input[i] = Sampler::generateGaussianSample(1.0);
      dft.pushSample(input[i]/n);
    }
  }  

  plot(input);

  std::vector<double> y = dft.abs();
  COUT("Mean = " << Stats::mean(y));
  COUT("Sum  = " << Stats::sum(y));
  COUT("y[0] = " << y[0]);
  COUT("y[1] = " << y[1]);

  plot(y, dft.getFrequencyResolution().Hz());

  
}

unsigned plot(std::vector<double>& input, double dx)
{
  std::vector<float> xvec;
  std::vector<float> yvec;

  xvec.resize(input.size());
  yvec.resize(input.size());

  for(unsigned i=0; i < input.size(); i++)
    yvec[i] = input[i];

  float xmin, xmax;
  float ymin, ymax, ymean;
  bool first = true;

  for(unsigned i=0; i < yvec.size(); i++) {

    xvec[i] = i * dx;;

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
  cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
  cpgline(yvec.size(), &xvec[0], &yvec[0]);

  return 0;
}
