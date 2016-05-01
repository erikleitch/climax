#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Chisq.h"
#include "gcp/util/ChisqVariate.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Sampler.h"

#include "gcp/pgutil/PgUtil.h"

#include <vector>

using namespace std;
using namespace gcp::util;
using namespace gcp::program;


KeyTabEntry Program::keywords[] = {
  { "ntrial",     "1000", "i", "Number of trials"},
  { "nsamp",      "10", "i", "Number of samples per trial"},
  { "chisq",      "42", "d", "Chisquared value to test"},
  { "ndof",       "51", "i", "Number of dof"},
  { "seed",       "1",  "i", "Seed"},
  { "nhist",      "10", "i", "number of histogram bins"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

void pte(double chisq, unsigned ndof)
{
  ChisqVariate cv;
  cv.setChisq(chisq, ndof);
  COUT("PTE = " << cv.samplingPte() << " pdf = " << cv.pdf() << " ln = " << cv.pdf().lnValue() << " " 
       << ((double)(ndof)/2-1)*log(chisq) - chisq/2);

  double val1 = (((double)(53312)/2-1)*log(54541.0) - 54541.0/2);
  double val2 = (((double)(53312)/2-1)*log(54540.0) - 54540.0/2);

  COUT("diff = " << val1 - val2);

}


int Program::main()
{
  unsigned nDof = Program::getIntegerParameter("ndof");
  
  double chisqmin = 0.1;
  double chisqmax = 3.0*nDof;
  double dchisq = (chisqmax - chisqmin) / 100;

  std::vector<double> xvals(100);
  std::vector<double> yvals(100);

  for(unsigned i=0; i < 100; i++) {
    double chisq = chisqmin + dchisq * i;
    xvals[i] = chisq;
    yvals[i] = Sampler::chisqPdf(chisq, nDof);
  }

  PgUtil::linePlot(xvals, yvals);

  return 0;
}
