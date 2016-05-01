/**.......................................................................
 * MATLAB Mex file for gridding a population
 * of sources
 *
 * Use like:
 *
 * d=gcpMatSampleTest(x(double[]), dydx(double[]), nSamp(int))
 *
 */
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"
#include "gcp/array/code/share/slalib/slalib.h"
#include "gcp/translator/DelayEngineNormal.h"
#include "gcp/array/code/share/include/rtcnetcoms.h"

#include "gcp/util/PtSrcGen.h"
#include "gcp/util/Sampler.h"

#include "mex.h"
#include "matrix.h"

#include <iostream.h>
#include <math.h>

using namespace std;
using namespace gcp::util;
using namespace gcp::matlab;

/**.......................................................................
 * Entry point from the matlab environment
 */
void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  gcp::util::Logger::installStdoutPrintFn(MexHandler::stdoutPrintFn);
  gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);
  
  gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
  gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
  gcp::util::ErrHandler::installLogFn(MexHandler::logFn);

  std::ostringstream usage;

  usage << "Usage: " << std::endl;
  usage << std::endl;
  usage << "   [ps] = gcpMatPoissPdf(ks, lambda)" << std::endl;
  usage << std::endl;
  usage << " where:" << std::endl;
  usage << std::endl;
  usage << "    ks     -- (double) arbitrary-size array of integer values" << std::endl;
  usage << "    lambda -- (double) the mean of the distribution" << std::endl;
  usage << std::endl;

  if(nrhs < 2) {
    ThrowError("Wrong number of arguments: " << usage.str());
  }

  // Push the top-level mxArray into a parser

  MexParser kParser(prhs[0]);
  MexParser lParser(prhs[1]);

  // Extract point source population parameters

  double l = *lParser.getDoubleData();
  double* kVals  = kParser.getDoubleData();

  double* pdfVals = MexHandler::createDoubleArray(&plhs[0], kParser);

  for(unsigned i=0; i < kParser.getNumberOfElements(); i++) {
    pdfVals[i] = Sampler::poissPdf((unsigned)kVals[i], l);
  }

  return;
}
