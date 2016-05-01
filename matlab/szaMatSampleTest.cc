/**.......................................................................
 * MATLAB Mex file for calculating UVW, given an array of antenna
 * locations, a source declination and a list of HA
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
  
  MexParser xParser(prhs[0]);
  MexParser yParser(prhs[1]);

  double* x = xParser.getDoubleData();
  double* y = yParser.getDoubleData();
  unsigned nSamp = (unsigned)(*MexParser::getDoubleData(prhs[2]));
    
  double* s = MexHandler::createDoubleArray(&plhs[0], nSamp);

  Sampler sampler;


  COUT("Here 0");
  sampler.setdYdX(xParser.getNumberOfElements(), x, y);

  COUT("Here 1");
  std::vector<double> samps = sampler.generateSamples(nSamp);

  COUT("Here 2");
  for(unsigned i=0; i < nSamp; i++)
    s[i] = samps[i];

  COUT("Here 3");
  return;
}
