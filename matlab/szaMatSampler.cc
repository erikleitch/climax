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
  
  MexParser nParser(prhs[0]);
  MexParser xParser(prhs[1]);
  MexParser yParser(prhs[2]);

  unsigned n   = (unsigned)(*nParser.getDoubleData());
  double*  x   = xParser.getDoubleData();
  double*  y   = yParser.getDoubleData();

  Sampler sampler;
  unsigned length = xParser.getNumberOfElements();

  COUT("Setting lengtyh to: " << length << " n = " << n);

  sampler.setdYdX(length, x, y);

  std::vector<double> samples = sampler.generateSamples(n);

  double* ret = MexHandler::createDoubleArray(&plhs[0], n);

  for(unsigned i=0; i < n; i++) {
    if(i < 10) {
      COUT(samples[i]);
    }
    ret[i] = samples[i];
  }

  return;
}
