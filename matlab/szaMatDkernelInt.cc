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
#include "gcp/matlab/MlPtsrc.h"

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
  
  MexParser     fluxParser(prhs[0]);
  MexParser    gammaParser(prhs[1]);
  MexParser    sigmaParser(prhs[2]);
  MexParser     nsigParser(prhs[3]);
  MexParser     nbinParser(prhs[4]);

  double* fluxes    = fluxParser.getDoubleData();
  unsigned nsrc     = fluxParser.getNumberOfElements();

  double gamma      = *(gammaParser.getDoubleData());
  double sigma      = *(sigmaParser.getDoubleData());
  double nsig       = *(nsigParser.getDoubleData());

  unsigned nbin     = (unsigned)(*(nbinParser.getDoubleData()));

  double* psi = MexHandler::createDoubleArray(&plhs[0], nsrc);

  MlPtsrc::dkernelint(nsrc, fluxes, psi, gamma, sigma, nsig, nbin);

  return;
}
