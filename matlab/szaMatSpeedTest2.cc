/**.......................................................................
 * MATLAB Mex file for calculating UVW, given an array of antenna
 * locations, a source declination and a list of HA
 *
 * Use like:
 *
 * d=gcpMatCalcUvw(lat, antLoc(ENU), refLoc(ENU), dec(rad), HA(hours))
 *
 */
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"
#include "gcp/array/code/share/slalib/slalib.h"
#include "gcp/util/Coordinates.h"

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

  double* xprof   = MexParser::getDoubleData(prhs[0]);
  double* yprof   = MexParser::getDoubleData(prhs[1]);
  double* length  = MexParser::getDoubleData(prhs[2]);
  double* minSep  = MexParser::getDoubleData(prhs[3]);

  double dxprof = (xprof[1] - xprof[0]);
  unsigned nprof = MexParser::getNumberOfElements(prhs[0]);
  
  // The theoretical uvr min/max
  
  double uvrGlobalMin = *minSep;
  double uvrGlobalMax = (*length) * sqrt(2.0);

  uvrGlobalMin = uvrGlobalMin < xprof[0] ? uvrGlobalMin : xprof[0];
  uvrGlobalMax = uvrGlobalMax > xprof[nprof-1] ? uvrGlobalMax: xprof[nprof-1];
  
  // Calculate additional elements we need to add to the profile
  // length to encompass the global min/max
  
  unsigned nxmin = (unsigned)((xprof[0] - uvrGlobalMin)/dxprof);
  unsigned nxmax = (unsigned)((uvrGlobalMax - xprof[nprof-1])/dxprof);
  
  unsigned nNewProf = nprof + nxmin + nxmax;
  
  const int dims[2] = {1, nNewProf};

  double* xNewProf = MexHandler::createDoubleArray(&plhs[0], 2, dims);
  double* yNewProf = MexHandler::createDoubleArray(&plhs[1], 2, dims);

  // And put the profile in the appropriate bins of the new array
  
  for(unsigned i=0; i < nxmin; i++) {
    xNewProf[i] = xprof[0] - (nxmin - i) * dxprof;
    yNewProf[i] = 0.0;
  }

  for(unsigned i=nxmin; i < nprof; i++) {
    xNewProf[i] = xprof[i-nxmin];
    yNewProf[i] = yprof[i-nxmin];
  }

  for(unsigned i=nprof; i < nNewProf; i++) {
    xNewProf[i] = xprof[nprof-1] + (i-nprof+1) * dxprof;
    yNewProf[i] = 0.0;
  }

  return;
}
