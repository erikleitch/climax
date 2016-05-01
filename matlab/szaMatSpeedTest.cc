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
  
  double* xyz     = MexParser::getDoubleData(prhs[0]);
  double* cosHa   = MexParser::getDoubleData(prhs[1]);
  double* sinHa   = MexParser::getDoubleData(prhs[2]);
  double* dec     = MexParser::getDoubleData(prhs[3]);

  unsigned nHa    = MexParser::getNumberOfElements(prhs[1]);
  const int* dims = MexParser::getDimensions(prhs[0]);
  unsigned nAnt   = dims[1];

  unsigned nBase = (nAnt * (nAnt-1))/2;

  const int uvDims[2] = {nBase, nHa};

  double* u     = MexHandler::createDoubleArray(&plhs[0], 2, uvDims);
  double* v     = MexHandler::createDoubleArray(&plhs[1], 2, uvDims);

  double x,y,z;

  double dx, dy, dz;
  unsigned ind1, ind2;
  double cd = cos((*dec)/180 * M_PI); 
  double sd = sin((*dec)/180 * M_PI); 

  unsigned uvInd;
  unsigned iBase = 0;

  for(unsigned iAnt1=0; iAnt1 < nAnt-1; iAnt1++)
    for(unsigned iAnt2=iAnt1+1; iAnt2 < nAnt; iAnt2++) {

      ind1 = iAnt1*3;
      ind2 = iAnt2*3;

      dx =   xyz[ind1] -   xyz[ind2];
      dy = xyz[ind1+1] - xyz[ind2+1];
      dz = xyz[ind1+2] - xyz[ind2+2];

      for(unsigned iHa=0; iHa < nHa; iHa++) {
	uvInd = iHa * nBase + iBase;
	u[uvInd] =       sinHa[iHa] * dx +      cosHa[iHa] * dy;
	v[uvInd] = -sd * cosHa[iHa] * dx + sd * sinHa[iHa] * dy + cd * dz;
      }

      ++iBase;
    }


  return;
}
