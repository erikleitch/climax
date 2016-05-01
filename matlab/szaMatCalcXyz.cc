/**.......................................................................
 * MATLAB Mex file for calculating TMS XYZ coordinates, given an array
 * of antenna locations, and a latitude
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
  
  // Get latitude

  Angle lat;

  if(MexParser::isString(prhs[0]))
    lat.setDegrees(MexParser::getString(prhs[0]));
  else
    lat.setRadians(*MexParser::getDoubleData(prhs[0]));

  // Get antenna locations

  double* enu     = MexParser::getDoubleData(prhs[1]);
  double east     = enu[0];
  double north    = enu[1];
  double up       = enu[2];

  // Create output TMS XYZ array
  
  double* xyzOut = MexHandler::createDoubleArray(&plhs[0], 3);

  // Altitude doesn't matter, since we're going to take relative
  // positions in the end anyway, so just pass 0.0

  Vector<double> xyz    = Coordinates::laAndUenToXyz(lat, 0.0, up, east, north, false);

  xyzOut[0] = xyz[0];
  xyzOut[1] = xyz[1];
  xyzOut[2] = xyz[2];

  return;
}
