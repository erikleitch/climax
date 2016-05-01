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

  // Check arguments

  MexHandler::checkArgs(nlhs, nrhs, 3, 5);

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

  // Get reference location

  double* enuRef  = MexParser::getDoubleData(prhs[2]);
  double eastRef  = enuRef[0];
  double northRef = enuRef[1];
  double upRef    = enuRef[2];

  // Get dec

  DecAngle dec;

  if(MexParser::isString(prhs[3]))
    dec.setDegrees(MexParser::getString(prhs[3]));
  else
    dec.setRadians(*MexParser::getDoubleData(prhs[3]));
 
  // Get hour angles
  
  double* has  = MexParser::getDoubleData(prhs[4]);
  unsigned nha = MexParser::getNumberOfElements(prhs[4]);

  // Create output lat/longitude arrays, in radians
  
  double* u = MexHandler::createDoubleArray(&plhs[0], nha);
  double* v = MexHandler::createDoubleArray(&plhs[1], nha);
  double* w = MexHandler::createDoubleArray(&plhs[2], nha);

  // Altitude doesn't matter, since we're going to take relative
  // positions in the end anyway, so just pass 0.0

  Vector<double> xyzRef = Coordinates::laAndUenToXyz(lat, 0.0, upRef, eastRef, northRef);
  Vector<double> xyz    = Coordinates::laAndUenToXyz(lat, 0.0, up, east, north);
  
  HourAngle ha;
  Vector<double> uvw;

  for(unsigned i=0; i < nha; i++) {
    ha.setHours(has[i]);

    uvw = Coordinates::haDecAndXyzToUvw(ha, dec, xyz[0]-xyzRef[0], xyz[1]-xyzRef[1], xyz[2]-xyzRef[2]);

    u[i] = uvw[0];
    v[i] = uvw[1];
    w[i] = uvw[2];
  }

  return;
}
