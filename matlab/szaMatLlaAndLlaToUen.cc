/**.......................................................................
 * MATLAB Mex file for converting from J2000 Ra/Dec to Galactic coordinates
 *
 * Use like:
 *
 * d=gcpMatEqGal(ra, dec)
 *
 */
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

#include "gcp/util/Coordinates.h"
#include "gcp/util/DecAngle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Vector.h"

#include "gcp/array/code/share/slalib/slalib.h"

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
  

  Angle lngFid, latFid, lng, lat;
  double altFid, alt;

  latFid.setDegrees(MexParser::getString(prhs[0]));
  lngFid.setDegrees(MexParser::getString(prhs[1]));
  altFid = *MexParser::getDoubleData(prhs[2]);

  lat.setDegrees(MexParser::getString(prhs[3]));
  lng.setDegrees(MexParser::getString(prhs[4]));
  alt = *MexParser::getDoubleData(prhs[5]);

  Vector<double> uen = Coordinates::llaAndLlaToUen(lngFid, latFid, altFid, lng, lat, alt);

  // Create ouput lat/longitude arrays, in radians
  
  double* uenPtr = MexHandler::createDoubleArray(&plhs[0], 3);

  uenPtr[0] = uen[0];
  uenPtr[1] = uen[1];
  uenPtr[2] = uen[2];

  return;
}
