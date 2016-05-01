/**.......................................................................
 * MATLAB Mex file for calculating shadowing, given a specified GCP
 * and CARMA array
 *
 * Use like:
 *
 * shadowed(bool[Nobs x Nant]) = gcpMatShadowCarma(az(double[Nobs x Nant]), el(double[Nobs x Nant]), 
 *                                                 gcpConf(string), carmaConf(string))
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

  std::ostringstream os;

  os << "Usage: " << std::endl
     << "   [az el] = gcpMatAzEl(double/string lat, double dec, double ha[NxM])" << std::endl
     << "Where:" << std::endl
     << " lat is the latitude, in degrees, or a sexagesimal degree string" << std::endl
     << " dec is the declination in degrees, or a sexagesimal degree strings" << std::endl
     << " ha is an arbitrarily-sized array of hour angles, in hours" << std::endl
     << std::endl;

  //------------------------------------------------------------
  // Input argument parsing
  //------------------------------------------------------------

  MexParser latParser(prhs[0]);
  MexParser decParser(prhs[1]);
  MexParser haParser(prhs[2]);

  Angle lat;
  if(latParser.isString()) {
    lat.setDegrees(latParser.getString());
  } else {
    lat.setDegrees(*latParser.getDoubleData());
  }

  DecAngle dec;
  if(decParser.isString()) {
    dec.setDegrees(decParser.getString());
  } else {
    dec.setDegrees(*decParser.getDoubleData());
  }

  double* haPtr   = haParser.getDoubleData();
  unsigned nHa = haParser.getNumberOfElements();

  //------------------------------------------------------------
  // Create output vector to have the same dimensions as the input
  //------------------------------------------------------------

  double* azPtr = MexHandler::createDoubleArray(&plhs[0], haParser.getNumberOfDimensions(), haParser.getDimensions());
  double* elPtr = MexHandler::createDoubleArray(&plhs[1], haParser.getNumberOfDimensions(), haParser.getDimensions());

  //------------------------------------------------------------  
  // Now Calculate the az and el
  //------------------------------------------------------------  
  
  Vector<Angle> angVec;
  HourAngle ha;
  for(unsigned i=0; i < nHa; i++) {
    ha.setHours(haPtr[i]);
    angVec = Coordinates::laAndHaDecToAzEl(lat, 2195.7, ha, dec, false);

    double azVal = angVec[0].degrees();

    while(azVal > 360)
      azVal -= 360;
    while(azVal < 0)
      azVal += 360;

    azPtr[i] = azVal;
    elPtr[i] = angVec[1].degrees();
  }

  return;
}
