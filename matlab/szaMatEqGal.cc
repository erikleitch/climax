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

#include "gcp/util/HourAngle.h"
#include "gcp/util/DecAngle.h"

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
  
  // Check input/output arguments
  
  if(nrhs < 2) {
    ThrowError("Usage: [l b] = gcpMatEqGal(ra(rad), dec(rad))" 
	       << std::endl << std::endl
	       << "Where l = Galactic longitude (rad)" << std::endl
	       << "and   b = Galactic latitude (rad)" << std::endl);
  }


  // Check input/output arguments
  
  if(!MexParser::dimensionsMatch(prhs[0], prhs[1])) {
    ThrowError("Dimensions of ra and dec arrays don't match");
  }

  MexParser raParser(prhs[0]);
  MexParser decParser(prhs[1]);

  // Create ouput lat/longitude arrays, in radians
  
  double* latPtr = MexHandler::createDoubleArray(&plhs[0], 
						 MexParser::getNumberOfDimensions(prhs[0]), 
						 MexParser::getDimensions(prhs[0]));

  double* lngPtr = MexHandler::createDoubleArray(&plhs[1], 
						 MexParser::getNumberOfDimensions(prhs[0]), 
						 MexParser::getDimensions(prhs[0]));

  if(raParser.isString() && decParser.isString()) {

    HourAngle ra;
    ra.setHours(raParser.getString());
    DecAngle dec;
    dec.setDegrees(decParser.getString());

    slaEqgal(ra.radians(), dec.radians(), &latPtr[0], &lngPtr[0]);

    Angle lat,lng;

    lat.setRadians(latPtr[0]);
    lng.setRadians(lngPtr[0]);

    COUT("RA = " << ra  << " Dec = " << dec << " l  = " << lat << " b   = " << lng);

  } else {

    // Get RA/DEC, assumed to be in radians

    double* raPtr  = raParser.getDoubleData();
    double* decPtr = decParser.getDoubleData();
    
    
    for(unsigned i=0; i < MexParser::getNumberOfElements(prhs[0]); i++) {
      slaEqgal(raPtr[i], decPtr[i], &latPtr[i], &lngPtr[i]);
    }

  }
}
