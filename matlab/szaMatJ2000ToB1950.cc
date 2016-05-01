/**.......................................................................
 * MATLAB Mex file for converting from J2000 coordinates to B1950
 * coordinates
 *
 */
#include <iostream.h>
#include <math.h>

#include "gcp/util/HourAngle.h"
#include "gcp/util/Astrometry.h"
#include "gcp/util/Debug.h"
#include "gcp/util/DecAngle.h"
#include "gcp/util/TimeVal.h"

#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

#include "mex.h"
#include "matrix.h"

using namespace gcp::util;
using namespace gcp::matlab;

void j2000ToB1950(const mxArray* prhs[], mxArray** ra, mxArray** dec);

void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  try {  

    gcp::util::Logger::installStdoutPrintFn(MexHandler::stdoutPrintFn);
    gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);

    gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
    gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
    gcp::util::ErrHandler::installLogFn(MexHandler::logFn);

    if(nrhs != 2) {
      std::cerr << "Wrong number of arguments. " << std::endl
		<< "Should be: raJ2000(hours) decJ2000(degrees)" 
		<< std::endl;
      return;
    }

    j2000ToB1950(prhs, &plhs[0], &plhs[1]);

  } catch(gcp::util::Exception& err) {
    mexErrMsgTxt(err.what());
  } catch(...) {
    mexErrMsgTxt("Caught an unknown error");
  }
}

void j2000ToB1950(const mxArray* prhs[], mxArray** ra, mxArray** dec)
{
  HourAngle raB1950, raJ2000;
  DecAngle decB1950, decJ2000;

  MexParser mp;

  mp.setTo(prhs[0]);
  if(mp.isString()) {
    raJ2000.setHours(mp.getString());
  } else {
    raJ2000.setHours(*mp.getDoubleData());
  }

  mp.setTo(prhs[1]);
  if(mp.isString()) {
    decJ2000.setDegrees(mp.getString());
  } else {
    decJ2000.setDegrees(*mp.getDoubleData());
  }

  gcp::util::Astrometry::j2000ToB1950(raJ2000, decJ2000, raB1950, decB1950);

  COUT("raB1950 = " << raB1950);
  COUT("decB1950 = " << decB1950);

  *ra  = MexHandler::createMatlabArray(1, DataType::DOUBLE);
  double* raPtr = (double*)mxGetData(*ra);
  *raPtr = raB1950.hours();

  *dec = MexHandler::createMatlabArray(1, DataType::DOUBLE);
  double* decPtr = (double*)mxGetData(*dec);
  *decPtr = decB1950.degrees();
} 

