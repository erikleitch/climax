/**.......................................................................
 * MATLAB Mex file for converting from B1950 coordinates to J2000
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

void b1950ToJ2000(const mxArray* prhs[], mxArray** ra, mxArray** dec);

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
		<< "Should be: raB1950(hours) decB1950(degrees)" 
		<< std::endl;
      return;
    }

    b1950ToJ2000(prhs, &plhs[0], &plhs[1]);

  } catch(gcp::util::Exception& err) {
    mexErrMsgTxt(err.what());
  } catch(...) {
    mexErrMsgTxt("Caught an unknown error");
  }
}

void b1950ToJ2000(const mxArray* prhs[], mxArray** ra, mxArray** dec)
{
  HourAngle raB1950, raJ2000;
  DecAngle decB1950, decJ2000;

  MexParser mp;

  mp.setTo(prhs[0]);

  if(mp.isString()) {
    raB1950.setHours(mp.getString());
  } else {
    raB1950.setHours(*mp.getDoubleData());
  }

  mp.setTo(prhs[1]);
  if(mp.isString()) {
    decB1950.setDegrees(mp.getString());
  } else {
    decB1950.setDegrees(*mp.getDoubleData());
  }

  gcp::util::Astrometry::b1950ToJ2000(raB1950, decB1950, raJ2000, decJ2000);

  COUT("raJ2000 = " << raJ2000);
  COUT("decJ2000 = " << decJ2000);

  *ra  = MexHandler::createMatlabArray(1, DataType::DOUBLE);
  double* raPtr = (double*)mxGetData(*ra);
  *raPtr = raJ2000.hours();

  *dec = MexHandler::createMatlabArray(1, DataType::DOUBLE);
  double* decPtr = (double*)mxGetData(*dec);
  *decPtr = decJ2000.degrees();
} 

