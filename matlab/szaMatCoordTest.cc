#include "gcp/util/Debug.h"

#include "gcp/util/Coordinates.h"
#include "gcp/util/Declination.h"
#include "gcp/util/HourAngle.h"
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

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
  

  HourAngle ra0;
  ra0.setHours(*(MexParser::getDoubleData(prhs[0])));

  Declination dec0;
  dec0.setDegrees(*(MexParser::getDoubleData(prhs[1])));

  double* thetaPtr = MexParser::getDoubleData(prhs[2]);
  double* rhoPtr   = MexParser::getDoubleData(prhs[3]);

  Angle theta;
  theta.setDegrees(*thetaPtr);

  Angle rho;
  rho.setDegrees(*rhoPtr);
  
  HourAngle ra;
  Declination dec;

  COUT("ra0   = " << ra0);
  COUT("dec0  = " << dec0);
  COUT("theta = " << theta.degrees());
  COUT("rho   = " << rho.degrees());

  Coordinates::raDecAndThetaRhoToRaDec(ra0, dec0, theta, rho, ra, dec);

  int mdims[1] = {2};
  plhs[0] = MexHandler::createMatlabArray(1, mdims, DataType::DOUBLE);
  plhs[1] = MexHandler::createMatlabArray(1, mdims, DataType::DOUBLE);
  double* raPtr = (double*)mxGetData(plhs[0]);
  double* decPtr = (double*)mxGetData(plhs[1]);

  *raPtr  = ra.hours();
  *decPtr = dec.degrees();

  return;
}

