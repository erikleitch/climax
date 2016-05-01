/**.......................................................................
 * MATLAB Mex file for converting from Apparent place to J2000
 *
 */
#include <iostream>
#include <math.h>

#include "gcp/util/HourAngle.h"
#include "gcp/util/Astrometry.h"
#include "gcp/util/Date.h"
#include "gcp/util/Debug.h"
#include "gcp/util/Declination.h"
#include "gcp/util/TimeVal.h"

#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

#include "mex.h"
#include "matrix.h"

using namespace gcp::util;
using namespace gcp::matlab;

void appToJ2000(const mxArray* prhs[], mxArray** ra, mxArray** dec);

void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  try {  

    gcp::util::Logger::installStdoutPrintFn(MexHandler::stdoutPrintFn);
    gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);

    gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
    gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
    gcp::util::ErrHandler::installLogFn(MexHandler::logFn);

    if(nrhs != 3) {
      std::cerr << "Wrong number of arguments. " << std::endl
		<< "Should be: obsRa(hours) obsDec(degrees) obsDate(mjd)" 
		<< std::endl;
      return;
    }

    appToJ2000(prhs, &plhs[0], &plhs[1]);

  } catch(gcp::util::Exception& err) {
    mexErrMsgTxt(err.what());
  } catch(...) {
    mexErrMsgTxt("Caught an unknown error");
  }
}

void appToJ2000(const mxArray* prhs[], mxArray** ra, mxArray** dec)
{
  HourAngle meanRa, obsRa;
  Declination meanDec, obsDec;
  TimeVal obsMjd;

  MexParser mp;

  mp.setTo(prhs[0]);
  if(mp.isString()) {
    obsRa.setHours(mp.getString());
  } else {
    double* raPtr = mp.getDoubleData();
    obsRa.setHours(*raPtr);
  }

  mp.setTo(prhs[1]);
  if(mp.isString()) {
    obsDec.setDegrees(mp.getString());
  } else {
    double* decPtr = mp.getDoubleData();
    obsDec.setDegrees(*decPtr);
  }

  mp.setTo(prhs[2]);
  if(mp.isString()) {
    Date date;
    date.setToDateAndTime(mp.getString());
    obsMjd.setMjd(date.getMjd());
    COUT("Date = " << date << " MJD = " << obsMjd.getMjd());
  } else {
    double* mjdPtr = mp.getDoubleData();
    obsMjd.setMjd(*mjdPtr);
  }

  // Need to convert from epoch of observation to standard epoch

  std::cout << "apparent Ra is: " << obsRa << std::endl;
  std::cout << "apparent Dec is: " << obsDec << std::endl;
  std::cout << "obsMjd is: " << obsMjd << std::endl;

  gcp::util::Astrometry::apparentToJ2000Place(obsRa, obsDec, obsMjd, 
					      meanRa, meanDec);

  std::cout << "J2000 meanRa is:  " << meanRa << std::endl;
  std::cout << "J2000 meanDec is: " << meanDec << std::endl;

  double* raPtr;
  double* decPtr;

  *ra  = MexHandler::createMatlabArray(1, DataType::DOUBLE);
  raPtr = (double*)mxGetData(*ra);
  *raPtr = meanRa.hours();

  *dec = MexHandler::createMatlabArray(1, DataType::DOUBLE);
  decPtr = (double*)mxGetData(*dec);
  *decPtr = meanDec.degrees();
} 

