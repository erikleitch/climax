/**.......................................................................
 * MATLAB Mex file for reading from GCP archive data files.
 *
 * Use like:
 *
 * d=gcpMatFastReadArc({'array.frame.record','corr.band0.usb[0][0]',
 *                     'antenna*.tracker.actual double', 'antenna*.tracker.source string'},
 *                     '06-jan-2005:15','06-jan-2005:16',
 *                     '/data/gcpdaq/arc','/home/gcpdaq/carma_unstable/gcp/array/conf/cal');
 *
 */
#include "gcp/util/Debug.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Date.h"
#include "gcp/util/Declination.h"
#include "gcp/util/TimeVal.h"

#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

#include "mex.h"
#include "matrix.h"

#include <iostream>
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
  
  HourAngle meanRa, obsRa;
  Declination meanDec, obsDec;
  TimeVal obsMjd;

  MexParser mp;

  mp.setTo(prhs[0]);
  if(mp.isString()) {
    meanRa.setHours(mp.getString());
    COUT("RA = " << mp.getString());
  } else {
    double* raPtr = mp.getDoubleData();
    meanRa.setHours(*raPtr);
  }

  mp.setTo(prhs[1]);
  if(mp.isString()) {
    meanDec.setDegrees(mp.getString());
  } else {
    double* decPtr = mp.getDoubleData();
    meanDec.setDegrees(*decPtr);
  }

  mp.setTo(prhs[2]);
  if(mp.isString()) {
    Date date;
    date.setToDateAndTime(mp.getString());
    obsMjd.setMjd(date.getMjd());
    COUT("Date = " << date << " MJD = " << obsMjd.getMjd());
  } else {
    COUT("Mjd is an array");
    double* mjdPtr = mp.getDoubleData();
    obsMjd.setMjd(*mjdPtr);
  }

  unsigned nMjd = mp.getNumberOfElements();

  return;
}
