/**.......................................................................
 * MATLAB Mex file for reading from UVF files
 *
 * Use like:
 *
 * d=gcpMatReadUvf({'array.frame.record','corr.band0.usb[0][0]',
 *                     'antenna*.tracker.actual double', 'antenna*.tracker.source string'},
 *                     '06-jan-2005:15','06-jan-2005:16',
 *                     '/data/gcpdaq/arc','/home/gcpdaq/carma_unstable/gcp/array/conf/cal');
 *
 */
#include "gcp/util/ChisqVariate.h"
#include "gcp/util/MiriadIo.h"
#include "gcp/util/FitsUvfReader.h"
#include "gcp/util/FitsBinTableReader.h"

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
  
  double* chisq = MexParser::getDoubleData(prhs[0]);
  double* ndof  = MexParser::getDoubleData(prhs[1]);

  double* retVal = MexHandler::createDoubleArray(plhs, 1);

  ChisqVariate chisqVar;
  chisqVar.setChisq(*chisq, (unsigned) *ndof);

  *retVal = chisqVar.pdf().value();

  return;
}
