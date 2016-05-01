/**.......................................................................
 * MATLAB Mex file for splitting two dates into a series of date
 * ranges, separated by a given time increment
 *
 * dates = gcpMatDates(start, stop);
 *
 */
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

#include "gcp/util/Date.h"
#include "gcp/util/DecAngle.h"
#include "gcp/util/HourAngle.h"

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
  
  if(nrhs < 3) {
    ThrowError("Usage: dates = gcpMatDates(start, stop, interval)"
	       << std::endl << std::endl
	       << "Where start = 'dd-mmm-yyyy:hh:mm:ss'" << std::endl
	       << "Where stop  = 'dd-mmm-yyyy:hh:mm:ss'" << std::endl
	       << "Where interval is the number of hours to break the data into" << std::endl);
  }


  MexParser startParser(prhs[0]);
  MexParser stopParser(prhs[1]);
  MexParser hourParser(prhs[2]);
  
  Date start, stop;
  
  start.setToDateAndTime(startParser.getString());
  stop.setToDateAndTime(stopParser.getString());

  double hourDiff = (stop.mjd() - start.mjd()) * 24;
  double hourIncr = *hourParser.getDoubleData();

  unsigned nIncr = (unsigned)(hourDiff/hourIncr);

  if(nIncr == 0)
    nIncr += 1;

  std::vector<int> dims(2);
  dims[0] = nIncr;
  dims[1] = 2;

  plhs[0] = mxCreateCellArray(2, &dims[0]);
  std::ostringstream os;

  Date currStart, currStop;
  unsigned i;
  for(i=0; i < nIncr-1; i++) {

    os.str("");
    currStart = start;
    currStart.addHours(i * hourIncr);

    currStop = start;
    currStop.addHours((i+1) * hourIncr);

    os << currStart.mjdToArcCal();
    mxSetCell(plhs[0], MexHandler::getMatlabIndex(nIncr, 2, i, 0), mxCreateString(os.str().c_str()));

    os.str("");
    os << currStop.mjdToArcCal();
    mxSetCell(plhs[0], MexHandler::getMatlabIndex(nIncr, 2, i, 1), mxCreateString(os.str().c_str()));
  }

  currStart = start;
  currStart.addHours(i * hourIncr);

  currStop = stop;

  os.str("");
  os << currStart.mjdToArcCal();
  mxSetCell(plhs[0], MexHandler::getMatlabIndex(nIncr, 2, i, 0), mxCreateString(os.str().c_str()));
  
  os.str("");
  os << currStop.mjdToArcCal();
  mxSetCell(plhs[0], MexHandler::getMatlabIndex(nIncr, 2, i, 1), mxCreateString(os.str().c_str()));
  
  return;
}
