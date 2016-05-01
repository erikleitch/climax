/**.......................................................................
 * MATLAB Mex file for accessing the uranus model
 *
 */
#include <iostream.h>
#include <math.h>

#include "gcp/util/Debug.h"
#include "gcp/matlab/UranusModel.h"
#include "gcp/matlab/MexHandler.h"

#include "mex.h"
#include "matrix.h"

using namespace gcp::util;
using namespace gcp::matlab;

void calcUranusFlux(const mxArray* prhs[], mxArray** flux, mxArray** error);

void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  gcp::util::Logger::installStdoutPrintFn(MexHandler::stdoutPrintFn);
  gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);
  
  gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
  gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
  gcp::util::ErrHandler::installLogFn(MexHandler::logFn);

  if(nrhs != 5) {
    std::cerr << "Wrong number of arguments. " << std::endl
	      << "Should be: directory fileName temp|flux|size mjd[N], frequency[N]" << std::endl;
    return;
  }
  
  int dims = 10;
  
  calcUranusFlux(prhs, &plhs[0], &plhs[1]);
}

void calcUranusFlux(const mxArray* prhs[], mxArray** flux, mxArray** error)
{
  UranusModel uranus;

  uranus.setDirectory(prhs[0]);
  uranus.setFileName(prhs[1]);
  uranus.setType(prhs[2]);
  uranus.setMjdArray(prhs[3]);
  uranus.setFrequencyArray(prhs[4]);

  uranus.getData(flux, error);
} 

