/**.......................................................................
 * MATLAB Mex file for accessing the mars model
 *
 */
#include <iostream.h>
#include <math.h>

#include "gcp/util/Debug.h"
#include "gcp/matlab/MarsModelNew.h"
#include "gcp/matlab/MexHandler.h"

#include "mex.h"
#include "matrix.h"

using namespace gcp::util;
using namespace gcp::matlab;

void calcMarsFlux(const mxArray* prhs[], mxArray** flux, mxArray** error);

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
  
  calcMarsFlux(prhs, &plhs[0], &plhs[1]);
}

void calcMarsFlux(const mxArray* prhs[], mxArray** flux, mxArray** error)
{
  MarsModelNew mars;

  mars.setDirectory(prhs[0]);
  mars.setFileName(prhs[1]);
  mars.setType(prhs[2]);
  mars.setMjdArray(prhs[3]);
  mars.setFrequencyArray(prhs[4]);

  mars.getData(flux, error);
} 

