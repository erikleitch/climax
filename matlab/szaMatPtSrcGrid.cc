/**.......................................................................
 * MATLAB Mex file for gridding a population
 * of sources
 *
 * Use like:
 *
 * d=gcpMatSampleTest(x(double[]), dydx(double[]), nSamp(int))
 *
 */
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"
#include "gcp/array/code/share/slalib/slalib.h"
#include "gcp/translator/DelayEngineNormal.h"
#include "gcp/array/code/share/include/rtcnetcoms.h"

#include "gcp/util/PtSrcGen.h"
#include "gcp/util/Sampler.h"

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
  
  MexParser fParser(prhs[0]);
  MexParser xParser(prhs[1]);
  MexParser yParser(prhs[2]);
  MexParser xwParser(prhs[3]);
  MexParser ywParser(prhs[4]);
  MexParser xNpixParser(prhs[5]);
  MexParser yNpixParser(prhs[6]);


  double* f    =  fParser.getDoubleData();
  double* x    =  xParser.getDoubleData();
  double* y    =  yParser.getDoubleData();
  double xw    = *(xwParser.getDoubleData());
  double yw    = *(ywParser.getDoubleData());
  unsigned xNpix = (unsigned)*(xNpixParser.getDoubleData());
  unsigned yNpix = (unsigned)*(yNpixParser.getDoubleData());

  unsigned nSrc = fParser.getNumberOfElements();

  double xPixSize = xw/xNpix;
  double yPixSize = yw/yNpix;

  std::vector<int> dims(2);
  dims[0] = xNpix;
  dims[1] = yNpix;
  double* fGridPtr = MexHandler::createDoubleArray(&plhs[0], 2, &dims[0]);
  double* nGridPtr = MexHandler::createDoubleArray(&plhs[1], 2, &dims[0]);

  unsigned nGrid = xNpix * yNpix;
  for(unsigned iGrid=0; iGrid < nGrid; iGrid++) {
    fGridPtr[iGrid] = 0.0;
    nGridPtr[iGrid] = 0.0;
  }

  unsigned ix, iy, iMat;

  for(unsigned iSrc=0; iSrc < nSrc; iSrc++) {

    ix = (unsigned)floor((x[iSrc] + xw/2)/xPixSize);
    iy = (unsigned)floor((y[iSrc] + yw/2)/yPixSize);

    iMat = MexHandler::getMatlabIndex(xNpix, yNpix, ix, iy);

    fGridPtr[iMat] += f[iSrc]; 
    nGridPtr[iMat] +=     1.0;

  }

  return;
}
