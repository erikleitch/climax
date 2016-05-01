/**.......................................................................
 * MATLAB Mex file for calculating UVW, given an array of antenna
 * locations, a source declination and a list of HA
 *
 * Use like:
 *
 * [map] = gcpMatMap(u, v, npix, size)
 *
 */
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"
#include "gcp/array/code/share/slalib/slalib.h"
#include "gcp/util/Coordinates.h"

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
  
  //------------------------------------------------------------
  // Deal with source info
  //------------------------------------------------------------

  MexParser   uParser(prhs[0]);
  MexParser   vParser(prhs[1]);
  MexParser   npixParser(prhs[2]);
  MexParser   sizeParser(prhs[3]);

  double* uPtr  = uParser.getDoubleData();
  double* vPtr  = vParser.getDoubleData();
  unsigned npix = (unsigned)*(npixParser.getDoubleData());
  double size   = *(sizeParser.getDoubleData());

  // First dimension of eastPtr should be the number of iterations

  unsigned nuv = uParser.getNumberOfElements();

  //------------------------------------------------------------
  // Output structure
  //------------------------------------------------------------

  std::vector<int> dims(2);
  dims[0] = 1;
  dims[1] = 1;

  plhs[0] = mxCreateStructArray(2, &dims[0], 0, NULL);

  dims[0] = npix;
  dims[1] = npix;

  double*      xOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],       "x", dims);
  double*      yOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],       "y", dims);
  double*    mapOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],     "map", dims);

  //------------------------------------------------------------
  // Loop over iterations
  //------------------------------------------------------------

  unsigned posInd, histInd, chisqInd;
  double chisqMin;

  size = size/180*M_PI;

  double xmin,xmax,ymin,ymax;

  xmin = -size/2;
  xmax = +size/2;
  ymin = -size/2;
  ymax = +size/2;

  double dx,dy;

  dx = size/npix;
  dy = size/npix;

  double x,y,arg;
  unsigned ind;

  COUT("Here 0");
  for(unsigned ix = 0; ix < npix; ix++) {
    x = xmin + dx*ix;
    for(unsigned iy = 0; iy < npix; iy++) {
      ind = ix * npix + iy;
      mapOutPtr[ind] = 0.0;

      xOutPtr[ind] = x * 180/M_PI;
      yOutPtr[ind] = y * 180/M_PI;
      
      y = ymin + dy*iy;
      for(unsigned iuv=0; iuv < nuv; iuv++) {
	arg = 2*M_PI*(uPtr[iuv] * x + vPtr[iuv] * y);
	mapOutPtr[ind] += cos(arg);
      }
    }
  }

  COUT("Here 1");
  return;
}
