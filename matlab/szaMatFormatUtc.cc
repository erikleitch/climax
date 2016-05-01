/**.......................................................................
 * MATLAB Mex file for formatting MJD dates as UTC strings
 *
 * Use like:
 *
 * d=gcpMatFormatUtc(date)
 */
#include <iostream.h>
#include <math.h>

#include "gcp/util/Debug.h"
#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/RegDate.h"
#include "gcp/util/TimeVal.h"

#include "gcp/matlab/MexParser.h"

#include "mex.h"
#include "matrix.h"

using namespace gcp::util;
using namespace gcp::matlab;

void outputUtc(mxArray* plhs[], const mxArray* prhs[]);

void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{

  try {  

    if(nrhs != 1) {
      std::cerr << "Wrong number of arguments. " << std::endl
		<< "Should be: mjd[1]" << std::endl;
      return;
    }
    outputUtc(plhs, prhs);
  } catch(gcp::util::Exception& err) {
    mexErrMsgTxt(err.what());
  } catch(...) {
    mexErrMsgTxt("Caught an unknown error");
  }
}

void outputUtc(mxArray* plhs[], const mxArray* prhs[])
{
  MexParser mp(prhs[0]);

  if(mp.getNumberOfDimensions() != 2)
    ThrowError("Date must be of size 1x1");

  double* dptr = 0;
  dptr = mp.getDoubleData();

  TimeVal timeVal;

  timeVal.setMjd(*dptr);
  RegDate regDate(timeVal);
  std::ostringstream os;

  os << regDate;

  plhs[0] = mxCreateString(os.str().c_str());

  return;
}

