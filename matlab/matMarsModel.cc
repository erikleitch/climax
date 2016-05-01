/**.......................................................................
 * MATLAB Mex file for reading and interpolating the Caltech Mars Model
 *
 */
#include <iostream.h>
#include <math.h>

#include "gcp/util/Debug.h"
#include "gcp/util/ModelReader.h"

#include "mex.h"
#include "matrix.h"

using namespace gcp::util;
using namespace gcp::matlab;

void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  try {  

    if(nrhs != 3) {
      std::cerr << "Wrong number of arguments. " << std::endl
		<< "Should be: fileName mjd[N] freq[N]" << std::endl;
      return;
    }

    returnPlanckTemps(plhs, prhs);

  } catch(gcp::util::Exception& err) {
    mexErrMsgTxt(err.what());
  } catch(...) {
    mexErrMsgTxt("Caught an unknown error");
  }
}

void returnPlanckTemps(mxArray* plhs[], const mxArray* prhs[])
{
  char* fileName;
  ModelReader mr(".", mxArrayToString(prhs[0]));

  MexParser mp;

  mp.setTo(prhs[1]);
  double mjds = mp.getDoubleData();

  mp.setTo(prhs[2]);
  double freqs = mp.getDoubleData();


}

