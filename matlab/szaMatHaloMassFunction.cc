/**.......................................................................
 * MATLAB Mex file for calculating halo mass functions
 *
 * Use like:
 *
 * d=gcpMatHaloMassFunction(mass(solar), sigma8, redshift)
 *
 */
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

#include "gcp/util/HaloMassFunction.h"

#include "gcp/array/code/share/slalib/slalib.h"
#include "gcp/sfd/subs_predict.h"

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
  try {

    if(nrhs != 3) {
      std::cerr << "Usage: d = gcpMatHaloMassFunction(mass(solar), sigma8(e.g., 0.7), redshift(e.g., 1.0))" << std::endl << std::endl;
      return;
    }

    gcp::util::Logger::installStdoutPrintFn(MexHandler::stdoutPrintFn);
    gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);
    gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
    gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
    gcp::util::ErrHandler::installLogFn(MexHandler::logFn);
  
    // Check input/output arguments
  
    if(!MexParser::dimensionsMatch(prhs[0], prhs[1])) {
      ThrowError("Dimensions of l and b arrays don't match");
    }
    
    if(!MexParser::dimensionsMatch(prhs[1], prhs[2])) {
      ThrowError("Dimensions of l/b and freq arrays don't match");
    }
    
    // Get l/b, assumed to be in radians
    
    double* mass     = MexParser::getDoubleData(prhs[0]);
    double* sigma8   = MexParser::getDoubleData(prhs[1]);
    double* redshift = MexParser::getDoubleData(prhs[2]);
    
    unsigned n = MexParser::getNumberOfElements(prhs[0]);
    
    double* dndm = MexHandler::createDoubleArray(&plhs[0], 
						 MexParser::getNumberOfDimensions(prhs[0]), 
						 MexParser::getDimensions(prhs[0]));
    
    for(unsigned i=0; i < n; i++) {
      dndm[i] = gcp::util::HaloMassFunction::haloMassFunction(HaloMassFunction::TYPE_TINKER, mass[i], sigma8[i], redshift[i]);
    }
    
  } catch(gcp::util::Exception& err) {
    mexErrMsgTxt(err.what());
  } catch(...) {
    mexErrMsgTxt("Caught an unknown error");
  }
}
