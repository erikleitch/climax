/**.......................................................................
 * MATLAB Mex file for calculating SFD dust model
  *
 * Use like:
 *
 * d=gcpMatEqGal(galLat, galLng, nu)
 *
 */
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

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

    if(nrhs != 4) {
      std::cerr << "Usage: d = gcpMatSfd(l(rad), b(rad), freq(GHz), units(string))" << std::endl << std::endl
		<< " Where l = Galactic longitude" << std::endl
		<< "       b = Galactic latitude"  << std::endl << std::endl
		<< "and where the last argument should be one of: \"MJy\", \"microK\", or \"thermo\"," << std::endl 
		<< "to get values in MJy/sr, brightness temperature micro-Kelvin, " << std::endl 
		<< "or thermodynamic micro-Kelvin, respectively." << std::endl;
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
    
    if(!MexParser::isString(prhs[3])) {
      ThrowError("last argument should be one of: \"MJy\", \"microK\", or \"thermo\", " << std::endl << "to get values in MJy/sr, brightness temperature micro-Kelvin, " << std::endl << "or thermodynamic micro-Kelvin, respectively");
    }
    
    // Get l/b, assumed to be in radians
    
    double* l    = MexParser::getDoubleData(prhs[0]);
    double* b    = MexParser::getDoubleData(prhs[1]);
    double* freq = MexParser::getDoubleData(prhs[2]);
    std::string units = MexParser::getString(prhs[3]);
    
    // Create output lat/longitude arrays, in radians
    
    double* Inu = MexHandler::createDoubleArray(&plhs[0], 
						MexParser::getNumberOfDimensions(prhs[0]), 
						MexParser::getDimensions(prhs[0]));
    
    unsigned nGal = MexParser::getNumberOfElements(prhs[0]);
    
    predict_thermal(Inu, nGal, l, b, freq, 
		    "/data/gcpdaq/dust/maps/", "I4096", (char*)units.c_str(), 8, 1, 1, 0);
    
  } catch(gcp::util::Exception& err) {
    mexErrMsgTxt(err.what());
  } catch(...) {
    mexErrMsgTxt("Caught an unknown error");
  }
}
