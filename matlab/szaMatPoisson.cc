/**.......................................................................
 * MATLAB Mex file for calculating UVW, given an array of antenna
 * locations, a source declination and a list of HA
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
  
  MexParser nMeanParser(prhs[0]);
  MexParser nSampParser(prhs[1]);

  unsigned nMean = (unsigned)*(nMeanParser.getDoubleData());
  unsigned nSamp = nSampParser.getNumberOfElements();
  
  double* samples = MexHandler::createDoubleArray(&plhs[0], 
						  nSampParser.getNumberOfDimensions(), 
						  nSampParser.getDimensions()); 

  double seed = -1;

  if(nrhs > 2) {
    MexParser seedParser(prhs[2]);
    seed = *(seedParser.getDoubleData());
  }

  COUT("Seed is: " << seed);

  if(seed > 0) {
    Sampler::seed((unsigned)seed);
  }


  if(nMean < 100) {
    
    // For small n, generate poisson samples
    
    std::vector<unsigned> n = 
      Sampler::generatePoissonSamples(nMean, nSamp);
    
    for(unsigned i=0; i < nSamp; i++) {
      samples[i] = (double)n[i];
    }

  } else {
    
    // In the limit of large N, poisson(k) --> gauss(k,sqrt(k))
    
    std::vector<double> n = 
      Sampler::generateGaussianSamples(sqrt((double)nMean), nSamp);

    for(unsigned i=0; i < nSamp; i++) {
      samples[i] = nearbyint((double)(n[i] + nMean));
    }

  }

  return;
}
