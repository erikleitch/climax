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

  std::ostringstream usage;

  usage << "Usage: " << std::endl;
  usage << std::endl;
  usage << "   [rt prt] = gcpMatProbN(p, r, rsigma)" << std::endl;
  usage << std::endl;
  usage << " where p should have fields:" << std::endl;
  usage << std::endl;
  usage << "    gamma -- power law spectral index" << std::endl;
  usage << "    kmin  -- the minimum number of sources to consider" << std::endl;
  usage << "    kmax  -- the maximum number of sources to consider" << std::endl;
  usage << std::endl;

  if(nrhs < 3) {
    ThrowError("Wrong number of arguments: " << usage);
  }

  // Push the top-level mxArray into a parser

  MexParser structParser(prhs[0]);
  MexParser rParser(prhs[1]);
  MexParser rsigmaParser(prhs[2]);

  double rObs      = *rParser.getDoubleData();
  double rsigmaObs = *rsigmaParser.getDoubleData();

  // Extract point source population parameters

  if(!structParser.fieldExists("gamma")) {
    ThrowError("No field gamma found: " << usage.str());
  }

  double gamma   = *(structParser.getFieldAsDouble("gamma"));

  if(!structParser.fieldExists("rmin")) {
    ThrowError("No field rmin found: " << usage.str());
  }

  double rmin    = *(structParser.getFieldAsDouble("rmin"));

  if(!structParser.fieldExists("rmax")) {
    ThrowError("No field rmax found: " << usage.str());
  }

  double rmax    = *(structParser.getFieldAsDouble("rmax"));

  if(!structParser.fieldExists("nr")) {
    ThrowError("No field nr found: " << usage.str());
  }

  unsigned nr    = (unsigned)*(structParser.getFieldAsDouble("nr"));

  double dr = (rmax - rmin)/(nr-1);

  std::vector<double> rs(nr);
  for(unsigned i=0; i < nr; i++) {
    rs[i] = rmin + i*dr;
  }

  // We assume that k is specified in #/sqarcmin, and
  // that the power law is specified in units of mJy

  if(!structParser.fieldExists("kmin")) {
    ThrowError("No field kmin found: " << usage.str());
  }

  unsigned kmin    = (unsigned)*(structParser.getFieldAsDouble("kmin"));

  if(!structParser.fieldExists("kmax")) {
    ThrowError("No field kmax found: " << usage.str());
  }

  unsigned kmax    = (unsigned)*(structParser.getFieldAsDouble("kmax"));

  double* rt   = MexHandler::createDoubleArray(&plhs[0], nr);
  double* prtp = MexHandler::createDoubleArray(&plhs[1], nr);
  double* prtm = MexHandler::createDoubleArray(&plhs[2], nr);

  //------------------------------------------------------------
  // Compute the expected mean for each value of rs
  //
  // In terms of the confusion flux S_c = (k\Omega)^{1\over gamma}
  //
  // the mean is given by dN/dS = k\Omega S^-\gamma = r^-\gamma,
  //.computed as ls below.
  //
  // In terms of discrete sources, this mean should be rounded to the
  // nearest integer at each flux level (lsr below).
  //
  // At each flux level r, we consider a fluctuation of k sources
  // about the (rounded) mean, or total number of sources = lsr + k
  // 
  //------------------------------------------------------------

  std::vector<double> ls(nr);
  std::vector<double> lsr(nr);

  for(unsigned ir=0; ir < nr; ir++) {
    ls[ir]  = pow(rs[ir], -gamma);
    lsr[ir] = round(ls[ir]);
    rt[ir]  = rs[ir];
    prtp[ir] = 0.0;
    prtm[ir] = 0.0;
  }

  double pkm, pkp, gkm, gkp, argp, argm;

  // Iterate over all k
  
  for(unsigned k=kmin; k <= kmax; k++) {
    
    // For each value of k, we have to integrate over the
    // probability of k sources of all fluxes (rs)
    
    for(unsigned ir=0; ir < nr; ir++) {
      
      // Calculate the Poisson probability of a fluctuation of k
      // sources over the mean of lsr, given a true mean ls
      
      pkp = Sampler::poissPdf((unsigned)(lsr[ir]+k), ls[ir]);

      if(k > lsr[ir]) {
	pkm = 0.0;
      } else {
	pkm = Sampler::poissPdf((unsigned)(lsr[ir]-k), ls[ir]);
      }
      
      // Calculate the probability of a noise fluctuation of rObs -
      // k*rs
      
      argp = (rObs - k*rs[ir])/rsigmaObs;
      argm = (rObs + k*rs[ir])/rsigmaObs;
      gkp = exp(-(argp*argp)/2);
      gkm = exp(-(argm*argm)/2);
      
      prtp[ir] += pkp * gkp;
      prtm[ir] += pkm * gkm;

    }

  }

  return;
}
