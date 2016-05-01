#include "gcp/matlab/MlPtsrc.h"

#include <cmath>

using namespace std;

using namespace gcp::matlab;

/**.......................................................................
 * Constructor.
 */
MlPtsrc::MlPtsrc() {}

/**.......................................................................
 * Destructor.
 */
MlPtsrc::~MlPtsrc() {}

void MlPtsrc::kernelparams(double& peak, double& width, double sobs, double gamma, double sigma)
{
  peak = 0.5*(sobs + sqrt(sobs*sobs - 4*gamma*sigma*sigma));
  double d2psids2 = gamma/(peak*peak) - 1.0/(sigma*sigma);
  width = 1.0/sqrt(-d2psids2);
}

void MlPtsrc::kernelint(unsigned nsrc, double* sobs, double* retVals, double gamma, double sigma, double nsig, unsigned nbin)
{
  // Iterate over sources, computing the integral for each one
 
  double pkpsi, sigpsi;

  for(unsigned isrc=0; isrc < nsrc; isrc++) {

    // Determine kernel parameters for this source

    kernelparams(pkpsi, sigpsi, sobs[isrc], gamma, sigma);

    double smin = pkpsi - nsig * sigpsi;
    double smax = pkpsi + nsig * sigpsi;

    if(smin <= 0)
      smin = sigpsi/10;

    double ds = (smax - smin)/(nbin-1);

    retVals[isrc] = 0.0;
    for(unsigned ibin=0; ibin < nbin; ibin++) {
      double s = smin + ds*ibin;
      retVals[isrc] += pow(s, -gamma) * exp(-(s - sobs[isrc])*(s - sobs[isrc])/(2*sigma*sigma)) * ds;
    }
  }
}

void MlPtsrc::dkernelint(unsigned nsrc, double* sobs, double* retVals, double gamma, double sigma, double nsig, unsigned nbin)
{
  // Iterate over sources, computing the integral for each one
 
  double pkpsi, sigpsi;

  for(unsigned isrc=0; isrc < nsrc; isrc++) {

    // Determine kernel parameters for this source

    kernelparams(pkpsi, sigpsi, sobs[isrc], gamma, sigma);

    double smin = pkpsi - nsig * sigpsi;
    double smax = pkpsi + nsig * sigpsi;

    if(smin <= 0)
      smin = sigpsi/10;

    double ds = (smax - smin)/(nbin-1);

    retVals[isrc] = 0.0;
    for(unsigned ibin=0; ibin < nbin; ibin++) {
      double s = smin + ds*ibin;
      retVals[isrc] += -log(s) * pow(s, -gamma) * exp(-(s - sobs[isrc])*(s - sobs[isrc])/(2*sigma*sigma)) * ds;
    }
  }
}


void MlPtsrc::d2kernelint(unsigned nsrc, double* sobs, double* retVals, double gamma, double sigma, double nsig, unsigned nbin)
{
  // Iterate over sources, computing the integral for each one
 
  double pkpsi, sigpsi;

  for(unsigned isrc=0; isrc < nsrc; isrc++) {

    // Determine kernel parameters for this source

    kernelparams(pkpsi, sigpsi, sobs[isrc], gamma, sigma);

    double smin = pkpsi - nsig * sigpsi;
    double smax = pkpsi + nsig * sigpsi;

    if(smin <= 0)
      smin = sigpsi/10;

    double ds = (smax - smin)/(nbin-1);

    retVals[isrc] = 0.0;
    for(unsigned ibin=0; ibin < nbin; ibin++) {
      double s = smin + ds*ibin;
      double ls = log(s);
      retVals[isrc] += -ls*ls * pow(s, -gamma) * exp(-(s - sobs[isrc])*(s - sobs[isrc])/(2*sigma*sigma)) * ds;
    }
  }
}

