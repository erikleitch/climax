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

  std::ostringstream usage;

  usage << "Usage: " << std::endl
	<< std::endl
	<< "       [f x y] = gcpMatPtSrcGen(k, gamma, fluxMin, xDeg, yDeg, doRand, seed, fluxMax)" << std::endl
	<< std::endl
	<< "Where: " << std::endl
	<< std::endl
	<< " k       = normalization, in #/(sqarcmin * mJy)" << std::endl
	<< " gamma   = power-law spectral index (dN/dS ~ S^-gamma)" << std::endl
	<< " fluxMin = min flux to generate, in mJy" << std::endl 
	<< " xDeg    = x-dimension of the field, in degrees" << std::endl 
	<< " yDeg    = y-dimension of the field, in degrees" << std::endl 
	<< " doRand  = if true  (1), generate a random number of sources, drawn from a Poisson distribution of the correct mean" << std::endl
	<< "           if false (0), generate exactly the mean number of sources" << std::endl
	<< " seed    = if < 0, do not reseed the random number generator, else use 'seed' as the starting point" << std::endl
	<< " fluxMax = if < 0, generate sources to arbitrary flux" << std::endl 
	<< "           if > 0, generate sources to fluxMax, in mJy" << std::endl << std::endl;

  if(nrhs < 8) {
    ThrowError(usage.str());
  }

  MexParser kParser(prhs[0]);
  MexParser gParser(prhs[1]);
  MexParser fminParser(prhs[2]);
  MexParser axParser(prhs[3]);
  MexParser ayParser(prhs[4]);
  MexParser doRandParser(prhs[5]);
  MexParser sParser(prhs[6]);
  MexParser fmaxParser(prhs[7]);

  double* k       = kParser.getDoubleData();
  double* gamma   = gParser.getDoubleData();
  double* seed    = sParser.getDoubleData();

  bool doRand = *(doRandParser.getDoubleData()) > 0;

  //  COUT("Seed is: " << *seed);

  if(*seed > 0) {
    Sampler::seed((unsigned)*seed);
  }

  // We assume that k is specified in #/sqarcmin, and
  // that the power law is specified in units of mJy

  Flux fu(Flux::MilliJansky(), 1.0);
  SolidAngle au(SolidAngle::SqArcMinutes(), 1.0);

  PtSrcGen gen;
  gen.setDnDs(*k, *gamma, fu, au);

  // Get the minimum flux, also in mJy

  Flux fluxMin;
  fluxMin.setMilliJy(*fminParser.getDoubleData());

  Flux fluxMax;
  fluxMax.setMilliJy(*fmaxParser.getDoubleData());

  bool havefMax = false;
  if(fluxMax.Jy() >= 0.0) {
    havefMax = true;
  }

  // Get the dimensions of each side of the field, in degrees

  Angle x, y;
  x.setDegrees(*axParser.getDoubleData());
  y.setDegrees(*ayParser.getDoubleData());

  SolidAngle area;
  area.setSr(x.radians() * y.radians());

  std::vector<Flux>  srcFlux;
  std::vector<Angle> srcX;
  std::vector<Angle> srcY;

  if(havefMax)
    gen.generateSources(fluxMin, fluxMax, x, y, srcFlux, srcX, srcY, doRand);
  else
    gen.generateSources(fluxMin, x, y, srcFlux, srcX, srcY, doRand);

  double* sfA = MexHandler::createDoubleArray(&plhs[0], srcFlux.size());
  double* sxA = MexHandler::createDoubleArray(&plhs[1], srcFlux.size());
  double* syA = MexHandler::createDoubleArray(&plhs[2], srcFlux.size());

  for(unsigned i=0; i < srcFlux.size(); i++) {
    sfA[i] = srcFlux[i].mJy();
    sxA[i] = srcX[i].degrees();
    syA[i] = srcY[i].degrees();
  }

  return;
}
