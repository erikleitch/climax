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
  usage << "   [p0 p1 p2] = gcpMatPtSrcProbCalc(p)" << std::endl;
  usage << std::endl;
  usage << " where p should have fields:" << std::endl;
  usage << std::endl;
  usage << "    k       -- n per mJy, per square arcminute" << std::endl;
  usage << "    gamma   -- power law spectral index" << std::endl;
  usage << "    xdeg    -- X Field size, degrees" << std::endl;
  usage << "    ydeg    -- Y Field size, degrees" << std::endl;
  usage << "    nx      -- Number of pixels to grid the field in, in X" << std::endl;
  usage << "    ny      -- Number of piyels to grid the field in, in Y" << std::endl;
  usage << "    rmin    -- Ratio of the minimum flux to the confusion flux" << std::endl;
  usage << "    rmax    -- Ratio of the maximum flux to the confusion flux" << std::endl;
  usage << "    rmingen -- minimum flux ratio to generate" << std::endl;
  usage << "    rmaxgen -- maximum flux ratio to generate" << std::endl;
  usage << "    rsigma  -- Ratio of the noise sigma to the confusion flux" << std::endl;
  usage << "    nr      -- The number of bins in r to consider" << std::endl;
  usage << "    niter   -- The number of iterations" << std::endl;
  usage << std::endl;

  if(nrhs < 1) {
    ThrowError("Wrong number of arguments: " << usage);
  }

  // Push the top-level mxArray into a parser

  MexParser structParser(prhs[0]);

  // Extract point source population parameters

  if(!structParser.fieldExists("k")) {
    ThrowError("No field k found: " << usage.str());
  }

  double   k     = *(structParser.getFieldAsDouble("k"));

  if(!structParser.fieldExists("gamma")) {
    ThrowError("No field gamma found: " << usage.str());
  }

  double   gamma = *(structParser.getFieldAsDouble("gamma"));

  // We assume that k is specified in #/sqarcmin, and
  // that the power law is specified in units of mJy

  Flux fu(Flux::MilliJansky(), 1.0);
  SolidAngle au(SolidAngle::SqArcMinutes(), 1.0);

  PtSrcGen gen;
  gen.setDnDs(k, gamma, fu, au);

  // Field size and resolution

  if(!structParser.fieldExists("xdeg")) {
    ThrowError("No field xdeg found: " << usage.str());
  }

  double   xdeg  = *(structParser.getFieldAsDouble("xdeg"));

  if(!structParser.fieldExists("ydeg")) {
    ThrowError("No field ydeg found: " << usage.str());
  }

  double   ydeg  = *(structParser.getFieldAsDouble("ydeg"));

  Angle x, y;
  x.setDegrees(xdeg);
  y.setDegrees(ydeg);

  if(!structParser.fieldExists("nx")) {
    ThrowError("No field nx found: " << usage.str());
  }

  unsigned nx    = (unsigned)*(structParser.getFieldAsDouble("nx"));

  if(!structParser.fieldExists("ny")) {
    ThrowError("No field ny found: " << usage.str());
  }

  unsigned ny    = (unsigned)*(structParser.getFieldAsDouble("ny"));

  double dx = xdeg / nx;
  double dy = ydeg / ny;

  // The range of r = (flux / confusion flux) we're interested in

  if(!structParser.fieldExists("rmin")) {
    ThrowError("No field rmin found: " << usage.str());
  }

  double   rmin  = *(structParser.getFieldAsDouble("rmin"));

  if(!structParser.fieldExists("rmax")) {
    ThrowError("No field rmax found: " << usage.str());
  }

  double   rmax  = *(structParser.getFieldAsDouble("rmax"));

  if(!structParser.fieldExists("nr")) {
    ThrowError("No field nr found: " << usage.str());
  }

  if(!structParser.fieldExists("rmingen")) {
    ThrowError("No field rmingen found: " << usage.str());
  }

  double   rmingen  = *(structParser.getFieldAsDouble("rmingen"));

  if(!structParser.fieldExists("rmaxgen")) {
    ThrowError("No field rmaxgen found: " << usage.str());
  }

  double   rmaxgen  = *(structParser.getFieldAsDouble("rmaxgen"));

  if(!structParser.fieldExists("nr")) {
    ThrowError("No field nr found: " << usage.str());
  }

  unsigned nr    = (unsigned)*(structParser.getFieldAsDouble("nr"));

  double   dr    = (rmax - rmin)/(nr-1);
  
  double*  n0Ptr = MexHandler::createDoubleArray(&plhs[0], nr);
  double*  n1Ptr = MexHandler::createDoubleArray(&plhs[1], nr);
  double*  n2Ptr = MexHandler::createDoubleArray(&plhs[2], nr);

  for(unsigned ir=0; ir < nr; ir++) {
    n0Ptr[ir] = 0;
    n1Ptr[ir] = 0;
    n2Ptr[ir] = 0;
  }

  // The number of iterations

  if(!structParser.fieldExists("niter")) {
    ThrowError("No field niter found: " << usage.str());
  }

  unsigned nIter = (unsigned)*(structParser.getFieldAsDouble("niter"));

  if(!structParser.fieldExists("rsigma")) {
    ThrowError("No field rsigma found: " << usage.str());
  }

  double rsigma   = *(structParser.getFieldAsDouble("rsigma"));

  // Determine whether we are starting from a random seed, or a
  // pre-determined one

  int startSeed = -1;

  if(structParser.fieldExists("seed")) {
    startSeed = (int)*(structParser.getFieldAsDouble("seed"));
  }

  if(startSeed > 0) {
    Sampler::seed((unsigned)startSeed);
  }

  // Determine whether or not to generate a random number of sources,
  // or the exact mean

  bool doRand = false;

  if(structParser.fieldExists("randMean")) {
    doRand = (bool)*(structParser.getFieldAsDouble("randMean"));
  }

  // Now calculate the confusion flux

  double sc = pow((k * dx * dy * (60*60)),(1.0/gamma));

  Flux fluxMin;
  fluxMin.setMilliJy(rmingen * sc);

  Flux fluxMax;
  fluxMax.setMilliJy(rmaxgen * sc);

  std::vector<Flux>  srcFlux;
  std::vector<Angle> srcX;
  std::vector<Angle> srcY;

  unsigned nGrid = nx * ny;
  std::vector<double>   rGrid(nGrid);
  std::vector<double>   gGrid(nGrid);
  std::vector<unsigned> nSrcGrid(nGrid);

  // Main loop -- generate sources over the specified field,
  // accumulate them onto a grid, and count the number of sources in
  // each flux bin

  for(unsigned iIter=0; iIter < nIter; iIter++) {

    if(iIter % 100 == 0) {
      COUT("iter = " << iIter);
    }

    // Generate sources for this iteration

    gen.generateSources(fluxMin, fluxMax, x, y, srcFlux, srcX, srcY, doRand);

    // Initialize the grid, generating gaussiain random noise for each
    // pixel.  Note that rsigma is also a ratio to the confusion flux,
    // so no need to divide by sc here


    gGrid = Sampler::generateGaussianSamples(rsigma, nGrid);

    for(unsigned iGrid=0; iGrid < nGrid; iGrid++) {
      rGrid[iGrid]    = gGrid[iGrid];
      nSrcGrid[iGrid] =            0;
    }

    // Now iterate over sources, and determine which pixel each belongs to

    unsigned nSrc = srcFlux.size();
    unsigned ix, iy, iPix;

    for(unsigned iSrc=0; iSrc < nSrc; iSrc++) {

      ix = (unsigned)floor((srcX[iSrc].degrees() + xdeg/2)/dx);
      iy = (unsigned)floor((srcY[iSrc].degrees() + ydeg/2)/dy);
      
      iPix = iy * nx + ix;
      
      rGrid[iPix]    += srcFlux[iSrc].mJy() / sc;
      nSrcGrid[iPix] += 1;
    }

    // Now iterate over the grid, and determine which bin in r it belongs to

    unsigned ir;

    for(unsigned iGrid=0; iGrid < nGrid; iGrid++) {

      ir = (unsigned)floor((rGrid[iGrid]-rmin)/dr);

      if(ir > 0 && ir < nr-1) {

	if(nSrcGrid[iGrid] == 0) {
	  ++n0Ptr[ir];
	}

	if(nSrcGrid[iGrid] == 1) {
	  ++n1Ptr[ir];
	}

	if(nSrcGrid[iGrid] == 2) {
	  ++n2Ptr[ir];
	}

      }
    }

  }

  return;
}
