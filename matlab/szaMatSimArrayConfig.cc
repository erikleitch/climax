/**.......................................................................
 * MATLAB Mex file for calculating UVW, given an array of antenna
 * locations, a source declination and a list of HA
 *
 * Use like:
 *
 * [uvw chisq] = gcpMatArrayConfig(lat, east(niter, nant), north(niter, nant), ha(nha), dec(ndec), wts(ndec),
 *                                 xprof(nprof), yprof(nprof), nangle, length, minSep, diam, relwts(nparam));
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

//#define DO_SB

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

  MexParser   latParser(prhs[0]);
  MexParser  eastParser(prhs[1]);
  MexParser northParser(prhs[2]);
  MexParser    haParser(prhs[3]);
  MexParser   decParser(prhs[4]);
  MexParser    wtParser(prhs[5]);

  double* latPtr = latParser.getDoubleData();
  double cLat    = cos((*latPtr) / 180 * M_PI); 
  double sLat    = sin((*latPtr) / 180 * M_PI); 

  double* eastPtr  = eastParser.getDoubleData(); 
  double* northPtr = northParser.getDoubleData(); 
  double* haPtr    = haParser.getDoubleData(); 
  double* decPtr   = decParser.getDoubleData();
  double* wtPtr    = wtParser.getDoubleData();
  
  // First dimension of eastPtr should be the number of iterations

  unsigned nIter = eastParser.getDimension(0);

  // Second should be the number of antennas

  unsigned nAnt  = eastParser.getDimension(1);
  unsigned nBase = (nAnt * (nAnt-1))/2;

  // Variables to deal with antenna positions

  std::vector<double> east(nAnt);
  std::vector<double> north(nAnt);
  std::vector<double> x(nAnt);
  std::vector<double> y(nAnt);
  std::vector<double> z(nAnt);
  double dx, dy, dz;

  // Variable to deal with the HA

  unsigned nHa = haParser.getNumberOfElements();
  std::vector<double> cHa(nHa);
  std::vector<double> sHa(nHa);

  for(unsigned iHa=0; iHa < nHa; iHa++) {
    cHa[iHa] = cos(haPtr[iHa] / 12 * M_PI);
    sHa[iHa] = sin(haPtr[iHa] / 12 * M_PI);
  }

  // Precompute cos()/sin() of the declinations

  unsigned nDec = decParser.getNumberOfElements();

  std::vector<double> cDec(nDec);
  std::vector<double> sDec(nDec);

  for(unsigned i=0; i < nDec; i++) {
    cDec[i] = cos(decPtr[i] / 180 * M_PI); 
    sDec[i] = sin(decPtr[i] / 180 * M_PI); 
  }

  // Variables to deal with UV

  std::vector<double>      uVec(nDec * nBase * nHa);
  std::vector<double>      vVec(nDec * nBase * nHa);
  std::vector<double>    uvrVec(nDec * nBase * nHa);
  std::vector<double> shadowVec(nDec * nBase * nHa);
  std::vector<double> shadowAntVec(nDec * nAnt * nHa);

  double*      uPtr =      &uVec[0];
  double*      vPtr =      &vVec[0];
  double*    uvrPtr =    &uvrVec[0];
  double* shadowPtr = &shadowVec[0];
  double* shadowAntPtr = &shadowAntVec[0];

  //------------------------------------------------------------
  // Deal with profile info
  //------------------------------------------------------------

  MexParser xprofParser(prhs[6]);
  double* xprofPtr = xprofParser.getDoubleData();
  double* yprofPtr = MexParser::getDoubleData(prhs[7]);
  unsigned nAngle  = (unsigned)*(MexParser::getDoubleData(prhs[8]));
  double dAng      = M_PI/nAngle;

  double* length   = MexParser::getDoubleData(prhs[9]);
  double* minSep   = MexParser::getDoubleData(prhs[10]);
  double  diam     = *(MexParser::getDoubleData(prhs[11]));

  double* relwts   = MexParser::getDoubleData(prhs[12]);
  unsigned nRelWt  = MexParser::getNumberOfElements(prhs[12]);

  bool doSb = relwts[5] > 0.0;

  double fwhm  = 1.2/diam;
  double sigma = fwhm/2.35 * 180/M_PI;

  double dxprof  = (xprofPtr[1] - xprofPtr[0]);
  unsigned nprof = xprofParser.getNumberOfElements();

  //------------------------------------------------------------
  // The theoretical uvr min/max
  //------------------------------------------------------------
  
  double uvrGlobalMin = *minSep;
  double uvrGlobalMax = (*length) * sqrt(2.0);

  uvrGlobalMin = uvrGlobalMin < xprofPtr[0] ? uvrGlobalMin : xprofPtr[0];
  uvrGlobalMax = uvrGlobalMax > xprofPtr[nprof-1] ? uvrGlobalMax: xprofPtr[nprof-1];
  
  //------------------------------------------------------------
  // Calculate additional elements we need to add to the profile
  // length to encompass the global min/max
  //------------------------------------------------------------
  
  unsigned nxmin = (unsigned)((xprofPtr[0] - uvrGlobalMin)/dxprof);
  unsigned nxmax = (unsigned)((uvrGlobalMax - xprofPtr[nprof-1])/dxprof);
  
  unsigned nNewProf = nprof + nxmin + nxmax;
  
  std::vector<double> xNewProf(nNewProf);
  std::vector<double> yNewProf(nNewProf);

  // Create an array of uvr histograms, corresponding to the number of
  // decs & angles into which we are dividing up the UV plane
  
  std::vector<std::vector<std::vector<double> > > uvrHists(nDec);
  std::vector<std::vector<double> > yHistSums(nDec);

  for(unsigned iDec=0; iDec < nDec; iDec++) {
    uvrHists[iDec].resize(nAngle);
    yHistSums[iDec].resize(nAngle);

    for(unsigned iAng=0; iAng < nAngle; iAng++) {
      uvrHists[iDec][iAng].resize(nNewProf);
      yHistSums[iDec][iAng] = 0.0;
    }
  }

  //------------------------------------------------------------
  // And put the profile in the appropriate bins of the new array
  //------------------------------------------------------------
  
  double xProfMin = xprofPtr[0] - nxmin * dxprof;

  double yProfSum=0.0, yHistSum=0.0;

  for(unsigned i=0; i < nxmin; i++) {
    xNewProf[i] = xprofPtr[0] - (nxmin - i) * dxprof;
    yNewProf[i] = 0.0;
  }

  for(unsigned i=nxmin; i < nprof; i++) {
    xNewProf[i] = xprofPtr[i-nxmin];
    yNewProf[i] = yprofPtr[i-nxmin];
    yProfSum += yNewProf[i];
  }

  for(unsigned i=nprof; i < nNewProf; i++) {
    xNewProf[i] = xprofPtr[nprof-1] + (i-nprof+1) * dxprof;
    yNewProf[i] = 0.0;
  }

  for(unsigned i=0; i < nNewProf; i++) {
    yNewProf[i] /= (yProfSum * dxprof);
  }

  //-----------------------------------------------------------------------
  // Create an array of synthesized beam profiles, corresponding to the
  // number of decs & angles into which we are dividing up the UV plane
  //-----------------------------------------------------------------------
  
  unsigned nPhi = 4;   // Number of points at which to evalulate the
		       // synthesized beam
  unsigned nTheta = 10; // Number of points at which to evalulate the
			 // synthesized beam
  double thetaMax = 0.2; // Degrees
  double dTheta = thetaMax / nTheta;

  std::vector<std::vector<std::vector<double> > > synthBeamCuts(nDec);
  std::vector<std::vector<double> > phis(nDec);
  std::vector<double> bmPa(nDec);
  std::vector<double> bmMin(nDec);
  std::vector<double> bmMaj(nDec);
  std::vector<double> theta(nTheta);

  for(unsigned iDec=0; iDec < nDec; iDec++) {

    phis[iDec].resize(nPhi);
    synthBeamCuts[iDec].resize(nPhi);

    for(unsigned iPhi=0; iPhi < nPhi; iPhi++) {
      synthBeamCuts[iDec][iPhi].resize(nTheta);
    }
  }

  //------------------------------------------------------------
  // Output structure
  //------------------------------------------------------------

  std::vector<int> dims(3);
  dims[0] = 1;
  dims[1] = 1;

  plhs[0] = mxCreateStructArray(2, &dims[0], 0, NULL);

  dims[0] = nDec;
  dims[1] = nBase;
  dims[2] = nHa;

  double*        uOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],       "u", dims);
  double*        vOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],       "v", dims);
  double*      uvrOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],     "uvr", dims);
  double*   shadowOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],  "shadow", dims);

  dims[0] = nDec;
  dims[1] = nAnt;
  dims[2] = nHa;

  double*   shadowAntOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],  "shadowAnt", dims);

  double*     nvisOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],    "nvis",    1);

  double*    uvminOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],   "uvmin",    1);
  double*    uvmaxOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],   "uvmax",    1);
  double*   uvrminOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],  "uvrmin",    1);
  double*   uvrmaxOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],  "uvrmax",    1);
  double*   dxprofOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],  "dxprof",    1);

  double*    xprofOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],   "xprof",   nNewProf);
  double*    yprofOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],   "yprof",   nNewProf);

  double*       haOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],      "ha", nHa);

  for(unsigned iHa=0; iHa < nHa; iHa++)
    haOutPtr[iHa] = haPtr[iHa];

  dims.resize(2);
  dims[0] = nDec;
  dims[1] = nPhi;

  double*  phiOutPtr        = MexHandler::addNamedDoubleStructField(plhs[0], "phi", dims);

  double*  thetaOutPtr      = MexHandler::addNamedDoubleStructField(plhs[0],  "theta",   nTheta);

  for(unsigned iTheta=0; iTheta < nTheta; iTheta++) {
    thetaOutPtr[iTheta] = (dTheta * iTheta)/180 * M_PI;
  }

  *nvisOutPtr   = nBase * nHa;
  *dxprofOutPtr = dxprof;
  
  for(unsigned i=0; i < nNewProf; i++) {
    xprofOutPtr[i] = xNewProf[i];
    yprofOutPtr[i] = yNewProf[i];
  }

  dims.resize(2);
  dims[0] = nIter;
  dims[1] = nDec;

  double* chisqUvOutPtr     = MexHandler::addNamedDoubleStructField(plhs[0], "chisqUv",     dims);
  double* chisqShadowOutPtr = MexHandler::addNamedDoubleStructField(plhs[0], "chisqShadow", dims);
  double* chisqMaxShadowOutPtr = MexHandler::addNamedDoubleStructField(plhs[0], "chisqMaxShadow", dims);
  double* chisqMaxTimeOutPtr = MexHandler::addNamedDoubleStructField(plhs[0], "chisqMaxTimeShadow", dims);
  double* chisqBmOutPtr     = MexHandler::addNamedDoubleStructField(plhs[0], "chisqBm",     dims);
  double* chisqBcOutPtr     = MexHandler::addNamedDoubleStructField(plhs[0], "chisqBc",     dims);

  double* bmPaOutPtr        = MexHandler::addNamedDoubleStructField(plhs[0], "bmPa",        dims);
  double* bmMinOutPtr       = MexHandler::addNamedDoubleStructField(plhs[0], "bmMin",       dims);
  double* bmMajOutPtr       = MexHandler::addNamedDoubleStructField(plhs[0], "bmMaj",       dims);

  double*  decOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],  "decs",   nDec);
  double*  wtOutPtr  = MexHandler::addNamedDoubleStructField(plhs[0],  "decWts",    nDec);
  double*  relwtOutPtr  = MexHandler::addNamedDoubleStructField(plhs[0],  "relWts",    nRelWt);

  for(unsigned iDec=0; iDec < nDec; iDec++) {
    decOutPtr[iDec] = decPtr[iDec];
    wtOutPtr[iDec]  = wtPtr[iDec];
  }

  for(unsigned iWt=0; iWt < nRelWt; iWt++) {
    relwtOutPtr[iWt]  = relwts[iWt];
  }

  mxArray* cmStructPtr = MexHandler::addNamedStructField(plhs[0],       "best");
  double* chisqMinOutPtr = MexHandler::addNamedDoubleStructField(cmStructPtr,"chisqmin",    1);

  dims.resize(3);
  dims[0] = nDec;
  dims[1] = nAngle;
  dims[2] = nNewProf;

  double*  yHistsOutPtr = MexHandler::addNamedDoubleStructField(cmStructPtr,   "yhists", dims);

  dims.resize(3);
  dims[0] = nDec;
  dims[1] = nPhi;
  dims[2] = nTheta;

  double*  synthBeamCutsPtr = MexHandler::addNamedDoubleStructField(cmStructPtr, "synthbeam", dims);

  dims.resize(2);
  dims[0] = 3;
  dims[1] = nAnt;

  double* uenOutPtr = MexHandler::addNamedDoubleStructField(cmStructPtr, "uen",        dims);

  double* bestBmPaOutPtr  = MexHandler::addNamedDoubleStructField(cmStructPtr,    "bmPa",     nDec);
  double* bestBmMinOutPtr = MexHandler::addNamedDoubleStructField(cmStructPtr,    "bmMin",    nDec);
  double* bestBmMajOutPtr = MexHandler::addNamedDoubleStructField(cmStructPtr,    "bmMaj",    nDec);

  double* bestChisqUvOutPtr     = MexHandler::addNamedDoubleStructField(cmStructPtr, "chisqUv",     nDec);
  double* bestChisqShadowOutPtr = MexHandler::addNamedDoubleStructField(cmStructPtr, "chisqShadow", nDec);
  double* bestChisqMaxShadowOutPtr = MexHandler::addNamedDoubleStructField(cmStructPtr, "chisqMaxShadow", nDec);
  double* bestChisqMaxTimeOutPtr = MexHandler::addNamedDoubleStructField(cmStructPtr, "chisqMaxTimeShadow", nDec);
  double* bestChisqBmOutPtr     = MexHandler::addNamedDoubleStructField(cmStructPtr, "chisqBm",     nDec);
  double* bestChisqBcOutPtr     = MexHandler::addNamedDoubleStructField(cmStructPtr, "chisqBc",     nDec);

  //------------------------------------------------------------
  // Loop over iterations
  //------------------------------------------------------------

  unsigned posInd, histInd, chisqInd, bmInd, phiInd;
  double u2mean, v2mean, uvmean;
  double bmfac;

  unsigned nmean;

  double chisqMin;
  unsigned dIter = nIter/10 > 1000 ? 1000 : nIter/10;

  if(dIter == 0)
    dIter = 100;

  double wtSum = 0.0;
  for(unsigned iDec = 0; iDec < nDec; iDec++)
    for(unsigned iWt = 0; iWt < nRelWt; iWt++) {
      wtSum += wtPtr[iDec] * relwts[iWt];
    }

  std::vector<double> chisqUvDec(nDec);
  std::vector<double> chisqShadowDec(nDec);
  std::vector<double> chisqMaxShadowDec(nDec);
  std::vector<double> chisqMaxTimeShadowDec(nDec);
  std::vector<double> chisqBmDec(nDec);
  std::vector<double> chisqBcDec(nDec);

  for(unsigned iIter=0; iIter < nIter; iIter++) {
    
    if(iIter % dIter == 0) {
      COUT(iIter);
    }

    // Convert antenna positions from UEN to XYZ for this iteration

    for(unsigned iAnt=0; iAnt < nAnt; iAnt++) {
      posInd = iAnt * nIter + iIter;

      x[iAnt] = -sLat * northPtr[posInd];
      y[iAnt] =          eastPtr[posInd];
      z[iAnt] =  cLat * northPtr[posInd];
    }

    // Loop over declinations

    double chisq = 0.0;
    double umin, umax;
    double vmin, vmax;
    double uvmin, uvmax;
    double uvrmin, uvrmax;

    for(unsigned iDec=0; iDec < nDec; iDec++) {

      chisqUvDec[iDec]            = 0.0;
      chisqShadowDec[iDec]        = 0.0;
      chisqMaxShadowDec[iDec]     = 0.0;
      chisqMaxTimeShadowDec[iDec] = 0.0;

      double shadowFrac    = 0.0;
      double maxShadow     = 0.0;
      double maxTimeShadow = 0.0;

      bmInd = iDec * nIter + iIter;

      // Now calculate uv
      
      unsigned iBase = 0;

      double u,v,uvr;
      unsigned uvInd;
      
      // Zero the UVR histograms for this source

      for(unsigned iAng=0; iAng < nAngle; iAng++) {
	yHistSums[iDec][iAng] = 0.0;
	for(unsigned iHist=0; iHist < nNewProf; iHist++) {
	  uvrHists[iDec][iAng][iHist] = 0.0;
	}
      }

      // Zero the beam cuts for this source

      for(unsigned iPhi=0; iPhi < nPhi; iPhi++) {
	for(unsigned iTheta=0; iTheta < nTheta; iTheta++) {
	  synthBeamCuts[iDec][iPhi][iTheta] = 0.0;
	}
      }

      u2mean=0.0;
      v2mean=0.0;
      uvmean=0.0;
      nmean=0;

      double w1,w2;
      unsigned iAntShadowed, antInd;
      for(unsigned iAnt1=0; iAnt1 < nAnt-1; iAnt1++)
	for(unsigned iAnt2=iAnt1+1; iAnt2 < nAnt; iAnt2++) {
	  
	  dx = x[iAnt1] - x[iAnt2];
	  dy = y[iAnt1] - y[iAnt2];
	  dz = z[iAnt1] - z[iAnt2];

	  double nShadow=0.0;
	  for(unsigned iHa=0; iHa < nHa; iHa++) {

	    uvInd = (iHa * nBase + iBase) * nDec + iDec;
	    
	    u =               sHa[iHa] * dx +              cHa[iHa] * dy;
	    v = -sDec[iDec] * cHa[iHa] * dx + sDec[iDec] * sHa[iHa] * dy + cDec[iDec] * dz;

	    //	    w1 = x[iAnt1]*cDec[iDec]*cHa[iHa] - y[iAnt1]*cDec[iDec]*sHa[iHa] + z[iAnt1]*sDec[iDec];
	    //	    w2 = x[iAnt2]*cDec[iDec]*cHa[iHa] - y[iAnt2]*cDec[iDec]*sHa[iHa] + z[iAnt2]*sDec[iDec];

	    uvr = sqrt(u * u + v * v);
	    
	    // Accumulate moments needed to calculate the orientation
	    // of the beam

	    u2mean += (u*u - u2mean)/(nmean+1);
	    v2mean += (v*v - v2mean)/(nmean+1);
	    uvmean += (u*v - uvmean)/(nmean+1);
	    ++nmean;

	    if(iDec==0 && iBase==0 && iHa==0) {

	      umin = u;
	      umax = u;
	      vmin = v;
	      vmax = v;
	      uvrmin = uvr;
	      uvrmax = uvr;

	    } else {

	      umin = umin < u ? umin : u;
	      umax = umax > u ? umax : u;

	      vmin = vmin < v ? vmin : v;
	      vmax = vmax > v ? vmax : v;

	      uvrmin = uvrmin < uvr ? uvrmin : uvr;
	      uvrmax = uvrmax > uvr ? uvrmax : uvr;

	    }

	    uvmin = umin < vmin ? umin : vmin;
	    uvmax = umax > vmax ? umax : vmax;

	    // Store the uv/r information
	    
	    uPtr[uvInd]   =   u;
	    vPtr[uvInd]   =   v;
	    uvrPtr[uvInd] = uvr;

	    // Store the percent shadowing for this point

	    shadowPtr[uvInd] = (uvr > diam ? 0 : (diam-uvr)/diam);
	    shadowFrac += shadowPtr[uvInd];

	    maxShadow = (shadowPtr[uvInd] > maxShadow) ? shadowPtr[uvInd] : maxShadow;

	    nShadow += (shadowPtr[uvInd] > 0.0 ? 1.0 : 0.0);

	    // If the current baseline is shadowed for this point,
	    // figure out which antenna is shadowed

	    //	    iAntShadowed = (abs(w1) > abs(w2)) ? iAnt1 : iAnt2;
	    //	    antInd = (iHa * nAnt + iAntShadowed) * nDec + iDec;
	    //	    shadowAntPtr[antInd] = shadowPtr[uvInd];

	    // Get the histogram index this point corresponds to

	    histInd = (unsigned)((uvr - xProfMin)/dxprof);
	    
	    // And which histogram it belongs to. Fold the data into
	    // the +ive uv-halfplane, so that we only need to test
	    // angles between 0 and 180

	    if(v < 0) {
	      u = -u;
	      v = -v;
	    }

	    double theta = atan2(v,u);
	    int iHist = (int)(theta/dAng);

	    if(iHist < 0 || iHist > nAngle-1)
	      iHist = 0;

	    if(histInd < 0)
	      histInd = 0;
	    
	    if(histInd > nNewProf-1)
	      histInd = nNewProf-1;

#if 1

	    uvrHists[iDec][iHist][histInd] += 1.0;
	    yHistSums[iDec][iHist]         += 1.0;

#endif

	  }

	  nShadow /= nHa;

	  maxTimeShadow = (maxTimeShadow > nShadow) ? maxTimeShadow : nShadow;

	  ++iBase;

	}

      //------------------------------------------------------------
      // Calculate beam symmetry
      //------------------------------------------------------------

      bmPa[iDec] = M_PI + 0.5 * atan2(2*uvmean, u2mean - v2mean);
      double ftmp = sqrt((u2mean - v2mean)*(u2mean - v2mean) + 4.0 * uvmean*uvmean);

      bmMin[iDec] = 0.7/(sqrt(2.0*(u2mean+v2mean) + 2.0*ftmp));
      bmMaj[iDec] = 0.7/(sqrt(2.0*(u2mean+v2mean) - 2.0*ftmp));
      
      bmPaOutPtr[bmInd]  = bmPa[iDec];
      bmMinOutPtr[bmInd] = bmMin[iDec];
      bmMajOutPtr[bmInd] = bmMaj[iDec];

      //------------------------------------------------------------
      // Calculate sidelobe levels
      //------------------------------------------------------------

      if(doSb) {

	double sigAv = sqrt(bmMin[iDec]*bmMaj[iDec])/2.35;
	chisqBcDec[iDec] = 0.0;

	for(unsigned iPhi=0; iPhi < nPhi; iPhi++) {

	  phiInd = iDec * nPhi + iPhi;
	  phiOutPtr[phiInd] = bmPa[iDec] + M_PI/4 * iPhi;
	  
	  double cp = cos(phiOutPtr[phiInd]);
	  double sp = sin(phiOutPtr[phiInd]);
	  
	  double sbNorm = 0.0, sb;
	  double chisqPhiBc = 0.0;
	  for(unsigned iTheta=0; iTheta < nTheta; iTheta++) {

	    for(unsigned iBase=0; iBase < nBase; iBase++) {
	      for(unsigned iHa=0; iHa < nHa; iHa++) {
		uvInd = (iHa * nBase + iBase) * nDec + iDec;
		u = uPtr[uvInd];
		v = vPtr[uvInd];
		sb = cos(2*M_PI*thetaOutPtr[iTheta]*(u*cp + v*sp))/(nBase*nHa);
		synthBeamCuts[iDec][iPhi][iTheta] += sb;
		sbNorm += sb;
	      }
	    }

	    if(thetaOutPtr[iTheta] > 3*sigAv) {
	      sb = synthBeamCuts[iDec][iPhi][iTheta];	    
	      chisqPhiBc += sb * sb;
	    }
	    
	  }

	  // Normalize so that integral of the synthesized beam is 1

	  chisqPhiBc /= (sbNorm * sbNorm);

	  // And give equal weight to each cut (nPhi)

	  chisqBcDec[iDec] += chisqPhiBc/nPhi;
	}
      }
	    
      //------------------------------------------------------------
      // Now we have the uvr histograms.  Calculate the delta from the
      // profile
      //------------------------------------------------------------

      double delta,d2;

      for(unsigned iAng=0; iAng < nAngle; iAng++) {
	std::vector<double>& hist = uvrHists[iDec][iAng];
	std::vector<double>& histSum = yHistSums[iDec];
	
	for(unsigned iBin=0; iBin < nNewProf; iBin++) {
	  hist[iBin] /= (histSum[iAng] * dxprof);
	  delta = (hist[iBin] - yNewProf[iBin]) * dxprof;
	  d2 = delta * delta;
	  chisqUvDec[iDec] += d2;
	}
      }

      // Store the contribution to the chi-squared from the uvr
      // profile condition, for this source only
      
      chisqInd = iDec * nIter + iIter;
      chisqUvOutPtr[chisqInd] = chisqUvDec[iDec];
      
      // Calculate the contribution to chi-squared from the percent
      // baseline shadowing condition

      chisqShadowDec[iDec] = shadowFrac / (nBase * nHa);
      chisqShadowOutPtr[chisqInd] = chisqShadowDec[iDec];

      // Calculate the contribution to chi-squared from the max shadow
      // condition

      chisqMaxShadowDec[iDec] = maxShadow;
      chisqMaxShadowOutPtr[chisqInd] = chisqMaxShadowDec[iDec];

      // Calculate the contribution to chi-squared from the max shadow
      // condition

      chisqMaxTimeShadowDec[iDec] = maxTimeShadow;
      chisqMaxTimeOutPtr[chisqInd] = chisqMaxTimeShadowDec[iDec];

      // Calculate the contribution to chi-squared from the beam
      // symmetry condition

      bmfac = (2.0/(1 + bmMin[iDec]/bmMaj[iDec]) - 1); // Tends to 1 if bmin >> bmaj
					               // or bmin >> bmaj, and 0 if
					               // bmin ~ bmaj

      chisqBmDec[iDec] = bmfac*bmfac;
      chisqBmOutPtr[chisqInd] = chisqBmDec[iDec];

      chisqBcOutPtr[chisqInd] = chisqBcDec[iDec];

      // And increment the global chi-squared, weighted appropriately
      // for this declination
      
      chisq += (
		wtPtr[iDec] * (relwts[0] * chisqUvDec[iDec] + 
			       relwts[1] * chisqShadowDec[iDec] + 
			       relwts[2] * chisqMaxShadowDec[iDec] + 
			       relwts[3] * chisqMaxTimeShadowDec[iDec] + 
			       relwts[4] * chisqBmDec[iDec] + 
			       relwts[5] * chisqBcDec[iDec])
		)/wtSum;

    } // End loop over decs


    if(iIter==0 || chisq < chisqMin) {
      
      for(unsigned i=0; i < nDec*nBase*nHa; i++ ) {
	     uOutPtr[i] =      uPtr[i];
	     vOutPtr[i] =      vPtr[i];
	   uvrOutPtr[i] =    uvrPtr[i];
	shadowOutPtr[i] = shadowPtr[i];
      }

      for(unsigned i=0; i < nDec*nAnt*nHa; i++ ) {
	shadowAntOutPtr[i] = shadowAntPtr[i];
      }

      for(unsigned iAnt=0; iAnt < nAnt; iAnt++) {
	posInd = iAnt * nIter + iIter;

	uenOutPtr[iAnt*3]   =              0.0;
	uenOutPtr[iAnt*3+1] =  eastPtr[posInd];
	uenOutPtr[iAnt*3+2] = northPtr[posInd];
      }

      unsigned ind;
      for(unsigned iDec=0; iDec < nDec; iDec++) {

	for(unsigned iAng=0; iAng < nAngle; iAng++) {
	  for(unsigned iBin=0; iBin < nNewProf; iBin++) {
	    ind = (iBin * nAngle + iAng) * nDec + iDec;
	    yHistsOutPtr[ind] = uvrHists[iDec][iAng][iBin];
	  }
	}

	for(unsigned iPhi=0; iPhi < nPhi; iPhi++) {
	  for(unsigned iTheta=0; iTheta < nTheta; iTheta++) {
	    ind = (iTheta * nPhi + iPhi) * nDec + iDec;
	    synthBeamCutsPtr[ind] = synthBeamCuts[iDec][iPhi][iTheta];
	  }
	}

	bestBmPaOutPtr[iDec]  =  bmPa[iDec];
	bestBmMinOutPtr[iDec] = bmMin[iDec];
	bestBmMajOutPtr[iDec] = bmMaj[iDec];

	bestChisqUvOutPtr[iDec]        =  chisqUvDec[iDec];
	bestChisqShadowOutPtr[iDec]    =  chisqShadowDec[iDec];
	bestChisqMaxShadowOutPtr[iDec] =  chisqMaxShadowDec[iDec];
	bestChisqMaxTimeOutPtr[iDec]   =  chisqMaxTimeShadowDec[iDec];
	bestChisqBmOutPtr[iDec]        =  chisqBmDec[iDec];
      	bestChisqBcOutPtr[iDec]        =  chisqBcDec[iDec];
      }

      *chisqMinOutPtr  = chisq;

      *uvminOutPtr  = uvmin;
      *uvmaxOutPtr  = uvmax;
      *uvrminOutPtr = uvrmin;
      *uvrmaxOutPtr = uvrmax;

      chisqMin = chisq;
    }

  } // End loop over iterations
  
  return;
}
