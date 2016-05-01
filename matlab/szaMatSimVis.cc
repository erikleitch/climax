/**.......................................................................
 * MATLAB Mex file for returning simulated visibilities, given a list
 * of (u, v, and frequencies)
 *
 * Use like:
 *
 * d=gcpMatSimVis(u, v, vis(complex), var(complex));
 *
 */
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

#include "gcp/array/code/share/slalib/slalib.h"

#include "gcp/fftutil/ImGen.h"
#include "gcp/fftutil/SimVis.h"

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

  std::ostringstream os;

  os << "Usage: gcpMatSimVis(ad, u[n x nBase x nFreq], v[n x nBase x nFreq], vis[n x nBase x nFreq], " 
     << std::endl
     << "var[n x nBase x nFreq])" << std::endl << std::endl
     << "Where ad should contain fields:" << std::endl << std::endl
     << "   nu(1xNBAND)" << std::endl
     << "   Field_size_deg(double)" << std::endl
     << "   N_pix(double)" << std::endl
     << "   beam(N_pix x N_pix x NBAND)" << std::endl << std::endl
     << " " << std::endl << std::endl
     << "and u,v,vis and var should all be the same length" << std::endl << std::endl
     << "Note that vis will be overwritten by this function" << std::endl << std::endl;

  //------------------------------------------------------------
  // Input argument parsing
  //------------------------------------------------------------

  MexParser adParser(prhs[0]);
  MexParser uParser(prhs[1]);
  MexParser vParser(prhs[2]);
  MexParser visParser(prhs[3]);

  // Check that the first argument is a struct

  if(!adParser.isStruct()) {
    ThrowError(os.str());
  }

  // Check that the u,v and vis array dimensions match

  if(!MexParser::dimensionsMatch(prhs[1], prhs[2]) || 
     !MexParser::dimensionsMatch(prhs[1], prhs[3])) {
    ThrowError("Dimensions of the u,v,vis arrays must match" << std::endl << os.str());
  }

  // Check that the vis argument is complex

  if(!visParser.isComplex()) {
    ThrowError("The vis array must be complex" << std::endl << os.str());
  }

  double* reVarPtr=0;
  double* imVarPtr=0;

  bool useVar = false;

  if(nrhs > 4) {

    MexParser varParser(prhs[4]);
    if(!varParser.isComplex()) {
      ThrowError("The var array must be complex" << std::endl << os.str());
    }

    reVarPtr = varParser.getDoubleData();
    imVarPtr = varParser.getImagDoubleData();

    COUT("Using the var array to generate noise");
    useVar = true;
  }
    
  //------------------------------------------------------------
  // Now get data
  //------------------------------------------------------------

  Angle size;
  size.setDegrees(*(adParser.getFieldAsDouble("Field_size_deg")));
  COUT("Generating images that are: " << size << " degrees on a side"); 

  unsigned npix = (unsigned)*(adParser.getFieldAsDouble("N_pix"));
  COUT("Generating images that are: " << npix << " x " << npix);

  MexParser nuParser(adParser.getField("nu"));
  MexParser beamParser(adParser.getField("beam"));

  if(beamParser.getDimension(0) != npix || beamParser.getDimension(1) != npix || 
     beamParser.getDimension(2) != nuParser.getDimension(1)) {
    beamParser.printDimensions();
    nuParser.printDimensions();
    ThrowError("Third dimension of field 'beam' must be the same as the second dimension of field 'nu'");
  }

  double* nuPtr = nuParser.getDoubleData();
  double* bmPtr = beamParser.getDoubleData();

  unsigned nData = uParser.getDimension(0);
  unsigned nBase = uParser.getDimension(1);

  double* uPtr  = uParser.getDoubleData(prhs[1]);
  double* vPtr  = vParser.getDoubleData(prhs[2]);

  double* rePtr = (double*)mxGetData(prhs[3]);
  double* imPtr = (double*)mxGetImagData(prhs[3]);

  //------------------------------------------------------------
  // Generate CMB image to the same size as the beam
  //------------------------------------------------------------

  gcp::util::ImGen imgen;

  imgen.setImageSize(size);
  imgen.setNpix(npix);

  //------------------------------------------------------------
  // Get the power spectrum normalization, in micro-Kelvin
  // matter
  //------------------------------------------------------------

  gcp::util::Temperature tnorm;

  if((adParser.getField("Tnorm_uK")) != 0) {
    tnorm.setK((*adParser.getFieldAsDouble("Tnorm_uK")) * 1e-6);
  } else {
    tnorm.setK(100e-6);
  }

  //------------------------------------------------------------
  // Get the l at which to normalize the power spectrum
  //------------------------------------------------------------

  double lnorm;

  if((adParser.getField("lnorm")) != 0) {
    lnorm = *adParser.getFieldAsDouble("lnorm");
  } else {
    lnorm = 200;
  }
  
  COUT("Normalizing power spectrum to: " << tnorm << " at l = " << lnorm);
  imgen.setPowerSpectrumNormalization(lnorm, tnorm);

  //------------------------------------------------------------
  // Get the type of power spectrum to use
  //------------------------------------------------------------

  std::string type;

  if((adParser.getField("Type")) != 0) {
    type = adParser.getFieldAsString("Type");
  } else {
    type = "pow";
  }
  
  if(type == "pow") {
    COUT("Generating power law spectrum");
    imgen.setPowerSpectrumType(ImGen::TYPE_POWER_LAW);
  } else if(type == "gauss") {
    COUT("Generating Gaussian spectrum");
    imgen.setPowerSpectrumType(ImGen::TYPE_GAUSSIAN);
  } else {
    ThrowError("Unrecognized power spectrum type: " << type);
  }

  //------------------------------------------------------------ 
  // In the case of power spectrum normalization, get the power
  // spectrum index to use
  //------------------------------------------------------------

  double powInd;

  if((adParser.getField("powInd")) != 0) {
    powInd = *adParser.getFieldAsDouble("powInd");
  } else {
    powInd = 0.0;
  }

  if(type == "pow") {
    COUT("Using power law index: " << powInd);
  }

  imgen.setPowInd(powInd);

  //------------------------------------------------------------ 
  // In the case of gaussian power spectrum, get the fwhm (in l)
  // to use
  //------------------------------------------------------------

  double fwhm;

  if((adParser.getField("Fwhm")) != 0) {
    fwhm = *adParser.getFieldAsDouble("Fwhm");
  } else {
    fwhm = 0.0;
  }

  if(type == "gauss") {
    COUT("Using Gaussian fwhm: " << fwhm);
  }

  imgen.setFwhm(fwhm);

  //------------------------------------------------------------
  // If we were told to seed the image with an explicit seed, do so now
  //------------------------------------------------------------

  if(adParser.fieldExists("imSeed")) {
    unsigned int s = (unsigned int)*adParser.getFieldAsDouble("imSeed");
    imgen.seed(s);
  } else {
    imgen.seedRandom();
  }

  //------------------------------------------------------------
  // Finally, generate the image
  //------------------------------------------------------------

  imgen.generateImages("I");

  //------------------------------------------------------------
  // Insert this in the simulation container as the CMB image
  //------------------------------------------------------------

  SimVis sim;
  sim.setImageSize(size);
  sim.setImageNpix(npix);

#if 0

  // Substitute a sinusoid for the actual image

  double dx = size.radians()/npix;
  double dy = size.radians()/npix;
  double x ,y,ufreq = 5000,vfreq = 1000;
  unsigned ind;
  for(unsigned j=0; j < npix; j++) {
    y = j * dy;
    for(unsigned i=0; i < npix; i++) {
      x = i * dx;
      ind = i + j*npix;
      imgen.images_[0].data_[ind] = sin(2*M_PI*ufreq*x) + sin(2*M_PI*vfreq*y);
    }
  }

#else
  unsigned ind;
  for(unsigned j=0; j < npix; j++) {
    for(unsigned i=0; i < npix; i++) {
      ind = i + j*npix;
      imgen.images_[0].data_[ind] = 0.0;
    }
  }
#endif

#if 1
  {
    // Write out this image
    
    std::vector<int> dims(2);
    dims[0] = npix;
    dims[1] = npix;
    
    double* imagePtr = MexHandler::createDoubleArray(plhs, 2, &dims[0]);
    
    unsigned imInd;
    for(unsigned j=0; j < npix; j++) {
      for(unsigned i=0; i < npix; i++) {
	imInd = i + j*npix;
	imagePtr[imInd] = imgen.images_[0].data_[imInd];
      }
    }
  }
#endif

  sim.installCmbImage(imgen.images_[0], Stokes::STOKES_I);

  //------------------------------------------------------------
  // Install frequencies
  //------------------------------------------------------------

  unsigned nFreq = nuParser.getDimension(1);
  std::vector<Frequency> freq(nFreq);

  for(unsigned i=0; i < nFreq; i++) {
    freq[i].setHz(nuPtr[i]);
  }

  sim.setFrequencies(freq);

  //------------------------------------------------------------
  // Now install the beams
  //------------------------------------------------------------

  SimVis::SimImage bm;

  unsigned ndata = npix*npix;
  for(unsigned iFreq=0; iFreq < nFreq; iFreq++) {
    bm.initialize(&bmPtr[iFreq*ndata], ndata, SimVis::SimImage::UNITS_JY);
    sim.installBeamImage(bm, iFreq);
  }

  //------------------------------------------------------------ 
  // Check if the user specified a visibility rms to use
  //------------------------------------------------------------

  double rms;

  if((adParser.getField("Rms_Jy")) != 0) {
    rms = *adParser.getFieldAsDouble("Rms_Jy");

    if(!useVar) {
      COUT("Using rms = " << rms << " to generate noise");
    }

  } else {
    rms = 0.0;
  }

  //------------------------------------------------------------
  // If we were told to seed the noise with an explicit seed, do so
  // now
  //------------------------------------------------------------

  if(adParser.fieldExists("noiseSeed")) {
    unsigned int s = (unsigned int)*adParser.getFieldAsDouble("noiseSeed");
    imgen.seed(s);
  } else {
    imgen.seedRandom();
  }

  //-----------------------------------------------------------------------
  // Now do the calculation, overwriting the visibilities
  //-----------------------------------------------------------------------

  SimVis::FT ft;

  bool bad;
  unsigned nPrint=10, iPrint=0;
  double reNoise, imNoise;
  for(unsigned iFreq=0; iFreq < nFreq; iFreq++) {

    // Compute the transform for this frequency

    COUT("Calculating transform for frequency: " << freq[iFreq]);

    sim.observeImage(iFreq, ft, Stokes::STOKES_I);

    double* u  = uPtr  + iFreq*nData*nBase;
    double* v  = vPtr  + iFreq*nData*nBase;
    double* re = rePtr + iFreq*nData*nBase;
    double* im = imPtr + iFreq*nData*nBase;

    double* reVar = (reVarPtr) ? reVarPtr + iFreq*nData*nBase : 0;
    double* imVar = (imVarPtr) ? imVarPtr + iFreq*nData*nBase : 0;

    // Now calculate visibilities from the transform we just computed,
    // overwriting the re and im arrays in the process

    for(unsigned iVis=0; iVis < nData*nBase; iVis++) {

      sim.getVis(ft, u[iVis], v[iVis], re[iVis], im[iVis], bad);

      // Add random noise if requested

      if(useVar || rms > 0.0) {
	double reRms = useVar ? sqrt(reVar[iVis]) : rms;
	double imRms = useVar ? sqrt(imVar[iVis]) : rms;

	reNoise = ImGen::gauss_rand(reRms);
	imNoise = ImGen::gauss_rand(imRms);

	re[iVis] += reNoise;
	im[iVis] += imNoise;

	if(iPrint < nPrint) {
	  COUT("Re noise was: " << reNoise << " Im noise was: " << imNoise);
	}

	++iPrint;
      }

      // Set to NaN (for matlab) if u,v value couldn't be interpolated

      if(bad) {
	re[iVis] = numeric_limits<double>::quiet_NaN();
	im[iVis] = numeric_limits<double>::quiet_NaN();
      }

    }
  }

  return;
}
