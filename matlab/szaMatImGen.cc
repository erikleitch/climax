/**.......................................................................
 * MATLAB Mex file for returning simulated images
 *
 * Use like:
 *
 * im = gcpMatImGen(ad)
 *
 */
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

#include "gcp/array/code/share/slalib/slalib.h"

#include "gcp/fftutil/ImGen.h"

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

  os << "Usage: im = gcpMatSimVis(ad) " 
     << std::endl
     << "Where ad should contain fields:" << std::endl << std::endl
     << "   Field_size_deg(double)" << std::endl
     << "   N_pix(double)" << std::endl
     << "   beam(N_pix x N_pix x NBAND)" << std::endl << std::endl
     << " " << std::endl << std::endl << std::endl;

  //------------------------------------------------------------
  // Input argument parsing
  //------------------------------------------------------------

  MexParser adParser(prhs[0]);

  // Check that the first argument is a struct

  if(!adParser.isStruct()) {
    ThrowError(os.str());
  }

  //------------------------------------------------------------
  // Now get data
  //------------------------------------------------------------

  Angle size;
  size.setDegrees(*(adParser.getFieldAsDouble("Field_size_deg")));
  COUT("Generating images that are: " << size << " degrees on a side"); 

  unsigned npix = (unsigned)*(adParser.getFieldAsDouble("N_pix"));
  COUT("Generating images that are: " << npix << " x " << npix);

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
  // Finally, generate the image
  //------------------------------------------------------------

  imgen.generateImages("I");

  // Write out this image

  std::vector<int> dims(2);
  dims[0] = npix;
  dims[1] = npix;

  double* imPtr = MexHandler::createDoubleArray(plhs, 2, &dims[0]);

  unsigned ind;
  for(unsigned j=0; j < npix; j++) {
    for(unsigned i=0; i < npix; i++) {
      ind = i + j*npix;
      imPtr[ind] = imgen.images_[0].data_[ind];
    }
  }
}
