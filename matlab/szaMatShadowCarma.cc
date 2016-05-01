/**.......................................................................
 * MATLAB Mex file for calculating shadowing, given a specified GCP
 * and CARMA array
 *
 * Use like:
 *
 * shadowed(bool[Nobs x Nant]) = gcpMatShadowCarma(az(double[Nobs x Nant]), el(double[Nobs x Nant]), 
 *                                                 gcpConf(string), carmaConf(string))
 *
 */
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

#include "gcp/array/code/share/slalib/slalib.h"

#include "gcp/util/CarmaConfig.h"

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

  os << "Usage: " << std::endl
     << "   shadowed = gcpMatShadowCarma(double az[Nobs x 8], double el[Nobs x 8], " << std::endl
     << "                                string gcpConf, string carmaConf)" << std::endl
     << "Where:" << std::endl
     << " az should be in degrees" << std::endl
     << " el should be in degrees" << std::endl
     << " gcpConf should be one of 'I' (imaging), 'H' (high-dec), or 'L' (low-dec)" << std::endl
     << " carmaConf should be one of 'A', 'B', 'C', 'D', or 'E'" << std::endl
     << std::endl;

  //------------------------------------------------------------
  // Input argument parsing
  //------------------------------------------------------------

  MexParser azParser(prhs[0]);
  MexParser elParser(prhs[1]);
  MexParser gcpConfParser(prhs[2]);
  MexParser carmaConfParser(prhs[3]);

  // Check that the az & el dimensions match

  if(!MexParser::dimensionsMatch(prhs[0], prhs[1])) {
    ThrowError("Dimensions of the az & el arrays must match" << std::endl << os.str());
  }

  // Check that the az/el arrays have 8 telescopes as their last dimension

  unsigned nDim = azParser.getNumberOfDimensions();
  if(azParser.getDimension(nDim-1) != 8) {
    ThrowError("The az & el arrays must be [Nobs x 8] for the GCP" << std::endl << os.str());
  }

  double* azPtrBase = azParser.getDoubleData();
  double* elPtrBase = elParser.getDoubleData();

  unsigned nObs = azParser.getNumberOfElements()/8;
  unsigned nSza = 8;

  //------------------------------------------------------------
  // Create output vector to have the same dimensions as the input
  //------------------------------------------------------------

  plhs[0] = MexHandler::createMatlabArray(nDim, azParser.getDimensions(), DataType::BOOL);
  unsigned char* shadowBasePtr = (unsigned char*)mxGetData(plhs[0]);

  //------------------------------------------------------------  
  // Now get the configurations
  //------------------------------------------------------------  

  CarmaConfig cc;
  std::vector<CarmaConfig::PadLocation> carma;

  if(carmaConfParser.getString() == "T") {
    cc.initializeCurrentConfiguration();
    cc.addPad(29, CarmaConfig::BIMA);

    carma = cc.getCurrentConfiguration();
  } else {
    carma = cc.getConfiguration(carmaConfParser.getString());
  }

#if 0
  COUT("Carma configuration is: " << std::endl);

  for(unsigned i=0; i < carma.size(); i++) {
    COUT(carma[i]);
  }

  COUT("");
#endif

  std::vector<CarmaConfig::PadLocation> gcp   = cc.getConfiguration(gcpConfParser.getString());
  COUT("Abotu to sort");
  gcp   = cc.sortByAntNumber(gcp);
  COUT("Abotu to sort: done");

#if 0
  COUT("GCP configuration is: " << std::endl);
  for(unsigned i=0; i < gcp.size(); i++) {
    COUT(gcp[i]);
  }
  COUT("");
#endif

  //------------------------------------------------------------  
  // Finally, loop through observations, checking shadowing
  //------------------------------------------------------------  

  Angle az, el;

  // Iterate over GCP antennas

  bool shadowed;
  for(unsigned iSza=0; iSza < nSza; iSza++) {

    double* azPtr = azPtrBase + nObs * iSza;
    double* elPtr = elPtrBase + nObs * iSza;
    unsigned char* shadowPtr = shadowBasePtr + nObs * iSza;
 
    CarmaConfig::PadLocation& gcpAnt = gcp[iSza];
    
    // Iterate over observations for this antenna

    for(unsigned iObs=0; iObs < nObs; iObs++) {

      az.setDegrees(azPtr[iObs]);
      el.setDegrees(elPtr[iObs]);

      // Iterate over CARMA pads for this antenna, this observation

      shadowed = false;
      for(unsigned iCarma=0; iCarma < carma.size(); iCarma++) {
	shadowed |= gcpAnt.isShadowed(az, el, carma[iCarma]);
      }

      shadowPtr[iObs] = shadowed;
    }
  }

  return;
}
