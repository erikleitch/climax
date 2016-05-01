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
     << "       " << std::endl
     << "   shadowed = gcpMatShadowCarma(double az[nObs x nArray1], double el[nObs x nArray1], " << std::endl
     << "                                double arrayPads1, double arrayAntTypes1, double arrayPads2, double arrayAntTypes2)" << std::endl
     << "Where:" << std::endl
     << "      " << std::endl
     << " az should be in degrees" << std::endl
     << " el should be in degrees" << std::endl
     << " gcpConf should be one of 'I' (imaging), 'H' (high-dec), or 'L' (low-dec)" << std::endl
     << " carmaConf should be one of 'A', 'B', 'C', 'D', or 'E'" << std::endl
     << std::endl;

  //------------------------------------------------------------
  // Input argument parsing
  //------------------------------------------------------------

  if(nrhs < 6) {
    ThrowError(os.str());
  }

  MexParser azParser(prhs[0]);
  MexParser elParser(prhs[1]);
  MexParser arrayPads1Parser(prhs[2]);
  MexParser arrayAntTypes1Parser(prhs[3]);
  MexParser arrayPads2Parser(prhs[4]);
  MexParser arrayAntTypes2Parser(prhs[5]);

  // Check that the az & el dimensions match

  if(!MexParser::dimensionsMatch(prhs[0], prhs[1])) {
    ThrowError("Dimensions of the az & el arrays must match" << std::endl << os.str());
  }

  double* arrayPads1Ptr     = arrayPads1Parser.getDoubleData();
  double* arrayAntTypes1Ptr = arrayAntTypes1Parser.getDoubleData();
  unsigned nArray1          = arrayPads1Parser.getNumberOfElements();

  // Check that the az/el arrays have nArray1 telescopes as their last dimension

  unsigned nDim = azParser.getNumberOfDimensions();
  if(azParser.getDimension(nDim-1) != nArray1) {
    ThrowError("The az & el arrays must be [nObs x nArray1]" << std::endl << os.str());
  }

  double* arrayPads2Ptr     = arrayPads2Parser.getDoubleData();
  double* arrayAntTypes2Ptr = arrayAntTypes2Parser.getDoubleData();
  unsigned nArray2          = arrayPads2Parser.getNumberOfElements();

  double* azPtrBase = azParser.getDoubleData();
  double* elPtrBase = elParser.getDoubleData();

  unsigned nObs = azParser.getNumberOfElements() / nArray1;

  //------------------------------------------------------------
  // Create output vector to have the same dimensions as the input
  //------------------------------------------------------------

  plhs[0] = MexHandler::createMatlabArray(nDim, azParser.getDimensions(), DataType::BOOL);
  unsigned char* shadowBasePtr = (unsigned char*)mxGetData(plhs[0]);

  //------------------------------------------------------------  
  // Now get the configurations
  //------------------------------------------------------------  


  CarmaConfig cc;
  std::vector<CarmaConfig::PadLocation> array1;
  std::vector<CarmaConfig::PadLocation> array2;
  
  for(unsigned i=0; i < nArray1; i++) {

    unsigned padNo   = (unsigned)arrayPads1Ptr[i];
    unsigned antType = (unsigned)arrayAntTypes1Ptr[i];

    if(antType == 0) {
      antType = CarmaConfig::GCP;
    } else if(antType == 1) {
      antType = CarmaConfig::BIMA;
    } else if(antType == 2) {
      antType = CarmaConfig::OVRO;
    }

    COUT("Adding pad: " << antType);
    cc.addPad(padNo, antType);
  }

  array1 = cc.getCurrentConfiguration();
  cc.initializeCurrentConfiguration();

  for(unsigned i=0; i < nArray2; i++) {

    unsigned padNo   = (unsigned)arrayPads2Ptr[i];
    unsigned antType = (unsigned)arrayAntTypes2Ptr[i];

    if(antType == 0) {
      antType = CarmaConfig::GCP;
    } else if(antType == 1) {
      antType = CarmaConfig::BIMA;
    } else if(antType == 2) {
      antType = CarmaConfig::OVRO;
    }

    COUT("Adding pad: " << antType);
    cc.addPad(padNo, antType);
  }

  array2 = cc.getCurrentConfiguration();

#if 1
  COUT("Array1 configuration is: " << std::endl);

  for(unsigned i=0; i < array1.size(); i++) {
    COUT(array1[i]);
  }

  COUT("");

  COUT("Array2 configuration is: " << std::endl);

  for(unsigned i=0; i < array2.size(); i++) {
    COUT(array2[i]);
  }

  COUT("");
#endif

  //------------------------------------------------------------  
  // Finally, loop through observations, checking shadowing
  //------------------------------------------------------------  

  Angle az, el;

  // Iterate over array1 antennas

  bool shadowed;
  for(unsigned iArr1=0; iArr1 < nArray1; iArr1++) {

    double* azPtr = azPtrBase + nObs * iArr1;
    double* elPtr = elPtrBase + nObs * iArr1;
    unsigned char* shadowPtr = shadowBasePtr + nObs * iArr1;
 
    CarmaConfig::PadLocation& arr1Ant = array1[iArr1];
    
    // Iterate over observations for this antenna

    for(unsigned iObs=0; iObs < nObs; iObs++) {

      if(iArr1 == 0 && iObs < 10) {
	COUT("AZ = " << azPtr[iObs]);
      }

      az.setDegrees(azPtr[iObs]);
      el.setDegrees(elPtr[iObs]);

      // Iterate over CARMA pads for this antenna, this observation

      shadowed = false;
      for(unsigned iArr2=0; iArr2 < array2.size(); iArr2++) {

	if(arr1Ant.padNumber_ != array2[iArr2].padNumber_) {
	  if(iArr1 == 0 && iObs < 10) {
	    COUT("Checking shadowing by: " << array2[iArr2]);
	  }
	  
	  shadowed |= arr1Ant.isShadowed(az, el, array2[iArr2]);
	}
      }

      shadowPtr[iObs] = shadowed;
    }
  }

  return;
}
