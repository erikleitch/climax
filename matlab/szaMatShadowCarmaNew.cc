/**.......................................................................
 * MATLAB Mex file for calculating shadowing, given a specified GCP
 * and CARMA array
 *
 * Use like:
 *
 * shadowed(bool[Nobs x Nant]) = gcpMatShadowCarma(az(double[Nobs x Nant]), 
 *                                                 el(double[Nobs x Nant]), 
 *                                                 d);
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

void setUpArray(CarmaConfig& cc, std::vector<CarmaConfig::PadLocation>& carma, 
		ArrayConfig::Type& currArray, ArrayConfig::Type& lastArray, 
		unsigned iObs, unsigned nObs, unsigned nAnt, 
		unsigned* carmaPadNumberPtr, unsigned* carmaAntTypePtr=0);

/**.......................................................................
 * Entry point from the matlab environment
 */
void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  // Redirect error handling

  gcp::util::Logger::installStdoutPrintFn(MexHandler::stdoutPrintFn);
  gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);
  
  gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
  gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
  gcp::util::ErrHandler::installLogFn(MexHandler::logFn);

  // Create Usage message to be displayed on input error

  std::ostringstream usage;
  
  usage << "Usage: " 
	<< std::endl
	<< std::endl
	<< "   shadowed = gcpMatShadowCarma(double az[Nobs x 8], " << std::endl
	<< "                                double el[Nobs x 8], " << std::endl
	<< "                                d)"
	<< std::endl << std::endl
	<< "where d must contain fields: " 
	<< std::endl << std::endl
	<< "   array.carma.config"  << std::endl
	<< "   array.carma.nPad"    << std::endl
	<< "   array.carma.pad"     << std::endl
	<< "   array.carma.antType" << std::endl
	<< std::endl
	<< "   array.sza.config"  << std::endl
	<< "   array.sza.pad"     << std::endl
	<< std::endl
	<< "Note: if these fields are not present in your data structure, " << std::endl
	<< "  it may mean that you are reducing data taken before the CARMA " << std::endl
	<< "  array configuration was stored in the archive.  " << std::endl
	<< std::endl
	<< "  If so, you should call gcpMatCarmaConfiguration(d) before calling" << std::endl
	<< "  this function to fill in the required information" << std::endl << std::endl;

  //------------------------------------------------------------
  // Input argument parsing
  //------------------------------------------------------------

  if(nrhs != 3) {
    ThrowError("Not enough input arguments: " << std::endl << usage.str());
  }

  // Check that the az & el dimensions match

  if(!MexParser::dimensionsMatch(prhs[0], prhs[1])) {
    ThrowError("Dimensions of the az & el arrays must match:" << std::endl << usage.str());
  }

  // Check that the az/el arrays have 8 telescopes as their last dimension
  
  MexParser azParser(prhs[0]);
  MexParser elParser(prhs[1]);

  unsigned nDim = azParser.getNumberOfDimensions();
  if(azParser.getDimension(nDim-1) != 8) {
    ThrowError("The az & el arrays must be [Nobs x 8] for the GCP:" << std::endl << usage.str());
  }

  double* azPtrBase = azParser.getDoubleData();
  double* elPtrBase = elParser.getDoubleData();

  unsigned nObs = azParser.getNumberOfElements()/8;
  unsigned nSza = 8;

  // Check that required members are present

  MexParser base(prhs[2]);

  if(!base.fieldExists("array.carma.config") || 
     !base.fieldExists("array.carma.nPad") || 
     !base.fieldExists("array.carma.pad") || 
     !base.fieldExists("array.carma.antType")) {
    ThrowError(usage.str());
  }

  // Check that the carma configuration array makes sense

  MexParser carmaConfigParser(base.getField("array.carma.config"));

  if(!carmaConfigParser.isUchar()) {
    ThrowError("The configuration array must an array of unsigned chars: " << usage.str());
  } 

  if(carmaConfigParser.getNumberOfElements() != nObs) {
    ThrowError("The configuration array must also be " << nObs << " long: " << usage.str());
  }

  unsigned char* carmaConfigPtr = carmaConfigParser.getUcharData();

  // Check that the gcp configuration array makes sense

  MexParser gcpConfigParser(base.getField("array.gcp.config"));

  if(!gcpConfigParser.isUchar()) {
    ThrowError("The configuration array must be an array of unsigned chars: " << usage.str());
  } 

  if(gcpConfigParser.getNumberOfElements() != nObs) {
    ThrowError("The configuration array must also be " << nObs << " long: " << usage.str());
  }

  unsigned char* gcpConfigPtr = gcpConfigParser.getUcharData();

  // Check that the number of antennas makes sense

  MexParser carmaNantParser(base.getField("array.carma.nPad"));

  if(!carmaNantParser.isUint()) {
    ThrowError("The array of CARMA pads must be an array of unsigned ints: " << usage.str());
  } 

  if(carmaNantParser.getNumberOfElements() != nObs) {
    ThrowError("The array of CARMA pads must also be " << nObs <<  " long: " << usage.str());
  }

  unsigned int* carmaNantPtr = carmaNantParser.getUintData();

  // Check that the CARMA pad numbers make sense

  MexParser carmaPadNumberParser(base.getField("array.carma.pad"));

  if(!carmaPadNumberParser.isUint()) {
    ThrowError("The pad number array must be an array of unsigned ints: " << usage.str());
  } 

  if(carmaPadNumberParser.getNumberOfDimensions() != 2) {
    ThrowError("The pad number array must be [" << nObs << " x nMaxAnt] long: " << usage.str());
  } 
  
  unsigned nMaxAnt = carmaPadNumberParser.getDimension(1);
  unsigned int* carmaPadNumberPtr = carmaPadNumberParser.getUintData();

  // Check that the GCP pad numbers make sense

  MexParser gcpPadNumberParser(base.getField("array.gcp.pad"));

  if(!gcpPadNumberParser.isUint()) {
    ThrowError("The pad number array must be an array of unsigned ints: " << usage.str());
  } 

  if(gcpPadNumberParser.getNumberOfDimensions() != 2) {
    ThrowError("The pad number array must be [" << nObs << " x nMaxAnt] long: " << usage.str());
  } 
  
  unsigned int* gcpPadNumberPtr = gcpPadNumberParser.getUintData();

  // Check that the antenna type arrays make sense

  MexParser carmaAntTypeParser(base.getField("array.carma.antType"));

  if(!carmaAntTypeParser.isUint()) {
     ThrowError("The antenna type array must be an array of unsigned ints: " << usage.str());
  }

  if(carmaAntTypeParser.getNumberOfDimensions() != 2 || carmaAntTypeParser.getDimension(1) != nMaxAnt) {
    ThrowError("The antenna type array must be [" << nObs << " x " << nMaxAnt << "] long: " << usage.str());
  }

  unsigned int* carmaAntTypePtr = carmaAntTypeParser.getUintData();
  
  //------------------------------------------------------------
  // Create output vector to have the same dimensions as the input
  // az/el arrays
  //------------------------------------------------------------

  plhs[0] = MexHandler::createMatlabArray(nDim, azParser.getDimensions(), DataType::BOOL);
  unsigned char* shadowBasePtr = (unsigned char*)mxGetData(plhs[0]);

  //------------------------------------------------------------  
  // Now loop through observations, checking shadowing
  //------------------------------------------------------------  

  Angle az, el;
  bool shadowed;

  CarmaConfig cc;

  std::vector<CarmaConfig::PadLocation> carma;
  std::vector<CarmaConfig::PadLocation> gcp;

  // Iterate over observations

  ArrayConfig::Type lastSzaArray = ArrayConfig::UNKNOWN;
  ArrayConfig::Type lastCarmaArray = ArrayConfig::UNKNOWN;
  ArrayConfig::Type currCarmaArray, currSzaArray;

  for(unsigned iObs=0; iObs < nObs; iObs++) {

    // Set up the CARMA array for this observation

    currCarmaArray = (ArrayConfig::Type)(*(carmaConfigPtr + iObs));
    setUpArray(cc, carma, currCarmaArray, lastCarmaArray, iObs, nObs, *(carmaNantPtr + iObs), 
	       carmaPadNumberPtr, carmaAntTypePtr);

    // Set up the GCP array for this observation

    currSzaArray = (ArrayConfig::Type)(*(gcpConfigPtr + iObs));
    setUpArray(cc, gcp, currSzaArray, lastSzaArray, iObs, nObs, 8,
	       gcpPadNumberPtr);

    // Iterate over GCP antennas

    for(unsigned iSza=0; iSza < nSza; iSza++) {

      double* azPtr = azPtrBase + nObs * iSza;
      double* elPtr = elPtrBase + nObs * iSza;
      unsigned char* shadowPtr = shadowBasePtr + nObs * iSza;
 
      CarmaConfig::PadLocation& gcpAnt = gcp[iSza];
    
      az.setDegrees(azPtr[iObs]);
      el.setDegrees(elPtr[iObs]);
      
      // Iterate over CARMA pads for this antenna, this observation

      shadowed = false;

      // Only check shadowing if there was a CARMA array to worry about

      if(currCarmaArray != ArrayConfig::NOCARMA) {

	for(unsigned iCarma=0; iCarma < carma.size(); iCarma++) {
	  shadowed |= gcpAnt.isShadowed(az, el, carma[iCarma]);
	}
      }

      // Store the result

      shadowPtr[iObs] = shadowed;
    }
  }

  return;
}

/**.......................................................................
 * Return array data for the current array configuration
 */
void setUpArray(CarmaConfig& cc, std::vector<CarmaConfig::PadLocation>& pads, 
		ArrayConfig::Type& currArray, ArrayConfig::Type& lastArray, 
		unsigned iObs, unsigned nObs, unsigned nAnt, 
		unsigned* padNumberPtr, unsigned* antTypePtr)
{
  // Only set up the array manually if the array configuration is
  // unknown

  if(currArray != ArrayConfig::UNKNOWN) {
    
    // If the current array is the same as the last array, do nothing.

    if(currArray == lastArray) {
      return;

      // Else set up for the new configuration

    } else if(currArray != ArrayConfig::NOCARMA) {
      cc.setCurrentConfiguration(currArray);
    } else {
      return;
    }

    // If the array is not recognized, set it up manually

  } else {
    cc.initializeCurrentConfiguration();

    unsigned padNo, antType;
    for(unsigned iAnt=0; iAnt < nAnt; iAnt++) {
      padNo   = *(padNumberPtr + iAnt * nObs + iObs);
      antType = antTypePtr ? *(antTypePtr   + iAnt * nObs + iObs) : CarmaConfig::GCP;
      cc.addPad(padNo, antType);
    }
  }

  // Store the new configuration in the pad array and update the
  // lastArray variable

  pads = cc.getCurrentConfiguration();
  lastArray = currArray;
}

