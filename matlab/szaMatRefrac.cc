/**.......................................................................
 * MATLAB Mex file for calculating refraction coefficients
 *
 */
#include <iostream.h>
#include <math.h>

#include "gcp/util/Debug.h"
#include "gcp/util/Atmosphere.h"

#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

#include "mex.h"
#include "matrix.h"

using namespace gcp::util;
using namespace gcp::matlab;

void calcRefraction(const mxArray* prhs[], mxArray* plhs[]);

void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  try {  

    if(nrhs != 7) {
      std::cerr << "Wrong number of arguments. " << std::endl
		<< "Should be: elevation(degrees) altitude(meters) latitude(dd:mm:ss) wavelength(cm) temperature[NxM](C) pressure[NxM](mBar) rel_humidity[NxM](%)" << std::endl;
      return;
    }

    gcp::util::Logger::installStdoutPrintFn(MexHandler::stdoutPrintFn);
    gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);

    gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
    gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
    gcp::util::ErrHandler::installLogFn(MexHandler::logFn);

    calcRefraction(prhs, plhs);

  } catch(gcp::util::Exception& err) {
    mexErrMsgTxt(err.what());
  } catch(...) {
    mexErrMsgTxt("Caught an unknown error");
  }
}

/**.......................................................................
 * Calculate the refraction coefficients
 */
void calcRefraction(const mxArray* prhs[], mxArray* plhs[])
{
  // Get the elevation

  Angle elevation(Angle::Degrees(), *MexParser::getDoubleData(prhs[0]));

  // Get the (fixed) altitude, latitude, and wavelength

  Length altitude(Length::Meters(), 
		  *MexParser::getDoubleData(prhs[1]));

  Angle latitude(MexParser::getString(prhs[2]));                        

  Wavelength radioWavelength(Length::Centimeters(), 
			     *MexParser::getDoubleData(prhs[3])); 

  Wavelength opticalWavelength = Atmosphere::opticalWavelength_;

  // Get symbolic pointers to the temperature, pressure and humidity arrays

  const mxArray* mxTemp  = prhs[4];
  const mxArray* mxPress = prhs[5];
  const mxArray* mxHumid = prhs[6];

  // Check that the array sizes match
  
  if(!MexParser::dimensionsMatch(mxTemp, mxPress) || 
     !MexParser::dimensionsMatch(mxTemp, mxHumid))
    ThrowError("Temperature/Pressure/Humidity array dimensions don't agree");

  // Get the double arrays corresponding to the input mxArrays

  double* tempC        = MexParser::getDoubleData(prhs[4]); // Temperature
  double* pressureMbar = MexParser::getDoubleData(prhs[5]); // Pressure
  double* humidPerc    = MexParser::getDoubleData(prhs[6]); // Humidity

  // Create return arrays

  int ndim        = MexParser::getNumberOfDimensions(mxTemp);
  const int* dims = MexParser::getDimensions(mxTemp);

  double* radOff = MexHandler::createDoubleArray(&plhs[0], ndim, dims);
  double* optOff = MexHandler::createDoubleArray(&plhs[1], ndim, dims);

  // Now loop over array elements, calculating refraction coeffs for each one

  Temperature temp;
  Atmosphere::RefractionCoefficients coeff;

  unsigned nel = MexParser::getNumberOfElements(mxTemp);

  // Finally, calculate over all elements of the input arrays

  for(unsigned iel=0; iel < nel; iel++) {

    temp.setC(tempC[iel]);
    coeff = Atmosphere::refractionCoefficients(altitude, temp, 
					       pressureMbar[iel], 
					       humidPerc[iel], 
					       radioWavelength, latitude);

    radOff[iel] = (Atmosphere::offset(elevation, coeff)).arcsec();

    coeff = Atmosphere::refractionCoefficients(altitude, temp, 
					       pressureMbar[iel], 
					       humidPerc[iel], 
					       opticalWavelength, latitude);

    optOff[iel] = (Atmosphere::offset(elevation, coeff)).arcsec();
  }
} 

