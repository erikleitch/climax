/**.......................................................................
 * MATLAB Mex file for writing out Miriad files
 *
 * Use like:
 *
 *    gcpMatMiriad2(d);
 *
 */
#include <iostream>
#include <cmath>

#include "gcp/util/Astrometry.h"
#include "gcp/util/Coordinates.h"
#include "gcp/util/Debug.h"
#include "gcp/util/DecAngle.h"
#include "gcp/util/DelayAntennaLocation.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/FitsIo.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/Location.h"
#include "gcp/util/MiriadIo.h"

#include "gcp/matlab/SzaUvfWriter.h"
#include "gcp/matlab/MexHandler.h"

#include "mex.h"
#include "matrix.h"

using namespace gcp::util;
using namespace gcp::matlab;

void writeMiriadFile(int nrhs, const mxArray* prhs[]);
void setFrequencyInfo(MiriadIo& io, const mxArray* freqPtr, 
		      double* visSpecRePtr, double* visSpecImPtr);
void setArrayInfo(MiriadIo& io, unsigned nFrame, unsigned nAnt, 
		  double* sitePtr);
void setAntennaInfo(MiriadIo& io, unsigned iFrame, unsigned nFrame, unsigned nAnt, double* xyzPtr);
void setIntTime(MiriadIo& io, unsigned iFrame, unsigned int* combPtr);
void setSiteInfo(MiriadIo& io, unsigned nFrame, double* sitePtr);
void setSourceInfo(MiriadIo& io, 
		   unsigned iFrame, unsigned nFrame, unsigned nAnt,
		   unsigned int* featPtr, const mxArray* srcNamePtr, 
		   double* equatGeocPtr,
		   double* mjdPtr);
void setPointingInfo(MiriadIo& io, unsigned iFrame, unsigned nFrame, 
		     unsigned nAnt, double* antAzElPtr);
void setWeatherInfo(MiriadIo& io, unsigned iFrame, double* airTempPtr, 
		    double* windSpeedPtr, double* windDirPtr, 
		    double* relHumidPtr, double* pressurePtr);

void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  try {  

    gcp::util::Logger::installStdoutPrintFn(MexHandler::stdoutPrintFn);
    gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);

    gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
    gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
    gcp::util::ErrHandler::installLogFn(MexHandler::logFn);

    writeMiriadFile(nrhs, prhs);

  } catch(gcp::util::Exception& err) {
    mexErrMsgTxt(err.what());
  } catch(...) {
    mexErrMsgTxt("Caught an unknown error");
  }
}

/**.......................................................................
 * Main function, to write a Miriad file
 */
void writeMiriadFile(int nrhs, const mxArray* prhs[])
{
  //------------------------------------------------------------
  // Sanity check arguments:
  //------------------------------------------------------------

  std::ostringstream usage;

  usage << "Usage: " << std::endl
	<< std::endl
	<< "  gcpMatMiriad2(d, m, fileName, appendMode)" << std::endl
	<< std::endl << std::endl
	<< "Where: " << std::endl
	<< std::endl
	<< "  d = struct returned by read_arc()" << std::endl
	<< "  m = struct returned by gcp2uv()" << std::endl
	<< "  fileName = the name of the miriad 'file' to create" << std::endl
	<< "  appendMode = either 'append' to append, or 'new' to create from scratch" << std::endl;
    
  MexParser regDataParser(prhs[0]);

  if(!regDataParser.isStruct()) {
    ThrowError(usage.str());
  }

  //------------------------------------------------------------
  // These will point to individual fields in the structs
  //------------------------------------------------------------

  unsigned int* featPtr = 0;
  unsigned int* combPtr = 0;
  double* mjdPtr        = 0;
  double* lstPtr        = 0;

  double* airTempPtr    = 0;
  double* relHumidPtr   = 0;
  double* pressurePtr   = 0;

  double* antAzElPtr    = 0;
  double* equatGeocPtr  = 0;  
  double* locationPtr   = 0;

  double* sitePtr       = 0;

  const mxArray* versionPtr = 0;
  const mxArray* srcNamePtr = 0;

  double* visWideRePtr  = 0;
  double* visWideImPtr  = 0;
  double* rmsPtr        = 0;
  double* uvwPtr        = 0;
  double* xyzPtr        = 0;

  double* visSpecRePtr  = 0;
  double* visSpecImPtr  = 0;

  //------------------------------------------------------------
  // Get pointers to individual fields in the parent structure
  //------------------------------------------------------------

  COUT("Here 0");
  if(regDataParser.fieldExists("array.frame.features")) {
    featPtr = regDataParser.getFieldAsUint("array.frame.features");
  }

  if(regDataParser.fieldExists("array.frame.nsnap")) {
    combPtr = regDataParser.getFieldAsUint("array.frame.nsnap");
  }

  if(regDataParser.fieldExists("antenna.tracker.lst")) {
    lstPtr = regDataParser.getFieldAsDouble("antenna.tracker.lst");
  }

  if(regDataParser.fieldExists("array.weather.airTemperature")) {
    airTempPtr = regDataParser.getFieldAsDouble("array.weather.airTemperature");
  }

  if(regDataParser.fieldExists("array.weather.relativeHumidity")) {
    relHumidPtr = regDataParser.getFieldAsDouble("array.weather.relativeHumidity");
  }

  if(regDataParser.fieldExists("array.weather.pressure")) {
    pressurePtr = regDataParser.getFieldAsDouble("array.weather.pressure");
  }

  if(regDataParser.fieldExists("antenna.tracker.actual")) {
    antAzElPtr = regDataParser.getFieldAsDouble("antenna.tracker.actual");
  }

  if(regDataParser.fieldExists("antenna.tracker.equat_geoc")) {
    equatGeocPtr = regDataParser.getFieldAsDouble("antenna.tracker.equat_geoc");
  }

  if(regDataParser.fieldExists("array.delay.location")) {
    locationPtr = regDataParser.getFieldAsDouble("array.delay.location");
  }

  if(regDataParser.fieldExists("array.delay.siteFiducial")) {
    sitePtr = regDataParser.getFieldAsDouble("array.delay.siteFiducial");
  }

  double* freqPtr = 0;
  if(regDataParser.fieldExists("array.realRF")) {
    freqPtr = regDataParser.getFieldAsDouble("array.realRF");
  }

  double* visRePtr=0;
  double* visImPtr=0;
  if(regDataParser.fieldExists("vis")) {
    visRePtr = regDataParser.getFieldAsDouble("vis");
    visImPtr = regDataParser.getImagFieldAsDouble("vis");
  }

  const mxArray* mxPtr = regDataParser.getField("array.frame.utc");
  unsigned nFrame = MexParser::getDimension(mxPtr, 0);
  mjdPtr = regDataParser.getFieldAsDouble("array.frame.utc");

  COUT("Here 1");

  //-----------------------------------------------------------------------
  // Create output members
  //-----------------------------------------------------------------------

  const mxArray* trackerPtr = regDataParser.getField("antenna.tracker");
  double* azelOutPtr  = MexHandler::copyNamedDoubleStructField((mxArray*)trackerPtr, "azel", "actual");
  double* delayOutPtr = MexHandler::copyNamedDoubleStructField((mxArray*)trackerPtr, "tropo", "lst");
  double* phaseOutPtr = MexHandler::copyNamedDoubleStructField((mxArray*)prhs[0], "tropoPhase", "vis");

  COUT("Here 2");

  //-----------------------------------------------------------------------
  // Now iterate through the data, calculating tropospheric delays
  //-----------------------------------------------------------------------

  gcp::util::DelayAntennaLocation antDLoc, refDLoc;

  COUT("Here 2a1");

  antDLoc.useGeometricDelay(true);
  COUT("Here 2a2");
  antDLoc.useAdjustableDelay(false);
  COUT("Here 2a3");
  antDLoc.useFixedDelay(false);
  COUT("Here 2a4");
  antDLoc.useTroposphericDelay(true);

  COUT("Here 2b");

  Angle lng, lat;
  double alt;
  Pressure press;
  Temperature airTemp;
  HourAngle ha;
  DecAngle dec;

  COUT("Here 3");
  unsigned uInd, eInd, nInd, raInd, decInd, lstInd, azInd, elInd, delayInd;

  // Store the fiducial site parameters -- these won't change 

  lng.setDegrees(sitePtr[MexHandler::getMatlabIndex(nFrame, 3, 0, 0)]);
  lat.setDegrees(sitePtr[MexHandler::getMatlabIndex(nFrame, 3, 0, 1)]);
  alt =          sitePtr[MexHandler::getMatlabIndex(nFrame, 3, 0, 2)];
  
  // Make the reference location be pad 32 (ie, 0,0,0 at fiducial
  // location, which is the lat/lng of pad 32)

  refDLoc.setFiducialSite(lng, lat, alt);
  refDLoc.setOffset(0.0, 0.0, 0.0);

  // Set the fiducial site to be the same for the arbitrary location

  antDLoc.setFiducialSite(lng, lat, alt);
  antDLoc.atmos().setRx(Rx::RX30GHZ);

  gcp::util::Vector<Angle> azel;
  Delay delay;

  Source src;
  HourAngle raApp, raMean;
  DecAngle decApp, decMean;

  src.setType(gcp::array::SRC_J2000);

  COUT("Here 4");

  for(unsigned iFrame=0; iFrame < nFrame; iFrame++) {

    COUT("iFrame: "  << iFrame);

    // The atmospheric parameters are not antenna-specific

    press.setMilliBar(pressurePtr[iFrame]);
    airTemp.setC(airTempPtr[iFrame]);
    
    refDLoc.atmos().setPressure(pressurePtr[iFrame]);
    refDLoc.atmos().setHumidity(relHumidPtr[iFrame]/100);
    refDLoc.atmos().setAirTemperature(airTemp);
    refDLoc.atmos().setRx(Rx::RX30GHZ);

    // Now iterate over antennas
    
    for(unsigned iAnt=0; iAnt < 8; iAnt++) {

      // Get indices into relevant matlab arrays

      uInd     = MexHandler::getMatlabIndex(nFrame, 8, 3, iFrame, iAnt, 0);
      eInd     = MexHandler::getMatlabIndex(nFrame, 8, 3, iFrame, iAnt, 1);
      nInd     = MexHandler::getMatlabIndex(nFrame, 8, 3, iFrame, iAnt, 2);

      raInd    = MexHandler::getMatlabIndex(nFrame, 3, 8, iFrame, 0, iAnt);
      decInd   = MexHandler::getMatlabIndex(nFrame, 3, 8, iFrame, 1, iAnt);
      lstInd   = MexHandler::getMatlabIndex(nFrame, 1, 8, iFrame, 0, iAnt);
      
      azInd    = MexHandler::getMatlabIndex(nFrame, 3, 8, iFrame, 0, iAnt);
      elInd    = MexHandler::getMatlabIndex(nFrame, 3, 8, iFrame, 1, iAnt);
      delayInd = MexHandler::getMatlabIndex(nFrame, 1, 8, iFrame, 0, iAnt);

      // Set the location of this antenna, relative to the fiducial
      // site we already specified

      antDLoc.setOffset(locationPtr[uInd], locationPtr[eInd], locationPtr[nInd]);

      // Set the atmospheric parameters we need

      antDLoc.atmos().setPressure(pressurePtr[iFrame]);
      antDLoc.atmos().setHumidity(relHumidPtr[iFrame]/100);
      antDLoc.atmos().setAirTemperature(airTemp);

      // Just as a check, calculate the az/el corresponding to this location

      double mjd = mjdPtr[iFrame];
      raApp.setHours(equatGeocPtr[raInd]);
      decApp.setDegrees(equatGeocPtr[decInd]);

      src.extend(mjd, raApp, decApp, 0.0);

      ha.setHours(lstPtr[lstInd] - equatGeocPtr[raInd]);
      azel = Coordinates::laAndHaDecToAzEl(antDLoc.latitude(false), antDLoc.altitude(false), ha, decApp, false);

      azelOutPtr[azInd] = azel[0].degrees();

      if(azelOutPtr[azInd] < 0.0)
	azelOutPtr[azInd] += 360;

      // Store the elevation, including refraction correction

      azelOutPtr[elInd] = azel[1].degrees() + antDLoc.atmos().offset(azel[1]).degrees();
      
      // Calculate the tropospheric delay, in meters

      delay = antDLoc.totalDelay(mjd, &src, &refDLoc, false);
      delayOutPtr[delayInd] = delay.meters();
    }
  }

  COUT("About to calculate phase corrections");

  // Now construct the baseline-based residual phases due to the
  // troposphere

  unsigned visInd, ant1Ind, ant2Ind, iBase;
  Wavelength wave;
  Frequency freq;

  for(unsigned iFrame=0; iFrame < nFrame; iFrame++) {

    iBase = 0;
    for(unsigned iAnt1=0; iAnt1 < 7; iAnt1++) {

      ant1Ind = MexHandler::getMatlabIndex(nFrame, 1, 8, iFrame, 0, iAnt1);
      for(unsigned iAnt2=iAnt1+1; iAnt2 < 8; iAnt2++) {

	ant2Ind = MexHandler::getMatlabIndex(nFrame, 1, 8, iFrame, 0, iAnt2);
 	for(unsigned iFreq=0; iFreq < 16; iFreq++) {

	  visInd  = MexHandler::getMatlabIndex(nFrame, 28, 16, iFrame, iBase, iFreq);

	  // Get the baseline phase, in radians corresponding to these antenna delays

	  freq.setGHz(freqPtr[iFreq]);
	  wave.setFrequency(freq);

	  double dPh = 2*M_PI*(delayOutPtr[ant1Ind] - delayOutPtr[ant2Ind]) / wave.meters();

	  phaseOutPtr[visInd] = dPh / M_PI * 180;

	  double re  = visRePtr[visInd];
	  double im  = visImPtr[visInd];

	  double cPh = cos(dPh);
	  double sPh = sin(dPh);

	  //------------------------------------------------------------
	  // And correct the visibility for this phase:
	  //
	  // C' = C * exp(-i Ph) = (R + i*I) * (cPh - i*sPh)
	  //                     = (R*cPh + I*sPh) + i*(I*cPh - R*sPh)
	  //------------------------------------------------------------

	  visRePtr[visInd] = re * cPh + im * sPh;
	  visImPtr[visInd] = im * cPh - re * sPh;
	  //	  visRePtr[visInd] = 1;
	  //	  visImPtr[visInd] = 2;

	  if(iAnt1 == 3 && iAnt2 == 4 && iFreq==6) {
	    COUT("dPh = " << dPh);
	    COUT("cPh = " << cPh);
	    COUT("sPh = " << sPh);
	    COUT("re = " << re);
	    COUT("im = " << im);
	    COUT("re' = " << visRePtr[visInd]);
	    COUT("im' = " << visImPtr[visInd]);
	  }


	}

	++iBase;
      }

    }
  }

  return;
}

