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
#include "gcp/util/Debug.h"
#include "gcp/util/DecAngle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/FitsIo.h"
#include "gcp/util/Frequency.h"
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
    
  if(nrhs != 4) {
    ThrowError(usage.str());
  }

  MexParser regDataParser(prhs[0]);
  MexParser uvDataParser(prhs[1]);
  MexParser fileParser(prhs[2]);
  MexParser modeParser(prhs[3]);

  if(!regDataParser.isStruct()) {
    ThrowError(usage.str());
  }

  if(!uvDataParser.isStruct()) {
    ThrowError(usage.str());
  }

  if(!fileParser.isString()) {
    ThrowError(usage.str());
  }

  if(!modeParser.isString()) {
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
  double* windDirPtr    = 0;
  double* windSpeedPtr  = 0;
  double* pressurePtr   = 0;

  double* antAzElPtr    = 0;
  double* equatGeocPtr  = 0;  
  double* locationPtr   = 0;

  double* sitePtr       = 0;

  const mxArray* versionPtr     = 0;
  const mxArray* srcNamePtr     = 0;
  const mxArray* goodSrcNamePtr = 0;

  double* freqReqPtr = 0;
  double* visWideRePtr  = 0;
  double* visWideImPtr  = 0;
  double* rmsPtr        = 0;
  double* uvwPtr        = 0;
  double* xyzPtr        = 0;
  bool* visFlagPtr      = 0;

  const mxArray* freqPtr = 0;

  double* visSpecRePtr  = 0;
  double* visSpecImPtr  = 0;

  //------------------------------------------------------------
  // Get pointers to individual fields in the parent structure
  //------------------------------------------------------------

  if(regDataParser.fieldExists("array.frame.features")) {
    featPtr = regDataParser.getFieldAsUint("array.frame.features");
  }

  if(regDataParser.fieldExists("array.frame.nsnap")) {
    combPtr = regDataParser.getFieldAsUint("array.frame.nsnap");
  }

  if(regDataParser.fieldExists("array.frame.utc")) {
    mjdPtr = regDataParser.getFieldAsDouble("array.frame.utc");
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

  if(regDataParser.fieldExists("array.weather.windSpeed")) {
    windSpeedPtr = regDataParser.getFieldAsDouble("array.weather.windSpeed");
  }

  if(regDataParser.fieldExists("array.weather.windDirection")) {
    windDirPtr = regDataParser.getFieldAsDouble("array.weather.windDirection");
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

  if(regDataParser.fieldExists("visPreAvg")) {
    visSpecRePtr = regDataParser.getFieldAsDouble("visPreAvg");
    visSpecImPtr = regDataParser.getImagFieldAsDouble("visPreAvg");
  }

  if(regDataParser.fieldExists("flags.visFlag")) {
    visFlagPtr = regDataParser.getFieldAsLogical("flags.visFlag");
  }

  //------------------------------------------------------------
  // Extract the uv data fields from the second argument
  //------------------------------------------------------------

  if(uvDataParser.fieldExists("vis")) {
    visWideRePtr = uvDataParser.getFieldAsDouble("vis");
    visWideImPtr = uvDataParser.getImagFieldAsDouble("vis");


#if 1
    for(unsigned i=1; i < 10; i++) {
      COUT("Re[" << i << "] = " << *(visWideRePtr + i)
	   << " Im[" << i << "] = " << *(visWideImPtr + i));
    }
#endif
  }

  if(uvDataParser.fieldExists("rms")) {
    rmsPtr = uvDataParser.getFieldAsDouble("rms");
  }

  if(uvDataParser.fieldExists("uvw")) {
    uvwPtr = uvDataParser.getFieldAsDouble("uvw");
  }

  if(uvDataParser.fieldExists("xyzMiriad")) {
    xyzPtr = uvDataParser.getFieldAsDouble("xyzMiriad");
  }

  if(uvDataParser.fieldExists("freq")) {
    freqPtr = uvDataParser.getField("freq");
  }

  //------------------------------------------------------------
  // Version
  //------------------------------------------------------------

  if(regDataParser.fieldExists("version")) {
    versionPtr = regDataParser.getField("version");
  }

  if(regDataParser.fieldExists("freq")) {
    freqReqPtr = regDataParser.getFieldAsDouble("freq");
  }

  //------------------------------------------------------------
  // Source name handling requires a bit more care.  This is a 
  // 'cell array' in matlab and must be handled differently from
  // other registers.
  //------------------------------------------------------------

  if(regDataParser.fieldExists("antenna.tracker.source")) {
    srcNamePtr = regDataParser.getField("antenna.tracker.source");
  }

  if(regDataParser.fieldExists("source")) {
    goodSrcNamePtr = regDataParser.getField("source");
  }

  //------------------------------------------------------------
  // Iterate over frames
  //------------------------------------------------------------

  unsigned nFrame = MexParser::getDimension(srcNamePtr, 0);
  unsigned nAnt   = MexParser::getDimension(srcNamePtr, 2);

  // Open the file with the requested mode

  MiriadIo io;
  io.openFile(fileParser.getString(), modeParser.getString());

  io.setVersion(MexParser::getString(versionPtr));
  io.setTelescopeName("GCP");
  io.setInstrument("Ka-Band Receivers");

  // Default to 30 GHz setup (conjugate phases)

  io.conjugateBaselines(true);

  // But if the requested frequency was 90, don't conjugate phases

  if(freqReqPtr && *freqReqPtr > 30) {
    io.conjugateBaselines(false);
  }

  // Set the vis, uvw, rms data up front, since these only need to be
  // set once

  io.setVisWide(visWideRePtr, visWideImPtr);
  io.setVisFlags(visFlagPtr);
  io.setRms(rmsPtr);
  io.setUvw(uvwPtr);
  io.setMjd(mjdPtr);
  io.setLst(lstPtr);

  // Set up IF frequency information

  setFrequencyInfo(io, freqPtr, visSpecRePtr, visSpecImPtr);

  // Set the observatory site once -- we don't expect this to change

  setArrayInfo(io, nFrame, nAnt, sitePtr);

  for(unsigned iFrame=0; iFrame < nFrame; iFrame++) {

    //------------------------------------------------------------
    // Set the integration time for this frame
    //------------------------------------------------------------

    setIntTime(io, iFrame, combPtr);

    //------------------------------------------------------------
    // Checking for noise source -- we need to create a fake source
    // called 'NOISE' so that CARMA pipeline doesn't choke
    //------------------------------------------------------------

    setSourceInfo(io, iFrame, nFrame, nAnt, featPtr, goodSrcNamePtr, equatGeocPtr, mjdPtr);

    //------------------------------------------------------------
    // Set up antenna locations for this frame
    //------------------------------------------------------------

    setAntennaInfo(io, iFrame, nFrame, nAnt, xyzPtr);

    //------------------------------------------------------------
    // Iterate over antennas to set the current az/el
    //------------------------------------------------------------

    setPointingInfo(io, iFrame, nFrame, nAnt, antAzElPtr);

    //------------------------------------------------------------
    // Set weather information
    //------------------------------------------------------------

    setWeatherInfo(io, iFrame, airTempPtr, 
		   windSpeedPtr, windDirPtr, 
		   relHumidPtr, pressurePtr);

    //------------------------------------------------------------
    // Finally, write this frame of data
    //------------------------------------------------------------

    io.writeCarmaFormatData(iFrame);
  }

  // Finally, close the file

  io.closeFile();
  io.reportVisStats();
}

/**.......................................................................
 * Set frequency information
 */
void setIntTime(MiriadIo& io, unsigned iFrame, unsigned int* combPtr)
{
  if(combPtr==0)
    return;

  Time intTime;
  intTime.setSeconds(combPtr[iFrame] * 0.5);
  io.setIntTime(intTime);
}

void setChannelInfo(MiriadIo& io, const mxArray* freqPtr)
{
  io.setNumberOfChannelsPerIf(17);
}

/**.......................................................................
 * Set frequency information
 */
void setFrequencyInfo(MiriadIo& io, const mxArray* freqPtr, double* visAvgRePtr, double* visAvgImPtr)
{
  if(freqPtr==0)
    return;

  unsigned nFreq = MexParser::getDimension(freqPtr, 1);
  double* freqData = MexParser::getDoubleData(freqPtr);

  std::vector<Frequency> freqs(nFreq);
  std::vector<Frequency> deltaFreqs(nFreq);
  
  for(unsigned i=0; i < nFreq; i++) {
    freqs[i].setHz(freqData[i]);
    deltaFreqs[i].setGHz(0.5);
  }
  
  io.setIfFrequencies(freqs);
  io.setDeltaIfFrequencies(deltaFreqs);

  // If we have pre-averaged visibilities, write them now

  if(visAvgRePtr && visAvgImPtr) {
    io.setNumberOfChannelsPerIf(17);
    Frequency freq;

    freq.setMHz(500.0 / 16);
    io.setDeltaChannelFrequency(freq);

    io.setVisSpec(visAvgRePtr, visAvgImPtr);
  }
}

/**.......................................................................
 * Set site information
 */
void setArrayInfo(MiriadIo& io, unsigned nFrame, unsigned nAnt, double* sitePtr)
{
  if(sitePtr == 0)
    return;

  Angle lng,lat;
  lng.setDegrees(*(sitePtr));
  lat.setDegrees(*(sitePtr+nFrame));

  io.setLatitude(lat);
  io.setLongitude(lng);
  io.setNumberOfTelescopes(nAnt);
  io.setNumberOfBaselines(28);
  io.setNumberOfFrames(nFrame);
}

void setAntennaInfo(MiriadIo& io, unsigned iFrame, unsigned nFrame, unsigned nAnt, double* xyzPtr)
{
  // Nnumber of telescopes

  io.setNumberOfTelescopes(nAnt);

  // Locations

  std::vector<std::vector<Length> > locs;

  locs.resize(nAnt);

  for(unsigned iAnt=0; iAnt < nAnt; iAnt++) {
    locs[iAnt].resize(3);
    for(unsigned iLoc=0; iLoc < 3; iLoc++) {
      unsigned ind = MexHandler::getMatlabIndex(nFrame, nAnt, 3, iFrame, iAnt, iLoc);
      locs[iAnt][iLoc].setMeters(xyzPtr[ind]);
    }
  }
  
  io.setTelescopeLocations(locs);

  // Fixed parameters

  Length diam;
  diam.setMeters(3.5);
  io.setTelescopeDiameter(diam);

  // A round number

  io.setTelescopeApertureEfficiency(0.6);
}

/**.......................................................................
 * Set source information for this frame
 */
void setSourceInfo(MiriadIo& io, unsigned iFrame, unsigned nFrame, unsigned nAnt,
		   unsigned int* featPtr, 
		   const mxArray* srcNamePtr, 
		   double* equatGeocPtr,
		   double* mjdPtr)
{
  static std::string lastSrc, currSrc;
  static HourAngle raApp;
  static DecAngle decApp;
  unsigned iRa, iDec;
  unsigned noiseBit  = 1<<10;
  unsigned targetBit = 1<<0;
  unsigned atmCalBit = 1<<1;
  unsigned ampCalBit = 1<<2;
  unsigned absCalBit = 1<<3;
  unsigned bpCalBit  = 1<<4;
  unsigned featureMask;

  // Get the current source

  currSrc = MexParser::getString(mxGetCell(srcNamePtr, iFrame));
  
  // Check the feature bit for this frame

  featureMask = *(featPtr+iFrame);
  if(featureMask & noiseBit) {
    io.setSourceName("NOISE");
  } else {
    io.setSourceName(currSrc);
  }
  
  // Set the apparent RA/Dec for this MJD
 
  iRa  = MexHandler::getMatlabIndex(nFrame, 3, nAnt, iFrame, 0, 0);
  iDec = MexHandler::getMatlabIndex(nFrame, 3, nAnt, iFrame, 1, 0);
  
  raApp.setHours(*(equatGeocPtr+iRa));
  decApp.setDegrees(*(equatGeocPtr+iDec));
  
  io.setRaApp(raApp);
  io.setDecApp(decApp);

  // Set te delta RA/DEC as well

  HourAngle dra;
  dra.setHours(0.0);

  DecAngle ddec;
  ddec.setDegrees(0.0);

  io.setDRaApp(dra);
  io.setDDecApp(ddec);
  
  // Need to convert from epoch of observation to standard epoch
  //  
  // If the current source is not the same as the last one, update the
  // J2000 RA/DEC info now.  This doesn't change, however, so only
  // update it when the source changes, since this is an expensive
  // calculation.
  
  if(currSrc != lastSrc) {

    HourAngle raMean;
    DecAngle decMean;

    TimeVal time;
    time.setMjd(*(mjdPtr+iFrame));

    gcp::util::Astrometry::apparentToJ2000Place(raApp, decApp, time,
					       raMean, decMean);

    // Need to convert from epoch of observation to standard epoch
  
    io.setRa(raMean);
    io.setDec(decMean);
  }

  // Check the feature mask for notable bits

  if(featureMask & atmCalBit) {
    io.setPurpose('A');
  } else if(featureMask & bpCalBit) {
    io.setPurpose('B');
  } else if(featureMask & targetBit) {
    io.setPurpose('S');
  } else if(featureMask & ampCalBit) {
    io.setPurpose('G');
  } else if(featureMask & absCalBit) {
    io.setPurpose('F');
  } else {
    io.setPurpose('O');
  }

  lastSrc = currSrc;
}

void setPointingInfo(MiriadIo& io, unsigned iFrame, unsigned nFrame, 
		     unsigned nAnt, 
		     double* antAzElPtr)
{
  static unsigned antMax = 8;
  static std::vector<Angle> az(antMax);
  static std::vector<Angle> el(antMax);
  unsigned iAz, iEl;

  if(nAnt != antMax) {
    az.resize(nAnt);
    el.resize(nAnt);
  }

  for(unsigned iAnt=0; iAnt < nAnt; iAnt++) {
    iAz = MexHandler::getMatlabIndex(nFrame, 3, nAnt, iFrame, 0, iAnt);
    iEl = MexHandler::getMatlabIndex(nFrame, 3, nAnt, iFrame, 1, iAnt);
    
    az[iAnt].setDegrees(*(antAzElPtr+iAz));
    el[iAnt].setDegrees(*(antAzElPtr+iEl));
  }

  io.setTelescopeAzimuth(az);
  io.setTelescopeElevation(el);
}

void setWeatherInfo(MiriadIo& io, unsigned iFrame, double* airTempPtr, 
		    double* windSpeedPtr, double* windDirPtr, 
		    double* relHumidPtr, double* pressurePtr)
{
  static Temperature airTemp;
  static Angle windDir;
  static Speed windSpeed;
  static Pressure pressure;
  static double relHumid;

  if(airTempPtr) {
    airTemp.setC(*(airTempPtr+iFrame));
    io.setAirTemperature(airTemp);
  }

  if(windSpeedPtr) {
    windSpeed.setMetersPerSec(*(windSpeedPtr+iFrame));
    io.setWindSpeed(windSpeed);
  }

  if(windDirPtr) {
    windDir.setDegrees(*(windDirPtr+iFrame));
    io.setWindDirection(windDir);
  }

  if(pressurePtr) {
    pressure.setMilliBar(*(pressurePtr+iFrame));
    io.setPressure(pressure);
  }

  if(relHumidPtr) {
    io.setRelativeHumidity(*(relHumidPtr+iFrame)/100);
  }
}

