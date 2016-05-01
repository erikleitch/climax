/**.......................................................................
 * MATLAB Mex file for calculating UVW, given an array of antenna
 * locations, a source declination and a list of HA
 *
 * Use like:
 *
 * d=gcpMatCalcUvw(latitude(rad), antLoc(UENxN), dec(rad), HA(hours))
 *
 */
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"
#include "gcp/array/code/share/slalib/slalib.h"
#include "gcp/translator/DelayEngineNormal.h"
#include "gcp/array/code/share/include/rtcnetcoms.h"

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
  
  // Get RA/DEC, assumed to be in radians

  Angle lat;

  if(MexParser::isString(prhs[0]))
    lat.setDegrees(MexParser::getString(prhs[0]));
  else
    lat.setRadians(*MexParser::getDoubleData(prhs[0]));

  double* uen   = MexParser::getDoubleData(prhs[1]);

  for(unsigned i=0; i < 8; i++) {
    COUT("Ant: " << i);
    for(unsigned j=0; j < 3; j++) {
      unsigned ind = i*3 + j;
      if(j==0)
	COUT("E = " << *(uen+ind));
      if(j==1)
	COUT("N = " << *(uen+ind));
      if(j==2)
	COUT("U = " << *(uen+ind));
    }
  }

  Declination dec;
  if(MexParser::isString(prhs[2]))
    dec.setDegrees(MexParser::getString(prhs[2]));
  else
    dec.setRadians(*MexParser::getDoubleData(prhs[2]));

  COUT("DEC = " << dec);

  HourAngle ha;
  if(MexParser::isString(prhs[3]))
    ha.setHours(MexParser::getString(prhs[3]));
  else
    ha.setHours(*MexParser::getDoubleData(prhs[3]));

  COUT("HA = " << ha);

  gcp::translator::DelayEngineNormal delayEngine;

  Angle lng;
  lng.setDegrees("-118:17:45.9");
  lat.setDegrees("37:13:57.5");
  double alt = 1208  ;

  delayEngine.setSite(lng, lat, alt);

  HourAngle ra;
  ra.setHours("12:00:00");
    
  delayEngine.trackSource(ra, dec);
  delayEngine.selectRx();

  for(unsigned i=0; i < 8; i++) {
    COUT("Ant: " << i);

    double e = uen[i*3 + 0];
    double n = uen[i*3 + 1];
    double u = uen[i*3 + 2];

    AntNum ant(i);
    delayEngine.setAntennaLocation(ant.getId(), e, n, u);
  }

  std::vector<Delay> delays;
  delays.resize(8);

  delayEngine.setUseDelay(gcp::array::GEOMETRIC, true);

  delayEngine.setRefAnt(AntNum::ANT2);

  //  delayEngine.calculateDelays(delays, ha, dec, AntNum::ANTALL);

  // Create output lat/longitude arrays, in radians
  
  double* dl = MexHandler::createDoubleArray(&plhs[0], 8);

  for(unsigned i=0; i < 8; i++)
    dl[i] = delays[i].nanoSeconds();

  return;
}
