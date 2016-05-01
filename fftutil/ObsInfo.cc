#include "gcp/fftutil/FitsIoHandler.h"
#include "gcp/fftutil/ObsInfo.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/util/Astrometry.h"
#include "gcp/util/Constants.h"
#include "gcp/util/Geoid.h"
#include "gcp/util/ObsParameter.h"
#include "gcp/util/Sampler.h"
#include "gcp/util/String.h"
#include "gcp/util/Timer.h"

using namespace std;
using namespace gcp::util;

#define SIM_TIMER_TEST

#ifdef SIM_TIMER_TEST
  Timer t1, t2, t3;
  double t1time = 0.0;
  double t2time = 0.0;
  double t3time = 0.0;
#endif

/**.......................................................................
 * Constructor.
 */
ObsInfo::ObsInfo() 
{
  infoMask_   = ObsParameter::INFO_NONE;
  noiseType_  = NOISE_NONE;
  noiseRms_.setVal(1.0, "mJy");

  obsRa_.hasValue_      = false;
  obsDec_.hasValue_     = false;

  obsEquinox_ = 2000;

  visibilities_.resize(0);
}

/**.......................................................................
 * Destructor.
 */
ObsInfo::~ObsInfo() 
{
  visibilities_.resize(0);
}

ObsInfo::ObsInfo(const ObsInfo& obs)
{
  *this = obs;
}

ObsInfo::ObsInfo(ObsInfo& obs)
{
  *this = obs;
}

/**.......................................................................
 * Set the number of antennas in this observation
 */
void ObsInfo::setNumberOfAntennas(unsigned nAnt)
{
  antennas_.resize(nAnt);
  nBaseline_ = (nAnt * (nAnt-1))/2;

  for(unsigned iAnt=0; iAnt < antennas_.size(); iAnt++) {
    antennas_[iAnt].antNo_ = iAnt + 1;
  }

  // If the location has been set, pass the array location down to any
  // antennas that have been added so far

  if(infoMask_ & ObsParameter::INFO_LOCATION_ARRAY) {
    for(unsigned iAnt=0; iAnt < antennas_.size(); iAnt++) {
      antennas_[iAnt].setReferenceLla(lla_);
    }
  }

  // Mark this information as received

  infoMask_ |= ObsParameter::INFO_NANT;
}

/**.......................................................................
 * Get the number of antennas in this observation
 */
unsigned ObsInfo::getNumberOfAntennas()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_NANT);
  return antennas_.size();
}

/**.......................................................................
 * Get the number of baselines in this observation
 */
unsigned ObsInfo::getNumberOfBaselines()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_NANT);
  return nBaseline_;
}

/**.......................................................................
 * Specify the type of an antenna
 */
void ObsInfo::setAntennaType(Antenna::AntennaType type, int index)
{
  if(antennas_.size() == 0) {
    ThrowError("The number of antennas is currently unknown.  Use setNumberOfAntennas()");
  }

  // Negative index means default to all

  if(index < 0) {
    for(unsigned i=0; i < antennas_.size(); i++) {
      antennas_[i].setType(type);
    }
  } else if(index > antennas_.size()-1) {
    ThrowError("Invalid index: " << index << " should be: 0 - " << antennas_.size()-1);
  } else {
    antennas_[index].setType(type);
  }
}

/**.......................................................................
 * Specify the type of an antenna. but only if it is currently unknown
 */
void ObsInfo::setAntennaTypeIfUnknown(Antenna::AntennaType type, int index)
{
  if(antennas_.size() == 0) {
    ThrowError("The number of antennas is currently unknown.  Use setNumberOfAntennas()");
  }

  // Negative index means default to all

  if(index < 0) {
    for(unsigned i=0; i < antennas_.size(); i++) {
      if(antennas_[i].type_ == Antenna::ANT_UNKNOWN)
	antennas_[i].setType(type);
    }
  } else if(index > antennas_.size()-1) {
    ThrowError("Invalid index: " << index << " should be: 0 - " << antennas_.size()-1);
  } else {
    if(antennas_[index].type_ == Antenna::ANT_UNKNOWN)
      antennas_[index].setType(type);
  }
}

/**.......................................................................
 * Specify the diameter of an antenna
 */
void ObsInfo::setAntennaDiameter(Length diameter, int index)
{
  if(antennas_.size() == 0) {
    ThrowError("The number of antennas is currently unknown.  Use setNumberOfAntennas()");
  }

  // Negative index means default to all

  if(index < 0) {
    for(unsigned i=0; i < antennas_.size(); i++) {
      antennas_[i].setDiameter(diameter);
    }
  } else if(index > antennas_.size()-1) {
    ThrowError("Invalid index: " << index << " should be: 0 - " << antennas_.size()-1);
  } else {
    antennas_[index].setDiameter(diameter);
  }
}

/**.......................................................................
 * Specify the aperture efficiency of an antenna
 */
void ObsInfo::setAntennaApertureEfficiency(Percent apEff, int index)
{
  if(antennas_.size() == 0) {
    ThrowError("The number of antennas is currently unknown.  Use setNumberOfAntennas()");
  }

  // Negative index means default to all

  if(index < 0) {
    for(unsigned i=0; i < antennas_.size(); i++) {
      antennas_[i].setApertureEfficiency(apEff);
    }
  } else if(index > antennas_.size()-1) {
    ThrowError("Invalid index: " << index << " should be: 0 - " << antennas_.size()-1);
  } else {
    antennas_[index].setApertureEfficiency(apEff);
  }
}

/**.......................................................................
 * Set the X coordinate of an antenna
 */
void ObsInfo::setAntennaX(Length X, unsigned index)
{
  if(antennas_.size() == 0) {
    ThrowError("The number of antennas is currently unknown.  Use setNumberOfAntennas()");
  }

  if(index < 0 || index > antennas_.size()-1) {
    ThrowError("Invalid index: " << index << " should be: 0 - " << antennas_.size()-1);
  } else {
    antennas_[index].setX(X);
  }

  checkAntennaLocations();
}

/**.......................................................................
 * Set the Y coordinate of an antenna
 */
void ObsInfo::setAntennaY(Length Y, unsigned index)
{
  if(antennas_.size() == 0) {
    ThrowError("The number of antennas is currently unknown.  Use setNumberOfAntennas()");
  }

  if(index < 0 || index > antennas_.size()-1) {
    ThrowError("Invalid index: " << index << " should be: 0 - " << antennas_.size()-1);
  } else {
    antennas_[index].setY(Y);
  }

  checkAntennaLocations();
}

/**.......................................................................
 * Set the Z coordinate of an antenna
 */
void ObsInfo::setAntennaZ(Length Z, unsigned index)
{
  if(antennas_.size() == 0) {
    ThrowError("The number of antennas is currently unknown.  Use setNumberOfAntennas()");
  }

  if(index < 0 || index > antennas_.size()-1) {
    ThrowError("Invalid index: " << index << " should be: 0 - " << antennas_.size()-1);
  } else {
    antennas_[index].setZ(Z);
  }

  checkAntennaLocations();
}

/**.......................................................................
 * Specify the location of an antenna (in ENU, relative to an array center, or absolute XYZ)
 */
void ObsInfo::setAntennaLocation(LengthTriplet location, int index)
{
  if(antennas_.size() == 0) {
    ThrowError("The number of antennas is currently unknown.  Use setNumberOfAntennas()");
  }

  if(index < 0 || index > antennas_.size()-1) {
    ThrowError("Invalid index: " << index << " should be: 0 - " << antennas_.size()-1);
  } else {
    antennas_[index].setLocation(location);
  }

  checkAntennaLocations();
}

/**.......................................................................
 * Check if all antennas have locations fully specified
 */
void ObsInfo::checkAntennaLocations()
{
  // Check if this update has now left us with locations for all
  // antennas

  bool antsHaveLocation = true;
  for(unsigned iAnt=0; iAnt < antennas_.size(); iAnt++) {
    if(!antennas_[iAnt].hasLocation()) {
      antsHaveLocation = false;
      break;
    }
  }

  // If all antennas now have locations, mark this information as
  // received

  if(antsHaveLocation) {
    infoMask_ |= ObsParameter::INFO_LOCATION_ANTENNA;
  } 
}

/**.......................................................................
 * Return true of antennas have locations specified
 */
bool ObsInfo::antsHaveLocations()
{
  return (infoMask_ & ObsParameter::INFO_LOCATION_ANTENNA);
}

/**.......................................................................
 * Plot antennas
 */
void ObsInfo::plotAntennas()
{
  if(getParameter("dev", false)->data_.hasValue()) {
    PgUtil::open(getStringVal("dev"));
  } else {
    PgUtil::open("/xs");
}

  try {
    plotAntennasEnu();
  } catch(Exception& err) {
    plotAntennasXyz();
  }
}

void ObsInfo::plotAntennasEnu()
{
  std::vector<float> x(antennas_.size());
  std::vector<float> y(antennas_.size());
  double xmn, xmx, ymn, ymx;
  double xmin, xmax, ymin, ymax;
  double radMeters;
  double xCenter=0.0, yCenter=0.0;

  //------------------------------------------------------------
  // First compute the center of the array
  //------------------------------------------------------------

  for(unsigned iAnt=0; iAnt < antennas_.size(); iAnt++) {
    xCenter += (antennas_[iAnt].getEnu().east_.meters()  - xCenter) / (iAnt + 1);
    yCenter += (antennas_[iAnt].getEnu().north_.meters() - yCenter) / (iAnt + 1);
  }

  //------------------------------------------------------------
  // Now calculate the boundaries of this array
  //------------------------------------------------------------

  for(unsigned iAnt=0; iAnt < antennas_.size(); iAnt++) {

    x[iAnt] = antennas_[iAnt].getEnu().east_.meters()  - xCenter;
    y[iAnt] = antennas_[iAnt].getEnu().north_.meters() - yCenter;

    radMeters = antennas_[iAnt].getRadius().meters();

    xmn = x[iAnt] - radMeters;
    xmx = x[iAnt] + radMeters;

    ymn = y[iAnt] - radMeters;
    ymx = y[iAnt] + radMeters;

    // Initialize min/max if this is the first antennas

    if(iAnt == 0) {
      xmin = xmn;
      xmax = xmx;
      ymin = ymn;
      ymax = ymx;
    }

    xmin = (xmin < xmn) ? xmin : xmn;
    xmax = (xmax > xmx) ? xmax : xmx;

    ymin = (ymin < ymn) ? ymin : ymn;
    ymax = (ymax > ymx) ? ymax : ymx;
  }

  double xrange = xmax - xmin;
  double yrange = ymax - ymin;

  xmin = xmin - xrange * 0.1;
  xmax = xmax + xrange * 0.1;

  ymin = ymin - yrange * 0.1;
  ymax = ymax + yrange * 0.1;

  PgUtil::setXmin(xmin);
  PgUtil::setXmax(xmax);

  PgUtil::setYmin(ymin);
  PgUtil::setYmax(ymax);

  PgUtil::setUsedefs(true);

  PgUtil::setInteractive(false);

  unsigned npt = 100;
  std::vector<double> xcirc(npt);
  std::vector<double> ycirc(npt);
  double dtheta = 2*M_PI / (npt-1), theta;

  for(unsigned iAnt=0; iAnt < antennas_.size(); iAnt++) {

    double radMeter = antennas_[iAnt].getRadius().meters();

    for(unsigned i=0; i < npt; i++) {
      theta = dtheta * i;
      xcirc[i] = x[iAnt] + radMeter * cos(theta);
      ycirc[i] = y[iAnt] + radMeter * sin(theta);
    }

    if(iAnt == 0) {

      PgUtil::setOverplot(false);

      PgUtil::setVp(true);
      PgUtil::setWin(true);
      PgUtil::setWnad(true);
      PgUtil::setBox(true);
      
      PgUtil::linePlot(xcirc, ycirc, "East (m)", "North (m)", "");
    } else {

      PgUtil::setOverplot(true);
      PgUtil::setVp(false);
      PgUtil::setWin(false);
      PgUtil::setBox(false);
      PgUtil::linePlot(xcirc, ycirc, "East (m)", "North (m)", "");

    }
  }

  // Unset any global defaults used in this functions

  PgUtil::setUsedefs(false);
  PgUtil::setVp(true);
  PgUtil::setWin(true);
  PgUtil::setWnad(false);
  PgUtil::setBox(true);
}

/**.......................................................................
 * Plot antennas
 */
void ObsInfo::plotAntennasXyz()
{
  std::vector<float> x(antennas_.size());
  std::vector<float> y(antennas_.size());
  double xmn, xmx, ymn, ymx;
  double xmin, xmax, ymin, ymax;
  double radMeters;
  double xCenter=0.0, yCenter=0.0;

  //------------------------------------------------------------
  // First compute the center of the array
  //------------------------------------------------------------

  for(unsigned iAnt=0; iAnt < antennas_.size(); iAnt++) {
    xCenter += (antennas_[iAnt].getXyz().X_.meters()  - xCenter) / (iAnt + 1);
    yCenter += (antennas_[iAnt].getXyz().Y_.meters() - yCenter) / (iAnt + 1);
  }

  //------------------------------------------------------------
  // Now calculate the boundaries of this array
  //------------------------------------------------------------

  for(unsigned iAnt=0; iAnt < antennas_.size(); iAnt++) {

    x[iAnt] = antennas_[iAnt].getXyz().X_.meters()  - xCenter;
    y[iAnt] = antennas_[iAnt].getXyz().Y_.meters() - yCenter;

    radMeters = antennas_[iAnt].getRadius().meters();

    xmn = x[iAnt] - radMeters;
    xmx = x[iAnt] + radMeters;

    ymn = y[iAnt] - radMeters;
    ymx = y[iAnt] + radMeters;

    // Initialize min/max if this is the first antennas

    if(iAnt == 0) {
      xmin = xmn;
      xmax = xmx;
      ymin = ymn;
      ymax = ymx;
    }

    xmin = (xmin < xmn) ? xmin : xmn;
    xmax = (xmax > xmx) ? xmax : xmx;

    ymin = (ymin < ymn) ? ymin : ymn;
    ymax = (ymax > ymx) ? ymax : ymx;
  }

  double xrange = xmax - xmin;
  double yrange = ymax - ymin;

  xmin = xmin - xrange * 0.1;
  xmax = xmax + xrange * 0.1;

  ymin = ymin - yrange * 0.1;
  ymax = ymax + yrange * 0.1;

  PgUtil::setXmin(xmin);
  PgUtil::setXmax(xmax);

  PgUtil::setYmin(ymin);
  PgUtil::setYmax(ymax);

  PgUtil::setUsedefs(true);

  PgUtil::setInteractive(false);

  unsigned npt = 100;
  std::vector<double> xcirc(npt);
  std::vector<double> ycirc(npt);
  double dtheta = 2*M_PI / (npt-1), theta;

  for(unsigned iAnt=0; iAnt < antennas_.size(); iAnt++) {

    double radMeter = antennas_[iAnt].getRadius().meters();

    for(unsigned i=0; i < npt; i++) {
      theta = dtheta * i;
      xcirc[i] = x[iAnt] + radMeter * cos(theta);
      ycirc[i] = y[iAnt] + radMeter * sin(theta);
    }

    if(iAnt == 0) {

      PgUtil::setOverplot(false);

      PgUtil::setVp(true);
      PgUtil::setWin(true);
      PgUtil::setWnad(true);
      PgUtil::setBox(true);
      
      PgUtil::linePlot(xcirc, ycirc, "X (m)", "Y (m)", "");
    } else {

      PgUtil::setOverplot(true);
      PgUtil::setVp(false);
      PgUtil::setWin(false);
      PgUtil::setBox(false);
      PgUtil::linePlot(xcirc, ycirc, "X (m)", "Y (m)", "");

    }
  }

  PgUtil::setUsedefs(false);
}

/**.......................................................................
 * Set the type of noise to generate
 */
void ObsInfo::setNoiseType(std::string type)
{
  if(type=="fixed") {
    setNoiseType(NOISE_FIXED);
  } else if(type=="realistic") {
    setNoiseType(NOISE_REALISTIC);
  } else if(type=="none") {
    setNoiseType(NOISE_NONE);
  } else {
    ThrowError("Unrecognized noise type: " << type);
  }
}

void ObsInfo::setNoiseType(NoiseType type)
{
  noiseType_ = type;
}

/**.......................................................................
 * If the noise type is fixed, this specifies the fixed noise rms to
 * use
 */
void ObsInfo::setFixedNoiseRms(double val, std::string units)
{
  noiseRms_.setVal(val, units);
  infoMask_ |= ObsParameter::INFO_NOISE_FIXED;
}

void ObsInfo::setSourceName(std::string name)
{
  sourceName_ = name;
  infoMask_ |= ObsParameter::INFO_SRC_NAME;
}

/**.......................................................................
 * Set the RA of observation
 */
void ObsInfo::setObsRa(HourAngle obsRa)
{
  obsRa_ = obsRa;
  infoMask_ |= ObsParameter::INFO_OBS_RA;
}

/**.......................................................................
 * Set the DEC of observation
 */
void ObsInfo::setObsDec(Declination dec)
{
  obsDec_ = dec;
  infoMask_ |= ObsParameter::INFO_OBS_DEC;
}

/**.......................................................................
 * Set the equinox of observation
 */
void ObsInfo::setObsEquinox(double equinox)
{
  obsEquinox_ = equinox;
  infoMask_ |= ObsParameter::INFO_OBS_EQN;
}

double ObsInfo::getObsEquinox()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_OBS_EQN);
  return obsEquinox_;
}

void ObsInfo::setObsHa(gcp::util::HourAngle startHa, gcp::util::HourAngle stopHa, 
		       gcp::util::HourAngle deltaHa)
{
  startHa_ = startHa;
  stopHa_  = stopHa;
  deltaHa_ = deltaHa;
  infoMask_ |= ObsParameter::INFO_OBS_HA;
}

/**.......................................................................
 * Set the array location
 */
void ObsInfo::setArrayLocation(Lla lla)
{
  lla_ = lla;
  infoMask_ |= ObsParameter::INFO_LOCATION_ARRAY;

  // Pass the array location down to any antennas that have been added so far

  for(unsigned iAnt=0; iAnt < antennas_.size(); iAnt++) {
    antennas_[iAnt].setReferenceLla(lla_);
  }

}

void ObsInfo::setTelescopeName(std::string name)
{
  telescopeName_ = name;
  infoMask_ |= ObsParameter::INFO_TELESCOPE;
}

void ObsInfo::setInstrumentName(std::string name)
{
  instrumentName_ = name;
  infoMask_ |= ObsParameter::INFO_INSTRUMENT;
}

void ObsInfo::simulateNoise(bool simNoise)
{
  simNoise_ = simNoise;
}

void ObsInfo::setAmbientTemperature(Temperature& tAmb)
{
  tAmb_ = tAmb;
  infoMask_ |= ObsParameter::INFO_NOISE_TAMB;
}

void ObsInfo::setOpacity(Percent& tau)
{
  tau_ = tau;
  infoMask_ |= ObsParameter::INFO_NOISE_TAU;
}

/**.......................................................................
 * Generate noise for this observation
 */
void ObsInfo::getNoiseRms(Flux& noiseRms,
			  HourAngle* ha, Declination* dec, Frequency* bandwidth, Time* dt,
			  Antenna* ant1, Antenna* ant2, PolarLengthVector* azel)
{
  static unsigned count=0;

  double noiseRmsJy = 0.0;

  //------------------------------------------------------------
  // Generate noise according to type
  //------------------------------------------------------------

  switch(noiseType_) {
    
  //------------------------------------------------------------
  // Fixed noise
  //------------------------------------------------------------

  case NOISE_FIXED:
    {
      if(ObsParameter::parametersAreMissing(infoMask_, ObsParameter::SIM_REQ_NOISE_FIXED))
	ThrowError("Not enough information has been given to calculate a fixed noise: this = " << this);

      noiseRmsJy = noiseRms_.getVal("Jy");
    }
    break;
    
    //------------------------------------------------------------
    // Realistic noise
    //------------------------------------------------------------
    
  case NOISE_REALISTIC:
    {
      if(ObsParameter::parametersAreMissing(infoMask_, ObsParameter::SIM_REQ_NOISE_REALISTIC)) {
	ObsParameter::printMissingParameters(infoMask_, ObsParameter::SIM_REQ_NOISE_REALISTIC);
	ThrowError("Not enough information has been given to calculate realistic noise");
      }
      
      if(ha==0 || dec==0 || ant1==0 || ant2==0)
	ThrowError("Null parameters passed to generateNoise(): Not enough information has"
		   "been given to calculate realistic noise");
      
      if(!(ant1->canComputeNoise() && ant2->canComputeNoise())) {
	ThrowError("Antenna parameters required for noise calculation have not been specified");
      }
      
#ifdef SIM_TIMER_TEST
      t1.start();
#endif

      PolarLengthVector azEl1;
      PolarLengthVector azEl2; 

      if(azel != 0) {
	azEl1 = *azel;
	azEl2 = *azel;
      } else {
	azEl1 = ant1->getAzEl(ha, dec);
	azEl2 = ant2->getAzEl(ha, dec);
      }
      
#ifdef SIM_TIMER_TEST
      t1.stop();
      t1time += t1.deltaInSeconds();
      t2.start();
#endif
      
      Temperature tatm1 = tAmb_ * (tau_.percentMax1() * 1.0/sin(azEl1.el_.radians()));
      Temperature tatm2 = tAmb_ * (tau_.percentMax1() * 1.0/sin(azEl2.el_.radians()));
      
#ifdef SIM_TIMER_TEST
      t2.stop();
      t2time += t2.deltaInSeconds();
      t3.start();
#endif
      
      Temperature t1 = ant1->getReceiverTemperature() + ant1->getGroundSpillover() + tatm1 + Constants::Tcmb_;
      Temperature t2 = ant2->getReceiverTemperature() + ant2->getGroundSpillover() + tatm2 + Constants::Tcmb_;
      
      double noiseRmsJy1 = (ant1->jyPerK().Jy() * t1.K()) / sqrt(dt->seconds() * bandwidth->Hz());
      double noiseRmsJy2 = (ant2->jyPerK().Jy() * t2.K()) / sqrt(dt->seconds() * bandwidth->Hz());

      //      COUT("Ant1 type = " << ant1->type_ << " t1 = " << t1 << " noiseRmsIj = " << noiseRmsJy1 << " dt = " << dt->seconds());
      //      COUT("Ant1 type = " << ant2->type_ << " t2 = " << t2 << " noiseRmsIj = " << noiseRmsJy2 << " bw = " << bandwidth->Hz());

      noiseRmsJy = sqrt(noiseRmsJy1 * noiseRmsJy2);

#ifdef SIM_TIMER_TEST
  t3.stop();
  t3time += t3.deltaInSeconds();
#endif
      
#if 0
      if(count < 10) {
	++count;
	
	COUT("ha = " << *ha << " dec = " << *dec << " El = " << azEl1.el_ << " ant1 lla = " << ant1->getReferenceLla());
	COUT("Calculated REAL Noise: from rms = " << noiseRmsJy << " Jy");
	COUT("tAmb_  = " << tAmb_);
	COUT("tau_  = " << tau_.percentMax1());
	COUT("trx   = " << ant1->getReceiverTemperature().K());
	COUT("spill = " << ant1->getGroundSpillover().K());
	
	COUT("tatm1 = " << tatm1 << " t1 = " << t1);
	COUT("tatm2 = " << tatm2 << " t2 = " << t2);
	COUT("ant1 jyperk = " << ant1->jyPerK().Jy());
	COUT("ant1 type   = " << ant1->typeStr());
	COUT("dt = " << dt->seconds() << " s");
	COUT("bw = " << bandwidth->GHz() << " GHz");
	
      }
#endif
    }
    break;
    
    //------------------------------------------------------------
    // No noise
    //------------------------------------------------------------
    
  default:
    noiseRmsJy = 0.0;
    break;
  }

  noiseRms.setJy(noiseRmsJy);
}

/**.......................................................................
 * Given a noise rms, generate random samples of noise consistent with it
 */
void ObsInfo::generateNoise(Flux& noiseRms, Flux& reNoise, Flux& imNoise, double& wt)
{
  double eps = 1e-12;
  double noiseRmsJy = noiseRms.Jy();

  if(noiseRmsJy < eps) {
    reNoise.setJy(0.0);
    imNoise.setJy(0.0);
    wt = 1.0;
  } else {

    // We assume that the noise rms being specified is the usual
    // Tsys/sqrt(t * bw).  For single-dish, this is the noise
    // associated with the intensity.  For an interferometer, the
    // variance of the Re and Im must sum to this variance, or
    // equivalently, the rms noise on the Re and Im is this rms,
    // divided by sqrt(2).  

    reNoise.setJy(Sampler::generateGaussianSample(noiseRmsJy/sqrt(2.0)));
    imNoise.setJy(Sampler::generateGaussianSample(noiseRmsJy/sqrt(2.0)));

    // The single weight that will be written with the data is the
    // inverse of the Re/Im variance.

    wt = 2.0/(noiseRmsJy*noiseRmsJy);
  }
}

/**.......................................................................
 * Generate noise for this observation
 */
void ObsInfo::generateNoise(Flux& reNoise, Flux& imNoise, double& wt, 
			    HourAngle* ha, Declination* dec, Frequency* bandwidth, Time* dt,
			    Antenna* ant1, Antenna* ant2)
{
  Flux noiseRms;
  getNoiseRms(noiseRms, ha, dec, bandwidth, dt, ant1, ant2);
  generateNoise(noiseRms, reNoise, imNoise, wt);
}

void ObsInfo::calculateStartJd()
{
  // Calculate the LST at the site, for the current time

  TimeVal time;
  time.setToCurrentTime();
  HourAngle lst = Astrometry::mjdUtcToLst(time.getMjd(), lla_.longitude_, 0.0, 0.0);

  // Now calculate the last absolute time corresponding to the requested startHa

  HourAngle deltaLst = startHa_ + obsRa_ - lst;

  double mjdStart = time.getMjd() + deltaLst.hours()/24;

  startJd_ = mjdStart + 2400000.5;
}

/**.......................................................................
 * Set frequency information
 */
void ObsInfo::setFrequencyInformation(std::vector<gcp::util::Frequency>& freqs,
				      std::vector<gcp::util::Frequency>& bws)
{
  if(freqs.size() != bws.size()) {
    ThrowError("Frequency and bandwidth arrays must be the same size");
  }

  frequencies_ = freqs;
  bandwidths_  = bws;

  nFreq_ = frequencies_.size();

  infoMask_ |= ObsParameter::INFO_FREQ;
}

/**.......................................................................
 * Set the number of groups in this observation
 */
void ObsInfo::setNumberOfGroups(unsigned nGroup)
{
  nGroup_ = nGroup;
  infoMask_ |= ObsParameter::INFO_NGROUP;
}

/**.......................................................................
 * Get the number of visibility groups in this observation
 */
unsigned ObsInfo::getNumberOfGroups()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_NGROUP);
  return nGroup_;
}

/**.......................................................................
 * Get the number of frequencies in this observation
 */
unsigned ObsInfo::getNumberOfFrequencies()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_FREQ);
  return nFreq_;
}

/**.......................................................................
 * Set the number of timestamps in this observation
 */
void ObsInfo::setNumberOfTimestamps(unsigned nTimestamp)
{
  // The number of groups is primary.  If we don't know the number of
  // antennas, then we can't convert from timestamps to number of
  // groups ( = ntimestamp x nbaseline)

  if(infoMask_ & ObsParameter::INFO_NANT) {
    nGroup_ = nTimestamp * nBaseline_;
    infoMask_ |= ObsParameter::INFO_NGROUP;
  } else {
    ThrowError("The number of antennas isn't know.  The number of groups cannot be calculated");
  }

}

/**.......................................................................
 * Get the number of timestamps in this observation
 */
unsigned ObsInfo::getNumberOfTimestamps()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_NGROUP | ObsParameter::INFO_NANT);
  return nGroup_ / nBaseline_;
}

/**.......................................................................
 * Set number of stokes parameters
 */
void ObsInfo::setNumberOfStokesParameters(unsigned nStokes)
{
  nStokes_ = nStokes;
  infoMask_ |= ObsParameter::INFO_NSTOKES;
}

/**.......................................................................
 * Get the number of stokes parameters in this observation
 */
unsigned ObsInfo::getNumberOfStokesParameters()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_NSTOKES);
  return nStokes_;
}

/**.......................................................................
 * Initialize internal arrays for holding simulated visibility data
 */
void ObsInfo::initializeSimulationVisibilityArray()
{
  ObsParameter::checkParameters(infoMask_, 
				ObsParameter::INFO_OBS_HA | 
				ObsParameter::INFO_NANT | 
				ObsParameter::INFO_FREQ | 
				ObsParameter::INFO_NSTOKES);

  unsigned nHa = (unsigned)((stopHa_ - startHa_) / deltaHa_);

  //------------------------------------------------------------
  // Reset size parameters for the underlaying visibility data set
  //------------------------------------------------------------

  setNumberOfTimestamps(nHa);

  //------------------------------------------------------------
  // Resize the visibility array to match the number of simulated
  // observations we will have
  //------------------------------------------------------------

  visibilities_.resize(nGroup_);

  //------------------------------------------------------------
  // Each vis group represents all IFs and Stokes parameters of a
  // single baseline at a specific timestamp
  //------------------------------------------------------------

  for(unsigned iVisGroup=0; iVisGroup < nGroup_; iVisGroup++) {

    ObsInfo::Vis& vis = visibilities_[iVisGroup];

    vis.dims_.resize(2);

    vis.dims_[0] = nStokes_;
    vis.dims_[1] = nFreq_;

    vis.re_.resize(nStokes_ * nFreq_);
    vis.im_.resize(nStokes_ * nFreq_);
    vis.wt_.resize(nStokes_ * nFreq_);
  }
}

float ObsInfo::getU(unsigned iFrame, unsigned iBaseline)
{
  unsigned groupIndex = iFrame * nBaseline_ + iBaseline;
  return (float) visibilities_[groupIndex].u_;
}

float ObsInfo::getV(unsigned iFrame, unsigned iBaseline)
{
  unsigned groupIndex = iFrame * nBaseline_ + iBaseline;
  return (float) visibilities_[groupIndex].v_;
}

float ObsInfo::getW(unsigned iFrame, unsigned iBaseline)
{
  unsigned groupIndex = iFrame * nBaseline_ + iBaseline;
  return (float) visibilities_[groupIndex].w_;
}

double ObsInfo::getJulianDate(unsigned iFrame)
{
  unsigned groupIndex = iFrame * nBaseline_;
  return visibilities_[groupIndex].jd_;
}

double ObsInfo::getJulianDate(unsigned iFrame, unsigned iBaseline)
{
  unsigned groupIndex = iFrame * nBaseline_ + iBaseline;
  return visibilities_[groupIndex].jd_;
}

float ObsInfo::getRe(unsigned iFrame, unsigned iBaseline, unsigned iIf, unsigned iStokes)
{
  unsigned groupIndex = iFrame * nBaseline_ + iBaseline;
  unsigned visIndex  = iStokes * nFreq_ + iIf;
  return (float) visibilities_[groupIndex].re_[visIndex];
}

float ObsInfo::getIm(unsigned iFrame, unsigned iBaseline, unsigned iIf, unsigned iStokes)
{
  unsigned groupIndex = iFrame * nBaseline_ + iBaseline;
  unsigned visIndex  = iStokes * nFreq_ + iIf;
  return (float) visibilities_[groupIndex].im_[visIndex];
}

float ObsInfo::getWt(unsigned iFrame, unsigned iBaseline, unsigned iIf, unsigned iStokes)
{
  unsigned groupIndex = iFrame * nBaseline_ + iBaseline;
  unsigned visIndex  = iStokes * nFreq_ + iIf;
  return (float) visibilities_[groupIndex].wt_[visIndex];
}

float ObsInfo::getAipsBaselineIndex(unsigned iFrame, unsigned iBaseline)
{
  unsigned groupIndex = iFrame * nBaseline_ + iBaseline;
  return (float) visibilities_[groupIndex].baseline_;
}

std::string ObsInfo::getTelescopeName()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_TELESCOPE);
  return telescopeName_;
}

std::string ObsInfo::getInstrumentName()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_INSTRUMENT);
  return instrumentName_;
}

Lla ObsInfo::getArrayLocation()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_LOCATION_ARRAY);
  return lla_;
}

HourAngle ObsInfo::getStartHa()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_OBS_HA);
  return startHa_;
}

HourAngle ObsInfo::getStopHa()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_OBS_HA);
  return stopHa_;
}

HourAngle ObsInfo::getDeltaHa()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_OBS_HA);
  return deltaHa_;
}

std::string ObsInfo::getSourceName()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_SRC_NAME);
  return sourceName_;
}

HourAngle ObsInfo::getObsRa()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_OBS_RA);
  return obsRa_;
}

Declination ObsInfo::getObsDec()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_OBS_DEC);
  return obsDec_;
}

std::vector<LengthTriplet> ObsInfo::getAntennaXyz()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_LOCATION_ANTENNA);

  std::vector<LengthTriplet> xyzs(antennas_.size());

  for(unsigned iAnt=0; iAnt < antennas_.size(); iAnt++) {
    xyzs[iAnt] = antennas_[iAnt].getXyz();
  }

  return xyzs;
}

std::vector<Frequency> ObsInfo::getFrequencies()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_FREQ);
  return frequencies_;
}

std::vector<Frequency> ObsInfo::getBandwidths()
{
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_FREQ);
  return bandwidths_;
}

//-----------------------------------------------------------------------
// Vis operators
//-----------------------------------------------------------------------

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, ObsInfo::Vis& vis)
{
  os << "u  = " << vis.u_ << std::endl;
  os << "v  = " << vis.v_ << std::endl;
  os << "w  = " << vis.w_ << std::endl;

  os << "freq[0] = " << vis.freq_[0] << std::endl;
  os << "freq[1] = " << vis.freq_[1] << std::endl;
  os << "if[0]   = " << vis.if_[0]   << std::endl;
  os << "if[1]   = " << vis.if_[1]   << std::endl;

  os << "baseline = " << vis.baseline_ << std::endl;
  os << "JD = " << std::setprecision(12) << vis.jd_ << std::endl;
  return os;
}

void ObsInfo::Vis::operator=(const FitsUvfReader::Vis& vis)
{
  *this = (FitsUvfReader::Vis&) vis;
}

void ObsInfo::Vis::operator=(FitsUvfReader::Vis& vis)
{
  dims_ = vis.dims_;

  u_ = vis.u_;
  v_ = vis.v_;
  w_ = vis.w_;

  baseline_ = vis.baseline_;

  jd_ = vis.jd_;
  
  index_  = vis.index_;
  re_     = vis.re_;
  im_     = vis.im_;
  wt_     = vis.wt_;
  stokes_ = vis.stokes_;
  freq_   = vis.freq_;
  if_     = vis.if_;
  ra_     = vis.ra_;
  dec_    = vis.dec_;
}

/**.......................................................................
 * Return true if we have enough information to simulate
 */
bool ObsInfo::canSimulate()
{
  switch (noiseType_) {
  case NOISE_FIXED:
    return !ObsParameter::parametersAreMissing(infoMask_, ObsParameter::SIM_REQ_INFO | ObsParameter::SIM_REQ_NOISE_FIXED);
    break;
  case NOISE_REALISTIC:
    return !ObsParameter::parametersAreMissing(infoMask_, ObsParameter::SIM_REQ_INFO | ObsParameter::SIM_REQ_NOISE_REALISTIC);
    break;
  default:
    return false;
    break;
  }
}

void ObsInfo::seed(int s)
{
  Sampler::seed(s);
}

void ObsInfo::writeUvfFile(std::string fileName)
{
  FitsIoHandler fitsio;
  fitsio.writeUvfFile(fileName, *this);
}

//-----------------------------------------------------------------------
// External calling interface
//-----------------------------------------------------------------------

/**.......................................................................
 * Initialize all parameters we will make available to the script
 * calling interface
 */
void ObsInfo::addParameters()
{
  addParameter("dev",            DataType::STRING, "Pgplot device to use when displaying this object");
  addParameter("instrumentName", DataType::STRING);
  addParameter("telescopeName",  DataType::STRING);

  addParameter("latitude",       DataType::STRING);
  addParameter("longitude",      DataType::STRING);
  addParameter("altitude",       DataType::FLOAT);

  addParameter("haStart",        DataType::STRING,  "The start HA to use when simulating data");
  addParameter("haStop",         DataType::STRING,  "The stop HA to use when simulating data");
  addParameter("haDelta",        DataType::STRING,  "The HA delta to use when simulating data");

  addParameter("sourceName",     DataType::STRING);
  addParameter("obsRa",          DataType::STRING);
  addParameter("obsDec",         DataType::STRING);
  addParameter("obsEquinox",     DataType::FLOAT);

  addParameter("nant",           DataType::UINT,    "Number of antennas in the array");

  addParameter("plot",           DataType::BOOL,    "True to plot the simulated array");
  addParameter("freqs",          DataType::STRING,  "A list of frequencies to simulate");
  addParameter("bws",            DataType::STRING,  "The bandwidths to use for each frequency channel");

  addParameter("noiseType",      DataType::STRING,  "The type of noise to simulate");
  addParameter("noiseRms",       DataType::DOUBLE,  "The RMS noise value to use if noiseType = fixed");
  addParameter("tambient",       DataType::DOUBLE,  "The ambient temperature to use if noiseType = real");
  addParameter("tau",            DataType::DOUBLE,  "The atmospheric opacity to use if noiseType = real");

  addParameter("seed",           DataType::UINT,    "An explicit seed for the random number generator");

  addParameter("nstokes",        DataType::UINT,    "The number of Stokes parameters to simulate");
}

/**.......................................................................
 * Add a parameter to our map.  If resizable=true, we create a new
 * vector of the specified type.
 */
void ObsInfo::setParameter(std::string name, std::string val, std::string units)
{
  String nameStr(name);

  //------------------------------------------------------------
  // Always call the underlying PM method:
  //------------------------------------------------------------

  ParameterManager::setParameter(name, val, units);

  //------------------------------------------------------------
  // Special keywords need additional action.  If 'nant' is specified,
  // resize our array of antennas, and add them to our parameter map,
  // so that they and their parameters can be accessed by name
  //------------------------------------------------------------

  if(name == "nant") {
    unsigned nant = getUintVal("nant");
    setNumberOfAntennas(nant);

    std::ostringstream os;
    for(unsigned iAnt=0; iAnt < antennas_.size(); iAnt++) {

      antennas_[iAnt].addParameters();

      os.str("");
      os << "ant" << iAnt+1;
      ParameterManager::addParameter(os.str(), antennas_[iAnt]);
    }

    if(nant > 0) {

      //------------------------------------------------------------
      // Also add generic versions (i.e., ant.apeff) so that the user
      // can short-hand parameters for all ants
      //------------------------------------------------------------

      ParameterManager::addParameter("ant", antennas_[0]);
    }

  } else if(name == "instrumentName") {
    setInstrumentName(getStringVal("instrumentName"));
  } else if(name == "telescopeName") {
    setTelescopeName(getStringVal("telescopeName"));
  } else if(name == "sourceName") {
    setSourceName(getStringVal("sourceName"));
  } else if(name == "latitude") {
    checkLocationParameters();
  } else if(name == "longitude") {
    checkLocationParameters();
  } else if(name == "altitude") {
    checkLocationParameters();
  } else if(name == "haStart") {
    checkHaParameters();
  } else if(name == "haStop") {
    checkHaParameters();
  } else if(name == "haDelta") {
    checkHaParameters();
  } else if(name == "obsRa") {
    HourAngle obsRa;
    obsRa.setHours(getStringVal("obsRa"));
    setObsRa(obsRa);
  } else if(name == "obsDec") {
    Declination obsDec;
    obsDec.setDegrees(getStringVal("obsDec"));
    setObsDec(obsDec);
  } else if(name == "obsEquinox") {
    setObsEquinox(getFloatVal("obsEquinox"));
  } else if(name == "freqs" || name == "bws") {
    checkFreqParameters();
  } else if(name == "noiseType") {
    setNoiseType(getStringVal("noiseType"));
  } else if(name == "noiseRms") {
    setFixedNoiseRms(getDoubleVal("noiseRms"), getParameter("noiseRms", true)->units_);
  } else if(name == "tambient") {
    Temperature temp;
    temp.setVal(getDoubleVal("tambient"), getParameter("tambient", true)->units_);
    setAmbientTemperature(temp);
  } else if(name == "tau") {
    Percent tau;
    tau.setPercentMax1(getDoubleVal("tau"));
    setOpacity(tau);
  } else if(name == "seed") {
    seed(getUintVal("seed"));
  } else if(name == "nstokes") {
    setNumberOfStokesParameters(getUintVal("nstokes"));
  } else if(name == "plot") {

    if(getBoolVal("plot")) {
      plotAntennas();
    }

    // Call down explicitly to the antenna named in the argument, or
    // its overloaded setParameter method won't get called

  } else if(nameStr.contains("ant") && nameStr.contains(".")) {

    String antNo   = nameStr.findNextInstanceOf("ant", true, ".", true, true);
    String varName = nameStr.remainder();

    // Try to set the named variable.  However, if antNo cannot be
    // converted to a valid int, then we assume the request is to set
    // the parameter for all antennas

    try {
      unsigned iAnt = antNo.toInt()-1;
      antennas_[iAnt].setParameter(varName.str(), val, units);
    } catch(...) {

      if(antennas_.size() == 0) {
	ThrowError("No antennas have currently been added to the array.  Use 'nant' to add them first");
      }

      for(unsigned iAnt=0; iAnt < antennas_.size(); iAnt++) {
	antennas_[iAnt].setParameter(varName.str(), val, units);
      }
    }
  }

  checkAntennaLocations();
}

void ObsInfo::checkLocationParameters()
{
  if(getParameter("latitude",  false)->data_.hasValue() && 
     getParameter("longitude", false)->data_.hasValue() && 
     getParameter("altitude",  false)->data_.hasValue()) {

    Lla lla;

    lla.latitude_.setVal(getStringVal("latitude"));
    lla.longitude_.setVal(getStringVal("longitude"));
    lla.altitude_.setVal(getFloatVal("altitude"), getParameter("altitude", true)->units_);

    setArrayLocation(lla);
  } 
}

/**.......................................................................
 * Method to check validity for externally specified HA range and delta
 */
void ObsInfo::checkHaParameters()
{
  if(getParameter("haStart",  false)->data_.hasValue() && 
     getParameter("haStop",   false)->data_.hasValue() && 
     getParameter("haDelta",  false)->data_.hasValue()) {

    HourAngle start, stop, delta;

    String startStr = getStringVal("haStart");
    String stopStr  = getStringVal("haStop");
    String deltaStr = getStringVal("haDelta");

    if(startStr.contains(":")) {
      start.setHours(startStr.str());
    }

    if(stopStr.contains(":")) {
      stop.setHours(stopStr.str());
    }

    if(deltaStr.contains(":")) {
      delta.setHours(deltaStr.str());
    }

    setObsHa(start, stop, delta);
  } 
}

/**.......................................................................
 * Method to check validity for externally specified frequencies and
 * bandwidths
 */
void ObsInfo::checkFreqParameters()
{
  if(getParameter("freqs", false)->data_.hasValue() && getParameter("bws", false)->data_.hasValue()) {

    String freqs = getStringVal("freqs");
    std::vector<double> freqVals = parseRange(freqs);

    String bws   = getStringVal("bws");
    std::vector<double> bwVals = parseRange(bws);

    if(bwVals.size() == 1) {
      double val = bwVals[0];
      bwVals.resize(freqVals.size());
      for(unsigned i=0; i < bwVals.size(); i++) {
	bwVals[i] = val;
      }
    }

    if(bwVals.size() != freqVals.size()) {
      ThrowSimpleColorError("You must specify a bandwidth for each value of the frequency array, or a single bandwidth that applies to all frequencies", "red");
    }
 
    std::vector<Frequency> freqArr(freqVals.size());
    std::vector<Frequency> bwArr(bwVals.size());
    
    for(unsigned i=0; i < freqVals.size(); i++) {
      freqArr[i].setVal(freqVals[i], getParameter("freqs", true)->units_);
      bwArr[i].setVal(bwVals[i], getParameter("bws", true)->units_);
    }

    setFrequencyInformation(freqArr, bwArr);
  }
}

void ObsInfo::operator=(const ObsInfo& obs)
{
  operator=((ObsInfo&) obs);
}

void ObsInfo::operator=(ObsInfo& obs)
{
  visibilities_ = obs.visibilities_;
  infoMask_     = obs.infoMask_;
  
  // The location of the array

  lla_ = obs.lla_;

  // The name of the array

  telescopeName_ = obs.telescopeName_;

  // The name of the instrument used for this observation

  instrumentName_ = obs.instrumentName_;

  //------------------------------------------------------------
  // Frequency information
  //------------------------------------------------------------

  frequencies_ = obs.frequencies_;
  bandwidths_  = obs.bandwidths_;

  //------------------------------------------------------------
  // The vector of antennas in the array
  //------------------------------------------------------------

  antennas_    = obs.antennas_;

  //------------------------------------------------------------
  // Observation parameters
  //------------------------------------------------------------

  startHa_     = obs.startHa_;
  stopHa_      = obs.stopHa_;
  deltaHa_     = obs.deltaHa_;

  sourceName_  = obs.sourceName_;
  obsRa_       = obs.obsRa_;
  obsDec_      = obs.obsDec_;
  obsEquinox_  = obs.obsEquinox_;

  //------------------------------------------------------------
  // Noise parameters
  //------------------------------------------------------------

  noiseType_   = obs.noiseType_;
  startJd_     = obs.startJd_;
  simNoise_    = obs.simNoise_;
  noiseRms_    = obs.noiseRms_;
  tAmb_        = obs.tAmb_;
  tau_         = obs.tau_;

  //------------------------------------------------------------
  // Size parameters
  //------------------------------------------------------------

  nGroup_    = obs.nGroup_;     // The number of visibility groups in this observation
  nBaseline_ = obs.nBaseline_;  // The number of baselines in this observation
  nFreq_     = obs.nFreq_;      // The number of frequencies in this observation
  nStokes_   = obs.nStokes_;    // The number of Stokes parameters in this observation
}

void ObsInfo::printTime()
{
#ifdef ObsParameter::SIM_TIMER_TEST
  COUT("T1 OBS= " << t1time);
  COUT("T2 OBS= " << t2time);
  COUT("T3 OBS= " << t3time);
#endif
}

double ObsInfo::getFixedNoiseRms(gcp::util::Unit::Units units)
{
  return getFixedNoiseRms(Unit::unitsToString(units));
}

double ObsInfo::getFixedNoiseRms(std::string units)
{
  try {
    if(noiseType_ == NOISE_FIXED)
      return noiseRms_.getVal(units);
    
    ThrowSimpleColorError("Noise type is not fixed.  Use '" << name_ << ".noiseType = fixed' to change this", "red");

  } catch(Exception& err) {
    XtermManip xtm;
    ostringstream os;
    os << err.what() << COLORIZE(xtm, "red", std::endl << std::endl << "(While retrieving fixed noise rms from " << noiseRms_.name_ << ")");
    ThrowError(os.str());
  }

  return 0.0;
}

void ObsInfo::initializeBaselineIndices()
{
  unsigned nAnt = antennas_.size();

  unsigned iBase=0;
  for(int iAnt1=0; iAnt1 < nAnt-1; iAnt1++) {
    for(int iAnt2=iAnt1+1; iAnt2 < nAnt; iAnt2++) {
      groupToAipsBaselineIndexMap_[iBase] = (iAnt1+1) * 256 + (iAnt2+1);
      ++iBase;
    }
  }

  nBaseline_ = (nAnt * (nAnt-1))/2;
}

/**.......................................................................
 * Return the unique merger of two sets of antennas.  We define this
 * as the resulting array should have the largest number of entries of
 * each antenna type that is encountered in both arrays
 */
std::vector<Antenna> ObsInfo::mergeAnts(ObsInfo& obs)
{
  std::map<Antenna*, unsigned> thisMap = getAntMap();
  std::map<Antenna*, unsigned> thatMap = obs.getAntMap();

  std::vector<Antenna> ants = antennas_;

  for(std::map<Antenna*, unsigned>::iterator iter=thatMap.begin(); iter != thatMap.end(); iter++) {
    Antenna* ant = iter->first;

    //------------------------------------------------------------
    // If this antenna doesn't exist in our map, add as many instances
    // of it as occur in that map
    //------------------------------------------------------------

    if(!existsInAntMap(ant, thisMap)) {

      unsigned nThat = thatMap[getAntMapKey(ant, thatMap)];

      for(unsigned iAdd=0; iAdd < nThat; iAdd++)
	ants.push_back(*ant);

      //------------------------------------------------------------
      // Else add as many more as we need to accomodate the number in
      // that object
      //------------------------------------------------------------

    } else {

      unsigned nThis = thisMap[getAntMapKey(ant, thisMap)];
      unsigned nThat = thatMap[getAntMapKey(ant, thatMap)];

      if(nThat > nThis) {
	for(unsigned iAdd=0; iAdd < (nThat-nThis); iAdd++)
	  ants.push_back(*ant);
      }
    }
  }

  return ants;
}

std::map<Antenna*, unsigned> ObsInfo::getAntMap()
{
  std::map<Antenna*, unsigned> antMap;

  for(unsigned iAnt=0; iAnt < antennas_.size(); iAnt++) {
    if(existsInAntMap(&antennas_[iAnt], antMap))
      antMap[getAntMapKey(&antennas_[iAnt], antMap)] += 1;
    else
      antMap[&antennas_[iAnt]] = 1;
  }

  return antMap;
}

/**.......................................................................
 * Return true if an antenna of this type (or size if type undetermined)
 * exists in our map
 */
bool ObsInfo::existsInAntMap(Antenna* ant, std::map<Antenna*, unsigned>& antMap)
{
  for(std::map<Antenna*, unsigned>::iterator iter=antMap.begin(); iter != antMap.end(); iter++)
    if(*(iter->first) == *ant)
      return true;

  return false;
}

/**.......................................................................
 * Return true if an antenna of this type (or size if type undetermined)
 * exists in our map
 */
Antenna* ObsInfo::getAntMapKey(Antenna* ant, std::map<Antenna*, unsigned>& antMap)
{
  for(std::map<Antenna*, unsigned>::iterator iter=antMap.begin(); iter != antMap.end(); iter++)
    if(*iter->first == *ant)
      return iter->first;

  ThrowError("No matching key found");
  return 0;
}

/**.......................................................................
 * Return true if an antenna of this type (or size if type undetermined)
 * exists in our map
 */
Antenna* ObsInfo::getAntMapKey(Antenna* ant, std::map<Antenna*, std::vector<Antenna*> >& antMap)
{
  for(std::map<Antenna*, std::vector<Antenna*> >::iterator iter=antMap.begin(); iter != antMap.end(); iter++)
    if(*iter->first == *ant)
      return iter->first;

  ThrowError("No matching key found");
  return 0;
}

std::map<Antenna*, std::vector<Antenna*> > ObsInfo::getAntVecMap()
{
  std::map<Antenna*, unsigned> antMap = getAntMap();
  std::map<Antenna*, std::vector<Antenna*> > vecMap;
  
  //------------------------------------------------------------
  // For each unique ant type, construct a vector of antennas in our
  // array that match the type
  //------------------------------------------------------------

  for(std::map<Antenna*, unsigned>::iterator iter=antMap.begin(); iter != antMap.end(); iter++) {
    Antenna* mapAnt = iter->first;

    for(unsigned iAnt=0; iAnt < antennas_.size(); iAnt++) {
      Antenna& vecAnt = antennas_[iAnt];

      if(*mapAnt == vecAnt)
	vecMap[mapAnt].push_back(&vecAnt);
    }
  }
  
  return vecMap;
}

void ObsInfo::setNumberOfBaselines(unsigned nBase)
{
  nBaseline_ = nBase;
}

void ObsInfo::markAntennaLocationsAsReceived()
{
  infoMask_ |= ObsParameter::INFO_LOCATION_ANTENNA;
}
