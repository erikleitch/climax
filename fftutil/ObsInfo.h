// $Id: ObsInfo.h,v 1.3 2012/05/29 16:20:28 eml Exp $

#ifndef GCP_UTIL_OBSINFO_H
#define GCP_UTIL_OBSINFO_H

/**
 * @file ObsInfo.h
 * 
 * Tagged: Thu Feb  3 13:57:05 PST 2011
 * 
 * @version: $Revision: 1.3 $, $Date: 2012/05/29 16:20:28 $
 * 
 * @author Erik Leitch
 */
#include "gcp/fftutil/Antenna.h"

#include "gcp/util/FitsUvfReader.h"
#include "gcp/util/Flux.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/GenericNormalization.h"
#include "gcp/util/Geoid.h"
#include "gcp/util/ParameterManager.h"
#include "gcp/util/Percent.h"

namespace gcp {
  namespace util {

    class ObsInfo : public ParameterManager {
    public:

      //------------------------------------------------------------
      // Struct defining a visibility for a single timestamp
      //------------------------------------------------------------

      struct Vis {

	std::vector<unsigned> dims_;

	double u_;
	double v_;
	double w_;
	unsigned int baseline_;
	double jd_;

	std::vector<unsigned> index_;
	std::vector<double> re_;
	std::vector<double> im_;
	std::vector<double> wt_;
	std::vector<double> stokes_;
	std::vector<double> freq_;
	std::vector<double> if_;
	std::vector<double> ra_;
	std::vector<double> dec_;

	void operator=(const gcp::util::FitsUvfReader::Vis& vis);
	void operator=(gcp::util::FitsUvfReader::Vis& vis);

	friend std::ostream& operator<<(std::ostream& os, Vis& vis);

	void initialize(unsigned nVis) {
	  re_.resize(nVis);
	  im_.resize(nVis);
	  wt_.resize(nVis);

	  u_  = 0.0;
	  v_  = 0.0;
	  w_  = 0.0;
	  jd_ = 0.0;
	  baseline_ = 0;

	  for(unsigned iVis=0; iVis < nVis; iVis++) {
	    re_[iVis] = 0.0;
	    im_[iVis] = 0.0;
	    wt_[iVis] = 0.0;
	  }
	}
      };

      enum NoiseType {
	NOISE_NONE,
	NOISE_FIXED,
	NOISE_REALISTIC
      };

      void printTime();

      //------------------------------------------------------------
      // Array details
      //------------------------------------------------------------

      void setTelescopeName(std::string telescope);
      std::string getTelescopeName();

      void setInstrumentName(std::string instrument);
      std::string getInstrumentName();

      //------------------------------------------------------------
      // Array location
      //------------------------------------------------------------

      void setArrayLocation(gcp::util::Lla lla);
      Lla getArrayLocation();

      //------------------------------------------------------------
      // Methods for setting observation parameters
      //------------------------------------------------------------

      void setObsHa(gcp::util::HourAngle startHa, 
		    gcp::util::HourAngle stopHa, 
		    gcp::util::HourAngle deltaHa);
      HourAngle getStartHa();
      HourAngle getStopHa();
      HourAngle getDeltaHa();

      void setSourceName(std::string name);
      std::string getSourceName();

      void setObsRa(gcp::util::HourAngle obsRa);
      HourAngle getObsRa();

      void setObsDec(gcp::util::Declination dec);
      Declination getObsDec();

      void setObsEquinox(double equinox);
      double getObsEquinox();

      //------------------------------------------------------------
      // Methods for setting noise parameters
      //------------------------------------------------------------

      void setNoiseType(std::string type);
      void setNoiseType(NoiseType type);
      void setFixedNoiseRms(double val, std::string units);

      void setAmbientTemperature(gcp::util::Temperature& tAmb);
      void setApertureEfficiency(gcp::util::Percent& apEff);
      void setOpacity(gcp::util::Percent& tau);

      // Re-seed the random number generator for noise generation

      void seed(int s);

      //------------------------------------------------------------
      // Methods for setting antenna parameters
      //------------------------------------------------------------

      virtual void setAntennaLocation(gcp::util::LengthTriplet enu, int index);
      virtual void setAntennaX(gcp::util::Length X, unsigned index);
      virtual void setAntennaY(gcp::util::Length Y, unsigned index);
      virtual void setAntennaZ(gcp::util::Length Z, unsigned index);

      void setAntennaType(gcp::util::Antenna::AntennaType type, int index=-1);
      void setAntennaTypeIfUnknown(gcp::util::Antenna::AntennaType type, int index=-1);
      void setAntennaDiameter(gcp::util::Length diameter, int index=-1);
      void setAntennaApertureEfficiency(gcp::util::Percent apEff, int index=-1);
      void checkAntennaLocations();

      std::vector<LengthTriplet> getAntennaXyz();

      //------------------------------------------------------------
      // Methods for setting frequency information
      //------------------------------------------------------------

      void setFrequencyInformation(std::vector<gcp::util::Frequency>& freqs,
				   std::vector<gcp::util::Frequency>& bws);

      std::vector<Frequency> getFrequencies();
      std::vector<Frequency> getBandwidths();

      unsigned getNumberOfFrequencies();

      void markAntennaLocationsAsReceived();

      //------------------------------------------------------------
      // Methods for setting sizes
      //------------------------------------------------------------

      void setNumberOfAntennas(unsigned nAnt);
      unsigned getNumberOfAntennas();

      void setNumberOfBaselines(unsigned nBase);
      unsigned getNumberOfBaselines();

      void setNumberOfStokesParameters(unsigned nStokes);
      unsigned getNumberOfStokesParameters();

      void setNumberOfTimestamps(unsigned nTimestamps);
      unsigned getNumberOfTimestamps();

      // Groups is primary information read from files -- number of
      // timestamps can only be inferred once the number of baselines
      // is known.

      void setNumberOfGroups(unsigned nGroups);
      unsigned getNumberOfGroups();

      //------------------------------------------------------------
      // Miscellaneous methods
      //------------------------------------------------------------

      // Return true if all antennas have locations

      bool antsHaveLocations();

      // Plot the array

      void plotAntennas();
      void plotAntennasEnu();
      void plotAntennasXyz();

      double getFixedNoiseRms(std::string units);
      double getFixedNoiseRms(gcp::util::Unit::Units units);

      void getNoiseRms(Flux& noiseRms,
		       HourAngle* ha, Declination* dec, Frequency* bandwidth, Time* dt,
		       Antenna* ant1, Antenna* ant2, PolarLengthVector* azel=0);

      void generateNoise(Flux& noiseRms, Flux& reNoise, Flux& imNoise, double& wt);

      void generateNoise(Flux& reNoise, Flux& imNoise, double& wt, 
			 HourAngle* ha=0, Declination* dec=0, Frequency* bandwidth=0, Time* dt=0,
			 Antenna* ant1=0, Antenna* ant2=0);

      // get() functions for use by external writers

      float getU(unsigned iFrame, unsigned iBaseline);
      float getV(unsigned iFrame, unsigned iBaseline);
      float getW(unsigned iFrame, unsigned iBaseline);
      double getJulianDate(unsigned iFrame);
      double getJulianDate(unsigned iFrame, unsigned iBaseline);
      float getRe(unsigned iFrame, unsigned iBaseline, unsigned iIf, unsigned iStokes);
      float getIm(unsigned iFrame, unsigned iBaseline, unsigned iIf, unsigned iStokes);
      float getWt(unsigned iFrame, unsigned iBaseline, unsigned iIf, unsigned iStokes);
      float getAipsBaselineIndex(unsigned iFrame, unsigned iBaseline);

      void initializeSimulationVisibilityArray();

      bool canSimulate();

      /**
       * Constructor.
       */
      ObsInfo();
      ObsInfo(const ObsInfo& obs);
      ObsInfo(ObsInfo& obs);

      /**
       * Destructor.
       */
      virtual ~ObsInfo();

      void operator=(const ObsInfo& obs); 
      void operator=(ObsInfo& obs); 

      void writeUvfFile(std::string fileName);

      std::string& name() {
	return name_;
      }

      //------------------------------------------------------------
      // Parsing interface
      //------------------------------------------------------------

      void addParameters();
      void setParameter(std::string name, std::string val, std::string units="");
      void checkLocationParameters();
      void checkHaParameters();
      void checkFreqParameters();

      void printMissingParameters(unsigned mask);

      void initializeBaselineIndices();

    public:

      //------------------------------------------------------------
      // The visibilities that comprise this observation
      //------------------------------------------------------------

      std::vector<Vis> visibilities_;

      //------------------------------------------------------------
      // A mask of informational parameters
      //------------------------------------------------------------

      unsigned infoMask_;

      //------------------------------------------------------------
      // Array parameters
      //------------------------------------------------------------

      // The location of the center of the array

      gcp::util::Lla lla_;

      // The name of the array

      std::string telescopeName_;

      // The name of the instrument used for this observation

      std::string instrumentName_;

      //------------------------------------------------------------
      // Frequency information
      //------------------------------------------------------------

      std::vector<Frequency> frequencies_;
      std::vector<Frequency> bandwidths_;

      //------------------------------------------------------------
      // The vector of antennas in the array
      //------------------------------------------------------------

      std::vector<gcp::util::Antenna> antennas_;
      std::map<unsigned, unsigned> groupToAipsBaselineIndexMap_;

      //------------------------------------------------------------
      // Observation parameters
      //------------------------------------------------------------

      gcp::util::HourAngle startHa_;
      gcp::util::HourAngle stopHa_;
      gcp::util::HourAngle deltaHa_;

      std::string sourceName_;
      HourAngle   obsRa_;
      Declination obsDec_;
      double      obsEquinox_;
      double      mjdMin_;
      double      mjdMax_;

      //------------------------------------------------------------
      // Noise parameters
      //------------------------------------------------------------

      NoiseType   noiseType_;
      double      startJd_;
      bool        simNoise_;
      GenericNormalization noiseRms_;
      Temperature tAmb_;
      Percent     tau_;

      //------------------------------------------------------------
      // Size parameters
      //------------------------------------------------------------

      unsigned nGroup_;     // The number of visibility groups in this observation
      unsigned nBaseline_;  // The number of baselines in this observation
      unsigned nFreq_;      // The number of frequencies in this observation
      unsigned nStokes_;    // The number of Stokes parameters in this observation

      //------------------------------------------------------------
      // Private methods of this class
      //------------------------------------------------------------

      void calculateStartJd();
      void simulateNoise(bool simNoise);

      std::map<gcp::util::Antenna*, unsigned> getAntMap();
      std::map<Antenna*, std::vector<Antenna*> > getAntVecMap();

      std::vector<gcp::util::Antenna> mergeAnts(ObsInfo& obs);
      bool existsInAntMap(Antenna* ant, std::map<gcp::util::Antenna*, unsigned>& antMap);
      Antenna* getAntMapKey(Antenna* ant, std::map<gcp::util::Antenna*, unsigned>& antMap);
      Antenna* getAntMapKey(Antenna* ant, std::map<Antenna*, std::vector<Antenna*> >& antMap);

    }; // End class ObsInfo

    std::ostream& operator<<(std::ostream& os, ObsInfo::Vis& vis);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_OBSINFO_H
