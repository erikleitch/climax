// $Id: VisIo.h,v 1.3 2012/05/09 21:17:48 eml Exp $

#ifndef GCP_UTIL_VISIO_H
#define GCP_UTIL_VISIO_H

/**
 * @file VisIo.h
 * 
 * Tagged: Mon Oct  3 15:32:27 PDT 2005
 * 
 * @version: $Revision: 1.3 $, $Date: 2012/05/09 21:17:48 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Angle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Delay.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/Geoid.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Length.h"
#include "gcp/util/Temperature.h"
#include "gcp/util/Pressure.h"

#include <vector>

namespace gcp {
  namespace util {

    class VisIo {
    public:

      /**
       * Constructor.
       */
      VisIo();

      /**
       * Destructor.
       */
      virtual ~VisIo();

      //------------------------------------------------------------
      // A struct for managing UV parameters
      //------------------------------------------------------------

      struct UvPar {
	bool isSet_;
	bool hasChangedSinceLastWrite_;
	
	UvPar() {
	  isSet_ = false;
	  hasChangedSinceLastWrite_ = false;
	}

	bool isSet() {
	  return isSet_;
	}

	bool changed() {
	  return hasChangedSinceLastWrite_;
	}

	void initialize() {
	  isSet_ = false;
	  hasChangedSinceLastWrite_ = false;
	};

	void update() {
	  isSet_ = true;
	  hasChangedSinceLastWrite_ = true;
	};

	void markAsWritten() {
	  hasChangedSinceLastWrite_ = false;
	};

      };

      //------------------------------------------------------------
      // Methods to set information needed to write the FITS file
      //------------------------------------------------------------

      void setSourceName(std::string srcName);
      std::string getSourceName();

      void setRa(HourAngle ra);
      HourAngle getRa();

      void setDec(Declination dec);
      Declination getDec();

      void setEquinox(double equinox);
      double getEquinox();

      void setRaApp(HourAngle ra);
      HourAngle getRaApp();

      void setDecApp(Declination dec);
      Declination getDecApp();

      void setRaRef(HourAngle ra);
      void setDecRef(Declination dec);
      void setDRaApp(HourAngle ra);
      void setDDecApp(Declination dec);

      //------------------------------------------------------------
      // Methods to set information about the array
      //------------------------------------------------------------

      void setInstrumentName(std::string instrument);
      std::string getInstrumentName();

      void setTelescopeName(std::string telescope);
      std::string getTelescopeName();

      void setLatitude(Angle& lat);
      Angle getLatitude();

      void setLongitude(Angle& lng);
      Angle getLongitude();

      //------------------------------------------------------------
      // Methods to set information about the telescopes
      //------------------------------------------------------------

      void setNumberOfTelescopes(unsigned nTel);
      unsigned getNumberOfTelescopes();

      void setFirstTelescopeNum(unsigned iTel);

      // A method to check consistency with the current number of
      // telescopes

      void checkNumberOfTelescopes(unsigned nTel);

      // Install a set of telescope locations
      
      void setTelescopeLocations(std::vector<LengthTriplet>& locations);
      std::vector<LengthTriplet> getTelescopeLocations();

      // Set the telescope diameters

      void setTelescopeDiameters(std::vector<Length>& diameters);
      void setTelescopeDiameter(Length& diameter);

      void setTelescopeAzimuth(std::vector<Angle>& az);
      void setTelescopeElevation(std::vector<Angle>& el);

      // Set the telescope aperture efficiencies

      void setTelescopeApertureEfficiencies(std::vector<float>& apeffs);
      void setTelescopeApertureEfficiency(float apeff);

      // Set the integration time

      void setIntTime(Time& time);

      //------------------------------------------------------------
      // Frequency information
      //------------------------------------------------------------

      // Install a set of IF frequencies

      void setNumberOfIfs(unsigned nIf);
      void setIfFrequencies(std::vector<Frequency>& frequencies);
      void setDeltaIfFrequencies(std::vector<Frequency>& frequencies);

      // An alternate method for specifying IF frequencies is to use a
      // starting frequency and a delta

      void setStartingIfFrequency(Frequency frequency);
      void setDeltaIfFrequency(Frequency& frequency);

      void setNumberOfChannelsPerIf(unsigned nChan);
      void setDeltaChannelFrequency(Frequency frequency);

      //------------------------------------------------------------
      // Methods to set weather information
      //------------------------------------------------------------

      void setAirTemperature(Temperature& temp);
      void setWindDirection(Angle& windDirection);
      void setWindSpeed(Speed& windSpeed);
      void setPressure(Pressure& pressure);
      void setRelativeHumidity(double relativeHumidity);

      //------------------------------------------------------------
      // Methods to set linelength information
      //------------------------------------------------------------

      void setLinelengthPhases(std::vector<Angle>& phase);
      void setLinelengthDelays(std::vector<Delay>& delay);

      //------------------------------------------------------------
      // Methods to set visibility data
      //------------------------------------------------------------

      void setUvw(double* uvw);
      void setVisWide(double* re, double* im);
      void setVisSpec(double* re, double* im);
      void setVisFlags(bool* visFlags);
      void setRms(double* rms);
      void setMjd(double* mjd);
      void setLst(double* lst);

      //-----------------------------------------------------------------------
      // Internal bookkeeping
      //-----------------------------------------------------------------------

      void setNumberOfFrames(unsigned nFrame) {
	nFrame_ = nFrame;
      }

      void setNumberOfBaselines(unsigned nBaseline) {
	nBaseline_ = nBaseline;
      }

      unsigned getNumberOfBaselines() {
        return nBaseline_;
      }
 
      void setBaselines(unsigned *baselines) {
        // Optional list of baseline numbers for
        // data sets that do not include a full set.
        baselines_ = baselines;
        baselinesPar_.update();
      }

      void setNumberOfStokesParameters(unsigned nStokes) {
	nStokes_ = nStokes;
      }

      //------------------------------------------------------------
      // Methods to write the file
      //------------------------------------------------------------

      virtual void checkParameters();
      virtual void openFile(std::string fileName) {};
      virtual void closeFile() {};

      void setPurpose(char purpose);      
      void setPurpose(std::string purpose);      

      //------------------------------------------------------------
      // Utility methods
      //------------------------------------------------------------

      // Return the antenna indices associated with a given visibility
      // index.

      void getTelescopeIndices(unsigned baslineIndex, 
			       unsigned* iRow, unsigned* iCol, 
			       unsigned nTel);

      float jyPerK(Length& diameter, float apeff);
      float jyPerK();

      void conjugateBaselines(bool conj);

    protected:

      bool doConj_;
      float jyPerK_;

      Time intTime_;
      UvPar intTimePar_;

      //------------------------------------------------------------
      // Information about the array
      //------------------------------------------------------------

      Angle latitude_;
      Angle longitude_;
      std::string telescope_;
      std::string instrument_;
      std::string date_;

      UvPar latitudePar_;
      UvPar longitudePar_;
      UvPar telescopePar_;
      UvPar instrumentPar_;

      //------------------------------------------------------------
      // Information about the source
      //------------------------------------------------------------

      std::string srcName_;

      // Mean coordinates

      HourAngle ra_;
      Declination dec_;

      // Equinox

      double equinox_;

      // Apparent

      HourAngle raApp_;
      Declination  decApp_;

      // Offsets

      HourAngle dRaApp_;
      Declination dDecApp_;

      HourAngle raRef_;
      Declination decRef_;

      UvPar srcNamePar_;

      UvPar raPar_;
      UvPar decPar_;

      UvPar equinoxPar_;

      UvPar raAppPar_;
      UvPar decAppPar_;

      UvPar raRefPar_;
      UvPar decRefPar_;

      UvPar dRaAppPar_;
      UvPar dDecAppPar_;

      //------------------------------------------------------------
      // Information about the frequency
      //------------------------------------------------------------

      // The number of IFs

      unsigned nIf_;
      Frequency startingIfFrequency_;
      Frequency deltaIfFrequency_;
      std::vector<Frequency> ifFrequencies_;
      std::vector<Frequency> deltaIfFrequencies_;
      Frequency ifCenterFrequency_;

      UvPar nIfPar_;
      UvPar ifFreqPar_;
      UvPar deltaIfFreqPar_;

      // The number of channels per IF

      unsigned nChannel_;
      Frequency startingFrequency_;
      Frequency deltaChannelFrequency_;

      UvPar nChannelPar_;
      UvPar deltaChannelFrequencyPar_;
      
      //------------------------------------------------------------
      // Information about telescopes
      //------------------------------------------------------------

      // The number of telescopes in our antenna table

      unsigned nTelescope_;      
      UvPar nTelescopePar_;

      std::vector<LengthTriplet> locations_;

      std::vector<Angle> az_;
      std::vector<Angle> el_;

      Length diameter_;
      std::vector<Length> diameters_;

      float apeff_;
      std::vector<float> apeffs_;

      UvPar locationsPar_;
      UvPar azPar_;
      UvPar elPar_;
      UvPar diameterPar_;
      UvPar apEffPar_;

      //------------------------------------------------------------
      // Weather information
      //------------------------------------------------------------

      Temperature airTemperature_;
      Angle windDirection_;
      Speed windSpeed_;
      float relativeHumidity_;
      Pressure pressure_;

      UvPar airTemperaturePar_;
      UvPar windDirectionPar_;
      UvPar windSpeedPar_;
      UvPar relativeHumidityPar_;
      UvPar pressurePar_;

      //------------------------------------------------------------
      //  Linelength information
      //------------------------------------------------------------

      std::vector<Angle> llPhases_;
      std::vector<Delay> llDelays_;

      UvPar llPhasesPar_;
      UvPar llDelaysPar_;

      //------------------------------------------------------------
      // Visibility data
      //------------------------------------------------------------

      double* uvw_;
      double* visWideRe_;
      double* visWideIm_;
      double* visSpecRe_;
      double* visSpecIm_;
      double* rms_;
      double* mjd_;
      double* lst_;
      bool* visFlags_;

      UvPar uvwPar_;
      UvPar visWidePar_;
      UvPar visSpecPar_;
      UvPar rmsPar_;
      UvPar mjdPar_;
      UvPar lstPar_;
      UvPar visFlagsPar_;

      //------------------------------------------------------------
      // Information about the observation
      //------------------------------------------------------------

      std::string purpose_;
      UvPar purposePar_;

      //------------------------------------------------------------
      // Internal bookkeeping
      //------------------------------------------------------------

      // The number of frames

      unsigned nFrame_;

      // The number of baselines

      unsigned nBaseline_;
      unsigned* baselines_;
      UvPar     baselinesPar_;
      unsigned  firstTelescopeNum_;

      // The number of Stokes parameters

      unsigned nStokes_;
      

    }; // End class VisIo

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_VISIO_H
