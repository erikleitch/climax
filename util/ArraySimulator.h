// $Id: ArraySimulator.h,v 1.2 2012/05/09 21:17:48 eml Exp $

#ifndef GCP_FFTUTIL_ARRAYSIMULATOR_H
#define GCP_FFTUTIL_ARRAYSIMULATOR_H

/**
 * @file ArraySimulator.h
 * 
 * Tagged: Thu Feb  3 13:57:05 PST 2011
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/09 21:17:48 $
 * 
 * @author tcsh: Erik Leitch
 */
#include "gcp/util/Geoid.h"

namespace gcp {
  namespace fftutil {

    class ArraySimulator {
    public:

      //------------------------------------------------------------
      // An enumeration of parameters needed for simulation
      //------------------------------------------------------------

      enum {
	SIMPAR_NONE              = 0x0,

	// Location parameters

	SIMPAR_LOCATION_ARRAY    = 0x1,
	SIMPAR_LOCATION_ANTENNA  = 0x2,
	SIMPAR_LOCATION_INFO     = SIMPAR_LOCATION_ARRAY | SIMPAR_LOCATION_ANTENNA,

	// Source parameters

	SIMPAR_SRC_RA            = 0x4,
	SIMPAR_SRC_DEC           = 0x8,
	SIMPAR_SRC_EQN           = 0x10,
	SIMPAR_SRC_INFO          = SIMPAR_SRC_RA | SIMPAR_SRC_DEC | SIMPAR_SRC_EQN,

	// Observation parameters

	SIMPAR_OBS_HA            = 0x20,
	SIMPAR_IMAGE             = 0x40,
	SIMPAR_OBS_INFO          = SIMPAR_OBS_HA | SIMPAR_IMAGE,

	// Noise parameters

	SIMPAR_NOISE_TFIXED      = 0x80,
	SIMPAR_NOISE_TAMB        = 0x100,
	SIMPAR_NOISE_TAU         = 0x200,
	SIMPAR_NOISE_BW          = 0x400,

	SIMPAR_NOISE_INFO        = SIMPAR_NOISE_TFIXED | SIMPAR_NOISE_TAMB | SIMPAR_NOISE_TAU | SIMPAR_NOISE_BW,

	// Define the set of all required parameters

	SIMPAR_REQUIRED          = SIMPAR_LOCATION_INFO | SIMPAR_SRC_INFO | SIMPAR_OBS_INFO,
      };

      enum NoiseType {
	NOISE_NONE,
	NOISE_FIXED,
	NOISE_REALISTIC
      };

      void setArrayLocation(gcp::util::Lla lla);

      //------------------------------------------------------------
      // Methods for setting observation parameters
      //------------------------------------------------------------

      void setObsHa(gcp::util::HourAngle startHa, 
		    gcp::util::HourAngle stopHa, 
		    gcp::util::HourAngle deltaHa);

      void setObsRa(gcp::util::HourAngle obsRa);
      void setObsDec(gcp::util::Declination dec);
      void setObsDec(gcp::util::Declination dec);
      void setEquinox(double equinox);

      //------------------------------------------------------------
      // Methods for setting noise parameters
      //------------------------------------------------------------

      void setNoiseType(NoiseType type);
      void setFixedNoiseRms(gcp::util::Flux noiseRms);

      void setTfixed(gcp::util::Temperature& tFixed);
      void setBandwidth(gcp::util::Frequency& bw);
      void setTambient(gcp::util::Temperature& tAmb);
      void setApertureEfficiency(gcp::util::Percent& apEff);
      void setOpacity(gcp::util::Percent& tau);

      //------------------------------------------------------------
      // Methods for setting antenna parameters
      //------------------------------------------------------------

      virtual void setNumberOfAntennas(unsigned nAnt);
      virtual void setAntennaLocation(gcp::util::LengthTriplet enu, int index);

      void setAntennaType(gcp::util::Antenna::AntennaType type, int index=-1);
      void setAntennaDiameter(gcp::util::Length diameter, int index=-1);
      void setAntennaApertureEfficiency(gcp::util::Percent apEff, int index=-1);
      void checkAntennaLocations();

      // Return true if all antennas have locations

      bool antsHaveLocations();

      /**
       * Constructor.
       */
      ArraySimulator();

      /**
       * Destructor.
       */
      virtual ~ArraySimulator();

    protected:

      //------------------------------------------------------------
      // A mask of parameters required for simulation
      //------------------------------------------------------------

      unsigned simulationParameterMask_;

      //------------------------------------------------------------
      // The location of the array
      //------------------------------------------------------------

      gcp::util::Lla lla_;

      //------------------------------------------------------------
      // The vector of antennas in the array
      //------------------------------------------------------------

      std::vector<gcp::util::Antenna> antennas_;

      //------------------------------------------------------------
      // Observation parameters
      //------------------------------------------------------------

      gcp::util::HourAngle startHa_;
      gcp::util::HourAngle stopHa_;
      gcp::util::HourAngle deltaHa_;

      //------------------------------------------------------------
      // Noise parameters
      //------------------------------------------------------------

      NoiseType noiseType_;
      double startJd_;

      bool                   simNoise_;
      gcp::util::Flux        noiseRms_;
      gcp::util::Temperature tFixed_;
      gcp::util::Temperature tAmb_;
      gcp::util::Percent     tau_;

      //------------------------------------------------------------
      // Private methods of this class
      //------------------------------------------------------------

      void reportMissingSimulationParameters();
      void calculateStartJd();

    }; // End class ArraySimulator

  } // End namespace fftutil
} // End namespace gcp



#endif // End #ifndef GCP_FFTUTIL_ARRAYSIMULATOR_H
