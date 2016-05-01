// $Id: $

#ifndef GCP_UTIL_OBSPARAMETER_H
#define GCP_UTIL_OBSPARAMETER_H

/**
 * @file ObsParameter.h
 * 
 * Tagged: Tue Mar  3 14:03:28 PST 2015
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include <iostream>

namespace gcp {
  namespace util {

    class ObsParameter {
    public:

      /**
       * Constructor.
       */
      ObsParameter();

      /**
       * Destructor.
       */
      virtual ~ObsParameter();

      //------------------------------------------------------------
      // An enumeration of parameters needed for simulation
      //------------------------------------------------------------

      enum {

	INFO_NONE              = 0x0,

	// Location parameters

	INFO_LOCATION_ARRAY    = 0x1,
	INFO_LOCATION_ANTENNA  = 0x2,
	SIM_REQ_LOCATION_INFO  = INFO_LOCATION_ARRAY | INFO_LOCATION_ANTENNA,

	// Source parameters

	INFO_OBS_RA            = 0x4,
	INFO_OBS_DEC           = 0x8,
	INFO_OBS_EQN           = 0x10,
	INFO_SRC_NAME          = 0x20,
	SIM_REQ_SRC_INFO       = INFO_OBS_RA | INFO_OBS_DEC | INFO_OBS_EQN,

	// Observation parameters

	INFO_OBS_HA            = 0x40,
	SIM_REQ_OBS_INFO       = INFO_OBS_HA,

	// Noise parameters

	INFO_NOISE_FIXED       = 0x80,
	SIM_REQ_NOISE_FIXED    = INFO_NOISE_FIXED,

	INFO_NOISE_TAMB        = 0x100,
	INFO_NOISE_TAU         = 0x200,
	SIM_REQ_NOISE_REALISTIC= INFO_NOISE_TAU | INFO_NOISE_TAMB,

	// Frequency information

	INFO_FREQ              = 0x800,
	SIM_REQ_FREQ_INFO      = INFO_FREQ,

	// Size information

	INFO_NANT              = 0x1000,
	INFO_NSTOKES           = 0x2000,
	INFO_NGROUP            = 0x4000,

	// Other information

	INFO_TELESCOPE         = 0x8000,
	INFO_INSTRUMENT        = 0x10000,
	INFO_DATE              = 0x20000,

	// Define the set of all required parameters (except for noise)

	SIM_REQ_INFO           = SIM_REQ_LOCATION_INFO | SIM_REQ_SRC_INFO | SIM_REQ_OBS_INFO | SIM_REQ_FREQ_INFO
      };

      static bool parametersAreMissing(unsigned mask, unsigned check);
      static void checkParameters(unsigned mask, unsigned check);
      static void reportMissingParameters(unsigned mask, unsigned check, std::ostringstream& os);
      static void printMissingParameters(unsigned mask, unsigned check);

    }; // End class ObsParameter

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_OBSPARAMETER_H
