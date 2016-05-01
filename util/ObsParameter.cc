#include "gcp/util/Exception.h"
#include "gcp/util/ObsParameter.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ObsParameter::ObsParameter() {}

/**.......................................................................
 * Destructor.
 */
ObsParameter::~ObsParameter() {}

void ObsParameter::checkParameters(unsigned mask, unsigned check)
{
  std::ostringstream os;
  if(parametersAreMissing(mask, check)) {
    reportMissingParameters(mask, check, os);
    ThrowSimpleError("Not enough information has been provided: " << os.str());
  }
}

bool ObsParameter::parametersAreMissing(unsigned mask, unsigned check)
{
  if((mask & check) != check) {
    return true;
  }

  return false;
}

/**.......................................................................
 * Notify the user of any required parameters that are missing
 */
void ObsParameter::reportMissingParameters(unsigned mask, unsigned check, std::ostringstream& os)
{
  if(check & INFO_LOCATION_ARRAY) {
    if(!(mask & INFO_LOCATION_ARRAY)) {
      os << "\nNo array location has been specified: use setArrayLocation()\n";
    }
  }

  if(check & INFO_LOCATION_ANTENNA) {
    if(!(mask & INFO_LOCATION_ANTENNA)) {
      os << "\nAntenna locations have not been fully specified: use setAntennaLocation()\n";
    }
  }

  if(check & INFO_OBS_RA) {
    if(!(mask & INFO_OBS_RA)) {
      os << "\nSource RA has not been specified: use setObsRa()\n";
    }
  }

  if(check & INFO_OBS_DEC) {
    if(!(mask & INFO_OBS_DEC)) {
      os << "\nSource DEC has not been specified: use setObsDec()\n";
    }
  }

  if(check & INFO_OBS_EQN) {
    if(!(mask & INFO_OBS_EQN)) {
      os << "\nSource EQN has not been specified: use setObsEquinox()\n";
    }
  }

  if(check & INFO_OBS_HA) {
    if(!(mask & INFO_OBS_HA)) {
      os << "\nStart/Stop/Delta HA of observation have not been specified: use setObsHa()\n";
    }
  }

  if(check & INFO_FREQ) {
    if(!(mask & INFO_FREQ)) {
      os << "\nNo Frequency array has been specified: use setFrequencyInformation()\n";
    }
  }

  if(check & INFO_NSTOKES) {
    if(!(mask & INFO_NSTOKES)) {
      os << "\nThe number of Stokes parameters has not been specified: use setNumberOfStokesParameters()\n";
    }
  }
}

void ObsParameter::printMissingParameters(unsigned mask, unsigned check)
{
  unsigned bitsNotSet = (mask & check) ^ check;

  std::ostringstream os;

  if(bitsNotSet & ObsParameter::INFO_LOCATION_ARRAY) 
    os << "Array reference location (latitude/longitude/altitude)" << std::endl;

  if(bitsNotSet & ObsParameter::INFO_LOCATION_ANTENNA) 
    os << "Antenna location (east/north/up or X/Y/Z)" << std::endl;

  if(bitsNotSet & ObsParameter::INFO_OBS_RA) 
    os << "Obs RA (obsRa)" << std::endl;

  if(bitsNotSet & ObsParameter::INFO_OBS_DEC) 
    os << "Obs DEC (obsDec)" << std::endl;

  if(bitsNotSet & ObsParameter::INFO_OBS_EQN) 
    os << "Obs Equinox (obsEquinox)" << std::endl;

  if(bitsNotSet & ObsParameter::INFO_SRC_NAME) 
    os << "Source name" << std::endl;

  if(bitsNotSet & ObsParameter::INFO_OBS_HA) 
    os << "Hour angle of observations (haStart/haStop/haDelta)" << std::endl;

  if(bitsNotSet & ObsParameter::INFO_NOISE_FIXED) 
    os << "Fixed noise rms (noiseRms)" << std::endl;

  if(bitsNotSet & ObsParameter::INFO_NOISE_TAMB) 
    os << "Ambient temperature (tambient)" << std::endl;

  if(bitsNotSet & ObsParameter::INFO_NOISE_TAU) 
    os << "Atmospheric opacity (tau)" << std::endl;

  if(bitsNotSet & ObsParameter::INFO_FREQ) 
    os << "Frequency information (freqs/bws)" << std::endl;

  if(bitsNotSet & ObsParameter::INFO_NANT) 
    os << "Number of antennas (nant)" << std::endl;

  if(bitsNotSet & ObsParameter::INFO_NSTOKES) 
    os << "Number of Stokes parameters (nstokes)" << std::endl;

  if(bitsNotSet & ObsParameter::INFO_NGROUP) 
    os << "Number of groups" << std::endl;

  if(bitsNotSet & ObsParameter::INFO_TELESCOPE) 
    os << "Telescope name (telescope)" << std::endl;

  if(bitsNotSet & ObsParameter::INFO_INSTRUMENT) 
    os << "Instrument name (instrument)" << std::endl;

  COUTCOLOR(std::endl << "The following information is missing: " << std::endl << std::endl << os.str() << std::endl, "red");
}

