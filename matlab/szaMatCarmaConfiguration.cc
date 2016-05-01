/**.......................................................................
 * MATLAB Mex file for calculating shadowing, given a specified GCP
 * and CARMA array
 *
 * Use like:
 *
 *    gcpMatShadowCarma(d)
 *
 * If d already contains fields array.carma.config
 *                                         .pad
 *                                         .nPad
 *                                         .antType
 *
 * then this function does nothing. If they do not exist, it creates
 * them and fills them with correct CARMA configuration information
 * for the timestamps given in array.frame.utc
 */
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

#include "gcp/array/code/share/slalib/slalib.h"

#include "gcp/util/CarmaConfig.h"
#include "gcp/util/Date.h"

#include "mex.h"
#include "matrix.h"

#include <iostream>
#include <iomanip>
#include <math.h>
#include <list>

using namespace std;
using namespace gcp::util;
using namespace gcp::matlab;

//=======================================================================
// Define a struct for handling a known CARMA array configuration
//=======================================================================

struct Config {

  Config(std::string conf, 
	 std::string startDate, std::string startTime, 
	 std::string stopDate,  std::string stopTime);

  Config(std::string conf, 
	 std::string startDate, std::string stopDate);

  void setTo(std::string conf, 
	     std::string startDate, std::string startTime, 
	     std::string stopDate,  std::string stopTime);
  
  void setTo(std::string conf, 
	     std::string startDate, std::string stopDate);

  friend std::ostream& operator<<(std::ostream& os, Config& config);

  Config(const Config& conf);
  Config(Config& conf);
  void operator=(const Config& conf);
  void operator=(Config& conf);

  // Convert the array specification to a CarmaConfig object

  void fillConfiguration();

  std::string config_;            // The nominal array configuration

  std::vector<unsigned> remPads_;     // An array of pads removed from the array
  std::vector<unsigned> addPads_;     // An array of pads added to the array
  std::vector<unsigned> addAntTypes_; // An array of antnena types
				      // added to the array

  Date startDate_;
  Date stopDate_;

  //------------------------------------------------------------
  // Utility members, derived from existing members of this object,
  // for efficiency.  These are set up when fillConfiguration() is
  // called
  //------------------------------------------------------------

  double startMjd_;
  double stopMjd_;

  CarmaConfig cc_; // The Configuration object that represents this
		   // configuration

  std::vector<CarmaConfig::PadLocation> pads_;
  ArrayConfig::Type type_;
};

Config::Config(std::string conf, 
	       std::string startDate, std::string startTime,
	       std::string stopDate, std::string stopTime)
{
  setTo(conf, startDate, startTime, stopDate, stopTime);
}

Config::Config(std::string conf, 
	       std::string startDate,
	       std::string stopDate)
{
  setTo(conf, startDate, stopDate);
}

Config::Config(const Config& conf)
{
  *this = (Config&) conf;
}

Config::Config(Config& conf)
{
  *this = conf;
}

void Config::operator=(const Config& conf)
{
  *this = (Config&)conf;
}

void Config::operator=(Config& conf)
{
  config_     = conf.config_;
  remPads_    = conf.remPads_;
  addPads_    = conf.addPads_;
  addAntTypes_= conf.addAntTypes_;
  startDate_  = conf.startDate_;
  stopDate_   = conf.stopDate_;
  startMjd_   = conf.startMjd_;
  stopMjd_    = conf.stopMjd_;
  cc_         = conf.cc_;
  pads_       = conf.pads_;
  type_       = conf.type_;
}

void Config::setTo(std::string conf, 
		   std::string startDate, std::string startTime,
		   std::string stopDate, std::string stopTime)
{
  config_ = conf;
  startDate_.setTo(startDate, startTime);
  stopDate_.setTo(stopDate, stopTime);
}

void Config::setTo(std::string conf, 
		   std::string startDate, std::string stopDate)
{
  config_ = conf;
  startDate_.setTo(startDate);
  stopDate_.setTo(stopDate);
}

void Config::fillConfiguration()
{
  // Set to the nominal configuration

  cc_.setCurrentConfiguration(config_);

  // Remove any pads that were specified

  for(unsigned i=0; i < remPads_.size(); i++) {
    cc_.removePad(remPads_[i]);
  }

  // Add any pads that were specified

  for(unsigned i=0; i < addPads_.size(); i++) {
    cc_.addPad(addPads_[i], addAntTypes_[i]);
  }

  pads_ = cc_.getCurrentConfiguration();

  type_ = cc_.confType();

  startMjd_ = startDate_.mjd();
  stopMjd_  = stopDate_.mjd();
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
operator<<(ostream& os, Config& config)
{
  os << setw(2) << config.config_  << " " << setw(2) << config.type_ 
     << " start = " << config.startDate_ 
     << " (" << std::setprecision(12) << config.startMjd_ << ")"
     << " stop  = " << config.stopDate_ 
     << " (" << std::setprecision(12) << config.stopMjd_ << ")";

  return os;
}

//=======================================================================
// Forward declarations of global functions used in this mex file
//=======================================================================

void fillConfig(unsigned nObs, unsigned iObs, 
		std::list<Config>::iterator& currConfig, 
		unsigned char* carmaConfPtr, 
		unsigned int* carmaPadNumberPtr, 
		unsigned int* carmaNpadPtr=0, 
		unsigned int* carmaAntTypePtr=0);

void checkConfig(unsigned nFrame, unsigned iFrame, bool reqFieldsExist,
		 std::list<Config>::iterator& currConfig, 
		 std::list<Config>::iterator& endConfig, 
		 double* mjdPtr, 
		 unsigned char* confPtr, 
		 unsigned int* padNumberPtr, 
		 unsigned int* nPadPtr=0, 
		 unsigned int* antTypePtr=0);

std::list<Config> initializeKnownCarmaConfigs();
std::list<Config> initializeKnownSzaConfigs();

bool lessThanPred(Config& c1, Config& c2);

//=======================================================================
// Main
//=======================================================================

/**.......................................................................
 * Entry point from the matlab environment
 */
void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  // Redirect error handling

  gcp::util::Logger::installStdoutPrintFn(MexHandler::stdoutPrintFn);
  gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);
  
  gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
  gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
  gcp::util::ErrHandler::installLogFn(MexHandler::logFn);

  // Create Usage message to be displayed on input error

  std::ostringstream usage;
  
  usage << "Usage: " 
	<< std::endl
	<< std::endl
	<< "   gcpMatCarmaConfiguration(d) " << std::endl
	<< std::endl
	<< " If d already contains register fields array.carma.*, "
	<< "this function does nothing" 
	<< std::endl << std::endl
	<< " If it does not, it creates them and attempts to fill them with "
	<< "known CARMA array configuration information"
	<< std::endl << std::endl
	<< " Note that d must contain at least the register: "
	<< std::endl << std::endl
	<< "      array.frame.utc "
	<< std::endl << std::endl
	<< "  formatted as a double mjd";

  // Check arguments

  if(nrhs != 1) {
    ThrowError(usage.str());
  }

  // Check that required members are present

  MexParser base(prhs[0]);

  if(!base.fieldExists("array.frame.utc")) {
    ThrowError(usage.str());
  }

  MexParser mjd(base.getField("array.frame.utc"));

  double* mjdPtr  = mjd.getDoubleData();
  unsigned nFrame = mjd.getNumberOfElements();

  // Get the array of know configurations

  std::list<Config> knownCarmaConfigs = initializeKnownCarmaConfigs();
  std::list<Config> knownSzaConfigs   = initializeKnownSzaConfigs();

  COUT(std::endl << "Known CARMA configurations are: " << std::endl);

  for(std::list<Config>::iterator conf=knownCarmaConfigs.begin(); 
      conf != knownCarmaConfigs.end(); conf++) {
    COUT(*conf);
  }

  COUT(std::endl << "Known GCP configurations are: " << std::endl);

  for(std::list<Config>::iterator conf=knownSzaConfigs.begin(); 
      conf != knownSzaConfigs.end(); conf++) {
    COUT(*conf);
  }

  //-----------------------------------------------------------------------
  // If the carma configuration arrays aren't present in the input
  // structure, add them now
  //-----------------------------------------------------------------------

  bool reqFieldsExist = true;

  if(!base.fieldExists("array.carma.config")) {
    reqFieldsExist = false;
    MexHandler::addRegisterField((mxArray*)prhs[0],"array.carma.config", nFrame);
  }

  if(!base.fieldExists("array.carma.pad")) {
    reqFieldsExist = false;
    MexHandler::addRegisterField((mxArray*)prhs[0],"array.carma.pad",    nFrame);
  }

  if(!base.fieldExists("array.carma.nPad")) {
    reqFieldsExist = false;
    MexHandler::addRegisterField((mxArray*)prhs[0],"array.carma.nPad",   nFrame);
  }

  if(!base.fieldExists("array.carma.antType")) {
    reqFieldsExist = false;
    MexHandler::addRegisterField((mxArray*)prhs[0],"array.carma.antType",nFrame);
  }

  //-----------------------------------------------------------------------
  // If the gcp configuration arrays aren't present in the input
  // structure, add them now
  //-----------------------------------------------------------------------

  if(!base.fieldExists("array.gcp.config")) {
    reqFieldsExist = false;
    MexHandler::addRegisterField((mxArray*)prhs[0],"array.gcp.config", nFrame);
  }

  if(!base.fieldExists("array.gcp.pad")) {
    reqFieldsExist = false;
    MexHandler::addRegisterField((mxArray*)prhs[0],"array.gcp.pad",    nFrame);
  }

  // If the required fields didn't exist, or the mjd spans a time for
  // which we have explicit configuration information, fill in the
  // array configuration info now

  double startMjd = *(mjdPtr);
  double stopMjd  = *(mjdPtr + nFrame);

  Config& lastCarmaConfig = knownCarmaConfigs.back();
  Config& lastSzaConfig   = knownSzaConfigs.back();

  if(!reqFieldsExist || 
     startMjd < lastSzaConfig.stopDate_.mjd() || 
     startMjd < lastCarmaConfig.stopDate_.mjd()) {

    // Now iterate through the mjds, setting up the appropriate array
    // for each timestamp
    
    unsigned char* carmaConfPtr      = base.getFieldAsUchar("array.carma.config");
    unsigned int*  carmaPadNumberPtr = base.getFieldAsUint( "array.carma.pad");
    unsigned int*  carmaNpadPtr      = base.getFieldAsUint( "array.carma.nPad");
    unsigned int*  carmaAntTypePtr   = base.getFieldAsUint( "array.carma.antType");

    unsigned char* gcpConfPtr        = base.getFieldAsUchar("array.gcp.config");
    unsigned int*  gcpPadNumberPtr   = base.getFieldAsUint( "array.gcp.pad");
    
    std::list<Config>::iterator currCarmaConfig = knownCarmaConfigs.begin();
    std::list<Config>::iterator currSzaConfig   = knownSzaConfigs.begin();
    std::list<Config>::iterator endCarmaConfig  = knownCarmaConfigs.end();
    std::list<Config>::iterator endSzaConfig    = knownSzaConfigs.end();
    
    // Now iterate over timeslots, filling the array configuration with
    // appropriate data for that timeslot, where those data are known
    
    for(unsigned iFrame=0; iFrame < nFrame; iFrame++) {

      checkConfig(nFrame, iFrame, reqFieldsExist, currCarmaConfig, endCarmaConfig, 
		  mjdPtr, carmaConfPtr, carmaPadNumberPtr, carmaNpadPtr, carmaAntTypePtr);

      checkConfig(nFrame, iFrame, reqFieldsExist, currSzaConfig,   endSzaConfig,   
		  mjdPtr, gcpConfPtr,   gcpPadNumberPtr);
    }
  }

  return;
}

/**.......................................................................
 * Check for known configurations matching the current date
 */
void checkConfig(unsigned nFrame, unsigned iFrame, bool reqFieldExist,
		 std::list<Config>::iterator& currConfig, 
		 std::list<Config>::iterator& endConfig, 
		 double* mjdPtr, unsigned char* confPtr, unsigned int* padNumberPtr, unsigned int* nPadPtr, 
		 unsigned int* antTypePtr)
{
  double currMjd = *(mjdPtr+iFrame);
      
  if(currMjd < currConfig->startMjd_) {

    // If the mjd is before the first known CARMA configuration,
    // then this data was taken before the move to the high site

    *(confPtr+iFrame) = ArrayConfig::NOCARMA;

  } else if(currMjd <= currConfig->stopMjd_) {

    // Else if the mjd is before the end time of the current known
    // configuration, fill the pointers with appropriate
    // configuration data

    fillConfig(nFrame, iFrame, currConfig, 
	       confPtr, padNumberPtr, 
	       nPadPtr, antTypePtr);

  } else if(currConfig == endConfig) {

    // If we are already on the last configuration that we know about,
    // then mark the configuration as unknown.  But only if the
    // configuration fields didn't already exist.

    if(!reqFieldExist) {
      COUT("Setting to UNKNOWN");
      *(confPtr+iFrame) = ArrayConfig::UNKNOWN;
    }

  } else {

    // Else try to find the next bracketing configuration

    do {
      ++currConfig;
    } while(currConfig != endConfig && currMjd > currConfig->stopMjd_);

    if(currMjd <= currConfig->stopMjd_) {
	  
      fillConfig(nFrame, iFrame, currConfig, 
		 confPtr, padNumberPtr, 
		 nPadPtr, antTypePtr);
    } else {

      if(!reqFieldExist) {
	*(confPtr+iFrame) = ArrayConfig::UNKNOWN;
      }

    }
  }
}

//=======================================================================
// Global functions used in this mex file
//=======================================================================

/**.......................................................................
 * Initialize known CARMA configurations
 */
std::list<Config> initializeKnownCarmaConfigs()
{
  std::list<Config> configs;

  // D array lasted until 07 Sep 2008

  configs.push_back(Config("DO", "15 Jul 2008", "08 Sep 2008"));
  
  // During 08-09 Sep 2008, CARMA was in an intermediate state

  configs.push_back(Config("DO", "08 Sep 2008", "10 Sep 2008"));

  // E array lasted from 10 Sep to 07 Oct 2008

  configs.push_back(Config("E", "10 Sep 2008", "08 Oct 2008"));

  // During 08-10 Oct 2008, CARMA was in an intermediate state

  configs.push_back(Config("E", "08 Oct 2008", "10 Oct 2008"));

  // From 10 Oct 2008 to 17 Nov 2008 CARMA was in C array

  configs.push_back(Config("C", "10 Oct 2008", "17 Nov 2008"));

  // From 17 Nov 2008 to 12 Jan 2009, CARMA is in B array 

  configs.push_back(Config("B", "17 Nov 2008", "12 Jan 2009"));

  // From 12 Jan 2009 to 12 Feb 2009, CARMA is in A array 
  // Timing of this switch is uncertain at the moment. 

  configs.push_back(Config("A", "12 Jan 2009", "12 Feb 2009"));

  // From 17 Nov 2008 onward, CARMA is in D array 
  // Actually a hybrid configuration for a while

  configs.push_back(Config("D", "12 Feb 2009", "09 Apr 2009"));

  // As of 09 Apr, CARMA is in C array 

  configs.push_back(Config("C", "09 Apr 2009", "04 Jun 2009"));

  // CARMA began the move to E array on the 4th, finished on the 6th

  configs.push_back(Config("E", "04 Jun 2009", "20 Jul 2009"));

  // CARMA began the move to D array on July 20th, finished on the 23rd

  configs.push_back(Config("D", "20 Jul 2009", "07 Sep 2009"));

  // cm-run configuration for month of September  Fake as D

  configs.push_back(Config("D", "07 Sep 2009", "06 Oct 2009"));

  // CARMA began the move to C array on October 6th

  configs.push_back(Config("C", "06 Oct 2009", "30 Nov 2009"));

  // CARMA began the move to PACS B array on 30 Nov 2009

  configs.push_back(Config("B", "30 Nov 2009", "01 Jan 2010"));

  // CARMA began the move to PACS A array on 01 Jan 2010

  configs.push_back(Config("A", "01 Jan 2010", "24 Feb 2010"));

  // CARMA began the move to C array on 24 Feb 2010

  configs.push_back(Config("C", "24 Feb 2010", "07 Apr 2010"));

  // CARMA began the move to D array on 07 Apr 2010

  configs.push_back(Config("D", "07 Apr 2010", "01 Jun 2010"));

  // Sort the configurations by starting date

  configs.sort(lessThanPred);

  // Now fill the CarmaConfig containers

  for(std::list<Config>::iterator conf=configs.begin(); conf != configs.end();
      conf++) {
    conf->fillConfiguration();
  }

  return configs;
}

/**.......................................................................
 * Initialize known GCP configurations
 */
std::list<Config> initializeKnownSzaConfigs()
{
  std::list<Config> configs;

  // I array lasted util 20 November.  Started moving GCP telescopes
  // on Wednesday, 18 Nov 2008.

  configs.push_back(Config("I", "15 Jul 2008", "20 Nov 2008"));
  
  // From 20 Nov 2008 to 06 Jan 2009, GCP is in the B-array paired antenna
  // configuration (BP)

  configs.push_back(Config("BP", "20 Nov 2008", "06 Jan 2009"));

  // From 07 Jan 2009 to 17 Feb 2009, GCP is in the A-array paired antenna
  // configuration (AP)

  configs.push_back(Config("AP", "07 Jan 2008", "17 Feb 2009"));

  // From 17 Feb 2009 onward, GCP is in L-array 

  configs.push_back(Config("L", "17 Feb 2009", "03 Dec 2009"));

  // From 3 Dec 2009, GCP is in PACS B array

  configs.push_back(Config("BP", "03 Dec 2009", "01 Jan 2010"));

  // GCP began the move to PACS A array on 01 Jan 2010

  configs.push_back(Config("AP", "01 Jan 2010", "24 Feb 2010"));

  // GCP began the move back to L array on 24 Feb 2010

  configs.push_back(Config("L", "24 Feb 2010", "25 Apr 2010"));

  // Now fill the CarmaConfig containers

  for(std::list<Config>::iterator conf=configs.begin(); conf != configs.end();
      conf++) {
    conf->fillConfiguration();
  }

  return configs;
}

/**.......................................................................
 * A predicate for sorting a list of configurations
 */
bool lessThanPred(Config& c1, Config& c2)
{
  return c1.startDate_.mjd() < c2.startDate_.mjd();
}

/**.......................................................................
 * Fill a timeslot with array information from the array of known
 * configurations
 */
void fillConfig(unsigned nObs, unsigned iObs, 
		std::list<Config>::iterator& currConfig, 
		unsigned char* confPtr, 
		unsigned int*  padNumberPtr, 
		unsigned int*  nPadPtr, 
		unsigned int*  antTypePtr)
{
  unsigned npad = currConfig->pads_.size();

  *(confPtr + iObs) = currConfig->type_;

  if(nPadPtr) {
    *(nPadPtr + iObs) = npad;
  }

  // Iterate over pads, storing the pad numbers and antenna types

  for(unsigned iPad=0; iPad < npad; iPad++) {

    *(padNumberPtr + iPad * nObs + iObs) = 
      currConfig->pads_[iPad].padNumber_;

    if(antTypePtr) {
      *(antTypePtr   + iPad * nObs + iObs) = 
	currConfig->pads_[iPad].ant_.getAntType();
    }

  }
}

