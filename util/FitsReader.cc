#include "gcp/util/FitsReader.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
FitsReader::FitsReader() 
{
  status_    = 0;
  fptr_      = 0;
  wasOpened_ = false;
  nHdu_      = 0;
  hduOffset_ = 0;
  nTable_    = 0;
  infoMask_  = 0x0;
}

/**.......................................................................
 * Destructor.
 */
FitsReader::~FitsReader() 
{
  FitsReader::close();
}

void FitsReader::close()
{
  // Only close the file if it was opened by this object
  
  if(wasOpened_) {
    if(fptr_) {
      //      if(ffclos(fptr_, &status_) > 0)
      //	ReportFitsError("Error closing file");
      ffclos(fptr_, &status_); // Ignore error status
      fptr_ = 0;
    }
  }
}

void FitsReader::open(std::string fileName) 
{
  fileName_ = fileName;

  if(ffopen(&fptr_, fileName.c_str(), 0, &status_) > 0) {
    ThrowSimpleFitsError("Unable to open file: " << fileName);
  }

  wasOpened_ = true;

  // Now initialize some info about the file

  initializeInfo();
}

void FitsReader::initializeInfo()
{
  countHdus();
  countTables();
}

void FitsReader::countHdus()
{
  if(ffthdu(fptr_, &nHdu_, &status_))
    ThrowFitsError("Error getting nHdu");
}

void FitsReader::countTables()
{
  int hdutype;
  char extname[80];
  char* cptr = extname;
  int nfound;

  for(unsigned iHdu=0; iHdu < nHdu_-1; iHdu++) {

    if(ffmrhd(fptr_, 1, &hdutype, &status_) > 0)
      ThrowFitsError("Problem moving to next HDU");
    
    // Only count tables
   
    if(hdutype == 2 || hdutype == 1) {
      ++nTable_;
    }

    // Read the extname of this HDU

    if(fits_read_keys_str(fptr_, "EXTNAME", 0, 1, &cptr, &nfound, &status_))
      ThrowFitsError("EXTNAME keyword not found in table header");
  }
      
  if(ffmrhd(fptr_, -(nHdu_-1), &hdutype, &status_) > 0)
    ThrowFitsError("Problem moving to the main HDU: " << -nHdu_);
}

bool FitsReader::hasTable(std::string tableName)
{
  int hdutype;
  char extname[80];
  char* cptr = extname;
  int nfound;
  bool found = false;
  int nHdu = 0;

  for(unsigned iHdu=0; iHdu < nHdu_-1; iHdu++) {

    if(ffmrhd(fptr_, 1, &hdutype, &status_) > 0)
      ThrowFitsError("Problem moving to next HDU");
    
    // Read the extname of this HDU

    if(fits_read_keys_str(fptr_, "EXTNAME", 0, 1, &cptr, &nfound, &status_))
      ThrowFitsError("EXTNAME keyword not found in table header");

    ++nHdu;

    if(extname == tableName) {
      found = true;
      break;
    }
  }
      
  if(ffmrhd(fptr_, -(nHdu-1), &hdutype, &status_) > 0)
    ThrowFitsError("Problem moving to the main HDU: " << -nHdu_);

  return found;
}

std::string FitsReader::getStringKey(std::string name)
{
  std::string keyword, comment;
  keyword.resize(80);
  comment.resize(80);

  if(ffgkys(fptr_, name.c_str(), &keyword[0], &comment[0], &status_))
    ThrowFitsError("Unable to get " << name << " keyword");

  return keyword;
}

long FitsReader::getLongKey(std::string name)
{
  std::string keyword, comment;
  comment.resize(80);
  long val;

  if(ffgkyj(fptr_, name.c_str(), &val, &comment[0], &status_) )
    ThrowFitsError("Unable to get " << name << " keyword");

  return val;
}

int FitsReader::getIntKey(std::string name)
{
  std::string keyword, comment;
  comment.resize(80);
  int val;

  if(ffgkyl(fptr_, name.c_str(), &val, &comment[0], &status_) )
    ThrowFitsError("Unable to get " << name << " keyword");

  return val;
}

/**.......................................................................
 * Move the the HDU offset by N from the current HDU
 */
void FitsReader::moveToMainHdu()
{
  moveToHduOffsetBy(-hduOffset_);
}

/**.......................................................................
 * Move the the HDU offset by N from the current HDU
 */
int FitsReader::moveToHduOffsetBy(unsigned offset) 
{
  if(hduOffset_+offset > nHdu()-1)
    ThrowError("Requested offset exceeds the number of HDUs (" << nHdu() << ") in this file");

  // Move forward by the requested number of hdus
 
  int hdutype;
  if(ffmrhd(fptr_, offset, &hdutype, &status_) > 0)
    ThrowFitsError("Problem moving to next HDU");
  
  // Increment the hdu offset
  
  hduOffset_ += offset;

  return hdutype;
}

/**.......................................................................
 * Method to move from the current location in the file to the next
 * HDU with the given extension name
 */
int FitsReader::moveForwardToTable(std::string extspec)
{
  // Iterate through the remaining HDUs, looking for a match

  unsigned nHdusLeft = nHdu()-hduOffset_-1;

  for(unsigned iHdu=0; iHdu < nHdusLeft; iHdu++) {

    int hdutype = moveToHduOffsetBy(1);

    // Check that this is a valid table

    if(hdutype == BINARY_TBL || hdutype == ASCII_TBL) {
      char extname[80];
      char* cptr = extname;
      int nfound;

      if(fits_read_keys_str(fptr_, "EXTNAME", 0, 1, &cptr, &nfound, &status_))
	ThrowFitsError("EXTNAME keyword not found in table header");

      if(strcmp(extname, extspec.c_str())==0) {
	tableName_ = extname;
	return hdutype;
      }

    }
  }

  ThrowError("Table: " << extspec << " not found");
  return 0;
}

/**.......................................................................
 * Method to move from the the first location in the file to the next
 * HDU with the given extension name
 */
int FitsReader::moveToTable(std::string extspec)
{
  moveToMainHdu();
  return moveForwardToTable(extspec);
}

/**.......................................................................
 * Move to the next valid table header
 */
int FitsReader::moveToNextTable()
{
  // Iterate through the remaining HDUs, looking for a match

  unsigned nHdusLeft = nHdu()-hduOffset_-1;

  for(unsigned iHdu=0; iHdu < nHdusLeft; iHdu++) {

    int hdutype = moveToHduOffsetBy(1);

    // Check that this is a valid table

    if(hdutype == BINARY_TBL || hdutype == ASCII_TBL) {

      char extname[80];
      char* cptr = extname;
      int nfound=0;

      if(fits_read_keys_str(fptr_, "EXTNAME", 0, 1, &cptr, &nfound, &status_))
	ThrowFitsError("EXTNAME keyword not found in table header");

      tableName_ = extname;
      return hdutype;
    }
  }

  ThrowError("No more tables found");
  return 0;
}

std::string FitsReader::getKeyword(std::string extname, std::string keyword)
{
  moveToTable(extname);

  char value[80];
  char* cptr = value;
  int nfound=0;
  
  if(fits_read_keys_str(fptr_, (char*)keyword.c_str(), 0, 1, &cptr, &nfound, &status_))
    ThrowFitsError("Keyword " << keyword << " not found in table header");
  
  return value;
}

std::string FitsReader::getKeyword(std::string keyword)
{
  moveToMainHdu();
  return getStringKey(keyword);
}

std::ostream& gcp::util::operator<<(std::ostream& os, const FitsReader::Axis& axis)
{
  os << "Axis: " << std::endl
     << "type        = " << axis.type_        << std::endl
     << "typeComment = " << axis.typeComment_ << std::endl
     << "unit        = " << axis.unit_        << std::endl
     << "n           = " << axis.n_           << std::endl
     << "refVal      = " << axis.refVal_      << std::endl
     << "refPix      = " << axis.refPix_      << std::endl
     << "delta       = " << axis.delta_       << std::endl << std::endl;

  return os;
}

/**.......................................................................
 * Read generic header information
 */
void FitsReader::getHeaderInfo()
{
  std::string comment;
  comment.resize(80);

  object_.resize(80);
  telescope_.resize(80);
  instrument_.resize(80);

  // Read the OBJECT keyword

  int nfound=0;

  if(ffgkys(fptr_, "OBJECT", &object_[0], &comment[0], &status_)) {
    ReportSimpleColorError("Unable to get OBJECT keyword", "yellow");
    status_ = 0;
  } else {
    infoMask_ |= ObsParameter::INFO_SRC_NAME;
  }

  if(ffgkys(fptr_, "TELESCOP", &telescope_[0], &comment[0], &status_)) {
    ReportSimpleColorError("Unable to get TELESCOP keyword", "yellow");
    status_ = 0;
  } else {
    infoMask_ |= ObsParameter::INFO_TELESCOPE;
  }

  if(ffgkys(fptr_, "INSTRUME", &instrument_[0], &comment[0], &status_)) {
    ReportSimpleColorError("Unable to get INSTRUME keyword", "yellow");
    status_ = 0;
  } else {
    infoMask_ |= ObsParameter::INFO_INSTRUMENT;
  }
  
  double ra;
  if(ffgkyd(fptr_, "OBSRA", &ra, &comment[0], &status_)) {
    ReportSimpleColorError("Unable to get OBSRA keyword", "yellow");
    status_ = 0;
  } else {
    obsRa_.setDegrees(ra);
    infoMask_ |= ObsParameter::INFO_OBS_RA;
  }

  double dec;
  if(ffgkyd(fptr_, "OBSDEC", &dec, &comment[0], &status_)) {
    ReportSimpleColorError("Unable to get OBSDEC keyword", "yellow");
    status_ = 0;
  } else {
    obsDec_.setDegrees(dec);
    infoMask_ |= ObsParameter::INFO_OBS_DEC;
  }

  dateObs_.resize(80);
  if(ffgkys(fptr_, "DATE-OBS", &dateObs_[0], &comment[0], &status_)) {
    ReportSimpleColorError("Warning: unable to get DATE-OBS keyword", "yellow");
    dateObs_ = "UNKNOWN";
    status_ = 0;
  } else {
    infoMask_ |= ObsParameter::INFO_DATE;
  }

  if(ffgkyd(fptr_, "EQUINOX", &equinox_, &comment[0], &status_)) {
    ReportSimpleColorError("Unable to get EQUINOX keyword: assuming these are J2000 coordinates" << std::endl, "yellow");
    equinox_ = 2000.0;
    status_ = 0;
  } else {
    infoMask_ |= ObsParameter::INFO_OBS_EQN;
  }
}

std::string FitsReader::telescope() {
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_TELESCOPE);
  return telescope_;
}
std::string FitsReader::instrument() {
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_INSTRUMENT);
  return instrument_;
}

std::string FitsReader::dateObs() {
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_DATE);
  return dateObs_;
}

std::string FitsReader::object() {
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_SRC_NAME);
  return object_;
}

HourAngle FitsReader::obsra() {
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_OBS_RA);
  return obsRa_;
}

Declination FitsReader::obsdec() {
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_OBS_DEC);
  return obsDec_;
}

double FitsReader::equinox() {
  ObsParameter::checkParameters(infoMask_, ObsParameter::INFO_OBS_EQN);
  return equinox_;
}

