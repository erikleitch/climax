#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/PtSrcReader.h"
#include "gcp/util/RegExpParser.h"

#include <iomanip>

using namespace std;

using namespace gcp::util;

// Static declarations for this class

Flux PtSrcReader::minFlux_ = 
Flux(Flux::Jansky(), 0.0);

Flux PtSrcReader::maxFlux_ = 
Flux(Flux::Jansky(), 1e10);

/**.......................................................................
 * Constructor.
 */
PtSrcReader::PtSrcReader(std::string catalogFile) : catalogFile_(catalogFile) 
{
  initialize();
}

/**.......................................................................
 * Constructor.
 */
PtSrcReader::PtSrcReader() 
{
  initialize();
}

/**.......................................................................
 * Initialize critical members
 */
void PtSrcReader::initialize() 
{
  raMin_.setHours(0.0); 
  raMax_.setHours(24.0); 
}

/**.......................................................................
 * Destructor.
 */
PtSrcReader::~PtSrcReader() {}

/**.......................................................................
 * Set the catalog file
 */
void PtSrcReader::setCatalogFile(std::string catalogFile) 
{
  catalogFile_ = catalogFile;
}

/**.......................................................................
 * Check if the second position is within radius of the first position
 */
bool PtSrcReader::checkAngle(PtSrcReader::Source& src, HourAngle& ra, Declination& dec, Angle& radius)
{
  double sd1 = sin(dec.radians());
  double cd1 = cos(dec.radians());

  double sd2 = sin(src.dec_.radians());
  double cd2 = cos(src.dec_.radians());

  double sr1 = sin(ra.radians());
  double cr1 = cos(ra.radians());

  double sr2 = sin(src.ra_.radians());
  double cr2 = cos(src.ra_.radians());

  double rad = acos(cr1*cd1*cr2*cd2 + sr1*cd1*sr2*cd2 + sd1*sd2);

  src.distance_.setRadians(rad);
  
  return rad <= radius.radians();
}

/**.......................................................................
 * Utility function to allow: std::cout << src
 */
void PtSrcReader::printHeader(std::ostream& os)
{
  os << std::right << setw(13) << "Field    ";
  os << std::right << setw(13) << "RA(J2000)";
  os << std::right << setw(13) << "DEC(J2000)";
  os << std::right << setw(10) << "dist";

  os << std::right << setw(5)  << "warn";

  os << std::right << setprecision(5) << setw(10) << "Peak";
  os << std::right << setprecision(5) <<  setw(9) << "Peak Err";

  os << std::right << setprecision(5) <<  setw(9) << "Rms";

  os << std::right << setprecision(5) << setw(10) << "Int Flux";
  os << std::right << setprecision(5) <<  setw(9) << "Int Err";

  os << std::right << setprecision(5) << setw(10) << "Maj Axis";
  os << std::right << setprecision(5) <<  setw(9) << "Maj Err";

  os << std::right << setprecision(5) << setw(10) << "Min Axis";
  os << std::right << setprecision(5) <<  setw(9) << "Min Err";

  os << std::right << setprecision(5) << setw(10) << "PA";
  os << std::right << setprecision(5) <<  setw(9) << "PA Err";

  // And units

  os << std::endl;
  os << std::right << setw(13) << "(hh:mm:ss.s)";
  os << std::right << setw(13) << "(dd:mm:ss.s)";
  os << std::right << setw(10) << "(\")";

  os << std::right << setw(5)  << "";

  os << std::right << setprecision(5) << setw(10) << "(mJy/bm)";
  os << std::right << setprecision(5) <<  setw(9) << "(mJy/bm)";

  os << std::right << setprecision(5) <<  setw(9) << "(mJy/bm)";

  os << std::right << setprecision(5) << setw(10) << "(mJy)";
  os << std::right << setprecision(5) <<  setw(9) << "(mJy)";

  os << std::right << setprecision(5) << setw(10) << "(\")";
  os << std::right << setprecision(5) <<  setw(9) << "(\")";

  os << std::right << setprecision(5) << setw(10) << "(\")";
  os << std::right << setprecision(5) <<  setw(9) << "(\")";

  os << std::right << setprecision(5) << setw(10) << "(deg)";
  os << std::right << setprecision(5) <<  setw(9) << "(deg)";
}

/**.......................................................................
 * Utility function to allow: std::cout << src
 */
std::ostream& gcp::util::operator<<(std::ostream& os, PtSrcReader::Source& src)
{
  os << std::right << setw(4) << src.survey_;
  os << std::right << setw(6) << src.name_;

  os << std::right << setw(13) << src.ra_;
  os << std::right << setw(13) << src.dec_;
  os << std::right << std::fixed << setprecision(2) << setw(10) << src.distance_.arcsec();

  os << std::right  << setw(5)  << src.warn_;

  os << std::right << std::fixed << setprecision(2) << setw(10) << src.peak_.mJy();
  os << std::right << std::fixed << setprecision(2) <<  setw(9) << src.peakErr_.mJy();

  os << std::right << std::fixed << setprecision(2) <<  setw(9) << src.rms_.mJy();

  os << std::right << std::fixed << setprecision(2) << setw(10) << src.int_.mJy();
  os << std::right << std::fixed << setprecision(2) <<  setw(9) << src.intErr_.mJy();

  char szStr[8];

  sprintf(szStr, "%c%5.2f", (src.decMajErr_.arcsec() < 0.0 ? '<' : ' '), src.decMaj_.arcsec());
  os << std::right << setw(10) << szStr;

  sprintf(szStr, "%3.2f", src.decMajErr_.arcsec());
  os << std::right << setw(9) << ((src.decMajErr_.arcsec() < 0.0) ? "---" : szStr);

  sprintf(szStr, "%c%5.2f", (src.decMinErr_.arcsec() < 0.0 ? '<' : ' '), src.decMin_.arcsec());
  os << std::right << setw(10) << szStr;

  sprintf(szStr, "%3.2f", src.decMinErr_.arcsec());
  os << std::right << setw(9) << ((src.decMinErr_.arcsec() < 0.0) ? "---" : szStr);

  os << std::right << std::fixed << setprecision(2) << setw(10) << src.decPa_.degrees();

  sprintf(szStr, "%3.2f", src.decPaErr_.degrees());
  os << std::right << setw(9) << ((src.decMajErr_.arcsec() < 0.0 || src.decMinErr_.arcsec() < 0.0) ? "---" : szStr);

  return os;
}

/**.......................................................................
 * Return a list of sources within radius of the requested position.
 *
 * Note that the FIRST catalog is in J2000 coordinates, so that the
 * requested positions should be too.
 */
std::vector<PtSrcReader::Source> 
PtSrcReader::findSources(HourAngle ra, Declination dec, Angle radius, Flux fMin, Flux fMax, bool doPrint)
{
  static PtSrcReader::Source src;
  std::vector<PtSrcReader::Source> srcs;
  bool first=true;

  // Calculate the RA range we need to search

  setRaRange(ra, dec, radius);

  // Make sure the catalog file is open

  openCatalogFile();

  // Iterate, reading from the catalog file until the EOF is reached

  while(!eof()) {
    src = readNextEntry();

    // If the source flux is within range, check the angle, but don't
    // do otherwise, since this is expensive.

    if(src.peak_ >= fMin && src.peak_ <= fMax) {
      if(checkAngle(src, ra, dec, radius)) {

	applyCorrections(src);

	if(first) {
	  std::ostringstream os;

	  if(doPrint) {
	    printHeader(os);
	    COUT(os.str());
	  }

	  first = false;
	}

	if(doPrint)
	  COUT(src);

	srcs.push_back(src);
      }
    }
  }

  // And close the file

  closeCatalogFile();

  // Return the source array

  return srcs;
}

/**.......................................................................
 * Return a list of sources that match the passed regexp string
 */
std::vector<PtSrcReader::Source> 
PtSrcReader::findSources(std::string regExpStr, bool exact, bool caseSensitive)
{
  static PtSrcReader::Source src;
  std::vector<PtSrcReader::Source> srcs;
  srcs.resize(0);
  bool first=true;
  String str;
  String regExpStrCopy = regExpStr;

  if(!caseSensitive) {
    regExpStrCopy = regExpStrCopy.toLower();
  }
  
  // Make sure the catalog file is open

  openCatalogFile();

  // Iterate, reading from the catalog file until the EOF is reached

  while(!eof()) {
    src = readNextEntry();

    str = src.name_;
    
    if(!caseSensitive) {
      str = str.toLower();
    }

    if((exact && regExpStrCopy == str) || (!exact && str.contains(regExpStrCopy.str()))) {
      srcs.push_back(src);
    }
  }

  // And close the file

  closeCatalogFile();

  // Return the source array

  return srcs;
}

/**.......................................................................
 * Return a list of sources within radius of the requested position.
 *
 * Note that the FIRST catalog is in J2000 coordinates, so that the
 * requested positions should be too.
 */
unsigned 
PtSrcReader::countSources(HourAngle ra, Declination dec, Angle radius, Flux fMin, Flux fMax)
{
  static PtSrcReader::Source src;
  unsigned nSrc = 0;

  // Calculate the RA range we need to search

  setRaRange(ra, dec, radius);

  // Make sure the catalog file is open

  openCatalogFile();

  // Iterate, reading from the catalog file until the EOF is reached

  while(!eof()) {

    src = readNextEntry();

    // If the source flux is within range, check the angle, but don't
    // do otherwwise, since this is expensive.

    if(src.peak_ >= fMin && src.peak_ <= fMax) {

      if(checkAngle(src, ra, dec, radius)) {

	++nSrc;

	if(nSrc==1) {
	  std::ostringstream os;
	  printHeader(os);
	  COUT(os.str());
	}

	applyCorrections(src);

	COUT(src);
      }
    }
  }

  // And close the file

  closeCatalogFile();

  // Return the number of sources

  return nSrc;
}

/**.......................................................................
 * Calculate the minimal RA range we need to search
 */
void PtSrcReader::setRaRange(HourAngle& ra, Declination& dec, Angle& radius) 
{
  // Ra limits only apply if the search radius is less than pi -
  // 2*dec.  Otherwise the search radius completely encloses the full
  // RA range for all decs > dec.

  if(radius.radians() < M_PI - 2*dec.radians()) {

    double cdec = cos(dec.radians());
    double dRa  = fabs(acos(1.0 - (1.0-cos(radius.radians()))/(cdec*cdec)));

    // Set the new RA limits.  Note that the HourAngle class will wrap
    // these into 0-2pi, so that raMin can end up > raMax.  If this is
    // the case, we will have to search any catalog sorted in order of
    // increasing RA in two ranges.

    raMin_.setRadians(ra.radians() - dRa);
    raMax_.setRadians(ra.radians() + dRa);
  }
}

void PtSrcReader::applyCorrections(PtSrcReader::Source& src) {}

/**.......................................................................
 * Report an error generated by the cfitsio library
 */
void PtSrcReader::throwCfitsioError(int status)
{
  char status_str[FLEN_STATUS];
  
  if (status) {
    fits_get_errstatus(status, status_str);  /* get the error description */
    ThrowSimpleError(status_str);
  }
}

/**.......................................................................
 * Index the list of sources
 */
void
PtSrcReader::indexSources()
{
  static PtSrcReader::Source src;
  unsigned nSrc = 0;

  // Initialize the number of sources

  nSrc_ = 0;

  // Initialize the index array

  for(unsigned iRa=0; iRa < 25; iRa++)
    sourceIndices_[iRa] = 0;

  // Make sure the catalog file is open

  openCatalogFile();

  // Iterate, reading from the catalog file until the EOF is reached

  while(!eof()) {

    // Read the next line of the catalog file

    src = readNextEntry();

    // Increment the total number of sources

    ++nSrc_;

    // And increment whichever index this source belongs to

    ++sourceIndices_[(unsigned) src.ra_.hours() + 1];
  }

  // Finally, correct from number of sources in each bin to absolute
  // starting source number 

  sourceIndices_[0] = 1;
  for(unsigned iRa=1; iRa < 24; iRa++) {
    sourceIndices_[iRa] += sourceIndices_[iRa-1];
  }

  // Close the file

  closeCatalogFile();
}
