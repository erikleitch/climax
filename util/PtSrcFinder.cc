#include "gcp/util/PtSrcFinder.h"

#include "gcp/util/Angle.h"
#include "gcp/util/Constants.h"
#include "gcp/util/Cosmology.h"
#include "gcp/util/Declination.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/PtSrcReader.h"
#include "gcp/util/FirstFitsReader.h"
#include "gcp/util/NvssReader.h"
#include "gcp/util/Flux.h"
#include "gcp/util/Mass.h"

#include <string.h>

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
PtSrcFinder::PtSrcFinder() {}

/**.......................................................................
 * Destructor.
 */
PtSrcFinder::~PtSrcFinder() {}

void PtSrcFinder::findSources(std::string dir, std::string cat, 
			      std::string raStr, std::string decStr, std::string radStr,
			      std::vector<HourAngle>& ras, std::vector<Declination>& decs)
{
  HourAngle ra;
  Declination dec;

  // Parse input arguments

  ra.setHours(raStr);
  dec.setDegrees(decStr);
  Angle radius;

  radius.setDegrees(radStr);

  findSources(dir, cat, ra, dec, radius, ras, decs);
}

void PtSrcFinder::findSources(std::string dir, std::string cat, 
			      HourAngle ra, Declination dec, Angle radius,
			      std::vector<HourAngle>& ras, std::vector<Declination>& decs)
{
  // First two arguments are the catalog file and catalog ID

  string catalogDir = dir;
  string catalog    = cat;

  NvssReader      nvss;
  FirstFitsReader first;

  PtSrcReader* reader=0;

  std::ostringstream os;
  if(strcmp(catalog.c_str(), "nvss")==0) {
    reader = &nvss;
    os << catalogDir << "/nvss.fits";
  } else {
    reader = &first;
    os << catalogDir << "/first.fits";
  }

  // Set the catalog file

  reader->setCatalogFile(os.str());

  // Default search criteria to min/max fluxes

  Flux fMin = PtSrcReader::minFlux_;
  Flux fMax = PtSrcReader::maxFlux_;

  // And search

  findSources(reader, ra, dec, radius, fMin, fMax, ras, decs);
}

void PtSrcFinder::findSources(PtSrcReader* reader, 
			      HourAngle& ra, Declination& dec, Angle& radius, 
			      Flux& fMin, Flux& fMax,
			      std::vector<HourAngle>& ras, std::vector<Declination>& decs)
{
  // First pass to determine how many sources

  // Now find sources and return the found values

  PtSrcReader::Source src;

  // Set the RA range to search

  COUT("Setting raRange to " << radius << " about " << ra << " " << dec);
  reader->setRaRange(ra, dec, radius);

  // Make sure the catalog file is open

  reader->openCatalogFile();

  // Iterate, reading from the catalog file until the EOF is reached

  while(!reader->eof()) {

    src = reader->readNextEntry();

    // If the source flux is within range, check the angle, but don't
    // do otherwise, since this is expensive.

    if(src.peak_ >= fMin && src.peak_ <= fMax) {
      if(reader->checkAngle(src, ra, dec, radius)) {

	reader->applyCorrections(src);

	ras.push_back(src.ra_);
	decs.push_back(src.dec_);
      }
    }
  }

  // And close the file

  reader->closeCatalogFile();
}

