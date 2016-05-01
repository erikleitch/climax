#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/NvssReader.h"
#include "gcp/util/FirstFitsReader.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "dir",            "./", "s", "directory where catalog resides"},
  { "cat",          "nvss", "s", "catalog to use (nvss | first)"},
  { "ra",       "12:00:00", "s", "RA to search"},
  { "dec",      "34:00:00", "s", "DEC to search"},
  { "rad",      "00:02:00", "s", "Radius to search"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

void findSources(PtSrcReader* reader, 
		 HourAngle& ra, Declination& dec, Angle& radius, 
		 Flux& fMin, Flux& fMax);

int Program::main()
{
  // First two arguments are the catalog file and catalog ID

  string catalogDir = Program::getStringParameter("dir");
  string catalog    = Program::getStringParameter("cat");

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

  HourAngle ra;
  Declination dec;

  // Parse input arguments

  ra.setHours(Program::getStringParameter("ra"));
  dec.setDegrees(Program::getStringParameter("dec"));
  Angle radius;

  radius.setDegrees(Program::getStringParameter("rad"));

  // Default search criteria to min/max fluxes

  Flux fMin = PtSrcReader::minFlux_;
  Flux fMax = PtSrcReader::maxFlux_;

  // And search

  findSources(reader, ra, dec, radius, fMin, fMax);
}

void findSources(PtSrcReader* reader, 
		 HourAngle& ra, Declination& dec, Angle& radius, 
		 Flux& fMin, Flux& fMax)
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

	COUT(src.ra_ << " " << src.dec_ << " dist = " << src.distance_.arcmin());
      }
    }
  }

  // And close the file

  reader->closeCatalogFile();
}
