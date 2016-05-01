/**.......................................................................
 * MATLAB Mex file for accessing point source catalogs
 *
 */
#include <math.h>

#include <string>
#include <iostream>
#include <sstream>

#include "gcp/util/Debug.h"
#include "gcp/util/FirstFitsReader.h"
#include "gcp/util/NvssReader.h"

#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

#include "mex.h"
#include "matrix.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::matlab;

/**.......................................................................
 * Return a list of sources within radius of the requested position.
 *
 * Note that the FIRST catalog is in J2000 coordinates, so that the
 * requested positions should be too.
 */
void findSources(mxArray* plhs[],
		 PtSrcReader* reader, 
		 HourAngle& ra, Declination& dec, Angle& radius, 
		 Flux& fMin, Flux& fMax)
{
  // First pass to determine how many sources

  unsigned nSrc = reader->countSources(ra, dec, radius, fMin, fMax);

  // Create an empty arraymap structure 

  int dims[2]={1,1};
  plhs[0] = mxCreateStructArray(2, dims, 0, NULL);

  double* raPtr        = MexHandler::addNamedDoubleStructField(plhs[0], "ra",        nSrc);
  double* raErrPtr     = MexHandler::addNamedDoubleStructField(plhs[0], "raErr",     nSrc);
  double* decPtr       = MexHandler::addNamedDoubleStructField(plhs[0], "dec",       nSrc);
  double* decErrPtr    = MexHandler::addNamedDoubleStructField(plhs[0], "decErr",    nSrc);

  double* peakPtr      = MexHandler::addNamedDoubleStructField(plhs[0], "peak",      nSrc);
  double* peakErrPtr   = MexHandler::addNamedDoubleStructField(plhs[0], "peakErr",   nSrc);

  double* intPtr       = MexHandler::addNamedDoubleStructField(plhs[0], "int",       nSrc);
  double* intErrPtr    = MexHandler::addNamedDoubleStructField(plhs[0], "intErr",    nSrc);

  double* rmsPtr       = MexHandler::addNamedDoubleStructField(plhs[0], "rms",       nSrc);

  double* decMajPtr    = MexHandler::addNamedDoubleStructField(plhs[0], "decMaj",    nSrc);
  double* decMajErrPtr = MexHandler::addNamedDoubleStructField(plhs[0], "decMajErr", nSrc);
  double* decMinPtr    = MexHandler::addNamedDoubleStructField(plhs[0], "decMin",    nSrc);
  double* decMinErrPtr = MexHandler::addNamedDoubleStructField(plhs[0], "decMinErr", nSrc);
  double* decPaPtr     = MexHandler::addNamedDoubleStructField(plhs[0], "decPa",     nSrc);
  double* decPaErrPtr  = MexHandler::addNamedDoubleStructField(plhs[0], "decPaErr",  nSrc);

  double* resMajPtr    = MexHandler::addNamedDoubleStructField(plhs[0], "resMaj",    nSrc);
  double* resMinPtr    = MexHandler::addNamedDoubleStructField(plhs[0], "resMin",    nSrc);

  double* warnPtr      = MexHandler::addNamedDoubleStructField(plhs[0], "warn",      nSrc);
  double* distPtr      = MexHandler::addNamedDoubleStructField(plhs[0], "dist",      nSrc);

  // Now find sources and return the found values

  PtSrcReader::Source src;

  // Set the RA range to search

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

	*(raPtr++)        = src.ra_.hours();
	*(raErrPtr++)     = src.raErr_.hours();
	*(decPtr++)       = src.dec_.degrees();
	*(decErrPtr++)    = src.decErr_.degrees();

	*(peakPtr++)      = src.peak_.mJy();
	*(peakErrPtr++)   = src.peakErr_.mJy();

	*(intPtr++)       = src.int_.mJy();
	*(intErrPtr++)    = src.intErr_.mJy();

	*(rmsPtr++)       = src.rms_.mJy();

	*(decMajPtr++)    = src.decMaj_.arcsec();
	*(decMajErrPtr++) = src.decMajErr_.arcsec();
	*(decMinPtr++)    = src.decMin_.arcsec();
	*(decMinErrPtr++) = src.decMinErr_.arcsec();
	*(decPaPtr++)     = src.decPa_.degrees();
	*(decPaErrPtr++)  = src.decPaErr_.degrees();

	*(resMajPtr++)    = src.resMaj_.arcsec();
	*(resMinPtr++)    = src.resMin_.arcsec();

	*(warnPtr++)      = (double)src.warn_;
	*(distPtr++)      = src.distance_.arcsec();
      }
    }
  }

  // And close the file

  reader->closeCatalogFile();
}

void getSources(mxArray* plhs[], const mxArray* prhs[], unsigned nrhs)
{
  // First two arguments are the catalog file and catalog ID

  string catalogDir = MexParser::getString(prhs[0]);
  string catalog    = MexParser::getString(prhs[1]);

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

  if(mxIsChar(prhs[2]))
    ra.setHours(MexParser::getString(prhs[2]));
  else
    ra.setHours(*MexParser::getDoubleData(prhs[2]));

  if(mxIsChar(prhs[3]))
    dec.setDegrees(MexParser::getString(prhs[3]));
  else
    dec.setDegrees(*MexParser::getDoubleData(prhs[3]));

  Angle radius;

  if(mxIsChar(prhs[4]))
    radius.setDegrees(MexParser::getString(prhs[4]));
  else
    radius.setDegrees(*MexParser::getDoubleData(prhs[4]));

  // Default search criteria to min/max fluxes

  Flux fMin = PtSrcReader::minFlux_;
  Flux fMax = PtSrcReader::maxFlux_;

  // If arguments were supplied, use them

  switch(nrhs) {
  case 7:
    fMax.setMilliJy(*(MexParser::getDoubleData(prhs[6])));
  case 6:
    fMin.setMilliJy(*(MexParser::getDoubleData(prhs[5])));
    break;
  default:
    break;
  };

  // And search

  findSources(plhs, reader, ra, dec, radius, fMin, fMax);
} 

/**.......................................................................
 * Main
 */
void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  gcp::util::Logger::installStdoutPrintFn(MexHandler::stdoutPrintFn);
  gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);
  
  gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
  gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
  gcp::util::ErrHandler::installLogFn(MexHandler::logFn);

  try {  

    if(nrhs < 5 || nrhs > 7) {

      CERR(std::endl << "Wrong number of arguments. Should be:" << std::endl << std::endl
	   << "  catalogFile catalog(nvss|first) ra(hh:mm:ss.ss) dec(dd:mm:ss.ss) radius(dd:mm:ss.ss) [fmin(mJy/beam) fmax(mJy/beam)]" 
	   << std::endl << std::endl 
	   << "where arguments in [] are optional" << std::endl);

      return;
    }

    getSources(plhs, prhs, nrhs);

  } catch(gcp::util::Exception& err) {
    mexErrMsgTxt(err.what());
  } catch(...) {
    mexErrMsgTxt("Caught an unknown error");
  }
}

