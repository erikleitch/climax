#include "gcp/util/Exception.h"
#include "gcp/util/FirstFitsReader.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/String.h"

#include <iomanip>

using namespace std;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
FirstFitsReader::FirstFitsReader(std::string catalogFile):
  PtSrcFitsReader(catalogFile) {}

/**.......................................................................
 * Constructor.
 */
FirstFitsReader::FirstFitsReader() : PtSrcFitsReader() {}

/**.......................................................................
 * Destructor.
 */
FirstFitsReader::~FirstFitsReader() {}

/**.......................................................................
 * Read the next entry from the catalog file
 */
PtSrcReader::Source FirstFitsReader::parseData()
{
  // Calculate the index into the current chunk containing this source

  unsigned iSrc = iRow_ % chunkSize_;
  
  // Stick the relevant data into a source structure and return it.
  
  // Store the flux, rms and axes - we will need these for
  // error calculations below

  PtSrcReader::Source src;

  src.ra_.setDegrees(ras_[iSrc]);
  src.dec_.setDegrees(decs_[iSrc]);

  src.rawPeak_.setJy(peakFluxes_[iSrc]);
  src.peak_.setJy(peakFluxes_[iSrc]);
  src.int_.setJy(intFluxes_[iSrc]);

  src.fitMaj_.setDegrees(majorAxes_[iSrc]);
  src.fitMin_.setDegrees(minorAxes_[iSrc]);
  src.fitPa_.setDegrees(positionAngles_[iSrc]);

  src.decMaj_.setDegrees(decMajorAxes_[iSrc]);
  src.decMin_.setDegrees(decMinorAxes_[iSrc]);
  src.decPa_.setDegrees(decPositionAngles_[iSrc]);

  src.rms_.setJy(rmsFluxes_[iSrc]);

  // NVSS catalog assumes every source is real

  src.warn_ = (bool) warns_[iSrc];

  // Return the source

  return src;
}

/**.......................................................................
 * Read the next chunk of data from the FITS file
 *
 * Data structure from the NVSS catalog.ps file:   
 *
 *   1 RA(2000)     J2000 Right Ascension (degrees)
 *   2 DEC(2000)    J2000 declination (degrees)
 *   3 WARN         Warning flag (T/F)
 *   4 PEAK INT     Peak Stokes I (Jy/beam)
 *   5 INT FLUX     Integrated FLux(Jy)
 *   6 I RMS        RMS in the I image (Jy)
 *   7 DEC MAJOR AX Deconvolved major axis (degrees)
 *   8 DEC MINOR AX Deconvolved minor axis (degrees)
 *   9 DEC PA       Deconvolved PA (degrees)
 *  10 FIT MAJOR AX Fitted major axis (degrees)
 *  11 FIT MINOR AX Fitted minor axis (degrees)
 *  12 FIT PA       Fitted PA (degrees)
 */
void FirstFitsReader::readFitsData(long startRow, long nElement)
{
  float floatnull = 0;
  static char strnull[10] = {" "};
  int anynull;

  // Now read the data

  fits_read_col(fitsFile_, TDOUBLE,  1, startRow, 1, nElement, &floatnull,  ras_,
		&anynull, &status_);
	
  fits_read_col(fitsFile_, TDOUBLE,  2, startRow, 1, nElement, &floatnull,  decs_,
		&anynull, &status_);

  fits_read_col(fitsFile_, TFLOAT,   3, startRow, 1, nElement, &floatnull,  warns_,
		&anynull, &status_);

  fits_read_col(fitsFile_, TFLOAT,   4, startRow, 1, nElement, &floatnull,  peakFluxes_,
		&anynull, &status_);

  fits_read_col(fitsFile_, TFLOAT,   5, startRow, 1, nElement, &floatnull,  intFluxes_,
		&anynull, &status_);

  fits_read_col(fitsFile_, TFLOAT,   6, startRow, 1, nElement, &floatnull,  rmsFluxes_,
		&anynull, &status_);

  fits_read_col(fitsFile_, TFLOAT,   7, startRow, 1, nElement, &floatnull,  decMajorAxes_,
		&anynull, &status_);

  fits_read_col(fitsFile_, TFLOAT,   8, startRow, 1, nElement, &floatnull,  decMinorAxes_,
		&anynull, &status_);

  fits_read_col(fitsFile_, TFLOAT,   9, startRow, 1, nElement, &floatnull,  decPositionAngles_,
		&anynull, &status_);

  fits_read_col(fitsFile_, TFLOAT,  10, startRow, 1, nElement, &floatnull,  majorAxes_,
		&anynull, &status_);

  fits_read_col(fitsFile_, TFLOAT,  11, startRow, 1, nElement, &floatnull,  minorAxes_,
		&anynull, &status_);

  fits_read_col(fitsFile_, TFLOAT,  12, startRow, 1, nElement, &floatnull,  positionAngles_,
		&anynull, &status_);
}

void FirstFitsReader::applyCorrections(PtSrcReader::Source& src)
{
  firstReader_.applyCorrections(src);
}
