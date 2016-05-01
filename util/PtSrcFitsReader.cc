#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/PtSrcFitsReader.h"

#include <iomanip>

using namespace std;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
PtSrcFitsReader::PtSrcFitsReader() : PtSrcReader() 
{
  initialize();
}

/**.......................................................................
 * Constructor.
 */
PtSrcFitsReader::PtSrcFitsReader(std::string catalogFile) :
  PtSrcReader(catalogFile)
{
  initialize();
}

/**.......................................................................
 * Initialize critical members
 */
void PtSrcFitsReader::initialize()
{
  status_   = 0;
  fitsFile_ = 0;
  initRange();

  for(unsigned i=0; i < chunkSize_; i++)
    sourceNames_[i] = 0;

  for(unsigned i=0; i < chunkSize_; i++)
    sourceNames_[i] = (char*)malloc(10);
}

/**.......................................................................
 * Destructor.
 */
PtSrcFitsReader::~PtSrcFitsReader() 
{
  for(unsigned i=0; i < chunkSize_; i++) {
    if(sourceNames_[i] != 0) {
      free(sourceNames_[i]);
      sourceNames_[i] = 0;
    }
  }
}

/**.......................................................................
 * Open the catalog file
 */
void PtSrcFitsReader::openCatalogFile()
{
  closeCatalogFile();

  // Now open it

  if(fits_open_file(&fitsFile_, catalogFile_.c_str(), READONLY, &status_))
      throwCfitsioError(status_);

  // And position the file at the first table HDU
    
  int hduType;
  if(fits_movabs_hdu(fitsFile_, 2, &hduType, &status_)) 
    throwCfitsioError(status_);

  // And check that this is a valid binary table

  if (hduType != BINARY_TBL)  
    ThrowError("Expected to find a binary table in this HDU");

  // Read the file indices

  int nfound;
  if (fits_read_keys_lng(fitsFile_, "INDEX", 0, 24, indices_, &nfound, &status_))
    throwCfitsioError(status_);

  if(nfound != 24)
    ThrowError("Found an unexpected number of indices:" << nfound);

  // Finally, read the number of rows in the table

  if(fits_get_num_rows(fitsFile_, &nRowTotal_, &status_))
    throwCfitsioError(status_);

  // Pad out the index to make the indexing logic cleaner

  indices_[24] = nRowTotal_ + 1;

  // And get the number of columns

  int nCol=0;
  if(fits_get_num_cols(fitsFile_, &nCol, &status_))
    throwCfitsioError(status_);

  // Step to the first range

  incrementRange();
}

/**.......................................................................
 * Close the catalog file
 */
void PtSrcFitsReader::closeCatalogFile()
{
  // Close the file if it is already open

  if(fitsFile_) {

    if(fits_close_file(fitsFile_, &status_)) 
      throwCfitsioError(status_);
    
    fitsFile_ = 0;
  }

  initRange();
}

/**.......................................................................
 * Read the next entry from the catalog file
 */
PtSrcReader::Source PtSrcFitsReader::readNextEntry()
{
  // Which chunk would the source be in?

  unsigned iChunk = iRow_/chunkSize_;

  // If the source is not located within the last read chunk, read the
  // next one.  Increment by one, so that iChunk_ = 0 in the test
  // below will correspond to no data having been read

  if(iChunk+1 > iChunk_)
    readNextChunk();

  // Increment the row counter

  ++iRow_;

  // Return the source

  return parseData();
}

/**.......................................................................
 * Read the next chunk of data from the FITS file
 *
 * Data structure from the NVSS catalog.ps file:   
 *
 *   1 RA(2000)     J2000 Right Ascension (degrees)
 *   2 DEC(2000)    J2000 declination (degrees)
 *   3 PEAK INT     Peak Stokes I (Jy/beam)
 *   4 MAJOR AX     Fitted major axis (degrees)
 *   5 MINOR AX     Fitted minior axis (degrees)
 *   6 POSANGLE     Fitted position angle of the major axis (degrees)
 *   7 Q CENTER     Interpolated Q value at position of I peak (Jy/beam)
 *   8 U CENTER     Interpolated U value at position of I peak (Jy/beam)
 *   9 P FLUX       Integrated polarized (linear) flux (Jy)
 *  10 I RMS        RMS noise in Stokes I image (Jy/beam)
 *  11 POL RMS      RMS noise in Stokes Q,U image (Jy/beam)
 *  12 RES RMS      RMS Stokes I residual after fit
 *  13 RES PEAK     Peak Stokes I residual (Jy/beam)
 *  14 RES FLUX     Integrated Stokes I residual (Jy)
 *  15 CENTER X     X pixel in FIELD of center (pixels)
 *  16 CENTER Y     Y pixel in FIELD of center (pixels)
 *  17 FIELD        Name of the 4x4 degree field (string)
 *  18 JD PROCESSED Julian date on which this entry was derived from an image (day)

 */
void PtSrcFitsReader::readNextChunk()
{
  float floatnull = 0;
  static char strnull[10] = {" "};
  int anynull;

  // If we have already read the last chunk of the current range,
  // increment to the next range

  if(iChunk_ == nChunk_) {
    incrementRange();
  }

  // Calculate the starting row to read from (cfitsio wants this
  // indexed from 1).

  long startRow = (iChunk_ * chunkSize_) + rowMin_;

  // Calculate the number of elements to read on this pass.  If this
  // is the last chunk in this range, read only the remaining
  // elements.

  long nElement = (iChunk_ == nChunk_ - 1) ?  (nRow_ % chunkSize_) : chunkSize_;

  // Parse the data we just read

  readFitsData(startRow, nElement);

  // And increment which chunk we are on

  ++iChunk_;
}

/**.......................................................................
 * Return true if the input stream for the catalog file is at the EOF
 * mark
 */
bool PtSrcFitsReader::eof()
{
  return (status_ == END_OF_FILE) || nRow_ == 0 || (iRange_ == nRange_ && iRow_ == nRow_-1);
}

/**.......................................................................
 * Overload the base-class method to set up RA ranges to search
 */
void PtSrcFitsReader::setRaRange(HourAngle& ra, Declination& dec, Angle& radius)
{
  PtSrcReader::setRaRange(ra, dec, radius);

  // Note that stop indices can range to 24, since I have added an
  // extra element to the index array

  // If raMin_ < raMax_, then we are not crossing the RA=0 boundary,
  // and only have to search one RA range
 
  if(raMin_ < raMax_) {

    nRange_ = 1;

    rangeStartInd_[0] = (unsigned int) raMin_.hours();
    rangeStopInd_[0]  = (unsigned int) raMax_.hours() + 1;

    // Else we are.  In this case, the catalog needs to be searched in
    // two ranges for maximum efficiency

  } else {

    nRange_ = 2;

    rangeStartInd_[0] = 0;
    rangeStopInd_[0]  = (unsigned int) raMax_.hours() + 1;

    rangeStartInd_[1] = (unsigned int) raMin_.hours();
    rangeStopInd_[1]  = 24;
  }
}

void PtSrcFitsReader::initRange()
{
  iChunk_   = 0;
  nChunk_   = 0;
  iRow_     = 0;
  iRange_   = 0;
  nRange_   = 1;
}

/**.......................................................................
 * Increment to the next range.
 */
void PtSrcFitsReader::incrementRange()
{
  // If we are already on the second range, we are done.

  if(iRange_ == nRange_)
    return;

  // The starting row for the next set of chunk reads

  rowMin_ = indices_[rangeStartInd_[iRange_]];

  // Set the total number of rows to read.  indices_[] contains the
  // starting row of the first source in the coresponding range, so
  // the number of sources to read is just the difference.

  nRow_   =  indices_[rangeStopInd_[iRange_]] - indices_[rangeStartInd_[iRange_]];

  // Set the number of chunks we will read.  

  nChunk_ = nRow_/chunkSize_ + 1;

  // And reset the chunk and row counters

  iRow_   = 0;
  iChunk_ = 0;

  // And increment the range counter

  ++iRange_;
}
