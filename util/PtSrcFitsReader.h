// $Id: PtSrcFitsReader.h,v 1.2 2007/06/23 03:56:15 gcpdaq Exp $

#ifndef GCP_UTIL_PTSRCFITSREADER_H
#define GCP_UTIL_PTSRCFITSREADER_H

/**
 * @file NvssReader.h
 * 
 * Tagged: Tue Aug 15 18:00:57 PDT 2006
 * 
 * @version: $Revision: 1.2 $, $Date: 2007/06/23 03:56:15 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/PtSrcReader.h"

#include "fitsio.h"

namespace gcp {
  namespace util {

    class PtSrcFitsReader : public PtSrcReader {
    public:

      /**
       * Constructor.
       */
      PtSrcFitsReader(std::string catalogFile);
      PtSrcFitsReader();

      /**
       * Destructor.
       */
      virtual ~PtSrcFitsReader();

      /**
       * Initialize critical members
       */
      void initialize();

      // Open the catalog file

      void openCatalogFile();

      // Close the catalog file

      void closeCatalogFile();

      // Read the next entry from the catalog file

      PtSrcReader::Source readNextEntry();

      // Read the next chunk of data from the FITS file

      void readNextChunk();

      // Return true if we are at the end of file

      bool eof();

      void setRaRange(HourAngle& ra, Declination& dec, Angle& radius);
      void initRange();

    protected:

      //------------------------------------------------------------
      // FITS-related members

      // A pointer to a FITS file

      fitsfile* fitsFile_;

      // A status flag used by cfitsio routines

      int status_;

      // An array of indices for telling which source number
      // corresponds to each new RA range

      long indices_[25];

      static const unsigned chunkSize_ = 10;
      unsigned nChunk_, iChunk_;

      long nRow_, nRowTotal_;
      unsigned iRow_;

      double ras_[chunkSize_];
      double decs_[chunkSize_];
      float peakFluxes_[chunkSize_];
      float rmsFluxes_[chunkSize_];
      float majorAxes_[chunkSize_];
      float minorAxes_[chunkSize_];
      float positionAngles_[chunkSize_];
      char* sourceNames_[chunkSize_];

      // The number of RA ranges we are searching

      unsigned int nRange_;
      unsigned int iRange_;
      unsigned rowMin_;

      // Up to two start indices

      unsigned int rangeStartInd_[2];

      // Up to two stop indices

      unsigned int rangeStopInd_[2];

      // Force inheritors to define the following functions

      virtual PtSrcReader::Source parseData() = 0;
      virtual void readFitsData(long startRow, long nElement) = 0;

      // Increment to the next range.

      void incrementRange();

    }; // End class PtSrcFitsReader

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PTSRCFITSREADER_H
