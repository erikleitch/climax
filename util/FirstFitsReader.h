// $Id: FirstFitsReader.h,v 1.1 2006/09/08 17:41:56 gcpdaq Exp $

#ifndef GCP_UTIL_FIRSTFITSREADER_H
#define GCP_UTIL_FIRSTFITSREADER_H

/**
 * @file FirstFitsReader.h
 * 
 * Tagged: Tue Aug  8 18:13:20 PDT 2006
 * 
 * @version: $Revision: 1.1 $, $Date: 2006/09/08 17:41:56 $
 * 
 * @author Erik Leitch
 */
#include <fstream>
#include <string>
#include <vector>

#include "gcp/util/Angle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/FirstReader.h"
#include "gcp/util/Flux.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/PtSrcFitsReader.h"
#include "gcp/util/String.h"

namespace gcp {
  namespace util {

    class FirstFitsReader : public PtSrcFitsReader {
    public:

      /**
       * Constructor.
       */
      FirstFitsReader(std::string catalogFile);
      FirstFitsReader();

      /**
       * Destructor.
       */
      virtual ~FirstFitsReader();

      //------------------------------------------------------------
      // Overloaded base-class methods from PtSrcFitsReader

      PtSrcReader::Source parseData();

      void readFitsData(long startRow, long nElement);

      void applyCorrections(PtSrcReader::Source& src);

    private:

      FirstReader firstReader_;
      float warns_[PtSrcFitsReader::chunkSize_];
      float intFluxes_[PtSrcFitsReader::chunkSize_];
      float decMajorAxes_[PtSrcFitsReader::chunkSize_];
      float decMinorAxes_[PtSrcFitsReader::chunkSize_];
      float decPositionAngles_[PtSrcFitsReader::chunkSize_];

    }; // End class FirstFitsReader

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_FIRSTFITSREADER_H
