// $Id: FirstReader.h,v 1.5 2006/09/08 17:38:28 gcpdaq Exp $

#ifndef GCP_UTIL_FIRSTREADER_H
#define GCP_UTIL_FIRSTREADER_H

/**
 * @file FirstReader.h
 * 
 * Tagged: Tue Aug  8 18:13:20 PDT 2006
 * 
 * @version: $Revision: 1.5 $, $Date: 2006/09/08 17:38:28 $
 * 
 * @author Erik Leitch
 */
#include <fstream>
#include <string>
#include <vector>

#include "gcp/util/Angle.h"
#include "gcp/util/Flux.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/PtSrcReader.h"
#include "gcp/util/String.h"

namespace gcp {
  namespace util {

    class FirstReader {
    public:

      /**
       * Constructor.
       */
      FirstReader();

      /**
       * Destructor.
       */
      virtual ~FirstReader();

      void applyCorrections(PtSrcReader::Source& src);

      void setRestoringBeam(PtSrcReader::Source& src);
      
      // Set the size error for an axis

      void setSizeError(PtSrcReader::Source& src, Angle& axis, Angle& error);
      
      // Get the SNR for this source

      double getSnr(PtSrcReader::Source& src);
	
      // Construct the average beam for uncertainty calculations
      
      Angle getAvBeam(PtSrcReader::Source& src);
      
      // Set position errors
      
      void setPositionErrors(PtSrcReader::Source& src);

    private:

      static Angle eps_;

      static Angle northResMaj_;
      static Angle northResMin_;

      static Angle medianResMaj_;
      static Angle medianResMin_;

      static Angle southResMaj_;
      static Angle southResMin_;

      static Angle medianDecLimit_;
      static Angle southDecLimit_;

    }; // End class FirstReader

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_FIRSTREADER_H
