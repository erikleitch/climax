// $Id: $

#ifndef GCP_UTIL_VOLUME_H
#define GCP_UTIL_VOLUME_H

/**
 * @file Volume.h
 * 
 * Tagged: Tue Jan  8 17:38:52 NZDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Area.h"
#include "gcp/util/ConformableQuantity.h"
#include "gcp/util/Length.h"
#include "gcp/util/NumberDensity.h"

namespace gcp {
  namespace util {

    class Area;
    class Length;
    class NumberDensity;

    class Volume : public ConformableQuantity {
    public:

      /**
       * Constructor.
       */
      Volume();

      /**
       * Destructor.
       */
      virtual ~Volume();

      /**
       * Set the volume of this object
       */
      void setCubicCentimeters(double cm3); 

      /**
       * Get the volume of this object
       */
      inline double cubicCentimeters() const {
	return val_;
      }

      // A volume divided by a length is an area

      Area   operator/(const Length& length);

      // A volume divided by an area is a length

      Length operator/(const Area& area);

      void operator=(const Volume& volume);
      void operator=(Volume& volume);

      friend std::ostream& operator<<(std::ostream& os, Volume& volume);
      friend std::ostream& operator<<(std::ostream& os, const Volume& volume);

    }; // End class Volume

    std::ostream& operator<<(std::ostream& os, Volume& volume);
    std::ostream& operator<<(std::ostream& os, const Volume& volume);

    // A factor divided by a volume is a number density

    NumberDensity operator/(double fac, const Volume& volume);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_VOLUME_H
