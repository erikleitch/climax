// $Id: $

#ifndef GCP_UTIL_AREA_H
#define GCP_UTIL_AREA_H

/**
 * @file Area.h
 * 
 * Tagged: Tue Jan  8 17:38:52 NZDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/ConformableQuantity.h"
#include "gcp/util/Length.h"
#include "gcp/util/Volume.h"

namespace gcp {
  namespace util {

    class Length;
    class Volume;

    class Area : public ConformableQuantity {
    public:

      class SquaredCentimeters {};

      /**
       * Constructor.
       */
      Area();
      Area(const SquaredCentimeters& units, double cm2);
      void initialize(); 

      /**
       * Destructor.
       */
      virtual ~Area();

      /**
       * Set the area of this object
       */
      void setSquaredCentimeters(double sqcm); 
      void setSquaredMpc(double sqMpc); 

      /**
       * Get the area of this object
       */
      inline double squaredCentimeters() const {
	return val_;
      }

      static const double cm2PerParsec2_;
      static const double pc2PerMpc2_;

      /**
       * Get the area of this object
       */
      double squaredMpc();

      // An area divided by an area is a number

      double operator/(const Area& area);

      // An area divided by a length is a length

      Length operator/(const Length& length);

      // An area multiplied by a length is a volume

      Volume operator*(const Length& length);

      void operator=(const Area& area);
      void operator=(Area& area);

      friend std::ostream& operator<<(std::ostream& os, Area& area);
      friend std::ostream& operator<<(std::ostream& os, const Area& area);

    }; // End class Area

    std::ostream& operator<<(std::ostream& os, Area& area);
    std::ostream& operator<<(std::ostream& os, const Area& area);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_AREA_H
