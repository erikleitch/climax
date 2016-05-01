#ifndef GCP_UTIL_COLUMNDENSITY_H
#define GCP_UTIL_COLUMNDENSITY_H

/**
 * @file ColumnDensity.h
 * 
 * Tagged: Tue Jan  8 17:38:52 NZDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Area.h"
#include "gcp/util/ConformableQuantity.h"
#include "gcp/util/Temperature.h"
#include "gcp/util/Volume.h"

namespace gcp {
  namespace util {

    class Area;
    class Temperature;
    class Volume;

    class ColumnDensity : public ConformableQuantity {
    public:

      /**
       * Constructor.
       */
      ColumnDensity();

      /**
       * Destructor.
       */
      virtual ~ColumnDensity();

      /**
       * Set the volume of this object
       */
      void setInverseSquaredCentimeters(double cm2); 
      void setInverseSquaredMeters(double m2); 

      /**
       * Get the density of this object
       */
      inline double inverseSquaredCentimeters() const {
	return val_;
      }

      inline double inverseSquaredMeters() const {
	double cmPerM = 100;
	return val_ * cmPerM * cmPerM;
      }

      // A Column density multiplied by an area is a number

      double operator*(const Area& area);

      void operator=(const ColumnDensity& density);
      void operator=(ColumnDensity& density);

      friend std::ostream& operator<<(std::ostream& os, ColumnDensity& cd);
      friend std::ostream& operator<<(std::ostream& os, const ColumnDensity& cd);

    }; // End class ColumnDensity

    std::ostream& operator<<(std::ostream& os, ColumnDensity& cd);
    std::ostream& operator<<(std::ostream& os, const ColumnDensity& cd);

    // A factor divided by a column density is an area

    Area operator/(double fac, const ColumnDensity& cd);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_COLUMNDENSITY_H
