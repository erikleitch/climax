#ifndef GCP_UTIL_NUMBERDENSITY_H
#define GCP_UTIL_NUMBERDENSITY_H

/**
 * @file NumberDensity.h
 * 
 * Tagged: Tue Jan  8 17:38:52 NZDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Area.h"
#include "gcp/util/ConformableQuantity.h"
#include "gcp/util/Pressure.h"
#include "gcp/util/Temperature.h"
#include "gcp/util/Volume.h"

namespace gcp {
  namespace util {

    class Area;
    class Pressure;
    class Temperature;
    class Volume;

    class NumberDensity : public ConformableQuantity {
    public:

      /**
       * Constructor.
       */
      NumberDensity();

      /**
       * Destructor.
       */
      virtual ~NumberDensity();

      /**
       * Set the volume of this object
       */
      void setInverseCubicCentimeters(double cm3); 
      void setInverseCubicMeters(double m3); 

      /**
       * Get the density of this object
       */
      inline double inverseCubicCentimeters() const {
	return val_;
      }

      inline double inverseCubicMeters() const {
	double cmPerM = 100;
	return val_ * cmPerM * cmPerM * cmPerM;
      }

      // A number density multiplied by a volume is a number

      double operator*(const Volume& volume);
      double operator*(Volume& volume);

      // A number density multiplied by a temperature implies a
      // pressure (P = nkT)

      Pressure operator*(const Temperature& temp);
      Pressure operator*(Temperature& temp);

      void operator=(const NumberDensity& density);
      void operator=(NumberDensity& density);

      friend std::ostream& operator<<(std::ostream& os, NumberDensity& nd);
      friend std::ostream& operator<<(std::ostream& os, const NumberDensity& nd);

    }; // End class NumberDensity

    std::ostream& operator<<(std::ostream& os, NumberDensity& nd);
    std::ostream& operator<<(std::ostream& os, const NumberDensity& nd);

    // A factor divided by a number density is a volume

    Volume operator/(double fac, const NumberDensity& nd);
    Volume operator/(double fac, NumberDensity& nd);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_NUMBERDENSITY_H
