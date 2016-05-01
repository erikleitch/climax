#ifndef GCP_UTIL_MASS_H
#define GCP_UTIL_MASS_H

/**
 * @file Mass.h
 * 
 * Tagged: Wed Dec  1 11:58:54 PST 2004
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Density.h"
#include "gcp/util/Energy.h"
#include "gcp/util/Temperature.h"
#include "gcp/util/Wavelength.h"

namespace gcp {
  namespace util {
    
    class Density;
    class Energy;
    class Frequency;
    class Temperature;
    class Wavelength;

    class Mass : public ConformableQuantity {
    public:
      
      class Gram{};

      /**
       * Constructor.
       */
      Mass();
      Mass(const Gram& units, double grams);
      
      /**
       * Destructor.
       */
      virtual ~Mass();

      void initialize(void);

      void setGrams(double g);
      void setKiloGrams(double kg);
      void setKg(double kg);
      void setSolarMass(double mSolar);

      double solarMass();
      double g();
      double kg();
      Energy energy();

      void operator=(const Energy& energy);
      void operator=(Energy& energy);

      void operator=(const Mass& mass);
      void operator=(Mass& mass);

      bool operator>(Mass& mass);
      bool operator>(const Mass& mass);

      bool operator>=(Mass& mass);
      bool operator>=(const Mass& mass);

      bool operator<(Mass& mass);
      bool operator<(const Mass& mass);

      bool operator<=(Mass& mass);
      bool operator<=(const Mass& mass);

      Mass operator*(double factor);
      void operator*=(double factor);

      Mass operator/(double factor);
      void operator/=(double factor);

      // A mass divided by a density is a volume

      Volume operator/(Density& density);

      static const double gPerKg_;

    }; // End class Mass
    
  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_MASS_H
