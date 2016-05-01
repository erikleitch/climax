#ifndef GCP_UTIL_ENERGY_H
#define GCP_UTIL_ENERGY_H

/**
 * @file Energy.h
 * 
 * Tagged: Wed Dec  1 11:58:54 PST 2004
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Length.h"
#include "gcp/util/Mass.h"
#include "gcp/util/Momentum.h"
#include "gcp/util/Speed.h"
#include "gcp/util/Temperature.h"
#include "gcp/util/Wavelength.h"

namespace gcp {
  namespace util {
    
    class Frequency;
    class Mass;
    class Pressure;
    class Temperature;
    class Wavelength;

    class Energy : public ConformableQuantity {
    public:
      
      /**
       * Constructor.
       */
      Energy();
      
      /**
       * Destructor.
       */
      virtual ~Energy();

      static const double joulePerErg_;
      static const double joulePerEv_;
      static const double joulePerBtu_;
      static const double evPerKev_;
      static const double ergPerEv_;
      static const double ergPerKev_;

      // Operators dealing implicitly with the energy of light

      void operator=(const Frequency& nu);
      void operator=(Frequency& nu);
      void operator=(const Temperature& temp);
      void operator=(Temperature& temp);
      void operator=(const Wavelength& wave);
      void operator=(Wavelength& wave);
      void operator=(const Mass& mass);
      void operator=(Mass& mass);

      double operator/(Energy& energy);

      // Energy divided by a volume is a pressure

      Pressure operator/(Volume& volume);

      // Energy divided by a speed is a momentum

      Momentum operator/(const Speed& speed);
      Momentum operator/(Speed& speed);

      void setErgs(double ergs);
      void setEv(double eV);
      void setKev(double keV);
      void setJoules(double joules);
      void setElectronVolts(double eV);
      void setBtu(double btu);

      double ergs();
      double joules();
      double eV();
      double keV();
      double btu();
      Mass mass();

      // Return the speed of an object of Mass m with this energy

      Speed speed(Mass& mass);
      double beta(Mass& mass);
      double gamma(Mass& mass);

      void operator=(const Energy& energy);
      void operator=(Energy& energy);

    }; // End class Energy
    
  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_ENERGY_H
