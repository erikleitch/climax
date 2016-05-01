#ifndef GCP_UTIL_TEMPERATURE_H
#define GCP_UTIL_TEMPERATURE_H

/**
 * @file Temperature.h
 * 
 * Tagged: Wed Dec  1 11:48:14 PST 2004
 * 
 * @author Erik Leitch
 */
#include "gcp/util/ConformableQuantity.h"
#include "gcp/util/Unit.h"
#include "gcp/util/Pressure.h"
#include "gcp/util/NumberDensity.h"

namespace gcp {
  namespace util {
    
    class Pressure;
    class NumberDensity;

    class Temperature : public ConformableQuantity {
    public:
      
      class Kelvin : public Unit {
      public:
	void addNames();
      };

      class Centigrade{};
      class Celsius   {};
      class Fahrenheit{};
      class MicroKelvin{};

      /**
       * Constructor.
       */
      Temperature();
      Temperature(const Kelvin& units,      double kelvins);
      Temperature(const MicroKelvin& units, double microKelvins);
      Temperature(const Centigrade& units,  double centigrade);
      Temperature(const Celsius& units,     double celsius);
      Temperature(const Fahrenheit& units,  double fahrenheit);
      
      /**
       * Destructor.
       */
      virtual ~Temperature();
      
      void setC(double centigrade);
      void setKeV(double keV);
      void setF(double fahrenheit);
      void setK(double kelvin);
      void setMicroK(double microKelvins);

      double C();
      double F();
      double K();
      double microK();
      double milliK();
      double keV();

      static const double kelvinZeroPointInC_;
      static const double kelvinPerKev_;

      void initialize();

      Temperature operator+(Temperature& temp);
      Temperature operator+(const Temperature& temp);

      Temperature operator*(double fac);
      Pressure operator*(const NumberDensity& nd);

      double operator/(Temperature& temp);
      double operator/(const Temperature& temp);

      void operator=(const Temperature& var);
      void operator=(Temperature& var);

      // Allows cout << Temperature

      friend std::ostream& operator<<(std::ostream& os, Temperature& temp);

    }; // End class Temperature
    
  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_TEMPERATURE_H
