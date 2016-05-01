#ifndef GCP_UTIL_FREQUENCY_H
#define GCP_UTIL_FREQUENCY_H

/**
 * @file Frequency.h
 * 
 * Tagged: Fri Aug 20 12:44:45 PDT 2004
 * 
 * @author Erik Leitch
 */
#include <iostream>

#include "gcp/util/Energy.h"
#include "gcp/util/Speed.h"

namespace gcp {
  namespace util {
    
    class Energy;
    class Rx;
    class Wavelength;

    class Frequency : public ConformableQuantity {
    public:
      
      // A few useful conversions

      static const double HzPerGHz_;
      static const double HzPerMHz_;

      class MegaHz {};
      class GigaHz {};

      /**
       * Constructor.
       */
      Frequency();
      Frequency(const MegaHz& units, double MHz);
      Frequency(const GigaHz& units, double GHz);
      Frequency(Wavelength& wavelength);
      
      /**
       * Destructor.
       */
      virtual ~Frequency();
      
      // Set the frequency, in GHz

      void setGHz(double GHz);

      // Set the frequency, in MHz

      void setMHz(double MHz);

      // Set the frequency, in MHz

      void setHz(double Hz);

      // Return the frequency, in GHz

      inline double GHz() {
	return val_ / HzPerGHz_;
      }

      // Return the frequency, in MHz

      inline double MHz() {
	return val_ / HzPerMHz_;
      }

      inline unsigned short yigUnits() {
	return (unsigned short)MHz();
      }

      // Return the frequency, in Hz

      inline double Hz() const {
	return val_;
      }

      double microns();
      double centimeters();
      double meters();

      Wavelength wavelength();

      /**
       * Allows cout << Length
       */
      friend std::ostream& operator<<(std::ostream& os, Frequency& frequency);

      Frequency operator-(Frequency& frequency);
      Frequency operator+(Frequency& frequency);
      double operator/(Frequency& frequency);

      bool operator<(Frequency& frequency);
      bool operator>(Frequency& frequency);

      void operator=(const Energy& energy);
      void operator=(Energy& energy);

      void operator=(const Frequency& frequency);
      void operator=(Frequency& frequency);

      void initialize();

    private:

      friend class Rx;
      friend class CorrelatorBand;

      // Constructor -- only Rx can call this constructor

      Frequency(double Hz);

    }; // End class Frequency
    
  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_FREQUENCY_H
