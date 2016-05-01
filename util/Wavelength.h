#ifndef GCP_UTIL_WAVELENGTH_H
#define GCP_UTIL_WAVELENGTH_H

/**
 * @file Wavelength.h
 * 
 * Tagged: Wed Dec  1 11:58:54 PST 2004
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Energy.h"
#include "gcp/util/Length.h"
#include "gcp/util/Speed.h"

namespace gcp {
  namespace util {
    
    class Energy;
    class Frequency;

    class Wavelength : public Length {
    public:
      
      class Angstroms {};

      static const double cmPerAngstrom_;
      static Speed lightSpeed_;

      /**
       * Constructor.
       */
      Wavelength();
      Wavelength(const Frequency& frequency);
      Wavelength(const Length::Centimeters& units, double cm);
      Wavelength(const Microns& units, double microns);
      
      /**
       * Destructor.
       */
      virtual ~Wavelength();
      
      void operator=(const Energy& energy);
      void operator=(Energy& energy);

      void operator=(const Wavelength& wavelength);
      void operator=(Wavelength& wavelength);

      void setFrequency(Frequency& freq);
      void setFrequency(const Frequency& freq);

      void setMicrons(double microns);
      
      void setAngstroms(double angstroms);
    
      double microns();
      double angstroms();
      double cm();
      double centimeters();
      double meters();

      Frequency frequency();

      void initialize();

    }; // End class Wavelength
    
  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_WAVELENGTH_H
