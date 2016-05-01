#ifndef GCP_UTIL_SPEED_H
#define GCP_UTIL_SPEED_H

/**
 * @file Speed.h
 * 
 * Tagged: Wed Dec  1 23:39:12 PST 2004
 * 
 * @author Erik Leitch
 */
#include "gcp/util/ConformableQuantity.h"
#include "gcp/util/HubbleConstant.h"
#include "gcp/util/Length.h"

namespace gcp {
  namespace util {
    
    class Speed : public ConformableQuantity {
    public:
      
      class CentimetersPerSec {};
      class KilometersPerSec {};
      class MetersPerSec {};

      /**
       * Constructor.
       */
      Speed();
      Speed(const CentimetersPerSec& units, double cmPerSec);
      Speed(const MetersPerSec& units, double mPerSec);
      
      Speed(const Speed& mom);
      Speed(Speed& mom);

      /**
       * Destructor.
       */
      virtual ~Speed();
      
      /**
       * Set a speed
       */
      void setCentimetersPerSec(double cmPerSec);
      void setCmPerSec(double cmPerSec);
      void setMetersPerSec(double mPerSec);
      void setMilesPerHour(double mph);

      /**
       * Get a speed
       */
      double centimetersPerSec();
      double cmPerSec();
      double metersPerSec();
      double mph();

      void initialize();

      static const unsigned secPerHour_;
      static const double metersPerMile_;

      Length operator/(HubbleConstant& H);
      double operator/(Speed& speed);

      void operator=(const Speed& s);
      void operator=(Speed& s);

    }; // End class Speed
    
  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_SPEED_H
