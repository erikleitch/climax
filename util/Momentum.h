// $Id: $

#ifndef GCP_UTIL_MOMENTUM_H
#define GCP_UTIL_MOMENTUM_H

/**
 * @file Momentum.h
 * 
 * Tagged: Wed Oct 23 17:23:54 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Mass.h"
#include "gcp/util/Speed.h"

namespace gcp {
  namespace util {

    class Momentum : public ConformableQuantity {
    public:

      /**
       * Constructor.
       */
      Momentum();

      Momentum(const Momentum& mom);
      Momentum(Momentum& mom);

      /**
       * Destructor.
       */
      virtual ~Momentum();

      void setGramsCmPerSec(double gCmPerSec);
      void setKgMetersPerSec(double kgMeterPerSec);

      double gramCmPerSec();
      double kgMeterPerSec();

      Speed operator/(Mass& mass);
      Mass  operator/(Speed& speed);

      Speed operator/(const Mass& mass);
      Mass  operator/(const Speed& speed);

      void operator=(const Momentum& mom);
      void operator=(Momentum& mom);

    }; // End class Momentum

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_MOMENTUM_H
