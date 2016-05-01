// $Id: Flux.h,v 1.2 2012/05/02 23:44:51 eml Exp $

#ifndef GCP_UTIL_FLUX_H
#define GCP_UTIL_FLUX_H

/**
 * @file Flux.h
 * 
 * Tagged: Wed Sep 14 17:14:39 PDT 2005
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/02 23:44:51 $
 * 
 * @author Erik Leitch
 */
#include <iostream>

#include "gcp/util/ConformableQuantity.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/Temperature.h"
#include "gcp/util/SolidAngle.h"

namespace gcp {
  namespace util {

    class Flux : public ConformableQuantity {
    public:

      class Jansky {};
      class MilliJansky {};
      class MegaJansky {};

      /**
       * Constructor.
       */
      Flux();
      Flux(const Jansky& units, double Jy);
      Flux(const MilliJansky& units, double mJy);
      Flux(const MegaJansky& units, double MJy);
      Flux(Frequency& freq, Temperature& temp, SolidAngle& omega);

      /**
       * Destructor.
       */
      virtual ~Flux();

      void setJy(double Jy);
      void setMilliJy(double mJy);
      void setMegaJy(double MJy);

      // Return the flux, in Jy

      inline double Jy() {
	return val_;
      }

      // Return the flux, in mJy

      inline double mJy() {
	return val_ * mJyPerJy_;
      }

      // Return the flux, in MJy

      inline double MJy() {
	return val_ / JyPerMJy_;
      }

      static const double mJyPerJy_;
      static const double JyPerMJy_;

      void initialize();

      friend std::ostream& operator<<(std::ostream& os, Flux& flux);

      bool operator>=(Flux& flux);
      bool operator<=(Flux& flux);

      void operator=(const Flux& flux);
      void operator=(Flux& flux);

    }; // End class Flux

    std::ostream& operator<<(std::ostream& os, Flux& flux);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_FLUX_H
