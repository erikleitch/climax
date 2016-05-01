// $Id: Intensity.h,v 1.2 2012/05/02 23:44:51 eml Exp $

#ifndef GCP_UTIL_INTENSITY_H
#define GCP_UTIL_INTENSITY_H

/**
 * @file Intensity.h
 * 
 * Tagged: Fri Aug  1 23:13:07 PDT 2008
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/02 23:44:51 $
 * 
 * @author tcsh: Erik Leitch.
 */
#include "gcp/util/ConformableQuantity.h"
#include "gcp/util/Flux.h"
#include "gcp/util/SolidAngle.h"

namespace gcp {
  namespace util {

    class Intensity : public ConformableQuantity {
    public:

      class JanskyPerSr {};
      class MegaJanskyPerSr {};

      /**
       * Constructor.
       */
      Intensity();
      Intensity(const JanskyPerSr& units, double JyPerSr);
      Intensity(const MegaJanskyPerSr& units, double MJyPerSr);

      /**
       * Destructor.
       */
      virtual ~Intensity();

      // Return the flux, in Jy

      inline double JyPerSr() {
	return val_;
      }

      inline double mJyPerSr() {
	return val_ / 1e-3;
      }

      inline double MJyPerSr() {
	return val_ / 1e6;
      }

      void setJyPerSr(double JyPerSr);
      void setMJyPerSr(double MJyPerSr);

      void operator=(const Intensity& intensity);
      void operator=(Intensity& intensity);

    private:

      void initialize();

    }; // End class Intensity

    Intensity operator/(const Flux& flux, const SolidAngle& omega);
    Intensity operator/(Flux& flux, SolidAngle& omega);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_INTENSITY_H
