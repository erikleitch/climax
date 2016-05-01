// $Id: Planck.h,v 1.1.1.1 2010/07/13 17:56:53 eml Exp $

#ifndef GCP_UTIL_PLANCK_H
#define GCP_UTIL_PLANCK_H

/**
 * @file Planck.h
 * 
 * Tagged: Thu Aug 21 11:15:49 PDT 2008
 * 
 * @version: $Revision: 1.1.1.1 $, $Date: 2010/07/13 17:56:53 $
 * 
 * @author Erik Leitch.
 */
#include "gcp/util/Frequency.h"
#include "gcp/util/Intensity.h"
#include "gcp/util/Temperature.h"

namespace gcp {
  namespace util {

    class Planck {
    public:

      /**
       * Constructor.
       */
      Planck();

      /**
       * Destructor.
       */
      virtual ~Planck();

      // Return the Planck intensity

      static Intensity IPlanck(Frequency nu, Temperature T);
      
      // Return the dimensionless Planck x-factor

      static double xPlanck(Frequency nu, Temperature T);

      // Return Planck dI/dT, in units of (Jy/sr) / K

      static double JyPerSrPerKPlanck(Frequency nu, Temperature T);

    private:
    }; // End class Planck

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PLANCK_H
