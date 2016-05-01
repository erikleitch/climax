#ifndef GCP_UTIL_INTEGRATOR_H
#define GCP_UTIL_INTEGRATOR_H

/**
 * @file Integrator.h
 * 
 * Tagged: Wed Oct 23 16:47:33 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gsl/gsl_integration.h"

#define INT_FN(fn) double (fn)(double x, void* params)

namespace gcp {
  namespace util {

    class Integrator {
    public:

      /**
       * Constructor.
       */
      Integrator();

      /**
       * Destructor.
       */
      virtual ~Integrator();

      void initializeGslMembers();

      void installGslFn(INT_FN(*gslIntFn), void* params);

      double integrateFromLowlimToHighlim(double lowLim, double highLim);
      double integrateFromNegInftyToPosInfty();
      double integrateFromLowlimToInfty(double lowLim);
      double integrateFromNegInftyToHighlim(double highLim);

      double integrateFromLowlimToHighlim(INT_FN(*gslIntFn), void* params, double lowLim, double highLim);
      double integrateFromNegInftyToPosInfty(INT_FN(*gslIntFn), void* params);
      double integrateFromLowlimToInfty(INT_FN(*gslIntFn), void* params, double lowLim);
      double integrateFromNegInftyToHighlim(INT_FN(*gslIntFn), void* params, double highLim);

      //------------------------------------------------------------
      // Parameters needed to pass to the GSL integration routines
      //------------------------------------------------------------

      double gslEpsAbs_;
      double gslEpsRel_;

      gsl_integration_workspace* gslWork_;
      size_t gslLimit_;

      gsl_function gslFn_;

    }; // End class Integrator

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_INTEGRATOR_H
