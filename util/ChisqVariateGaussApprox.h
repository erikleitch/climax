// $Id: $

#ifndef GCP_UTIL_CHISQVARIATEGAUSSAPPROX_H
#define GCP_UTIL_CHISQVARIATEGAUSSAPPROX_H

/**
 * @file ChisqVariateGaussApprox.h
 * 
 * Tagged: Wed May  1 10:48:59 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/ChisqVariate.h"

namespace gcp {
  namespace util {

    class ChisqVariateGaussApprox : public ChisqVariate {
    public:

      /**
       * Constructor.
       */
      ChisqVariateGaussApprox();

      /**
       * Destructor.
       */
      virtual ~ChisqVariateGaussApprox();

      Probability pdf();

    private:
    }; // End class ChisqVariateGaussApprox

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_CHISQVARIATEGAUSSAPPROX_H
