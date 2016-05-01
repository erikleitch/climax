// $Id: $

#ifndef GCP_UTIL_JOINTGAUSSIANVARIATE_H
#define GCP_UTIL_JOINTGAUSSIANVARIATE_H

/**
 * @file JointGaussianVariate.h
 * 
 * Tagged: Wed Apr 24 10:43:02 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/ChisqVariate.h"

namespace gcp {
  namespace util {

    class JointGaussianVariate : public ChisqVariate {
    public:

      /**
       * Constructor.
       */
      JointGaussianVariate();

      /**
       * Destructor.
       */
      virtual ~JointGaussianVariate();

      Probability pdf();

    private:
    }; // End class JointGaussianVariate

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_JOINTGAUSSIANVARIATE_H
