#ifndef GCP_UTIL_QUADRATICINTERPOLATORNORMAL_H
#define GCP_UTIL_QUADRATICINTERPOLATORNORMAL_H

/**
 * @file QuadraticInterpolatorNormal.h
 * 
 * Tagged: Tue Mar 16 17:16:24 UTC 2004
 * 
 * @author Erik Leitch
 */
#include "gcp/util/QuadraticInterpolator.h"

namespace gcp {
  namespace util {
    
    /**
     * Class for interpolating normal (non-angle) ordinates
     */
    class QuadraticInterpolatorNormal : public QuadraticInterpolator {
    public:
      
      /**
       * Constructor.
       */
      QuadraticInterpolatorNormal(double emptyValue=0.0, bool relativeX=false);

    }; // End class QuadraticInterpolatorNormal
    
  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_QUADRATICINTERPOLATORNORMAL_H
