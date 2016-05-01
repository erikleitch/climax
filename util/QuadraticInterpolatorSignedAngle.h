#ifndef GCP_UTIL_QUADRATICINTERPOLATORSIGNEDANGLE_H
#define GCP_UTIL_QUADRATICINTERPOLATORSIGNEDANGLE_H

/**
 * @file QuadraticInterpolatorSignedAngle.h
 * 
 * Tagged: Tue Mar 16 16:59:18 UTC 2004
 * 
 * @author Erik Leitch
 */
#include "gcp/util/QuadraticInterpolator.h"

namespace gcp {
  namespace util {
    
    /**
     * Class for interpolating signed (-pi <= v < pi) angles
     */
    class QuadraticInterpolatorSignedAngle : public QuadraticInterpolator {
    public:

      QuadraticInterpolatorSignedAngle(double emptyValue, bool relativeX=false);

    private:

      double fixAngle(double angle);

    }; // End class QuadraticInterpolatorSignedAngle
    
  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_QUADRATICINTERPOLATORSIGNEDANGLE_H
