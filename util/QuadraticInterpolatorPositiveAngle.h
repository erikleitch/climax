#ifndef GCP_UTIL_QUADRATICINTERPOLATORPOSITIVEANGLE_H
#define GCP_UTIL_QUADRATICINTERPOLATORPOSITIVEANGLE_H

/**
 * @file QuadraticInterpolatorPositiveAngle.h
 * 
 * Tagged: Tue Mar 16 17:02:03 UTC 2004
 * 
 * @author Erik Leitch
 */
#include "gcp/util/QuadraticInterpolator.h"

namespace gcp {
  namespace util {
    
    /**
     * Class for interpolating signed (0 <= v < 2.pi) angles
     */
    class QuadraticInterpolatorPositiveAngle : public QuadraticInterpolator {
    public:

      QuadraticInterpolatorPositiveAngle(double emptyValue, bool relativeX=false);

    private:

      double fixAngle(double angle);

    }; // End class QuadraticInterpolatorPositiveAngle
    
  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_QUADRATICINTERPOLATORPOSITIVEANGLE_H
