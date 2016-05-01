// $Id: $

#ifndef GCP_UTIL_GENERIC2DSINUSOID_H
#define GCP_UTIL_GENERIC2DSINUSOID_H

/**
 * @file Generic2DSinusoid.h
 * 
 * Tagged: Sun Jan 18 14:52:00 NZDT 2015
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include "gcp/fftutil/Generic2DAngularModel.h"
#include "gcp/fftutil/Image.h"

#include "gcp/pgutil/PgModel.h"

#include "gcp/util/Angle.h"

namespace gcp {
  namespace util {

    class Generic2DSinusoid : public gcp::util::Generic2DAngularModel {
    public:

      /**
       * Constructor.
       */
      Generic2DSinusoid();

      /**
       * Destructor.
       */
      virtual ~Generic2DSinusoid();

      double envelope(unsigned type, double xRad, double yRad);

      gcp::util::Angle wavelength_;
      gcp::util::Angle phase_;

    }; // End class Generic2DSinusoid

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_GENERIC2DSINUSOID_H
