// $Id: $

#ifndef GCP_UTIL_SPECTRALTYPE_H
#define GCP_UTIL_SPECTRALTYPE_H

/**
 * @file SpectralType.h
 * 
 * Tagged: Mon Jul 16 11:33:24 PDT 2012
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/EnumeratedVariate.h"

namespace gcp {
  namespace util {

    class SpectralType : public EnumeratedVariate {
    public:

      enum {
	SPEC_NONE    = 0x0,
	SPEC_ALPHA   = 0x2,
	SPEC_SZ      = 0x4,
	SPEC_ITOH    = 0x8,
      };

      /**
       * Constructor.
       */
      SpectralType();

      /**
       * Destructor.
       */
      virtual ~SpectralType();

      // Overload the base-class string setVal method

      void initializeMaps();

    }; // End class SpectralType

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_SPECTRALTYPE_H
