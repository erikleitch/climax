// $Id: Stokes.h,v 1.1.1.1 2010/07/13 17:56:32 eml Exp $

#ifndef GCP_UTIL_STOKES_H
#define GCP_UTIL_STOKES_H

/**
 * @file Stokes.h
 * 
 * Tagged: Wed Aug 20 16:28:18 PDT 2008
 * 
 * @version: $Revision: 1.1.1.1 $, $Date: 2010/07/13 17:56:32 $
 * 
 * @author Erik Leitch.
 */
#include <iostream>

namespace gcp {
  namespace util {

    class Stokes {
    public:

      /**
       * Constructor.
       */
      Stokes();

      /**
       * Destructor.
       */
      virtual ~Stokes();

      enum Param {
	STOKES_NONE =  0x0,
	STOKES_I    =  0x1,
	STOKES_Q    =  0x2,
	STOKES_U    =  0x4,
	STOKES_T    =  0x8,
	STOKES_E    = 0x10,
	STOKES_B    = 0x20,
      };

      friend std::ostream& operator<<(std::ostream& os, const Stokes::Param& param);

    }; // End class Stokes

    std::ostream& operator<<(std::ostream& os, const Stokes::Param& param);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_STOKES_H
