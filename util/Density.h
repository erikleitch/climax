// $Id: $

#ifndef GCP_UTIL_DENSITY_H
#define GCP_UTIL_DENSITY_H

/**
 * @file Density.h
 * 
 * Tagged: Wed Jun 19 11:36:50 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/ConformableQuantity.h"
#include "gcp/util/Volume.h"

namespace gcp {
  namespace util {

    class Mass;

    class Density : public ConformableQuantity {
    public:

      /**
       * Constructor.
       */
      Density();

      /**
       * Destructor.
       */
      virtual ~Density();

      // Set the density of this object

      void setGPerCm3(double gPerCm3);

      // Get the density of this object

      double gPerCm3();

      // A density times a volume is a mass

      Mass operator*(const Volume& vol);

      void operator=(const Density& density);
      void operator=(Density& density);

    private:
    }; // End class Density

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_DENSITY_H
