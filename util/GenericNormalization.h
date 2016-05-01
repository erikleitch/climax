// $Id: $

#ifndef GCP_UTIL_GENERICNORMALIZATION_H
#define GCP_UTIL_GENERICNORMALIZATION_H

/**
 * @file GenericNormalization.h
 * 
 * Tagged: Thu Oct  3 13:50:21 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "util/VariableUnitQuantity.h"
#include "util/Unit.h"

namespace gcp {
  namespace util {

    class GenericNormalization : public VariableUnitQuantity {
    public:

      /**
       * Constructor.
       */
      GenericNormalization();

      /**
       * Destructor.
       */
      virtual ~GenericNormalization();

      void operator=(GenericNormalization& gn);
      void operator=(const GenericNormalization& gn);
      GenericNormalization(GenericNormalization& gn);
      GenericNormalization(const GenericNormalization& gn);

      void setVal(double val, std::string units);
      double getVal(std::string units);
      double getVal(gcp::util::Unit::Units units);

    private:

      ConformableQuantity* currentQuantity_;
      std::vector<ConformableQuantity*> knownQuantities_;

    }; // End class GenericNormalization

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_GENERICNORMALIZATION_H
