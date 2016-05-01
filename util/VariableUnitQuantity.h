// $Id: VariableUnitQuantity.h,v 1.2 2012/05/30 16:53:11 eml Exp $

#ifndef GCP_UTIL_VARIABLEUNITQUANTITY_H
#define GCP_UTIL_VARIABLEUNITQUANTITY_H

/**
 * @file VariableUnitQuantity.h
 * 
 * Tagged: Tue May 15 16:22:29 PDT 2012
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/30 16:53:11 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/ConformableQuantity.h"

namespace gcp {
  namespace util {

    class VariableUnitQuantity : public ConformableQuantity {
    public:

      /**
       * Constructor.
       */
      VariableUnitQuantity();

      /**
       * Destructor.
       */
      virtual ~VariableUnitQuantity();

      void allowUnitless(bool allow);

      bool unitlessAllowed();

      virtual void setVal(double val, std::string units);
      void setVal(std::string val);

      virtual double getVal(double val, std::string units);
      
      virtual void setUnits(std::string units);

      void operator=(const VariableUnitQuantity& vuq);
      void operator=(VariableUnitQuantity& vuq);

    private:

      bool stringIsEmpty(std::string& str);

      bool unitlessAllowed_;

    }; // End class VariableUnitQuantity

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_VARIABLEUNITQUANTITY_H
