#ifndef GCP_UTIL_CONFORMABLEQUANTITY_H
#define GCP_UTIL_CONFORMABLEQUANTITY_H

/**
 * @file ConformableQuantity.h
 * 
 * Tagged: Wed Dec  1 11:44:54 PST 2004
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Exception.h"
#include "gcp/util/Variate.h"

#include <map>

namespace gcp {
  namespace util {
    
    // A pure interface for unit'ed quantities

    class ConformableQuantity : public Variate {
    public:
      
      /**
       * Constructor.
       */
      ConformableQuantity();
      
      /**
       * Destructor.
       */
      virtual ~ConformableQuantity();
      
      virtual void initialize() {};

      bool isFinite() {
	return finite_;
      }

#if 0
      ConformableQuantity(const ConformableQuantity& cq);
      ConformableQuantity(ConformableQuantity& cq);
      void operator=(const ConformableQuantity& cq);
      void operator=(ConformableQuantity& cq);
#endif

    public:

      void printUnits();
      std::string getUnitsString();

      Variate::Conversion addConversion(std::string units, double scaleFactor=1.0, double offset=0.0);
      void printConversions();

      virtual void setVal(double val, std::string units);

      virtual void setVal(std::string val);

      virtual double getVal(double val, std::string units);

      double getUnitVal() {
	return Variate::getUnitVal();
      }

      double getUnitVal(std::string units);
      Variate::Conversion* findConversion(std::string units);

      virtual void setUnits(std::string units);

      virtual bool unitlessAllowed() {
	return false;
      }

      void setFinite(bool finite) {
	finite_ = finite;
      }

      bool finite_;

      std::map<std::string, Variate::Conversion> unitConversions_;

      bool isValidUnit(std::string units);

      std::string typeName_;

    }; // End class ConformableQuantity
    
  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_CONFORMABLEQUANTITY_H
