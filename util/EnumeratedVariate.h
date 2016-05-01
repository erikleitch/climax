// $Id: $

#ifndef GCP_UTIL_ENUMERATEDVARIATE_H
#define GCP_UTIL_ENUMERATEDVARIATE_H

/**
 * @file EnumeratedVariate.h
 * 
 * Tagged: Thu Mar 28 10:35:59 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Variate.h"

#include <map>
#include <sstream>

namespace gcp {
  namespace util {

    class EnumeratedVariate : public Variate {
    public:

      /**
       * Constructor.
       */
      EnumeratedVariate();

      /**
       * Destructor.
       */
      virtual ~EnumeratedVariate();

      // Overloaded base-class methods

      void setVal(std::string val);
      std::string getStringVal();

      // Other methods

      void listValidOptions(std::ostringstream& os);
      virtual void initializeMaps() = 0;

      unsigned type_;

    protected:

      std::map<unsigned int, std::string>  idToNameMap_;
      std::map<unsigned int, std::string>  explMap_;
      std::map<std::string,  unsigned int> nameToIdMap_;

    }; // End class EnumeratedVariate

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_ENUMERATEDVARIATE_H
