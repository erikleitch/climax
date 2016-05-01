// $Id: $

#ifndef GCP_UTIL_STRINGFACTORY_H
#define GCP_UTIL_STRINGFACTORY_H

/**
 * @file StringFactory.h
 * 
 * Tagged: Wed Feb  4 15:00:19 PST 2015
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include <string>
#include <vector>

namespace gcp {
  namespace util {

    class StringFactory {
    public:

      /**
       * Destructor.
       */
      virtual ~StringFactory();

      static char* getString(std::string str);

    private:

      // Constructor.

      StringFactory();

      static std::vector<std::string> strings_;

    }; // End class StringFactory

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_STRINGFACTORY_H
