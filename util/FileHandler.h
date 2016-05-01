// $Id: $

#ifndef GCP_UTIL_FILEHANDLER_H
#define GCP_UTIL_FILEHANDLER_H

/**
 * @file FileHandler.h
 * 
 * Tagged: Thu Oct  3 08:24:19 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include <string>

namespace gcp {
  namespace util {

    class FileHandler {
    public:

      /**
       * Constructor.
       */
      FileHandler();

      /**
       * Destructor.
       */
      virtual ~FileHandler();

      static bool fileExists(std::string name);

    private:
    }; // End class FileHandler

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_FILEHANDLER_H
