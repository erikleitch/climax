// $Id: MiriadReader.h,v 1.1 2011/01/20 19:40:44 eml Exp $

#ifndef GCP_UTIL_MIRIADREADER_H
#define GCP_UTIL_MIRIADREADER_H

/**
 * @file MiriadReader.h
 * 
 * Tagged: Mon Nov 22 17:22:22 PST 2010
 * 
 * @version: $Revision: 1.1 $, $Date: 2011/01/20 19:40:44 $
 * 
 * @author tcsh: Erik Leitch
 */
#include <string>

namespace gcp {
  namespace util {

    class MiriadReader {
    public:

      /**
       * Constructor.
       */
      MiriadReader();

      /**
       * Destructor.
       */
      virtual ~MiriadReader();

      void readMiriadDir(std::string dir);

    private:
    }; // End class MiriadReader

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_MIRIADREADER_H
