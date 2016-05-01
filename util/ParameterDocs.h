// $Id: $

#ifndef GCP_UTIL_PARAMETERDOCS_H
#define GCP_UTIL_PARAMETERDOCS_H

/**
 * @file ParameterDocs.h
 * 
 * Tagged: Wed Jun 12 09:54:45 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/ParameterManager.h"

namespace gcp {
  namespace util {

    class ParameterDocs : public ParameterManager {
    public:

      /**
       * Constructor.
       */
      ParameterDocs();

      /**
       * Destructor.
       */
      virtual ~ParameterDocs();

    private:
    }; // End class ParameterDocs

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PARAMETERDOCS_H
