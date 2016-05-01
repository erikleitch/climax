// $Id: $

#ifndef GCP_UTIL_OBSMANAGER_H
#define GCP_UTIL_OBSMANAGER_H

/**
 * @file ObsManager.h
 * 
 * Tagged: Thu Jul 26 11:03:42 PDT 2012
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "ObsInfo.h"

namespace gcp {
  namespace util {

    class ObsManager {
    public:

      /**
       * Constructor.
       */
      ObsManager();

      /**
       * Destructor.
       */
      virtual ~ObsManager();

      //------------------------------------------------------------
      // Method to allocate an obsinfo object
      //------------------------------------------------------------

      gcp::util::ObsInfo* addObs(std::string modelName);

      //------------------------------------------------------------
      // Method to return an allocated obsinfo object by name
      //------------------------------------------------------------

      gcp::util::ObsInfo* getObs(std::string name);

      //------------------------------------------------------------
      // Internal map of known obs info objects
      //------------------------------------------------------------

      std::map<std::string, gcp::util::ObsInfo*> obsMap_;

    }; // End class ObsManager

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_OBSMANAGER_H
