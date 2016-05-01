// $Id: $

#ifndef GCP_UTIL_PGMODELMANAGER_H
#define GCP_UTIL_PGMODELMANAGER_H

/**
 * @file PgModelManager.h
 * 
 * Tagged: Mon Oct 21 10:20:11 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/pgutil/PgModel.h"
#include "gcp/pgutil/Trans.h"

#include <string>

namespace gcp {
  namespace util {

    class PgModelManager {
    public:

      /**
       * Constructor.
       */
      PgModelManager();

      /**
       * Destructor.
       */
      virtual ~PgModelManager();

      void display();
      void removeModel(float x, float y);
      void getModel(float xstart, float ystart, bool& read, char& key, std::string unit, Trans& trans);
      void printModels(std::string& unit, Trans& trans);
      
      void addModel(PgModel& model);

      std::vector<PgModel> models_;
      bool display_;

    }; // End class PgModelManager

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PGMODELMANAGER_H
