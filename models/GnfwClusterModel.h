// $Id: GnfwClusterModel.h,v 1.2 2012/05/08 21:58:00 eml Exp $

#ifndef GCP_MODELS_GNFWCLUSTERMODEL_H
#define GCP_MODELS_GNFWCLUSTERMODEL_H

/**
 * @file GnfwClusterModel.h
 * 
 * Tagged: Tue Apr  3 14:11:51 PDT 2012
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/08 21:58:00 $
 * 
 * @author Erik Leitch
 */
#include "gcp/fftutil/Model.h"

#include "gcp/models/ClusterModel.h"

namespace gcp {
  namespace models {

    class GnfwClusterModel : ClusterModel {
    public:

      /**
       * Constructor.
       */
      GnfwClusterModel();

      /**
       * Destructor.
       */
      virtual ~GnfwClusterModel();

    }; // End class GnfwClusterModel

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_GNFWCLUSTERMODEL_H
