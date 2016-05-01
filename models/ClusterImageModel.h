// $Id: $

#ifndef GCP_MODELS_CLUSTERIMAGEMODEL_H
#define GCP_MODELS_CLUSTERIMAGEMODEL_H

/**
 * @file ClusterImageModel.h
 * 
 * Tagged: Mon Feb  9 09:34:50 PST 2015
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include "fftutil/ImageManager.h"

#include "models/ClusterModel.h"

namespace gcp {
  namespace models {

    class ClusterImageModel : public ClusterModel {
    public:

      /**
       * Constructor.
       */
      ClusterImageModel();

      /**
       * Destructor.
       */
      virtual ~ClusterImageModel();

    private:

      double radioImageEnvelope(double xRad, double yRad);
      double xrayImageEnvelope(double xRad, double yRad);

      gcp::util::ImageManager imageManager_;

    }; // End class ClusterImageModel

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_CLUSTERIMAGEMODEL_H
