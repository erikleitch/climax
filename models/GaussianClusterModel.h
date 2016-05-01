// $Id: GaussianClusterModel.h,v 1.5 2012/05/16 18:00:50 eml Exp $

#ifndef GCP_MODELS_GAUSSIANCLUSTERMODEL_H
#define GCP_MODELS_GAUSSIANCLUSTERMODEL_H

/**
 * @file GaussianClusterModel.h
 * 
 * Tagged: Fri Sep 17 15:47:31 PDT 2010
 * 
 * @version: $Revision: 1.5 $, $Date: 2012/05/16 18:00:50 $
 * 
 * @author tcsh: Erik Leitch
 */
#include "gcp/fftutil/Image.h"
#include "gcp/fftutil/Model.h"

#include "gcp/models/ClusterModel.h"

#include "gcp/util/Temperature.h"
#include "gcp/util/Frequency.h"

namespace gcp {
  namespace models {

    class GaussianClusterModel : public ClusterModel {
      //class GaussianClusterModel : public gcp::util::Generic2DAngularModel {
    public:

      /**
       * Constructor.
       */
      GaussianClusterModel();

      /**
       * Destructor.
       */
      virtual ~GaussianClusterModel();

      // Set the major/minor axis widths

      void setMajSigma(gcp::util::Angle majSigma);
      void setMinSigma(gcp::util::Angle minSigma);
      void setSigma(gcp::util::Angle sigma);

      void setMajFwhm(gcp::util::Angle majFwhm);
      void setMinFwhm(gcp::util::Angle minFwhm);
      void setFwhm(gcp::util::Angle fwhm);

      double radioEnvelope(double xRad, double yRad);

      void fillSzImage(gcp::util::Image& image, gcp::util::Frequency& frequency);

    public:
      
      gcp::util::Angle majSigma_;
      gcp::util::Angle minSigma_;

    }; // End class GaussianClusterModel
    
  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_GAUSSIANCLUSTERMODEL_H
