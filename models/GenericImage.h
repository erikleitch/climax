// $Id: $

#ifndef GCP_MODELS_GENERICIMAGE_H
#define GCP_MODELS_GENERICIMAGE_H

/**
 * @file GenericImage.h
 * 
 * Tagged: Thu Nov  1 18:51:32 PDT 2012
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/fftutil/Generic2DAngularModel.h"
#include "gcp/fftutil/Image.h"
#include "gcp/fftutil/ImageManager.h"

namespace gcp {
  namespace models {

    class GenericImage : public gcp::util::Generic2DAngularModel {
    public:

      /**
       * Constructor.
       */
      GenericImage();

      /**
       * Destructor.
       */
      virtual ~GenericImage();

      // Evaluate the dimensionless shape of this profile

      double envelope(double xRad, double yRad);
      double radioEnvelope(double xRad, double yRad);
      double xrayImageEnvelope(double xRad, double yRad);
      double genericEnvelope(double xRad, double yRad);

      void initializeImage();

      void display();
      void difmapDisplay();

      void checkSetup();
      void displayIfRequested();

      gcp::util::Image image_;
      gcp::util::ImageManager imageManager_;

    private:

      // Utility objects used for calculation

      gcp::util::Angle xOffCalc_;
      gcp::util::Angle yOffCalc_;

      bool imageInitialized_;
      bool near_;
      bool setupChecked_;

    }; // End class GenericImage

  } // End namespace 
} // End namespace gcp



#endif // End #ifndef GCP__GENERICIMAGE_H
