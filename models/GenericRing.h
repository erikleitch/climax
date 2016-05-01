// $Id: GenericRing.h,v 1.2 2012/05/11 20:21:40 eml Exp $

#ifndef GCP_MODELS_GENERICRING_H
#define GCP_MODELS_GENERICRING_H

/**
 * @file GenericRing.h
 * 
 * Tagged: Mon May  7 17:27:29 PDT 2012
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/11 20:21:40 $
 * 
 * @author Erik Leitch
 */
#include "gcp/fftutil/Generic2DAngularModel.h"
#include "gcp/fftutil/Image.h"

#include "gcp/pgutil/PgModel.h"

#include "gcp/util/Angle.h"

namespace gcp {
  namespace models {

    class GenericRing : public gcp::util::Generic2DAngularModel {
    public:

      class GenericRingEvalData : public Generic2DAngularModel::EvalData {
      public:

	double majRad_;
	double minRad_;
	double width_;
      };

      /**
       * Constructor.
       */
      GenericRing();

      /**
       * Destructor.
       */
      virtual ~GenericRing();

      // Set the major/minor axis widths

      void setRadius(gcp::util::Angle radius);
      void setWidth(gcp::util::Angle width);
      void setAxialRatio(double ratio);

      double envelope(unsigned type, double xRad, double yRad);

      // For multi-threaded execution

      double envelope(void* evalData, unsigned type, double xRad, double yRad);
      void*  allocateEvalData();
      void   initializeEvalData(void* evalData);
      
      // Check this model's setup for sense

      void checkSetup();

      gcp::util::PgModel pgModel();

    private:

      gcp::util::VariableUnitQuantity axialRatio_;
      gcp::util::Angle radius_;
      gcp::util::Angle width_;

    }; // End class GenericRing

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_GENERICRING_H
