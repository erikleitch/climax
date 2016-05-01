// $Id: GenericLine.h,v 1.2 2012/05/11 20:21:40 eml Exp $

#ifndef GCP_MODELS_GENERICLINE_H
#define GCP_MODELS_GENERICLINE_H

/**
 * @file GenericLine.h
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

    class GenericLine : public gcp::util::Generic2DAngularModel {
    public:

      class GenericLineEvalData : public Generic2DAngularModel::EvalData {
      public:

	double lengthRad_;
	double widthRad_;
	double arg_;
      };

      /**
       * Constructor.
       */
      GenericLine();

      /**
       * Destructor.
       */
      virtual ~GenericLine();

      // Set the length

      void setLength(gcp::util::Angle length);
      void setWidth(gcp::util::Angle width);

      double envelope(unsigned type, double xRad, double yRad);

      // For multi-threaded execution

      double envelope(void* evalData, unsigned type, double xRad, double yRad);
      void*  allocateEvalData();
      void   initializeEvalData(void* evalData);
      
      // Check this model's setup for sense

      void checkSetup();

      gcp::util::PgModel pgModel();

    private:

      gcp::util::Angle length_;
      gcp::util::Angle width_;

    }; // End class GenericLine

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_GENERICLINE_H
