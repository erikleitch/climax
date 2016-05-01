// $Id: Generic2DGaussian.h,v 1.2 2012/05/11 20:21:40 eml Exp $

#ifndef GCP_MODELS_GENERIC2DGAUSSIAN_H
#define GCP_MODELS_GENERIC2DGAUSSIAN_H

/**
 * @file Generic2DGaussian.h
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

    class Generic2DGaussian : public gcp::util::Generic2DAngularModel {
    public:

      class Generic2DGaussianEvalData : public Generic2DAngularModel::EvalData {
      public:

	double majSigRad_;
	double minSigRad_;
	double arg_;
      };

      /**
       * Constructor.
       */
      Generic2DGaussian();

      /**
       * Destructor.
       */
      virtual ~Generic2DGaussian();

      // Set the major/minor axis widths

      void setMajSigma(gcp::util::Angle majSigma);
      void setAxialRatio(double ratio);
      void setSigma(gcp::util::Angle sigma);

      void setMajFwhm(gcp::util::Angle majFwhm);
      void setFwhm(gcp::util::Angle fwhm);

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
      gcp::util::Angle majSigma_;

    }; // End class Generic2DGaussian

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_GENERIC2DGAUSSIAN_H
