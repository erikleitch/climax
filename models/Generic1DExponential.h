// $Id: Generic1DExponential.h,v 1.4 2012/05/12 00:45:52 eml Exp $

#ifndef GCP_MODELS_GENERIC1DEXPONENTIAL_H
#define GCP_MODELS_GENERIC1DEXPONENTIAL_H

/**
 * @file Generic1DExponential.h
 * 
 * Tagged: Sun Apr 29 09:45:28 PDT 2012
 * 
 * @version: $Revision: 1.4 $, $Date: 2012/05/12 00:45:52 $
 * 
 * @author Erik Leitch
 */
#include "gcp/fftutil/Generic1DModel.h"

#include "gcp/util/Variate.h"

namespace gcp {
  namespace models {

    class Generic1DExponential : public gcp::util::Generic1DModel {
    public:

      class Generic1DExponentialEvalData {
      public:

	double dx_;
	double lr_;
	double ld_;
	double peak_;
	double xPeak_;
      };

      /**
       * Constructor.
       */
      Generic1DExponential();

      /**
       * Destructor.
       */
      virtual ~Generic1DExponential();

      void setPeak(double peak);
      void setXPeak(double xPeak);
      void setLambdaRise(double lambda);
      void setLambdaDecay(double lambda);

      double eval(double x);
      double eval(void* evalData, double x);

      void* allocateEvalData();
      void initializeEvalData(void* evalData);

      void checkSetup();

    private:

      gcp::util::Variate peak_;
      gcp::util::Variate xPeak_;
      gcp::util::Variate lambdaRise_;
      gcp::util::Variate lambdaDecay_;

    }; // End class Generic1DExponential

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_GENERIC1DEXPONENTIAL_H
