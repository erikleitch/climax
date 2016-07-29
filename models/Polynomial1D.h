// $Id: $

#ifndef GCP_MODELS_POLYNOMIAL1D_H
#define GCP_MODELS_POLYNOMIAL1D_H

/**
 * @file Polynomial1D.h
 * 
 * Tagged: Wed May  7 09:22:30 PDT 2014
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include "gcp/fftutil/Generic1DModel.h"

#include "gcp/util/VariableUnitQuantity.h"

namespace gcp {
  namespace models {

    class Polynomial1D : public gcp::util::Generic1DModel {
    public:

      /**
       * Constructor.
       */
      Polynomial1D();

      /**
       * Destructor.
       */
      virtual ~Polynomial1D();

      void setNumberOfCoefficients(unsigned n);
      void setParameter(std::string name, std::string val, std::string units, bool external);

      std::vector<gcp::util::VariableUnitQuantity> coeffs_;

      double eval(double x);

    }; // End class Polynomial1D

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_POLYNOMIAL1D_H
