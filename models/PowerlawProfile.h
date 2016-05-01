#ifndef GCP_MODELS_POWERLAWPROFILE_H
#define GCP_MODELS_POWERLAWPROFILE_H

#include "gcp/models/GenericRadiallySymmetric3DModel.h"

#include <vector>

//-----------------------------------------------------------------------
// Piecewise power-law model
//-----------------------------------------------------------------------

namespace gcp {
  namespace models {

    class PowerlawProfile : public GenericRadiallySymmetric3DModel {
    public:

      /**
       * Constructor.
       */
      PowerlawProfile();

      /**
       * Destructor.
       */
      virtual ~PowerlawProfile();

      void initialize();

      // Define different types of radial models

      double radialModel(unsigned type, double r, void* params);
#if 0
      virtual double radialRadioModelAM(double r, void* params);
#endif
      virtual double radialRadioModelEml(double r, void* params);
      virtual double radialXrayModel(double r, void* params);

      //------------------------------------------------------------
      // Return true if shape parameters are fixed
      //------------------------------------------------------------

      bool shapeParametersAreFixed();
      void checkSetup();
      void setParameter(std::string name, std::string val, std::string units);
      void setNumberOfSegments(unsigned n);
      double profileFn(double x, double alpha, double beta);

    public:

      //------------------------------------------------------------
      // Within each segment, the model is currently defined as:
      //    
      //                   1
      //     p(r) = ---------------
      //             /        a \ b
      //             | fac + r  |
      //             \          /
      //    
      // Radius is in units of thetaCore
      //
      // Exponents for index i hold from radius_[i] to radius_[i+1]
      // Exponents for index 0 hold from r=0 to radius_[1]
      // Last exponent holds to r = infinity
      //
      // The default value for fac is 1.0.  Although I have left fac
      // controllable, note that fac = 0.0 is
      // not integrable, and fac > 1.0 is not normalized to 1.0 at r = 0
      //------------------------------------------------------------

      gcp::util::VariableUnitQuantity fac_;

      double norm_;
      std::vector<double> piecewiseNorm_;

      std::vector<gcp::util::VariableUnitQuantity> radius_;
      std::vector<gcp::util::VariableUnitQuantity> alpha_;
      std::vector<gcp::util::VariableUnitQuantity> beta_;

      void binSearchForRadius(double rad, unsigned& iLow, unsigned& iHigh);

    }; // End class PowerlawProfile

  } // End namespace models
} // End namespace gcp


#endif // End #ifndef GCP_MODELS_POWERLAWPROFILE_H
