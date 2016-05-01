// $Id: $

#ifndef GCP_MODELS_GNFWMODEL_H
#define GCP_MODELS_GNFWMODEL_H

/**
 * @file GnfwModel.h
 * 
 * Tagged: Mon Jul 23 21:25:09 PDT 2012
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/models/GenericRadiallySymmetric3DModel.h"

//-----------------------------------------------------------------------
// Generalized form of the pressure model from Nagai et al 2007, ApJ
// 668, 1
//-----------------------------------------------------------------------

namespace gcp {
  namespace models {

    class GnfwModel : public GenericRadiallySymmetric3DModel {
    public:

      /**
       * Constructor.
       */
      GnfwModel();

      /**
       * Destructor.
       */
      virtual ~GnfwModel();

      void initialize();

      // Define different types of radial models

      double radialModel(unsigned type, double r, void* params);
      double radialGenericModel(double r, void* params);
      virtual double radialRadioModel(double r, void* params);
      virtual double radialXrayModel(double r, void* params);

      //------------------------------------------------------------
      // Return true if shape parameters are fixed
      //------------------------------------------------------------

      bool shapeParametersAreFixed();

      void checkSetup();

      double gnfwX(double x);

      //------------------------------------------------------------
      // Methods for deriving variates
      //------------------------------------------------------------

      static VAR_DERIVE_FN(deriveMgas500);
      static VAR_DERIVE_FN(deriveMtot500);
      static VAR_DERIVE_FN(deriveR500);
      static VAR_DERIVE_FN(deriveM500);
      static VAR_DERIVE_FN(deriveT500);
      static VAR_DERIVE_FN(deriveYsph500);
      static VAR_DERIVE_FN(deriveThetaCore);
      static VAR_DERIVE_FN(deriveRadioNormalization);
      static VAR_DERIVE_FN(deriveVolumeIntegralFactorR500);
      static VAR_DERIVE_FN(deriveScaleFactorMassT500);

      void fillVolumeIntegralFactorR500();
      void fillScaleFactorMassT500();
      void fillMgas500();
      void fillMtot500();
      void fillR500();
      void fillYsph500();
      void fillR500FromThetaCore();
      void fillR500FromRhoCrit();
      void fillT500();

      gcp::util::Temperature getTemperature(gcp::util::Mass mDelta, double delta);
      gcp::util::Length getRDeltaFromRhoCrit(gcp::util::Mass mDelta, double delta);

      virtual void fillM500();
      virtual void fillRadioNormalization();
      void fillThetaCore();

      gcp::util::Mass scaleFactorMassT500();
      void integratedGasMassT500(gcp::util::Mass& mass, double& volumeIntegralFactor);
      void integratedGasMassT500(gcp::util::Mass& mass, gcp::util::Angle& theta);
      void integratedTotalMassT500(gcp::util::Mass& mass, double& volumeIntegralFactor);
      void integratedTotalMassT500(gcp::util::Mass& mass, gcp::util::Angle& theta);

      bool integratedMassRequired();
      void interpolateForThetaCore();
      void interpolateForRescale();

      void calculateThetaCoreFromRhoCrit();

      virtual double getGnfwNormalization();

    public:

      //------------------------------------------------------------
      // Variates common to this class of models
      //------------------------------------------------------------

      gcp::util::VariableUnitQuantity alpha_;
      gcp::util::VariableUnitQuantity beta_;
      gcp::util::VariableUnitQuantity gamma_;
      gcp::util::VariableUnitQuantity fac_;
      gcp::util::VariableUnitQuantity c_;
      gcp::util::VariableUnitQuantity p0_;

      gcp::util::VariableUnitQuantity rescale_;

      //------------------------------------------------------------
      // Derived variates common to this class of models
      //------------------------------------------------------------

      gcp::util::Length r500_;
      gcp::util::Mass   mGas500_;
      gcp::util::Mass   mTot500_;
      gcp::util::Area   ySph500_;
      gcp::util::Mass   m500_;
      gcp::util::Temperature t500_;
      
      //------------------------------------------------------------
      // Convenience variables for computation
      //------------------------------------------------------------

      gcp::util::Angle thetaR500_;
      gcp::util::VariableUnitQuantity volumeIntegralFactorR500_;
      gcp::util::Mass scaleFactorMassT500_;

    }; // End class GnfwModel

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_GNFWMODEL_H
