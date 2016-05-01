// $Id: $

#ifndef GCP_MODELS_GENERICRADIALLYSYMMETRIC3DMODEL_H
#define GCP_MODELS_GENERICRADIALLYSYMMETRIC3DMODEL_H

/**
 * @file GenericRadiallySymmetric3DModel.h
 * 
 * Tagged: Mon Jul 23 14:17:45 PDT 2012
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/models/ClusterModel.h"
#include "gcp/models/CosmologyModel.h"

#include "gcp/util/QuadraticInterpolatorNormal.h"

#include "gsl/gsl_integration.h"

#define GSL_HANDLER_FN(fn) void (fn)(const char* reason, const char* file, int line, int gsl_errno)

namespace gcp {
  namespace models {

    class GenericRadiallySymmetric3DModel : public ClusterModel {
    public:

      /**
       * Constructor.
       */
      GenericRadiallySymmetric3DModel();

      /**
       * Destructor.
       */
      virtual ~GenericRadiallySymmetric3DModel();

      void initialize();
      void setGslLimit(size_t limit);
      void initializeGslMembers();

      //------------------------------------------------------------
      // Virtual interface for inheritors to define
      //------------------------------------------------------------

      // The (normalized to unity) 3D radial model

      virtual double radialModel(unsigned type, double x, void* params);

      // True if shape parameters of this model are fixed

      virtual bool shapeParametersAreFixed();

      //------------------------------------------------------------
      // End virtual interface
      //------------------------------------------------------------

      // Evaluate the line integral of the radial model at the
      // specified dimensionless cylindrical radius

      double lineIntegral(unsigned type, double xSky, void* params=0);
      double volumeIntegral(unsigned type, double xSky, void* params=0);

      // Overloaded envelope functions from the Generic2DAngularModel base-class

      double genericEnvelope(double xRad, double yRad);
      double radioEnvelope(double xRad, double yRad);
      double xrayImageEnvelope(double xRad, double yRad);

      // Overloaded sample() function from the base-class.  Sets a
      // flag whenever a new sample is generated

      void sample();

      // Overloaded fillImage() function from the base-class.
      // Interpolates the line integral the first time it's called after a new
      // sample is generated.

      virtual void fillImage(unsigned type,  gcp::util::Image& image,           
			     void* params=0);

      virtual void debugPrint();

      // Return the integral of this shape function, out to the
      // specified radius

      virtual gcp::util::SolidAngle solidAngleIntegral(unsigned type, gcp::util::Angle radius);
      virtual void solidAngleIntegral(unsigned type, gcp::util::Angle& radius, gcp::util::SolidAngle& sa);
      double volumeIntegral(unsigned type, gcp::util::Angle& radius);

      virtual void initializeDerivedVariates(Model* caller);

      static VAR_DERIVE_FN(deriveRadioLineIntegralNormalization);
      static VAR_DERIVE_FN(deriveXrayLineIntegralNormalization);
      static VAR_DERIVE_FN(deriveVolumeIntegralFactor);
      static VAR_DERIVE_FN(deriveYsph);
      static VAR_DERIVE_FN(deriveMtot);
      static VAR_DERIVE_FN(deriveMgas);
      static VAR_DERIVE_FN(derivePe0);
      static VAR_DERIVE_FN(deriveScaleFactorY);
      static VAR_DERIVE_FN(deriveScaleFactorMass);
      static VAR_DERIVE_FN(deriveScaleFactorArea);
      static VAR_DERIVE_FN(deriveScaleFactorPressure);
      static VAR_DERIVE_FN(deriveInnerRadius);
      static VAR_DERIVE_FN(deriveOuterRadius);
      static VAR_DERIVE_FN(deriveThetaInner);
      static VAR_DERIVE_FN(deriveThetaOuter);
      static VAR_DERIVE_FN(deriveRestMassThermalEnergyRatio);

      void fillLineIntegralNormalization();
      void fillVolumeIntegralFactor();
      void fillYsph();
      void fillMtot();
      void fillMgas();
      void fillPe0();
      void fillScaleFactorY();
      void fillScaleFactorMass();
      void fillScaleFactorArea();
      void fillScaleFactorPressure();
      void fillInnerRadius();
      void fillOuterRadius();
      void fillThetaInner();
      void fillThetaOuter();
      void fillRestMassThermalEnergyRatio();

      void computeLineIntegralNormalization(unsigned type);
      double computeVolumeIntegralFactor(gcp::util::Angle& thetaInner, gcp::util::Angle& thetaOuter);

      //------------------------------------------------------------
      // Limiting x at which we will truncate the line integral
      //------------------------------------------------------------

      double xLim_;

      //------------------------------------------------------------
      // Members needed for derived variate computation
      //------------------------------------------------------------

      gcp::util::HubbleConstant H_;
      gcp::util::HubbleConstant H0_;

      gcp::util::VariableUnitQuantity scaleFactorY_;

      gcp::util::Mass scaleFactorMass_;
      gcp::util::Pressure scaleFactorPressure_;
      gcp::util::Area yInt_;
      gcp::util::Area scaleFactorArea_;
      gcp::util::Area dA2_;

      gcp::util::VariableUnitQuantity volumeIntegralFactor_;
      gcp::util::VariableUnitQuantity restMassThermalEnergyRatio_;

      gcp::util::Length innerRadius_;
      gcp::util::Length outerRadius_;
      gcp::util::Angle thetaInner_;
      gcp::util::Angle thetaOuter_;
      gcp::util::Mass mass_;
      gcp::util::Pressure pressure_;

      //------------------------------------------------------------
      // Derived variates
      //------------------------------------------------------------

      gcp::util::Pressure pe0_;
      gcp::util::Area     ySph_;
      gcp::util::Mass     mGas_;
      gcp::util::Mass     mTot_;

    public:

      //------------------------------------------------------------
      // Calculate values used for subsequent interpolation of the line integral
      //------------------------------------------------------------

      void calculateInterpolationValues(unsigned type, gcp::util::Image& image);
      void calculateInterpolationValuesForAllScaleRadii(unsigned type, gcp::util::Image& image);
      void calculateInterpolationValuesForCurrentScaleRadius(unsigned type, gcp::util::Image& image);

      //------------------------------------------------------------
      // Interpolate the line integral of the radialModel()
      //------------------------------------------------------------

      double interpolateLineIntegral(unsigned typ, double x);
      double interpolateArealIntegral(unsigned type, double x);
      double interpolateVolumeIntegral(unsigned type, double x);

      //------------------------------------------------------------
      // A flag set to true whenever a new sample has been generated
      //------------------------------------------------------------

      std::map<unsigned, bool> newSample_;

      //------------------------------------------------------------
      // True if this type of model must be reinterpolated for each new sample
      //------------------------------------------------------------

      std::map<unsigned, bool> needsRecomputing_;

    public:

      //------------------------------------------------------------
      // Evaluate interpolation values for the passed type, for the
      // current nInterp_ and deltaInterp_ values
      //------------------------------------------------------------

      virtual void calculateInterpolationValues(unsigned type);
      virtual void calculateAllInterpolationValues(unsigned type);
      virtual void calculateVolumeInterpolationValues(unsigned type);

      //------------------------------------------------------------
      // Evaluate the kernel of the line and volume integrals at the
      // specified x = theta/thetaCore_
      //------------------------------------------------------------

      static double evaluateLineIntegralKernel(double x, void* params=0);
      static double evaluateLineIntegralKernel2(double x, void* params=0);
      static double evaluateVolumeIntegralKernel(double x, void* params=0);

      void integratedYsph(gcp::util::Area& yInt, double& volumeIntegralFactor);
      void integratedGasMass(gcp::util::Mass& mass, double& volumeIntegralFactor);
      void integratedGasMass(gcp::util::Mass& mass, gcp::util::Angle& theta);
      void integratedTotalMass(gcp::util::Mass& mass, double& volumeIntegralFactor);
      void integratedTotalMass(gcp::util::Mass& mass, gcp::util::Angle& theta);
      void pressure(gcp::util::Pressure& pressure);

      //------------------------------------------------------------
      // Parameters needed to pass to the GSL integration routines
      //------------------------------------------------------------

      gsl_integration_workspace* gslWork_;
      size_t gslLimit_;

      gsl_function gslLineIntegralFn_;
      gsl_function gslVolumeIntegralFn_;

      void* params_;
      unsigned currentIntegrationType_;

      //------------------------------------------------------------
      // The (dimensionless) cylindrical radius at which to evaluate
      // the line integral, ie., xSky_ = theta_sky / theta_c
      //------------------------------------------------------------

      double xSky_;

      static GSL_HANDLER_FN(intErrHandler);

    protected:
      
      virtual void checkSetup();
      unsigned getMaxNSigma();

    public:

      //------------------------------------------------------------
      // Arrays used for interpolation and integration
      //------------------------------------------------------------

      std::map<unsigned, std::vector<double> > interpolatedLineIntegral_;
      std::map<unsigned, std::vector<double> > arealIntegral_;
      std::map<unsigned, std::vector<double> > volumeIntegral_;
      std::map<unsigned, gcp::util::VariableUnitQuantity> lineIntegralNormalization_;

      double   deltaInterp_;
      unsigned nInterp_;

      gcp::util::QuadraticInterpolatorNormal quadInterp_;

    }; // End class GenericRadiallySymmetric3DModel

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_GENERICRADIALLYSYMMETRIC3DMODEL_H
