// $Id: $

#ifndef GCP_UTIL_SCATTERING_H
#define GCP_UTIL_SCATTERING_H

/**
 * @file Scattering.h
 * 
 * Tagged: Thu Oct 24 11:39:40 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Integrator.h"
#include "gcp/util/Intensity.h"
#include "gcp/util/Temperature.h"

#define REDIST_FN(fn) double (fn)(double t, void* params)
#define ETE_FN(fn) double (fn)(void* params)

namespace gcp {
  namespace util {

    class Scattering {
    public:

      /**
       * Constructor.
       */
      Scattering();

      /**
       * Destructor.
       */
      virtual ~Scattering();

      class GridParameter {
      public:

	void load(int fd);
	void findBracketing(double val, int& i1, int& i2);
	double distance(double val, unsigned ind);

	unsigned npt_;
	double min_;
	double max_;
	double delta_;

	friend std::ostream& operator<<(std::ostream& os, const GridParameter& param);
      };

      class PowerlawGridder {
      public:

	PowerlawGridder() {};

	void initialize(std::string fileName);

	double interpolate(double alpha, double p1, double p2, double x);

	double getVal(int iAlpha, int iP1, int iP2, int iX);

	std::vector<GridParameter> parameters_;
	std::vector<std::vector<double> > vals_;
      };

      //------------------------------------------------------------
      // Integratable functions
      //------------------------------------------------------------

      // The mono-energetic scattering kernel

      static INT_FN(intScatteringKernel);

      // The scattered Planck spectrum kernel

      static INT_FN(intScatteredSpectrumKernel);

      // The equivalent thermal energy kernel

      static INT_FN(equivalentThermalEnergyKernel);

      // Momentum distributions

      static INT_FN(intThermalMomentumDistribution);
      static INT_FN(intPowerlawMomentumDistribution);
      static INT_FN(intThermalTailDistribution);

      static REDIST_FN(photonRedistributionFunctionFinite);
      static REDIST_FN(photonRedistributionFunctionInfinite);
      static REDIST_FN(photonRedistributionFunctionPowerlaw);

      static ETE_FN(equivalentThermalEnergyFunctionFinite);
      static ETE_FN(equivalentThermalEnergyFunctionInfinite);

      double photonRedistributionFunction(double t);
      double photonRedistributionFunctionMono(double t, double p);
      double powerlawPhotonRedistributionFunction(double t, double alpha, double p1, double p2);

      double thermalMomentumDistribution(double p, double eta, double norm);
      double powerlawMomentumDistribution(double p, double alpha, double p1, double p2);
      double thermalTailDistribution(double p, double eta, double norm, double alpha, double p1, double p2);

      void initializeThermalDistribution(Temperature& Te);
      void initializePowerlawDistribution(double alpha, double p1, double p2);
      void initializePowerlawDistribution2(double alpha, double p1, double p2);
      void initializeMonoEnergeticDistribution(double p);

      void initializeThermalTailDistribution(Temperature& Te, double alpha, double p1, double p2);

      double scatteredSpectralShape(double x);
      double planckSpectralShape(double x);
      double kompaneetsSpectralShape(double x);
      double h(double x);

      Intensity planck(Frequency& freq, Temperature& temp);
      Intensity scatteredPlanckSpectrum(Frequency& freq, Temperature& temp);

      double incompleteBetaBx(double a, double b, double x);
      double incompleteBetaIx(double a, double b, double x);

      double powerlawEvalFn(double t, double p, double alpha);

      void comptonYToDeltaI(Frequency& freq, Intensity& YtoI);
      void comptonYToDeltaT(Frequency& freq, Temperature& YtoT);

      double g(double x);
      double jmi(double x);

      void computePowerlawGrid(double alphaMin, double alphaMax, unsigned nAlpha,
			       double p1Min,    double p1Max,    unsigned nP1,
			       double p2Min,    double p2Max,    unsigned nP2,
			       double xMin,     double xMax,     unsigned nx,
			       std::string fileName);

      void loadPowerlawGrid(std::string fileName);

    public:

      INT_FN(*momentumDistribution_);
      REDIST_FN(*photonRedistributionFunction_);
      ETE_FN(*equivalentThermalEnergyFunction_);
      
      Temperature electronTemperature_;
      double eta_;
      double norm_;
      Integrator integratorDist_;
      Integrator integratorSpec_;
      double t_;
      double x_;
      double alpha_;
      double p1_;
      double p2_;
      double lowLim_;
      double highLim_;
      double equivalentThermalEnergyPerRestmass_;
      Intensity planckNormalization_;

      PowerlawGridder powerlawGridder_;

      bool debug_;
      void setDebug(bool debug);

    }; // End class Scattering

    std::ostream& operator<<(std::ostream& os, const Scattering::GridParameter& param);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_SCATTERING_H
