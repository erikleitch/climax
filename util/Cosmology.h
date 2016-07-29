// $Id: Cosmology.h,v 1.1 2012/05/02 23:44:51 eml Exp $

#ifndef GCP_UTIL_COSMOLOGY_H
#define GCP_UTIL_COSMOLOGY_H

/**
 * @file Cosmology.h
 * 
 * Tagged: Thu Apr  5 08:51:54 PDT 2012
 * 
 * @version: $Revision: 1.1 $, $Date: 2012/05/02 23:44:51 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Density.h"
#include "gcp/util/HubbleConstant.h"
#include "gcp/util/Length.h"
#include "gcp/util/ParameterManager.h"

#include "gsl/gsl_integration.h"

namespace gcp {
  namespace util {

    class Cosmology : public ParameterManager {
    public:

      /**
       * Constructor.
       */
      Cosmology();

      /**
       * Destructor.
       */
      virtual ~Cosmology();

      //------------------------------------------------------------
      // Set cosmological parameters
      //------------------------------------------------------------

      void setH0(HubbleConstant H0);
      void setOmegaM(double omegaM);
      void setOmegaL(double omegaL);
      void setOmegaK(double omegaK);
      void setRedshift(double z);
      
      //------------------------------------------------------------
      // Query cosmological parameters
      //------------------------------------------------------------
      
      double dimensionlessHubbleConstant();
      double dimensionlessHubbleConstant(double z);
      double h();
      double h(double z);

      HubbleConstant H0();
      HubbleConstant H(double z);
      HubbleConstant H();
      bool isFlat();
      double omegaK();
      double omegaL();
      double omegaM();
      double redshift();

      // Return dH -- the Hubble distance

      Length hubbleDistance();

      // Dimensionless E(z) and E^2(z) for given cosmology

      double E();
      double E2();
      Length comovingDistance();
      Length transverseComovingDistance();
      Length angularDiameterDistance();
      Length luminosityDistance();
      Length lightTravelDistance();

      double E(double z);
      double E2(double z);
      Length comovingDistance(double z);
      Length transverseComovingDistance(double z);
      Length angularDiameterDistance(double z);
      Length luminosityDistance(double z);
      Length lightTravelDistance(double z);

      Density criticalDensity();
      Density criticalDensity(double z);

      void setParameter(std::string name, std::string val, std::string units=" ", bool external=true);

      void checkParameter(std::string par);
      void checkParametersOr(std::string par1, std::string par2);
      void checkParametersOr(std::string par1, std::string par2, std::string par3);
      void checkParametersAnd(std::string par1, std::string par2);

    private:

      //------------------------------------------------------------
      // Fundamental parameters of this class
      //------------------------------------------------------------

      HubbleConstant H0_;
      double z_;
      double omegaM_;
      double omegaL_;

      //------------------------------------------------------------
      // Needed for numerical integration
      //------------------------------------------------------------

      void setGslLimit(size_t limit);
      void initializeGslMembers();

      //------------------------------------------------------------
      // Parameters needed to pass to the GSL integration routines
      //------------------------------------------------------------

      static const double eps_;
      gsl_integration_workspace* gslWork_;
      size_t gslLimit_;
      double gslEpsAbs_;
      double gslEpsRel_;
      gsl_function gslComovingFn_;
      gsl_function gslLightTravelFn_;

      static double evaluateComovingDistanceKernel(double z, void* params=0);
      static double evaluateLightTravelDistanceKernel(double z, void* params=0);
      
    }; // End class Cosmology

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_COSMOLOGY_H
