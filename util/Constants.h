#ifndef GCP_UTIL_CONSTANTS_H
#define GCP_UTIL_CONSTANTS_H

/**
 * @file Constants.h
 * 
 * Tagged: Tue Aug 10 13:17:59 PDT 2004
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Length.h"
#include "gcp/util/Mass.h"
#include "gcp/util/Speed.h"
#include "gcp/util/Temperature.h"

namespace gcp {
  namespace util {
    
    class Constants {
    public:
      
      static Temperature Tcmb_;
      static const double hPlanckCgs_;
      static const double kBoltzCgs_;
      static const double kBoltzSi_;
      static const double JyPerCgs_;
      static const double electronMassCgs_;
      static const double protonMassCgs_;
      static const double utSecPerSiderealSec_;
      static const double gravitationalConstantCgs_;
      static const double gravitationalConstantSi_;

      static Area   sigmaT_;
      static Speed  lightSpeed_;
      static Length au_;
      static Length defaultEarthRadius_;
      static Mass   electronMass_;
      static Mass   protonMass_;
      static Mass   solarMass_;

      static const double milesPerDegree_;

      /**
       * Constructor.
       */
      Constants();
      
      /**
       * Destructor.
       */
      virtual ~Constants();
      
      virtual double cgs() {
	return cgs_;
      }

      virtual double si() {
	return si_;
      }

    protected:

      void setGgs(double cgs) {
	cgs = cgs_;
      }

      void setSi(double si) {
	si = si_;
      }
      
    private:

      double cgs_;
      double si_;

    }; // End class Constants
    
  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_CONSTANTS_H
