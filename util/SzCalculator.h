// $Id: SzCalculator.h,v 1.1.1.1 2010/07/13 17:56:53 eml Exp $

#ifndef GCP_UTIL_SZCALCULATOR_H
#define GCP_UTIL_SZCALCULATOR_H

/**
 * @file SzCalculator.h
 * 
 * Tagged: Fri Aug  1 22:12:12 PDT 2008
 * 
 * @version: $Revision: 1.1.1.1 $, $Date: 2010/07/13 17:56:53 $
 * 
 * @author tcsh: Erik Leitch.
 */
#include "gcp/util/Flux.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/Intensity.h"
#include "gcp/util/Temperature.h"

namespace gcp {
  namespace util {

    class SzCalculator {
    public:

      /**
       * Constructor.
       */
      SzCalculator();

      /**
       * Destructor.
       */
      virtual ~SzCalculator();

      // Calculate the factor by which comptonY should be multiplied
      // to convert to CMB temperature decrement/increment

      static Temperature comptonYToDeltaT(Frequency& freq);

      // Calculate the factor by which comptonY should be multiplied
      // to convert to Flux/sr

      static Intensity comptonYToDeltaI(Frequency& freq);

      // Calculate the dimensionless x factor that enters into the
      // Planck function

      static double planckX(Frequency& freq, Temperature& temp);

      // Calculate the frequency given a dimensionless x

      static Frequency planckFreq(double x, Temperature& temp);

      // Evaluate the Planck function and its derivative wrt to T at
      // given T and freq.

      static Intensity dPlanck(Frequency& freq, Temperature& temp);
      static Intensity planck(Frequency& freq, Temperature& temp);

      //-----------------------------------------------------------------------
      // Versions of the above that involve no intermediate object
      // creation, for efficient computation
      //-----------------------------------------------------------------------

      static void comptonYToDeltaT(Frequency& freq, Temperature& temp);
      static void comptonYToDeltaI(Frequency& freq, Temperature& temp, Intensity& intensity);
      static void dPlanck(Frequency& freq, Temperature& temp, Intensity& intensity);
      static void planck(Frequency& freq, Temperature& temp, Intensity& intensity);

      static void comptonYToDeltaTItoh(Temperature& electronTemperature, Frequency& freq, Temperature& YtoT);
      static void comptonYToDeltaIItoh(Temperature& Te, Frequency& freq, Temperature& temp, Intensity& intensity);
      static double itohFn(Temperature& electronTemperature, Frequency& freq);
      static double kompFn(Temperature& electronTemperature, Frequency& freq);

      static void itohX(double x, double& XI, double& SI);
      static std::vector<std::vector<double> > initializeItohY0Coeffs();
      static std::vector<std::vector<double> > initializeItohY1Coeffs();
      static std::vector<std::vector<double> > initializeItohY2Coeffs();
      static std::vector<std::vector<double> > initializeItohY3Coeffs();
      static std::vector<std::vector<double> > initializeItohY4Coeffs();

      static double computeItohY(double XI, double SI, 
				 std::vector<std::vector<double> >& coeffArr);
      static double itohY4(double XI, double SI);
      static double itohY3(double XI, double SI);
      static double itohY1(double XI, double SI);
      static double itohY2(double XI, double SI);
      static double itohY0(double XI, double SI);

      static double itohY1Comp(double XI, double SI);

    private:

      static std::vector<std::vector<double> > itohY0Coeffs_;
      static std::vector<std::vector<double> > itohY1Coeffs_;
      static std::vector<std::vector<double> > itohY2Coeffs_;
      static std::vector<std::vector<double> > itohY3Coeffs_;
      static std::vector<std::vector<double> > itohY4Coeffs_;

    }; // End class SzCalculator

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_SZCALCULATOR_H
