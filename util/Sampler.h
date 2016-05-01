#ifndef GCP_UTIL_SAMPLER_H
#define GCP_UTIL_SAMPLER_H

/**
 * @file Sampler.h
 * 
 * Tagged: Fri Mar 27 15:41:06 PDT 2009
 * 
 * 
 * @author Erik Leitch
 */

#define SAMPLER_FN(fn) double (fn)(double x, double* args)

#include <vector>

#include "gcp/util/Matrix.h"
#include "gcp/util/Vector.h"

namespace gcp {
  namespace util {

    class Sampler {
    public:

      /**
       * Constructor.
       */
      Sampler();

      /**
       * Destructor.
       */
      virtual ~Sampler();

      //------------------------------------------------------------
      // Methods for specifying the sampling function
      //------------------------------------------------------------

      // Methods for specifying the sampling function, as two arrays
      // specifying dy/dx
      
      void setdYdX(unsigned n, double* x, double* y);
      void setdYdX(std::vector<double>& x, std::vector<double>& y);

      // Method for specifying the sampling dy/dx as a function

      void setdYdX(SAMPLER_FN(fn), double xMin, double xMax, double dx, double* args);

      // Method for specifying the integral of dy/dx (Y(x' > x))
      // directly

      void setYX(SAMPLER_FN(fn));

      //------------------------------------------------------------
      // Generate samples according to the user-specified sampling function
      //------------------------------------------------------------

      double generateSample();

      std::vector<double> generateSamples(unsigned nSamp);

      //------------------------------------------------------------
      // Generate Poisson samples
      //------------------------------------------------------------

      static std::vector<unsigned> 
	generatePoissonSamples(double mean, unsigned nSamp);

      static unsigned
	generatePoissonSample(double mean);

      //------------------------------------------------------------
      // Generate univariate Gaussian samples
      //------------------------------------------------------------

      static double
	generateGaussianSample(double sigma);

      double
	generateGaussianSampleNonStatic(double sigma);

      //------------------------------------------------------------
      // Generate truncated univariate Gaussian samples
      //------------------------------------------------------------

      static double
	generateTruncatedGaussianSample(double sigma, double mean, double xmin, double xmax);

      //------------------------------------------------------------
      // Generate multi-variate Gaussian samples
      //------------------------------------------------------------

      static std::vector<double> 
	generateGaussianSamples(double sigma, unsigned nSamp);

      static std::vector<double> 
	generateTruncatedGaussianSamples(double mean, double sigma, double xmin, double xmax, unsigned nSamp);

      static Vector<double> 
	generateMultiVariateGaussianSample(Vector<double>& initialVal, 
					   Vector<double>& mean, 
					   Vector<double>& sigma,
					   Matrix<double>& correl);

      static Vector<double> 
	generateMultiVariateGaussianSample(Vector<double>& initialVal, 
					   Vector<double>& mean, 
					   Matrix<double>& invCov);

      static void 
	multiVariateSampleIterator(unsigned k,
				   Vector<double>& val, 
				   Vector<double>& mean, 
				   Matrix<double>& invCov);

      //------------------------------------------------------------
      // Generate uniform samples xmin <= x <= xmax
      //------------------------------------------------------------

      static double
	generateUniformSample(double xmin, double xmax);

      static std::vector<double>
	generateUniformSamples(double xmin, double xmax, unsigned nSamp);

      //------------------------------------------------------------
      // Seed the random number generator
      //------------------------------------------------------------

      static void seed(unsigned int s);
      static void seedRandom();

      //------------------------------------------------------------
      // Utility functions
      //------------------------------------------------------------

      // Return the natural log of the gamma function (series
      // approximation)

      static double lnGamma(double x);

      // Return the natural log of the factorial of n.  Exact for
      // small n, usese series approximation to the Gamma functions
      // for large n

      static double lnFactrl(unsigned n);

      //============================================================
      // Pdfs, Cdfs and Ptes for common variates
      //============================================================

      //------------------------------------------------------------
      // User-specified
      //------------------------------------------------------------

      double lnUserPdf(double x);
      double userPdf(double x);
      double userCdf(double x);
      double userPte(double x);

      //------------------------------------------------------------
      // Chi-square
      //------------------------------------------------------------

      static double lnChisqPdf(double chisq, unsigned nDof);
      static double chisqPdf(double chisq, unsigned nDof);
      static double chisqCdf(double chisq, unsigned nDof);
      static double chisqPte(double chisq, unsigned nDof);

      //------------------------------------------------------------
      // Gaussian
      //------------------------------------------------------------

      static double lnGaussPdf(double x, double mean, double sigma);
      static double gaussPdf(double x, double mean, double sigma);
      static double gaussCdf(double x, double mean, double sigma);
      static double gaussPte(double x, double mean, double sigma);

      //------------------------------------------------------------
      // Truncated gaussian
      //------------------------------------------------------------

      static double lnTruncatedGaussPdf(double x, double mean, double sigma, double xmin, double xmax, double& norm, bool& normIsCalculated);
      static double truncatedGaussPdf(double x, double mean, double sigma, double xmin, double xmax, double& norm, bool& normIsCalculated);
      static double truncatedGaussCdf(double x, double mean, double sigma, double xmin, double xmax, double& norm, bool& normIsCalculated);
      static double truncatedGaussPte(double x, double mean, double sigma, double xmin, double xmax, double& norm, bool& normIsCalculated);
      static double truncatedGaussIntegral(double mean, double sigma, double xmin, double xmax, double& norm, bool& normIsCalculated);

      //------------------------------------------------------------
      // Multivariate gaussian
      //------------------------------------------------------------

      static double lnMultiVariateGaussPdf(Vector<double>& x, 
					   Vector<double>& mean, 
					   Matrix<double>& invCov,
					   double detC);

      static double lnMultiVariateGaussPdf(Vector<double>& x, 
					   Vector<double>& mean, 
					   Vector<double>& sigma, 
					   Matrix<double>& correl);

      static double multiVariateGaussPdf(Vector<double>& x, 
					 Vector<double>& mean, 
					 Matrix<double>& invCov,
					 double detC);

      static double multiVariateGaussPdf(Vector<double>& x, 
					 Vector<double>& mean, 
					 Vector<double>& sigma,
					 Matrix<double>& correl);

      //------------------------------------------------------------
      // Poisson
      //------------------------------------------------------------

      static double lnPoissPdf(unsigned k, double lambda);
      static double poissPdf(unsigned k, double lambda);
      static double poissCdf(unsigned k, double lambda);
      static double poissPte(unsigned k, double lambda);

      //------------------------------------------------------------
      // Uniform
      //------------------------------------------------------------

      static double lnUniformPdf(double x, double xMin, double xMax);
      static double uniformPdf(double x, double xMin, double xMax);
      static double uniformCdf(double x, double xMin, double xMax);
      static double uniformPte(double x, double xMin, double xMax);

      //------------------------------------------------------------
      // Utility sampler functions
      //------------------------------------------------------------

      static SAMPLER_FN(gaussian);

      //------------------------------------------------------------
      // Utility functions
      //------------------------------------------------------------

      static int gammp(double a, double x, double *val);
      static double erfcc(double x);

      static void histogram(std::vector<double>& vals, unsigned nbin, std::vector<float>& x, std::vector<float>& n);

    private:

      static const double posInf_;
      static const double negInf_;
      static const double sqrt2pi_;

      // Binary search for samples

      double binSearchForSample();

      // True if we have a function to integrate

      bool haveFn_;

      // True if the integrated function is specified numerically, or
      // as a functional form

      bool isFn_;

      // The number of points in our integrated function

      unsigned nPt_;

      // Differential pdf

      std::vector<double> y_;      
      std::vector<double> x_;

      // Integral version of the above

      std::vector<double> yInt_;      
      std::vector<double> xInt_;

    }; // End class Sampler

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_SAMPLER_H
