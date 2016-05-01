// $Id: Distribution.h,v 1.1 2012/05/02 23:44:51 eml Exp $

#ifndef GCP_UTIL_DISTRIBUTION_H
#define GCP_UTIL_DISTRIBUTION_H

/**
 * @file Distribution.h
 * 
 * Tagged: Wed Apr 11 17:17:50 PDT 2012
 * 
 * @version: $Revision: 1.1 $, $Date: 2012/05/02 23:44:51 $
 * 
 * @author Erik Leitch
 */

#define EVAL_FN(fn) Probability (fn)(double x, Distribution* dist)
#define SAMP_FN(fn) double (fn)(Distribution* dist)

#include <vector>

#include "gcp/util/Probability.h"
#include "gcp/util/Sampler.h"

namespace gcp {
  namespace util {

    class Variate;

    class Distribution {
    public:

      // Enumerate known types of distributions

      enum Type {
	DIST_UNKNOWN     = 0x0,
	DIST_USER        = 0x1,
	DIST_GAUSS       = 0x2,
	DIST_CHISQ       = 0x4,
	DIST_POISS       = 0x8,
	DIST_UNIFORM     = 0x10,
	DIST_JOINT_GAUSS = 0x20,
	DIST_TRUNC_GAUSS = DIST_GAUSS | DIST_UNIFORM,
	DIST_UNSPEC      = 0x40
      };

      // Enumerate parameters that must be specified

      enum {
	PARAM_NONE         = 0x0,
	PARAM_UNKNOWN      = 0x0,

	PARAM_GAUSS_MEAN   = 0x1,
	PARAM_GAUSS_SIGMA  = 0x2,

	PARAM_CHISQ_NDOF   = 0x4,

	PARAM_UNIFORM_XMIN = 0x8,
	PARAM_UNIFORM_XMAX = 0x10,

	PARAM_USER_DYDX    = 0x20,

	PARAM_NORM         = 0x40,

	PARAM_POISS_MEAN   = 0x80,
      };

      /**
       * Constructor.
       */
      Distribution();

      /**
       * Destructor.
       */
      virtual ~Distribution();

      //------------------------------------------------------------
      // Set the type of this distribution
      //------------------------------------------------------------

      void setType(Type type);
      Type getType();

      //------------------------------------------------------------
      // Parameters for Poisson distributions
      //------------------------------------------------------------

      void setPoissMean(double mean);
      double getPoissMean();

      //------------------------------------------------------------
      // Parameters for gaussian distributions
      //------------------------------------------------------------

      void setGaussMean(double mean);
      void setGaussSigma(double sigma);

      double getGaussMean();
      double getGaussSigma();

      void setGaussMean(Variate& mean);
      void setGaussSigma(Variate& sigma);

      void invalidateTruncatedGaussianParameters();

      //------------------------------------------------------------
      // Parameters for uniform distributions
      //------------------------------------------------------------

      void setUniformXMin(double xMin);
      void setUniformXMax(double xMax);

      void setUniformXMin(Variate& xMin);
      void setUniformXMax(Variate& xMax);

      double getUniformXMin();
      double getUniformXMax();

      //------------------------------------------------------------
      // Parameters for chisq distributions
      //------------------------------------------------------------

      void setChisqNdof(unsigned nDof);
      unsigned chisqNdof();
      
      //------------------------------------------------------------
      // Parameters for user-specified distributions
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

      //-----------------------------------------------------------------------
      // Generate a sample from this distribution
      //-----------------------------------------------------------------------

      double sample();

      //-----------------------------------------------------------------------
      // Return the pdf of the passed value, for this distribution
      //-----------------------------------------------------------------------

      Probability pdf(double x);
      Probability cdf(double x);
      Probability pte(double x);
      double generateSample();

      friend std::ostream& operator<<(std::ostream& os, const Distribution& distribution);
      friend std::ostream& operator<<(std::ostream& os, Distribution& distribution);

    private:

      // Check the type validity

      void checkType(Type type);
      void checkTypes(Type type1, Type type2);
      void checkParams();

      static EVAL_FN(unspecPdf);
      static EVAL_FN(userPdf);
      static EVAL_FN(chisqPdf);
      static EVAL_FN(gaussPdf);
      static EVAL_FN(truncGaussPdf);
      static EVAL_FN(poissPdf);
      static EVAL_FN(uniformPdf);
      static EVAL_FN(jointGaussPdf);

      static EVAL_FN(unspecCdf);
      static EVAL_FN(userCdf);
      static EVAL_FN(chisqCdf);
      static EVAL_FN(gaussCdf);
      static EVAL_FN(truncGaussCdf);
      static EVAL_FN(poissCdf);
      static EVAL_FN(uniformCdf);
      static EVAL_FN(jointGaussCdf);

      static SAMP_FN(unspecSample);
      static SAMP_FN(userSample);
      static SAMP_FN(chisqSample);
      static SAMP_FN(gaussSample);
      static SAMP_FN(truncGaussSample);
      static SAMP_FN(poissSample);
      static SAMP_FN(uniformSample);
      static SAMP_FN(jointGaussSample);

      // The functions to evaluate when pdf() or cdf() are called

      EVAL_FN(*pdfFn_);
      EVAL_FN(*cdfFn_);
      SAMP_FN(*sampleFn_);

      // The type of this distribution

      Type type_;

      // Required parameters for this type

      unsigned requiredMask_;

      // Parameters that have been specified

      unsigned paramMask_;

      // Parameters for gaussian distributions

      struct {
	double mean_;
	double sigma_;
	double norm_;            // Used for truncated gaussians only
	bool normIsCalculated_;  // Used for truncated gaussians only
      } gauss_;

      // Parameters for poisson distributions

      struct {
	double mean_;
      } poiss_;

    public:

      // Parameters for chisq distributions

      struct {
	unsigned nDof_;
      } chisq_;

      // Parameters for uniform distributions

      struct {
	double xMin_;
	double xMax_;
      } uniform_;

      // Parameters for arbitrary user-specified distributions

      struct {
	Sampler sampler_;
      } user_;

    }; // End class Distribution

    std::ostream& operator<<(std::ostream& os, const Distribution& distribution);
    std::ostream& operator<<(std::ostream& os, Distribution& distribution);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_DISTRIBUTION_H
