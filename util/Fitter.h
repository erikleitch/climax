#ifndef GCP_UTIL_FITTER_H
#define GCP_UTIL_FITTER_H

/**
 * @file Fitter.h
 * 
 * Tagged: Mon Nov 25 13:45:47 PST 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include <vector>

#define FITTER_EVAL_FN(fn) void (fn)(float x, float a[], float *y, float dyda[], int npar)
#define FITTER_EVAL_FN2D(fn) void (fn)(float x1, float x2, float a[], float *y, float dyda[], int npar)

namespace gcp {
  namespace util {

    class Fitter {
    public:

      enum  Func_type {F_GAUSS, F_GAUSSB, F_POLY, F_POWSPEC, F_POWER, F_GAUSS2D, F_GAUSSTEST};

      struct Fit_scan {
	std::vector<double> x;
	std::vector<double> x2; // If fitting 2D function
	std::vector<double> y;
	std::vector<double> sd;
	int npt;

	Fit_scan(int npt);
	~Fit_scan();
      };
    
      // Constructor.
    
      Fitter();
    
      // Destructor.

      virtual ~Fitter();
    
      static FITTER_EVAL_FN(poly);
      static FITTER_EVAL_FN(gauss_nobase);
      static FITTER_EVAL_FN(gausstest);
      static FITTER_EVAL_FN(powspec);
      static FITTER_EVAL_FN(power);
      static FITTER_EVAL_FN2D(gauss2d);

      void covsrt(float **covar, int ma, int ia[], int mfit);
      float** matrix(long nrow, long ncol);
      float** free_matrix(float **m);

      int mrqmin2d(std::vector<double>& x1, std::vector<double>& x2, std::vector<double>& y, 
		   std::vector<double>& sig, float a[], 
		   int ia[], int ma, float **covar, float **alpha, float *chisq,
		   FITTER_EVAL_FN2D(*funcs),
		   float *alamda);

      int mrqcof2d(std::vector<double>& x1, std::vector<double>& x2, std::vector<double>& y, 
		   std::vector<double>& sig, float a[], 
		   int ia[], int ma, float **alpha, float beta[], float *chisq,
		   FITTER_EVAL_FN2D(*funcs));

      int mrqmin(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, float a[], 
		 int ia[], int ma, float **covar, float **alpha, float *chisq,
		 FITTER_EVAL_FN(*funcs),
		 float *alamda);

      int mrqcof(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, float a[], 
		 int ia[], int ma, float **alpha, float beta[], float *chisq,
		 FITTER_EVAL_FN(*funcs));

      int fit(Fit_scan *fitscan, float apar[], int ma, Fitter::Func_type type, 
	      float *chisq_fin, int printchisq);

      Fit_scan getPowerSpectrumFitData(std::vector<double>& yarr, std::vector<unsigned>* multiplicity, 
				       double& p0, double& ks, double& alpha, unsigned& nSamp);

      Fit_scan getPowerSpectrumFitDataNew(std::vector<double>& yarr, std::vector<unsigned>* multiplicity, unsigned nFit,
					  double& p0, double& ks, double& alpha, unsigned& nSamp);

      void fitPowerSpectrum(std::vector<double>& yarr, std::vector<unsigned>* multiplicity, unsigned nFit,
                            double& p0, double& ks, double& alpha, unsigned& nSamp);

      void plotPowerSpectrum(std::vector<double>& yarr, std::vector<unsigned>* multiplicity, unsigned nFit,
                             double& p0, double& ks, double& alpha);

      void getIndices(unsigned n, unsigned iInterval, double frac, int& iStart, int& iStop);

      int gaussj(float **a, int n, float **b, int m);

    }; // End class Fitter

  } // End namespace util
} // End namespace gcp

#endif // End #ifndef GCP_UTIL_FITTER_H
