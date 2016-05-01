// $Id: $

#ifndef GCP_UTIL_STATS_H
#define GCP_UTIL_STATS_H

/**
 * @file Stats.h
 * 
 * Tagged: Tue Apr 23 09:57:58 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include <vector>

namespace gcp {
  namespace util {

    class Stats {
    public:

      /**
       * Constructor.
       */
      Stats();

      /**
       * Destructor.
       */
      virtual ~Stats();


      static void checkIndices(unsigned n, int& iStart, int& iStop);

      static double min(std::vector<double>& vec, int iStart=0, int iStop=0);
      static double max(std::vector<double>& vec, int iStart=0, int iStop=0);
      static double max(std::vector<double>& vec, unsigned& ind, int iStart=0, int iStop=0);
      static void minmax(std::vector<double>& vec, double& min, unsigned& iMin, double& max, unsigned& iMax, int iStart=0, int iStop=0);
      static double mean(std::vector<double>& vec, std::vector<unsigned>* multiplicity=0, int iStart=0, int iStop=0);
      static double mean(std::vector<unsigned>& vec, std::vector<unsigned>* multiplicity=0, int iStart=0, int iStop=0);
      static double sum(std::vector<double>& vec, std::vector<unsigned>* multiplicity=0, int iStart=0, int iStop=0);
      static double sum(std::vector<unsigned>& vec, std::vector<unsigned>* multiplicity=0, int iStart=0, int iStop=0);
      static double mode(std::vector<double>& vec, unsigned nbin, std::vector<unsigned>* multiplicity=0, int iStart=0, int iStop=0);
      static double variance(std::vector<double>& vec, std::vector<unsigned>* multiplicity=0, int iStart=0, int iStop=0);
      static double variance(std::vector<unsigned>& vec, std::vector<unsigned>* multiplicity=0, int iStart=0, int iStop=0);
      static double rms(std::vector<double>& vec, std::vector<unsigned>* multiplicity=0, int iStart=0, int iStop=0);
      static double rms(std::vector<unsigned>& vec, std::vector<unsigned>* multiplicity=0, int iStart=0, int iStop=0);
      static void histogram(std::vector<double>& vec, unsigned nbin, std::vector<double>& x, std::vector<double>& y, std::vector<unsigned>* multiplicity=0, int iStart=0, int iStop=0);
      static void confidenceIntervalN(std::vector<double>& vec, unsigned nbin, double nSigma, double refVal, double& lowVal, double& highVal, std::vector<unsigned>& multiplicity, int iStart=0, int iStop=0);
      static void confidenceIntervalN(std::vector<double>& vec, unsigned nbin, double nSigma, double refVal, double& lowVal, double& highVal, std::vector<unsigned>* multiplicity=0, int iStart=0, int iStop=0);

    private:
    }; // End class Stats

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_STATS_H
