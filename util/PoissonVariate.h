// $Id: $

#ifndef GCP_UTIL_POISSONVARIATE_H
#define GCP_UTIL_POISSONVARIATE_H

/**
 * @file PoissonVariate.h
 * 
 * Tagged: Tue Mar 25 14:10:14 PDT 2014
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include <iostream>
#include "gcp/util/Variate.h"

namespace gcp {
  namespace util {

    class PoissonVariate : public Variate {
    public:

      /**
       * Constructor.
       */
      PoissonVariate();

      /**
       * Destructor.
       */
      virtual ~PoissonVariate();

      void initialize();
      void setValue(unsigned k, double mean);
      void setMean(double mean);
      double mean();

#if 0
      virtual Probability pdf();
      Probability likelihood();
#endif
      void plotPdf(double min, double max, unsigned n);

    private:

    }; // End class PoissonVariate

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_POISSONVARIATE_H
