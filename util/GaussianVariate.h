// $Id: GaussianVariate.h,v 1.1 2012/05/02 23:44:51 eml Exp $

#ifndef GCP_UTIL_GAUSSIANVARIATE_H
#define GCP_UTIL_GAUSSIANVARIATE_H

/**
 * @file GaussianVariate.h
 * 
 * Tagged: Fri Nov 19 11:41:00 PST 2010
 * 
 * @version: $Revision: 1.1 $, $Date: 2012/05/02 23:44:51 $
 * 
 * @author tcsh: Erik Leitch
 */
#include <iostream>
#include "gcp/util/Variate.h"

namespace gcp {
  namespace util {

    class GaussianVariate : public Variate {
    public:

      /**
       * Constructor.
       */
      GaussianVariate();

      /**
       * Destructor.
       */
      virtual ~GaussianVariate();

      void initialize();

      void setMean(double mean);
      void setSigma(double sigma);
      void setVal(double val);

      double val();

      // Allows cout << Variate

      friend std::ostream& operator<<(std::ostream& os, const GaussianVariate& gvar);
      friend std::ostream& operator<<(std::ostream& os, GaussianVariate& gvar);

    }; // End class GaussianVariate

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_GAUSSIANVARIATE_H
