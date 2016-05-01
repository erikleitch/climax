// $Id: UniformVariate.h,v 1.1 2012/05/02 23:44:52 eml Exp $

#ifndef GCP_UTIL_UNIFORMVARIATE_H
#define GCP_UTIL_UNIFORMVARIATE_H

/**
 * @file UniformVariate.h
 * 
 * Tagged: Mon Apr 16 10:45:12 PDT 2012
 * 
 * @version: $Revision: 1.1 $, $Date: 2012/05/02 23:44:52 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Variate.h"

namespace gcp {
  namespace util {

    class UniformVariate : public Variate {
    public:

      /**
       * Constructor.
       */
      UniformVariate();
      UniformVariate(Variate& xMin, Variate& xMax);
      UniformVariate(double xMin, double xMax);

      /**
       * Destructor.
       */
      virtual ~UniformVariate();

    private:

      void initialize();

    }; // End class UniformVariate

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_UNIFORMVARIATE_H
