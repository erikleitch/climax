// $Id: ChisqVariate.h,v 1.3 2012/05/08 21:58:00 eml Exp $

#ifndef GCP_UTIL_CHISQVARIATE_H
#define GCP_UTIL_CHISQVARIATE_H

/**
 * @file ChisqVariate.h
 * 
 * Tagged: Fri Nov 19 11:41:00 PST 2010
 * 
 * @version: $Revision: 1.3 $, $Date: 2012/05/08 21:58:00 $
 * 
 * @author tcsh: Erik Leitch
 */
#include <iostream>
#include "gcp/util/Variate.h"

namespace gcp {
  namespace util {

    class ChisqVariate : public Variate {
    public:

      /**
       * Constructor.
       */
      ChisqVariate();

      /**
       * Destructor.
       */
      virtual ~ChisqVariate();

      void initialize();
      void setChisq(double chisq, unsigned nDof);
      void setNdof(unsigned nDof);

      virtual Probability pdf();
      Probability likelihood();

      //------------------------------------------------------------
      // Access methods
      //------------------------------------------------------------

      double chisq();
      double reducedChisq();
      unsigned nDof();

      // Add a chisqVal to this object

      inline void operator+=(double chisqVal) {
	unsigned nDofNew = nDof() + 1;
	val_ += (chisqVal - val_) / nDofNew;
	setNdof(nDofNew);
      };
      
      inline void directAdd(double chisqVal, unsigned ndof) {
	unsigned nDofNew = nDof() + ndof;
	val_ = val_ + (chisqVal - ndof*val_) / nDofNew;
	setNdof(nDofNew);
      };

      // Add two variates

      void operator+=(const ChisqVariate& chisq);
      void operator+=(ChisqVariate& chisq);

      // Compare two variates

      bool operator<(const ChisqVariate& chisq);
      bool operator<(ChisqVariate& chisq);

      void operator=(const ChisqVariate& chisq);
      void operator=(ChisqVariate& chisq);

      // Allows cout << Chisq

      friend std::ostream& operator<<(std::ostream& os, const ChisqVariate& chisq);
      friend std::ostream& operator<<(std::ostream& os, ChisqVariate& chisq);

    }; // End class ChisqVariate

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_CHISQVARIATE_H
