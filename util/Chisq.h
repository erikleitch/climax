// $Id: Chisq.h,v 1.1 2011/01/20 19:40:44 eml Exp $

#ifndef GCP_UTIL_CHISQ_H
#define GCP_UTIL_CHISQ_H

/**
 * @file Chisq.h
 * 
 * Tagged: Thu Nov 11 23:35:43 PST 2010
 * 
 * @version: $Revision: 1.1 $, $Date: 2011/01/20 19:40:44 $
 * 
 * @author tcsh: Erik Leitch
 */
#include <iostream>

namespace gcp {
  namespace util {

    class Chisq {
    public:

      /**
       * Constructor.
       */
      Chisq();

      /**
       * Destructor.
       */
      virtual ~Chisq();

      double chisq();
      double reducedChisq();
      unsigned nChisq();

      // Initialize this object to zero

      void initialize();

      // Add a chisqVal to this object

      void operator+=(double chisqVal);

      // Add two chis-squares

      void operator+=(const Chisq& chisq);
      void operator+=(Chisq& chisq);

      // Allows cout << Chisq

      friend std::ostream& operator<<(std::ostream& os, const Chisq& chisq);
      friend std::ostream& operator<<(std::ostream& os, Chisq& chisq);

      double probabilityToExceed();
      void setChisq(double chisq, unsigned nChisq);

    private:

      double reducedChisq_;
      unsigned nChisq_;

    }; // End class Chisq

    std::ostream& operator<<(std::ostream& os, const Chisq& chisq);
    std::ostream& operator<<(std::ostream& os, Chisq& chisq);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_CHISQ_H
