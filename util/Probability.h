// $Id: Probability.h,v 1.1 2012/05/02 23:44:51 eml Exp $

#ifndef GCP_UTIL_PROBABILITY_H
#define GCP_UTIL_PROBABILITY_H

/**
 * @file Probability.h
 * 
 * Tagged: Tue Apr 17 22:51:50 PDT 2012
 * 
 * @version: $Revision: 1.1 $, $Date: 2012/05/02 23:44:51 $
 * 
 * @author Erik Leitch
 */
#include <sstream>

namespace gcp {
  namespace util {

    class Probability {
    public:

      /**
       * Constructor.
       */
      Probability();

      /**
       * Destructor.
       */
      virtual ~Probability();

      void setValue(double value);
      void setLnValue(double lnValue);

      double value();
      double lnValue();

      void operator*=(const Probability& prob);
      void operator*=(Probability& prob);
      Probability operator*(const Probability& prob);
      Probability operator*(Probability& prob);

      void operator/=(const Probability& prob);
      void operator/=(Probability& prob);
      Probability operator/(const Probability& prob);
      Probability operator/(Probability& prob);

      bool operator<(const Probability& prob);
      bool operator>(const Probability& prob);

      bool operator<(Probability& prob);
      bool operator>(Probability& prob);
      bool operator<(double val);
      bool operator>(double val);

      bool isValid();
      bool isValid(double val);

      friend std::ostream& operator<<(std::ostream& os, const Probability& probability);
      friend std::ostream& operator<<(std::ostream& os, Probability& probability);

      bool lessThan(double val1, double val2);
      bool greaterThan(double val1, double val2);

    private:

      // Privately, probabilities are stored as ln(prob)

      double lnProb_;

    }; // End class Probability

    std::ostream& operator<<(std::ostream& os, const Probability& probability);
    std::ostream& operator<<(std::ostream& os, Probability& probability);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PROBABILITY_H
