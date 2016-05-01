// $Id: Pressure.h,v 1.2 2012/05/02 23:44:51 eml Exp $

#ifndef GCP_UTIL_PRESSURE_H
#define GCP_UTIL_PRESSURE_H

/**
 * @file Pressure.h
 * 
 * Tagged: Fri Nov 14 17:07:29 PST 2008
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/02 23:44:51 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/ConformableQuantity.h"
#include "gcp/util/NumberDensity.h"
#include "gcp/util/Temperature.h"

namespace gcp {
  namespace util {

    class NumberDensity;
    class Temperature;

    class Pressure : public ConformableQuantity {
    public:
      
      class MilliBar {};

      /**
       * Constructor.
       */
      Pressure();
      Pressure(const MilliBar& unit, double mBar);

      /**
       * Destructor.
       */
      virtual ~Pressure();

      static const double barPerPascal_;
      static const double milliBarPerBar_;
      static const double pascalPerTorr_;
      static const double dynPerCm2PerPascal_;
      static const double keVPerCm3PerDynPerCm2_;
      static const double keVPerCm3PerMilliBar_;
      static const double ergPerCm3PerMilliBar_;

      void initialize();

      void setMilliBar(double mBar);
      void setBar(double bar);
      void setPascal(double pascal);
      void setDynPerCm2(double dynPerCm2);
      void setKeVPerCm3(double keVPerCm3);
      void setErgPerCm3(double ergPerCm3);

      double milliBar();
      double bar();
      double pascals();
      double torr();
      double mmHg();
      double keVPerCm3();
      double ergPerCm3();
      double dynPerCm2();

      double operator/(const Pressure& press);
      double operator/(Pressure& press);

      Pressure operator*(double fac);
      void operator/=(double fac);

      //------------------------------------------------------------
      // A pressure divided by a temperature implies a
      // number densty (P/kT = n)
      //------------------------------------------------------------

      NumberDensity operator/(const Temperature& temp);
      NumberDensity operator/(Temperature& temp);

      void operator=(const Pressure& pressure);
      void operator=(Pressure& pressure);

      friend std::ostream& operator<<(std::ostream& os, Pressure& pressure);

    }; // End class Pressure

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PRESSURE_H
