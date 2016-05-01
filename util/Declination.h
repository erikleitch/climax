#ifndef GCP_UTIL_DECLINATION_H
#define GCP_UTIL_DECLINATION_H

/**
 * @file Declination.h
 * 
 * Tagged: Tue Aug 10 14:03:21 PDT 2004
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Angle.h"

namespace gcp {
  namespace util {
    
    class Declination : public Angle {
    public:
      
      /**
       * Constructor.
       */
      Declination() {
	initialize();
      };
      
      /**
       * Destructor.
       */
      virtual ~Declination() {};

      inline double arcMinutes() {
	return val_ * arcMinPerRad_;
      }

      inline double arcSeconds() {
	return val_ * arcSecPerRad_;
      }

      /** 
       * Add two Declinations
       */
      Declination operator+(Declination& angle) {
	Declination sum;
	sum.setRadians(val_);
	sum.addRadians(angle.radians());
	return sum;
      }

      /** 
       * Add an angle to this Declination
       */
      Declination operator+(Angle& angle) {
	Declination sum;
	sum.setRadians(val_);
	sum.addRadians(angle.radians());
	return sum;
      }
      
      /** 
       * Subtract two Declinations
       */
      Declination operator-(Declination& angle) {
	Declination diff;
	diff.setRadians(val_);
	diff.addRadians(-angle.radians());
	return diff;
      }

      /** 
       * Subtract an angle from this Declination
       */
      Declination operator-(Angle& angle) {
	Declination diff;
	diff.setRadians(val_);
	diff.addRadians(-angle.radians());
	return diff;
      }

      inline int getIntegerDegrees() {
	bool neg = val_ < 0.0;
	double arad = fabs(val_);
	unsigned degs = (unsigned)(arad * degPerRad_);

	return (neg ? -1 : 1)*degs;
      }

      inline unsigned getIntegerArcMinutes() {
	double arad = fabs(val_);
	unsigned degs  = (unsigned)(arad * degPerRad_);
	unsigned mins  = (unsigned)((arad * degPerRad_ - degs)*60);
	return mins;
      }

      inline unsigned getIntegerArcSeconds() {
	double arad = fabs(val_);
	unsigned degs  = (unsigned)(arad   * degPerRad_);
	unsigned mins  = (unsigned)((arad  * degPerRad_ - degs)*60);
	unsigned secs  = (unsigned)(((arad * degPerRad_ - degs)*60 - mins)*60);
	return secs;
      }

      inline unsigned getIntegerMilliArcSeconds() {
	double arad = fabs(val_);
	unsigned degs  = (unsigned)(arad   * degPerRad_);
	unsigned mins  = (unsigned)((arad  * degPerRad_ - degs)*60);
	unsigned secs  = (unsigned)(((arad * degPerRad_ - degs)*60 - mins)*60);
	return (unsigned)((((arad * degPerRad_ - degs)*60 - mins)*60 - secs) * 1000);
      }

    private:

      static const double arcSecPerRad_;
      static const double arcMinPerRad_;

    }; // End class Declination
    
  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_DECLINATION_H
