#ifndef GCP_UTIL_ASTROMETRY_H
#define GCP_UTIL_ASTROMETRY_H

/**
 * @file Astrometry.h
 * 
 * Tagged: Mon Aug  9 17:32:59 UTC 2004
 * 
 * @author 
 */
#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"

namespace gcp {
  namespace util {

    class QuadraticInterpolator;

    class Astrometry {
    public:
      
      struct Date {
	int year_;
	int month_;
	int day_;
	int hour_;
	int min_;
	int sec_;
	int nsec_;
      };

      /**
       * Constructor.
       */
      Astrometry();
      
      /**
       * Destructor.
       */
      virtual ~Astrometry();
      
      /**
       * Extend the quadratic interpolation table of ut1 - utc
       * versus MJD UTC.
       *
       * @throws Exception
       */
      void extendUt1Utc(double mjd, double ut1utc);
      
      /**
       * Extend the quadratic interpolation table of the equation
       * of the equinoxes versus Terrestrial Time (as a Modified
       * Julian Date).
       *
       * @throws Exception
       */
      void extendEqnEqx(double tt, double eqneqx);
      
      /**
       * Get the value of UT1-UTC for a given UTC.
       *
       * @throws Exception
       */
      double getUt1Utc(double mjd);
      
      /**
       * Get the value of the equation of the equinoxes for a
       * given terrestrial time.
       *
       * @throws Exception
       */
      double getEqnEqx(double tt);
      
      /**
       * Return true if ephemeris parameters can be interpolated for
       * this timestamp
       */
      bool canBracket(double mjdUtc);

      /**
       * Return the local sidereal time for a given site and UTC.
       *
       * Input:
       *
       *  utc      double    The current date and time (UTC), expressed as a
       *                     Modified Julian Date.
       *  longitude double   The longitude at which the lst is desired
       *  ut1Utc    double   The current value of UT1-UTC. If you don't need
       *                     more than one second of accuracy, this can be
       *                     given as zero.
       *
       *  eqnEqx   double    The current value of the equation of the
       *                     equinoxes if you want apparent sidereal time,
       *                     or 0 if you can make do with mean sidereal time.
       *                     The equation of the equinoxes is a slowly varying
       *                     number that is computationally intensive to 
       *                     calculate so it doesn't make sense to calculate 
       *                     it anew on each call.
       * Output:
       *
       *  return   double    The local sidereal time, expressed in radians.
       */
      static HourAngle mjdUtcToLst(double mjdUtc, Angle longitude, 
				   double ut1Utc, double eqnEqx);
      
      /**
       * Same as above, using internal ephemerides
       */
      HourAngle mjdUtcToLst(double mjdUtc, Angle longitude);
      
      /**
       * Return the Terestrial time (aka Ephemeris Time), corresponding to a
       * given UTC (expressed as a Modified Julian date).
       *
       * Input:
       *
       *  mjd      double   The Modified Julian date.
       *
       * Output:
       *
       *  return   double   The corresponding Terestrial Time
       */
      static double mjdUtcToMjdTt(double mjdUtc);
      static Date mjdUtcToCalendarDate(double mjdUtc);

      static void meanToApparentPlace(HourAngle& meanRa, Declination& meanDec, TimeVal& date, double equinox,
						  HourAngle& apparentRa, Declination& apparentDec);

      static void j2000ToApparentPlace(HourAngle& meanRa, Declination& meanDec, TimeVal& date, 
						   HourAngle& apparentRa, Declination& apparentDec);

      static void apparentToMeanPlace(HourAngle& apparentRa, Declination& apparentDec, TimeVal& date, 
						  double equinox,
						  HourAngle& meanRa, Declination& meanDec);

      static void apparentToJ2000Place(HourAngle& apparentRa, Declination& apparentDec, TimeVal& date, 
						   HourAngle& meanRa, Declination& meanDec);

      static void b1950ToJ2000(HourAngle& raB1950, Declination& decB1950, HourAngle& raJ2000, Declination& decJ2000);

      static void j2000ToB1950(HourAngle& raJ2000, Declination& decJ2000, HourAngle& raB1950, Declination& decB1950);

      static Angle angularSeparation(HourAngle ra1, Declination dec1, HourAngle ra2, Declination dec2);
      
      static void flatSkyApproximationSeparations(Angle& xSep, Angle& ySep,
						  HourAngle& ra1, Declination& dec1,
						  HourAngle& ra2, Declination& dec2);

    private:
      
      static const double secondsPerDay_;
      static const double pi_;
      static const double twopi_;
      
      /**
       * An interpolators for the UT1-UTC correction
       */
      QuadraticInterpolator* ut1Utc_;
      
      /**
       * An interpolator for the Equation of the equinox
       */
      QuadraticInterpolator* eqnEqx_;
      
    }; // End class Astrometry
    
  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_ASTROMETRY_H
