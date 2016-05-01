// $Id: Time.h,v 1.2 2012/05/02 23:44:52 eml Exp $

#ifndef GCP_UTIL_TIME_H
#define GCP_UTIL_TIME_H

/**
 * @file Time.h
 * 
 * Tagged: Mon Oct  3 16:09:14 PDT 2005
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/02 23:44:52 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/ConformableQuantity.h"

namespace gcp {
  namespace util {

    class Time : public ConformableQuantity {
    public:

      class Seconds {};
      class NanoSeconds {};

      // Scale factors used by this class

      static const double nsPerSecond_;
      static const double secPerMinute_;
      static const double secPerHour_;
      static const double secPerDay_;
      static const double dayPerYear_;

      /**
       * Constructor.
       */
      Time();
      Time(const Time& time);
      Time(const Seconds& units, double s);
      Time(const NanoSeconds& units, double ns);

      /**
       * Destructor.
       */
      virtual ~Time();

      // Set the time of this object

      void setSeconds(double s)
	{
	  val_ = s;
	}

      void setNanoSeconds(double ns)
	{
	  val_ = ns / nsPerSecond_;
	}

      inline double years() const {
	return val_ / (secPerDay_ * dayPerYear_);
      }

      inline double days() const {
	return val_ / secPerDay_;
      }

      inline double hours() const {
	return val_ / secPerHour_;
      }

      inline double minutes() const {
	return val_ / secPerMinute_;
      }

      inline double seconds() const {
	return val_;
      }

      inline double nanoSeconds() const {
	return val_ * nsPerSecond_;
      }

      void initialize();

      friend std::ostream& operator<<(std::ostream& os, Time& time);

    }; // End class Time

    std::ostream& operator<<(std::ostream& os, Time& time);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_TIME_H
