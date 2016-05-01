#ifndef GCP_UTIL_LENGTH_H
#define GCP_UTIL_LENGTH_H

/**
 * @file Length.h
 * 
 * Tagged: Wed Aug 25 02:58:04 PDT 2004
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Area.h"
#include "gcp/util/ConformableQuantity.h"
#include "gcp/util/Matrix.h"
#include "gcp/util/Time.h"
#include "gcp/util/Volume.h"

#include <iostream>
#include <cmath>

namespace gcp {
  namespace util {
    
    class Area;
    class Volume;

    class Length : public ConformableQuantity {
    public:
      
      class Centimeters {};
      class Meters {};
      class Kilometers {};
      class Microns {};
      class Feet {};
      class Inches {};
      class Parsec {};
      class KiloParsec {};
      class MegaParsec {};
      class GigaParsec {};

      // Scale factors used by this class

      static const double cmPerM_;
      static const double mmPerCm_;
      static const double cmPerKm_;
      static const double micronsPerCm_;
      static const double cmPerInch_;
      static const double inchesPerFoot_;
      static const double cmPerParsec_;
      static const double cmPerAu_;
      static const double pcPerKpc_;
      static const double pcPerMpc_;
      static const double mpcPerGpc_;
      static const double pcPerLightYear_;
      static const double kmPerMpc_;
      static const double milesPerKm_;

      /**
       * Constructor.
       */
      Length();
      Length(const Centimeters& units, double cm);
      Length(const Meters& units, double m);
      Length(const Kilometers& units, double km);
      Length(const Feet& units, double feet);
      Length(const Inches& units, double inches);
      Length(const Parsec& units, double pc);
      Length(const KiloParsec& units, double kpc);
      Length(const MegaParsec& units, double mpc);
      Length(const GigaParsec& units, double gpc);
      
      /**
       * Copy constructor
       */
      Length(const Length& length);

      /**
       * Destructor.
       */
      virtual ~Length();
      
      /**
       * Set the length of this object
       */
      void setCentimeters(double cm); 

      void setMeters(double m) 
      {
	setCentimeters(m * cmPerM_);
      }

      void setMillimeters(double mm) 
      {
	setCentimeters(mm / mmPerCm_);
      }

      void setKilometers(double km) 
      {
	setCentimeters(km * cmPerKm_);
      }

      void setFeet(double feet)
      {
	setCentimeters(feet * inchesPerFoot_ * cmPerInch_);
      }

      void setInches(double inches)
      {
	setCentimeters(inches * cmPerInch_);
      }

      void setParsec(double parsec)
      {
	setCentimeters(parsec * cmPerParsec_);
      }

      void setGpc(double gpc)
      {
	setGigaParsec(gpc);
      }

      void setGigaParsec(double gpc)
      {
	setMegaParsec(gpc * mpcPerGpc_);
      }

      void setKpc(double kpc)
      {
	setKiloParsec(kpc);
      }

      void setKiloParsec(double kpc)
      {
	setParsec(kpc * pcPerKpc_);
      }
	
      void setMpc(double mpc)
      {
	setMegaParsec(mpc);
      }

      void setMegaParsec(double mpc)
      {
	setParsec(mpc * pcPerMpc_);
      }

      void setAu(double au)
      {
	setCentimeters(au * cmPerAu_);
      }

      // Set the length of this object as a light travel time

      void setLightTravelTime(Time time);

      /**
       * Get the length of this object
       */
      inline double centimeters() const {
	return val_;
      }

      inline double cm() const {
	return val_;
      }

      inline double meters() const {
	return val_ / cmPerM_;
      }

      inline double kilometers() const {
	return val_ / cmPerKm_;
      }

      inline double microns() const {
	return val_ * micronsPerCm_;
      }

      inline double feet() const {
	return val_ / (inchesPerFoot_ * cmPerInch_);
      }

      inline double inches() const {
	return val_ / cmPerInch_;
      }

      inline double parsec() const {
	return val_ / cmPerParsec_;
      }

      inline double kiloParsec() const {
	return parsec() / pcPerKpc_;
      }

      inline double Kpc() const {
	return kiloParsec();
      }

      inline double megaParsec() const {
	return parsec() / pcPerMpc_;
      }

      inline double Mpc() const {
	return megaParsec();
      }

      inline double gigaParsec() const {
	return megaParsec() / mpcPerGpc_;
      }

      inline double Gpc() const {
	return gigaParsec();
      }

      inline double lightYear() const {
	return parsec() / pcPerLightYear_;
      }

      Time time();

      // Operators

      /** 
       * Assignment
       */
      void operator=(Length& length) {
	val_ = length.val_;
      }

      void operator=(const Length& length) {
	val_ = length.val_;
      }

      /** 
       * Add two Lengths
       */
      Length operator+(Length& length);
      Length operator+(const Length& length);
      
      /** 
       * Subtract two Lengths
       */
      Length operator-(Length& length);
      Length operator-(const Length& length);

      /** 
       * Multiply a length by a constant
       */
      Length operator*(double multFac);
      void operator*=(double multFac);

      /** 
       * Divide two Lengths
       */
      double operator/(const Length& length);
      double operator/(Length& length);

      /** 
       * Multiply two Lengths
       */
      Area operator*(const Length& length);
      Area operator*(Length& length);

      Volume operator*(const Area& area);
      Volume operator*(Area& area);

      Length operator/(double fac);
      void operator/=(double fac);

      void operator+=(const Length& length);
      void operator+=(Length& length);

      void operator-=(const Length& length);
      void operator-=(Length& length);

      bool operator==(const Length& length);
      bool operator==(Length& length);

      bool operator>(const Length& length);
      bool operator>(Length& length);

      bool operator<(const Length& length);
      bool operator<(Length& length);

      /**
       * Allows cout << Length
       */
      friend std::ostream& operator<<(std::ostream& os, Length& length);
      friend std::ostream& operator<<(std::ostream& os, const Length& length);

      void initialize();

    }; // End class Length

    // Define a right matrix multiplier operator

    Vector<Length> operator*(Matrix<double>& mat, Vector<Length>& vec);

    void operator/=(Vector<Length>& vec, double fac);

    Vector<Length> operator*(Vector<double>& vec, Length& fac);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_LENGTH_H
