// $Id: SolidAngle.h,v 1.2 2012/05/02 23:44:51 eml Exp $

#ifndef GCP_UTIL_SOLIDANGLE_H
#define GCP_UTIL_SOLIDANGLE_H

/**
 * @file SolidAngle.h
 * 
 * Tagged: Wed Sep 14 17:52:22 PDT 2005
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/02 23:44:51 $
 * 
 * @author Erik Leitch
 */
#include <iostream>

#include "gcp/util/ConformableQuantity.h"
#include "gcp/util/Angle.h"

namespace gcp {
  namespace util {

    class SolidAngle : public ConformableQuantity {
    public:

      class Steradians {};
      class SqArcMinutes {};
      class SqDegrees {};

      /**
       * Constructor.
       */
      SolidAngle();
      SolidAngle(const Steradians& units,   double sr);
      SolidAngle(const SqArcMinutes& units, double sqarcmin);
      SolidAngle(const SqDegrees& units, double sqdeg);
      SolidAngle(Angle& fwhm);
      SolidAngle(Angle& fwhma, Angle& fwhmb);

      /**
       * Destructor.
       */
      virtual ~SolidAngle();

      void initialize();

      void setSr(double sr);
      void setSqDegrees(double sqdeg);
      void setSqArcMin(double sqarcmin);

      inline double sqArcMin() {
	return val_ * Angle::arcMinPerRad_ * Angle::arcMinPerRad_;
      }

      inline double sqDegrees() {
	return val_ * Angle::degPerRad_ * Angle::degPerRad_;
      }

      inline double sr() {
	return val_;
      }

      void operator=(const SolidAngle& angle);
      void operator=(SolidAngle& angle);

    }; // End class SolidAngle

    SolidAngle operator*(const Angle& a1, const Angle& a2);
    SolidAngle operator*(Angle& a1, Angle& a2);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_SOLIDANGLE_H
