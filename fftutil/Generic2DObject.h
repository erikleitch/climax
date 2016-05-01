// $Id: $

#ifndef GCP_UTIL_GENERIC2DOBJECT_H
#define GCP_UTIL_GENERIC2DOBJECT_H

/**
 * @file Generic2DObject.h
 * 
 * Tagged: Wed Mar 20 09:48:38 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Declination.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/SolidAngle.h"
#include "gcp/util/Unit.h"

#include "gcp/fftutil/ImageAxis.h"

namespace gcp {
  namespace util {

    class Generic2DObject {
    public:

      /**
       * Constructor.
       */
      Generic2DObject();

      /**
       * Destructor.
       */
      virtual ~Generic2DObject();

      virtual void operator*=(double mult) = 0;

      bool hasData_;

      bool hasAbsolutePosition_; // True if this object has an absolute
                                 // position specified

      HourAngle ra_;             // If specified, the RA associated with the reference position of this object
      Declination dec_;          // If specified, the DEC associated with the reference position of this object

      Unit::Units units_;        // The native units of this image
      bool hasUnits_;

      Frequency frequency_;
      bool hasFrequency_;

      Unit::Units getUnits();

      double nativeToJy(Frequency& nu, SolidAngle& beam);
      void convertToJy(Frequency& nu, SolidAngle& beam);
      void setUnits(Unit::Units units);
      void setUnits(std::string units);

      virtual void setRaDec(HourAngle& ra, Declination& dec);
      virtual void setFrequency(Frequency& freq);

      virtual HourAngle&   getRa();
      virtual Declination& getDec();

      virtual ImageAxis& xImageAxis() = 0;
      virtual ImageAxis& yImageAxis() = 0;

      // Return true if this object has been initialized with data

      bool hasData();
      void setHasData(bool hasData);

    }; // End class Generic2DObject

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_GENERIC2DOBJECT_H
