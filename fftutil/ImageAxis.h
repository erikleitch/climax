// $Id: ImageAxis.h,v 1.5 2012/05/29 16:20:28 eml Exp $

#ifndef GCP_UTIL_IMAGEAXIS_H
#define GCP_UTIL_IMAGEAXIS_H

/**
 * @file Axis.h
 * 
 * Tagged: Wed Jun 23 15:11:11 PDT 2010
 * 
 * @version: $Revision: 1.5 $, $Date: 2012/05/29 16:20:28 $
 * 
 * @author Erik Leitch.
 */
#include "gcp/util/Angle.h"

namespace gcp {
  namespace util {

    class ImageAxis {
    public:

      enum {
	AXIS_NONE = 0x0,
	AXIS_X    = 0x1,
	AXIS_Y    = 0x2,
	AXIS_BOTH = AXIS_X | AXIS_Y
      };

      /**
       * Constructor.
       */
      ImageAxis();

      /**
       * Destructor.
       */
      virtual ~ImageAxis();

      // Set methods
      
      virtual void setNpix(unsigned n);
      virtual void setAngularSize(Angle delta);
      virtual void setSpatialFrequencyResolution(double inverseRadians);
      virtual void setProjection(std::string projection);
      virtual void setSense(int sense);

      // Check methods
      
      bool hasNpix();
      bool hasAngularSize();
      
      // Get methods
      
      unsigned getNpix();
      Angle& getAngularSize();
      Angle& getAngularResolution();
      std::string getProjection();
      int getSense();

      virtual double getSpatialFrequencyResolution();
      double getMinimumSpatialFrequency();
      double getMaximumSpatialFrequency();
      
      void setAxisType(unsigned type);

      bool operator==(const ImageAxis& axis);
      bool operator==(ImageAxis& axis);

      ImageAxis(const ImageAxis& axis);
      ImageAxis(ImageAxis& axis);
      void operator=(const ImageAxis& axis);
      void operator=(ImageAxis& axis);

      virtual void setScale(double scale, std::string name);

      std::string name();
      double scale();

      double scale_;
      std::string name_;
      Angle size_;
      Angle resolution_;
      unsigned n_;
      
      bool sizeIsSpecified_;
      bool nIsSpecified_;

      // Sense of this axis.  
      //
      // If sense =  1, the axis values increase from 0 --> n
      // If sense = -1, the axis values decrease from 0 --> n

      int sense_;

      unsigned type_;

      std::string projection_;

    }; // End class ImageAxis

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_IMAGEAXIS_H
