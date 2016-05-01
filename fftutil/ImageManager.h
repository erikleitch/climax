// $Id: $

#ifndef GCP_UTIL_IMAGEMANAGER_H
#define GCP_UTIL_IMAGEMANAGER_H

/**
 * @file ImageManager.h
 * 
 * Tagged: Mon Feb  9 09:55:24 PST 2015
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include "fftutil/Image.h"

#include "util/ParameterManager.h"

namespace gcp {
  namespace util {

    class ObsInfo;

    class ImageManager : public ParameterManager {
    public:

      /**
       * Constructor.
       */
      ImageManager();

      /**
       * Destructor.
       */
      virtual ~ImageManager();

      Image initializeImage(std::string fileName, ObsInfo* obs=0);
      Image initializeImage(std::string fileName, Image& image, ObsInfo* obs=0);
      Image initializeImageParameters(Image& image, bool rebin=true);
      bool isImage(std::string fileName);
      bool isEvent(std::string fileName);
      void getNpix(unsigned& nx, unsigned& ny);
      double getVal(Image& image, std::string str);

      void parseCoordinate(String& valStr, double& x, double& y, String& unit);
      void parseVal(String& valStr, double& val, String& unit);
      void parseRange(String& valStr, double& valMin, double& valMax, String& unit);
      void parseRanges(String& valStr, double& xMin, double& xMax, 
		       double& yMin, double& yMax, String& unit);

      Image extractRegion(String& region, Image& image);
      void convertRange(String& range, Image& image, 
			int& ixmin, int& ixmax, int& iymin, int& iymax,
			bool isRel=false, double xCenterVal=0.0, double yCenterVal=0.0, String centerUnit=String(""));


    }; // End class ImageManager

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_IMAGEMANAGER_H
