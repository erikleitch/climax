// $Id: $

#ifndef GCP_UTIL_XRAYIMAGE_H
#define GCP_UTIL_XRAYIMAGE_H

/**
 * @file XrayImage.h
 * 
 * Tagged: Mon Feb  9 09:43:16 PST 2015
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include "gcp/fftutil/Image.h"

namespace gcp {
  namespace util {

    class XrayImage : public Image {
    public:

      /**
       * Constructor.
       */
      XrayImage();

      /**
       * Destructor.
       */
      virtual ~XrayImage();

    private:
    }; // End class XrayImage

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_XRAYIMAGE_H
