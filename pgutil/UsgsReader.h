// $Id: $

#ifndef GCP_UTIL_USGSREADER_H
#define GCP_UTIL_USGSREADER_H

/**
 * @file UsgsReader.h
 * 
 * Tagged: Tue Sep  1 14:36:14 PDT 2009
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author tcsh: username: Command not found.
 */
#include <string>
#include <vector>

#include "gcp/util/Angle.h"

#include "gcp/fftutil/Image.h"

namespace gcp {
  namespace util {

    class UsgsReader {
    public:

      enum Type {
	TYPE_UNKNOWN = 0x0,
	TYPE_FLT     = 0x1,
	TYPE_DEM     = 0x2,
      };

      /**
       * Constructor.
       */
      UsgsReader();

      /**
       * Destructor.
       */
      virtual ~UsgsReader();

      enum ByteOrder {
	LSB_FIRST,
	USB_FIRST
      };

      void setTo(std::string directory, std::string prefix);
      void read(std::string directory, std::string prefix);
      void read();

      bool initialized();
      void parseHeaderInfo();
      void parseFltHeaderInfo();
      void parseDemHeaderInfo();
      void readImageData();
      void display();

      void getNextKeywordPair(String& str, String& tok, String& val);
      void reverseBytes(short& s);

      std::string getDataFile(std::string directory, std::string prefix);
      bool needsSwap();

    public:

      bool isInitialized_;
      bool imageReceived_;

      ByteOrder byteOrder_;

      Image image_;

      unsigned nCols_;
      unsigned nRows_;

      Angle xllCorner_;
      Angle yllCorner_;
      Angle cellSize_;

      double noDataVal_;
      unsigned nBitPerPixel_;

      Type type_;

      std::string directoryPrefix_;
      std::string hdrFileName_;
      std::string imageFileName_;

      std::vector<float> imageData_;

    }; // End class UsgsReader

  } // End namespace util
} // End namespace gcp

#endif // End #ifndef GCP_UTIL_USGSREADER_H
