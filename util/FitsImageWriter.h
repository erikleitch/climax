// $Id: FitsImageWriter.h,v 1.1 2012/05/08 21:58:00 eml Exp $

#ifndef GCP_UTIL_FITSIMAGEWRITER_H
#define GCP_UTIL_FITSIMAGEWRITER_H

/**
 * @file FitsImageWriter.h
 * 
 * Tagged: Mon May  7 11:11:21 PDT 2012
 * 
 * @version: $Revision: 1.1 $, $Date: 2012/05/08 21:58:00 $
 * 
 * @author Erik Leitch
 */
#include "fitsio.h"

#include <string>
#include <valarray>
#include <vector>

namespace gcp {
  namespace util {

    class FitsImageWriter {
    public:

      /**
       * Constructors
       */
      FitsImageWriter();
      FitsImageWriter(std::string fileName);

      /**
       * Destructor.
       */
      virtual ~FitsImageWriter();

      void initialize();

      void open(std::string fileName);
      void close();

      // Write data for this image

      void writeData(std::vector<float>& data);
      void writeData(std::valarray<float>& data);
      void writeStandardKeys();
      void writeStandardKeys(int bitsPerPixel, std::vector<long>& axisLengths);
      void writeImageKeys(int bitsPerPixel, std::vector<long>& axisLengths);

      void putKey(std::string name, std::string val, std::string comment);
      void putKey(std::string name, bool val, std::string comment);
      void putKey(std::string name, char val, std::string comment);
      void putKey(std::string name, short val, std::string comment);
      void putKey(std::string name, int val, std::string comment);
      void putKey(std::string name, long val, std::string comment);
      void putKey(std::string name, float val, std::string comment);
      void putKey(std::string name, double val, std::string comment);

      // Methods for setting required header items

      void setBitsPerPixel(unsigned bitPix);
      void setAxes(std::vector<long>& axisLength);

    private:

      std::string fileName_;
      fitsfile* fptr_;
      int status_;
      bool wasOpened_;

      int bitPix_;
      long nAxis_;
      std::vector<long> axisLengths_;

      bool hasBitPix_;
      bool hasAxes_;

    }; // End class FitsImageWriter

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_FITSIMAGEWRITER_H
