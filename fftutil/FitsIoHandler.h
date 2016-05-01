#ifndef GCP_UTIL_FITSIOHANDLER_H
#define GCP_UTIL_FITSIOHANDLER_H

/**
 * @file FitsIoHandler.h
 * 
 * Tagged: Mon May  2 19:04:15 UTC 2005
 * 
 * @author 
 */
#include <string>

#include "gcp/fftutil/ObsInfo.h"

#include "gcp/util/DataType.h"
#include "gcp/util/Declination.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/VisIo.h"

// Functions of the following type are called by put_phdu() to
// convert a passed string value to the appropriate type.

#define WRITEFN(fn) void (fn)(char *buf, const char *val);

// Functions of the following type are called by read_phdu() to
// convert a passed string value to the appropriate type.

#define READFN(fn)  void (fn)(FitsIoHandler::FitsDataType& data, char *str);

namespace gcp {
  namespace util {
    
    class FitsIoHandler : public VisIo {
    public:
      
      /**
       * Constructor.
       */
      FitsIoHandler();
      
      /**
       * Destructor.
       */
      virtual ~FitsIoHandler();
      
      //------------------------------------------------------------
      // Some const members of this class
      //------------------------------------------------------------

      // Length of an ASCII header

      static const unsigned nBytePerHeader_=80;
      
      // Number of headers per data record

      static const unsigned nHeaderPerDataRecord_ = 36;

      // Length of a data record

      static const unsigned nBytePerDataRecord_ = 
	nBytePerHeader_ * nHeaderPerDataRecord_;

      // Length of the header keyword.

      static const unsigned nBytePerKeyword_=8;

      // A logical value 'T' or 'F' for logical keywords must occur in
      // this position in a string:

      static const unsigned logPos_=30;      

      // The length of a single record of the frequency table

      static const unsigned nBytePerFrequencyTableEntry_ = 32;

      // The number of seconds per day

      static const unsigned secondsPerDay_ = 86400;

      // Enumerate recognized Fits axis types.

      enum FitsAxis {
	AX_DEG,        /* Degrees -- type of RA/DEC */
	AX_RAD,
	AX_UV,
	AX_UNKNOWN
      };
      
      // And a container for associating a header string with this
      // unit.
      
      struct FitsAxisCard {
	std::string str_;
	std::string label_;
	FitsIoHandler::FitsAxis axis_;
      };
      
      // A static list of recognized FITS axis types

      static FitsAxisCard fitsAxes_[];
      static unsigned nFitsAxes_;

      // Enumerate recognized Fits bunits.

      enum FitsBunit {
	BU_MJYSR,     /* MJy/sr */
	BU_JYBEAM,    /* Jy/Beam */
	BU_MUK,       /* MicroKelvin */
	BU_UNKNOWN
      };
      
      // And a container for associating a header string with this
      // bunit.

      struct FitsBunitCard {
	std::string str_;
	std::string label_;
	FitsBunit bunit_;
      };

      // A struct for handling all possible FITS data types

      struct FitsDataType {
	gcp::util::DataType stdVal_;
	std::string stringVal_;
	FitsIoHandler::FitsBunit bunitVal_;
	FitsIoHandler::FitsAxis axisVal_;
      };

      // A static list of recognized FITS axis types

      static FitsBunitCard fitsUnits_[];
      static unsigned nFitsUnits_;

      // Declare a container for a single header card.

      struct Phdu {
	char name_[FitsIoHandler::nBytePerKeyword_+1];
	READFN(*readfn_);
	WRITEFN(*writefn_);
	int required_;
      };

      static FitsIoHandler::Phdu phdus_[];
      static unsigned nPhdus_;

      // Write a header key to a file.

      void putPhdu(const char *name, const char *val, const char *comment, 
		   FILE *fp=0);

      // Methods for formatting various FITS keyword types

      static void nullStr(char *buf, const char *val);
      static void  logStr(char *buf, const char *val);
      static void  intStr(char *buf, const char *val);
      static void  fltStr(char *buf, const char *val);
      static void  strStr(char *buf, const char *val);

      // Methods for reading various FITS keyword types

      static void  rdNull(FitsDataType& data, char *str);
      static void   rdLog(FitsDataType& data, char *str);
      static void   rdStr(FitsDataType& data, char *str);
      static void   rdInt(FitsDataType& data, char *str);
      static void   rdFlt(FitsDataType& data, char *str);
      static void rdBunit(FitsDataType& data, char *str);
      static void  rdAxis(FitsDataType& data, char *str);

      //------------------------------------------------------------
      // Methods to write headers
      //------------------------------------------------------------

      void setDate();

      // All headers will be initialized and finished the same way

      void initHeader();
      void finishHeader(FILE* fp=0);

      // Methods to write the bodies for different types of header

      void writeFileHeaderBody(FILE* fp=0);
      void writeAntennaTableHeaderBody(FILE* fp=0);
      void writeFrequencyTableHeaderBody(FILE* fp=0);

      // Methods to write the whole header, for different header types

      void writeFileHeader(FILE* fp=0);
      void writeAntennaTableHeader(FILE* fp=0);
      void writeFrequencyTableHeader(FILE* fp=0);

      //------------------------------------------------------------
      // Methods to write the data segment
      //------------------------------------------------------------

      void initVisibilityData();
      void installVisibilityData(double* vis, double* date, double* uvw);

      void writeVisibilityDataBody(FILE* fp=0);
      void writeVisibilityDataBody(ObsInfo& obs);
      void writeVisibilityData(FILE* fp=0);
      void writeVisibilityData(ObsInfo& obs);

      void writeFakeVisibilityData(double* vis, double* date, double* uvw, FILE* fp=0);
      void writeFakeVisibilityDataBody(double* vis, double* date, double* uvw, FILE* fp=0);

      void finishVisibilityData(FILE* fp=0);

      //------------------------------------------------------------
      // Methods to write the antenna table
      //------------------------------------------------------------

      void initAntennaTableDataBody();
      void writeAntennaTableDataBody(FILE* fp=0, ObsInfo* obs=0);
      void finishAntennaTableDataBody(FILE* fp=0);

      // Write the data for an antenna table

      void writeAntennaTableData(FILE* fp=0, ObsInfo* obs=0);
      void writeAntennaTableEntry(int i, double X, double Y, double Z, FILE* fp=0);
      void writeAntennaTableEntry(int i, std::string anname, double X, double Y, double Z, FILE* fp=0);

      // Write the whole antenna table

      void writeAntennaTable(FILE* fp=0);
      void writeAntennaTable(ObsInfo& obs);

      //------------------------------------------------------------
      // Methods to write the frequency table
      //------------------------------------------------------------

      void initFrequencyTableDataBody();
      void writeFrequencyTableDataBody(FILE* fp=0);
      void finishFrequencyTableDataBody(FILE* fp=0);

      void writeFrequencyTableData(FILE* fp=0);
      void writeFrequencyTable(FILE* fp=0);

      //-----------------------------------------------------------------------
      // Write a UVF file
      //-----------------------------------------------------------------------

      void writeUvfFile(FILE* fp=0);
      void writeUvfFile(ObsInfo& obs);

      // Write a file based on an external observation info object

      void writeUvfFile(std::string fileName, ObsInfo& obs);

      void writeFakeUvfFile(double* data, double* date, double* uvw, FILE* fp=0);

      // A method to return whichever is the current file descriptor

      FILE* getFptr(FILE* fp);

      void openFile(std::string fileName);
      void closeFile();

      // Inheritors should define these accessor methods for writing UVF files

    protected:

      virtual float getU(unsigned iFrame, unsigned iBaseline);
      virtual float getV(unsigned iFrame, unsigned iBaseline);
      virtual float getW(unsigned iFrame, unsigned iBaseline);
      virtual double getJulianDate(unsigned iFrame);
      virtual float getRe(unsigned iFrame, unsigned iBaseline, unsigned iIf, unsigned iStokes);
      virtual float getIm(unsigned iFrame, unsigned iBaseline, unsigned iIf, unsigned iStokes);
      virtual float getWt(unsigned iFrame, unsigned iBaseline, unsigned iIf, unsigned iStokes);
      virtual float getAipsBaselineIndex(unsigned iFrame, unsigned iBaseline);

    private:

      // If true, swap bytes to write from native Endianness to FITS bi-Endianness

      bool swapBytes_;

      // Internal counters we will use when writing headers and data

      unsigned nHdu_;
      unsigned nDataByte_;

      //------------------------------------------------------------
      // Some information needed to write a FITS file
      //------------------------------------------------------------

      FILE* fp_;
 
      // Temporary handles uesd by the base-class

      double* visPtr_;
      double* uvwPtr_;
      double* datePtr_;

      // Utility functions for converting from one Endianness to the
      // other

      void cp4r4(unsigned char *dest, unsigned char *orig, size_t nitem);
      
      // Orig: 8-byte datatype.  Dest: 8-byte datatype with byte
      // order reversed.
      
      void cp8r8(unsigned char *dest, unsigned char *orig, size_t nitem);

      void fwrite(unsigned int* ptr, size_t nel, FILE* stream);
      void fwrite(float* ptr,        size_t nel, FILE* stream);
      void fwrite(double* ptr,       size_t nel, FILE* stream);

    }; // End class FitsIoHandler
    
  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_FITSIOHANDLER_H
