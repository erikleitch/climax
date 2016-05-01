// $Id: FitsReader.h,v 1.2 2010/07/13 21:17:14 eml Exp $

#ifndef GCP_UTIL_FITSREADER_H
#define GCP_UTIL_FITSREADER_H

/**
 * @file FitsReader.h
 * 
 * Tagged: Mon Jun  4 15:49:15 PDT 2007
 * 
 * @version: $Revision: 1.2 $, $Date: 2010/07/13 21:17:14 $
 * 
 * @author Erik Leitch.
 */
#ifndef ThrowFitsError
#define ThrowFitsError(arg) \
{\
  char errtext[31];\
  ffgerr(status_, errtext);\
  ThrowError(arg << " (" << errtext << " (" << status_ << "))");\
}
#endif

#ifndef ThrowSimpleFitsError
#define ThrowSimpleFitsError(arg) \
{\
  char errtext[31];\
  ffgerr(status_, errtext);\
  ThrowSimpleError(arg << " (" << errtext << " (" << status_ << "))");\
}
#endif

#ifndef ThrowSimpleFitsErrorColor
#define ThrowSimpleFitsErrorColor(arg, color)		\
{\
  char errtext[31];\
  ffgerr(status_, errtext);\
  ThrowSimpleColorError(arg << " (" << errtext << " (" << status_ << "))", color); \
}
#endif

#ifndef ReportFitsError
#define ReportFitsError(arg) \
{\
  char errtext[31];\
  ffgerr(status_, errtext);\
  ReportError(arg << " (" << errtext << " (" << status_ << "))");\
}
#endif

#ifndef ReportSimpleFitsError
#define ReportSimpleFitsError(arg) \
{\
  char errtext[31];\
  ffgerr(status_, errtext);\
  ReportSimpleError(arg << " (" << errtext << " (" << status_ << "))");\
}
#endif

#ifndef ReportSimpleFitsErrorColor
#define ReportSimpleFitsErrorColor(arg, color)		\
{\
  char errtext[31];\
  ffgerr(status_, errtext);\
  ReportSimpleErrorColor(arg << " (" << errtext << " (" << status_ << "))", color); \
}
#endif

#include <string>

#include "gcp/util/Declination.h"
#include "gcp/util/Exception.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/ObsParameter.h"

#include "fitsio.h"

namespace gcp {
  namespace util {

    class FitsReader {
    public:

      struct Axis {
	long n_;
	std::string type_;
	std::string typeComment_;
	std::string unit_;
	double refVal_;
	double refPix_;
	double delta_;
	unsigned iCol_;
	bool isPresent_;

	friend std::ostream& operator<<(std::ostream& os, const FitsReader::Axis& axis);
      };

      /**
       * Constructor.
       */
      FitsReader();

      /**
       * Destructor.
       */
      virtual ~FitsReader();

      unsigned nHdu() {
	return nHdu_;
      }

      unsigned nTable() {
	return nTable_;
      }

      void open(std::string file);
      void close();

      std::string getStringKey(std::string name);
      long getLongKey(std::string name);
      int  getIntKey(std::string name);

      int moveToHduOffsetBy(unsigned offset);
      void moveToMainHdu();
      int moveForwardToTable(std::string extension);
      int moveToTable(std::string extension);
      int moveToNextTable();
      bool hasTable(std::string tableName);

      virtual void getHeaderInfo();

      std::string tableName() {
	return tableName_;
      }

      std::string getKeyword(std::string extname, std::string keyword);
      std::string getKeyword(std::string keyword);

      std::string telescope();
      std::string instrument();
      std::string dateObs();
      std::string object();
      HourAngle obsra();
      Declination obsdec();
      double equinox();

    protected:

      std::string fileName_;
      fitsfile*   fptr_;
      int         nHdu_;
      unsigned    hduOffset_; // Offset of the current HDU from the main HDU;
      unsigned    nTable_;
      bool        wasOpened_;
      int         status_;
      std::string tableName_;

      Declination obsDec_;
      HourAngle obsRa_;
      double equinox_;
      std::string instrument_;
      std::string telescope_;
      std::string object_;
      std::string dateObs_;

      unsigned infoMask_;

    private:

      void initializeInfo();
      void countHdus();
      void countTables();

    }; // End class FitsReader

    std::ostream& operator<<(std::ostream& os, const FitsReader::Axis& axis);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_FITSREADER_H
