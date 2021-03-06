// $Id: MiriadIo.h,v 1.2 2011/11/01 21:57:53 eml Exp $

#ifndef GCP_UTIL_MIRIADIO_H
#define GCP_UTIL_MIRIADIO_H

/**
 * @file MiriadIo.h
 * 
 * Tagged: Mon Oct  3 15:31:19 PDT 2005
 * 
 * @version: $Revision: 1.2 $, $Date: 2011/11/01 21:57:53 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Angle.h"
#include "gcp/util/Pressure.h"
#include "gcp/util/Speed.h"
#include "gcp/util/Temperature.h"
#include "gcp/util/VisIo.h"

namespace gcp {
  namespace util {

    class MiriadIo : public VisIo {
    public:

      /**
       * Constructor.
       */
      MiriadIo();

      /**
       * Destructor.
       */
      virtual ~MiriadIo();

      // Versions of the base-class routines which handle Miriad files

      void openFile(std::string name, std::string openMode);
      void openFileForRead(std::string name);
      void closeFile();
      void writeFile(double* data, double* date, double* uvw, double* rms);
      void writeFakeFile(double* data, double* date, double* uvw, double* rms);
      void getStats(std::string name, unsigned& nRecord, unsigned& minVisPerRecord, unsigned& maxVisPerRecord);

      // Method to write a single frame of data, including uv
      // variables written by CARMA
      
      // Write antenna pointing information

      void writeAntennaPointing();
      
      // Write weather parameters

      void writeWeatherParameters();
      void writeTimeParameters(unsigned iFrame);
      void writeCarmaFormatData(unsigned iFrame);

      void resetVisStats();
      void reportVisStats();

      void setVersion(std::string version);

      void setFloatOutputFormat();

      //------------------------------------------------------------
      // Methods for reading miriad files
      //------------------------------------------------------------

      void readMiriadFile(std::string name);
      void readNextRecord();
      void parseHeader();

      void parseArrayInfo();
      void parseAntennaLocations();
      void parseSourceInfo();

      void accessVisData();

      // Return true when at the end of the file

      bool atEnd();

      // Variables items

      std::string readStringVar(std::string name);

      std::vector<float> readFloatVar(std::string name, unsigned n);
      float readFloatVar(std::string name);

      std::vector<double> readDoubleVar(std::string name, unsigned n);
      double readDoubleVar(std::string name);

      std::vector<int> readIntVar(std::string name, unsigned n);
      int readIntVar(std::string name);

      void getTypeAndLength(std::string name, char& type, int& len);
      bool variableExists(std::string name);

      // Header items

      std::string readStringHdrItem(std::string name);

      float readFloatHdrItem(std::string name);

      double readDoubleHdrItem(std::string name);

      int readIntHdrItem(std::string name);

    public:
      
      std::string version_;

      // A UV handle used by Miriad writing routines

      int uvh_;

      unsigned nCarma_;

      unsigned iFirstGoodChannel_;
      unsigned nGoodChannel_;

      void writeFixedParameters();
      void writeVisibilityData(double* data, double* date, double* uvw, 
			       double* rms);

      void writeVisibilityData(unsigned iFrame);
      void writeWidebandVisibilityData(unsigned iFrame, unsigned iBaseline);
      void writeSpectralVisibilityData(unsigned iFrame, unsigned iBaseline);

      void writeFakeVisibilityData(double* data, double* date, double* uvw, 
				   double* rms);

      void writeSourceParameters();
      void writeSiteParameters();
      void writeAntennaParameters();
      void writeArrayParameters();
      void writeLinelengthParameters();
      void writeSysTemps(double* rmsInJy, unsigned iFrame);

      bool haveVisWideData();
      bool haveVisSpecData();

      unsigned nVis_;
      unsigned goodData_;
      unsigned badData_;
      unsigned badFlagData_;
      unsigned badWtData_;

      // Variables used for reading in data

      std::vector<double> preamble_;
      std::vector<float>  readVisData_;
      std::vector<int>    readFlags_;

    public:

      int recordLen_;
      bool recordLenIsInitialized_;

      int recordNo_;
      int nVisPerRecord_;
      int nVisLastRead_;
      bool noDataReadYet_;
      unsigned maxVisPerRecord_;
      unsigned minVisPerRecord_;

    }; // End class MiriadIo

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_MIRIADIO_H
