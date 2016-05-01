// $Id: FitsBinTableReader.h,v 1.1.1.1 2010/07/13 17:56:52 eml Exp $

#ifndef GCP_UTIL_FITSBINTABLEREADER_H
#define GCP_UTIL_FITSBINTABLEREADER_H

/**
 * @file FitsBinTableReader.h
 * 
 * Tagged: Wed May  9 15:07:15 PDT 2007
 * 
 * @version: $Revision: 1.1.1.1 $, $Date: 2010/07/13 17:56:52 $
 * 
 * @author Erik Leitch.
 */

#include <iostream>
#include <vector>
#include <string>

#include "gcp/util/DataType.h"
#include "gcp/util/FitsReader.h"

#define FITS_BIN_MAX_DIM 4

namespace gcp {
  namespace util {

    class FitsBinTableReader : public FitsReader {
    public:

      struct Column {
	std::string name_;
	std::string unit_;
	int typecode_;
	long repeat_;
	long width_;
	std::vector<long> naxes_;

	friend std::ostream& operator<<(std::ostream& os, Column& col);

      };

      /**
       * Constructor.
       */
      FitsBinTableReader(std::string file);
      FitsBinTableReader(std::string file, std::string extension);
      FitsBinTableReader(std::string file, int extension);

      FitsBinTableReader(fitsfile* fptr);
      FitsBinTableReader(fitsfile* fptr, std::string extension);

      void setTo(fitsfile* fptr, std::string extension);

      void initialize();
      void getTableInfo(int hdutype);
      void getNextTableInfo();
      void getTableInfoByExtnum(int extnum);
      void getTableInfo(std::string extension);

      void countTables();

      /**
       * Destructor.
       */
      virtual ~FitsBinTableReader();

      unsigned nRow() {
	return (unsigned)nRow_;
      }

      unsigned nCol() {
	return (unsigned)nCol_;
      }

      std::vector<double> getData(unsigned colNo);
      std::vector<double> getData(std::string colName);

      std::vector<bool> getBoolData(unsigned colNo);
      std::vector<unsigned char> getUcharData(unsigned colNo);
      std::vector<char> getCharData(unsigned colNo);
      std::vector<unsigned short> getUshortData(unsigned colNo);
      std::vector<short> getShortData(unsigned colNo);
      std::vector<unsigned int> getUintData(unsigned colNo);
      std::vector<int> getIntData(unsigned colNo);
      std::vector<unsigned long> getUlongData(unsigned colNo);
      std::vector<long> getLongData(unsigned colNo);
      std::vector<float> getFloatData(unsigned colNo);
      std::vector<double> getDoubleData(unsigned colNo);
      std::vector<std::string> getStringData(unsigned colNo);

      std::vector<bool> getBoolData(std::string colName);
      std::vector<unsigned char> getUcharData(std::string colName);
      std::vector<char> getCharData(std::string colName);
      std::vector<unsigned short> getUshortData(std::string colName);
      std::vector<short> getShortData(std::string colName);
      std::vector<unsigned int> getUintData(std::string colName);
      std::vector<int> getIntData(std::string colName);
      std::vector<unsigned long> getUlongData(std::string colName);
      std::vector<long> getLongData(std::string colName);
      std::vector<float> getFloatData(std::string colName);
      std::vector<double> getDoubleData(std::string colName);
      std::vector<std::string> getStringData(std::string colName);

      std::vector<double> getDataAsDouble(unsigned colNo);

      void getData(unsigned colNo, void* array, 
		   void* nulval, int* anynul);

      long colRepeat(unsigned colNo) {
	return columns_[colNo].repeat_;
      }

      int colFitsDataType(unsigned colNo) {
	return columns_[colNo].typecode_;
      }

      std::string colName(unsigned colNo) {
	return columns_[colNo].name_;
      }

      std::string colUnit(unsigned colNo) {
	return columns_[colNo].unit_;
      }

      DataType::Type colDataType(unsigned colNo);

      static DataType::Type colDataType(int colNo);

      FitsBinTableReader::Column* getColumn(std::string name);

      std::vector<long> getColNaxes(std::string name);

      void printColumns();

      unsigned getColNo(std::string name);
      std::string getUnit(std::string name);

      FitsReader::Axis getAxis(std::string valContains);

    private:

      int nKeyword_;
      long int nRow_;
      int nCol_;
      char** tType_;
      char** tForm_;
      char** tUnit_;
      char binName_[100];

      std::vector<Column> columns_;

    }; // End class FitsBinTableReader

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_FITSBINTABLEREADER_H
