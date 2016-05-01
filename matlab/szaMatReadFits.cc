/**.......................................................................
 * MATLAB Mex file for reading from UVF files
  *
 * Use like:
 *
 * d=gcpMatReadUvf({'array.frame.record','corr.band0.usb[0][0]',
 *                     'antenna*.tracker.actual double', 'antenna*.tracker.source string'},
 *                     '06-jan-2005:15','06-jan-2005:16',
 *                     '/data/gcpdaq/arc','/home/gcpdaq/carma_unstable/gcp/array/conf/cal');
 *
 */
#include "gcp/util/FitsImageReader.h"
#include "gcp/util/FitsBinTableReader.h"

#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

#include "mex.h"
#include "matrix.h"

#include <iostream.h>
#include <math.h>

using namespace std;
using namespace gcp::util;
using namespace gcp::matlab;

#define ASSIGN_VALS(cast) \
  {\
    cast* mptr = (cast*)data;\
    unsigned iM, iC;\
    unsigned nEl0 = reader.nRow();\
    unsigned nEl1 = reader.colRepeat(iCol);\
    for(unsigned iEl1=0; iEl1 < nEl1; iEl1++)\
      for(unsigned iEl0=0; iEl0 < nEl0; iEl0++) {\
         iC = iEl0*nEl1 + iEl1;\
         iM = iEl1*nEl0 + iEl0;\
         *(mptr + iM) = vals[iC];\
      }\
  }

mxArray* readFitsFile(const mxArray* prhs[]);
void addAxes(mxArray* hdr, FitsImageReader& reader);
mxArray* addHeader(mxArray* hdr, std::string& file, FitsImageReader& reader);
void addTables(mxArray* hdr, std::string fileName);
void addTable(mxArray* tables, FitsBinTableReader& tableReader, std::string extname, unsigned index);
void assignTableData(FitsBinTableReader& reader, unsigned iCol, void* data);
void addData(mxArray* data, FitsImageReader& reader);

/**.......................................................................
 * Entry point from the matlab environment
 */
void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  gcp::util::Logger::installStdoutPrintFn(MexHandler::stdoutPrintFn);
  gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);
  
  gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
  gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
  gcp::util::ErrHandler::installLogFn(MexHandler::logFn);
  
  // Check input/output arguments
  
  //  MexHandler::checkArgs(nlhs, nrhs, 1, 5);
  
  plhs[0] = readFitsFile(prhs);
}
   
/**.......................................................................
 * Read the FITS image file
 */
mxArray* readFitsFile(const mxArray* prhs[])
{
  mxArray* ptr=0;

  // Open the file

  char* cptr = mxArrayToString(prhs[0]);
  std::string fName(cptr);

  gcp::util::FitsImageReader reader(fName);

  // Create a top-level struct
  
  int dims[2]={1,1};
  ptr = mxCreateStructArray(2, dims, 0, NULL);

  mxArray* hdr = addHeader(ptr, fName, reader);

  // Add axis struct array to the header

  addAxes(hdr, reader);

  // Finally, read the data!

  addData(ptr, reader);

  return ptr;
}

/**.......................................................................
 * Create two substructures
 */
mxArray* addHeader(mxArray* ptr, std::string& file, gcp::util::FitsImageReader& reader)
{
  // Create two substructures

  mxArray* hdr  = MexHandler::addNamedStructField(ptr, "header");

  // Add the file to it

  MexHandler::addNamedStringStructField(hdr, "file", file);

  return hdr;
}

/**.......................................................................
 * Add axis information to the returned structure
 */
void addAxes(mxArray* hdr, FitsImageReader& reader)
{
  unsigned* naxis = MexHandler::addNamedUintStructField(hdr,   "naxis", 1);
  *naxis = reader.nAxis();

  mxArray* axis = 
    MexHandler::addNamedStructField(hdr,   "axes", reader.nAxis());

  // For each element of the array, add pertinent fields

  unsigned imAxis=0;
  for(int iAxis=reader.nAxis()-1; iAxis >= 0; iAxis--, imAxis++) {

    MexHandler::addNamedStringStructField(axis,       "name", 
					  reader.axisName(iAxis), imAxis);

    MexHandler::addNamedStringStructField(axis,       "comment", 
					  reader.axisComment(iAxis), imAxis);

    unsigned* n = MexHandler::addNamedUintStructField(axis,   "n", 1, imAxis);
    *n = reader.axisSize(iAxis);

    double* refVal = MexHandler::addNamedDoubleStructField(axis,   "refVal", 1, imAxis);
    *refVal = reader.axisRefVal(iAxis);

    double* refPix = MexHandler::addNamedDoubleStructField(axis,   "refPix", 1, imAxis);
    *refPix = reader.axisRefPix(iAxis);

    double* delta = MexHandler::addNamedDoubleStructField(axis,   "delta", 1, imAxis);
    *delta = reader.axisDelta(iAxis);

    double* data = MexHandler::addNamedDoubleStructField(axis,   "data", reader.axisSize(iAxis), imAxis);
    reader.getAxisData(data, iAxis);
  }    
}

/**.......................................................................
 * Add table information to the returned structure
 */
void addTables(mxArray* hdr, std::string fileName)
{
  FitsBinTableReader tableReader(fileName);

  // Create an array of tables

  unsigned* ntable = MexHandler::addNamedUintStructField(hdr,   "ntable", 1);
  *ntable = 2;

  mxArray* tables = 
    MexHandler::addNamedStructField(hdr, "tables", *ntable);

  // Move to the next table
  
  tableReader.getTableInfo("AIPS AN");
  addTable(tables, tableReader, "AIPS AN", 0);

  tableReader.getTableInfo("AIPS FQ");
  addTable(tables, tableReader, "AIPS FQ", 1);
}

/**.......................................................................
 * Add table information to the returned structure
 */
void addTable(mxArray* tables, FitsBinTableReader& tableReader, 
	      std::string extname, unsigned index)
{
  // Each table will have a name

  MexHandler::addNamedStringStructField(tables, "name", extname, index);

  // Each table will record the number of columns

  unsigned* ncol = MexHandler::addNamedUintStructField(tables,   "ncol", 1, index);
  *ncol = tableReader.nCol();

  // Each table will have an array of columns
  
  mxArray* cols = 
    MexHandler::addNamedStructField(tables, "cols", tableReader.nCol(), index);
  
  // Each column will have the same fields

  for(int iCol=tableReader.nCol()-1; iCol >= 0; iCol--) {
    
    MexHandler::addNamedStringStructField(cols,       "name", 
					  tableReader.colName(iCol), iCol);
    
    unsigned* n = MexHandler::addNamedUintStructField(cols,   "n", 1, iCol);
    *n = tableReader.colRepeat(iCol);
    
    // Add the data field

    int anynul=0;
    void* data = MexHandler::addNamedStructField(cols, "data", &tableReader, iCol, iCol);

    // Now fill it with data from the FITS file.  Ignore fields with 0
    // repeat count -- I don't know what this means, or what to do
    // with it.

    if(tableReader.colRepeat(iCol) == 0)
      continue;

    if(tableReader.colDataType(iCol) != DataType::STRING) {

      // If this is numeric data, just create numeric arrays

      assignTableData(tableReader, iCol, data);

    } else {

      // If these are strings, create a cell array

      std::vector<string> vals = tableReader.getStringData(iCol);

      mxArray* array = mxGetField(cols, iCol, "data");

      for(unsigned i=0; i < vals.size(); i++)
	mxSetCell(array, i, mxCreateString(vals[i].c_str()));

    }
  }  
}

/**.......................................................................
 * Crap function to assign matlab arrays from returned FITS table arrays
 */
void assignTableData(FitsBinTableReader& reader, unsigned iCol, void* data)
{
  switch (reader.colDataType(iCol)) {
  case DataType::BOOL:
    {
      std::vector<bool> vals = reader.getBoolData(iCol);
      ASSIGN_VALS(bool);
    }
    break;
  case DataType::CHAR:
    {
      std::vector<char> vals = reader.getCharData(iCol);
      ASSIGN_VALS(char);
    }
    break;
  case DataType::UCHAR:
    {
      std::vector<unsigned char> vals = reader.getUcharData(iCol);
      ASSIGN_VALS(unsigned char);
    }
    break;
  case DataType::SHORT:
    {
      std::vector<short> vals = reader.getShortData(iCol);
      ASSIGN_VALS(short);
    }
    break;
  case DataType::USHORT:
    {
      std::vector<unsigned short> vals = reader.getUshortData(iCol);
      ASSIGN_VALS(unsigned short);
    }
    break;
  case DataType::INT:
    {
      std::vector<int> vals = reader.getIntData(iCol);
      ASSIGN_VALS(int);
    }
    break;
  case DataType::UINT:
    {
      std::vector<unsigned int> vals = reader.getUintData(iCol);
      ASSIGN_VALS(unsigned int);
    }
    break;
  case DataType::LONG:
    {
      std::vector<long> vals = reader.getLongData(iCol);
      ASSIGN_VALS(long);
    }
    break;
  case DataType::ULONG:
    {
      std::vector<unsigned long> vals = reader.getUlongData(iCol);
      ASSIGN_VALS(unsigned long);
    }
    break;
  case DataType::FLOAT:
    {
      std::vector<float> vals = reader.getFloatData(iCol);
      ASSIGN_VALS(float);
    }
    break;
  case DataType::DOUBLE:
    {
      std::vector<double> vals = reader.getDoubleData(iCol);
      ASSIGN_VALS(double);
    }
    break;
  default:
    break;
  }
}

void addData(mxArray* data, FitsImageReader& reader)
{
  std::vector<unsigned> mindices;

  // Read the data

  FitsImageReader::Image image;
  reader.readData(image);
  
  // Add the image data field
  
  std::vector<int> dims = reader.axisDims();
  MexHandler::getIndicesMatlab(mindices, dims);

  float* imArr = MexHandler::addNamedFloatStructField(data, "image", dims);

  for(unsigned i=0; i < mindices.size(); i++) {
    imArr[i] = image.data_[mindices.size()-mindices[i]];
    //   imArr[i] = image.data_[mindices.size()-i];
  }
}
