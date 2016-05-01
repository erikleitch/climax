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
#include "gcp/util/FitsUvfReader.h"
#include "gcp/util/FitsBinTableReader.h"

#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

#include "mex.h"
#include "matrix.h"

#include <iostream>
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

mxArray* readUvfFile(const mxArray* prhs[]);
void addAxes(mxArray* hdr, FitsUvfReader& reader);
void addTables(mxArray* hdr, std::string fileName);
void addTable(mxArray* tables, FitsBinTableReader& tableReader, std::string extname, unsigned index);
void assignTableData(FitsBinTableReader& reader, unsigned iCol, void* data);
void addData(mxArray* data, FitsUvfReader& reader);

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
  
  plhs[0] = readUvfFile(prhs);
}
   
/**.......................................................................
 * Read the UVF file
 */
mxArray* readUvfFile(const mxArray* prhs[])
{
  mxArray* ptr=0;

  // Open the file

  char* cptr = mxArrayToString(prhs[0]);
  std::string fName(cptr);

  gcp::util::FitsUvfReader reader(fName);

  // Create a top-level struct
  
  int dims[2]={1,1};
  ptr = mxCreateStructArray(2, dims, 0, NULL);

  // Create two substructures

  mxArray* hdr  = MexHandler::addNamedStructField(ptr, "header");

  // Add header fields

  MexHandler::addNamedStringStructField(hdr,       "file",     
					fName);

  MexHandler::addNamedStringStructField(hdr,       "object",     
					reader.object());

  double* ra = MexHandler::addNamedDoubleStructField(hdr,       "obsra", 1);
  *ra = reader.obsra().hours();

  double* dec = MexHandler::addNamedDoubleStructField(hdr,       "obsdec", 1);
  *dec = reader.obsdec().degrees();

  MexHandler::addNamedStringStructField(hdr,       "obsdate",     
					reader.dateObs());

  MexHandler::addNamedStringStructField(hdr,       "telescope",  
					reader.telescope());

  MexHandler::addNamedStringStructField(hdr,       "instrument", 
					reader.instrument());
  
  unsigned int* ngroup = MexHandler::addNamedUintStructField(hdr,       "ngroup", 1);
  *ngroup = reader.nGroup();

  // Add axis struct array to the header

  addAxes(hdr, reader);
  addTables(hdr, fName);

  // Finally, read the data!

  addData(ptr, reader);

  return ptr;
}

/**.......................................................................
 * Add axis information to the returned structure
 */
void addAxes(mxArray* hdr, FitsUvfReader& reader)
{
  unsigned* naxis = MexHandler::addNamedUintStructField(hdr,   "naxis", 1);
  *naxis = reader.nAxis()-1;

  mxArray* axis = 
    MexHandler::addNamedStructField(hdr,   "axes", reader.nAxis()-1);

  // For each element of the array, add pertinent fields

  unsigned imAxis=0;
  for(int iAxis=reader.nAxis()-1; iAxis > 0; iAxis--, imAxis++) {

    MexHandler::addNamedStringStructField(axis,       "name", 
					  reader.axisName(iAxis), imAxis);

    MexHandler::addNamedStringStructField(axis,       "comment", 
					  reader.axisComment(iAxis), imAxis);

    unsigned* n = MexHandler::addNamedUintStructField(axis,   "n", 1, imAxis);
    *n = reader.axisSize(iAxis);

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
  *ntable = tableReader.nTable();

  mxArray* tables = 
    MexHandler::addNamedStructField(hdr, "tables", *ntable);

  bool stop = false;

  unsigned iTable=0;
  while(!stop) {
    try {
      tableReader.getNextTableInfo();
      addTable(tables, tableReader, tableReader.tableName(), iTable);
      ++iTable;
    } catch(Exception& err) {
      stop = true;
    } catch(...) {
      stop = true;
    }
  }
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

    if(tableReader.colRepeat(iCol) == 0) {
      continue;
    }

  if(tableReader.colDataType((unsigned)iCol) != DataType::STRING) {

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

void addData(mxArray* data, FitsUvfReader& reader)
{
  FitsUvfReader::Vis vis;
  std::vector<unsigned> dims = reader.groupAxisDimsMatlabOrder();
  std::vector<unsigned> mindices;
  bool first=true;

  // Add the groups data field

  mxArray*  groups = MexHandler::addNamedStructField(data, "groups", reader.nGroup());

  // Now iterate over all groups in the file, reading each one into
  // the corresponding element of the matlab array

  for(unsigned iGroup=0; iGroup < reader.nGroup(); iGroup++) {

    // Add singleton data fields for each of the following

    float*  u  = MexHandler::addNamedFloatStructField(groups,      "u", 1, iGroup);
    float*  v  = MexHandler::addNamedFloatStructField(groups,      "v", 1, iGroup);
    float*  w  = MexHandler::addNamedFloatStructField(groups,      "w", 1, iGroup);
    unsigned* a1     = MexHandler::addNamedUintStructField(groups,     "ant1", 1, iGroup);
    unsigned* a2     = MexHandler::addNamedUintStructField(groups,     "ant2", 1, iGroup);
    unsigned* base   = MexHandler::addNamedUintStructField(groups,     "base", 1, iGroup);
    unsigned* source = MexHandler::addNamedUintStructField(groups,     "source", 1, iGroup);
    float* jd  = MexHandler::addNamedFloatStructField(groups,     "jd", 1, iGroup);


    // Add an array with the same dimensionality as the group array in
    // the FITS file

    mxArray* visArr = MexHandler::addNamedStructField(groups, "vis", dims.size(), (const int*)&dims[0], DataType::COMPLEX_FLOAT, iGroup);
    mxArray* wtArr  = MexHandler::addNamedStructField(groups,  "wt", dims.size(), (const int*)&dims[0], DataType::FLOAT, iGroup);

    float* rePtr = (float*)mxGetData(visArr);
    float* imPtr = (float*)mxGetImagData(visArr);
    float* wtPtr = (float*)mxGetData(wtArr);

    // Now read the data

    reader.readData(iGroup, vis);

    // Set the group data

    *u      = vis.u_;
    *v      = vis.v_;
    *w      = vis.w_;
    *jd     = vis.jd_;
    *base   = vis.baseline_;
    *source = vis.source_;

    *a1 = (*base) / 256;
    *a2 = (*base) - (*a1)*256;

    // Now set the visibility data.  If this is the first call,
    // establish the matlab order of the indices

    if(first) {
      MexHandler::getIndicesMatlab(mindices, dims);
      first = false;
    }

    for(unsigned i=0; i < mindices.size(); i++) {
      rePtr[i] = vis.re_[mindices[i]];
      imPtr[i] = vis.im_[mindices[i]];
      wtPtr[i] = vis.wt_[mindices[i]];
    }

  }
}
