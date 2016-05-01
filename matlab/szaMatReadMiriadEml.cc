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
  
  plhs[0] = readMiriadFile(prhs);
}
   
/**.......................................................................
 * Read the Miriad file
 */
mxArray* readMiriadFile(const mxArray* prhs[])
{
  mxArray* ptr=0;

  // Open the file

  char* cptr = mxArrayToString(prhs[0]);
  std::string fName(cptr);

  gcp::util::MiriadIo io;

  // Get stats we need to read the file

  unsigned nGroup, minVisPerRecord, maxVisPerRecord;
  io.getStats(nGroup, minVisPerRecord, maxVisPerRecord);

  // Create a top-level struct
  
  int dims[2]={1,1};
  ptr = mxCreateStructArray(2, dims, 0, NULL);

  // Create two substructures

  mxArray* hdr  = MexHandler::addNamedStructField(ptr, "header");

  //------------------------------------------------------------
  // Add header fields
  //------------------------------------------------------------

  MexHandler::addNamedStringStructField(hdr,       "file", fName);

  unsigned int* ngroup = MexHandler::addNamedUintStructField(hdr,       "ngroup", 1);
  *ngroup = nGroup;

  // Finally, read the data!

  addData(ptr, fName);

  return ptr;
}

void addData(mxArray* data, std::string fName, unsigned maxVisPerRecord)
{
  FitsUvfReader::Vis vis;
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
    unsigned* base  = MexHandler::addNamedUintStructField(groups,     "base", 1, iGroup);
    float* jd  = MexHandler::addNamedFloatStructField(groups,     "jd", 1, iGroup);

    // Add an array with the same dimensionality as the group array in
    // the Miriad file

    mxArray* visArr = MexHandler::addNamedStructField(groups, "vis", dims.size(), (const int*)&dims[0], DataType::COMPLEX_FLOAT, iGroup);
    mxArray* wtArr  = MexHandler::addNamedStructField(groups,  "wt", dims.size(), (const int*)&dims[0], DataType::FLOAT, iGroup);

    float* rePtr = (float*)mxGetData(visArr);
    float* imPtr = (float*)mxGetImagData(visArr);
    float* wtPtr = (float*)mxGetData(wtArr);

    // Now read the data

    reader.readData(iGroup, vis);

    // Set the group data

    *u  = vis.u_;
    *v  = vis.v_;
    *w  = vis.w_;
    *jd = vis.jd_;
    *base = vis.baseline_;

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
