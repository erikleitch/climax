/**.......................................................................
 * MATLAB Mex file for reading from GCP archive data files.
 *
 * Use like:
 *
 * d=gcpMatFastReadArc({'array.frame.record','corr.band0.usb[0][0]',
 *                     'antenna*.tracker.actual double', 
 *                     'antenna*.tracker.source string'},
 *                     '06-jan-2005:15','06-jan-2005:16',
 *                     '/data/gcpdaq/arc',
 *                     '/home/gcpdaq/carma_unstable/gcp/array/conf/cal');
 *
 */
#include "gcp/util/ArchiveReader.h"
#include "gcp/util/RegDescription.h"
#include "gcp/util/RegParser.h"

#include "gcp/program/Program.h"

#include "gcp/matlab/MatArchiveConvFn.h"
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

#include "mex.h"
#include "matrix.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>

using namespace std;
using namespace gcp::util;
using namespace gcp::matlab;

struct RegInfo {
  RegDescription desc_;
  DataType::Type outputType_;
  RegAspect aspect_;
  unsigned nEl_;
  unsigned nFrame_;
  std::valarray<unsigned int> cInds_;
  std::valarray<unsigned int> mInds_;
  ArchiveReader::ArchiveTransposeType transpose_;

  RegInfo() {
    outputType_ = DataType::UNKNOWN;
    aspect_     = REG_PLAIN;
    nEl_        = 0;
    nFrame_     = 0;
    transpose_  = ArchiveReader::NONE;
  }

  RegInfo(const RegInfo& info) {
    *this = (RegInfo&) info;
  }

  RegInfo(RegInfo& info) {
    *this = info;
  }

  void operator=(const RegInfo& info) {
    *this = (RegInfo&) info;
  }

  void operator=(RegInfo& info) {
    desc_       = info.desc_;
    outputType_ = info.outputType_;
    aspect_     = info.aspect_;
    nFrame_     = info.nFrame_;
    transpose_  = info.transpose_;

    initialize();
  }

  void initialize() {

    nEl_        = desc_.nEl();

    cInds_.resize(nEl_);
    mInds_.resize(nEl_);
    
    CoordAxes axes = desc_.axes();
    unsigned nEl0  = axes.nEl(0);
    unsigned nEl1  = axes.nAxis()==2 ? axes.nEl(1) : 1;
    unsigned iM, iC;

    if(transpose_ == ArchiveReader::NONE) {

      for(unsigned iEl1=0; iEl1 < nEl1; iEl1++) {
	for(unsigned iEl0=0; iEl0 < nEl0; iEl0++) {
	  iC = iEl0*nEl1 + iEl1;
	  cInds_[iC] = iEl0*nEl1 + iEl1;
	  mInds_[iC] = iEl1*nFrame_*nEl0 + iEl0*nFrame_; // Add iFrame
							 // to get iM
	}
      }

    } else if(transpose_ == ArchiveReader::LAST) {

      for(unsigned iEl1=0; iEl1 < nEl1; iEl1++) {
	for(unsigned iEl0=0; iEl0 < nEl0; iEl0++) {
	  iC = iEl0*nEl1 + iEl1;
	  cInds_[iC] = iEl0*nEl1 + iEl1;
	  mInds_[iC] = iEl1*nEl0 + iEl0;                // Add
							// iFrame*nEl0*nEl1
							// to get iM
	}

      }

    } else if(transpose_ == ArchiveReader::FIRST) {

      for(unsigned iEl1=0; iEl1 < nEl1; iEl1++) {
	for(unsigned iEl0=0; iEl0 < nEl0; iEl0++) {
	  iC = iEl0*nEl1 + iEl1;
	  cInds_[iC] = iEl0*nEl1 + iEl1;
	  mInds_[iC] = iEl0*nFrame_*nEl1 + iEl1;        // Add
							// iFrame*nEl1
							// to get iM
	}

      }
    }
  }
};

// Each time the register changes during a read, we must call this
// function to update register information

std::vector<RegInfo> cacheRegInfo(ArchiveReader& reader, unsigned nFrame);

// Set pointers to matlab arrays as internal memory to the
// ArchiveReader

void setMemory(ArchiveReader& reader, std::vector<MexHandler::MxArray>& vArray, 
	       std::vector<RegInfo>& regInfo);

/**.......................................................................
 * Create the hierarchical board/register structure to be returned to
 * the matlab environment
 */
void createOutputValueStructure(mxArray* plhs[], const mxArray* prhs[],  
				std::vector<RegInfo>& regs,
				unsigned nFrame,
				std::vector<MexHandler::MxArray>& vArray,
				bool rename);

/**
 * Create empty matlab arrays to be filled in by the monitor stream
 */
mxArray* createRegValArray(mxArray* ptr, 
			   RegDescription& reg,
			   DataType::Type formatType,
			   unsigned nFrame,
			   ArchiveReader::ArchiveTransposeType transpose,
			   bool rename,
			   std::string name);

/**.......................................................................
 * Dispatch function for a string read from the archive
 */
static STRING_DISPATCH_FN(dispatchMexString)
{
  mxArray* array = (mxArray*) args;
  mxSetCell(array, ind, mxCreateString((const char*)carr));
}

/**.......................................................................
 * Cache register information that won't change throughout the read
 */
std::vector<RegInfo> cacheRegInfo(ArchiveReader& reader, unsigned nFrame)
{
  std::vector<RegInfo> regInfo;
  RegInfo info;

  for(unsigned iReg=0; iReg < reader.regs_.size(); iReg++) {
    ArchiveReader::ArchiveRegister& reg = reader.regs_[iReg];

    info.desc_       = reader.regDescs_[iReg];
    info.outputType_ = reg.outputType_;
    info.aspect_     = reg.aspect_;
    info.nFrame_     = nFrame;
    info.transpose_  = reg.transpose_;

    info.initialize();

    regInfo.push_back(info);
  }

  return regInfo;
}

/**.......................................................................
 * Create the hierarchical board/register structure to be returned to
 * the matlab environment
 */
void createOutputValueStructure(mxArray* plhs[], const mxArray* prhs[],  
				std::vector<RegInfo>& regs,
				unsigned nFrame,
				std::vector<MexHandler::MxArray>& vArray,
				bool rename)
{
  int dims[2]={1,1};
  
  // Create an empty arraymap structure 
  
  plhs[0] = mxCreateStructArray(2, dims, 0, NULL);
  
  // Get the number of registers requested?
  
  unsigned nreg = regs.size();

  // For each register, add a register field

  std::ostringstream os;

  for(unsigned iReg=0; iReg < regs.size(); iReg++) {
    mxArray* boardPtr = 0;

    if(!rename) {
      boardPtr = MexHandler::addRegisterField(plhs[0], regs[iReg].desc_);
    } else {
      os.str("");
      os << "reg" << iReg;
      boardPtr = MexHandler::addHierNamedStructField(plhs[0], os.str());
    }

    vArray[iReg] = createRegValArray(boardPtr, regs[iReg].desc_, 
				     regs[iReg].outputType_, nFrame, 
				     regs[iReg].transpose_, rename, os.str());
  }
}

/**.......................................................................
 * Create empty matlab arrays to be filled in by the monitor stream
 */
mxArray* createRegValArray(mxArray* ptr, 
			   RegDescription& reg,
			   DataType::Type formatType,
			   unsigned nFrame,
			   ArchiveReader::ArchiveTransposeType transpose,
			   bool rename,
			   std::string name)
{
  int dims[4] = {0,0,0,0};
  unsigned nDim;
  CoordAxes axes = reg.axes();

  // Create and return the array in which values will be stored. We
  // want to format the array to have the same number of dimensions as
  // the C array.

  // If this is a string array, the last dimension is the string
  // length, and we don't count it in the dimensions of the output
  // array

  if(formatType == DataType::STRING) {

    nDim = 1 + axes.nAxis()-1;
    dims[0] = nFrame;

    for(unsigned iDim=1; iDim < axes.nAxis()-1; iDim++) {
      dims[iDim] = axes.nEl(iDim);
    }
    
  } else {

    // Else if not transposing the array, size the output array to be
    // of size nFrame x dimensions of the register

    if(transpose == ArchiveReader::NONE) {

      nDim = axes.nAxis()+1;
      dims[0] = nFrame;

      for(unsigned iDim=0; iDim < axes.nAxis(); iDim++) {
	dims[iDim+1] = axes.nEl(iDim);
      }

      // Else if we are tansposing LAST, the last dimension will be nFrame
      // x dimension of the last axis

    } else if(transpose == ArchiveReader::LAST) {      

      nDim = 1 + axes.nAxis()-1;
      for(unsigned iDim=0; iDim < axes.nAxis()-1; iDim++) {
	dims[iDim] = axes.nEl(iDim);
      }

      dims[axes.nAxis()-1] = nFrame * axes.nEl(axes.nAxis()-1);

      // Else if we are tansposing FIRST, the first dimension will be nFrame
      // x dimension of the last axis

    } else if(transpose == ArchiveReader::FIRST) {      

      nDim = 1 + axes.nAxis()-1;
      for(unsigned iDim=0; iDim < axes.nAxis()-1; iDim++) {
	dims[iDim+1] = axes.nEl(iDim);
      }

      dims[0] = nFrame * axes.nEl(axes.nAxis()-1);
    }

  }

  mxArray* array = MexHandler::createMatlabArray(nDim, dims, formatType);
  
  if(!rename) {
    mxSetField(ptr, 0, reg.aspect()==REG_PLAIN ? reg.blockName().c_str() :
	       reg.aspectName().c_str(), array);
  } else {
    mxSetField(ptr, 0, name.c_str(), array);
  }
  
  return array;
}

/**.......................................................................
 * Set the matlab pointer for each register as the memory into which the
 * reader will write its calibrated register values
 */
void setMemory(ArchiveReader& reader, 
	       std::vector<MexHandler::MxArray>& vArray, 
	       std::vector<RegInfo>& regInfo)
{
  for(unsigned i=0; i < regInfo.size(); i++) {
    ArchiveReader::ArchiveRegister& reg = reader.regs_[i];

    reg.setInputIndices(regInfo[i].cInds_);
    reg.setOutputIndices(regInfo[i].mInds_);
    reg.setExternalMemory(vArray[i].vPtr_, mxGetImagData(vArray[i].array_), vArray[i].array_);
  }
}

void gcp::program::Program::initializeUsage() {};

int gcp::program::Program::main() {
  return 0;
}

/**.......................................................................
 * Entry point from the matlab environment
 */
void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  std::ostringstream usage;

  usage << "Usage: d = " 
	<< "gcpMatReadArcOpt({'receiver.bolometers.adc'}, '24-jan-2009:01:00:00', "
	<< "'24-jan-2009:02:00:00', '/data/sptdaq/arc', "
	<< "'/home/sptdaq/gcp/array/conf/cal')";

  // Reassign output functions

#if 1
  gcp::util::Logger::installStdoutPrintFn(&MexHandler::stdoutPrintFn);
  gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);

  gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
  gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
  gcp::util::ErrHandler::installLogFn(MexHandler::logFn);
#endif

  // Set the string dispatch function appropriate for this environment

  ArchiveConvFn::setStringDispatchFn(dispatchMexString);

  // Check input/output arguments

  if(nrhs < 5) {
    ThrowError(usage.str());
  }

  //-----------------------------------------------------------------------
  // Intialize the archive reader
  //-----------------------------------------------------------------------

  std::string directory;  // The directory in which to look for files 
  std::string calfile;    // The name of the calibration file 
  std::string start_date; // Start time as date:time string 
  std::string end_date;   // End time as date:time string 

  // Get start/end date:time strings 
  
  MexParser startDate(prhs[1]);
  start_date = startDate.getString();

  MexParser endDate(prhs[2]);
  end_date = endDate.getString();

  MexParser dir(prhs[3]);
  directory  = dir.getString();
  
  MexParser cal(prhs[4]);
  calfile    = cal.getString();

  bool rename = false;
  if(nrhs==6) {
    MexParser parser(prhs[5]);
    rename = *parser.getLogicalData();
  }

  ArchiveReader reader(directory, calfile, start_date, end_date);
  reader.getFileList();

  //-----------------------------------------------------------------------
  // Add default registers, plus any registers the user requested
  //-----------------------------------------------------------------------
  
  // And add the user-requested regs

  unsigned nreg = mxGetNumberOfElements(prhs[0]);

  for(unsigned iReg=0; iReg < nreg; iReg++) 
    reader.addRegister(mxArrayToString(mxGetCell(prhs[0],iReg)));

  //-----------------------------------------------------------------------
  // Force reading of the first array map, just to check the register
  // selection
  //-----------------------------------------------------------------------

  reader.readFirstArrayMap();

  //-----------------------------------------------------------------------
  // Count the number of frames
  //-----------------------------------------------------------------------

  unsigned nFrame = reader.countFrames();
  std::vector<RegInfo> regInfo = cacheRegInfo(reader, nFrame);

  //-----------------------------------------------------------------------
  // Create the hierarchical matlab structure corresponding to the
  // arraymap
  //-----------------------------------------------------------------------

  std::vector<MexHandler::MxArray> vArray(regInfo.size());
  createOutputValueStructure(plhs, prhs, regInfo, nFrame, vArray, rename);

  //-----------------------------------------------------------------------
  // Loop, reading frames
  //-----------------------------------------------------------------------

  reader.resetToBeginning();

  //-----------------------------------------------------------------------
  // Now set our newly-created memory as memory into which the reader
  // will copy its calibrated values
  //-----------------------------------------------------------------------

  setMemory(reader, vArray, regInfo);

  while(reader.readNextFrame()) {
    reader.readRegs();
  }

  std::cout << "\rReading..." << std::fixed << std::setprecision(0) << "100%";
  fflush(stdout);

  COUT("");
}

