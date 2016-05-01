/**.......................................................................
 * MATLAB Mex file for reading from GCP archive data files.
 *
 * Use like:
 *
 * d=gcpMatFastReadArc({'array.frame.record','corr.band0.usb[0][0]',
 *                     'antenna*.tracker.actual double', 'antenna*.tracker.source string'},
 *                     '06-jan-2005:15','06-jan-2005:16',
 *                     '/data/gcpdaq/arc','/home/gcpdaq/carma_unstable/gcp/array/conf/cal');
 *
 */
#include "gcp/array/code/unix/libmonitor_src/monitor_stream.h"

#include "gcp/util/Debug.h"
#include "gcp/util/Monitor.h"
#include "gcp/util/RegDescription.h"
#include "gcp/util/RegParser.h"

#include "gcp/matlab/MexHandler.h"

#include "mex.h"
#include "matrix.h"

#include <iostream.h>
#include <math.h>

using namespace std;
using namespace gcp::util;
using namespace gcp::matlab;

// Macro to assign values from a real C array (double* dPtr) to a
// matlab array (MexHandler::MxArray* array).  Matlab arrays index
// backwards from C arrays; i.e., when encoding a multi-dimensional
// array into 1-dimension, the index which changes most rapidly in C
// is the slowest changing index in matlab

#define ASSIGN_VALS(cast) \
  {\
    cast* ptr = (cast*)array.vPtr_;\
    switch(axes.nAxis()) {\
    case 2:\
      {\
        unsigned nEl0 = axes.nEl(0);\
        unsigned nEl1 = axes.nEl(1);\
        unsigned iM, iC;\
        for(unsigned iEl1=0; iEl1 < nEl1; iEl1++)\
          for(unsigned iEl0=0; iEl0 < nEl0; iEl0++) {\
            iC = iEl0*nEl1 + iEl1;\
            iM = iEl1*nFrame*nEl0 + iEl0*nFrame + iFrame;\
            *(ptr + iM) = (cast)*(dPtr + iC);\
          }\
      }\
      break;\
    case 1:\
      {\
        unsigned nEl0 = axes.nEl(0);\
        unsigned iM, iC;\
        for(unsigned iEl0=0; iEl0 < nEl0; iEl0++) {\
          iC = iEl0;\
          iM = iEl0 * nFrame + iFrame;\
          *(ptr + iM) = (cast)*(dPtr + iC);\
        }\
      }\
      break;\
    default:\
      break;\
   }\
  }

// Macro to assign values from a complex C array (double* dPtr) to a
// matlab array (MexHandler::MxArray* array).  Matlab arrays index
// backwards from C arrays; i.e., when encoding a multi-dimensional
// array into 1-dimension, the index which changes most rapidly in C
// is the slowest changing inde in matlab

#define ASSIGN_COMPLEX_VALS(cast) \
{\
  static Complex<float> cmplx;\
  unsigned nEl0, nEl1;\
  unsigned iM, iC;\
  cast* fPtr     = (cast*)array.vPtr_;\
  cast* realPtr  = (cast*)mxGetData(array.array_);\
  cast* imagPtr  = (cast*)mxGetImagData(array.array_);\
  switch(axes.nAxis()) {\
  case 2:\
    nEl0 = axes.nEl(0);\
    nEl1 = axes.nEl(1);\
    break;\
  case 1:\
    nEl0 = axes.nEl(0);\
    nEl1 = 1;\
    break;\
  };\
  for(unsigned iEl1=0; iEl1 < nEl1; iEl1++)\
    for(unsigned iEl0=0; iEl0 < nEl0; iEl0++) {\
      iC = iEl0*nEl1 + iEl1;\
      iM = iEl1*nFrame*nEl0 + iEl0*nFrame + iFrame;\
      cmplx = *((Complex<cast>::Data*)(dPtr + iC));\
      switch (reg.aspect()) {\
      case REG_PLAIN:\
	*(realPtr + iM) = (cast)cmplx.real();\
	*(imagPtr + iM) = (cast)cmplx.imag();\
	break;\
      case REG_REAL:\
	*(fPtr + iM) = (cast)cmplx.real();\
	break;\
      case REG_IMAG:\
	*(fPtr + iM) = (cast)cmplx.imag();\
	break;\
      case REG_AMP:\
	*(fPtr + iM) = (cast)cmplx.amp();\
	break;\
      case REG_PHASE:\
	*(fPtr + iM) = (cast)cmplx.phaseInDegrees();\
	break;\
      }\
    }\
  } 

unsigned countFrames(const mxArray* prhs[], 
		     std::vector<gcp::util::RegDescription>& regs,
		     std::vector<gcp::util::MonitorDataType::FormatType>& formats);

void readFrames(const mxArray* prhs[], std::vector<MexHandler::MxArray>& vArray, unsigned nFrame);

/**
 * Each time the register changes during a read, we must call this
 * function to update register information
 */
void cacheRegInfo(Monitor* monitor, std::vector<RegDescription>& regs, 
		  std::vector<int>& startSlots, 
		  std::vector<gcp::util::MonitorDataType::FormatType>& nativeFormats, 
		  std::vector<gcp::util::MonitorDataType::FormatType>& selectedFormats, 
		  double** slotPtr);

void assignOutputValues(mxArray* plhs[], const mxArray* prhs[],  
			std::vector<gcp::util::RegDescription>& regs,
			std::vector<std::vector<std::vector<MonitorDataType> > >& regVals);

/**.......................................................................
 * Create the hierarchical board/register structure to be returned to
 * the matlab environment
 */
void createOutputValueStructure(mxArray* plhs[], const mxArray* prhs[],  
				std::vector<gcp::util::RegDescription>& regs,
				std::vector<gcp::util::MonitorDataType::FormatType>& formats,
				unsigned nFrame,
				std::vector<MexHandler::MxArray>& vArray);

/**
 * Create empty matlab arrays to be filled in by the monitor stream
 */
mxArray* createRegValArray(mxArray* ptr, 
			   RegDescription& reg,
			   MonitorDataType::FormatType formatType,
			   unsigned nFrame);

/**
 * Assign values returned from the monitor stream into the
 * corresponding matlab array
 */
void assignVals(MexHandler::MxArray& array, 
		RegDescription& reg,
		MonitorDataType::FormatType nativeFormat,
		MonitorDataType::FormatType selectedFormat,
		double* dPtr,
		unsigned iFrame,
		unsigned nFrame);

void assignStringVals(mxArray* array, double* dPtr, unsigned len, unsigned iFrame);
void assignDateVals(gcp::matlab::MexHandler::MxArray& array, double* dPtr, 
		    gcp::util::MonitorDataType::FormatType format, unsigned iFrame);

/**.......................................................................
 * Entry point from the matlab environment
 */
void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  int nreg;           // Number of registers being read 
  int i;

  try {  

    gcp::util::Logger::installStdoutPrintFn(MexHandler::stdoutPrintFn);
    gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);

    gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
    gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
    gcp::util::ErrHandler::installLogFn(MexHandler::logFn);

    // Check input/output arguments

    MexHandler::checkArgs(nlhs, nrhs, 1, 5);

    // Read register values now
    
    std::vector<gcp::util::RegDescription> regs;
    std::vector<gcp::util::MonitorDataType::FormatType> formats;

    // Count how many frames are in the requested monitor stream

    unsigned nFrame = countFrames(prhs, regs, formats);

    // Create a matlab array the same size as the number of registers

    std::vector<MexHandler::MxArray> vArray(regs.size());

    // Create the hierarchical matlab structure corresponding to the arraymap

    createOutputValueStructure(plhs, prhs, regs, formats, nFrame, vArray);

    // Read frames

    readFrames(prhs, vArray, nFrame);

  } catch(gcp::util::Exception& err) {
    mexErrMsgTxt(err.what());
  } catch(...) {
    return;
  }
}

/**.......................................................................
 * Count the number of frames in the monitor stream
 */
unsigned countFrames(const mxArray* prhs[], 
		     std::vector<gcp::util::RegDescription>& regs,
		     std::vector<gcp::util::MonitorDataType::FormatType>& formats)

{
  char* directory;    // The directory in which to look for files 
  char* calfile;      // The name of the calibration file 
  char* start_date;   // Start time as date:time string 
  char* end_date;     // End time as date:time string 

  // How many registers requested?
  
  unsigned nreg = mxGetNumberOfElements(prhs[0]);
  
  // Get start/end date:time strings 
  
  start_date = mxArrayToString(prhs[1]);
  end_date   = mxArrayToString(prhs[2]);
  
  // Get arcdir and calfile locations 
  
  directory  = mxArrayToString(prhs[3]);
  calfile    = mxArrayToString(prhs[4]);
  
  Monitor monitor(directory, calfile, start_date, end_date);

  // Always add the UTC register
  
  monitor.addRegister("array.frame.utc double");
  
  // And add the user-requested regs

  for(unsigned iReg=0; iReg < nreg; iReg++) {

    try {
      monitor.addRegister(mxArrayToString(mxGetCell(prhs[0],iReg)));
    } catch(Exception& err) {
      COUT(err.what());
    } catch(...) {
    }

  }
  
  // Return the vector of register descriptions
  
  unsigned nFrame = monitor.countFrames();

  regs    = monitor.selectedRegs();
  formats = monitor.selectedFormats();
 
  return nFrame;
}

/**.......................................................................
 * Read frames from the monitor stream
 */
void readFrames(const mxArray* prhs[], std::vector<MexHandler::MxArray>& vArray, unsigned nFrame) 
{
  char* directory;    // The directory in which to look for files 
  char* calfile;      // The name of the calibration file 
  char* start_date;   // Start time as date:time string 
  char* end_date;     // End time as date:time string 

  // How many registers requested?
  
  unsigned nreg = mxGetNumberOfElements(prhs[0]);
  
  // Get start/end date:time strings 
  
  start_date = mxArrayToString(prhs[1]);
  end_date   = mxArrayToString(prhs[2]);
  
  // Get arcdir and calfile locations 
  
  directory  = mxArrayToString(prhs[3]);
  calfile    = mxArrayToString(prhs[4]);
  
  Monitor monitor(directory, calfile, start_date, end_date);

  // Always add the UTC register

  monitor.addRegister("array.frame.utc double");
  
  // And add the user-requested regs
  
  for(unsigned iReg=0; iReg < nreg; iReg++) {

    try {
      monitor.addRegister(mxArrayToString(mxGetCell(prhs[0],iReg)));
    } catch(...) {
    }

  }
  
  monitor.reinitialize();

  std::vector<gcp::util::RegDescription> regs;
  std::vector<gcp::util::MonitorDataType::FormatType> nativeFormats;
  std::vector<gcp::util::MonitorDataType::FormatType> selectedFormats;
  std::vector<int> startSlots;
  double* slotPtr=0;
  unsigned iFrame=0;

  gcp::array::MsReadState state;

  cacheRegInfo(&monitor, regs, startSlots, nativeFormats, selectedFormats, 
	       &slotPtr);

  // Loop, reading frames

  while((state=monitor.readNextFrame()) != gcp::array::MS_READ_ENDED) {

    // If the register map changed, re-cache the register information

    if(state==gcp::array::MS_READ_REGMAP) {
      cacheRegInfo(&monitor, regs, startSlots, nativeFormats, selectedFormats, &slotPtr);
    }

    // Loop over all requested registers from this frame

    if(state == gcp::array::MS_READ_DONE) {
      for(unsigned iReg=0; iReg < regs.size(); iReg++) {
	assignVals(vArray[iReg], regs[iReg], nativeFormats[iReg], selectedFormats[iReg],
		   (slotPtr+startSlots[iReg]), iFrame, nFrame);
      }
      ++iFrame;
    }
  }
}

/**.......................................................................
 * Each time the register changes during a read, we must call this
 * function to update register information
 */
void cacheRegInfo(Monitor* monitor, std::vector<RegDescription>& regs, 
		  std::vector<int>& startSlots, 
		  std::vector<gcp::util::MonitorDataType::FormatType>& nativeFormats, 
		  std::vector<gcp::util::MonitorDataType::FormatType>& selectedFormats, 
		  double** slotPtr)
{
  // Cache the register information
  
  regs            = monitor->selectedRegs();
  nativeFormats   = monitor->nativeFormats();
  selectedFormats = monitor->selectedFormats();

  startSlots.resize(regs.size());

  for(unsigned iReg=0; iReg < regs.size(); iReg++) 
    startSlots[iReg] = regs[iReg].startSlot();

  *slotPtr = monitor->getCalSlotPtr();
}

/**.......................................................................
 * Create the hierarchical board/register structure to be returned to
 * the matlab environment
 */
void createOutputValueStructure(mxArray* plhs[], const mxArray* prhs[],  
				std::vector<gcp::util::RegDescription>& regs,
				std::vector<gcp::util::MonitorDataType::FormatType>& formats,
				unsigned nFrame,
				std::vector<MexHandler::MxArray>& vArray)
{
  int dims[2]={1,1};
  
  // Create an empty arraymap structure 
  
  plhs[0] = mxCreateStructArray(2, dims, 0, NULL);
  
  // Get the number of registers requested?
  
  unsigned nreg = mxGetNumberOfElements(prhs[0]);

  // For each register, add a register field

  for(unsigned iReg=0; iReg < regs.size(); iReg++) {
    mxArray* boardPtr = MexHandler::addRegisterField(plhs[0], regs[iReg]);
    vArray[iReg]      = createRegValArray(boardPtr, regs[iReg], formats[iReg], nFrame);
  }
}

/**.......................................................................
 * Create empty matlab arrays to be filled in by the monitor stream
 */
mxArray* createRegValArray(mxArray* ptr, 
			   RegDescription& reg,
			   MonitorDataType::FormatType formatType,
			   unsigned nFrame)
{
  int dims[4];
  unsigned nDim;
  CoordAxes axes = reg.axes();

  // Create and return the array in which values will be stored. We
  // want to format the array to have the same number of dimensions as
  // the C array.

  dims[0] = nFrame;
  if(formatType != MonitorDataType::FM_STRING) {
    nDim = axes.nAxis();
    for(unsigned iDim=0; iDim < nDim; iDim++)
      dims[iDim+1] = axes.nEl(iDim);

  } else {
    dims[1] = 1;
    nDim = 1;
  }

  mxArray* array = MexHandler::createMatlabArray(nDim+1, dims, formatType);
  
  mxSetField(ptr, 0, reg.aspect()==REG_PLAIN ? reg.blockName().c_str() :
  	     reg.aspectName().c_str(), array);
  
  return array;
}

/**.......................................................................
 * Assign values returned from the monitor stream into the
 * corresponding matlab array
 */
void assignVals(MexHandler::MxArray& array, 
		RegDescription& reg,
		MonitorDataType::FormatType nativeFormat,
		MonitorDataType::FormatType selectedFormat,
		double* dPtr,
		unsigned iFrame,
		unsigned nFrame)
{
  CoordAxes& axes = reg.axes_;

  if(selectedFormat == MonitorDataType::FM_STRING) {
    assignStringVals(array.array_, dPtr, reg.nEl(), iFrame);
  } else if(nativeFormat == MonitorDataType::FM_DATE) {
    assignDateVals(array, dPtr, selectedFormat, iFrame);
  } else {
 
   switch (selectedFormat) {

    case MonitorDataType::FM_BOOL:
    case MonitorDataType::FM_UCHAR:
      ASSIGN_VALS(unsigned char);
      break;
    case MonitorDataType::FM_CHAR:
      ASSIGN_VALS(char);
      break;
    case MonitorDataType::FM_USHORT:
      ASSIGN_VALS(unsigned short);
      break;
    case MonitorDataType::FM_SHORT:
      ASSIGN_VALS(short);
      break;
    case MonitorDataType::FM_UINT:
      ASSIGN_VALS(unsigned int);
      break;
    case MonitorDataType::FM_INT:
      ASSIGN_VALS(int);
      break;
    case MonitorDataType::FM_ULONG:
      ASSIGN_VALS(unsigned long);
      break;
    case MonitorDataType::FM_LONG:
      ASSIGN_VALS(long);
      break;
    case MonitorDataType::FM_FLOAT:
      ASSIGN_VALS(float);
      break;
    case MonitorDataType::FM_COMPLEX_FLOAT:
      ASSIGN_COMPLEX_VALS(float);
      break;
    case MonitorDataType::FM_DOUBLE:
      ASSIGN_VALS(double);
      break;
    }
  }
}

/**.......................................................................
 * Set a string value in a matlab cell array
 */
void assignStringVals(mxArray* array, double* dPtr, unsigned len, unsigned iFrame)
{
  static std::ostringstream os;

  os.str("");
  for(unsigned i=0; i < len; i++)
    os << (char)*(dPtr+i);
  
  mxSetCell(array, iFrame, mxCreateString(os.str().c_str()));
}

/**.......................................................................
 * Convert a date to a string value in a matlab cell array
 */
void assignDateVals(MexHandler::MxArray& array, double* dPtr, MonitorDataType::FormatType format, unsigned iFrame)
{
  RegDate date(*((RegDate::Data*)dPtr));
  static int count=0;

  ++count;
  if(count % 100 == 0)
    COUT("Reading: " << date);

  if(format == MonitorDataType::FM_DATE) 
    mxSetCell(array.array_, iFrame, mxCreateString(date.str().c_str()));
  else 
    *((double*)array.vPtr_ + iFrame) = date.mjd();
}
