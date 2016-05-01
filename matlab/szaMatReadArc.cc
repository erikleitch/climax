/**.......................................................................
 * MATLAB Mex file for reading from GCP archive data files.
 *
 * Use like:
 *
 * d=gcpMatReadArc({'array.frame.record','corr.band0.usb[0][0]',
 *                  'antenna*.tracker.actual double', 'antenna*.tracker.source string'},
 *                  '06-jan-2005:15','06-jan-2005:16',
 *                  '/data/gcpdaq/arc','/home/gcpdaq/carma_unstable/gcp/array/conf/cal');
 *
 * >> d
 * d = 
 * 
 *        array: [1x1 struct]
 *         corr: [1x1 struct]
 *     antenna0: [1x1 struct]
 *     antenna1: [1x1 struct]
 *     antenna2: [1x1 struct]
 *     antenna3: [1x1 struct]
 *     antenna4: [1x1 struct]
 *     antenna5: [1x1 struct]
 *     antenna6: [1x1 struct]
 *     antenna7: [1x1 struct]
 * 
 * >> d.antenna0.tracker.source(1)
 * d.antenna0.tracker.source(1)
 * 
 * ans = 
 * 
 *     '(none)'
 * 
 * >> d.array.frame.record(1:5)
 * d.array.frame.record(1:5)
 * 
 * ans =
 * 
 *       289539
 *       289559
 *       289579
 *       289599
 *       289619
 * 
 * >> d.antenna0.tracker.actual(1,:)
 * d.antenna0.tracker.actual(1,:)
 * 
 * ans =
 * 
 *   -97.4541   10.8045         0
 * 
 * >> d.antenna0.tracker.actual(1:5,1)
 * d.antenna0.tracker.actual(1:5,1)
 * 
 * ans =
 * 
 *   -97.4541
 *   -97.4541
 *   -97.4541
 *   -97.4541
 *   -97.4541
 */

#include "gcp/util/Monitor.h"

#include <iostream.h>
#include <math.h>

#include "gcp/util/Debug.h"
#include "gcp/util/RegDescription.h"
#include "gcp/util/RegParser.h"

#include "mex.h"
#include "matrix.h"

#define ASSIGN_VALS(cast, member) \
  {\
    cast* ptr = (cast*)vPtr;\
    switch(axes.nAxis()) {\
    case 2:\
      {\
        unsigned nEl0 = axes.nEl(0);\
        unsigned nEl1 = axes.nEl(1);\
        unsigned iEl;\
        unsigned i=0;\
        for(unsigned iEl1=0; iEl1 < nEl1; iEl1++)\
          for(unsigned iEl0=0; iEl0 < nEl0; iEl0++)\
            for(unsigned iFrame=0; iFrame < nFrame; iFrame++, i++) {\
              iEl = iEl0 * nEl1 + iEl1;\
              *(ptr + i) = (cast)regVals[iFrame][iReg][iEl].val.data_.member;\
            }\
      }\
      break;\
    case 1:\
      {\
        unsigned nEl0 = axes.nEl(0);\
        unsigned iEl;\
        unsigned i=0;\
       for(unsigned iEl0=0; iEl0 < nEl0; iEl0++)\
          for(unsigned iFrame=0; iFrame < nFrame; iFrame++, i++) {\
            iEl = iEl0;\
            *(ptr + i) = (cast)regVals[iFrame][iReg][iEl].val.data_.member;\
          }\
      }\
      break;\
    default:\
      break;\
   }\
  }

#define ASSIGN_COMPLEX_VALS(cast) \
  {\
    cast* realPtr  = (cast*)mxGetData(array);\
    cast* imagPtr  = (cast*)mxGetImagData(array);\
    switch(axes.nAxis()) {\
    case 2:\
      {\
        unsigned nEl0 = axes.nEl(0);\
        unsigned nEl1 = axes.nEl(1);\
        unsigned iEl;\
        unsigned i=0;\
        for(unsigned iEl1=0; iEl1 < nEl1; iEl1++)\
          for(unsigned iEl0=0; iEl0 < nEl0; iEl0++)\
            for(unsigned iFrame=0; iFrame < nFrame; iFrame++, i++) {\
              iEl = iEl0 * nEl1 + iEl1;\
              *(realPtr + i) = (cast)regVals[iFrame][iReg][iEl].val.data_.cf.real_;\
              *(imagPtr + i) = (cast)regVals[iFrame][iReg][iEl].val.data_.cf.imag_;\
            }\
      }\
      break;\
    case 1:\
      {\
        unsigned nEl0 = axes.nEl(0);\
        unsigned iEl;\
        unsigned i=0;\
        for(unsigned iEl0=0; iEl0 < nEl0; iEl0++)\
          for(unsigned iFrame=0; iFrame < nFrame; iFrame++, i++) {\
            iEl = iEl0;\
            *(realPtr + i) = (cast)regVals[iFrame][iReg][iEl].val.data_.cf.real_;\
            *(imagPtr + i) = (cast)regVals[iFrame][iReg][iEl].val.data_.cf.imag_;\
          }\
      }\
      break;\
    default:\
      break;\
   }\
  }

using namespace std;
using namespace gcp::util;

void checkArgs(int nlhs, int nrhs);

void readRegs(const mxArray* prhs[], 
	      std::vector<gcp::util::RegDescription>& regs,
	      std::vector<std::vector<std::vector<MonitorDataType> > >& regVals);

void assignOutputValues(mxArray* plhs[], const mxArray* prhs[],  
			std::vector<gcp::util::RegDescription>& regs,
			std::vector<std::vector<std::vector<MonitorDataType> > >& regVals);

mxArray* addRegisterField(mxArray* ptr, RegDescription& reg);
mxArray* addNamedStructField(mxArray* parentPtr, std::string fieldName);
mxArray* addNamedCellField(mxArray* parentPtr, std::string fieldName, unsigned iDim);

mxClassID matlabTypeOf(MonitorDataType& dataType);
mxClassID matlabTypeOf(MonitorDataType::FormatType formatType);

mxComplexity matlabComplexityOf(MonitorDataType& dataType);
mxComplexity matlabComplexityOf(MonitorDataType::FormatType formatType);

void assignRegVals(mxArray* ptr, 
		   RegDescription& reg,
		   unsigned iReg, 
		   std::vector<std::vector<std::vector<MonitorDataType> > >& regVals);

void* createRegValArray(mxArray* ptr, 
			RegDescription& reg,
			MonitorDataType::FormatType formatType,
			unsigned nFrame);

void assignStringVals(mxArray* array,
		      unsigned iReg,
		      std::vector<std::vector<std::vector<MonitorDataType> > >& regVals);

void assignDateVals(mxArray* array,
		     unsigned iReg,
		      std::vector<std::vector<std::vector<MonitorDataType> > >& regVals);

mxArray* createMatlabArray(int ndim, int* dims, MonitorDataType& dataType);
mxArray* createMatlabArray(int ndim, int* dims, MonitorDataType::FormatType formatType);

void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  int nreg;           // Number of registers being read 
  int i;

  try {  
  // The Array of register values returned by the Monitor object
  
  std::vector<std::vector<std::vector<MonitorDataType> > > regVals;
  
  // Read register values now
  
  std::vector<gcp::util::RegDescription> regs;

  readRegs(prhs, regs, regVals);
  
  // Now assign output values
  
  unsigned nFrame = regVals.size();

  if(nFrame == 0)
    return;
  
  unsigned nReg = regVals[0].size();
  
  if(nReg == 0)
    return;

  assignOutputValues(plhs, prhs, regs, regVals);
  } catch(...) {}
}

/**.......................................................................
 * Add a register block to the struct
 */
mxArray* addRegisterField(mxArray* ptr, RegDescription& reg)
{
  mxArray* regMapPtr = addNamedStructField(ptr,       reg.regMapName());
  mxArray* boardPtr  = addNamedStructField(regMapPtr, reg.boardName());
  mxArray* regPtr    = addNamedStructField(boardPtr,  reg.blockName());

  if(reg.aspect() != REG_PLAIN) {
    (void)addNamedStructField(regPtr,  reg.aspectName());
    return regPtr;
  } else {
    return boardPtr;
  }
}

/**.......................................................................
 * Add a named struct field to the passed structutr
 */
mxArray* addNamedStructField(mxArray* parentPtr, std::string fieldName)
{
  int dims[2]={1,1};
  
  // Now add a field for the board
  
  mxArray* childPtr=NULL;
  
  if((childPtr=mxGetField(parentPtr, 0, fieldName.c_str()))==NULL) {
    
    // Add it
    
    mxAddField(parentPtr, fieldName.c_str());
    
    // create empty register substructure
    
    childPtr = mxCreateStructArray(2, dims, 0, NULL);
    
    // link it to board structure
    
    mxSetField(parentPtr, 0, fieldName.c_str(), childPtr);
  }
  
  return childPtr;
}

/**.......................................................................
 * Add a named struct field to the passed structutr
 */
mxArray* addNamedCellField(mxArray* parentPtr, std::string fieldName, unsigned iDim)
{
  int dims[2]={1, iDim};
  
  // Now add a field for the board
  
  mxArray* childPtr=NULL;
  
  if((childPtr=mxGetField(parentPtr, 0, fieldName.c_str()))==NULL) {
    
    // Add it
    
    mxAddField(parentPtr, fieldName.c_str());
    
    // create empty register substructure
    
    childPtr = mxCreateCellArray(2, dims);
    
    // link it to board structure
    
    mxSetField(parentPtr, 0, fieldName.c_str(), childPtr);
  }
  
  return childPtr;
}

void checkArgs(int nlhs, int nrhs)
{
  // Get the command-line arguments from matlab 
  
  if(nrhs != 5)
    mexErrMsgTxt("Must have 5 input arguments");
  
  if(nlhs != 1)
    mexErrMsgTxt("Must have 1 output arguments");
} 

void readRegs(const mxArray* prhs[], 
	      std::vector<gcp::util::RegDescription>& regs,
	      std::vector<std::vector<std::vector<MonitorDataType> > >& regVals)
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
  
  Debug::setLevel(0);
  
  try {
    
    Monitor monitor(directory, calfile, start_date, end_date);

    // Always add the UTC register

    monitor.addRegister("array.frame.utc double");

    // And add the user-requested regs
    
    for(unsigned iReg=0; iReg < nreg; iReg++) {
      monitor.addRegister(mxArrayToString(mxGetCell(prhs[0],iReg)));
    }

    // Finally, run this object    

    regVals = monitor.readRegsAsDataTypes();

    // Return the vector of register descriptions
    
    regs = monitor.selectedRegs();

  } catch(Exception& err) {
    mexPrintf("%s\n", err.what());
    return;
  }
}

void assignOutputValues(mxArray* plhs[], const mxArray* prhs[],  
			std::vector<gcp::util::RegDescription>& regs,
			std::vector<std::vector<std::vector<MonitorDataType> > >& regVals)
{
  int dims[2]={1,1};
  unsigned nFrame = regVals.size();
  
  // Create empty arraymap structure 
  
  plhs[0] = mxCreateStructArray(2, dims, 0, NULL);
  
  // How many registers requested?
  
  unsigned nreg = mxGetNumberOfElements(prhs[0]);
  
  for(unsigned iReg=0; iReg < regs.size(); iReg++) {
    mxArray* boardPtr = addRegisterField(plhs[0], regs[iReg]);

    assignRegVals(boardPtr, regs[iReg], iReg, regVals);
  }
}

/**.......................................................................
 * Create the hierarchical board/register structure to be returned to
 * the matlab environment
 */
void createOutputValueStructure(mxArray* plhs[], const mxArray* prhs[],  
				std::vector<gcp::util::RegDescription>& regs,
				std::vector<gcp::util::MonitorDataType::FormatType>& formats,
				unsigned nFrame,
				std::vector<void*>& registerData)
{
  int dims[2]={1,1};

  // Create empty arraymap structure 
  
  plhs[0] = mxCreateStructArray(2, dims, 0, NULL);
  
  // How many registers requested?
  
  unsigned nreg = mxGetNumberOfElements(prhs[0]);
  
  for(unsigned iReg=0; iReg < regs.size(); iReg++) {
    mxArray* boardPtr = addRegisterField(plhs[0], regs[iReg]);
    registerData[iReg] = createRegValArray(boardPtr, regs[iReg], formats[iReg], 
					   nFrame);
  }
}

/**.......................................................................
 * Set a string value in a matlab cell array
 */
void assignStringVals(mxArray* array,
		     unsigned iReg,
		     std::vector<std::vector<std::vector<MonitorDataType> > >& regVals)
{
  for(unsigned iFrame=0; iFrame < regVals.size(); iFrame++)
    mxSetCell(array, iFrame, mxCreateString(regVals[iFrame][iReg][0].stringVal.c_str()));
}

/**.......................................................................
 * Convert a date to a string value in a matlab cell array
 */
void assignDateVals(mxArray* array,
		    unsigned iReg,
		    std::vector<std::vector<std::vector<MonitorDataType> > >& regVals)
{
  for(unsigned iFrame=0; iFrame < regVals.size(); iFrame++) {
    RegDate date(regVals[iFrame][iReg][0].val.data_.date);
    mxSetCell(array, iFrame, mxCreateString(date.str().c_str()));
  }
}

/**.......................................................................
 * Create empty matlab arrays to be filled in by the monitor stream
 */
void* createRegValArray(mxArray* ptr, 
			RegDescription& reg,
			MonitorDataType::FormatType formatType,
			unsigned nFrame)
{
  int dims[4];
  int nDim = 2; // Initialize to 2

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

  mxArray* array = createMatlabArray(nDim+1, dims, formatType);

  mxSetField(ptr, 0, reg.aspect()==REG_PLAIN ? reg.blockName().c_str() :
  	     reg.aspectName().c_str(), array);

  return mxGetData(array);
}

/**.......................................................................
 * Assign values returned from the monitor stream into the
 * corresponding matlab array
 */
void assignRegVals(mxArray* ptr, 
		   RegDescription& reg,
		   unsigned iReg, 
		   std::vector<std::vector<std::vector<MonitorDataType> > >& regVals)
{
  unsigned nFrame = regVals.size();
  unsigned nEl    = regVals[0][iReg].size();
  int dims[4]     = {nFrame, nEl};
  int nDim = 2; // Initialize to 2

  MonitorDataType& dataType = regVals[0][iReg][0];

  CoordAxes axes = reg.axes();

  // Create and return the array in which values will be stored We
  // want to format the array to have the same number of dimensions as
  // the C array.

  dims[0] = nFrame;
  if(dataType.selectedFormat != MonitorDataType::FM_STRING) {
    nDim = axes.nAxis();
    for(unsigned iDim=0; iDim < nDim; iDim++)
      dims[iDim+1] = axes.nEl(iDim);

  } else {
    dims[1] = 1;
    nDim = 1;
  }

  mxArray* array = createMatlabArray(nDim+1, dims, dataType);
  mxSetField(ptr, 0, reg.aspect()==REG_PLAIN ? reg.blockName().c_str() :
  	     reg.aspectName().c_str(), array);

  void* vPtr = mxGetData(array);

  switch (dataType.selectedFormat) {
  case MonitorDataType::FM_BOOL:
  case MonitorDataType::FM_UCHAR:
    ASSIGN_VALS(unsigned char, uc);
    break;
  case MonitorDataType::FM_CHAR:
    ASSIGN_VALS(char, c);
    break;
  case MonitorDataType::FM_UINT:
    ASSIGN_VALS(unsigned int, ui);
    break;
  case MonitorDataType::FM_INT:
    ASSIGN_VALS(int, i);
    break;
  case MonitorDataType::FM_ULONG:
    ASSIGN_VALS(unsigned long, ul);
    break;
  case MonitorDataType::FM_LONG:
    ASSIGN_VALS(long, l);
    break;
  case MonitorDataType::FM_FLOAT:
    ASSIGN_VALS(float, f);
    break;
  case MonitorDataType::FM_COMPLEX_FLOAT:
    ASSIGN_COMPLEX_VALS(float);
    break;
  case MonitorDataType::FM_DOUBLE:
    ASSIGN_VALS(double, d);
    break;
  case MonitorDataType::FM_STRING:
    assignStringVals(array, iReg, regVals);
    break;
  case MonitorDataType::FM_DATE:
    assignDateVals(array, iReg, regVals);
    break;
  }
}

mxArray* createMatlabArray(int ndim, int* dims, MonitorDataType& dataType)
{
  return createMatlabArray(ndim, dims, dataType.selectedFormat);
}

mxArray* createMatlabArray(int ndim, int* dims, 
			   MonitorDataType::FormatType formatType)
{
  switch (formatType) {
  case MonitorDataType::FM_STRING:
  case MonitorDataType::FM_DATE:
    return mxCreateCellArray(ndim, dims);
  default:
    return mxCreateNumericArray(ndim, dims, 
				matlabTypeOf(formatType), 
				matlabComplexityOf(formatType));
    break;
  }
}


mxClassID matlabTypeOf(MonitorDataType& dataType)
{
  return matlabTypeOf(dataType.selectedFormat);
}

mxClassID matlabTypeOf(MonitorDataType::FormatType formatType)
{
  switch (formatType) {
  case MonitorDataType::FM_BOOL:
  case MonitorDataType::FM_UCHAR:
    return mxUINT8_CLASS;
    break;
  case MonitorDataType::FM_CHAR:
    return mxINT8_CLASS;
    break;
  case MonitorDataType::FM_USHORT:
    return mxUINT16_CLASS;
    break;
  case MonitorDataType::FM_SHORT:
    return mxINT16_CLASS;
    break;
  case MonitorDataType::FM_UINT:
    return mxUINT32_CLASS;
    break;
  case MonitorDataType::FM_INT:
    return mxINT32_CLASS;
    break;
  case MonitorDataType::FM_ULONG:
    return mxUINT64_CLASS;
    break;
  case MonitorDataType::FM_LONG:
    return mxINT64_CLASS;
    break;
  case MonitorDataType::FM_FLOAT:
    return mxSINGLE_CLASS;
    break;
  case MonitorDataType::FM_COMPLEX_FLOAT:
    return mxSINGLE_CLASS;
    break;
  case MonitorDataType::FM_DOUBLE:
    return mxDOUBLE_CLASS;
    break;
  }
}

mxComplexity matlabComplexityOf(MonitorDataType& dataType)
{
  return matlabComplexityOf(dataType.selectedFormat);
}

mxComplexity matlabComplexityOf(MonitorDataType::FormatType formatType)
{
  switch (formatType) {
  case MonitorDataType::FM_COMPLEX_FLOAT:
    return mxCOMPLEX;
    break;
  default:
    return mxREAL;
    break;
  }
}
