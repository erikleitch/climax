/*
 * CLP 050109
 * Start on Matlab mex funtion to read from GCP data archive and return Matlab
 * data structure.
 * Based on:
 * /usr/local/matlab6p1/extern/examples/mex/mexcpp.cpp
 *
 * Makefile needs to do something like:
 
 gcpReadArcc.mexglx: gcpReadArcc.cc
 mex gcpReadArcc.cc -f /usr/local/matlab6p1/bin/cxxopts.sh \
 -I/home/gcpdaq/carma \
 -L../lib -lSzaUtil -lSzaArrayShare -lSzaMonitor -lSzaSla -lSzaSrc \
 -lrt -lreadline -ltermcap
 
 * Currently just provokes gcpMonitor style print to screen -
 * going further requires recomp libSzaMonitor.so which I can't
 * do on my laptop right now...
 *
 * To make even this run requires a file temp.mon containing "read"
 *
 */

#include "gcp/util/Monitor.h"


#include <iostream.h>
#include <math.h>

#include "gcp/util/Debug.h"
#include "gcp/util/RegDescription.h"
#include "gcp/util/RegParser.h"

#include "mex.h"
#include "matrix.h"

using namespace std;
using namespace gcp::util;

void checkArgs(int nlhs, int nrhs);

void readRegs(const mxArray* prhs[], 
	      std::vector<std::vector<std::vector<Monitor::DataType> > >& regVals);

void assignOutputValues(mxArray* plhs[], const mxArray* prhs[],  
			std::vector<std::vector<std::vector<Monitor::DataType> > >& regVals);

mxArray* addRegisterField(mxArray* ptr, RegDescription& reg, unsigned iDim);
mxArray* addNamedStructField(mxArray* parentPtr, std::string fieldName);
mxArray* addNamedCellField(mxArray* parentPtr, std::string fieldName, unsigned iDim);
mxClassID matlabTypeOf(Monitor::DataType& dataType);
mxComplexity matlabComplexityOf(Monitor::DataType& dataType);
void assignRegVals(mxArray* ptr, 
		   RegDescription& reg,
		   unsigned iReg, 
		   std::vector<std::vector<std::vector<Monitor::DataType> > >& regVals);

void assignStringVals(mxArray* array,
		     unsigned iReg,
		      std::vector<std::vector<std::vector<Monitor::DataType> > >& regVals);

void assignComplexFloatVals(mxArray* array,
			    unsigned iReg,
			    std::vector<std::vector<std::vector<Monitor::DataType> > >& regVals);

void assignDateVals(mxArray* array,
		     unsigned iReg,
		      std::vector<std::vector<std::vector<Monitor::DataType> > >& regVals);

mxArray* createMatlabArray(int dims[2], Monitor::DataType& dataType);

void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  int nreg;           // Number of registers being read 
  int i;
  
  // The Array of register values returned by the Montior object
  
  std::vector<std::vector<std::vector<Monitor::DataType> > > regVals;
  
  // Read register values now
  
  readRegs(prhs, regVals);
  
  // Now assign output values
  
  unsigned nFrame = regVals.size();
  
  if(nFrame == 0)
    return;
  
  unsigned nReg = regVals[0].size();
  
  if(nReg == 0)
    return;
  
  assignOutputValues(plhs, prhs, regVals);
}

/**.......................................................................
 * Add a register block to the struct
 */
mxArray* addRegisterField(mxArray* ptr, RegDescription& reg, unsigned iDim)
{
  mxArray* regMapPtr = addNamedStructField(ptr,       reg.regMapName());
  mxArray* boardPtr  = addNamedStructField(regMapPtr, reg.boardName());
  mxArray* regPtr    = addNamedStructField(boardPtr,  reg.blockName());

  if(reg.aspect() != REG_PLAIN) {
    (void)addNamedStructField(regPtr,  reg.aspectName());
    return regPtr;
  } else
    return boardPtr;
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
	      std::vector<std::vector<std::vector<Monitor::DataType> > >& regVals)
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
    
    // Add the requested regs
    
    for(unsigned iReg=0; iReg < nreg; iReg++) 
      monitor.addRegister(mxArrayToString(mxGetCell(prhs[0],iReg)));
    
    // Finally, run this object    
    
    regVals = monitor.readRegsAsDataTypes();
    
  } catch(Exception& err) {
    mexPrintf("%s\n", err.what());
    return;
  }
}

void assignOutputValues(mxArray* plhs[], const mxArray* prhs[],  
			std::vector<std::vector<std::vector<Monitor::DataType> > >& regVals)
{
  int dims[2]={1,1};
  unsigned nFrame = regVals.size();
  
  // Create empty arraymap structure 
  
  plhs[0] = mxCreateStructArray(2, dims, 0, NULL);
  
  // How many registers requested?
  
  unsigned nreg = mxGetNumberOfElements(prhs[0]);
  
  gcp::util::RegParser regParser;
  
  for(unsigned iRegDesc=0, iReg=0; iRegDesc < nreg; iRegDesc++) {
    std::vector<RegDescription> 
      regs=regParser.inputRegs(mxArrayToString(mxGetCell(prhs[0],iRegDesc)));

    for(unsigned iSubReg=0; iSubReg < regs.size(); iSubReg++, iReg++) {
      mxArray* boardPtr = addRegisterField(plhs[0], regs[iSubReg], 1);
      assignRegVals(boardPtr, regs[iSubReg], iReg, regVals);
    }
  }
}

void assignStringVals(mxArray* array,
		     unsigned iReg,
		     std::vector<std::vector<std::vector<Monitor::DataType> > >& regVals)
 {
   for(unsigned iFrame=0; iFrame < regVals.size(); iFrame++)
     mxSetCell(array, iFrame, mxCreateString(regVals[iFrame][iReg][0].stringVal.c_str()));
 }

void assignComplexFloatVals(mxArray* array,
			    unsigned iReg,
			    std::vector<std::vector<std::vector<Monitor::DataType> > >& regVals)
 {
   float* realPtr  = (float*)mxGetData(array);
   float* imagPtr  = (float*)mxGetImagData(array);
   unsigned nEl    = regVals[0][iReg].size();
   unsigned nFrame = regVals.size();

   for(unsigned iFrame=0; iFrame < regVals.size(); iFrame++)
     for(unsigned iEl=0; iEl < nEl; iEl++) {
       *(realPtr + (iEl*nFrame) + iFrame) = (float)regVals[iFrame][iReg][iEl].val.data_.cf.real_;
       *(imagPtr + (iEl*nFrame) + iFrame) = (float)regVals[iFrame][iReg][iEl].val.data_.cf.imag_;
     }
 }

void assignDateVals(mxArray* array,
		     unsigned iReg,
		     std::vector<std::vector<std::vector<Monitor::DataType> > >& regVals)
 {
   for(unsigned iFrame=0; iFrame < regVals.size(); iFrame++) {
     RegDate date(regVals[iFrame][iReg][0].val.data_.date);
     mxSetCell(array, iFrame, mxCreateString(date.str().c_str()));
   }
 }

void assignRegVals(mxArray* ptr, 
		   RegDescription& reg,
		   unsigned iReg, 
		   std::vector<std::vector<std::vector<Monitor::DataType> > >& regVals)
{
#define ASSIGN_VALS(cast, member) \
 {\
      cast* ptr = (cast*)vPtr;\
      for(unsigned iFrame=0; iFrame < nFrame; iFrame++)\
	for(unsigned iEl=0; iEl < nEl; iEl++) {\
	  *(ptr + (iEl*nFrame) + iFrame) = (cast)regVals[iFrame][iReg][iEl].val.data_.member;\
        }\
 }

  unsigned nFrame = regVals.size();
  unsigned nEl    = regVals[0][iReg].size();
  int dims[2]     = {nFrame, nEl};

  Monitor::DataType& dataType = regVals[0][iReg][0];

  // Create and return the array in which values will be stored

  mxArray* array = createMatlabArray(dims, dataType);
  mxSetField(ptr, 0, reg.aspect()==REG_PLAIN ? reg.blockName().c_str() :
	     reg.aspectName().c_str(), array);

  void* vPtr = mxGetData(array);

  switch (dataType.selectedFormat) {
  case Monitor::FM_BOOL:
  case Monitor::FM_UCHAR:
    ASSIGN_VALS(unsigned char, c);
    break;
  case Monitor::FM_CHAR:
    ASSIGN_VALS(char, c);
    break;
  case Monitor::FM_UINT:
    ASSIGN_VALS(unsigned int, ui);
    break;
  case Monitor::FM_INT:
    ASSIGN_VALS(int, i);
    break;
  case Monitor::FM_ULONG:
    ASSIGN_VALS(unsigned long, ul);
    break;
  case Monitor::FM_LONG:
    ASSIGN_VALS(long, l);
    break;
  case Monitor::FM_FLOAT:
    ASSIGN_VALS(float, f);
    break;
  case Monitor::FM_COMPLEX_FLOAT:
    assignComplexFloatVals(array, iReg, regVals);
    break;
  case Monitor::FM_DOUBLE:
    ASSIGN_VALS(double, d);
  case Monitor::FM_STRING:
    assignStringVals(array, iReg, regVals);
    break;
  case Monitor::FM_DATE:
    assignDateVals(array, iReg, regVals);
    break;
  }
}

mxArray* createMatlabArray(int dims[2], Monitor::DataType& dataType)
{
  switch (dataType.selectedFormat) {
  case Monitor::FM_STRING:
  case Monitor::FM_DATE:
    return mxCreateCellArray(2, dims);
  default:
    return mxCreateNumericArray(2, dims, matlabTypeOf(dataType), matlabComplexityOf(dataType));
    break;
  }
}

mxClassID matlabTypeOf(Monitor::DataType& dataType)
{
  switch (dataType.selectedFormat) {
  case Monitor::FM_BOOL:
  case Monitor::FM_UCHAR:
    return mxUINT8_CLASS;
    break;
  case Monitor::FM_CHAR:
    return mxINT8_CLASS;
    break;
  case Monitor::FM_UINT:
    return mxUINT32_CLASS;
    break;
  case Monitor::FM_INT:
    return mxINT32_CLASS;
    break;
  case Monitor::FM_ULONG:
    return mxUINT64_CLASS;
    break;
  case Monitor::FM_LONG:
    return mxINT64_CLASS;
    break;
  case Monitor::FM_FLOAT:
    return mxSINGLE_CLASS;
    break;
  case Monitor::FM_COMPLEX_FLOAT:
    return mxSINGLE_CLASS;
    break;
  case Monitor::FM_DOUBLE:
    return mxDOUBLE_CLASS;
    break;
  }
}

mxComplexity matlabComplexityOf(Monitor::DataType& dataType)
{
  switch (dataType.selectedFormat) {
  case Monitor::FM_COMPLEX_FLOAT:
    return mxCOMPLEX;
    break;
  default:
    return mxREAL;
    break;
  }
}
