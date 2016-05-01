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

#include <iostream.h>
#include <math.h>

#include "gcp/util/Debug.h"
#include "gcp/util/Directives.h"
#include "gcp/util/Monitor.h"

#include "mex.h"
#include "matrix.h"

using namespace std;
using namespace gcp::util;

#include "gcp/util/Test/Program.h"

void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  char* directory;    // The directory in which to look for files 
  char* calfile;      // The name of the calibration file 
  char* start_date;   // Start time as date:time string 
  char* end_date;     // End time as date:time string 
  int nreg;           // Number of registers being read 
  int i;

  // Get the command-line arguments from matlab 

  if(nrhs!=5)
    mexErrMsgTxt("Must have 5 input arguments");

  if(nlhs!=1)
    mexErrMsgTxt("Must have 1 output arguments");

  // How many registers requested?

  nreg = mxGetNumberOfElements(prhs[0]);
  
  // Get start/end date:time strings 

  start_date = mxArrayToString(prhs[1]);
  end_date   = mxArrayToString(prhs[2]);

  // Get arcdir and calfile locations 
 
  directory  = mxArrayToString(prhs[3]);
  calfile    = mxArrayToString(prhs[4]);

  Debug::setLevel(0);

  std::vector<std::vector<double> > regVals;

  try {

    Monitor monitor(directory, calfile, start_date, end_date);

    // Add the requested regs

    for(i=0;i<nreg;i++) 
      monitor.addRegister(mxArrayToString(mxGetCell(prhs[0],i)));

    // Finally, run this object    

    regVals = monitor.readRegsAsDoubles();

  } catch(Exception& err) {
    mexPrintf("%s\n", err.what());
    return;
  }

  unsigned nRow = regVals.size();

  if(nRow == 0)
    return;

  unsigned nCol = regVals[0].size();

  if(nCol == 0)
    return;

  plhs[0] = mxCreateDoubleMatrix(nRow, nCol, mxREAL);

  double* outputPtr = mxGetPr(plhs[0]);

  for(unsigned iRow=0; iRow < nRow; iRow++)
    for(unsigned iCol=0; iCol < nCol; iCol++)
      outputPtr[iCol * nRow + iRow] = regVals[iRow][iCol];

  return;
}
