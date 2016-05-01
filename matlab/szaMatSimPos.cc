/**.......................................................................
 * MATLAB Mex file for generating non-overlapping array positions
 * locations, a source declination and a list of HA
 *
 * Use like:
 *
 * [e n] = gcpMatSimPos(niter, nant, length, minSep, [e, n, fix])
 *
 */
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"
#include "gcp/array/code/share/slalib/slalib.h"
#include "gcp/util/Coordinates.h"

#include "mex.h"
#include "matrix.h"

#include <iostream.h>
#include <math.h>

using namespace std;
using namespace gcp::util;
using namespace gcp::matlab;

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

  //------------------------------------------------------------
  // Parse input argument
  //------------------------------------------------------------

  MexParser parser(prhs[0]);

  if(!parser.isStruct()) {
    std::ostringstream os;
    os << "Usage: sp = gcpMatSimPos(p)\n\n"
       << "Where p is a struct, with fields: " << std::endl << std::endl
       << " niter       = number of iterations\n"
       << " nant        = number of antennas\n"
       << " length      = length of the area on a side (meters)\n"
       << " minSep(nloc)= minimum allowed separation (meters)\n"
       << " dofix       = true to fix antenna locations\n"
       << " fix(nloc)   = true to fix the location of an antenna to a certain position\n"
       << " e(nloc)     = a list of east locations\n"
       << " n(nloc)     = a list of north locations\n"
       << " dosnap      = true to snap to locations\n"
       << " snap(nloc)  = if true, snap to this location when closer than minSep\n"
       << "               if false, avoid this location when closer than minSep\n"
       << " nSnapMin    = if dosnap=true, only return configurations where nsnap >= nSnapMin"
       << " nSnapMax    = if dosnap=true, only return configurations where nsnap <= nSnapMax";

    ThrowError(os.str());
  }

  //------------------------------------------------------------
  // Deal with input parsing
  //------------------------------------------------------------

  unsigned nIter = (unsigned)(*parser.getFieldAsDouble("niter"));
  unsigned nAnt  = (unsigned)(*parser.getFieldAsDouble("nant"));
  double length  =            *parser.getFieldAsDouble("length");
  double* minSepPtr=           parser.getFieldAsDouble("minsep");
  bool dofix     =     (bool)(*parser.getFieldAsDouble("dofix"));
  bool dosnap    =     (bool)(*parser.getFieldAsDouble("dosnap"));

  unsigned nFixedAnt = 0;
  std::vector<unsigned> fixedAntLocations(nAnt);

  double* eastPtr    = 0;
  double* northPtr   = 0;

  unsigned nLoc=0, nSnap=0;

  if(dofix || dosnap) {

    eastPtr  = parser.getFieldAsDouble("e");
    northPtr = parser.getFieldAsDouble("n");
    nLoc     = parser.getNumberOfElements("e");

    if(nLoc != parser.getNumberOfElements("n")) {
      ThrowError("e and n arrays must be of the same length");
    }
  }

  std::vector<unsigned> occupied(nLoc);
  std::vector<unsigned> snapVec(nLoc);

  //------------------------------------------------------------
  // Get the fix array information
  //------------------------------------------------------------

  if(dofix) {
    double* dPtr  = parser.getFieldAsDouble("fix");

    if(nLoc != parser.getNumberOfElements("fix")) {
      ThrowError("e, n and fix arrays must be of the same length");
    }

    nFixedAnt = 0;
    for(unsigned iLoc=0; iLoc < nLoc; iLoc++) {
      if(dPtr[iLoc] > 0) {
	fixedAntLocations[nFixedAnt++] = iLoc;
      }
    }
  }

  //------------------------------------------------------------
  // Get the snap array information
  //------------------------------------------------------------

  unsigned* snapPtr  = 0;  
  unsigned nSnapMin = 0, nSnapMax = 0;
  unsigned nSnapFound = 0;

  if(dosnap) {
    double* dPtr  = parser.getFieldAsDouble("snap");

    if(nLoc != parser.getNumberOfElements("snap")) {
      ThrowError("e, n and snap arrays must be of the same length");
    }

    snapPtr  = &snapVec[0];

    for(unsigned iLoc=0; iLoc < nLoc; iLoc++) {
      snapPtr[iLoc] = (dPtr[iLoc] > 0);
    }

    nSnap = nLoc;

    nSnapMin = (unsigned)(*(parser.getFieldAsDouble("nSnapMin")));
    nSnapMax = (unsigned)(*(parser.getFieldAsDouble("nSnapMax")));
  }

  //------------------------------------------------------------
  // Create output structure
  //------------------------------------------------------------

  std::vector<int> dims(3);
  dims[0] = 1;
  dims[1] = 1;

  plhs[0] = mxCreateStructArray(2, &dims[0], 0, NULL);

  dims.resize(2);
  dims[0] = nIter;
  dims[1] = nAnt;

  double* eOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],"e", dims);
  double* nOutPtr = MexHandler::addNamedDoubleStructField(plhs[0],"n", dims);

  //------------------------------------------------------------
  // Variable we need for testing separations
  //------------------------------------------------------------

  double cmpVal;

  //------------------------------------------------------------
  // Main loop.  Iterate until we have found niter valid sets of
  // positions
  //------------------------------------------------------------

  double e, n, de, dn;
  unsigned iIter, iLoc, locInd;
  unsigned antInd, iAnt, outputInd, iPosAnt, outputPosInd;
  unsigned antPosInd;

  for(unsigned iIter=0; iIter < nIter; iIter++) {

    //------------------------------------------------------------
    // Initialize all known locations to be unoccupied
    //------------------------------------------------------------

    for(iLoc=0; iLoc < nSnap; iLoc++) 
      occupied[iLoc] = false;

    //------------------------------------------------------------
    // Start by assigning any fixed antennas for this iteration
    //------------------------------------------------------------

    for(iAnt=0; iAnt < nFixedAnt; iAnt++) {
      
      outputInd = (nIter * iAnt) + iIter;

      locInd = fixedAntLocations[iAnt];

      eOutPtr[outputInd] = eastPtr[locInd];
      nOutPtr[outputInd] = northPtr[locInd];
      occupied[locInd]   = true;
    }

    //------------------------------------------------------------
    // Now generate positions just for the free antennas
    //------------------------------------------------------------

    nSnapFound = 0;
    for(iAnt=nFixedAnt; iAnt < nAnt; iAnt++) {

      outputInd = (nIter * iAnt) + iIter;

      bool valid = false;

      while(!valid) {

	bool stop  = false;

	e = (double)(rand())/RAND_MAX * length - length/2;
	n = (double)(rand())/RAND_MAX * length - length/2;
	
	valid = true;

	//------------------------------------------------------------
	// Check against all specified locations
	//------------------------------------------------------------
	
	for(iLoc=0; iLoc < nSnap && !stop; iLoc++) {
	  
	  de = e - eastPtr[iLoc];
	  dn = n - northPtr[iLoc];
	  
	  // If this location is within minSep of the current location
	  
	  cmpVal = minSepPtr[iLoc]*minSepPtr[iLoc];

	  if(de*de + dn*dn < cmpVal) {
	    
	    // Then snap to this point if not already occupied
	    
	    if(snapPtr[iLoc] && !occupied[iLoc]) {
	      
	      // If there is no limit on the max number of matching
	      // positions, or the number already found is less than
	      // the max number allowed, assign this position to a
	      // known antenna pad

	      if(nSnapMax == 0 || nSnapFound < nSnapMax) {
		eOutPtr[outputInd] = eastPtr[iLoc];
		nOutPtr[outputInd] = northPtr[iLoc];
		occupied[iLoc] = true;
		valid = true;
		++nSnapFound;
	      } else {
		valid = false;
	      }

	    } else {
	      valid = false;
	    }

	    // If we found a position within minSep, then we have
	    // either included this position, or excluded it, so we
	    // should stop checking, and either generate a new
	    // position for this antenna, or continue on to the next
	    
	    stop = true; 
	  }
	}

	//------------------------------------------------------------
	// If the position is still valid after vetting against all
	// known locations, check against the other antennas that
	// already have positions
	//------------------------------------------------------------

	if(valid && !stop) {

	  for(iPosAnt=0; iPosAnt < iAnt && !stop; iPosAnt++) {
      
	    outputPosInd = (nIter * iPosAnt) + iIter;
	    
	    de = e - eOutPtr[outputPosInd];
	    dn = n - nOutPtr[outputPosInd];

	    // If this location is too close to an existing antenna,
	    // flag the position as invalid (and break out of this
	    // loop)

	    if(de*de + dn*dn < cmpVal) {
	      valid = false;
	      stop  = true;
	    }
	  }
	}

	// If this position is still valid, assign it to the output
	// array

	if(valid && !stop) {

	  int nLeft = nAnt - iAnt - 1;
	  int nSnapLeft = nSnapMin - nSnapFound;

	  // If we are not requiring any positions to match the input
	  // array, or we have already satisfied the minimum matching
	  // condition, or the number of antennas we still have to
	  // find positions for is >= the number of matched positions
	  // we still need to find
#if 1
	  if(nSnapMin == 0 || nSnapFound >= nSnapMin || nLeft >= nSnapLeft) {
	    eOutPtr[outputInd] = e;
	    nOutPtr[outputInd] = n;
	  } else {
	    valid = false;
	  }
#else
	  eOutPtr[outputInd] = e;
	  nOutPtr[outputInd] = n;
#endif

	}

      } // End while(!valid)

    } // End loop over free antennas

  } // End loop over iterations
  
  return;
}
