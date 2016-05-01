#ifndef GCP_MATLAB_MEXHANDLER_H
#define GCP_MATLAB_MEXHANDLER_H

/**
 * @file MexHandler.h
 * 
 * Tagged: Thu Sep 15 09:37:48 PDT 2005
 * 
 * 
 * @author Erik Leitch
 */
#include <iostream>
#include <cmath>

#include "gcp/util/Debug.h"
#include "gcp/util/DataType.h"
#include "gcp/util/ErrHandler.h"
#include "gcp/util/Logger.h"

#include "gcp/matlab/MexParser.h"

#include "mex.h"
#include "matrix.h"

#include <sstream>
#include <vector>

#define MXPRINT(statement) \
{\
  std::ostringstream os;\
  os.str("");\
  os << statement << std::endl;\
  mexPrintf(os.str().c_str());\
}

#define MXPRINTN(statement) \
{\
  std::ostringstream os;\
  os.str("");\
  os << statement; \
  mexPrintf(os.str().c_str());\
}


namespace gcp {
  
  namespace util {
    class FitsBinTableReader;
  }
  
  namespace matlab {
    
    class MexHandler {
    public:
      
      // Structure to encapsulate a mex array and the data pointer
      // associated with it.
      
      struct MxArray {
	mxArray* array_;
	void* vPtr_;
	
	MxArray() {
	  array_ = 0;
	  vPtr_  = 0;
	};
	
	// Constructor calls the asignment operator
	
	MxArray(mxArray* array) {
	  *this = array;
	};
	
	// Assignment operator
	
	void operator=(mxArray* array) {
	  array_ = array;
	  vPtr_  = mxGetData(array);
	};
      };
      
      /**
       * Constructor.
       */
      MexHandler();
      
      /**
       * Destructor.
       */
      virtual ~MexHandler();
      
      static mxArray* createMatlabArray(int ndim, std::vector<int> dims, gcp::util::DataType::Type dataType);
      static mxArray* createMatlabArray(int ndim, const int* dims, gcp::util::DataType::Type dataType);
      static mxArray* createMatlabArray(int length, gcp::util::DataType::Type dataType);

      // Assign a return mxArray and return a double pointer to it
      
      static double* createDoubleArray(mxArray** mxArray, int ndim, const int* dims);
      static double* createDoubleArray(mxArray** mxArray, int ndim, std::vector<int> dims);
      static double* createDoubleArray(mxArray** mxArray, unsigned len);
      static float*  createFloatArray(mxArray** mxArray, int ndim, const int* dims);
      static double* createDoubleArray(mxArray** mxArray, MexParser& parser);
      
      static mxClassID matlabTypeOf(gcp::util::DataType::Type dataType);
      
      static mxComplexity matlabComplexityOf(gcp::util::DataType::Type dataType);
      
      /**.......................................................................
       * Create empty matlab arrays to be filled in by the monitor stream
       */
      static mxArray* addHierNamedStructField(mxArray* parentPtr, 
					      std::string fieldName,
					      gcp::util::DataType::Type=gcp::util::DataType::NONE,
					      unsigned nDim=0,
					      int* dims=0);

      static mxArray* addNamedStructField(mxArray* parentPtr, std::string fieldName, unsigned nElement=1, unsigned index=0);
      static mxArray* addNamedCellField(mxArray* parentPtr, std::string fieldName, unsigned iDim, unsigned index=0);
      
      static char* addNamedStringStructField(mxArray* parentPtr, std::string fieldName,
							 unsigned len, unsigned index=0);
      
      static char* addNamedStringStructField(mxArray* parentPtr, std::string fieldName,
							 std::string str, unsigned index=0);
      
      static mxArray* addNamedStructField(mxArray* parentPtr, std::string fieldName, int ndim, const int* dims, gcp::util::DataType::Type dataType, unsigned index=0);
      static mxArray* addNamedStructField(mxArray* parentPtr, std::string fieldName, int ndim, std::vector<int>& dims, gcp::util::DataType::Type dataType, unsigned index=0);

      static double* addNamedDoubleStructField(mxArray* parentPtr, 
					       std::string fieldName, 
					       unsigned nElement, unsigned index=0);
      
      static double* addNamedDoubleStructField(mxArray* parentPtr, 
					       std::string fieldName, 
					       std::vector<int> dims, unsigned index=0);
      
      static float* addNamedFloatStructField(mxArray* parentPtr, std::string fieldName,
					     std::vector<int> dims, unsigned index=0);

      static float* addNamedFloatStructField(mxArray* parentPtr, 
					     std::string fieldName, 
					     unsigned nElement, unsigned index=0);
      
      static unsigned* addNamedUintStructField(mxArray* parentPtr, 
					       std::string fieldName, 
					       unsigned nElement, unsigned index=0);
      
      // Add a named double struct field that has the same dimensions
      // as the field named in copyName

      static double* copyNamedDoubleStructField(mxArray* parentPtr, std::string fieldName, std::string copyName);

      // Create a matlab array appropriate for the numbered column
      // from a FITS table
      
      // Add a named field to the passed structure, and create an
      // properly-typed array for it
      
      static void* addNamedStructField(mxArray* parentPtr, std::string fieldName,
				       gcp::util::FitsBinTableReader* reader, unsigned iCol, unsigned index);
      
      static mxArray* createMatlabArray(gcp::util::FitsBinTableReader* reader, 
					unsigned colNo);
      
      static void checkArgs(int nlhsExpected, int nrhsExpected, 
			    int nlhsActual,   int nrhsActual);
      
      static LOG_HANDLER_FN(stdoutPrintFn);
      static LOG_HANDLER_FN(stderrPrintFn);
      
      static ERR_HANDLER_FN(throwFn);
      static ERR_HANDLER_FN(reportFn);
      static ERR_HANDLER_FN(logFn);

      // Utility functions for converting from a "multi-dimensional"
      // array coordinate into a one-dimensional representation.

      static unsigned indexStartingWithFastest(std::vector<unsigned>& coord, 
					       std::vector<unsigned>& dims);

      static unsigned indexStartingWithSlowest(std::vector<unsigned>& coord, 
					       std::vector<unsigned>& dims);

      // Return an array of C-order idices for the input dimensions.
      // Assumes dims is in order from fastest to slowest changing
      // indices
      
      static void getIndicesC(std::vector<unsigned>& cVec, 
			      std::vector<unsigned>& cDims);


      static unsigned getMatlabIndex(unsigned n1, unsigned n2,
				     unsigned i1, unsigned i2);

      static unsigned getMatlabIndex(unsigned n1, unsigned n2, unsigned n3,
				     unsigned i1, unsigned i2, unsigned i3);

      // Return an array of C-order idices for the input dimensions.
      // rAssumes dims is in order from fastest to slowest changing
      // indices
      
      static void getIndicesMatlab(std::vector<unsigned>& mVec, 
				   std::vector<unsigned>& mDims);
	
      static void getIndicesMatlab(std::vector<unsigned>& mVec, 
				   std::vector<int>& mDims);

      // Recursive function to return an array of C-order and
      // Matlab-order idices for the input dimensions.  Note that the
      // input dimension array must be in C-style (fastest to slowest)
      // dimension order
      
      static void getIndices(std::vector<unsigned>& cVec, std::vector<unsigned>& mVec,
			     std::vector<unsigned>& dims, 
			     int iDim=0, unsigned indLast=0, 
			     unsigned nLast=0, unsigned baseLast=0);

    private:

      static std::vector<mwSize> convertDims(int ndim, const int* dims);
	
    }; // End class MexHandler
    
  } // End namespace matlab
} // End namespace gcp



#endif // End #ifndef GCP_MATLAB_MEXHANDLER_H
