#ifndef GCP_UTIL_MEXPARSER_H
#define GCP_UTIL_MEXPARSER_H

/**
 * @file MexParser.h
 * 
 * Tagged: Wed May  4 23:26:43 PDT 2005
 * 
 * @author Erik Leitch
 */
#include <mex.h>
#include <matrix.h>

#include <string>
#include <vector>


// 64-bit versions of matlab define a mwSize typedef to obscure wide
// representations of dimensions, but 32-bit versions define the
// equivalent API with ints.

#ifdef MX_COMPAT_32
typedef int mwSize;
#else
typedef size_t mwSize;
#endif

namespace gcp {
  namespace matlab {
    
    class MexParser {
    public:
      
      /**
       * Constructor.
       */
      MexParser();
      MexParser(const mxArray*);

      void setTo(const mxArray*);
      void setTo(const mxArray**);
      
      /**
       * Destructor.
       */
      virtual ~MexParser();

      bool isDouble();
      static bool isDouble(mxArray* arr);
      static bool isDouble(const mxArray* arr);

      bool isLogical();
      static bool isLogical(mxArray* arr);
      static bool isLogical(const mxArray* arr);

      bool isUchar();
      static bool isUchar(mxArray* arr);
      static bool isUchar(const mxArray* arr);

      bool isUint();
      static bool isUint(mxArray* arr);
      static bool isUint(const mxArray* arr);
     
      bool isString();
      static bool isString(mxArray* arr);
      static bool isString(const mxArray* arr);

      bool isComplex();
      static bool isComplex(mxArray* arr);
      static bool isComplex(const mxArray* arr);

      // Return the number of dimensions in an mxArray

      static int getNumberOfDimensions(mxArray* arr);
      static int getNumberOfDimensions(const mxArray* arr);
      int getNumberOfDimensions();
      
      // Return the vector of dimensions of an mxArray

      static std::vector<int> getDimensions(mxArray* arr);
      static std::vector<int> getDimensions(const mxArray* arr);
      std::vector<int> getDimensions();

      static int getDimension(mxArray* arr, unsigned iDim);
      static int getDimension(const mxArray* arr, unsigned iDim);
      int getDimension(unsigned iDim);

      static unsigned int getNumberOfElements(mxArray* arr);
      static unsigned int getNumberOfElements(const mxArray* arr);
      unsigned int getNumberOfElements();

      int getMaxDimension();

      // Struct handling

      static bool isStruct(mxArray* arr);
      static bool isStruct(const mxArray* arr);
      bool isStruct();

      static unsigned int getNumberOfFields(mxArray* arr);
      static unsigned int getNumberOfFields(const mxArray* arr);
      unsigned int getNumberOfFields();

      unsigned getNumberOfElements(mxArray* arr, std::string fieldName);
      unsigned getNumberOfElements(const mxArray* arr, std::string fieldName);
      unsigned getNumberOfElements(std::string fieldName);

      static const mxArray* getField(const mxArray* arr, std::string fieldName);
      const mxArray* getField(std::string fieldName);

      static bool fieldExists(const mxArray* arr, std::string fieldName);
      bool fieldExists(std::string fieldName);

      static double* getFieldAsDouble(mxArray* arr, std::string fieldName);
      static double* getFieldAsDouble(const mxArray* arr, std::string fieldName);
      double* getFieldAsDouble(std::string fieldName);

      static double* getImagFieldAsDouble(mxArray* arr, std::string fieldName);
      static double* getImagFieldAsDouble(const mxArray* arr, std::string fieldName);
      double* getImagFieldAsDouble(std::string fieldName);

      static unsigned* getFieldAsUint(mxArray* arr, std::string fieldName);
      static unsigned* getFieldAsUint(const mxArray* arr, std::string fieldName);
      unsigned* getFieldAsUint(std::string fieldName);

      static unsigned char* getFieldAsUchar(mxArray* arr, std::string fieldName);
      static unsigned char* getFieldAsUchar(const mxArray* arr, std::string fieldName);
      unsigned char* getFieldAsUchar(std::string fieldName);

      static bool* getFieldAsLogical(mxArray* arr, std::string fieldName);
      static bool* getFieldAsLogical(const mxArray* arr, std::string fieldName);
      bool* getFieldAsLogical(std::string fieldName);

      static std::string getFieldAsString(mxArray* arr, std::string fieldName);
      static std::string getFieldAsString(const mxArray* arr, std::string fieldName);
      std::string getFieldAsString(std::string fieldName);

      void printDimensions();

      // Return pointers to appropriate data types

      static double* getDoubleData(mxArray* arr);
      static double* getDoubleData(const mxArray* arr);
      double* getDoubleData();
      double getDoubleVal(unsigned i1, unsigned i2);
      double getDoubleVal(unsigned i1, unsigned i2, unsigned i3);

      static float* getFloatData(mxArray* arr);
      static float* getFloatData(const mxArray* arr);
      float* getFloatData();

      static bool* getLogicalData(mxArray* arr);
      static bool* getLogicalData(const mxArray* arr);
      bool* getLogicalData();

      static unsigned char* getUcharData(mxArray* arr);
      static unsigned char* getUcharData(const mxArray* arr);
      unsigned char* getUcharData();

      static unsigned* getUintData(mxArray* arr);
      static unsigned* getUintData(const mxArray* arr);
      unsigned* getUintData();
      unsigned int getUintVal(unsigned i1, unsigned i2);
      unsigned int getUintVal(unsigned i1, unsigned i2, unsigned i3);

      static double* getImagDoubleData(mxArray* arr);
      static double* getImagDoubleData(const mxArray* arr);
      double* getImagDoubleData();

      // Return the string corresponding to this mxArray

      static std::string getString(mxArray* arr);
      static std::string getString(const mxArray* arr);
      std::string getString();

      // Return true if the dimensions of two mxArrays match

      static bool dimensionsMatch(mxArray* arr1, mxArray* arr2);
      static bool dimensionsMatch(const mxArray* arr1, const mxArray* arr2);

      bool operator==(MexParser& mp);
      bool operator!=(MexParser& mp);

      // Convert from Matlab mwSize type to int

      static std::vector<int> convertDims(int ndim, const mwSize* dims);

    private:
      
      const mxArray* array_;

    }; // End class MexParser
    
  } // End namespace matlab
} // End namespace gcp

#endif // End #ifndef GCP_UTIL_MEXPARSER_H
