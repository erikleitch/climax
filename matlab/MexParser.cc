#include <iostream>

#include "gcp/matlab/MexParser.h"

#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/String.h"

using namespace std;

using namespace gcp::matlab;
using namespace gcp::util;

#define CHECK_ARRAY(arr) \
  {\
    if(arr == 0)\
      ThrowError("Null array received");\
  }

/**.......................................................................
 * Constructor.
 */
MexParser::MexParser(const mxArray* arr) 
{
  setTo(arr);
}

MexParser::MexParser()
{
  array_ = 0;
}

void MexParser::setTo(const mxArray** arr) 
{
  CHECK_ARRAY(*arr);
  array_ = *arr;
}

void MexParser::setTo(const mxArray* arr) 
{
  CHECK_ARRAY(arr);
  array_ = arr;
}

/**.......................................................................
 * Destructor.
 */
MexParser::~MexParser() {}


/**.......................................................................
 * Return the dimensionality of an mxArray
 */
int MexParser::getNumberOfDimensions(mxArray* arr)
{
  CHECK_ARRAY(arr);

  return mxGetNumberOfDimensions(arr);
}

int MexParser::getNumberOfDimensions(const mxArray* arr)
{
  return getNumberOfDimensions((mxArray*) arr);
}

int MexParser::getNumberOfDimensions()
{
  getNumberOfDimensions(array_);
}

/**.......................................................................
 * Return the vector of dimensions for this array
 */
std::vector<int> MexParser::getDimensions(mxArray* arr)
{
  CHECK_ARRAY(arr);
  return convertDims(mxGetNumberOfDimensions(arr), mxGetDimensions(arr));
}

std::vector<int> MexParser::getDimensions(const mxArray* arr)
{
  return getDimensions((mxArray*) arr);
}

std::vector<int> MexParser::getDimensions()
{
  return getDimensions(array_);
}

/**.......................................................................
 * Return the number of elements in a dimension
 */
int MexParser::getDimension(unsigned iDim)
{
  return getDimension(array_, iDim);
}

int MexParser::getDimension(mxArray* arr, unsigned iDim)
{
  if(iDim > getNumberOfDimensions(arr)-1)
    ThrowError("Dimension: " << iDim << " too large for this array");

  std::vector<int> dims = getDimensions(arr);

  return dims[iDim];
}

int MexParser::getDimension(const mxArray* arr, unsigned iDim)
{
  return getDimension((mxArray*) arr, iDim);}

/**.......................................................................
 * Return true if the dimensions of two mxArrays match 
 */
bool MexParser::dimensionsMatch(mxArray* arr1, mxArray* arr2)
{
  // If the number of dimensions don't agree, they don't match

  int ndim;
  if((ndim=getNumberOfDimensions(arr1)) != getNumberOfDimensions(arr2))
    return false;

  // Else if any one dimension has a different length, they don't match

  std::vector<int> dims1 = getDimensions(arr1);
  std::vector<int> dims2 = getDimensions(arr2);

  for(unsigned idim = 0; idim < ndim; idim++)
    if(dims1[idim] != dims2[idim])
      return false;

  // Else they match

  return true;
}

bool MexParser::dimensionsMatch(const mxArray* arr1, const mxArray* arr2)
{
  return dimensionsMatch((mxArray*) arr1, (mxArray*) arr2);
}

/**.......................................................................
 * Return the dimensionality of an mxArray
 */
unsigned int MexParser::getNumberOfElements()
{
  return getNumberOfElements(array_);
}

unsigned int MexParser::getNumberOfElements(mxArray* arr)
{
  CHECK_ARRAY(arr);

  bool first=true;
  unsigned nEl=0;

  for(unsigned iDim=0; iDim < getNumberOfDimensions(arr); iDim++) {

    if(first) {
      nEl = getDimension(arr, iDim);
      first = false;
    } else {
      nEl *= getDimension(arr, iDim);
    };
  }

  return nEl;
}

unsigned int MexParser::getNumberOfElements(const mxArray* arr)
{
  return getNumberOfElements((mxArray*) arr);
}

/**.......................................................................
 * Return true if the array is a struct
 */
bool MexParser::isStruct(const mxArray* arr)
{
  CHECK_ARRAY(arr);
  return mxIsStruct(arr);
}

bool MexParser::isStruct(mxArray* arr)
{
  CHECK_ARRAY(arr);
  return mxIsStruct((const mxArray*) arr);
}

bool MexParser::isStruct()
{
  return isStruct(array_);
}

/**.......................................................................
 * Return the dimensionality of an mxArray
 */
unsigned int MexParser::getNumberOfFields()
{
  return getNumberOfFields(array_);
}

/**.......................................................................
 * Return the dimensionality of an mxArray
 */
unsigned int MexParser::getNumberOfFields(mxArray* arr)
{
  CHECK_ARRAY(arr);

  bool first=true;
  unsigned nEl=0;

  for(unsigned iDim=0; iDim < getNumberOfDimensions(arr); iDim++) {

    if(first) {
      nEl = getDimension(arr, iDim);
      first = false;
    } else {
      nEl *= getDimension(arr, iDim);
    };
  }

  return nEl;
}

unsigned int MexParser::getNumberOfFields(const mxArray* arr)
{
  return getNumberOfFields((mxArray*) arr);
}

unsigned MexParser::getNumberOfElements(std::string fieldName)
{
  return getNumberOfElements(mxGetField(array_, 0, fieldName.c_str()));
}

unsigned MexParser::getNumberOfElements(mxArray* arr, std::string fieldName) 
{
  return getNumberOfElements((const mxArray*) arr, fieldName);
}

unsigned MexParser::getNumberOfElements(const mxArray* arr, std::string fieldName)
{
  return getNumberOfElements(mxGetField(arr, 0, fieldName.c_str()));
}

const mxArray* MexParser::getField(string fieldName)
{
  return getField(array_, fieldName);
}

const mxArray* MexParser::getField(const mxArray* arr, std::string fieldName)
{
  String name(fieldName);
  String next;
  const mxArray* arrPtr = arr;

  do {
    next = name.findNextStringSeparatedByChars(".");

    if(!next.isEmpty()) {

      arrPtr = mxGetField(arrPtr, 0, next.str().c_str());

      // If the array is empty, this field doesn't exist
      
      if(arrPtr == 0)
	return 0;
    }

  } while(!next.isEmpty());

  // Else, return true

  return arrPtr;
}

bool MexParser::fieldExists(string fieldName)
{
  return fieldExists(array_, fieldName);
}

bool MexParser::fieldExists(const mxArray* arr, std::string fieldName)
{
  String name(fieldName);
  String next;
  const mxArray* arrPtr = arr;

  do {
    next = name.findNextStringSeparatedByChars(".");

    if(!next.isEmpty()) {
      arrPtr = getField(arrPtr, next.str());

      // If the array is empty, this field doesn't exist
      
      if(arrPtr == 0)
	return false;
    }

  } while(!next.isEmpty());

  // Else, return true

  return true;
}

unsigned char* MexParser::getFieldAsUchar(mxArray* arr, std::string fieldName)
{
  return getFieldAsUchar((const mxArray*) arr, fieldName);
}

unsigned char* MexParser::getFieldAsUchar(const mxArray* arr, std::string fieldName)
{
  return getUcharData(getField(arr, fieldName));
}

unsigned char* MexParser::getFieldAsUchar(std::string fieldName)
{
  return getFieldAsUchar(array_, fieldName);
}

bool* MexParser::getFieldAsLogical(std::string fieldName)
{
  return getFieldAsLogical(array_, fieldName);
}

bool* MexParser::getFieldAsLogical(mxArray* arr, std::string fieldName)
{
  return getFieldAsLogical((const mxArray*) arr, fieldName);
}

bool* MexParser::getFieldAsLogical(const mxArray* arr, std::string fieldName)
{
  return getLogicalData(getField(arr, fieldName));
}

double* MexParser::getFieldAsDouble(mxArray* arr, std::string fieldName)
{
  return getFieldAsDouble((const mxArray*)arr, fieldName);
}

double* MexParser::getFieldAsDouble(const mxArray* arr, std::string fieldName)
{
  return getDoubleData(getField(arr, fieldName));
}

double* MexParser::getFieldAsDouble(std::string fieldName)
{
  return getFieldAsDouble(array_, fieldName);
}

double* MexParser::getImagFieldAsDouble(mxArray* arr, std::string fieldName)
{
  return getImagFieldAsDouble((const mxArray*)arr, fieldName);
}

double* MexParser::getImagFieldAsDouble(const mxArray* arr, std::string fieldName)
{
  return getImagDoubleData(getField(arr, fieldName));
}

double* MexParser::getImagFieldAsDouble(std::string fieldName)
{
  return getImagFieldAsDouble(array_, fieldName);
}

unsigned* MexParser::getFieldAsUint(mxArray* arr, std::string fieldName)
{
  return getFieldAsUint((const mxArray*) arr, fieldName);
}

unsigned* MexParser::getFieldAsUint(const mxArray* arr, std::string fieldName)
{
  return getUintData(getField(arr, fieldName));
}

unsigned* MexParser::getFieldAsUint(std::string fieldName)
{
  return getFieldAsUint(array_, fieldName);
}

std::string MexParser::getFieldAsString(mxArray* arr, std::string fieldName)
{
  return getFieldAsString((const mxArray*) arr, fieldName);
}

std::string MexParser::getFieldAsString(const mxArray* arr, std::string fieldName)
{
  return getString(getField(arr, fieldName));
}

std::string MexParser::getFieldAsString(std::string fieldName)
{
  return getFieldAsString(array_, fieldName);
}

/**.......................................................................
 * Printe the dimensionality of this array
 */
void MexParser::printDimensions()
{
  std::vector<int> dims = getDimensions();

  std::cout << "Dimensions are: ";
  for(unsigned idim=0; idim < dims.size(); idim++)
    std::cout << "[" << dims[idim] << "]";
  std::cout << std::endl;
}

/**.......................................................................
 * Return if this object represents string data
 */
bool MexParser::isString()
{
  return isString(array_);
}

bool MexParser::isString(const mxArray* arr)
{
  return isString((mxArray*) arr);
}

bool MexParser::isString(mxArray* arr)
{
  return mxIsChar(arr);
}

/**.......................................................................
 * Return if this object represents complex data
 */
bool MexParser::isComplex()
{
  return isComplex(array_);
}

bool MexParser::isComplex(const mxArray* arr)
{
  return isComplex((mxArray*) arr);
}

bool MexParser::isComplex(mxArray* arr)
{
  return mxIsComplex(arr);
}

/**.......................................................................
 * Return if this object represents double data
 */
bool MexParser::isDouble()
{
  return isDouble(array_);
}

bool MexParser::isDouble(const mxArray* arr)
{
  return isDouble((mxArray*) arr);
}

bool MexParser::isDouble(mxArray* arr)
{
  return mxIsDouble(arr);
}

/**.......................................................................
 * Get a pointer to imaginary data, as a double array
 */
double* MexParser::getImagDoubleData()
{
  return getImagDoubleData(array_);
}

double* MexParser::getImagDoubleData(mxArray* arr)
{
  CHECK_ARRAY(arr);

  if(mxIsDouble(arr)) {

    if(mxIsComplex(arr)) {
      return (double*)mxGetImagData(arr);
    } else {
      ThrowError("Array does not represent complex data");
    }
  } else {
    ThrowError("Array does not represent double precision data");
  }
}

double* MexParser::getImagDoubleData(const mxArray* arr)
{
  return getImagDoubleData((mxArray*)arr);
}

/**.......................................................................
 * Get a pointer to double data
 */
double* MexParser::getDoubleData()
{
  return getDoubleData(array_);
}

double* MexParser::getDoubleData(mxArray* arr)
{
  CHECK_ARRAY(arr);

  if(mxIsDouble(arr))
    return (double*)mxGetData(arr);

  ThrowError("Array does not represent double precision data");
}

double* MexParser::getDoubleData(const mxArray* arr)
{
  return getDoubleData((mxArray*)arr);
}

/**.......................................................................
 * Return if this object represents unsigned data
 */
bool MexParser::isLogical()
{
  return isLogical(array_);
}

bool MexParser::isLogical(mxArray* arr)
{
  return isLogical((const mxArray*) arr);
}

bool MexParser::isLogical(const mxArray* arr)
{
  return mxIsLogical(arr);
}

/**.......................................................................
 * Return if this object represents unsigned data
 */
bool MexParser::isUchar()
{
  return isUchar(array_);
}

bool MexParser::isUchar(mxArray* arr)
{
  return isUchar((const mxArray*) arr);
}

bool MexParser::isUchar(const mxArray* arr)
{
  return mxIsUint8(arr);
}

/**.......................................................................
 * Return if this object represents unsigned data
 */
bool MexParser::isUint()
{
  return isUint(array_);
}

bool MexParser::isUint(mxArray* arr)
{
  return isUint((const mxArray*) arr);
}

bool MexParser::isUint(const mxArray* arr)
{
  return mxIsUint32(arr);
}

/**.......................................................................
 * Get a pointer to unsigned data
 */
unsigned char* MexParser::getUcharData()
{
  return getUcharData(array_);
}

unsigned char* MexParser::getUcharData(mxArray* arr)
{
  CHECK_ARRAY(arr);

  if(isUchar(arr))
    return (unsigned char*)mxGetData(arr);

  ThrowError("Array does not represent unsigned data");
}

unsigned char* MexParser::getUcharData(const mxArray* arr)
{
  return getUcharData((mxArray*)arr);
}

/**.......................................................................
 * Get a pointer to logical data
 */
bool* MexParser::getLogicalData()
{
  return getLogicalData(array_);
}

bool* MexParser::getLogicalData(mxArray* arr)
{
  CHECK_ARRAY(arr);

  if(isLogical(arr))
    return (bool*)mxGetData(arr);

  ThrowError("Array does not represent unsigned data");
}

bool* MexParser::getLogicalData(const mxArray* arr)
{
  return getLogicalData((mxArray*)arr);
}

/**.......................................................................
 * Get a pointer to unsigned data
 */
unsigned* MexParser::getUintData()
{
  CHECK_ARRAY(array_);

  if(isUint(array_))
    return (unsigned*)mxGetData(array_);

  ThrowError("Array does not represent unsigned data");
}

unsigned* MexParser::getUintData(mxArray* arr)
{
  CHECK_ARRAY(arr);

  if(isUint(arr))
    return (unsigned*)mxGetData(arr);

  ThrowError("Array does not represent unsigned data");
}

unsigned* MexParser::getUintData(const mxArray* arr)
{
  return getUintData((mxArray*)arr);
}

/**.......................................................................
 * Get a pointer to float data
 */
float* MexParser::getFloatData()
{
  return getFloatData(array_);
}

float* MexParser::getFloatData(mxArray* arr)
{
  CHECK_ARRAY(arr);

  if(mxIsSingle(arr))
    return (float*)mxGetData(arr);

  ThrowError("Array does not represent float precision data");
}

float* MexParser::getFloatData(const mxArray* arr)
{
  return getFloatData((mxArray*)arr);
}

/**.......................................................................
 * Get a pointer to string data
 */
std::string MexParser::getString()
{
  return getString(array_);
}

std::string MexParser::getString(mxArray* arr)
{
  CHECK_ARRAY(arr);

  if(mxIsChar(arr)) {
    std::string retVal(mxArrayToString(arr));
    return retVal;
  }

  ThrowError("Array does not represent a string");
}

std::string MexParser::getString(const mxArray* arr)
{
  return getString((mxArray*) arr);
}

int MexParser::getMaxDimension()
{
  int max=0;

  for(unsigned i=0; i < getNumberOfDimensions(); i++) 
    if(getDimension(i) > max)
      max = getDimension(i);

  return max;
}

bool MexParser::operator==(MexParser& mp)
{
  if(getNumberOfDimensions() != mp.getNumberOfDimensions())
    return false;

  for(unsigned i=0; i < getNumberOfDimensions(); i++) {
    if(getDimension(i) != mp.getDimension(i))
      return false;
  }

  return true;
}

bool MexParser::operator!=(MexParser& mp)
{
  if(getNumberOfDimensions() != mp.getNumberOfDimensions())
    return true;

  for(unsigned i=0; i < getNumberOfDimensions(); i++) {
    if(getDimension(i) != mp.getDimension(i))
      return true;
  }

  return false;
}

std::vector<int> MexParser::convertDims(int ndim, const mwSize* dims)
{
  std::vector<int> mtDims;
  mtDims.resize(ndim);

  for(unsigned i=0; i < ndim; i++) {
    mtDims[i] = dims[i];
  }

  return mtDims;
}
