#include <iostream>

#include "gcp/python/PyParser.h"
#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/String.h"

#if DIR_HAVE_NUMPY
#include "arrayobject.h"
#endif

using namespace std;

using namespace gcp::python;
using namespace gcp::util;

#define CHECK_ARRAY(arr) \
  {\
    if(arr == 0)\
      ThrowError("Null array received");\
  }

#define CHECK_DICT(arr) \
  {\
    CHECK_ARRAY(arr); \
    if(!isDict(arr)) {\
      ThrowError("Object does not represent a dictionary");\
    }\
  }
    
/**.......................................................................
 * Constructor.
 */
PyParser::PyParser(const PyObject* arr) 
{
  setTo(arr);
}

PyParser::PyParser()
{
  array_ = 0;
}

void PyParser::setTo(const PyObject** arr) 
{
  CHECK_ARRAY(*arr);
  array_ = *arr;
}

void PyParser::setTo(const PyObject* arr) 
{
  CHECK_ARRAY(arr);
  array_ = arr;
}

/**.......................................................................
 * Destructor.
 */
PyParser::~PyParser() {}

/**.......................................................................
 * Return the dimensionality of an PyObject
 */
int PyParser::getNumberOfDimensions(PyObject* arr)
{
  CHECK_ARRAY(arr);

  // If the object is a tuple, or a list or a sequence, it can have
  // multiple dimensions

  if(isArray(arr)) {
    return PyObject_Size(arr);
  } else {
    return 1;
  }
}

int PyParser::getNumberOfDimensions(const PyObject* arr)
{
  return getNumberOfDimensions((PyObject*) arr);
}

int PyParser::getNumberOfDimensions()
{
  getNumberOfDimensions(array_);
}

PyObject* PyParser::getArrayItem(PyObject* arr, Py_ssize_t index)
{
  CHECK_ARRAY(arr);

  if(index+1 > getSize(arr)) {
    ThrowError("Request for element " << index << " from an array of size " << getSize(arr));
  }

  if(isArray(arr)) {

    if(isList(arr)) {
      return PyList_GetItem(arr, index); // No new reference
    } else if(isSequence(arr)) { 
      PyObject* obj = PySequence_GetItem(arr, index); // Creates a new reference
      Py_DECREF(obj); // Remove the new reference.
      return obj;
    } else if(isTuple(arr)) {
      return PyTuple_GetItem(arr, index); // No new reference.
    }

  } else {
    return arr;
  }

}

/**.......................................................................
 * Return the vector of dimensions for this array
 */
std::vector<int> PyParser::getDimensions(PyObject* obj)
{
  CHECK_ARRAY(obj);

  unsigned nDim = getNumberOfDimensions(obj);

  std::vector<int> dims;
  dims.resize(nDim);

  if(isArray(obj)) {
    for(unsigned iDim=0; iDim < nDim; iDim++) {
      PyObject* obj = getArrayItem(obj, iDim);
      dims[iDim] = getSize(obj);
    }
  } else {
    dims[0] = getSize(obj);
  }

  return dims;
}

std::vector<int> PyParser::getDimensions(const PyObject* arr)
{
  return getDimensions((PyObject*) arr);
}

std::vector<int> PyParser::getDimensions()
{
  return getDimensions(array_);
}

/**.......................................................................
 * Return the number of elements in a dimension
 */
int PyParser::getDimension(unsigned iDim)
{
  return getDimension(array_, iDim);
}

int PyParser::getDimension(PyObject* arr, unsigned iDim)
{
  if(iDim > getNumberOfDimensions(arr)-1)
    ThrowError("Dimension: " << iDim << " too large for this array");

  std::vector<int> dims = getDimensions(arr);

  return dims[iDim];
}

int PyParser::getDimension(const PyObject* arr, unsigned iDim)
{
  return getDimension((PyObject*) arr, iDim);
}

unsigned PyParser::getSize()
{
  return getSize(array_);
}

unsigned PyParser::getSize(const PyObject* obj)
{
  return getSize((PyObject*)obj);
}

unsigned PyParser::getSize(PyObject* obj)
{
  return PyObject_Size(obj);
}

/**.......................................................................
 * Return true if the dimensions of two PyObjects match 
 */
bool PyParser::dimensionsMatch(PyObject* arr1, PyObject* arr2)
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

bool PyParser::dimensionsMatch(const PyObject* arr1, const PyObject* arr2)
{
  return dimensionsMatch((PyObject*) arr1, (PyObject*) arr2);
}

/**.......................................................................
 * Return the dimensionality of an PyObject
 */
unsigned int PyParser::getNumberOfElements()
{
  return getNumberOfElements(array_);
}

unsigned int PyParser::getNumberOfElements(const PyObject* arr)
{
  return getNumberOfElements((PyObject*) arr);
}

unsigned int PyParser::getNumberOfElements(PyObject* arr)
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

/**.......................................................................
 * Return the number of dictionary entries
 */
unsigned int PyParser::getNumberOfDictEntries(const PyObject* arr)
{
  return getNumberOfDictEntries((PyObject*)arr);
}

unsigned int PyParser::getNumberOfDictEntries(PyObject* arr)
{
  CHECK_DICT(arr);
  return getSize(arr);
}

/**.......................................................................
 * Print the dimensionality of this array
 */
void PyParser::printDimensions()
{
  std::vector<int> dims;

  dims = getDimensions();

  std::cout << "Dimensions are: ";
  for(unsigned idim=0; idim < dims.size(); idim++)
    std::cout << "[" << dims[idim] << "]";
  std::cout << std::endl;
}

/**.......................................................................
 * Return true if this object represents an array
 */
bool PyParser::isArray()
{
  return isArray(array_);
}

bool PyParser::isArray(const PyObject* arr)
{
  return isArray((PyObject*) arr);
}

/**.......................................................................
 * An array (in my parlance) is a sequence of objects that is not a
 * string -- that is handled separately
 */
bool PyParser::isArray(PyObject* arr)
{
  return (isTuple(arr) || isList(arr) || isSequence(arr)) && !isString(arr);
}

/**.......................................................................
 * Return if this object represents tuple data
 */
bool PyParser::isTuple()
{
  return isTuple(array_);
}

bool PyParser::isTuple(const PyObject* arr)
{
  return isTuple((PyObject*) arr);
}

bool PyParser::isTuple(PyObject* arr)
{
  return PyTuple_Check(arr);
}

/**.......................................................................
 * Return if this object represents list data
 */
bool PyParser::isList()
{
  return isList(array_);
}

bool PyParser::isList(const PyObject* arr)
{
  return isList((PyObject*) arr);
}

bool PyParser::isList(PyObject* arr)
{
  return PyList_Check(arr);
}

/**.......................................................................
 * Return if this object represents sequence data
 */
bool PyParser::isSequence()
{
  return isSequence(array_);
}

bool PyParser::isSequence(const PyObject* arr)
{
  return isSequence((PyObject*) arr);
}

bool PyParser::isSequence(PyObject* arr)
{
  return PySequence_Check(arr);
}

/**.......................................................................
 * Return if this object represents dict data
 */
bool PyParser::isDict()
{
  return isDict(array_);
}

bool PyParser::isDict(const PyObject* arr)
{
  return isDict((PyObject*) arr);
}

bool PyParser::isDict(PyObject* arr)
{
  return PyDict_Check(arr);
}

/**.......................................................................
 * Return if this object represents string data
 */
bool PyParser::isString()
{
  return isString(array_);
}

bool PyParser::isString(const PyObject* arr)
{
  return isString((PyObject*) arr);
}

bool PyParser::isString(PyObject* arr)
{
  return PyString_Check(arr);
}

/**.......................................................................
 * Return if this object represents double data
 */
bool PyParser::isFloat()
{
  return isFloat(array_);
}

bool PyParser::isFloat(const PyObject* arr)
{
  return isFloat((PyObject*) arr);
}

bool PyParser::isFloat(PyObject* arr)
{
  return PyFloat_Check(arr);
}

/**.......................................................................
 * Return if this object represents unsigned data
 */
bool PyParser::isBool()
{
  return isBool(array_);
}

bool PyParser::isBool(PyObject* arr)
{
  return isBool((const PyObject*) arr);
}

bool PyParser::isBool(const PyObject* arr)
{
  return PyBool_Check(arr);
}

/**.......................................................................
 * Return if this object represents integer data
 */
bool PyParser::isInt()
{
  return isInt(array_);
}

bool PyParser::isInt(PyObject* arr)
{
  return isInt((const PyObject*) arr);
}

bool PyParser::isInt(const PyObject* arr)
{
  return PyInt_Check(arr);
}

/**.......................................................................
 * Get a pointer to string data
 */
std::string PyParser::getString(unsigned iEl)
{
  return getString(array_, iEl);
}

std::string PyParser::getString(const PyObject* arr, unsigned iEl)
{
  return getString((PyObject*) arr, iEl);
}

std::string PyParser::getString(PyObject* arr, unsigned iEl)
{
  PyObject* obj = getArrayItem(arr, iEl);

  if(!isString(obj)) {
    ThrowError("Object does not represent a string");
  }

  std::string retval(PyString_AsString(obj));

  return retval;
}

float PyParser::getFloat(unsigned iEl)
{
  return getFloat(array_, iEl);
}

float PyParser::getFloat(const PyObject* arr, unsigned iEl)
{
  return getFloat((PyObject*) arr, iEl);
}

float PyParser::getFloat(PyObject* arr, unsigned iEl)
{
  if(!isFloat(arr)) {
    ThrowError("Object does not represent a float");
  }

  float retval = PyFloat_AS_DOUBLE(arr);

  return retval;
}

int PyParser::getMaxDimension()
{
  int max=0;

  for(unsigned i=0; i < getNumberOfDimensions(); i++) 
    if(getDimension(i) > max)
      max = getDimension(i);

  return max;
}

bool PyParser::operator==(PyParser& mp)
{
  if(getNumberOfDimensions() != mp.getNumberOfDimensions())
    return false;

  for(unsigned i=0; i < getNumberOfDimensions(); i++) {
    if(getDimension(i) != mp.getDimension(i))
      return false;
  }

  return true;
}

bool PyParser::operator!=(PyParser& mp)
{
  if(getNumberOfDimensions() != mp.getNumberOfDimensions())
    return true;

  for(unsigned i=0; i < getNumberOfDimensions(); i++) {
    if(getDimension(i) != mp.getDimension(i))
      return true;
  }

  return false;
}


unsigned PyParser::getUintVal()
{
  return getUintVal(array_);
}

unsigned PyParser::getUintVal(PyObject* arr)
{
  long val = PyInt_AsLong(arr);

  if(!isInt(arr)) {
    ThrowError("Object does not represent integer data");
  }

  if(PyErr_Occurred() && val < 0) {
    ThrowError("Error occurred in conversion");
  }

  if(val < 0) {
    ThrowError("Value cannot be represented as an unsigned: " << val);
  }

  if(val > UINT_MAX) {
    ThrowError("Value is larger than the maximum that can be represented as an unsigned: " << val);
  }

  return (unsigned)val;
}

unsigned PyParser::getUintVal(const PyObject* arr)
{
  return getUintVal((PyObject*)arr);
}

bool PyParser::getBoolVal()
{
  return getBoolVal(array_);
}

bool PyParser::getBoolVal(const PyObject* arr)
{
  return getBoolVal((PyObject*)arr);
}

bool PyParser::getBoolVal(PyObject* arr)
{
  if(!isBool(arr)) {
    ThrowError("Object does not represent boolean data");
  }

  return (bool)getUintVal(arr);
}

PyParser::PythonType PyParser::pythonTypeOf(gcp::util::DataType::Type dataType)
{
  switch (dataType) {
  case DataType::BOOL:
  case DataType::UCHAR:
  case DataType::CHAR:
  case DataType::USHORT:
  case DataType::SHORT:
  case DataType::UINT:
  case DataType::INT:
  case DataType::ULONG:
  case DataType::LONG:
    return PY_INT;
    break;
  case DataType::FLOAT:
  case DataType::DOUBLE:
  case DataType::DATE:
    return PY_FLOAT;
    break;
  case DataType::STRING:
    return PY_STRING;
    break;
  case DataType::COMPLEX_FLOAT:
  case DataType::COMPLEX_DOUBLE:
    return PY_COMPLEX;
    break;
  default:
    return PY_UNKNOWN;
    break;
  }
}

double* PyParser::getNumpyDoublePtr(PyObject* arr, bool directData)
{
#if DIR_HAVE_NUMPY
  if(directData) {
    return (double*)((PyArrayObject*)arr)->data;
  } else {
    std::vector<int> dims = getNumpyDimensions(arr);
    std::vector<npy_intp> npyInd(dims.size());
    
    for(unsigned iDim=0; iDim < dims.size(); iDim++) {
      npyInd[iDim]  = 0;
    }

    return (double*)PyArray_GetPtr((PyArrayObject*)arr, &npyInd[0]);
  }
#else
  return 0;
#endif
}

/**.......................................................................
 * Get the number of dimensions for a numpy object
 */
int PyParser::getNumpyLength(PyObject* arr, unsigned index)
{
  int n = 0;

#if DIR_HAVE_NUMPY
  if(index > getNumpyNDimensions(arr)-1) {
    ThrowError("Request for element " << index << " from an array of size " << getNumpyNDimensions(arr));
  }

  npy_intp* d = ((PyArrayObject*)(arr))->dimensions;
  n = d[index];
#endif

  return n;
}

/**.......................................................................
 * Get the number of dimensions for a numpy object
 */
std::vector<int> PyParser::getNumpyDimensions(PyObject* arr)
{
  std::vector<int> dims;

#if DIR_HAVE_NUMPY
  dims.resize(getNumpyNDimensions(arr));
  npy_intp* d = ((PyArrayObject*)(arr))->dimensions;
  for(unsigned i=0; i < dims.size(); i++) 
    dims[i] = d[i];
#endif

  return dims;
}

/**.......................................................................
 * Get the number of dimensions for a numpy object
 */
unsigned PyParser::getNumpyNDimensions(PyObject* arr)
{
#if DIR_HAVE_NUMPY
  return ((PyArrayObject*)(arr))->nd;
#else
  return 0;
#endif
}
