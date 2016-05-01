#include "gcp/python/PyHandler.h"
#include "gcp/python/PyParser.h"

#include <sstream>

using namespace std;

using namespace gcp::python;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
PyHandler::PyHandler() {}

/**.......................................................................
 * Destructor.
 */
PyHandler::~PyHandler() {}

/**.......................................................................
 * Get a list of the specified length
 */
PyObject* PyHandler::getObject(gcp::util::DataType::Type dataType)
{
  PyParser::PythonType pythonType = PyParser::pythonTypeOf(dataType);

  switch (pythonType) {
  case PyParser::PY_INT:
    return PyInt_FromLong(0);
    break;
  case PyParser::PY_FLOAT:
    return PyFloat_FromDouble(0.0);
    break;
  case PyParser::PY_STRING:
    return PyString_FromString("");
    break;
  case PyParser::PY_COMPLEX:
    return PyComplex_FromDoubles(0.0, 0.0);
    break;
  default:    
    return 0;
    break;
  }
}

PyObject* PyHandler::getFlatList(unsigned len,
				 gcp::util::DataType::Type dataType,
				 char** data)
{
  std::vector<int> dims(1);
  dims[0] = len;
  return getFlatList(dims, dataType, data);
}

/**.......................................................................
 * Get a list of the specified length
 */
PyObject* PyHandler::getFlatList(std::vector<unsigned int> dims,
				 gcp::util::DataType::Type dataType,
				 char** data)
{
  // For now, convert multi-dimensional arrays to 1-D arrays

  unsigned nEl=dims[0];
  for(unsigned iDim=1; iDim < dims.size(); iDim++) {
    nEl *= dims[iDim];
  }

  // Here we just create a List that is nEl long, with typed objects
  // to be filled in later

  PyObject* obj = PyList_New(nEl);

  if(obj == 0) {
    ThrowError("Unable to allocate new list of length: " << nEl);
  }

  if(data != 0)
    *data = (char*)obj;

  return obj;
}

/**.......................................................................
 * Other calling sequences
 */
PyObject* PyHandler::getFlatList(std::vector<int> dims,
				 gcp::util::DataType::Type dataType,
				 char** data)
{
  return getFlatList(convertToVec(dims), dataType, data);
}

PyObject* PyHandler::getFlatList(int nDim, const int* dims, 
				 gcp::util::DataType::Type dataType,
				 char** data)
{
  return getFlatList(convertToVec(nDim, dims), dataType, data);
}

/**.......................................................................
 * Add a tuple of the specified length to each element of the passed
 * object (if the object is a tuple)
 */
void PyHandler::addTuple(PyObject** obj, std::vector<unsigned> dims, unsigned iDim, gcp::util::DataType::Type dataType)
{
  // If we've exhausted the dimensionality of this array, return

  if(iDim >= dims.size()) {

    // If no datatype was specified, just return with a blank tuple

    if(dataType != DataType::NONE) {

      if(*obj == 0) {
	ThrowError("addTuple got passed a null object");
      }

      unsigned nItem = PyTuple_Size(*obj);
      for(unsigned iItem=0; iItem < nItem; iItem++) {
	PyObject* item = getObject(dataType);
	PyTuple_SetItem(*obj, iItem, item);
      }
    }

    return;
  }

  // Else capture the size of the current dimension

  unsigned nEl = dims[iDim];

  // If this is the first call, allocate a new tuple and assign it to
  // the passed object

  if(iDim==0) {
    *obj = PyTuple_New(dims[iDim]);
    addTuple(obj, dims, iDim+1, dataType);
  } else {
    
    // Iterate over the current dimension of the passed object, adding a
    // tuple of length nel to each element of it
    
    unsigned nItem = PyTuple_Size(*obj);
    
    for(unsigned iItem=0; iItem < nItem; iItem++) {
      PyObject* tuple = PyTuple_New(nEl);
      addTuple(&tuple, dims, iDim+1, dataType);
      PyTuple_SetItem(*obj, iItem, tuple);
    }
  }
}

/**.......................................................................
 * Return a tuple sized to match the passed dimensions
 */
PyObject* PyHandler::getTuple(std::vector<unsigned> dims, gcp::util::DataType::Type dataType)
{
  PyObject* obj = 0;
  addTuple(&obj, dims, 0, dataType);
  return obj;
}

/**.......................................................................
 * Return a tuple sized to match the passed dimensions
 */
PyObject* PyHandler::getTuple(int nDim, const int* dims, gcp::util::DataType::Type dataType)
{
  return getTuple(convertToVec(nDim, dims), dataType);
}

/**.......................................................................
 * Add a list of the specified length to each element of the passed
 * object (if the object is a list)
 */
void PyHandler::addList(PyObject** obj, std::vector<unsigned> dims, unsigned iDim, gcp::util::DataType::Type dataType)
{
  // If we've exhausted the dimensionality of this array, return

  if(iDim >= dims.size()) {

    // If no datatype was specified, just return with a blank list

    if(dataType != DataType::NONE) {

      if(*obj == 0) {
	ThrowError("addList got passed a null object");
      }

      unsigned nItem = PyList_Size(*obj);
      for(unsigned iItem=0; iItem < nItem; iItem++) {
	PyObject* item = getObject(dataType);
	PyList_SetItem(*obj, iItem, item);
      }
    }

    return;
  }

  // Else capture the size of the current dimension

  unsigned nEl = dims[iDim];

  // If this is the first call, allocate a new list and assign it to
  // the passed object

  if(iDim==0) {
    *obj = PyList_New(dims[iDim]);
    addList(obj, dims, iDim+1, dataType);
  } else {
    
    // Iterate over the current dimension of the passed object, adding a
    // list of length nel to each element of it
    
    unsigned nItem = PyList_Size(*obj);
    
    for(unsigned iItem=0; iItem < nItem; iItem++) {
      PyObject* listObj = PyList_New(nEl);
      addList(&listObj, dims, iDim+1, dataType);
      PyList_SetItem(*obj, iItem, listObj);
    }
  }
}

/**.......................................................................
 * Return a list sized to match the passed dimensions
 */
PyObject* PyHandler::getList(std::vector<unsigned> dims, gcp::util::DataType::Type dataType)
{
  PyObject* obj = 0;
  addList(&obj, dims, 0, dataType);
  return obj;
}

/**.......................................................................
 * Return a list sized to match the passed dimensions
 */
PyObject* PyHandler::getList(int nDim, const int* dims, gcp::util::DataType::Type dataType)
{
  return getList(convertToVec(nDim, dims), dataType);
}

//=======================================================================
// Conversion functions
//=======================================================================

/**.......................................................................
 * Conversion functions
 */
std::vector<unsigned> PyHandler::convertToVec(int nDim, const int* dims)
{
  std::vector<unsigned> indVec(nDim);

  for(unsigned iDim=0; iDim < nDim; iDim++) {
    indVec[iDim] = dims[iDim];
  }

  return indVec;
}

/**.......................................................................
 * Conversion functions
 */
std::vector<unsigned> PyHandler::convertToVec(CoordAxes& axes)
{
  std::vector<unsigned> indVec(axes.nAxis());

  for(unsigned iAxis=0; iAxis < axes.nAxis(); iAxis++) {
    indVec[iAxis] = axes.nEl(iAxis);
  }

  return indVec;
}

/**.......................................................................
 * Convert from int vec to uint vec
 */
std::vector<unsigned> PyHandler::convertToVec(std::vector<int> dims)
{
  return convertToVec(dims.size(), &dims[0]);
}

#if DIR_HAVE_NUMPY
int PyHandler::numpyTypeOf(gcp::util::DataType::Type dataType)
{
  switch (dataType) {
  case DataType::BOOL:
    return PyArray_BOOL;
    break;
  case DataType::UCHAR:
    return PyArray_UBYTE;
    break;
  case DataType::CHAR:
    return PyArray_BYTE;
    break;
  case DataType::USHORT:
    return PyArray_USHORT;
    break;
  case DataType::SHORT:
    return PyArray_SHORT;
    break;
  case DataType::UINT:
    return PyArray_UINT;
    break;
  case DataType::INT:
    return PyArray_INT;
    break;
  case DataType::ULONG:
    return PyArray_ULONG;
    break;
  case DataType::LONG:
    return PyArray_LONG;
    break;
  case DataType::FLOAT:
    return PyArray_FLOAT;
    break;
  case DataType::DOUBLE:
    return PyArray_DOUBLE;
    break;
  case DataType::STRING:
    return PyArray_STRING;
    break;
  case DataType::DATE:
    return PyArray_DOUBLE;
    break;
  case DataType::COMPLEX_FLOAT:
    return PyArray_CFLOAT;
    break;
  case DataType::COMPLEX_DOUBLE:
    return PyArray_CDOUBLE;
    break;
  default:
    return PyArray_UBYTE;
    break;
  }
}

PyObject* PyHandler::getNumPyArray(unsigned len, gcp::util::DataType::Type dataType, char** data)
{
  std::vector<int> dims(1);
  dims[0] = len;
  return getNumPyArray(dims, dataType, data);
}

PyObject* PyHandler::getNumPyArray(std::vector<int> dims, gcp::util::DataType::Type dataType, char** data)
{
  COUT("Here 0");
  int type = numpyTypeOf(dataType);
  COUT("Here 1");
  std::vector<npy_intp> npyInds = getNumPyInds(dims);
  COUT("Here 2: dims.size() = " << dims.size());
  for(unsigned i=0; i < dims.size(); i++) {
    COUT("Dims[" << i << "] = " << dims[i]);
  }

  COUT("Calling pyarray from dims with dims.size() = " << dims.size() << " dims[0] = " << dims[0] << " type = " << type);
  PyObject* obj = PyArray_FromDims(dims.size(), &dims[0], type);

  COUT("Here 3");
  if(data != 0) {
  COUT("Here 4");
    *data = (char*)PyArray_GetPtr((PyArrayObject*)obj, &npyInds[0]);
      COUT("Here 5");
  }

  return obj;
}

/**.......................................................................
 * Return the dimensions of this register block
 */
std::vector<npy_intp> PyHandler::getNumPyInds(std::vector<int>& dims)
{
  std::vector<npy_intp> npyInds(dims.size());

  for(unsigned iDim=0; iDim < dims.size(); iDim++)
    npyInds[iDim]  = 0;

  return npyInds;
}
#endif

PyObject* PyHandler::getArray(unsigned len, DataType::Type type, char** data)
{
  PyObject* obj = 0;

#if DIR_HAVE_NUMPY
  obj = getNumPyArray(len, type, data);
#else
  obj = PyHandler::getFlatList(len, type, data);
#endif

  return obj;
}
