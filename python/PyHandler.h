// $Id: PyHandler.h,v 1.10 2008/11/22 02:11:21 eml Exp $

#ifndef GCP_PYTHON_PYHANDLER_H
#define GCP_PYTHON_PYHANDLER_H

/**
 * @file PyHandler.h
 * 
 * Tagged: Thu Sep 15 09:37:48 PDT 2005
 * 
 * @version: $Revision: 1.10 $, $Date: 2008/11/22 02:11:21 $
 * 
 * @author Erik Leitch
 */
#include <iostream>
#include <cmath>

#include "gcp/util/CoordAxes.h"
#include "gcp/util/Debug.h"
#include "gcp/util/Directives.h"
#include "gcp/util/DataType.h"
#include "gcp/util/ErrHandler.h"
#include "gcp/util/Logger.h"

#include <Python.h>

#if DIR_HAVE_NUMPY
#include "arrayobject.h"
#endif

#include <sstream>

namespace gcp {
  namespace python {
    
    class PyHandler {
    public:
      
      /**
       * Constructor.
       */
      PyHandler();
      
      /**
       * Destructor.
       */
      virtual ~PyHandler();
      
      static PyObject* getObject(gcp::util::DataType::Type dataType);

      static PyObject* getFlatList(unsigned len,
				   gcp::util::DataType::Type dataType=gcp::util::DataType::NONE,
				   char** data=0);

      static PyObject* getFlatList(std::vector<unsigned int> dims,
				   gcp::util::DataType::Type dataType=gcp::util::DataType::NONE,
				   char** data=0);

      static PyObject* getFlatList(std::vector<int> dims,
				   gcp::util::DataType::Type dataType=gcp::util::DataType::NONE,
				   char** data=0);

      static PyObject* getFlatList(int ndim, const int* dims, 
				   gcp::util::DataType::Type dataType=gcp::util::DataType::NONE,
				   char** data=0);

      // Methods for creating tuples

      static void addTuple(PyObject** obj, std::vector<unsigned> dims, unsigned iDim, 
			   gcp::util::DataType::Type dataType=gcp::util::DataType::NONE);

      static PyObject* getTuple(std::vector<unsigned> dims, gcp::util::DataType::Type dataType=gcp::util::DataType::NONE);
      static PyObject* getTuple(int nDim, const int* dims,  gcp::util::DataType::Type dataType=gcp::util::DataType::NONE);

      static void addList(PyObject** obj, std::vector<unsigned> dims, unsigned iDim, 
			   gcp::util::DataType::Type dataType=gcp::util::DataType::NONE);

      static PyObject* getList(std::vector<unsigned> dims, gcp::util::DataType::Type dataType=gcp::util::DataType::NONE);
      static PyObject* getList(int nDim, const int* dims,  gcp::util::DataType::Type dataType=gcp::util::DataType::NONE);
      
      static std::vector<unsigned> convertToVec(int nDim, const int* dims);
      static std::vector<unsigned> convertToVec(std::vector<int> dims);
      static std::vector<unsigned> convertToVec(gcp::util::CoordAxes& axes);

      static PyObject* getArray(unsigned len, gcp::util::DataType::Type type, char** data=0);

#if DIR_HAVE_NUMPY
      static int numpyTypeOf(gcp::util::DataType::Type dataType);
      static PyObject* getNumPyArray(std::vector<int> dims, gcp::util::DataType::Type dataType, char** data=0);
      static PyObject* getNumPyArray(unsigned len, gcp::util::DataType::Type dataType, char** data=0);
      static std::vector<npy_intp> getNumPyInds(std::vector<int>& dims);
#endif

    }; // End class PyHandler
    
  } // End namespace python
} // End namespace gcp



#endif // End #ifndef GCP_PYTHON_PYHANDLER_H
