// $Id: Generic1DModel.h,v 1.2 2012/05/12 00:45:52 eml Exp $

#ifndef GCP_UTIL_GENERIC1DMODEL_H
#define GCP_UTIL_GENERIC1DMODEL_H

/**
 * @file Generic1DModel.h
 * 
 * Tagged: Fri May  4 10:51:29 PDT 2012
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/12 00:45:52 $
 * 
 * @author Erik Leitch
 */
#include "gcp/fftutil/Model.h"

namespace gcp {
  namespace util {

    class Generic1DModel : public gcp::util::Model {
    public:

      class ExecData {
      public:
	
	Generic1DModel* model_;
	std::vector<double>* x_;
	std::vector<double>* y_;
	unsigned iSegment_;
	unsigned nSegment_;
	unsigned iStart_;
	unsigned iStop_;
	unsigned i_;

	void* evalData_;

	ExecData(Generic1DModel* model) {
	  initialize();
	  model_    = model;

	  // And allocate eval data now

	  evalData_ = model->allocateEvalData();
	};

	void initialize()
	{
	  model_    = 0;
	  x_        = 0;
	  y_        = 0;
	  iSegment_ = 0;
	  nSegment_ = 0;
	  iStart_   = 0;
	  iStop_    = 0;
	  evalData_ = 0;
	};

      };

      /**
       * Constructor.
       */
      Generic1DModel();

      /**
       * Destructor.
       */
      virtual ~Generic1DModel();

      // Return the value of a 1D model at the specified point

      virtual double eval(double x);
      virtual double eval(void* evalData, double x);

      // Method to fill a 1D array with this model

      void fill1DArray(std::vector<double>& x, std::vector<double>& y);

      void fill1DArraySingleThread(std::vector<double>& x, std::vector<double>& y);
      void fill1DArrayMultiThread(ExecData* ed);
      unsigned initializeExecData(std::vector<double>& x, std::vector<double>& y);
      static EXECUTE_FN(execFill1DArrayMultiThread);

      // Method to write fake 1D data from this model to a file

      void generateFake1DData(std::string fileName, std::vector<double>& x, std::vector<double>& y, double sigma, std::string units);

      // Inherited method to initialize our thread pool

      void setThreadPool(ThreadPool* pool);

      // If inheritors support multi-threaded execution, they should
      // define this method to return an object containing the stack needed to
      // evaluat the envelope function

      virtual void* allocateEvalData();
      virtual void initializeEvalData(void* evalData);

      std::vector<ExecData*> execData_;

    }; // End class Generic1DModel

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_GENERIC1DMODEL_H
