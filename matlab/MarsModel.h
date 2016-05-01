// $Id: MarsModel.h,v 1.1.1.1 2010/07/13 17:56:46 eml Exp $

#ifndef GCP_MATLAB_MARSMODEL_H
#define GCP_MATLAB_MARSMODEL_H

/**
 * @file MarsModel.h
 * 
 * Tagged: Thu Sep 15 10:07:49 PDT 2005
 * 
 * @version: $Revision: 1.1.1.1 $, $Date: 2010/07/13 17:56:46 $
 * 
 * @author Erik Leitch
 */
#include "mex.h"
#include "matrix.h"

#include "gcp/matlab/MexParser.h"
#include "gcp/util/DataType.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/TimeVal.h"

namespace gcp {
  namespace matlab {

    class MarsModel {
    public:

      /**
       * Constructor.
       */
      MarsModel();

      /**
       * Destructor.
       */
      virtual ~MarsModel();

      void setFileName(const mxArray* array);
      void setType(const mxArray* array);
      void setDirectory(const mxArray* array);
      void setMjdArray(const mxArray* array);
      void setFrequencyArray(const mxArray* array);
      void getData(mxArray** mxData, mxArray** mxErrorCode);

      enum Type {
	TEMP,
	FLUX,
	SIZE,
	EDIAM,
	PDIAM
      };

    private:

      MexParser mpFlux_;
      MexParser mpMjd_;
      std::string fileName_;
      std::string dir_;

      Type type_;

      double* freq_;
      unsigned int nFreq_;

      double* mjd_;
      unsigned int nMjd_;

    }; // End class MarsModel

  } // End namespace matlab
} // End namespace gcp



#endif // End #ifndef GCP_MATLAB_MARSMODEL_H
