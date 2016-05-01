// $Id: ModelType.h,v 1.2 2012/05/31 16:22:56 eml Exp $

#ifndef GCP_UTIL_MODELTYPE_H
#define GCP_UTIL_MODELTYPE_H

/**
 * @file ModelType.h
 * 
 * Tagged: Thu Apr 26 13:25:22 PDT 2012
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/31 16:22:56 $
 * 
 * @author Erik Leitch
 */
namespace gcp {
  namespace util {

    class ModelType {
    public:

      enum {
	MODEL_UNKNOWN        = 0x0,
	MODEL_ADDITIVE       = 0x1,  // Additive model
	MODEL_MULTIPLICATIVE = 0x2,  // Multiplicative model
      };

    }; // End class ModelType

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_MODELTYPE_H
