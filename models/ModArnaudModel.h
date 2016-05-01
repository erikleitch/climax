// $Id: $

#ifndef GCP_MODELS_MODARNAUDMODEL_H
#define GCP_MODELS_MODARNAUDMODEL_H

/**
 * @file ModArnaudModel.h
 * 
 * Tagged: Tue Jul 24 15:06:55 PDT 2012
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/models/GnfwModel.h"

#include "gcp/pgutil/PgModel.h"

#include "gcp/util/Mass.h"

//-----------------------------------------------------------------------
// Define the model from Arnaud et al 2010, A&A 517, A92
//-----------------------------------------------------------------------

namespace gcp {
  namespace models {

    class ModArnaudModel : public GnfwModel {
    public:

      /**
       * Constructor.
       */
      ModArnaudModel();

      /**
       * Destructor.
       */
      virtual ~ModArnaudModel();

      void fillM500();
      void fillRadioNormalization();

      double massPrefactorFirstOrder(double h70);
      double massExponent(double xg);

      double getGnfwNormalization();

      virtual gcp::util::PgModel pgModel();

      virtual double radialRadioModel(double r, void* params);
      virtual double radialXrayModel(double r, void* params);

      gcp::util::VariableUnitQuantity rfrac_;

      double radioRat_;
      double xrayRat_;

    }; // End class ModArnaudModel

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_MODARNAUDMODEL_H
