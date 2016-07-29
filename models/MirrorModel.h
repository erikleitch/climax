// $Id: $

#ifndef GCP_MODELS_MIRRORMODEL_H
#define GCP_MODELS_MIRRORMODEL_H

/**
 * @file MirrorModel.h
 * 
 * Tagged: Thu Jan 15 11:14:53 NZDT 2015
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include "gcp/fftutil/Generic2DAngularModel.h"

#include "gcp/util/Angle.h"
#include "gcp/util/Length.h"
#include "gcp/util/VariableUnitQuantity.h"

namespace gcp {
  namespace models {

    class MirrorModel : public gcp::util::Generic2DAngularModel {
    public:

      /**
       * Constructor.
       */
      MirrorModel();

      /**
       * Destructor.
       */
      virtual ~MirrorModel();

      void setNumberOfReferencePoints(unsigned n);
      void setNumberOfFiducialPoints(unsigned n);
      void setParameter(std::string name, std::string val, std::string units, bool external=true);

      void fillArray(unsigned type, gcp::util::Angle& axisUnits, 
		     std::valarray<double>& x, std::valarray<double>& y, 
		     std::valarray<double>& d, 
		     void* params);

    private:

      // Positions of the fiducial balls, in fixed coordinate system

      std::vector<gcp::util::VariableUnitQuantity> xF_;
      std::vector<gcp::util::VariableUnitQuantity> yF_;
      std::vector<gcp::util::VariableUnitQuantity> zF_;

      // Nominal positions of the reference points, in same coordinate system

      std::vector<gcp::util::VariableUnitQuantity> xR_;
      std::vector<gcp::util::VariableUnitQuantity> yR_;
      std::vector<gcp::util::VariableUnitQuantity> zR_;

      // Model positions of the reference points, in the same coordinate system

      std::vector<gcp::util::VariableUnitQuantity> xRm_;
      std::vector<gcp::util::VariableUnitQuantity> yRm_;
      std::vector<gcp::util::VariableUnitQuantity> zRm_;

      // Position of the CR in coordinate system

      gcp::util::VariableUnitQuantity xCr_;
      gcp::util::VariableUnitQuantity yCr_;
      gcp::util::VariableUnitQuantity zCr_;

      // Translation/rotations of the mirror coordinate system about
      // the nominal CR position

      gcp::util::VariableUnitQuantity xTrans_;
      gcp::util::VariableUnitQuantity yTrans_;
      gcp::util::VariableUnitQuantity zTrans_;

      gcp::util::Angle xRot_;
      gcp::util::Angle yRot_;
      gcp::util::Angle zRot_;

    }; // End class MirrorModel

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_MIRRORMODEL_H
