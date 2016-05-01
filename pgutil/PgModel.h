// $Id: $

#ifndef GCP_UTIL_PGMODEL_H
#define GCP_UTIL_PGMODEL_H

/**
 * @file PgModel.h
 * 
 * Tagged: Mon Oct 21 10:23:01 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Angle.h"

#include "gcp/pgutil/Trans.h"

namespace gcp {
  namespace util {

    class PgModel {
    public:

      enum Type {
	TYPE_NONE,
	TYPE_DELTA,  // Delta-function model
	TYPE_BETA,   // Beta-model
	TYPE_ELBETA, // Elliptical Beta-model
	TYPE_ARNAUD, // Arnaud-model
	TYPE_GAUSS,  // Gaussian model
	TYPE_BOX,    // Box
	TYPE_BEAM,   // Beam
      };

      /**
       * Constructor.
       */
      PgModel();
      PgModel(PgModel& mod);
      PgModel(const PgModel& mod);

      /**
       * Destructor.
       */
      virtual ~PgModel();

      void operator=(const PgModel& mod);
      void operator=(PgModel& mod);

      void rectify();
      void precalculateShape();
      void draw();
      bool hasTwoAxes();
      int getColor();

      void print(std::string& unit, unsigned iMod, Trans& trans);
      std::string formatPosition(std::string varName, double pos);

      Type type_;

      float peak_;

      float xMid_;
      float yMid_;

      float xRad1_;
      float yRad1_;

      float xRad2_;
      float yRad2_;

      float rad1_;
      float rad2_;

      float radMin_;
      float radMax_;

      bool drawCenter_;
      bool fill_;

      float xMin_;
      float xMax_;
      float yMin_;
      float yMax_;

      gcp::util::Angle rot_;
      gcp::util::Angle angle_;

      std::vector<float> x_;
      std::vector<float> y_;

      static Type keyToType(char key);

    private:

      void render();

    }; // End class PgModel

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PGMODEL_H
