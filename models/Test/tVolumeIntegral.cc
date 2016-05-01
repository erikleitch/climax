#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Timer.h"

#include "gcp/models/ArnaudModel.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "d",        "0", "d", "distance at which the 3D model is centered"},
  { "r",        "0", "d", "cylindrical radius at which to evaluate the line integral"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};


namespace gcp {
  namespace models {

    class TestModel : public GenericRadiallySymmetric3DModel {
    public:

      /**
       * Constructor.
       */
      TestModel() {};

      /**
       * Destructor.
       */
      virtual ~TestModel() {};

      // Define the radial model to be a solid sphere (constant value at all r)

      double radialModel(unsigned type, double x, void* params) {
	return 1.0;
      }

      // Return true if shape parameters are fixed

      bool shapeParametersAreFixed() {
	return true;
      }

    }; // End class TestModel

  } // End namespace models
} // End namespace gcp

int Program::main()
{
  TestModel model;

  model.nInterp_     = 100;
  model.deltaInterp_ = 0.01;

  model.calculateVolumeInterpolationValues(DataSetType::DATASET_RADIO);

  COUT("integral = " << model.interpolateVolumeIntegral(DataSetType::DATASET_RADIO, 0.1));
  COUT("integral = " << model.interpolateVolumeIntegral(DataSetType::DATASET_RADIO, 0.2));
  COUT("integral = " << model.interpolateVolumeIntegral(DataSetType::DATASET_RADIO, 0.4));
  COUT("integral = " << model.interpolateVolumeIntegral(DataSetType::DATASET_RADIO, 0.6));
  COUT("integral = " << model.interpolateVolumeIntegral(DataSetType::DATASET_RADIO, 0.8));
  COUT("integral = " << model.interpolateVolumeIntegral(DataSetType::DATASET_RADIO, 1.0));

  return 0;
}
