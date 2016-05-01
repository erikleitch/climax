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

      // Define the radial model

      double radialModel(unsigned type, double x, void* params) {
	double val = 1.0/(1.0 + x*x);
	return val*val;
      }

      // Return true if shape parameters are fixed

      bool shapeParametersAreFixed() {
	return true;
      }

      /**.......................................................................
       * Use the current nInterp_ and deltaInterp_ to compute interpolation
       * values for the passed type
       */
      void calculateInterpolationValues(unsigned type)
      {
	if(nInterp_ != interpolatedLineIntegral_[type].size()) {
	  interpolatedLineIntegral_[type].resize(nInterp_);
	  arealIntegral_[type].resize(nInterp_);
	}
	
	double x;
	double norm = lineIntegral(type, 0.0);
	
	for(unsigned i=0; i < nInterp_; i++) {
	  x = i * deltaInterp_;
	  
	  //------------------------------------------------------------
	  // Get the line integral at the specified point
	  //------------------------------------------------------------

	  interpolatedLineIntegral_[type][i] = 1.0/((1.0 + x*x)*(1.0 + x*x));
	  
	  //------------------------------------------------------------
	  // Since the profiles for this class by definition don't depend on
	  // angle (but only radius), we can calculate the areal integration
	  // of this function at the same time by performing a simple line
	  // integral of the function we've just evaluated.  
	  //
	  // Given f(r), we want the integral:
	  //
	  //          1     /x
	  // F(x) = ------  |  f(r) 2pi*r dr
	  //        pi*x^2  /0
	  //
	  //------------------------------------------------------------
	  
	  if(i == 0) {
	    
	    //------------------------------------------------------------
	    //  Integral at zero radius is zero
	    //------------------------------------------------------------
	    
	    arealIntegral_[type][i] = 0.0;
	    
	  } else {
	    
	    //------------------------------------------------------------
	    // Take the value of the function at the midpoint
	    //------------------------------------------------------------
	    
	    unsigned i2 = i;
	    unsigned i1 = i-1;
	    x -= deltaInterp_/2;
	    
	    //------------------------------------------------------------
	    // For the current value of i, we take f(r) * 2*pi*r * dr, with r a
	    // the midpoint of the last two samples
	    //------------------------------------------------------------
	    
	    double val = (interpolatedLineIntegral_[type][i2] + interpolatedLineIntegral_[type][i1])/2 * 2 * M_PI * x * deltaInterp_;
	    arealIntegral_[type][i] = val + arealIntegral_[type][i-1];
	  }
	}
	
	//------------------------------------------------------------
	// Now we normalize the integral by the area at each radius
	//------------------------------------------------------------
	
	for(unsigned i=0; i < nInterp_; i++) {
	  x = i * deltaInterp_;
	  arealIntegral_[type][i] = arealIntegral_[type][i]/(M_PI*x*x);
	}
      };

    }; // End class TestModel

  } // End namespace models
} // End namespace gcp

int Program::main()
{
  TestModel model;
  model.getVar("thetaCore")->setVal(1.0, "'");
  model.getVar("Sradio")->setVal(1.0, "Jy");
  model.getVar("thetaCore")->wasSpecified_ = true;
  model.getVar("Sradio")->wasSpecified_ = true;

  Angle size;
  size.setDegrees(0.2);
  Image image;
  image.xAxis().setNpix(512);
  image.xAxis().setAngularSize(size);
  image.yAxis().setNpix(512);
  image.yAxis().setAngularSize(size);
  
  COUT("Here 0");
  model.calculateInterpolationValuesForCurrentScaleRadius(DataSetType::DATASET_RADIO, image);
  COUT("Here 1");
  model.fillImage(DataSetType::DATASET_RADIO, image);
  COUT("Here 2");

  image.display();

  COUT("Here 3");

  COUT("integral = " << model.interpolateArealIntegral(DataSetType::DATASET_RADIO, 1.0));
  COUT("integral = " << model.interpolateArealIntegral(DataSetType::DATASET_RADIO, 2.0));
  COUT("integral = " << model.interpolateArealIntegral(DataSetType::DATASET_RADIO, 3.0));
  COUT("integral = " << model.interpolateArealIntegral(DataSetType::DATASET_RADIO, 5.0));
  COUT("integral = " << model.interpolateArealIntegral(DataSetType::DATASET_RADIO, 10.0));

  COUT("Here 4");

  return 0;
}
