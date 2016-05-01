#include "gcp/models/PowerlawProfile.h"

#include "gcp/fftutil/DataSetType.h"

#include "gcp/util/Constants.h"
#include "gcp/util/Exception.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

using namespace std;
using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
PowerlawProfile::PowerlawProfile() 
{
  addParameter("n", DataType::UINT, "The number of radial inflection segments");

  fac_.allowUnitless(true);
  fac_.setVal(1.0, "");
  addComponent(fac_);
  addComponentName(fac_, "fac", "Additive factor");

  initializeComponentsToFixed();
  norm_ = 1.0;

  general_.addParameter("Description", DataType::STRING, 
			"This model is a piecewise-continuous radial profile, "
			"where each segment is described by a power-law of the form:\n\n"
			"   /                  \\                1              \n"
			" p | r  <=  r  < r    |   =   -------------------      \n"
			"   \\  i           i+1 /                         bi    \n"
			"                              /           ai  \\       \n"
			"                              |  fac  +  r    |        \n"
			"                              \\               /        \n\n"

			"Thus fac = 0, bi = 0 corresponds to a simple power law, and fac = 1, ai = 2 corresponds to a beta model. "
			"The number of inflection points is specified using the parameter 'n'.  Note that 'n = 1' means that there "
			"will be two segments: one from r = 0 --> r1 and one from r = r1 --> infinity. \n\n"
			"Only the inflection points "
			"can be specified; the origin is always r=0.  The zeroth parameters a0 and b0 therefore always refer to the "
			"exponents of the segment from r=0 to r1. This profile shape is multiplied by the normalization "
			"(Sradio or Sxray, etc); the shape is therefore always normalized to 1 at r=0.  This and the requirement that it "
			"be integrable mean that if fac=0, the profile is held constant at its r1 value if r < r1.");
}

/**.......................................................................
 * Destructor.
 */
PowerlawProfile::~PowerlawProfile() {}

/**.......................................................................
 * Check setup of this model
 */
void PowerlawProfile::checkSetup()
{
  //------------------------------------------------------------
  // Now check that exponents have been specified for this model
  //------------------------------------------------------------

  std::ostringstream os;

  if(fac_.val_ < 0.0) {
    ThrowSimpleColorError("Additive factor must be >= 0.0", "red");
  } else if(!(radius_.size() > 1 || fac_.val_ > 0.0)) {
    ThrowSimpleColorError("Additive factor must be > 0 if no inflection points are specified", "red");
  }

  for(unsigned i=0; i < radius_.size(); i++) {

    os.str("");
    os << "r" << i;

    if(i > 0)
      checkVar(os.str());

    //------------------------------------------------------------
    // Require all exponents if fac > 0, or exponents for i > 0 if
    // not, since in this case index 0 exponents are never used
    //------------------------------------------------------------

    if(fac_.val_ > 0.0 || (!(fac_.val_ > 0.0) && i > 0)) {

      os.str("");
      os << "a" << i;
      
      checkVar(os.str());
      
      os.str("");
      os << "b" << i;
      
      checkVar(os.str());
    }

    //------------------------------------------------------------
    // Ensure that our radius array is monotonic
    //------------------------------------------------------------

    if(i > 0 && radius_[i].val_ <= radius_[i-1].val_) {
      ThrowSimpleColorError("Inflection points must be specified in monotonically increasing order", "red");
    }
  }

  //------------------------------------------------------------
  // Precompute the relevant norms
  //------------------------------------------------------------

  norm_ = profileFn(radius_[0].val_, alpha_[0].val_, beta_[0].val_);

  unsigned iStart=1;
  piecewiseNorm_[0] = 1.0;

  if(!(fac_.val_ > 0.0)) {
    iStart=2;
    piecewiseNorm_[1] = 1.0;
  }

  for(unsigned i=iStart; i < radius_.size(); i++) {
    piecewiseNorm_[i] = piecewiseNorm_[i-1] * 
      profileFn(radius_[i].val_, alpha_[i-1].val_, beta_[i-1].val_) / profileFn(radius_[i].val_, alpha_[i].val_, beta_[i].val_);
  }

  //------------------------------------------------------------
  // Finally, perform any additional base-class checks
  //------------------------------------------------------------

  GenericRadiallySymmetric3DModel::checkSetup();
}

double PowerlawProfile::radialModel(unsigned type, double x, void* params)
{
  switch (type) {
  case DataSetType::DATASET_RADIO:
    return radialRadioModelEml(x, params);
    break;
  case DataSetType::DATASET_XRAY_IMAGE:
    return radialXrayModel(x, params);
    break;
  case DataSetType::DATASET_GENERIC:
    return radialRadioModelEml(x, params);
    break;
  default:
    ThrowColorError("Unsupported dataset type: " << type, "red");
    return 0.0;
    break;
  }
}

double PowerlawProfile::radialRadioModelEml(double x, void* params)
{
  unsigned iLow=0, iHigh=0;
  binSearchForRadius(x, iLow, iHigh);
  return piecewiseNorm_[iLow] * profileFn(x, alpha_[iLow].val_, beta_[iLow].val_);
}

double PowerlawProfile::radialXrayModel(double x, void* params)
{
  ThrowSimpleColorError("No radialXrayModel() function has been defined for this model", "red");
  return 0.0;
}

/**.......................................................................
 * Return true if all parameters affecting the shape (note: not the
 * scale) of the model are fixed.  However, if thetaCore has no prior
 * specified, this means that we can't precompute its possible range,
 * so return false so that the interpolation grid is recomputed on
 * each trial of thetaCore
 */
bool PowerlawProfile::shapeParametersAreFixed()
{
  if(thetaCore_.prior().getType() == Distribution::DIST_UNSPEC) 
    return false;

  for(unsigned i=0; i< alpha_.size(); i++)
    if(alpha_[i].isVariable() || radius_[i].isVariable() || beta_[i].isVariable())
      return false;

  return true;
}

void PowerlawProfile::setParameter(std::string name, std::string val, std::string units)
{
  String nameStr(name);

  //------------------------------------------------------------
  // Always call the underlying PM method:
  //------------------------------------------------------------

  ParameterManager::setParameter(name, val, units);

  if(name == "n") {
    unsigned n = getUintVal("n");
    setNumberOfSegments(n);
  }
}

/**.......................................................................
 * Change the number of segments in the model
 */
void PowerlawProfile::setNumberOfSegments(unsigned n)
{
  int nPrev = piecewiseNorm_.size();
  
  piecewiseNorm_.resize(n+1);
  radius_.resize(n+1);
  alpha_.resize(n+1);
  beta_.resize(n+1);

  std::ostringstream os, expos;

  for(unsigned i=0; i < piecewiseNorm_.size(); i++) {

    //------------------------------------------------------------
    // rad0 is not available to the user -- this will always be set to
    // zero
    //------------------------------------------------------------

    if(i == 0) {
      radius_[0].val_ = 0.0;
    } else {
      os.str("");
      expos.str("");
      
      os << "r" << i;
      expos << "The " << i << String::numericalSuffix(i) << " inflection radius (in units of theteCore)";
      
      if(nPrev == 0 || i > nPrev-1) {
	addComponent(radius_[i]);
	addComponentName(radius_[i], os.str(), expos.str());
      }
    }

    os.str("");
    expos.str("");

    os << "a" << i;
    if(i==0) {
      expos << "The 'a' power-law exponent outward of radius zero";
    } else {
      expos << "The " << i << String::numericalSuffix(i) << " 'a' power-law exponent outward of the corresponding radius";
    }

    if(nPrev == 0 || i > nPrev-1) {
      addComponent(alpha_[i]);
      addComponentName(alpha_[i], os.str(), expos.str());
    }

    os.str("");
    expos.str("");

    os << "b" << i;
    if(i==0) {
      expos << "The 'b' power-law exponent outward of radius zero";
    } else {
      expos << "The " << i << String::numericalSuffix(i) << " 'b' power-law exponent outward of the corresponding radius";
    }

    if(nPrev == 0 || i > nPrev-1) {
      addComponent(beta_[i]);
      addComponentName(beta_[i], os.str(), expos.str());
    }

    //------------------------------------------------------------
    // Initialize these components to fixed
    //------------------------------------------------------------

    radius_[i].isVariable() = false;
    radius_[i].allowUnitless(true);
    radius_[i].val_ = 0.0;

    alpha_[i].isVariable() = false;
    alpha_[i].allowUnitless(true);
    alpha_[i].val_ = 0.0;

    beta_[i].isVariable() = false;
    beta_[i].allowUnitless(true);
    beta_[i].setVal(1.0, "");
  }
}

void PowerlawProfile::binSearchForRadius(double rad, unsigned& iLow, 
					 unsigned& iHigh)
{
  double yHigh, yLow, yMid;
  unsigned n = radius_.size();
  iLow=0, iHigh=radius_.size()-1;
  unsigned mid;

  if(rad > radius_[n-1].val_) {
    iLow = iHigh = n-1;
    return;
  }

  if(rad < radius_[0].val_) {
    iLow = 0;
    iHigh = n > 1 ? 1 : 0;
    return;
  }
    
  while(iHigh - iLow > 1) {

    mid = (iHigh+iLow)/2;

    yHigh = radius_[iHigh].val_;
    yLow  = radius_[iLow].val_;
    yMid  = radius_[mid].val_;

    if(rad > yMid) {
      iLow = mid;
    } else if(rad < yMid) {
      iHigh = mid;
    } else {
      iLow = iHigh = mid;
    }
  }
}

double PowerlawProfile::profileFn(double x, double alpha, double beta)
{
  //------------------------------------------------------------
  // If the function as-defined is not integrable, fix the value to
  // the first inflection point for r < r1
  //------------------------------------------------------------

  if(!(fac_.val_ > 0.0)) {

    double r1 = radius_[1].val_;
    double a1 = alpha_[1].val_;
    double b1 = beta_[1].val_;

    return ((x >= r1) ? 1.0/pow(fac_.val_ + pow(x, alpha), beta) : 1.0/pow(fac_.val_ + pow(r1, a1), b1))/norm_;

    //------------------------------------------------------------
    // Else evaluate as requested
    //------------------------------------------------------------
    
  } else {
    return (1.0/pow(fac_.val_ + pow(x, alpha), beta))/norm_;
  }
}

