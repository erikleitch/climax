#include "gcp/models/MirrorModel.h"

#include "gcp/util/DataType.h"

using namespace std;

using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
MirrorModel::MirrorModel() 
{
  addParameter("nf", DataType::UINT);
  addParameter("nr", DataType::UINT);

  addComponent(xTrans_);
  addComponent(yTrans_);
  addComponent(zTrans_);

  xTrans_.allowUnitless(true);
  yTrans_.allowUnitless(true);
  zTrans_.allowUnitless(true);

  addComponentName(xTrans_, "xtrans", "x-translation");
  addComponentName(yTrans_, "ytrans", "y-translation");
  addComponentName(zTrans_, "ztrans", "z-translation");

  addComponent(xCr_);
  addComponent(yCr_);
  addComponent(zCr_);

  xCr_.allowUnitless(true);
  yCr_.allowUnitless(true);
  zCr_.allowUnitless(true);

  addComponentName(xCr_, "xcr", "x-translation");
  addComponentName(yCr_, "ycr", "y-translation");
  addComponentName(zCr_, "zcr", "z-translation");

  addComponent(xRot_);
  addComponent(yRot_);
  addComponent(zRot_);

  addComponentName(xRot_, "xrot", "x-rotation");
  addComponentName(yRot_, "yrot", "y-rotation");
  addComponentName(zRot_, "zrot", "z-rotation");
  
  initializeComponentsToFixed();
}

MirrorModel::~MirrorModel() 
{}

void MirrorModel::setParameter(std::string name, std::string val, std::string units, bool external)
{
  String nameStr(name);

  //------------------------------------------------------------
  // Always call the underlying PM method:
  //------------------------------------------------------------

  ParameterManager::setParameter(name, val, units, external);

  if(name == "nr") {
    unsigned nr = getUintVal("nr");
    setNumberOfReferencePoints(nr);
  }

  if(name == "nf") {
    unsigned nf = getUintVal("nf");
    setNumberOfFiducialPoints(nf);
  }
}

void MirrorModel::setNumberOfFiducialPoints(unsigned n)
{
  xF_.resize(n);
  yF_.resize(n);
  zF_.resize(n);

  std::ostringstream os, expos;

  for(unsigned i=0; i < n; i++) {
    addComponent(xF_[i]);
    addComponent(yF_[i]);
    addComponent(zF_[i]);

    xF_[i].allowUnitless(true);
    yF_[i].allowUnitless(true);
    zF_[i].allowUnitless(true);

    os.str("");
    expos.str("");

    os << "xf" << i;
    expos << "The " << i << "th fiducial x-coordinate";
    addComponentName(xF_[i], os.str(), expos.str());

    os.str("");
    expos.str("");

    os << "yf" << i;
    expos << "The " << i << "th fiducial y-coordinate";
    addComponentName(yF_[i], os.str(), expos.str());

    os.str("");
    expos.str("");

    os << "zf" << i;
    expos << "The " << i << "th fiducial z-coordinate";
    addComponentName(zF_[i], os.str(), expos.str());

    // Initialize this component to fixed

    xF_[i].isVariable() = false;
    xF_[i].val_ = 0.0;

    yF_[i].isVariable() = false;
    yF_[i].val_ = 0.0;

    zF_[i].isVariable() = false;
    zF_[i].val_ = 0.0;
  }
}

void MirrorModel::setNumberOfReferencePoints(unsigned n)
{
  xR_.resize(n);
  yR_.resize(n);
  zR_.resize(n);

  // And create arrays for the model positions

  xRm_.resize(n);
  yRm_.resize(n);
  zRm_.resize(n);

  std::ostringstream os, expos;

  for(unsigned i=0; i < n; i++) {
    addComponent(xR_[i]);
    addComponent(yR_[i]);
    addComponent(zR_[i]);

    xR_[i].allowUnitless(true);
    yR_[i].allowUnitless(true);
    zR_[i].allowUnitless(true);

    os.str("");
    expos.str("");

    os << "xr" << i;
    expos << "The " << i << "th reference x-coordinate";
    addComponentName(xR_[i], os.str(), expos.str());

    os.str("");
    expos.str("");

    os << "yr" << i;
    expos << "The " << i << "th fiducial y-coordinate";
    addComponentName(yR_[i], os.str(), expos.str());

    os.str("");
    expos.str("");

    os << "zr" << i;
    expos << "The " << i << "th fiducial z-coordinate";
    addComponentName(zR_[i], os.str(), expos.str());

    // Initialize this component to fixed

    xR_[i].isVariable() = false;
    xR_[i].val_ = 0.0;

    yR_[i].isVariable() = false;
    yR_[i].val_ = 0.0;

    zR_[i].isVariable() = false;
    zR_[i].val_ = 0.0;
  }
}

void MirrorModel::fillArray(unsigned type, Angle& axisUnits, 
			    std::valarray<double>& x, std::valarray<double>& y, 
			    std::valarray<double>& d, 
			    void* params)
{
  // Calculate the positions of the reference points, given the
  // current rotation/translation about the CR position

  double cx = cos(xRot_.radians());
  double sx = sin(xRot_.radians());

  double cy = cos(yRot_.radians());
  double sy = sin(yRot_.radians());

  double cz = cos(zRot_.radians());
  double sz = sin(zRot_.radians());

  // Current rotation matrix

  Matrix<double> m(3,3);

  m[0][0] =             cy*cz; m[0][1] =            -cy*sz; m[0][2] =     sy;
  m[1][0] =  sx*sy*cz + cx*sz; m[1][1] = -sx*sy*sz + cx*cz; m[1][2] = -sx*cy;
  m[2][0] = -cx*sy*cz + sx*sz; m[2][1] =  cx*sy*sz + sx*cz; m[2][2] =  cx*cy;

  //  COUT("M = " << std::endl << m);

  Vector<double> xp(3);
  Vector<double> xm(3);
  Vector<double> cr(3);
  Vector<double> crd(3);
  Vector<double> trans(3);

  // Get the coordinates of the CR

  cr[0] = xCr_.val_;
  cr[1] = yCr_.val_;
  cr[2] = zCr_.val_;

  //  COUT("CR = " << cr);

  // Get the delta translation

  crd[0] = xTrans_.val_;
  crd[1] = yTrans_.val_;
  crd[2] = zTrans_.val_;

  //  COUT("CRD = " << crd);
  
  for(unsigned i=0; i < xR_.size(); i++) {
    
    xp[0] = xR_[i].val_;
    xp[1] = yR_[i].val_;
    xp[2] = zR_[i].val_;

    //    COUT("XP" << i << " = " << std::endl << xp);

    // Get the coordinates of this point relative to the CR

    xp = xp - cr;

    //    COUT("XP" << i << " = " << std::endl << xp);

    // Now rotate it by the current rotation matrix

    xm = m * xp;

    //    COUT("XM" << i << " = " << std::endl << xp);

    // And translate it back by the CR + deltas

    xm = xm + cr;
    xm = xm + crd;

    //    COUT("XM" << i << " = " << std::endl << xm);

    xRm_[i].val_ = xm[0];
    yRm_[i].val_ = xm[1];
    zRm_[i].val_ = xm[2];
  }

  // Now fill the array with model values

  Length dtmp;
  for(unsigned i=0; i < x.size(); i++) {
    unsigned iF = (unsigned)x[i];
    unsigned iR = (unsigned)y[i];

    double dx = xF_[iF].val_ - xRm_[iR].val_;
    double dy = yF_[iF].val_ - yRm_[iR].val_;
    double dz = zF_[iF].val_ - zRm_[iR].val_;

    // Set the model distance to be the coordinate separation, minus
    // the ball diameter

    dtmp.setMillimeters(sqrt(dx*dx + dy*dy + dz*dz));
    dtmp.setInches(dtmp.inches() - 0.75);

    d[i] = dtmp.inches() * (1.0 + 21.0/1e6 * 50);

    //    COUT("Computing chisq for iF = " << iF << " iR = " << iR << " d = " << d[i]);
  } 
}
