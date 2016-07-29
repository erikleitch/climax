#include "gcp/util/Constants.h"
#include "gcp/util/Cosmology.h"
#include "gcp/util/Density.h"
#include "gcp/util/HubbleConstant.h"

#include <cmath>

using namespace std;

using namespace gcp::util;

const double Cosmology::eps_ = 1e-12;

/**.......................................................................
 * Constructor.
 */
Cosmology::Cosmology() 
{
  addParameter("z",      DataType::DOUBLE);
  addParameter("h70",    DataType::DOUBLE);
  addParameter("h100",   DataType::DOUBLE);
  addParameter("H0",     DataType::DOUBLE);
  addParameter("omegaM", DataType::DOUBLE);
  addParameter("omegaL", DataType::DOUBLE);

  setParameter("h70",    "1.0");
  setParameter("omegaM", "0.3");
  setParameter("omegaL", "0.7");

  initializeGslMembers();
}

/**.......................................................................
 * Initialize members needed for numerical integration
 */
void Cosmology::initializeGslMembers()
{
  setGslLimit(100);
  gslEpsAbs_ = 0;
  gslEpsRel_ = 1e-7;

  gslWork_ = 0;
  gslWork_ = gsl_integration_workspace_alloc(gslLimit_);

  if(!gslWork_) {
    ThrowError("Unable to allocate GSL workspace");
  }

  gslComovingFn_.function = &Cosmology::evaluateComovingDistanceKernel;
  gslComovingFn_.params   = (void*)this;

  gslLightTravelFn_.function = &Cosmology::evaluateLightTravelDistanceKernel;
  gslLightTravelFn_.params   = (void*)this;
}

/**.......................................................................
 * Set the limit for GSL integration
 */
void Cosmology::setGslLimit(size_t limit)
{
  gslLimit_ = limit;
}

/**.......................................................................
 * Destructor.
 */
Cosmology::~Cosmology() 
{
  if(gslWork_) {
    gsl_integration_workspace_free(gslWork_);
    gslWork_ = 0;
  }
}

//------------------------------------------------------------
// Set methods
//------------------------------------------------------------

void Cosmology::setH0(HubbleConstant H0)
{
  std::ostringstream os;
  os << H0.kmPerSecPerMpc();
  setParameter("H0", os.str(), "km/s/Mpc");
}

void Cosmology::setOmegaM(double omegaM)
{
  std::ostringstream os;
  os << omegaM;
  setParameter("omegaM", os.str());
}

void Cosmology::setOmegaL(double omegaL)
{
  std::ostringstream os;
  os << omegaL;
  setParameter("omegaL", os.str());
}

void Cosmology::setRedshift(double z)
{
  std::ostringstream os;
  os << z;
  setParameter("z", os.str());
}

//------------------------------------------------------------
// Query methods
//------------------------------------------------------------

HubbleConstant Cosmology::H0()
{
  return H0_;
}

double Cosmology::omegaL()
{
  checkParameter("omegaL");
  return getDoubleVal("omegaL");
}

double Cosmology::omegaM() 
{
  checkParameter("omegaM");
  return getDoubleVal("omegaM");
}

double Cosmology::omegaK() 
{
  return 1.0 - (omegaM() + omegaL());
}

double Cosmology::redshift()
{
  checkParameter("z");
  return getDoubleVal("z");
}

/**.......................................................................
 * If this is a flat cosmology, then omegaK should be zero:
 * (omegaK = 1.0 - omegaM - omegaL)
 */
bool Cosmology::isFlat()
{
  return fabs(omegaK()) < eps_;
}

/**.......................................................................
 * Dimensionless E(z) for given cosmology
 */
double Cosmology::E()
{
  checkParameter("z");
  return E(z_);
}

double Cosmology::E(double z)
{
  return sqrt(E2(z));
}

/**.......................................................................
 * Dimensionless E^2(z) for given cosmology
 */
double Cosmology::E2()
{
  checkParameter("z");
  return E2(z_);
}

double Cosmology::E2(double z)
{
  double z1 = 1 + z;

  // Valid for all cosmologies
  //  COUT("OM = " << omegaM() << " OK = " << omegaK() << " OL = " << omegaL());

  return ((omegaM() * z1) + omegaK() ) * z1 * z1 + omegaL();
}

/**.......................................................................
 * Hubble constant at our current redshift
 */
HubbleConstant Cosmology::H()
{
  checkParameter("z");

  return H(z_);
}

/**.......................................................................
 * Hubble constant as a function of redshift
 */
HubbleConstant Cosmology::H(double z)
{
  checkParametersOr("H0", "h70", "h100");

  return H0_ * E(z);
}

double Cosmology::dimensionlessHubbleConstant()
{
  return H() / H(0);
}

double Cosmology::h()
{
  return dimensionlessHubbleConstant();
}

double Cosmology::dimensionlessHubbleConstant(double z)
{
  return H(z) / H(0);
}

double Cosmology::h(double z)
{
  return dimensionlessHubbleConstant(z);
}

/**.......................................................................
 * Return the Hubble distance, dH = c/H0
 */
Length Cosmology::hubbleDistance()
{
  checkParametersOr("H0", "h70", "h100");

  return Constants::lightSpeed_ / H0_;
}

/**.......................................................................
 * Comoving distance:
 *
 *              z
 *            /     dz
 *  d   = d   |    ----
 *   C     H  /    E(z)
 *              0
 */ 
Length Cosmology::comovingDistance()
{
  checkParameter("z");
  return comovingDistance(z_);
}

Length Cosmology::comovingDistance(double z)
{
  double result, abserr;

  gsl_integration_qags(&gslComovingFn_, 0, z, 
		       gslEpsAbs_, gslEpsRel_, gslLimit_, 
		       gslWork_, &result, &abserr);

  return hubbleDistance() * result;
}

/**.......................................................................
 * Return the transverse comoving distance
 */
Length Cosmology::transverseComovingDistance()
{
  checkParameter("z");
  return transverseComovingDistance(z_);
}

Length Cosmology::transverseComovingDistance(double z)
{
  double ok = omegaK();
  Length dc = comovingDistance(z);
  Length dh = hubbleDistance();

  if(isFlat()) {
    return dc;
  } else if(ok > 0.0) {
    double sqok = sqrt(ok);
    return dh * (sinh(sqok * (dc/dh)) / sqok);
  } else {
    double sqok = sqrt(-ok);
    return dh * (sin(sqok * (dc/dh)) / sqok);
  }
}

/**.......................................................................
 * Return the angular diameter distance
 */
Length Cosmology::angularDiameterDistance()
{
  checkParameter("z");

  return angularDiameterDistance(z_);
}

Length Cosmology::angularDiameterDistance(double z)
{
  return transverseComovingDistance(z) / (1.0 + z);
}

/**.......................................................................
 * Return the luminosity distance
 */
Length Cosmology::luminosityDistance()
{
  checkParameter("z");

  return luminosityDistance(z_);
}

Length Cosmology::luminosityDistance(double z)
{
  return transverseComovingDistance(z) * (1.0 + z);
}

/**.......................................................................
 * Return the light travel distance
 */
Length Cosmology::lightTravelDistance()
{
  checkParameter("z");

  return lightTravelDistance(z_);
}

Length Cosmology::lightTravelDistance(double z)
{
  double result, abserr;

  gsl_integration_qags(&gslLightTravelFn_, 0, z, 
		       gslEpsAbs_, gslEpsRel_, gslLimit_, 
		       gslWork_, &result, &abserr);

  return hubbleDistance() * result;
}

/**.......................................................................
 * Evaluate the kernel needed for calculating the comoving distance
 */
double Cosmology::evaluateComovingDistanceKernel(double z, void* params)
{
  Cosmology* cosmo = (Cosmology*) params;

  // And evaluate the kernel:

  return 1.0/cosmo->E(z);
}

/**.......................................................................
 * Evaluate the kernel needed for calculating the light-travel distance
 */
double Cosmology::evaluateLightTravelDistanceKernel(double z, void* params)
{
  Cosmology* cosmo = (Cosmology*) params;

  // And evaluate the kernel:

  return 1.0/((1.0 + z) * cosmo->E(z));
}

/**.......................................................................
 * Method to set a parameter for this object
 */
void Cosmology::setParameter(std::string name, std::string val, std::string units, bool external)
{
  String nameStr(name);

  // Always call the underlying PM method:
  
  ParameterManager::setParameter(name, val, units, external);
  
  if(name == "H0") {
    H0_.setVal(getDoubleVal("H0"), getParameter("H0", true)->units_);
  } else if(name == "h100") {
    H0_.setKmPerSecPerMpc(getDoubleVal("h100") * 100);
  } else if(name == "h70") {
    H0_.setKmPerSecPerMpc(getDoubleVal("h70") * 70);
  } else if(name == "z") {
    z_ = getDoubleVal("z");
  } else if(name == "omegaL") {
    omegaL_ = getDoubleVal("omegaL");
  } else if(name == "omegaM") {
    omegaM_ = getDoubleVal("omegaM");
  }
}

void Cosmology::checkParametersAnd(std::string par1, std::string par2)
{
  if(!(getParameter(par1, false)->data_.hasValue() && getParameter(par2, false)->data_.hasValue())) {
    ThrowError("No value has been set for " << par1 << " and " << par2 << " " << getParameter(par1, false)->owner_->name_ << "." << par1);
  }
}

void Cosmology::checkParametersOr(std::string par1, std::string par2)
{
  if(!(getParameter(par1, false)->data_.hasValue() || getParameter(par2, false)->data_.hasValue())) {
    ThrowError("No value has been set for " << par1 << " or " << par2);
  }
}

void Cosmology::checkParametersOr(std::string par1, std::string par2, std::string par3)
{
  if(!(getParameter(par1, false)->data_.hasValue() 
       || getParameter(par2, false)->data_.hasValue()
       || getParameter(par3, false)->data_.hasValue())) {
    ThrowError("No value has been set for " << par1 << " or " << par2 << " or " << par3);
  }
}

void Cosmology::checkParameter(std::string par)
{
  if(!(getParameter(par, false)->data_.hasValue())) {
    ThrowError("No value has been set for " << par << " " << getParameter(par, false)->owner_->name_ << "." << par);
  }
}

Density Cosmology::criticalDensity()
{
  return criticalDensity(z_);
}

Density Cosmology::criticalDensity(double z)
{
  HubbleConstant Hz = H(z);
  Density rho;

  rho.setGPerCm3(3.0/(8*M_PI*Constants::gravitationalConstantCgs_) * Hz.inverseSeconds() * Hz.inverseSeconds());

  return rho;
}
