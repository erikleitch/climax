#include "gcp/models/Polynomial1D.h"

#include "gcp/util/DataType.h"
#include "gcp/util/Exception.h"
#include "gcp/util/String.h"

using namespace std;

using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Polynomial1D::Polynomial1D() 
{
  addParameter("n", DataType::UINT);
}

/**.......................................................................
 * Destructor.
 */
Polynomial1D::~Polynomial1D() {}

void Polynomial1D::setParameter(std::string name, std::string val, std::string units, bool external)
{
  String nameStr(name);

  //------------------------------------------------------------
  // Always call the underlying PM method:
  //------------------------------------------------------------

  Model::setParameter(name, val, units, external);

  if(name == "n") {
    unsigned n = getUintVal("n");
    setNumberOfCoefficients(n);
  }
}

void Polynomial1D::setNumberOfCoefficients(unsigned n)
{
  coeffs_.resize(n);
  std::ostringstream os, expos;

  for(unsigned i=0; i < n; i++) {
    addComponent(coeffs_[i]);

    os.str("");
    expos.str("");

    os << "a" << i;
    expos << "The " << i << "th coefficient of the polynomial";

    addComponentName(coeffs_[i], os.str(), expos.str());

    // Initialize this component to fixed

    coeffs_[i].isVariable() = false;
    coeffs_[i].allowUnitless(true);
    coeffs_[i].val_ = 0.0;
  }
}

double Polynomial1D::eval(double x)
{
  double val = 0.0;
  for(int i=coeffs_.size()-1; i >=0; i--) {
    val = val*x + coeffs_[i].val_;
  }

  return val;
}
