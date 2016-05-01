#include "gcp/util/Exception.h"
#include "gcp/util/Probability.h"

#include <cmath>

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Probability::Probability() 
{
  lnProb_ = 0.0;
}

/**.......................................................................
 * Destructor.
 */
Probability::~Probability() {}

void Probability::setValue(double value)
{
  lnProb_ = log(value);
}

void Probability::setLnValue(double lnValue)
{
  lnProb_ = lnValue;
}

double Probability::value()
{
  if(isnan(lnProb_))
    ThrowError("Probability is NaN");

  return exp(lnProb_);
}

double Probability::lnValue()
{
  return lnProb_;
}

void Probability::operator*=(const Probability& prob)
{
  operator*=((Probability&) prob);
}

void Probability::operator*=(Probability& prob)
{
  lnProb_ += prob.lnProb_;
}

Probability Probability::operator*(const Probability& prob)
{
  return operator*((Probability&) prob);
}

Probability Probability::operator*(Probability& prob)
{
  Probability ret;
  ret.lnProb_ = lnProb_ + prob.lnProb_;
  return ret;
}

void Probability::operator/=(const Probability& prob)
{
  operator/=((Probability&) prob);
}

void Probability::operator/=(Probability& prob)
{
  lnProb_ -= prob.lnProb_;
}

Probability Probability::operator/(const Probability& prob)
{
  return operator/((Probability&) prob);
}

Probability Probability::operator/(Probability& prob)
{
  Probability ret;
  ret.lnProb_ = lnProb_ - prob.lnProb_;
  return ret;
}

bool Probability::operator<(const Probability& prob)
{
  return operator<((Probability&) prob);
}

bool Probability::operator>(const Probability& prob)
{
  return operator>((Probability&) prob);
}

bool Probability::operator<(double val)
{
  double lnVal = log(val);
  return lessThan(lnProb_, lnVal);
}

bool Probability::operator>(double val)
{
  double lnVal = log(val);
  return greaterThan(lnProb_, lnVal);
}

/**.......................................................................
 * Check if one probability is less than another.  If either
 * probability is invalid, we return false
 */
bool Probability::operator<(Probability& prob)
{
  return lessThan(lnProb_, prob.lnProb_);
}

/**.......................................................................
 * Check is one probability is greater than another.  If either
 * probability is invalid, we return false
 */
bool Probability::operator>(Probability& prob)
{
  return greaterThan(lnProb_, prob.lnProb_);
}

bool Probability::lessThan(double val1, double val2)
{
  if(isValid(val1) && isValid(val2)) {
    return val1 < val2;
  } else {
    return false;
  }
}

bool Probability::greaterThan(double val1, double val2)
{
  if(isValid(val1) && isValid(val2)) {
    return val1 > val2;
  } else {
    return false;
  }
}

/**.......................................................................
 * Return true if this contains a valid probability.  
 *
 * Note that above operators can cause:
 *
 * lnProb_ = -inf (prob = zero)
 * lnProb_ = +inf (result of dividing prob/0.0)
 * lnProb_ =  nan (result of dividing 0.0/0.0)
 *
 * I consider only the last to be invalid, since probability
 * comparisons with other values will return valid results.
 */
bool Probability::isValid()
{
  return isValid(lnProb_);
}

bool Probability::isValid(double val)
{
  return !isnan(val);
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, const Probability& probability)
{
  return operator<<(os, (Probability&)probability);
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, Probability& probability)
{
  os << probability.value();
  return os;
}
