#include "gcp/util/HubbleConstant.h"

#include <iomanip>

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
HubbleConstant::HubbleConstant() 
{
  addConversion("km/s/Mpc",  1.0);
}

/**.......................................................................
 * Destructor.
 */
HubbleConstant::~HubbleConstant() {}

HubbleConstant HubbleConstant::operator*(double fac)
{
  HubbleConstant hc;
  hc.setKmPerSecPerMpc(val_ * fac);
  return hc;
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, const HubbleConstant& hc)
{
  return operator<<(os, (HubbleConstant&)hc);
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, HubbleConstant& hc)
{
  os << std::setw(14) << hc.kmPerSecPerMpc() << " km/s/Mpc";
  return os;
}

double HubbleConstant::operator/(const HubbleConstant& hc)
{
  return operator/((HubbleConstant&) hc);
}

double HubbleConstant::operator/(HubbleConstant& hc)
{
  return kmPerSecPerMpc() / hc.kmPerSecPerMpc();
}

void HubbleConstant::operator=(const HubbleConstant& var)
{
  operator=((HubbleConstant&) var);
}

void HubbleConstant::operator=(HubbleConstant& var)
{
  val_ = var.val_;
}
