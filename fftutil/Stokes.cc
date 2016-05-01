#include "gcp/fftutil/Stokes.h"

/**.......................................................................
 * Write the contents of this object to an ostream
 */
std::ostream& 
gcp::util::operator<<(std::ostream& os, const Stokes::Param& param)
{
  switch(param) {
  case Stokes::STOKES_I:
    os << "STOKES I";
    break;
  case Stokes::STOKES_Q:
    os << "STOKES Q";
    break;
  case Stokes::STOKES_U:
    os << "STOKES U";
    break;
  case Stokes::STOKES_T:
    os << "STOKES T";
    break;
  case Stokes::STOKES_E:
    os << "STOKES E";
    break;
  case Stokes::STOKES_B:
    os << "STOKES B";
    break;
  default:
    os << "STOKES UNKNOWN";
    break;
  }
}

