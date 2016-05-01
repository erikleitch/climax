#include "gcp/util/Declination.h"
#include "gcp/util/Exception.h"

using namespace std;
using namespace gcp::util;

const double Declination::arcSecPerRad_ = 206265;
const double Declination::arcMinPerRad_ = 206265.0 / 60;
