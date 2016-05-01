#include "gcp/util/SpectralType.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
SpectralType::SpectralType() 
{
  initializeMaps();
}

/**.......................................................................
 * Destructor.
 */
SpectralType::~SpectralType() {}

/**.......................................................................
 * Define what spectral type this object encapsulates
 */
void SpectralType::initializeMaps()
{
  idToNameMap_[SPEC_NONE]  =  "none";
  idToNameMap_[SPEC_ALPHA] = "alpha";
  idToNameMap_[SPEC_SZ]    =    "sz";
  idToNameMap_[SPEC_ITOH]  =  "itoh";

  nameToIdMap_["none"]  = SPEC_NONE;
  nameToIdMap_["alpha"] = SPEC_ALPHA;
  nameToIdMap_["sz"]    = SPEC_SZ;
  nameToIdMap_["itoh"]  = SPEC_ITOH;

  explMap_[SPEC_NONE]  = "no spectral shape";
  explMap_[SPEC_ALPHA] = "spectral index";
  explMap_[SPEC_SZ]    = "SZ spectral shape";
  explMap_[SPEC_ITOH]  = "Itoh SZ spectral shape";
}
