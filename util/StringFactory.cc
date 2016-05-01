#include "gcp/util/Exception.h"
#include "gcp/util/StringFactory.h"

using namespace std;

using namespace gcp::util;

std::vector<std::string> StringFactory::strings_;

/**.......................................................................
 * Constructor.
 */
StringFactory::StringFactory() {}

/**.......................................................................
 * Destructor.
 */
StringFactory::~StringFactory()
{
  COUT("Inside StringFactory destructor");
}

char* StringFactory::getString(std::string str)
{
  strings_.push_back(str);
  return (char*)strings_[strings_.size()-1].c_str();
}


