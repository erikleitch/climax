#include "gcp/util/EnumeratedVariate.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
EnumeratedVariate::EnumeratedVariate() 
{
  isEnumerated_ = true;
}

/**.......................................................................
 * Destructor.
 */
EnumeratedVariate::~EnumeratedVariate() {}

std::string EnumeratedVariate::getStringVal() 
{
  if(idToNameMap_.find(type_) == idToNameMap_.end()) {

    std::ostringstream os;

    os << std::endl << "Unrecognized type: " << type_ << " for variate " << name_ << std::endl
       << "Should be one of the following: " << std::endl << std::endl;

    listValidOptions(os);

    ThrowColorError(os.str(), "red");
  }

  return idToNameMap_[type_];
}

/**.......................................................................
 * Define what spectral type this object encapsulates
 */
void EnumeratedVariate::setVal(std::string val)
{
  if(nameToIdMap_.find(val) == nameToIdMap_.end()) {

    std::ostringstream os;
    
    os << std::endl << "Unrecognized value: '" << val << "' for variate " << name_ << std::endl
       << "Should be one of the following: " << std::endl << std::endl;
    
    listValidOptions(os);
    
    ThrowColorError(os.str(), "red");
  }
  
  type_ = nameToIdMap_[val];
}

void EnumeratedVariate::listValidOptions(std::ostringstream& os)
{
  unsigned nameMaxLen=0, explMaxLen=0;
  for(std::map<unsigned, std::string>::iterator iter = idToNameMap_.begin(); iter != idToNameMap_.end(); iter++) {
    std::string name = iter->second;
    std::string expl = explMap_[iter->first];
    nameMaxLen = name.size() > nameMaxLen ? name.size() : nameMaxLen;
    explMaxLen = expl.size() > explMaxLen ? expl.size() : explMaxLen;
  }

  for(std::map<unsigned, std::string>::iterator iter = idToNameMap_.begin(); iter != idToNameMap_.end(); iter++)
    os << "  " << std::setw(nameMaxLen+1) << std::left << iter->second << " -- " 
       << std::setw(explMaxLen+1) << std::left << explMap_[iter->first] << std::endl;
}
