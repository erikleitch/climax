#include "gcp/util/Exception.h"
#include "gcp/util/FileHandler.h"

#include <sstream>
#include <fstream>

using namespace std;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
FileHandler::FileHandler() {}

/**.......................................................................
 * Destructor.
 */
FileHandler::~FileHandler() {}

bool FileHandler::fileExists(std::string name)
{
  std::ifstream fin;
  fin.open(name.c_str(), ios::in);

  if(fin) {
    fin.close();
    return true;
  } else {
    return false;
  }
}
