#include "gcp/util/Exception.h"
#include "gcp/util/Single.h"

using namespace std;

using namespace gcp::util;

//bool Single::init_ = Single::createInstance();
int Single::count_ = 0;

/**.......................................................................
 * Constructor.
 */
Single::Single() 
{
  count_ = 1;
}

/**.......................................................................
 * Destructor.
 */
Single::~Single() 
{
  COUT("Inside destructor");
}

bool Single::createInstance()
{
  static bool init = true;
  if(init) {
    static Single s;
    init = false;
  }
  return true;
}

void Single::printHello()
{
  COUT("Count = " << count_);
  ++count_;
}
