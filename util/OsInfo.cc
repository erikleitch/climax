#include "gcp/util/Exception.h"
#include "gcp/util/OsInfo.h"

#include <iostream>
#include <unistd.h>

using namespace std;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
OsInfo::OsInfo() {}

/**.......................................................................
 * Destructor.
 */
OsInfo::~OsInfo() {}

bool OsInfo::isBigEndian()
{
  return !isLittleEndian();
}

bool OsInfo::isLittleEndian()
{
  unsigned int intTest = 1;
  return (intTest & 0x1) == 1;
}

void OsInfo::printBits(unsigned int iVal)
{
  std::cout << iVal << " = " << std::hex << iVal << " = ";
  for(unsigned i=0; i < 32; i++) {
    if(i%4 == 0) {
      std::cout << " ";
    }
    std::cout << ((iVal >> i) & 0x1);
  }
  std::cout << std::endl;
}

int OsInfo::getNumberOfCpus(void)
{
  long nprocs       = -1;
  long nprocs_max   = -1;
  
  nprocs = sysconf(_SC_NPROCESSORS_ONLN);

  if(nprocs < 1) {
    ThrowSysError("Could not determine number of CPUs online");
    return 0;
  }

  nprocs_max = sysconf( _SC_NPROCESSORS_CONF );

  if(nprocs_max < 1) {
    ThrowSysError("Could not determine number of CPUs in host");
    return 0;
  }

  //  COUT("nprocs = " << nprocs << " of " << nprocs_max << " online");

  return nprocs; 
}
