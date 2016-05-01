#include <iostream>

#include "gcp/util/Debug.h"
#include "gcp/util/IoLock.h"

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
IoLock::IoLock() {}

// Create the static instance

//IoLock::IoLock();

// Initialize static variables

Mutex IoLock::coutMutex_;
Mutex IoLock::cerrMutex_;

/**.......................................................................
 * Destructor.
 */
IoLock::~IoLock() {}

void IoLock::lockCout()
{
  coutMutex_.lock();
}

void IoLock::unlockCout()
{
  coutMutex_.unlock();
}

void IoLock::lockCerr()
{
  cerrMutex_.lock();
}

void IoLock::unlockCerr()
{
  cerrMutex_.unlock();
}
