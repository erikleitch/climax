#include "gcp/util/ThreadSynchronizer.h"
#include "gcp/util/Exception.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ThreadSynchronizer::ThreadSynchronizer() {}

/**.......................................................................
 * Destructor.
 */
ThreadSynchronizer::~ThreadSynchronizer() {}

void ThreadSynchronizer::resize(unsigned nBit)
{
  doneMask_.resize(nBit);
}

void ThreadSynchronizer::reset(unsigned nBit)
{
  resize(nBit);
  reset();
  initWait();
}

void ThreadSynchronizer::reset()
{
  doneMask_.setAllBitsLow();
  //  COUT("Just reset: doneMask_= " << std::endl << doneMask_);
}

void ThreadSynchronizer::registerPending(unsigned iBit)
{
  doneMaskGuard_.lock();
  doneMask_.setBitLow(iBit);
  doneMaskGuard_.unlock();
}

void ThreadSynchronizer::registerDone(unsigned iBit, unsigned nBit)
{
  doneMaskGuard_.lock();

  doneMask_.setBitHigh(iBit);

  //  CTOUT(pthread_self() << " Set bit " << iBit << " out of " << nBit << " high");

  try {
    if(doneMask_.firstNBitsAreHigh(nBit)) {
      //      CTOUT("Broadcasting done (iBit = " << iBit << ") doneMask_ = " << std::endl << doneMask_);
      doneVar_.broadcastReady();
    }
  } catch(...) {
  }

  doneMaskGuard_.unlock();
}

void ThreadSynchronizer::initWait()
{ 
  doneVar_.lockMutex();
}

void ThreadSynchronizer::wait()
{ 
  doneVar_.waitUntilReadyPrelock();
}
