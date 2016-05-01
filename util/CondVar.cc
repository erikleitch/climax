#include "gcp/util/CondVar.h"
#include "gcp/util/Exception.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
CondVar::CondVar() 
{
  // Initialize the mutex and cond var

  pthread_mutex_init(&guard_, NULL);
  pthread_cond_init(&cond_, NULL);
}

/**.......................................................................
 * Destructor.
 */
CondVar::~CondVar() 
{
  pthread_mutex_destroy(&guard_);
  pthread_cond_destroy(&cond_);
}

void CondVar::lockMutex()
{
  if(pthread_mutex_lock(&guard_) != 0)
    ThrowSysError("pthread_mutex_lock");
}

/**.......................................................................
 * Version of waitUntilReady that assumes the mutex has already been
 * locked
 */
void CondVar::waitUntilReadyPrelock()
{
  // Now call cond_wait.  If we fail, unlock the mutex before exiting

  try {

    if(pthread_cond_wait(&cond_, &guard_) != 0)
      ThrowSysError("pthread_cond_wait");

  } catch(Exception& err) {
    pthread_mutex_unlock(&guard_);
    throw err;
  }

  if(pthread_mutex_unlock(&guard_) != 0) 
    ThrowSysError("pthread_mutex_unlock");
}

void CondVar::waitUntilReady()
{
  // Lock the mutex before calling cond_wait

  if(pthread_mutex_lock(&guard_) != 0)
    ThrowSysError("pthread_mutex_lock");

  // Now call cond_wait.  If we fail, unlock the mutex before exiting

  try {

    if(pthread_cond_wait(&cond_, &guard_) != 0)
      ThrowSysError("pthread_cond_wait");

  } catch(Exception& err) {
    pthread_mutex_unlock(&guard_);
    throw err;
  }

  if(pthread_mutex_unlock(&guard_) != 0) 
    ThrowSysError("pthread_mutex_unlock");
}

void CondVar::broadcastReady()
{
  if(pthread_mutex_lock(&guard_) != 0)
    ThrowSysError("pthread_mutex_lock");

  // The meat of the function -- just broadcast to threads waiting on
  // this condition variable.
   
  try {

    if(pthread_cond_broadcast(&cond_) != 0)
      ThrowSysError("pthread_cond_broadcast");   
   
  } catch(Exception& err) {
    pthread_mutex_unlock(&guard_);
    throw err;
  }

  // Unlock the mutex.

  if(pthread_mutex_unlock(&guard_) != 0)
    ThrowSysError("pthread_mutex_unlock");
}
