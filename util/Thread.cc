#include "gcp/util/Exception.h"
#include "gcp/util/Debug.h"
#include "gcp/util/Directives.h"
#include "gcp/util/IoLock.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/Thread.h"

#include <signal.h>
#include <errno.h>

using namespace std;
using namespace gcp::util;

/**.......................................................................
 * Thread constructor function
 */
Thread::Thread(THREAD_START(*startFn), THREAD_CLEAN(*cleanFn),
	       THREAD_PING(*pingFn), string name, 
	       unsigned startOrder, unsigned cancelOrder,
	       int priority, int schedPolicy)
{
  privateConstructor(startFn, cleanFn, pingFn, name, startOrder, cancelOrder, 
		     priority, schedPolicy);
}

/**.......................................................................
 * Thread constructor function
 */
Thread::Thread(THREAD_START(*startFn), THREAD_CLEAN(*cleanFn),
	       THREAD_PING(*pingFn), string name, 
	       unsigned startOrder, unsigned cancelOrder)
{
  privateConstructor(startFn, cleanFn, pingFn, name, startOrder, cancelOrder, 
		     0, SCHED_OTHER);
}

/**.......................................................................
 * Thread constructor function
 */
void Thread::privateConstructor(THREAD_START(*startFn), THREAD_CLEAN(*cleanFn),
				THREAD_PING(*pingFn), string name, 
				unsigned startOrder, unsigned cancelOrder,
				int priority, int schedPolicy)
{
  wasError_    = false;

  startOrder_  = startOrder;
  cancelOrder_ = cancelOrder;

  //------------------------------------------------------------
  // Sanity check arguments.
  //------------------------------------------------------------

  if(startFn == 0)
    throw Error("Thread::Thread: No startup function supplied.\n");

  //------------------------------------------------------------
  // Initialize the startup and shutdown functions of this threa
  //------------------------------------------------------------d

  userStartFn_  = startFn;
  userCleanFn_  = cleanFn;
  userPingFn_   = pingFn;

  startupArg_ = 0;

  //------------------------------------------------------------
  // Initialize the name of this thread
  //------------------------------------------------------------

  name_ = name;

  //------------------------------------------------------------
  // Initialize the guard mutexes of this thread
  //------------------------------------------------------------

  pthread_mutex_init(&ready_.guard,  NULL);
  pthread_mutex_init(&done_.guard,   NULL);
  pthread_mutex_init(&runningGuard_, NULL);

  //------------------------------------------------------------
  // Initialize the ready_ condition variables for this thread
  //------------------------------------------------------------

  pthread_cond_init(&ready_.cond, NULL);
  pthread_cond_init(&done_.cond,  NULL);

  //------------------------------------------------------------
  // We initialize this thread to block all (blockable) signals.  I
  // have confirmed that sigfillset() adds all signals up to SIGRTMAX
  // into the signal mask.
  //------------------------------------------------------------

  //  sigset_t allSignals;
  //  sigfillset(&allSignals);
  //  pthread_sigmask(SIG_BLOCK, &allSignals, NULL);
    
  //------------------------------------------------------------
  // Initialize the attributes
  //------------------------------------------------------------

  pthread_attr_init(&attr_);

  //------------------------------------------------------------
  // Set the inheritance policy
  //------------------------------------------------------------

  if(pthread_attr_setinheritsched(&attr_, PTHREAD_EXPLICIT_SCHED) != 0)
    ThrowSysError("pthread_attr_setinheritsched()");

  //------------------------------------------------------------
  // Set scheduling policy
  //------------------------------------------------------------

  if(pthread_attr_setschedpolicy(&attr_, schedPolicy) != 0)
    ThrowSysError("pthread_attr_setschedpolicy()");

  //------------------------------------------------------------
  // Set the priority
  //------------------------------------------------------------

  struct sched_param params;
  params.sched_priority = priority;
  uid_t euid = getuid();

  if(priority != 0 && euid != 0)
    ThrowError("You must be setuid root to execute a thread ("
	       << name_ << ") with priority: " << priority);

  if(pthread_attr_setschedparam(&attr_, &params) != 0)
    ThrowSysError("pthread_attr_setschedparam()");

  //------------------------------------------------------------
  // Initialize the state of this thread to not running: this will be
  // set to true when the resources for this thread have been
  // allocated and the thread's message service is running
  //------------------------------------------------------------

  setRunStatePrivate(false);
}

/**.......................................................................
 * Thread destructor
 */
Thread::~Thread()
{
  //------------------------------------------------------------
  // Destroy the guard mutex and condition variables
  //------------------------------------------------------------

  pthread_mutex_destroy(&ready_.guard);
  pthread_mutex_destroy(&runningGuard_);

  pthread_cond_destroy(&ready_.cond);

  //------------------------------------------------------------
  // Don't destroy the done condition variable here.  The cleanup
  //  handler uses it to signal to any waiting threads that this
  //  thread has now gone away.  If we destroy it here, then the
  //  cleanup handler will throw an exception when it tries to signal
  //  this variable.  Instead we will destroy it in the cleanup
  //  handler.
  //------------------------------------------------------------

  // pthread_mutex_destroy(&done_.guard);
  // pthread_cond_destroy(&done_.cond);
}

/**.......................................................................
 * A public interface to setRunStatePrivate()
 */
void Thread::setRunState(bool state)
{
  if(setRunStatePrivate(state))
    throw Error("Thread::setRunState: Error in setRunStatePrivate()");
}

/**.......................................................................
 * Change the boolean variable that reflects the running state of this
 * thread.
 */
bool Thread::setRunStatePrivate(bool state)
{
  int oldtype;
  LogStream ls;

  // Before locking the mutex, we want to set pthread_mutex_unlock()
  // as a cleanup handler.  This guarantees that the mutex will be
  // unlocked even if this thread is cancelled before we reach the
  // next pthread_mutex_unlock() statement below.
  //
  // Note that if this thread was created with canceltype
  // PTHREAD_CANCEL_ASYNCHRONOUS, we need to switch to
  // PTHREAD_CANCEL_DEFERRED so that a cancellation cannot occur
  // between cleanup_push() and mutex_lock() below, resulting in the
  // cleanup handler trying to unlock a mutex that isn't locked by the
  // current thread

  INSTALL_MUTEX_CLEANUP(runningGuard_, ls);

  // The calling thread must lock the mutex before altering the
  // thread's running state; we don't want multiple threads to try and
  // change the state at the same time.

  if(pthread_mutex_lock(&runningGuard_) != 0)
    ls.appendSysError(true, "pthread_mutex_lock");

  // Now change the state of this Thread object

  running_ = state;

  // Release the guard mutex

  if(pthread_mutex_unlock(&runningGuard_) != 0)
    ls.appendSysError(true, "pthread_mutex_unlock");

  // Now remove the handler, but don't execute it, so we can test the
  // return value, below

  UNINSTALL_MUTEX_CLEANUP(ls);

  return ls.isError();
}

/**.......................................................................
 * Thread run function
 */
void Thread::start(void* arg) 
{
  bool waserr=false;
  int oldtype;
  LogStream ls;

  startupArg_ = arg;

  // We call pthread_create with a pointer to each thread's startup
  // function.  Each startup function initializes the subsystem and
  // sets the appropriate pointer of ant pointing to its resources.

  // The calling thread must lock the mutex before pthread_cond_wait()
  // can be called

  INSTALL_MUTEX_CLEANUP(ready_.guard, ls);

  // The calling thread must lock the mutex before altering the
  // thread's running state; we don't want multiple threads to try and
  // change the state at the same time.

  if(pthread_mutex_lock(&ready_.guard) != 0)
    ls.appendSysError(true, "pthread_mutex_lock");
  
  // Create the thread, with startup function pointed to by this
  // thread's startfn

  if(pthread_create(&id_, &attr_, &startThread, this) != 0)
    ls.appendSysError(true, "pthread_create");

  DBPRINT(false, Debug::DEBUG7, "Created " << name_ << " thread ("
	  << id_ << ")");

  // Suspend execution of the calling thread until the ready variable
  // is signalled.  This is kind of dangerous, in that it requires the
  // user startup function explicitly to call broadcastReady(), but it
  // allows the user to synchronize startup precisely.

  if(pthread_cond_wait(&ready_.cond, &ready_.guard) != 0)
    ls.appendSysError(true, "pthread_cond_wait");

  // Set the state of this thread to running
  
   if(setRunStatePrivate(true))
     ls.appendMessage(true, "Error in setRunStatePrivate");

  // And unlock the guard mutex

   if(pthread_mutex_unlock(&ready_.guard) != 0)
     ls.appendSysError(true, "pthread_mutex_unlock");

  // Now remove the handler, but don't execute it, so we can test the
  // return value of pthread_mutex_unlock, below

  UNINSTALL_MUTEX_CLEANUP(ls)

  if(ls.isError()) 
    throw Error(ls);
};

/**.......................................................................
 * Thread cancellation function.  It is expected that the exiting
 * thread will call a function established by pthread_cleanup_push()
 * which will signal ready when it is done.
 */
void Thread::cancel()
{
  bool waserr=false;
  int oldtype;
  LogStream ls;

  // Only proceed if the thread is currently running.  isRunning()
  // will return false if cancel has already been called on this
  // thread

  if(!isRunning())
    return;

  // Before locking the mutex, we want to set pthread_mutex_unlock()
  // as a cleanup handler.  This guarantees that the mutex will be
  // unlocked even if this thread is cancelled before we reach the
  // next pthread_mutex_unlock() statement below.
  //
  // Note that if this thread was created with canceltype
  // PTHREAD_CANCEL_ASYNCHRONOUS, we need to switch to
  // PTHREAD_CANCEL_DEFERRED so that a cancellation cannot occur
  // between cleanup_push() and mutex_lock() below, resulting in the
  // cleanup handler trying to unlock a mutex that isn't locked by the
  // current thread

  // Push cleanup handler onto the stack to unlock the mutex in case
  // we are cancelled while the mutex is locked.

  INSTALL_MUTEX_CLEANUP(done_.guard, ls);

  if(pthread_mutex_lock(&done_.guard) != 0)
    ls.appendSysError(true, "pthread_mutex_lock");

  DBPRINT(true, Debug::DEBUG7, "Cancelling thread id: " << id_ 
	  << " (" << name_ << ")");

  if(pthread_cancel(id_) != 0)
    ls.appendSysError(true, "pthread_cancel");

  DBPRINT(false, Debug::DEBUG7, "About to wait on done cond var");

  // Suspend execution of the calling thread until it is signalled
  // that the cancellation function is done.

  if(pthread_cond_wait(&done_.cond, &done_.guard) != 0)
    ls.appendSysError(true, "pthread_cond_wait");

  DBPRINT(true, Debug::DEBUG7, "Signalled done");

  // Set our running state to false

  if(setRunStatePrivate(false))
     ls.appendMessage(true, "Error in setRunStatePrivate");

  // Unlock the guard mutex

  if(pthread_mutex_unlock(&done_.guard) != 0)
    ls.appendSysError(true, "pthread_mutex_unlock");

  // And remove the handler, but don't execute it!

  UNINSTALL_MUTEX_CLEANUP(ls);

  if(ls.isError())
    throw Error(ls);
}

/**.......................................................................
 * Thread ping function
 */
void Thread::ping(void* arg) 
{
  if(isPingable())  // Don't ping if no ping method was passed
    (*userPingFn_)(arg);
}

/**.......................................................................
 * Return true if this thread is running
 */
bool Thread::isRunning()
{
  return running_;
}

/**.......................................................................
 * Broadcast to any waiting threads that we are ready
 */
void Thread::broadcastReady()
{
  bool waserr=false;
  int oldtype;
  LogStream ls;

  // Before locking the mutex, we want to set pthread_mutex_unlock()
  // as a cleanup handler.  This guarantees that the mutex will be
  // unlocked even if this thread is cancelled before we reach the
  // next pthread_mutex_unlock() statement below.
  //
  // Note that if this thread was created with canceltype
  // PTHREAD_CANCEL_ASYNCHRONOUS, we need to switch to
  // PTHREAD_CANCEL_DEFERRED so that a cancellation cannot occur
  // between cleanup_push() and mutex_lock() below, resulting in the
  // cleanup handler trying to unlock a mutex that isn't locked by the
  // current thread

   INSTALL_MUTEX_CLEANUP(ready_.guard, ls);

   // Now lock the mutex.  We have to explicitly lock here, to avoid a
   // race condition in which we broadcast ready before another thread
   // begins to wait, which will cause the other thread to wait
   // forever.  The model here is that a thread waiting on this
   // condition variable locks the mutex first, and doesn't release it
   // until pthread_cond_wait() gets called.  This means that the
   // following line will cause us to block until the other thread is
   // actually waiting.

   if(pthread_mutex_lock(&ready_.guard) != 0)
     ls.appendSysError(true, "pthread_mutex_lock");

   // The meat of the function -- just broadcast to threads waiting on
   // this condition variable.
   
   if(pthread_cond_broadcast(&ready_.cond) != 0)
     ls.appendSysError(true, "pthread_cond_broadcast");   

   // Unlock the mutex.
   
   if(pthread_mutex_unlock(&ready_.guard) != 0)
      ls.appendSysError(true, "pthread_mutex_unlock");

   // And remove the handler, but don't execute it!
   
   UNINSTALL_MUTEX_CLEANUP(ls);
   
   if(ls.isError())
     throw Error(ls);
}

/**.......................................................................
 * Broadcast to any waiting threads that we are done
 */
void Thread::broadcastDone()
{
  errno = pthread_cond_broadcast(&done_.cond);
  if(errno != 0) {
    ThrowSysError(pthread_self() << " Error in pthread_cond_broadcast()");
  }
}

/**.......................................................................
 * Return true if the passed name matches
 */
bool Thread::matchName(string compname)
{
  return (name_ == compname);
}

/**.......................................................................
 * Return the name string of this thread
 */
string Thread::strName()
{
  return name_;
}

//-----------------------------------------------------------------------
// Static definitions

/**.......................................................................
 * Define a wrapper around pthread_mutex_unlock() suitable for passing to
 * pthread_cleanup_push(), which expects a void (*fn)(void* arg)
 */
THREAD_CLEAN(Thread::unlockMutex)
{
  pthread_mutex_t* mut = (pthread_mutex_t*) arg;

  // Ignore the return value for now -- Anyway, can we throw an
  // exception in a cleanup handler??

  pthread_mutex_unlock(mut);
}

/**.......................................................................
 * Return true if a ping function has been installed for this thread.
 */
bool Thread::isPingable()
{
  return (userPingFn_ != 0);
}

/**.......................................................................
 * Raise a signal to this thread.
 */
void Thread::raise(int sigNo)
{
  pthread_kill(id_, sigNo);
}

/**.......................................................................
 * Wrapper around the user-supplied startup function.  This is the
 * method which is actually passed to pthread_create().  
*/
THREAD_START(Thread::startThread)
{
  Thread* thread = (Thread*) arg;
  LogStream errStr;

  // Set up this thread to block all signals

  sigset_t allSignals;
  sigfillset(&allSignals);
  pthread_sigmask(SIG_BLOCK, &allSignals, NULL);

  // Set the cancel type to asynchronous, which means that the thread
  // will immediately call its cleanup handlers and exit when
  // cancelled.

  pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);

  // Install a cleanup handler to be called on thread cancellation.
  //
  // We have to enclose pthread_cleanup_push() and
  // pthread_cleanup_pop() in the same block, at the same level of
  // nesting, since they are macros which introduce and close,
  // respectively, a pair of {} braces.

  pthread_cleanup_push(&cleanThread, arg);

  // Now call the user-supplied startup function, catching any
  // exceptions that may be thrown
  
  try {
    thread->userStartFn_(thread->startupArg_);
  } catch(Exception& err) {
    ReportError("Caught an exception in thread: "
		<< thread->name_ << " (" 
		<< thread->id_ << ")): "
		<< err.what());
    thread->wasError_ = true;
  }
  
  thread->running_ = false;

  // Remove and execute our cleanup handler.
    
  pthread_cleanup_pop(1);

  return 0;
}

/**.......................................................................
 * A wrapper around the user-supplied cleanup function.  This is the
 * method which is actually passed to pthread_cleanup_push().  
*/
THREAD_CLEAN(Thread::cleanThread)
{

  Thread* thread = (Thread*) arg;

  try {

    errno = 0;
    if(thread->userCleanFn_ != 0)
      thread->userCleanFn_(thread->startupArg_);

    DBPRINT(false, Debug::DEBUG7, "Done calling clean fn for: " << thread->name_);

    // Broadcast to other threads that we are done

    errno = 0;
    thread->broadcastDone();

  } catch(Exception& err) {

    CERR("Caught an exception in Thread::cleanThread "
	 << "(thread: " << thread->name_ << ", (" 
	 << thread->id_ << ") " << err.what());

  } catch(...) {

    CERR("Caught an exception in Thread::cleanThread "
	 << "(thread: " << thread->name_ << ", (" 
	 << thread->id_ << ").");
  }

  pthread_cond_destroy(&thread->done_.cond);
}

void Thread::waitUntilReady()
{
  // Suspend execution of the calling thread until the ready variable
  // is signalled.  
  
  if(pthread_cond_wait(&ready_.cond, &ready_.guard) != 0)
    ThrowError("A test");
}

void Thread::setAffinity(unsigned cpu)
{

  // Affinity specification not supported on MAC OS X

  #if MAC_OSX
  ThrowSimpleColorError("Mac OS X does not permit thread affinities to be specified from user space -- ignoring setAffinity() request", "red");
  #else

  // Initialize a CPU set to zero

  cpu_set_t cpuSet;
  CPU_ZERO(&cpuSet);

  // Set the requested CPU

  CPU_SET(cpu, &cpuSet);

  // Set the affinity in the attribute mask

  if(pthread_attr_setaffinity_np(&attr_, sizeof(cpu_set_t), &cpuSet) != 0)
    ThrowSysError("pthread_attr_setaffinity_np()");

  #endif
}
