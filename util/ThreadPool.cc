#include "gcp/util/ThreadPool.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ThreadPool::ThreadPool(unsigned nThread) :
  SpawnableTask<ThreadPoolMsg>(true) 
{
  threads_.resize(nThread);

  for(unsigned iThread=0; iThread < nThread; iThread++) {

    ExecuteThread* thread = new ExecuteThread();

    thread->spawn();
    threads_[iThread] = thread;
    
    // And make this thread available for work

    availableThreads_.push_back(thread);
  }
}

/**.......................................................................
 * Constructor with vector of CPUs to bind to
 */
ThreadPool::ThreadPool(unsigned nThread, std::vector<unsigned>& cpus) :
  SpawnableTask<ThreadPoolMsg>(true) 
{
  unsigned nCpu = cpus.size();
  unsigned nThreadPerCpu = nCpu > nThread ? 1 : nThread/nCpu + 1;

  threads_.resize(nThread);
  bool first = true;

  for(unsigned iThread=0; iThread < nThread; iThread++) {

    unsigned iCpu = iThread / nThreadPerCpu;

    if(iCpu > nCpu)
      iCpu = nCpu-1;

    ExecuteThread* thread = new ExecuteThread(cpus[iCpu]);

    thread->spawn();
    threads_[iThread] = thread;
    
    // And make this thread available for work

    availableThreads_.push_back(thread);
  }
}

/**.......................................................................
 * Destructor.
 */
ThreadPool::~ThreadPool() 
{
  for(unsigned iThread=0; iThread < threads_.size(); iThread++) {
    if(threads_[iThread]) {
      delete threads_[iThread];
      threads_[iThread] = 0;
    }
  }
}

/**.......................................................................
 * Process a message for this object
 */
void ThreadPool::processMsg(ThreadPoolMsg* msg)
{
  switch(msg->type_) {
  case ThreadPoolMsg::EXECUTE:
    {
      // Allocate a new request at the point of request.  This will be
      // deleted by the worker thread when the execution has finished

      Request* req = new Request(msg->body_.execute_.fn_, msg->body_.execute_.args_, this);

      // And push it onto the queue of pending requests

      pendingRequests_.push_back(req);
    }
    break;
  case ThreadPoolMsg::THREAD_DONE:
    availableThreads_.push_back(msg->body_.threadDone_.thread_);
    break;
  default:
    ThrowError("Unrecognized message");
    break;
  }

  checkPendingRequests();
}

/**.......................................................................
 * Check for pending requests, servicing as threads become available
 */
void ThreadPool::checkPendingRequests()
{
  // Execute the minimum of the number of available threads and pending
  // requests

  unsigned nAvailThreads    = availableThreads_.size();
  unsigned nPendingRequests = pendingRequests_.size();

  unsigned nMin = nAvailThreads < nPendingRequests ? nAvailThreads : nPendingRequests;

  // Iterate over any pending requests that we can currently service

  for(unsigned iReq=0; iReq < nMin; iReq++) {

    Request*       req    = pendingRequests_.front();
    ExecuteThread* thread = availableThreads_.front();

    req->thread_ = thread;

    thread->execute(&executeCallback, (void*)req);

    pendingRequests_.pop_front();
    availableThreads_.pop_front();
  }
}

/**.......................................................................
 * A function to be called by a worker thread to execute a work request
 */
EXECUTE_FN(ThreadPool::executeCallback)
{
  Request* req = (Request*)args;

  // Execute the call

  req->fn_(req->args_);

  // Now notify the originator of this request that we are finished

  req->originator_->registerDone(req->thread_);

  // And delete the request (which was allocated by the originator)

  delete req;
}

void ThreadPool::registerDone(ExecuteThread* thread)
{
  ThreadPoolMsg msg;
  msg.genericMsgType_ = GenericTaskMsg::TASK_SPECIFIC;

  msg.type_ = ThreadPoolMsg::THREAD_DONE;

  msg.body_.threadDone_.thread_ = thread;

  sendTaskMsg(&msg);
}

void ThreadPool::execute(EXECUTE_FN(*fn), void* args)
{
  ThreadPoolMsg msg;
  msg.genericMsgType_ = GenericTaskMsg::TASK_SPECIFIC;

  msg.type_ = ThreadPoolMsg::EXECUTE;

  msg.body_.execute_.fn_   = fn;
  msg.body_.execute_.args_ = args;

  sendTaskMsg(&msg);
}

unsigned ThreadPool::nThread()
{
  return threads_.size();
}
