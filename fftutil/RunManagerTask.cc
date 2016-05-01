#include "gcp/fftutil/RunManagerTask.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................                                                                                    
 * Initialize the static RunManagerTask pointer to NULL.  This will be                                                                                         
 * used inside static functions (ie, signal handlers)                                                                                                         
 */
RunManagerTask* RunManagerTask::master_ = 0;

/**.......................................................................
 * Constructor.
 */
RunManagerTask::RunManagerTask(std::string file) 
{
  //------------------------------------------------------------
  // Keep a static pointer to ourselves; we will need this for access
  // to the RunManagerTask resources in signal handlers
  //------------------------------------------------------------

  master_ = this;
  cancelled_ = false;

  //------------------------------------------------------------
  // Set the run file
  //------------------------------------------------------------

  rm_.setRunFile(file);

  //------------------------------------------------------------
  // Create the signal task thread
  //------------------------------------------------------------

  threads_.push_back(new Thread(&startSignal,     &cleanSignal,     NULL, "Signal",     0));
  threads_.push_back(new Thread(&startRunManager, &cleanRunManager, NULL, "RunManager", 1));

  //------------------------------------------------------------
  // Install signals
  //------------------------------------------------------------

  sendInstallSignalMsg(SIGINT, &sendSignalRcvdMsg);

  //------------------------------------------------------------
  // Start up threads
  //------------------------------------------------------------

  startThreads(this);
		     
  //------------------------------------------------------------
  // Finally, start our message queue
  //------------------------------------------------------------

  run();

  for(int ithread=0; ithread < threads_.size(); ithread++) 
    threads_[ithread]->setRunState(false);

  if(cancelled_) {
    ReportSimpleColorError(std::endl << "Program was cancelled", "red");
  }
}

/**.......................................................................
 * Destructor.
 */
RunManagerTask::~RunManagerTask() 
{
}

/**.......................................................................
 * Thread in which the RunManager will execute
 */
THREAD_START(RunManagerTask::startRunManager)
{
  RunManagerTask* parent = (RunManagerTask*) arg;
  Thread* thread = 0;

  try {

    thread = parent->getThread("RunManager");
    thread->broadcastReady();

    if(!parent->cancelled_) {
      parent->rm_.run();
    }

    if(!parent->cancelled_) {
      parent->rm_.listModel();
    }

  } catch(Exception& err) {

    if(!parent->cancelled_) {
      ReportError(err.what());
    }
  }

  parent->sendStopMsg();

  return 0;
}

/**.......................................................................
 * Thread in which the RunManager will execute
 */
THREAD_START(RunManagerTask::startSignal)
{
  RunManagerTask* parent = (RunManagerTask*) arg;
  Thread* thread = 0;

  try {

    parent->signal_ = new gcp::util::SignalTask();

    thread = parent->getThread("Signal");
    thread->broadcastReady();

    parent->signal_->run();

  } catch(Exception& err) {

    ReportError(err.what());

    if(parent->signal_) {
      delete parent->signal_;
      parent->signal_ = 0;
    }
  }

  return 0;
}

THREAD_CLEAN(RunManagerTask::cleanSignal)
{
  RunManagerTask* master = (RunManagerTask*) arg;
  Thread* thread = 0;
  
#if 1
  if(master->signal_ != 0) {
    delete master->signal_;
    master->signal_ = 0;
  }
#endif

  thread = master->getThread("Signal");
}

THREAD_CLEAN(RunManagerTask::cleanRunManager)
{
  RunManagerTask* master = (RunManagerTask*) arg;
  Thread* thread = 0;
  
  thread = master->getThread("RunManager");
}

/**
 * Send a message that a signal was received
 */
SIGNALTASK_HANDLER_FN(RunManagerTask::sendSignalRcvdMsg)
{
  master_->cancelled_ = true;
  master_->sendStopMsg();
}



