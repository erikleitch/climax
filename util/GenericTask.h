#ifndef GENERICTASK_H
#define GENERICTASK_H

/**
 * @file GenericTask.h
 * 
 * Tagged: Fri Nov 14 12:39:33 UTC 2003
 * 
 * @author Erik Leitch
 */
#include <vector>
#include <unistd.h>
#include <sys/socket.h> // shutdown()
#include <sys/time.h>
#include <sys/types.h>

#include "gcp/util/Command.h"
#include "gcp/util/Debug.h"     
#include "gcp/util/IoLock.h"
#include "gcp/util/Exception.h" // Definition of Error macro
#include "gcp/util/FdSet.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/PipeQ.h"
#include "gcp/util/Thread.h"

namespace gcp {
  namespace util {
    
    /**
     * All tasks will have the following functionality:
     */
    template<class Msg>
      class GenericTask {
      
      public:
      
      /**
       * Method to send a stop message to this task via its
       * message queue.
       */
      void sendRestartMsg();
      
      /**
       * Method to send a stop message to this task via its
       * message queue.
       */
      void sendStopMsg();
      
      /**
       * Method to send a heartbeat message to this task via its
       * message queue.
       */
      void sendHeartBeatMsg();
      
      /**
       * Forward a message to this task via its message queue.  
       *
       * A public method to forward messages to this task via its
       * message queue.  We make this virtual so that inheriting
       * tasks can control precisely how messages from the outside
       * world are routed.
       *
       * protected method sendTaskMsg() (below) is available to
       * inheritors to send messages to themselves.
       *
       * @param msg Pointer to the message to be sent to
       * this task.
       */
      virtual void fwdTaskMsg(Msg* msg);
      
      protected:
      
      /**
       * Protected constructor ensures that the base class cannot
       * be instantiated.
       */
      GenericTask();
      
      /**
       * Constructor which initializes thread_, below.
       */
      GenericTask(Thread* thread);
      
      /**
       * Making the destructor virtual ensures that the right
       * destructor will be called for classes which inherit from
       * GenericTask.
       */
      virtual ~GenericTask();
      
      /**
       * Send a message to this task via its message queue.
       *
       * @param msg Pointer to the message to be sent to this
       * task.
       */
      void sendTaskMsg(Msg* msg);
      
      /**
       * If this GenericTask object was instantiated by another
       * thread, keep a pointer to it here.  If non-NULL, this
       * will be used to broadcast to other threads that this task
       * is shutting down (see ~GenericTask(), below).
       */
      Thread* thread_;
      
      /**
       * A vector of Thread objects managed by this task. 
       */
      std::vector<Thread*> threads_;
      
      /**
       * A method to start all threads.  This calls the run()
       * method of all Threads managed by this task, which calls
       * pthread_create() for each thread.
       *
       * @param arg The argument to be passed to the thread
       * startup function specified in the Thread constructor.
       */
      void startThreads(void* arg);
      
      /**
       * Start the next thread with the specified order
       */
      void startThread(void* arg, unsigned order);
      
      /**
       * Method to return the minimum start order for threads which
       * are not yet running.
       */
      unsigned getMinStartOrder();
      
      /**
       * Return true if there are still unstarted threads
       */
      bool threadsNeedStarting();
      
      /**
       * A method to cancel threads managed by this task.
       */
      void cancelThreads();
      
      /**
       * Cancel the next thread with the specified order
       */
      void cancelThread(unsigned order);
      
      /**
       * Method to return the minimum cancel order for threads which
       * are not yet running.
       */
      unsigned getMinCancelOrder();
      
      /**
       * Return true if there are still uncancelled threads
       */
      bool threadsNeedCancelling();
      
      /**
       * A method to ping all pingable threads managed by this task.
       *
       * @param arg The argument to be passed to the ping function
       * specified in the Thread constructor.
       */
      void pingThreads(void* arg);
      
      /**
       * Raise a signal to a named thread.
       */
      void raise(std::string name, int sigNo);
      
      /**
       * A method to start all threads managed by this task
       * running.
       */
      Thread* getThread(std::string name);
      
      /**
       * A method to test if the threads managed by this task are
       * running.
       */
      bool threadsAreRunning();
      
      /**
       * A message queue, implemented as a pipe, by which we can
       * communicate with this task.
       */
      PipeQ<Msg> msgq_;
      
      /**
       * A set of file descriptors associated with this task.
       */
      gcp::util::FdSet fdSet_;
      
      // Fd utilities
      
      /**
       * Shutdown a connection.  This method calls shutdown() and
       * close() on the file descriptor, and clears it from the set of
       * desfcriptors to be watched by this task.
       */
      void shutdownConnection(int fd);
      
      /**
       * This routine will simply block, servicing messages on the
       * message queue.  This is declared virtual so that
       * inheritors can overload with their own methods.
       *
       * @throws Exception
       */
      virtual void serviceMsgQ(void);
      
      /**
       * Restart this thread.
       */
      virtual void restart(void);
      
      /**
       * Force inheritors to define a run method.
       */
      virtual void run(void);
      
      /**
       * Process a message received on our message queue
       *
       * @throws Exception
       */
      virtual void processTaskMsg(bool* stop);
      
      /**
       * This method should be defined by each inheriting task to
       * process its own task-specific messages.
       *
       * @throws Exception
       */
      virtual void processMsg(Msg* msg);
      
      /**
       * Respond to a heartbeat message.
       *
       * @throws Exception
       */
      virtual void respondToHeartBeat();
      
      //------------------------------------------------------------
      // Signals and timers
      //------------------------------------------------------------
      
      /**
       * Respond to a message to install a timer.
       */
      virtual void installTimer(Msg* msg);
      
      /**
       * Respond to a message to install a signal.
       */
      virtual void installSignal(Msg* msg);
      
      /**
       * Respond to a message to enable/disable a timer.
       */
      virtual void enableTimer(Msg* msg);
      
      /**
       * Respond to a message to add/remove a handler
       */
      virtual void addHandler(Msg* msg);
      
      //------------------------------------------------------------
      // A vector of commands used by this task
      //------------------------------------------------------------

      std::vector<Command*> commands_;

      struct timeval* timeOut_;

      TimeVal commandTimeOut_;

      void stepCommands();
      
    }; // End class GenericTask
    
    // Template functions must be defined in the header file.
    
    /**
     * Constructor
     */
    template<class Msg>
      GenericTask<Msg>::GenericTask() 
      {
	thread_  = 0;
	timeOut_ = 0;
	
	// Register the message queue fd to be watched for readability
	
	fdSet_.zeroReadFdSet();
	fdSet_.registerReadFd(msgq_.fd());
      }
    
    /**
     * Constructor
     */
    template<class Msg>
      GenericTask<Msg>::GenericTask(Thread* thread) 
      {
	thread_  = thread;
	timeOut_ = 0;
	
	// Register the message queue fd to be watched for readability
	
	fdSet_.zeroReadFdSet();
	fdSet_.registerReadFd(msgq_.fd());
      }
    
    /**
     * Destructor
     */
    template<class Msg>
      GenericTask<Msg>::~GenericTask()
      {
	// Broadcast to other threads that this task is shutting
	// down.
	
	if(thread_ != 0)
	  thread_->setRunState(false); 
	
	// Cancel all threads managed by this task.
	
	cancelThreads();
	
	// And delete any allocated memory.  NB: the default
	// destructor for vector<Thread*> will not delete memory
	// pointed to by its elements.
	
	for(unsigned ithread=0; ithread < threads_.size(); ithread++)
	  if(threads_[ithread] != 0) {
	    delete threads_[ithread];
	  }

	// Delete any commands too

	for(unsigned iCommand=0; iCommand < commands_.size(); iCommand++)
	  if(commands_[iCommand] != 0)
	    delete commands_[iCommand];
      }
    
    /**
     * Restart method stub
     */
    template<class Msg>
      void GenericTask<Msg>::restart(void) {}
    
    /**
     * Run method defaults to calling the serviceMsgQ() method
     * below.
     */
    template<class Msg>
      void GenericTask<Msg>::run(void)
      {
	serviceMsgQ();
      }
    
    /**
     * Main Task event loop: when this is called, the task blocks forever
     * in select(), or until a stop message is received.
     *yes
     * @throws Exception (also via msgq_.rfds())
     */
    template<class Msg>
      void GenericTask<Msg>::serviceMsgQ(void) 
      {
	bool stop=false;
	int nready; // number of file descriptors ready for reading
	
	if(msgq_.fd() < 0) {
	  ThrowError("Received NULL file descriptor");
	}

	// Loop, checking the message queue file descriptor for readability
	
	while(!stop && (nready=select(fdSet_.size(), fdSet_.readFdSet(), 
				      NULL, NULL, NULL)) > 0) {
	  
	  // If no file descriptors were ready, throw an exception
	  
	  if(nready != 1)
	    ThrowError("Error");
	  
	  DBPRINT(true, Debug::DEBUG2, "About to call processTaskMsg: "
		  << "nready = "  << nready);
	  
	  
	  processTaskMsg(&stop);

	};

      };
    
    /**
     *  Method to process a message received on the Task message * queue
     */
    template<class Msg>
      void GenericTask<Msg>::processTaskMsg(bool* stop)
      {
	Msg msg;
	
	msgq_.readMsg(&msg);
	
	DBPRINT(true, Debug::DEBUG2, "Inside processTaskMsg: " 
		<< msg.genericMsgType_);
	
	switch (msg.genericMsgType_) {
	case Msg::HEARTBEAT: // Is this a heartbeat request?
	  respondToHeartBeat();
	  break;
	case Msg::RESTART:  // Is this a request to restart?
	  restart();
	  break;
	case Msg::STOP:     // Did we receive a request to shut
	  // down?
	  *stop = true;
	  break;
	default: // Else forward this message to the task-specific
	  // process method
	  processMsg(&msg);
	  break;
	}

      };
    
    /**
     *  Method to process a message received on the Task message * queue
     */
    template<class Msg>
      void GenericTask<Msg>::processMsg(Msg* msg) {};
    
    /**
     * Each thread should respond to a heartbeat request by re-setting its
     * running state to true
     */
    template<class Msg>
      void GenericTask<Msg>::respondToHeartBeat()
      {
	if(thread_ != 0)
	  thread_->setRunState(true);
      }
    
    /**
     * Send a restart request to this thread
     */
    template<class Msg>
      void GenericTask<Msg>::sendRestartMsg() 
      {
	Msg msg;
	msg.genericMsgType_ = Msg::RESTART;
	
	// We use fwdTaskMsg() here instead of sendTaskMsg() so
	// that inheriting tasks can control their own messages
	// routing.
	
	fwdTaskMsg(&msg);
      }
    
    /**
     * Send a shutdown request to this thread
     */
    template<class Msg>
      void GenericTask<Msg>::sendStopMsg() 
      {
	Msg msg;
	msg.genericMsgType_ = Msg::STOP;
	
	// We use fwdTaskMsg() here instead of sendTaskMsg() so
	// that inheriting tasks can control their own messages
	// routing.
	
	fwdTaskMsg(&msg);
      };
    
    /**
     * Send a heartbeat request to this task.
     */
    template<class Msg>
      void GenericTask<Msg>::sendHeartBeatMsg() 
      {
	Msg msg;
	msg.genericMsgType_ = Msg::HEARTBEAT;
	
	// We use fwdTaskMsg() here instead of sendTaskMsg() so
	// that inheriting tasks can control their own message
	// routing.
	
	fwdTaskMsg(&msg);
      };
    
    /**
     * A task-independent method for sending a message to a task via its
     * message queue
     */
    template<class Msg>
      void GenericTask<Msg>::sendTaskMsg(Msg* msg)
      { 
	//	DBPRINT(true, Debug::DEBUG10, "Inside sendTaskMsg for thread: "
	//		<< thread_->strName());
	
	try {
	  DBPRINT(true, Debug::DEBUG7, "Sending a message " << pthread_self());
	  msgq_.sendMsg(msg); 
	  DBPRINT(true, Debug::DEBUG7, "Sending a message: done " << pthread_self());
	} catch(Exception& err) {
	  COUT("Caught an exception in task instantiated by thread: " << thread_->name());
	  
	  // Rethrow
	  
	  throw err;
	}
      }
    
    /**
     * A task-independent method for forwarding a message to a task
     * via its message queue.
     */
    template<class Msg>
      void GenericTask<Msg>::fwdTaskMsg(Msg* msg)
      { 
	DBPRINT(true, Debug::DEBUG2, "Inside fwdTaskMsg");

	try {
	  msgq_.sendMsg(msg); 
	} catch(Exception& err) {
	  COUT("Caught an exception in task instantiated by thread: " << thread_->name());
	  
	  // Rethrow
	  
	  throw err;
	}
      }
    
    /**
     * Method to start the next thread with a given start priority
     */
    template<class Msg>
      void GenericTask<Msg>::startThreads(void* arg)
      { 
	// Start threads in priority order
	
	while(threadsNeedStarting()) {
	  unsigned order = getMinStartOrder();
	  startThread(arg, order);
	}
      }
    
    /**
     * Method to start the next thread with a given start priority
     */
    template<class Msg>
      void GenericTask<Msg>::startThread(void* arg, unsigned order)
      { 
	// Search through the list in first-in, first-started order
	
	for(int ithread=0; ithread < threads_.size(); ithread++) 
	  if(!threads_[ithread]->isRunning() && threads_[ithread]->startOrder()==order) {
	    threads_[ithread]->start(arg);
	    
	    // Check the error code from the start up function
	    
	    if(threads_[ithread]->wasError_) {
	      gcp::util::LogStream errStr;
	      errStr.appendMessage(true, "Error in thread startup function\n");
	      throw Error(errStr);
	    }
	  }
      }
    
    /**
     * Return true if there are still unstarted threads
     */
    template<class Msg>
      bool GenericTask<Msg>::threadsNeedStarting()
      {
	for(int ithread=0; ithread < threads_.size(); ithread++) 
	  if(!threads_[ithread]->isRunning())
	    return true;
	return false;
      }
    
    /**
     * Return true if there are still running threads
     */
    template<class Msg>
      bool GenericTask<Msg>::threadsNeedCancelling()
      {
	for(int ithread=0; ithread < threads_.size(); ithread++) 
	  if(threads_[ithread]->isRunning())
	    return true;
	return false;
      }
    
    /**
     * Method to return the minimum start order for threads which are not yet running.
     */
    template<class Msg>
      unsigned GenericTask<Msg>::getMinStartOrder()
      {
	bool first=true;
	unsigned minOrder;
	
	for(int ithread=0; ithread < threads_.size(); ithread++) {
	  if(!threads_[ithread]->isRunning()) {
	    
	    if(first) {
	      minOrder = threads_[ithread]->startOrder();
	      first = false;
	    } else {
	      minOrder = threads_[ithread]->startOrder() < minOrder ? 
		threads_[ithread]->startOrder() : minOrder;
	    }

	  }
	}   
	
	return minOrder;
      }
    
    /**
     * Method to return the minimum stop order for threads which are still running.
     */
    template<class Msg>
      unsigned GenericTask<Msg>::getMinCancelOrder()
      {
	bool first=true;
	unsigned minOrder;
	
	for(int ithread=0; ithread < threads_.size(); ithread++) {
	  if(threads_[ithread]->isRunning()) {
	    
	    if(first) {
	      minOrder = threads_[ithread]->cancelOrder();
	      first = false;
	    } else {
	      minOrder = threads_[ithread]->cancelOrder() < minOrder ? 
		threads_[ithread]->cancelOrder() : minOrder;
	    }

          }
	}
	
	return minOrder;
      }
    
    /**
     * Method to cancel all threads managed by this task.
     *
     * This will wait until all spawned threads have shut down,
     * since each call to Thread::cancel() waits on the 'done'
     * condition variable, which is signalled by each thread in
     * its cancellation cleanup handler.
     */
    template<class Msg>
      void GenericTask<Msg>::cancelThreads()
      { 
	// Cancel threads in first-started, last-canceled order.
	
	DBPRINT(true, Debug::DEBUG7, "Managing " << threads_.size() 
		<< " threads");

	while(threadsNeedCancelling()) {
	  unsigned order = getMinCancelOrder();
	  cancelThread(order);
	}
	
	DBPRINT(true, Debug::DEBUG7, "Leaving cancelThreads");
      }
    
    /**
     * Method to cancel all threads managed by this task.
     *
     * This will wait until all spawned threads have shut down,
     * since each call to Thread::cancel() waits on the 'done'
     * condition variable, which is signalled by each thread in
     * its cancellation cleanup handler.
     */
    template<class Msg>
      void GenericTask<Msg>::cancelThread(unsigned order)
      { 
	// Cancel threads in first-started, last-canceled order.
	
	for(int ithread=threads_.size()-1; ithread >=0; ithread--) {
	  std::cout << " thread : " << ithread << std::endl;
	  if(threads_[ithread]->isRunning() && threads_[ithread]->cancelOrder() == order) {
	    std::cout << " canceling thread : " << ithread << std::endl;
	    threads_[ithread]->cancel();
	  }
	}
      }
    
    /**
     * Method to ping all threads managed by this task.
     */
    template<class Msg>
      void GenericTask<Msg>::pingThreads(void* arg)
      { 
	for(int ithread=0; ithread < threads_.size(); ithread++) 
	  if(threads_[ithread]->isPingable()) {
	    threads_[ithread]->setRunState(false);
	    threads_[ithread]->ping(arg);    
	  }
      }
    
    /**
     * Return a pointer to the thread of this name
     *
     * Returns:
     *  The requested pointer, or NULL on error
     */
    template<class Msg>
      Thread* GenericTask<Msg>::getThread(std::string name)
      {
	for(int ithread=0; ithread < threads_.size(); ithread++) {
	  if(threads_[ithread]->matchName(name))
	    return threads_[ithread];
	}
	
	ThrowError("No matching thread found");
        return threads_[0];
      };
    
    
    // Raise a signal to a thread.
    
    template<class Msg>
      void GenericTask<Msg>::raise(std::string name, int sigNo)
      {
	for(int ithread=0; ithread < threads_.size(); ithread++) {
	  if(threads_[ithread]->matchName(name)) {
	    threads_[ithread]->raise(sigNo);
	    return;
	  }
	}
	
	ReportError("No matching thread found");
      };
    
    
    // Check the running status of all threads
    
    template<class Msg>
      bool GenericTask<Msg>::threadsAreRunning()
      {
	for(int ithread=0; ithread < threads_.size(); ithread++) {
	  if(!threads_[ithread]->isRunning()) {
	    ReportError("Thread " 
			<< threads_[ithread]->strName() 
			<< " is not running!");
	    
	    return false;
	  } else {
	    //	    COUT("Thread " 
	    //		 << threads_[ithread]->strName() 
	    //		 << " is running.");
	  }
	}
	return true;
      }
    
    //------------------------------------------------------------
    // Signals and timers
    //------------------------------------------------------------
    
    /**
     * Respond to a message to install a timer.
     */
    template<class Msg>
      void GenericTask<Msg>::installTimer(Msg* msg) {};
    
    /**
     * Respond to a message to install a signal.
     */
    template<class Msg>
      void GenericTask<Msg>::installSignal(Msg* msg) {};
    
    /**
     * Send a message to the master thread to install a timer.
     */
    template<class Msg>
      void GenericTask<Msg>::enableTimer(Msg* msg) {};
    
    /**
     * Send a message to the signal thread to add a handler.
     */
    template<class Msg>
      void GenericTask<Msg>::addHandler(Msg* msg) {};
    
    // Utilities
    
    /**
     * Shutdown a connection. 
     */
    template<class Msg>
      void GenericTask<Msg>::shutdownConnection(int fd)
      {
	if(fd >= 0) {
	  fdSet_.clear(fd);
	  ::shutdown(fd, 2);
	  ::close(fd);
	}
	
      }
    
    /**.......................................................................
     * Step through any commands we are running
     */
    template<class Msg>
      void GenericTask<Msg>::stepCommands() {

      // On entry, set the timeout to NULL.  If any command is still
      // active, we will set the timeout pointing to the timeout
      // returned by the command.  This is so that inheritors can use
      // timeOut_ without worrying about setting it back to NULL when
      // all commands are finished

      timeOut_ = 0;
      bool first=true;
      for(unsigned iCommand=0; iCommand < commands_.size(); iCommand++) {

	// If any command is still active, execute its next
	// instruction, and record the timeout.

	if(commands_[iCommand]->active()) {
	  commands_[iCommand]->executeNextInstruction(commandTimeOut_, first);
	  timeOut_ = commandTimeOut_.timeVal();

	  if(first)
	    first = false;
	}
      }
    }
    
  } // End namespace util
} // End namespace gcp


#endif // End #ifndef 
