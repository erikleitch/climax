// $Id: ThreadPool.h,v 1.3 2012/05/30 16:53:11 eml Exp $

#ifndef GCP_UTIL_THREADPOOL_H
#define GCP_UTIL_THREADPOOL_H

/**
 * @file ThreadPool.h
 * 
 * Tagged: Wed Mar 14 14:30:10 PDT 2012
 * 
 * @version: $Revision: 1.3 $, $Date: 2012/05/30 16:53:11 $
 * 
 * @author tcsh: Erik Leitch
 */
#include "gcp/util/ExecuteThread.h"
#include "gcp/util/SpawnableTask.h"

#include <vector>
#include <deque>

namespace gcp {
  namespace util {

    //------------------------------------------------------------
    // A utility class for sending messages to this object
    //------------------------------------------------------------

    class ThreadPoolMsg : public GenericTaskMsg {
    public:

      enum MsgType {
	EXECUTE,
	THREAD_DONE,
      };
      
      union {

	// A request to execute some work

	struct {
	  EXECUTE_FN(*fn_);
	  void*args_;
	} execute_;

	// A message from a worker thread that it has completed some work

	struct {
	  ExecuteThread* thread_;
	} threadDone_;

      } body_;

      // A type for this message

      MsgType type_;
    };

    //-----------------------------------------------------------------------
    // Gegneric class for executing in a separate thread
    //-----------------------------------------------------------------------

    class ThreadPool  : public SpawnableTask<ThreadPoolMsg> {
    public:

      // A request to execute some work

      struct Request {

	EXECUTE_FN(*fn_);
	void* args_;
	ThreadPool* originator_;
	ExecuteThread* thread_;

	Request(EXECUTE_FN(*fn), void* args, ThreadPool* originator) {
	  fn_         = fn;
	  args_       = args;
	  originator_ = originator;
	  thread_     = 0;
	};

	~Request() {};
      };

      // Constructor.

      ThreadPool(unsigned nThread);

      // Constructor with vector of allowable CPUs

      ThreadPool(unsigned nThread, std::vector<unsigned>& cpus);

      /**
       * Destructor.
       */
      virtual ~ThreadPool();

      // Public methods of this class

      void execute(EXECUTE_FN(*fn), void* args=0);
      void registerDone(ExecuteThread* thread);

      void checkPendingRequests();
      void processMsg(ThreadPoolMsg* msg);

      static EXECUTE_FN(executeCallback);

      unsigned nThread();

    private:

      std::vector<ExecuteThread*> threads_;

      std::deque<Request*>       pendingRequests_;
      std::deque<ExecuteThread*> availableThreads_;

    }; // End class ThreadPool

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_THREADPOOL_H
