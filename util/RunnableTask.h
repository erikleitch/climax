// $Id: RunnableTask.h,v 1.1 2012/03/14 20:35:38 eml Exp $

#ifndef GCP_UTIL_RUNNABLETASK_H
#define GCP_UTIL_RUNNABLETASK_H

/**
 * @file RunnableTask.h
 * 
 * Tagged: Fri Jan 26 17:59:28 NZDT 2007
 * 
 * @version: $Revision: 1.1 $, $Date: 2012/03/14 20:35:38 $
 * 
 * @author Erik Leitch
 */

#include "gcp/util/Runnable.h"

#include <iostream>

namespace gcp {
  namespace util {

    class RunnableTask : public Runnable {
    public:

      /**
       * Constructor.
       */
      RunnableTask(bool spawnThread) :
	Runnable(spawnThread, runFn)
	{
	  // CAUTION: spawn() calls the startup function passed to
	  // pthread_create() when the thread corresponding to this
	  // object is created in Runnable().  This startup function
	  // calls our static runFn() method, which in turn calls the
	  // run() function of this object.  
	  //
	  // If we call spawn() in the constructor, there's no
	  // guarantee that the virtual table for this object has been
	  // constructed by the time run() is called, which means that
	  // sometimes the call to run() will correctly call the
	  // overloaded version in the inherited class, but sometimes
	  // it will call the virtual base-class method below.
	  //
	  // Therefore, this construct is not safe as-is, and spawn()
	  // should NOT be called here.  You should instead construct
	  // the object in external code, then call spawn() on it
	  // directly.

	  //	  spawn(this);

	}

      /**
       * Destructor.
       */
      virtual ~RunnableTask() {};

      virtual void run() {
      };

      /**
       * A run method to be called from pthread_start()
       */
      static RUN_FN(runFn) {
	RunnableTask* runnable = (RunnableTask*) arg;
	runnable->run();
      }

    }; // End class RunnableTask

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_RUNNABLETASK_H
