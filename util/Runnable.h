// $Id: Runnable.h,v 1.1 2012/03/14 20:35:38 eml Exp $

#ifndef GCP_UTIL_RUNNABLE_H
#define GCP_UTIL_RUNNABLE_H

/**
 * @file Runnable.h
 * 
 * Tagged: Tue Dec 21 19:19:59 CST 2004
 * 
 * @version: $Revision: 1.1 $, $Date: 2012/03/14 20:35:38 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Thread.h"
#include "gcp/util/Exception.h"

#define RUN_FN(fn) void (fn)(void* arg)

namespace gcp {
  namespace util {
    
    class Runnable {
    public:
      
      /**
       * Constructor.
       */
      Runnable(bool spawnThread, RUN_FN(*runFn));
      
      Runnable(bool spawnThread, RUN_FN(*runFn), int priority);
      
      /**
       * Destructor.
       */
      virtual ~Runnable();
      
      /**
       * A startup function for the spawned thread.
       */
      static THREAD_START(startUp);
      
      void blockForever();
      
    protected:
      
      virtual void broadcastReady() 
	{
	  spawnedThread_->broadcastReady();
	}

      /**
       * If this object is spawned in its own thread, we will use a
       * Thread container to manage it.
       */
      Thread* spawnedThread_;
      
      /**
       * True if this object is spawned in a separate thread.
       */
      bool spawned_; 
      
      // A pointer to a function which will be called on startup
      
      RUN_FN(*runFn_);
      
      void spawn(void* arg);
      
    public:
      
      void spawn();
      
    }; // End class Runnable
    
  } // End namespace util
} // End namespace gcp




#endif // End #ifndef GCP_UTIL_RUNNABLE_H
