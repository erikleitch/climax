// $Id: CondVar.h,v 1.2 2012/05/11 20:21:40 eml Exp $

#ifndef GCP_UTIL_CONDVAR_H
#define GCP_UTIL_CONDVAR_H

/**
 * @file CondVar.h
 * 
 * Tagged: Wed Mar 14 17:46:09 PDT 2012
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/11 20:21:40 $
 * 
 * @author Erik Leitch
 */
#include <pthread.h>

namespace gcp {
  namespace util {

    class CondVar {
    public:

      /**
       * Constructor.
       */
      CondVar();

      /**
       * Destructor.
       */
      virtual ~CondVar();

      void lockMutex();
      void waitUntilReadyPrelock();

      void waitUntilReady();
      void broadcastReady();

    private:

      // A guard mutex for the variable

      pthread_mutex_t guard_;
	
      // A variable which we will use to signal other threads

      pthread_cond_t cond_;

    }; // End class CondVar

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_CONDVAR_H
