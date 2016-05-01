// $Id: ThreadSynchronizer.h,v 1.1 2012/05/16 21:26:44 eml Exp $

#ifndef GCP_UTIL_THREADSYNCHRONIZER_H
#define GCP_UTIL_THREADSYNCHRONIZER_H

/**
 * @file ThreadSynchronizer.h
 * 
 * Tagged: Wed May 16 13:50:51 PDT 2012
 * 
 * @version: $Revision: 1.1 $, $Date: 2012/05/16 21:26:44 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/BitMask.h"
#include "gcp/util/Mutex.h"
#include "gcp/util/CondVar.h"

namespace gcp {
  namespace util {

    class ThreadSynchronizer {
    public:

      /**
       * Constructor.
       */
      ThreadSynchronizer();

      /**
       * Destructor.
       */
      virtual ~ThreadSynchronizer();

      void initWait();
      void wait();
      void registerDone(unsigned iBit, unsigned nBit);
      void registerPending(unsigned iBit);
      void resize(unsigned nBit);
      void reset();
      void reset(unsigned nBit);

    private:

      BitMask doneMask_;
      Mutex doneMaskGuard_;
      CondVar doneVar_;

    }; // End class ThreadSynchronizer

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_THREADSYNCHRONIZER_H
