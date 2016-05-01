#ifndef GCP_UTIL_IOLOCK_H
#define GCP_UTIL_IOLOCK_H

/**
 * @file IoLock.h
 * 
 * Tagged: Sat May  8 08:22:36 PDT 2004
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Mutex.h"

#define IOCOUT(statement) \
{\
    gcp::util::IoLock::lockCout(); \
    std::cout << statement << std::endl; \
    gcp::util::IoLock::unlockCout(); \
}

#define IOCTOUT(statement) \
{\
    gcp::util::TimeVal timeVal;\
    timeVal.setToCurrentTime();\
    gcp::util::IoLock::lockCout(); \
    std::cout << timeVal << ": " << statement << std::endl; \
    gcp::util::IoLock::unlockCout(); \
}

#define IOCERR(statement) \
{\
    gcp::util::IoLock::lockCerr(); \
    std::cerr << statement << std::endl; \
    gcp::util::IoLock::unlockCerr(); \
}

namespace gcp {
  namespace util {
    
    class IoLock {
    public:
      
      /**
       * Destructor.
       */
      virtual ~IoLock();
      
      static void lockCout();
      static void unlockCout();
      static void lockCerr();
      static void unlockCerr();

    private:

      static Mutex coutMutex_;
      static Mutex cerrMutex_;

      /**
       * Private constructor prevents instantiation
       */
      IoLock();

    }; // End class IoLock
    
  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_IOLOCK_H
