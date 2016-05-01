// $Id: Timer.h,v 1.1 2012/03/14 20:35:38 eml Exp $

#ifndef GCP_UTIL_TIMER_H
#define GCP_UTIL_TIMER_H

/**
 * @file Timer.h
 * 
 * Tagged: Fri Aug 14 10:54:50 PDT 2009
 * 
 * @version: $Revision: 1.1 $, $Date: 2012/03/14 20:35:38 $
 * 
 * @author tcsh: Erik Leitch
 */
#include "gcp/util/TimeVal.h"

namespace gcp {
  namespace util {

    class Timer {
    public:

      /**
       * Constructor.
       */
      Timer();

      /**
       * Destructor.
       */
      virtual ~Timer();

      void start();
      void stop();
      double deltaInSeconds();
      double integratedElapsedSeconds();

      friend std::ostream& operator<<(std::ostream& os, const Timer& timer);
      friend std::ostream& operator<<(std::ostream& os, Timer& timer);

    private:

      TimeVal start_;
      TimeVal stop_;
      TimeVal diff_;
      double integratedSec_;

    }; // End class Timer

    std::ostream& operator<<(std::ostream& os, const Timer& timer);
    std::ostream& operator<<(std::ostream& os, Timer& timer);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_TIMER_H
