// $Id: ExecuteThread.h,v 1.2 2012/05/30 16:53:11 eml Exp $

#ifndef GCP_UTIL_EXECUTETHREAD_H
#define GCP_UTIL_EXECUTETHREAD_H

/**
 * @file ExecuteThread.h
 * 
 * Tagged: Wed Mar 14 13:39:47 PDT 2012
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/30 16:53:11 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/GenericTask.h"
#include "gcp/util/GenericTaskMsg.h"
#include "gcp/util/Runnable.h"
#include "gcp/util/SpawnableTask.h"

#define EXECUTE_FN(fn) void (fn)(void* args)

namespace gcp {
  namespace util {

    //------------------------------------------------------------
    // A utility class for sending messages to this object
    //------------------------------------------------------------

    class ExecuteThreadMsg : public GenericTaskMsg {
    public:

      enum MsgType {
	EXECUTE
      };
      
      union {

	struct {
	  EXECUTE_FN(*fn_);
	  void*args_;
	} execute_;

      } body_;

      // A type for this message

      MsgType type_;
    };

    //-----------------------------------------------------------------------
    // Gegneric class for executing in a separate thread
    //-----------------------------------------------------------------------

    class ExecuteThread : public SpawnableTask<ExecuteThreadMsg> {
    public:

      /**
       * Constructor.
       */
      ExecuteThread();
      ExecuteThread(unsigned cpu);

      /**
       * Destructor.
       */
      virtual ~ExecuteThread();

      void execute(EXECUTE_FN(*fn), void* args=0);

      void processMsg(ExecuteThreadMsg* msg);

    }; // End class ExecuteThread

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_EXECUTETHREAD_H
