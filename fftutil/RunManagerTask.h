// $Id: $

#ifndef GCP_UTIL_RUNMANAGERTASK_H
#define GCP_UTIL_RUNMANAGERTASK_H

/**
 * @file RunManagerTask.h
 * 
 * Tagged: Tue Dec 23 10:00:20 PST 2014
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include "gcp/util/GenericMasterTask.h"

#include "gcp/fftutil/RunManager.h"

namespace gcp {
  namespace util {

    class RunManagerTaskMsg : 
    public gcp::util::GenericMasterTaskMsg {

    public:

      enum MsgType {
	RUNMANAGER_DONE_MSG,
	SIGNAL_RCVD_MSG,
      };

      MsgType type_;

      inline void packRunManagerDoneMsg()
      {
	genericMsgType_ =
	  gcp::util::GenericTaskMsg::TASK_SPECIFIC;

	type_ = RUNMANAGER_DONE_MSG;
      }

      inline void packSignalRcvdMsg()
      {
	genericMsgType_ =
	  gcp::util::GenericTaskMsg::TASK_SPECIFIC;

	type_ = SIGNAL_RCVD_MSG;
      }

    };

    class RunManagerTask : 
    public gcp::util::GenericMasterTask<RunManagerTaskMsg> {

    public:

      /**
       * Constructor.
       */
      RunManagerTask(std::string file);

      /**
       * Destructor.
       */
      virtual ~RunManagerTask();

      RunManager rm_;

      // Send a message to the RunManagerTask that a signal has been received

      static SIGNALTASK_HANDLER_FN(sendSignalRcvdMsg);

      static THREAD_START(startSignal);
      static THREAD_CLEAN(cleanSignal);

      static THREAD_START(startRunManager);
      static THREAD_CLEAN(cleanRunManager);

      // A static pointer to ourselves for use in static functions

      static RunManagerTask* master_;

      bool cancelled_;

    }; // End class RunManagerTask

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_RUNMANAGERTASK_H
