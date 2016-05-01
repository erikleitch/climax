#ifndef GCP_UTIL_GENERICMASTERTASK_H
#define GCP_UTIL_GENERICMASTERTASK_H

/**
 * @file GenericMasterTask.h
 * 
 * Tagged: Tue Jul  6 17:41:36 PDT 2004
 * 
 * @author Erik Leitch
 */
#include <string>

#include "gcp/util/GenericTask.h"
#include "gcp/util/GenericMasterTaskMsg.h"
#include "gcp/util/SignalTask.h"
#include "gcp/util/TimeVal.h"

namespace gcp {
  namespace util {
    
    template<class Msg>
      class GenericMasterTask : 
      public gcp::util::GenericTask<Msg> {
      
      public:
      
      /**
       * Constructor.
       */
      GenericMasterTask();
      
      /**
       * Destructor.
       */
      virtual ~GenericMasterTask();
      
      protected:
      
      /**
       * Handles all signals for this process
       */
      gcp::util::SignalTask* signal_; 
      
      /**
       *  Method to process a message received on the Task message * queue
       */
      void processTaskMsg(bool* stop);
      
      /**
       * Send a message to the master thread to install a timer.  The
       * initial delay is set to the interval in this version.
       */
      void sendInstallTimerMsg(std::string name, int sigNo, 
			       unsigned long intervalSec, 
			       unsigned long intervalNsec,
			       SIGNALTASK_HANDLER_FN(*handler));
      
      /**
       * Send the signal task a message to install a timer.
       */
      void sendInstallTimerMsg(std::string name, int sigNo, 
			       unsigned long initSec, 
			       unsigned long initNsec,
			       unsigned long intervalSec, 
			       unsigned long intervalNsec,
			       SIGNALTASK_HANDLER_FN(*handler));
      
      /**
       * Respond to a message to install a timer.
       */
      void installTimer(Msg* msg);
      
      /**
       * Install timers of interest to us.
       */
      virtual void installTimers();
      
      /**
       * A method to install a signal
       */
      void sendInstallSignalMsg(int sigNo, SIGNALTASK_HANDLER_FN(*handler));
      
      /**
       * Respond to a message to install a signal.
       */
      void installSignal(Msg* msg);
      
      /**
       * Install timers of interest to us.
       */
      virtual void installSignals();
      
      /**
       * A message to enable/disable a timer.
       */
      void sendEnableTimerMsg(std::string name, bool enable);
      
      /**
       * Respond to a message to enable/disable a timer.
       */
      void enableTimer(Msg* msg);
      
      /**
       * A message to add/remove a signal handler
       */
      void sendAddHandlerMsg(std::string name, 
			     SIGNALTASK_HANDLER_FN(*handler),
			     bool add);
      
      /**
       * Respond to a message to add/remove a handler
       */
      void addHandler(Msg* msg);
      
    }; // End class GenericMasterTask
    
    // Template functions must be defined in the header file
    
    /*
     * Constructor
     */
    template<class Msg>
      GenericMasterTask<Msg>::GenericMasterTask() : 
      GenericTask<Msg>::GenericTask() 
      {
	DBPRINT(true, Debug::DEBUG7, "Inside MasterTask constructor");
	signal_ = 0;
      };
    
    /*
     * Detructor
     */
    template<class Msg>
      GenericMasterTask<Msg>::~GenericMasterTask() {};
    
    /**
     *  Method to process a message received on the Task message * queue
     */
    template<class Msg>
      void GenericMasterTask<Msg>::processTaskMsg(bool* stop)
      {
	Msg msg;
	
	GenericTask<Msg>::msgq_.readMsg(&msg);
	
	switch (msg.genericMsgType_) {
	case Msg::HEARTBEAT: // Is this a heartbeat request?
	  GenericTask<Msg>::respondToHeartBeat();
	  break;
	case Msg::RESTART:  // Is this a request to restart?
	  GenericTask<Msg>::restart();
	  break;
	case Msg::STOP:     // Did we receive a request to shut
	  // down?
	  *stop = true;
	  break;
	default:
	  {
	    GenericMasterTaskMsg::GenericMasterMsgType type = 
	      (GenericMasterTaskMsg::GenericMasterMsgType)msg.genericMsgType_;
	    switch(type) {
	    case Msg::INSTALL_TIMER:
	      installTimer(&msg);
	      break;
	    case Msg::INSTALL_SIGNAL:
	      DBPRINT(true, Debug::DEBUG31, "Got an install signal message");
	      installSignal(&msg);
	      break;
	    case Msg::ADD_HANDLER:
	      addHandler(&msg);
	      break;
	    case Msg::ENABLE_TIMER:
	      DBPRINT(true, Debug::DEBUG10, "Got an enable timer message");
	      enableTimer(&msg);
	      break;
	    default: // Else forward this message to the task-specific
	      // process method
              GenericMasterTask<Msg>::processMsg(&msg);
	      break;
	    }
	  }
	}
      };
    
    /**.......................................................................
     * Send a message to install a timer.
     */
    template<class Msg>
      void GenericMasterTask<Msg>::sendInstallTimerMsg(std::string name, int sigNo, 
						       unsigned long initSec, 
						       unsigned long initNsec,
						       unsigned long intervalSec, 
						       unsigned long intervalNsec,
						       SIGNALTASK_HANDLER_FN(*handler))
      {
	DBPRINT(true, DEBUG_SIGNAL, "Inside");
	
	Msg msg;
	
	msg.packInstallTimerMsg(name, sigNo, 
				initSec, initNsec, 
				intervalSec, intervalNsec, 
				handler);
	
	sendTaskMsg(&msg);
      }
    
    /**
     * Send a message to the master thread to install a timer.  The
     * initial delay is set to the interval in this version.
     */
    template<class Msg>
      void GenericMasterTask<Msg>::sendInstallTimerMsg(std::string name, int sigNo, 
						       unsigned long intervalSec, 
						       unsigned long intervalNsec,
						       SIGNALTASK_HANDLER_FN(*handler))
      {
	DBPRINT(true, DEBUG_SIGNAL, "Inside");
	
	Msg msg;
	
	msg.packInstallTimerMsg(name, sigNo, intervalSec, intervalNsec, handler);
	
	sendTaskMsg(&msg);
      }
    
    /**
     * Respond to a message to install a timer.
     */
    template<class Msg>
      void GenericMasterTask<Msg>::installTimer(Msg* msg)
      {
	// Install the timer.
	
	DBPRINT(true, Debug::DEBUG2, "Signal_ is: "
		<< (signal_ == 0 ? "NULL" : "not NULL"));
	
	if(signal_ != 0)
	  signal_->sendInstallTimerMsg(msg->genericMasterBody.installTimer.name,
				       msg->genericMasterBody.installTimer.sigNo,
				       msg->genericMasterBody.installTimer.initSec,
				       msg->genericMasterBody.installTimer.initNsec,
				       msg->genericMasterBody.installTimer.intervalSec,
				       msg->genericMasterBody.installTimer.intervalNsec,
				       msg->genericMasterBody.installTimer.handler);
      }
    
    /**
     * Install timers of interest to us.
     */
    template<class Msg>
      void GenericMasterTask<Msg>::installTimers() {};
    
    /**
     * Send a message to the master thread to install a signal.
     */
    template<class Msg>
      void GenericMasterTask<Msg>::sendInstallSignalMsg(int sigNo, 
							SIGNALTASK_HANDLER_FN(*handler))
      {
	Msg msg;
	
	DBPRINT(true, Debug::DEBUG2, "Signal_ is: "
		<< (signal_ == 0 ? "NULL" : "not NULL"));
	
	msg.packInstallSignalMsg(sigNo, handler);
	
	GenericMasterTask<Msg>::sendTaskMsg(&msg);
      }
    
    /**
     * Respond to a message to install a signal.
     */
    template<class Msg>
      void GenericMasterTask<Msg>::installSignal(Msg* msg)
      {
	// Install the signal.
	
	DBPRINT(true, Debug::DEBUG2, "Signal_ is: "
		<< (signal_ == 0 ? "NULL" : "not NULL"));
	
	
	if(signal_)
	  signal_->sendInstallSignalMsg(msg->genericMasterBody.installSignal.sigNo,
					msg->genericMasterBody.installSignal.handler);
      }
    
    /**
     * Install signals of interest to us.
     */
    template<class Msg>
      void GenericMasterTask<Msg>::installSignals() {};
    
    /**
     * Send a message to the master thread to install a timer.
     */
    template<class Msg>
      void GenericMasterTask<Msg>::sendEnableTimerMsg(std::string name, bool enable)
      {
	DBPRINT(true, Debug::DEBUG2, "Inside");
	
	Msg msg;
	
	msg.packEnableTimerMsg(name, enable);
	
	sendTaskMsg(&msg);
      }
    
    /**
     * Send a message to the master thread to install a timer.
     */
    template<class Msg>
      void GenericMasterTask<Msg>::sendAddHandlerMsg(std::string name, 
						     SIGNALTASK_HANDLER_FN(*handler), 
						     bool add)
      {
	DBPRINT(true, Debug::DEBUG2, "Inside, name = " << name
		<< " add = " << add);
	Msg msg;
	
	msg.packAddHandlerMsg(name, handler, add);
	
	sendTaskMsg(&msg);
      }
    
    /**
     * Send a message to the master thread to install a timer.
     */
    template<class Msg>
      void GenericMasterTask<Msg>::enableTimer(Msg* msg)
      {
	DBPRINT(true, Debug::DEBUG2, "Signal_ is: "
		<< (signal_ == 0 ? "NULL" : "not NULL"));
	
	if(signal_ != 0)
	  signal_->sendEnableTimerMsg(msg->genericMasterBody.enableTimer.name,
				      msg->genericMasterBody.enableTimer.enable);
      }
    
    /**
     * Send a message to the signal thread to add a handler.
     */
    template<class Msg>
      void GenericMasterTask<Msg>::addHandler(Msg* msg)
      {
	if(Debug::debugging(DEBUG_DELAY)) {
	  gcp::util::TimeVal timeVal;
	  timeVal.setToCurrentTime();
	  DBPRINT(true, DEBUG_DELAY, 
		  "Adding a handler to: " << msg->genericMasterBody.addHandler.name << " "
		  << timeVal.getSeconds());
	}
	
	DBPRINT(true, Debug::DEBUG2, "Signal_ is: "
		<< (signal_ == 0 ? "NULL" : "not NULL"));
	
	if(signal_ != 0)
	  signal_->sendAddHandlerMsg(msg->genericMasterBody.addHandler.name,
				     msg->genericMasterBody.addHandler.handler,
				     msg->genericMasterBody.addHandler.add);
      }
    
  } // End namespace util
} // End namespace gcp




#endif // End #ifndef GCP_UTIL_GENERICMASTERTASK_H
