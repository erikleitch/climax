#ifndef GCP_UTIL_GENERICTASKMSG_H
#define GCP_UTIL_GENERICTASKMSG_H

/**
 * @file GenericTaskMsg.h
 * 
 * Tagged: Fri Nov 14 12:39:34 UTC 2003
 * 
 * @author Erik Leitch
 */
#define SIGNALTASK_HANDLER_FN(fn) void (fn)(int sigNo, void* args)
#define SIGNAL_NAME_LEN 10

namespace gcp {
  namespace util {
    
    /**
     * A class to encapsulate message types for a generic task.
     *
     * Classes which extend from this class should simply add
     * whatever members are required to process additional
     * task-specific messages, for instance a union of task-specific
     * messages.
     *
     * NB: There is no explicit constructor for this class, since
     * the compiler won't allow classes with constructors to be
     * included as members of unions.  This means that we cannot
     * construct unions using objects which inherit from
     * GenericTaskMsg unless they also don't have constructors.
     */
    class GenericTaskMsg {
      
    public:
      
      /**
       * Enumerate supported generic message types.
       */
      enum GenericMsgType {
	HEARTBEAT,      // A request to respond to the heartbeat
	STOP,           // A message to shutdown
	RESTART,        // A message to restart threads managed by this
	// task.
	TASK_SPECIFIC,  // A task-specific message
	LAST            // This should always come last!
      };
      
      /**
       * A type for this message
       */
      GenericMsgType genericMsgType_;
      
    }; // End class GenericTaskMsg
    
  } // End namespace util
} // End namespace gcp


#endif // End #ifndef 
