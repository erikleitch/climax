// $Id: Command.h,v 1.1 2012/03/15 01:19:32 eml Exp $

#ifndef GCP_UTIL_COMMAND_H
#define GCP_UTIL_COMMAND_H

/**
 * @file Command.h
 * 
 * Tagged: Mon Oct 17 10:54:36 PDT 2005
 * 
 * @version: $Revision: 1.1 $, $Date: 2012/03/15 01:19:32 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Instruction.h"

#include <string>
#include <vector>

// A function to be called when the command has completed

#define COMMAND_DONE_HANDLER(fn) void (fn)(void* args)

// A function to be called if the command has failed

#define COMMAND_FAILED_HANDLER(fn) void (fn)(void* args)

namespace gcp {
  namespace util {

    // A class for managing canned commands that can be executed by
    // tasks via their message queues.

    class Command {
    public:

      /**
       * Constructor.
       */
      Command();

      /**
       * Destructor.
       */
      virtual ~Command();

      /**
       * Install a handler to be called when this command is complete
       */
      void installDoneHandler(COMMAND_DONE_HANDLER(*handler), 
			      void* args=0);
      
      /**
       * Install a handler to be called if this command fails
       */
      void installFailedHandler(COMMAND_FAILED_HANDLER(*handler), 
				void* args=0);

      // Insert another instruction into our list

      Instruction* insert(Instruction instruction);

      // Insert another command into our list

      void insert(Command& command);

      void run();

      TimeVal& timeOut();

      bool active();

      void restart(struct timeval* timeOut);

      // Execute the next instruction of this command
      
      void executeNextInstruction(TimeVal& timeOut, bool setToValue);
      
    protected:

      bool active_;

      COMMAND_DONE_HANDLER(*doneHandler_);
      void* doneArgs_;
      
      COMMAND_FAILED_HANDLER(*failedHandler_);
      void* failedArgs_;

      // A name for this command

      std::string name_;

      // A vector of instructions

      std::vector<Instruction> instructions_;

      // The current instruction
      
      std::vector<Instruction>::iterator nextInstruction_;
      
      TimeVal lastTime_;
      TimeVal currTime_;
      TimeVal diff_;
      TimeVal timeOut_;

      void initialize();
      void reset();
      bool isComplete();
      
      // Register that this command has completed

      void registerCompletion();

      // Register that this command has failed

      void registerFailure();

    }; // End class Command

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_COMMAND_H
