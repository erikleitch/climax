// $Id: Instruction.h,v 1.1 2012/03/15 01:19:32 eml Exp $

#ifndef GCP_UTIL_INSTRUCTION_H
#define GCP_UTIL_INSTRUCTION_H

/**
 * @file Instruction.h
 * 
 * Tagged: Mon Oct 17 12:17:35 PDT 2005
 * 
 * @version: $Revision: 1.1 $, $Date: 2012/03/15 01:19:32 $
 * 
 * @author BICEP Software Development
 */
#include "gcp/util/TimeVal.h"

#define INSTRUCTION(fn) gcp::util::Instruction::State (fn)(void* args)

namespace gcp {
  namespace util {

    // A class for managing a single instruction of a command

    class Instruction {
    public:

      enum State {
	DONE,
	AGAIN,
	FAILED
      };

      /**
       * Constructor.
       */
      Instruction(INSTRUCTION(*fn), void* args, TimeVal timeToNext, TimeVal timeToRetry);
      Instruction(INSTRUCTION(*fn), void* args, TimeVal timeToNext);
      Instruction(INSTRUCTION(*fn), void* args=0);
      Instruction();

      /**
       * Destructor.
       */
      virtual ~Instruction();

      // Execute this instruction.  Method should return true if the
      // instruction has successfully executed.

      virtual State execute(TimeVal& timeOut);
      
      void initialize();

      void setTimeToNext(TimeVal& tVal);
      void setTimeToRetry(TimeVal& tVal);

    private:
      
      // A function to call 

      INSTRUCTION(*fn_);

      // The arguments to call it with

      void* args_;

      // The time we should wait before executing another instruction
      // after this one

      TimeVal timeToNext_;
      TimeVal timeToRetry_;

    }; // End class Instruction

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_INSTRUCTION_H
