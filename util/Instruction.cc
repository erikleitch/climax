#include "gcp/util/Instruction.h"
#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Instruction::Instruction(INSTRUCTION(*fn), void* args, TimeVal timeToNext, TimeVal timeToRetry) 
{
  initialize();
  fn_ = fn;
  args_ = args;

  timeToNext_  = timeToNext;
  timeToRetry_ = timeToRetry;
}

/**.......................................................................
 * Constructor.
 */
Instruction::Instruction(INSTRUCTION(*fn), void* args, TimeVal timeToNext)
{
  initialize();
  fn_ = fn;
  args_ = args;

  timeToNext_  = timeToNext;

  // Default retry to 0.5 seconds

  timeToRetry_.setTime(0, 500000, 0);
}

/**.......................................................................
 * Constructor.
 */
Instruction::Instruction(INSTRUCTION(*fn), void* args)
{
  initialize();
  fn_ = fn;
  args_ = args;

  // Default time to next command to 0.5 seconds

  timeToNext_.setTime(0, 500000, 0);

  // Default retry to 0.5 seconds

  timeToRetry_.setTime(0, 500000, 0);
}

/**.......................................................................
 * Constructor.
 */
Instruction::Instruction()
{
  initialize();
}

void Instruction::initialize()
{
  fn_ = 0;
  args_ = 0;
  timeToNext_.setTime(0,0,0);
  timeToRetry_.setTime(0,0,0);
}

/**.......................................................................
 * Destructor.
 */
Instruction::~Instruction() {}

/**.......................................................................
 * Execute this instruction
 */
Instruction::State Instruction::execute(TimeVal& timeOut)
{
  // Execute the command.

  State state = fn_(args_);

  // If we should advance to the next instruction, set the timeout to
  // the time until the next command should be tried.
  //
  // If we should not, set the timeout to the time until this command
  // should be retried.

  if(state==DONE) 
    timeOut = timeToNext_;
  else if(state==AGAIN)
    timeOut = timeToRetry_;

  return state;
}

void Instruction::setTimeToNext(TimeVal& tVal)
{
  timeToNext_ = tVal;
}

void Instruction::setTimeToRetry(TimeVal& tVal)
{
  timeToRetry_ = tVal;
}
