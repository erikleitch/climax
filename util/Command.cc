#include "gcp/util/Command.h"
#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Command::Command() 
{
  initialize();

  name_   = "";
  active_ = false;
  
  doneHandler_ = 0;
  doneArgs_ = 0;

  failedHandler_ = 0;
  failedArgs_ = 0;

  timeOut_.setTime(0,0);
  lastTime_.setToCurrentTime();
  currTime_ = lastTime_;
}

/**.......................................................................
 * Destructor.
 */
Command::~Command() {}

/**.......................................................................
 * (Re)Initialize this command
 */
void Command::initialize()
{
  instructions_.resize(0);
  reset();
}

/**.......................................................................
 * (Re)Initialize this command
 */
void Command::reset()
{
  active_  = false;
  nextInstruction_ = instructions_.begin();
}

/**.......................................................................
 * (Re)Initialize this command
 */
void Command::restart(struct timeval* timeOut)
{
  // Set the timeout to 1 second before this command will be run

  timeOut_.setSeconds(1);
  timeOut = timeOut_.timeVal();
  
  reset();
  active_  = true;
}

/**.......................................................................
 * Insert another instruction into our list
 */
Instruction* Command::insert(Instruction instruction)
{
  instructions_.push_back(instruction);
  reset();
  return &instructions_[instructions_.size()-1];
}

/**.......................................................................
 * Insert another command
 */
void Command::insert(Command& command)
{
  for(unsigned iInstr=0; iInstr < command.instructions_.size(); iInstr++)
    instructions_.push_back(command.instructions_[iInstr]);
}

/**.......................................................................
 * Execute the next instruction
 */
void Command::executeNextInstruction(TimeVal& timeOut, bool setToValue)
{
  // If this command is now complete, call our wrap-up function
  
  if(isComplete())
    registerCompletion();

  // Else execute the next instruction, if this command is still active.  

  else if(active_) {

    // Check if enough time has elapsed to execute the next instruction

    currTime_.setToCurrentTime();
    diff_ = currTime_ - lastTime_;

    // If enough time has elapsed, execute the next command

    if(diff_ >= timeOut_) {

      COUT("currTime = " << currTime_ << ", lastTime = " << lastTime_ 
      	   << ", difference (" << diff_ << ") is >= timeout (" << timeOut_ << "): executing next command");
      
      Instruction::State state = nextInstruction_->execute(timeOut_);

      if(state==Instruction::DONE) 
	nextInstruction_++;
      else if(state==Instruction::FAILED)
	registerFailure();

      // Else set the timeout for this command to the requested
      // timeout, less the amount of time that has actually elapsed


    } else {

      COUT("currTime = " << currTime_ << ", lastTime = " << lastTime_ 
		 << ", difference (" << diff_ << ") is NOT >= timeout (" << timeOut_ << "): NOT executing next command");

      timeOut_ -= diff_;
    }
  
    // And set the time to the current time
    
    lastTime_.setToCurrentTime();

    // Finally, set the passed timeout to the smaller of the two timeouts

    if(setToValue)
      timeOut = timeOut_;
    else {
      if(timeOut_ < timeOut)
	timeOut = timeOut_;
    }

    CTOUT("Setting timeout to: " << timeOut);
  }
}

/**.......................................................................
 * Run this command
 */
void Command::run()
{
  reset();
  executeNextInstruction(timeOut_, true);
}

/**.......................................................................
 * Install a handler to be called when this command is complete
 */
void Command::
installDoneHandler(COMMAND_DONE_HANDLER(*handler), void* args)
{
  doneHandler_ = handler;
  doneArgs_ = args;
}

/**.......................................................................
 * Install a handler to be called if this command fails
 */
void Command::
installFailedHandler(COMMAND_FAILED_HANDLER(*handler), void* args)
{
  failedHandler_ = handler;
  failedArgs_ = args;
}

/**.......................................................................
 * Return true if this command is complete
 */
bool Command::isComplete()
{
  return nextInstruction_ == instructions_.end();
}

/**.......................................................................
 * Register that this command has completed
 */
void Command::registerCompletion()
{
  // Mark this command as inactive

  active_ = false;

  if(doneHandler_ != 0)
    doneHandler_(doneArgs_);
}

/**.......................................................................
 * Register that this command has failed
 */
void Command::registerFailure()
{
  // Mark this command as inactive

  active_ = false;

  if(failedHandler_ != 0)
    failedHandler_(failedArgs_);
}

TimeVal& Command::timeOut()
{
  return timeOut_;
}

bool Command::active()
{
  return active_;
}
