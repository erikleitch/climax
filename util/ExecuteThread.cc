#include "gcp/util/ExecuteThread.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ExecuteThread::ExecuteThread() :
  SpawnableTask<ExecuteThreadMsg>(true) 
{}

/**.......................................................................
 * Constructor with explicit cpu affinity specified (not supported on
 * certain platforms)
 */
ExecuteThread::ExecuteThread(unsigned cpu) :
  SpawnableTask<ExecuteThreadMsg>(true) 
{
  static bool first = true;

  try {
    spawnedThread_->setAffinity(cpu);
  } catch(Exception& err) {
    if(first) {
      COUT(err.what());
      first = false;
    }
  }
}

/**.......................................................................
 * Destructor.
 */
ExecuteThread::~ExecuteThread() {}

/**.......................................................................
 * Respond to a message to execute a command
 */
void ExecuteThread::processMsg(ExecuteThreadMsg* msg)
{
  switch (msg->type_) {
  case ExecuteThreadMsg::EXECUTE:
    msg->body_.execute_.fn_(msg->body_.execute_.args_);
    break;
  default:
    ThrowError("Unrecognized message type: " << msg->type_);
    break;
  }
}

void ExecuteThread::execute(EXECUTE_FN(*fn), void* args)
{
  ExecuteThreadMsg msg;
  msg.genericMsgType_ = GenericTaskMsg::TASK_SPECIFIC;

  msg.type_ = ExecuteThreadMsg::EXECUTE;

  msg.body_.execute_.fn_   = fn;
  msg.body_.execute_.args_ = args;

  sendTaskMsg(&msg);
}
