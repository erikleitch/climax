#include "gcp/util/Debug.h"

using namespace gcp::util;

// Create the static instance

//Debug::Debug();

// Initialize static variables

Debug::Level Debug::level_ = DEBUGANY;
Mutex Debug::mutex_;

/**.......................................................................
 * Constructor.
 */
Debug::Debug() {};

/**.......................................................................
 * Public method to change the debugging level
 */
void Debug::setLevel(Level level)
{
  level_ = level;
}

/**.......................................................................
 * Public method to change the debugging level
 */
void Debug::addLevel(Level level)
{
  unsigned lmask = static_cast<unsigned>(level_);
  unsigned rmask = static_cast<unsigned>(level);
  level_ = static_cast<Level>(lmask | rmask);
}

/**.......................................................................
 * Public method to change the debugging level
 */
void Debug::remLevel(Level level)
{
  unsigned lmask = static_cast<unsigned>(level_);
  unsigned rmask = static_cast<unsigned>(level);
  level_ = static_cast<Level>(lmask & ~rmask);
}

/**.......................................................................
 * Public method to change the debugging level.
 */
void Debug::setLevel(unsigned int level)
{
  level_ = (Level) level;
}

/**.......................................................................
 * Public method to query the debugging state
 */
bool Debug::debugging(Level level)
{
  return level_ & level;
}

void Debug::lock()
{
  mutex_.lock();
}

void Debug::unlock()
{
  mutex_.unlock();
}

