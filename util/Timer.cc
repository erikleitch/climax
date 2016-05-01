#include "gcp/util/Timer.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Timer::Timer() 
{
  integratedSec_ = 0.0;
}

/**.......................................................................
 * Destructor.
 */
Timer::~Timer() {}

void Timer::start()
{
  start_.setToCurrentTime();
}

void Timer::stop()
{
  stop_.setToCurrentTime();
  diff_ = stop_ - start_;
  integratedSec_ += deltaInSeconds();
}

double Timer::integratedElapsedSeconds()
{
  return integratedSec_;
}

double Timer::deltaInSeconds()
{
  return diff_.getTimeInSeconds();
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, const Timer& timer)
{
  return operator<<(os, (Timer&)timer);
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, Timer& timer)
{
  os << "Elapsed time = " << timer.deltaInSeconds() << " s";
  return os;
}
