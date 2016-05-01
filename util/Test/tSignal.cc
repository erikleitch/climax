#include <iostream>
#include <iomanip>

#include <cmath>
#include <pthread.h>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Signal.h"
#include "gcp/util/SignalTask.h"
#include "gcp/util/GenericTaskMsg.h"
#include "gcp/util/Thread.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "rad",        "0", "d", "Radians"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

static SIGNALTASK_HANDLER_FN(taskHandler)
{
  std::cout << pthread_self() << " Caught a signal" << std::endl;
}

static void handler(int sigNo)
{
  std::cout << pthread_self() << " Caught a signal" << std::endl;
}

static THREAD_START(startThread)
{
  cout << pthread_self() << " Inside startThread" << endl;

  //  signal(SIGINT, handler);

  while(true)
    select(0, NULL, NULL, NULL, NULL);

  cout << pthread_self() << " Leaving startThread" << endl;
}

pthread_t createThread()
{
  pthread_t id;
  pthread_create(&id, NULL, &startThread, NULL);
  return id;
}

int Program::main()
{
  //  Signal::blockAllSignals();
  //  signal(SIGINT, handler);

  //  COUT("Creating thread: " << createThread());

  SignalTask signalTask(true);

  signalTask.sendInstallSignalMsg(SIGINT, taskHandler);

  COUT("About to sleep");
  unsigned int ret = sleep(5);

  COUT(pthread_self() << " Exiting with ret = " << ret);

  return 0;
}

