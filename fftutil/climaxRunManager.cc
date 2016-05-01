#include <iostream>
#include <fstream>
#include <math.h>
#include <signal.h>
#include <sstream>

#include "gcp/fftutil/RunManager.h"

#include "gcp/program/Program.h"

#include "gcp/util/Exception.h"
#include "gcp/util/Logo.h"
#include "gcp/util/SignalTask.h"
#include "gcp/util/XtermManip.h"

using namespace gcp::program;
using namespace gcp::util;

//------------------------------------------------------------
// The only keywords we support are the file to run and an 
// option to turn the logo off
//------------------------------------------------------------

KeyTabEntry Program::keywords[] = {
  { "file",           "",  "s", "Run file to load"},
  { "logo",           "t", "b", "True to display logo"},
  {END_OF_KEYWORDS},
};

/**.......................................................................
 * Print usage information about this program
 */
void Program::initializeUsage() 
{
  std::ostringstream os;

  FORMATCOLOR(os, std::endl << "This program is controlled by a text-based run file, specified with the command-line 'file' option."
	      << std::endl << std::endl << "No other command-line options are allowed, so that runs with the same run file will be exactly reproducible." 
	      << std::endl << std::endl << "The code is intended to be self-documenting.  For more information, create a run file containing"
	      << std::endl << "the single line 'help;', and run this program on it, as in "
	      << std::endl << std::endl
	      << "\t\t\tclimaxRunManager file=myFile" << std::endl, "green");

  usage_ = os.str();
};

//=======================================================================
// Global variables used in signal handler
//=======================================================================

Logo logo_;
pthread_t mainId_;

/**.......................................................................
 * Signal handler for user interrupt
 */
static SIGNALTASK_HANDLER_FN(handler)
{
  //------------------------------------------------------------
  // Try and exit cleanly -- kill the main thread, and restore the
  // terminal if we were in the middle of displaying our logo
  //------------------------------------------------------------

  pthread_cancel(mainId_);

  if(!logo_.isDone())
    logo_.reset();

  //------------------------------------------------------------
  // Print message and exit
  //------------------------------------------------------------

  COUTCOLOR(std::endl << "User interrupt received -- CLIMAX exiting", "red");

  exit(1);
}

/**.......................................................................
 * Main -- parse a config file and execute its instructions
 */
int Program::main()
{
  mainId_ = pthread_self();

  //------------------------------------------------------------
  // Spawn a thread that will handle all signals, and install a
  // handler for SIGINT
  //------------------------------------------------------------

  SignalTask signalTask(true);
  signalTask.sendInstallSignalMsg(SIGINT, handler);

  //------------------------------------------------------------
  // The user must specify a file to run this program
  //------------------------------------------------------------

  if(!Program::hasValue("file")) {
    COUT(usage_);
    return 1;
  }

  //------------------------------------------------------------
  // Run climax
  //------------------------------------------------------------

  try {

    //------------------------------------------------------------
    // Display logo if requested
    //------------------------------------------------------------

    if(Program::getBooleanParameter("logo")) {
      logo_.display2();
    }

    //------------------------------------------------------------    
    // And run the file
    //------------------------------------------------------------    

    std::string file  = Program::getStringParameter("file");

    RunManager rm;
    rm.setRunFile(file);
    rm.run();
    rm.listModel();

  } catch(Exception& err) {
    COUT(err.what());
    return 1;
  }

  return 0;
}
