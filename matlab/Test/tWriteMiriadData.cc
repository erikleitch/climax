#include <iostream>
#include <sstream>
#include <cmath>

#include "gcp/program/Program.h"
#include "gcp/util/Exception.h"

#include "engine.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "start",        "", "s", "Start date: dd-mmm-yy:hh:mm:ss"},
  { "stop",         "", "s", "Stop  date: dd-mmm-yy:hh:mm:ss"},
  { "file",         "", "s", "Output file name"},
  { "freq",         "", "i", "Frequency (30 or 90)"},
  { "interval",     "1.0", "f", "Interval, in hours, into which to break the data for processing"},
  { "startup",      "addpath /home/gcpdaq/gcp_analysis;startup;", "s", "Startup string for matlab"},
  {END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{

  if(Program::isDefault("start") 
     || Program::isDefault("stop") 
     || Program::isDefault("file")) {
    ThrowError("You must specify a start and stop time, as well as an output file name");
  }

  Engine *ep;
  mxArray *T = NULL, *result = NULL;
  double time[10] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

  /*
   * Start the MATLAB engine locally by executing the string
   * "matlab"
   *
   * To start the session on a remote host, use the name of
   * the host as the string rather than \0
   *
   * For more complicated cases, use any string with whitespace,
   * and that string will be executed literally to start MATLAB
   */
  if (!(ep = engOpen("\0"))) {
    ThrowError("Can't start MATLAB engine");
  }

  engSetVisible(ep, 0);

  unsigned bufsize = 100000;
  char buffer[bufsize + 1];

  engEvalString(ep, Program::getStringParameter("startup").c_str());

  std::ostringstream os;
  os << "write2miriad(" 
     << "'" << Program::getStringParameter("start") << "', "
     << "'" << Program::getStringParameter("stop") << "', "
     << "'" << Program::getStringParameter("file") << "', "
     << "'new', "
     << Program::getDoubleParameter("interval") << ", "
     << Program::getDoubleParameter("freq") << ")";

  COUT("Evaluating string: " << os.str());

  // Evaluate it

  buffer[bufsize] = '\0';
  engOutputBuffer(ep, buffer, bufsize);
  engEvalString(ep, os.str().c_str());

  printf("%s", buffer+2);

  // We're done! Free memory, close MATLAB engine and exit.

  printf("Done!\n");
  engClose(ep);
	
  return 0;
}
