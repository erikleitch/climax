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
  { "startup",      "addpath /home/gcpdaq/gcp_analysis;startup;", "s", "Startup string for matlab"},
  { "command",      "",            "s", "Matlab command to run"},
  {END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{

  if(Program::isDefault("command")) {
    ThrowError("You must specify a command to run");
  }

  Engine *ep;
  mxArray *T = NULL, *result = NULL;

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
  if (!(ep = engOpen("matlab -nodisplay -nosplash -nodesktop"))) {
    ThrowError("Can't start MATLAB engine");
  }

  engSetVisible(ep, 0);

  unsigned bufsize = 100000;
  char buffer[bufsize + 1];

  engEvalString(ep, Program::getStringParameter("startup").c_str());

  std::ostringstream os;
  os << Program::getStringParameter("command");
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
