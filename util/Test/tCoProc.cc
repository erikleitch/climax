#include <iostream>

#include "gcp/program/Program.h"
#include "gcp/util/Exception.h"
#include "gcp/util/CoProc.h"

using namespace std;
using namespace gcp::program;
using namespace gcp::util;

void Program::initializeUsage() {};

KeyTabEntry Program::keywords[] = {
  { "command", "\0",  "s", "A test command"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

int Program::main()
{
#if 1
  FILE* stdIn=0;
  FILE* stdOut=0;
  FILE* stdErr=0;

  CoProc proc(Program::getStringParameter("command"));

  stdIn  = proc.stdIn()->writeFp();
  stdOut = proc.stdOut()->readFp();
  stdErr = proc.stdErr()->readFp();

  //  fprintf(stdIn, "hello\nworld\nablus\n");
  fclose(stdIn);

  char buff[80];

  COUT("Got this from stdout:");

  int c;
  std::ostringstream os;
  while((c=fgetc(stdOut)) != EOF) {
    os << (char)c;
    COUT((char)c << " (" << c << ")");
  }

  COUT(os.str());

  COUT("Got this from stderr:");
  os.str("");
  while((c=fgetc(stdErr)) != EOF) {
    os << (char)c;
    COUT((char)c << " (" << c << ")");
  }

  COUT(os.str());
  
#else

  std::vector<std::string> args 
    = CoProc::split(Program::getParameter("command"));

  for(unsigned i=0; i < args.size(); i++)
    COUT("args[" << i << "] = " << args[i]);
#endif

  return 0;
}
