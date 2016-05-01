#include <iostream>

#include "gcp/program/Program.h"
#include "gcp/util/PythonGenerator.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

void Program::initializeUsage() {};

KeyTabEntry Program::keywords[] = {
  {    "dir",      "",     "s", "Input/output directory"},
  {   "file",      "",     "s", "Input file"},
  { "prefix",      "",     "s", "Prefix of the output cc/dlm file"},
  { "suffix", "pycc",     "s", "Suffix of the output cc file"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

int Program::main()
{
  PythonGenerator pythonGen(Program::getStringParameter("file"), Program::getStringParameter("dir"));

  if(Program::hasValue("prefix"))
    pythonGen.setOutputPrefix(Program::getStringParameter("prefix"));

  if(Program::hasValue("suffix"))
    pythonGen.setOutputCcSuffix(Program::getStringParameter("suffix"));

  pythonGen.outputCcFile();

  return 0;
}
