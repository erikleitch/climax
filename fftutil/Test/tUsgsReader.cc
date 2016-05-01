#include <iostream>

#include "gcp/program/Program.h"
#include "gcp/pgutil/UsgsReader.h"

#include "gcp/util/Exception.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

void Program::initializeUsage() {};

KeyTabEntry Program::keywords[] = {
  { "dir",    "/Users/eml/projects/usgs/USGS/", "s", "Directory prefix"},
  { "prefix", "29124458",                       "s", "Prefix to read"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

int Program::main()
{
  UsgsReader reader;

  reader.setTo(Program::getStringParameter("dir"), Program::getStringParameter("prefix"));
  reader.parseHeaderInfo();
  reader.readImageData();

  reader.display();

  return 0;
}

