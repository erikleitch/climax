#include "gcp/util/Sort.h"
#include "gcp/util/Exception.h"
#include "gcp/util/CoProc.h"

#include<iostream>

#include <stdio.h>

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Sort::Sort() {}


/**.......................................................................
 * Destructor.
 */
Sort::~Sort() {}

std::vector<string> Sort::sort(vector<string>& entries)
{
  CoProc proc("sort");

  FILE* stdIn  = proc.stdIn()->writeFp();
  FILE* stdOut = proc.stdOut()->readFp();

  unsigned size=0;
  for(unsigned i=0; i < entries.size(); i++) {
    fprintf(stdIn, "%s\n", entries[i].c_str());
    size = entries[i].size() > size ? entries[i].size() : size;
  }

  fclose(stdIn);

  vector<string> sortedEntries;

  char c;
  std::ostringstream os;
  while((c = (char)fgetc(stdOut)) != EOF) {
    if(c == '\n') {
      sortedEntries.push_back(os.str());
      os.str("");
    } else {
      os << c;
    }
  }

  return sortedEntries;
}
