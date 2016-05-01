#include "gcp/util/Directives.h"
#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"

using namespace std;
using namespace gcp::util;

/**.......................................................................
 * Constructs an Error with a detailed message.
 * 
 * @param str string containing message.
 * @param filename where exception originated.
 * @param lineNumber where exception occurred.
 */
Exception::Exception(string str, const char* fileName, const int lineNumber, bool doReport, bool printHelp)
{
  initialize(str, doReport, printHelp);
}

void Exception::initialize(std::string message, bool doReport, bool printHelp)
{
  printHelp_ = printHelp;
  message_ = message;

  if(doReport) 
    report();
} 

/**.......................................................................
 * Constructs an Error with a detailed message.
 * 
 * @param os ostringstream containing message
 * @param filename where exception originated.
 * @param lineNumber exception occurred on.
 */
Exception::Exception(ostringstream& os, const char * fileName, 
		     const int lineNumber, bool doReport, bool printHelp)
{
  initialize(os.str(), false, printHelp);
}

/**.......................................................................
 * Constructs an Error with a detailed message.  
 * 
 * @param ls LogStream containing message 
 * @param fileName where exception originated.  
 * @param lineNumber where exception occurred.
 */
Exception::Exception(LogStream& ls, 
		     const char* fileName, const int lineNumber, bool doReport, bool printHelp)
{
  initialize(ls.getMessage(), false, printHelp);;
}

/** 
 * Constructor with a log stream.
 */
Exception::Exception(LogStream* ls, 
	  const char* fileName, const int lineNumber, bool doReport, bool printHelp)
{
  initialize(ls->getMessage(), false, printHelp);
}

/**.......................................................................
 * Destructor
 */
Exception::~Exception() {}

bool Exception::printHelp() 
{
  return printHelp_;
}

void Exception::setPrintHelp(bool print)
{
  printHelp_ = print;
}

