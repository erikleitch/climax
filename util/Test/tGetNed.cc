#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/CurlUtils.h"
#include "gcp/util/Exception.h"
#include "gcp/util/String.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "src",        "", "s", "Source to fetch"},
  { "print",        "f", "b", "True to print record"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  std::string src = Program::getStringParameter("src");

  std::ostringstream os;
  os << "http://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=" << src << "&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=NO";

  CurlUtils cu;
  std::string str = cu.getUrl(os.str());

  if(Program::getBooleanParameter("print")) {
    COUT("Got record: " << std::endl << str);
    return 0;
  }

  String nedStr(str);
  String test=nedStr.findNextInstanceOf(" ", false, "SOURCE LIST", "true");

  if(test.isEmpty()) {
    ThrowColorError("Source " << src << " not found", "red");
    return 1;
  }

  String srcStr = nedStr.findNextInstanceOf("/A>", true, "&nbsp", "true");
  srcStr.findNextInstanceOf(" ", false, ".", true, true);
  srcStr.advanceToNextNonWhitespaceChar();
  srcStr.findNextString();
  srcStr.advanceToNextNonWhitespaceChar();
  srcStr.findNextString();
  srcStr.advanceToNextNonWhitespaceChar();
  srcStr.findNextString();
  srcStr.advanceToNextNonWhitespaceChar();
  srcStr.findNextString();
  srcStr.advanceToNextNonWhitespaceChar();

  COUT("z = " << srcStr.remainder().toDouble());

  return 0;
}
