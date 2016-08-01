#include <iostream>
#include <iomanip>

#include <cmath>

#include <errno.h>
#include <unistd.h>

#include "gcp/program/Program.h"

#include "gcp/util/String.h"
#include "gcp/util/Exception.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "debuglevel", "0", "i", "Debug level"},
  { "val", "0", "i", "Value to print"},
  { "n",   "0", "i", "nchar"},
  { "str", " ", "s", "String to print"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  String var("my home directory is ${HOME} and my climax tools dir = ${CLIMAX_TOOLS}/${EMLTEST}");

  var.expandVar();

  COUT("var is now: " << var);

  return 0;
  
  String uname("~/this/is/a/test/~eml");

  uname.expandTilde();

  COUT("uname is now: " << uname);
  COUT("Getloging returns: " << getlogin());
  return 0;
}

int oldTests(Program* prog)
{
    String tStr("-max");

  COUT((tStr == "-max"));
  return 0;

  tStr = "a";
  try {
    tStr.toInt();
  } catch(Exception& err) {
    COUT(err.what());
  }

  try {
    tStr.toDouble();
  } catch(Exception& err) {
    COUT(err.what());
  }

  tStr = "1.02";

  try {
    COUT("Val = " << tStr.toInt());
  } catch(Exception& err) {
    COUT(err.what());
  }

  try {
    COUT("Val = " << tStr.toDouble());
  } catch(Exception& err) {
    COUT(err.what());
  }

  tStr = "1.02x";

  try {
    COUT("Val = " << tStr.toInt());
  } catch(Exception& err) {
    COUT(err.what());
  }

  try {
    COUT("Val = " << tStr.toDouble());
  } catch(Exception& err) {
    COUT(err.what());
  }

  return 0;

  std::vector<double> vals = tStr.parseRange();
  for(unsigned i=0; i < vals.size(); i++) {
    COUT("val = " << vals[i]);
  }
  return 0;

  COUT("str = '" << tStr.getNextNchars(prog->getIntegerParameter("n")));
  return 0;

  tStr = prog->getStringParameter("str");
  COUT("tStr = " << tStr);
  COUT("nnn = '" << tStr.findNextNonNumericString() << "'");
  return 0;

  COUT(String::formatHumanReadableInteger(prog->getIntegerParameter("val")));
  return 0;

  String testStr3("10^-5");
  testStr3.superscriptNumbersForPgplot();
  COUT(testStr3);

  return 0;

  String testStr2("10^2 Mpc x 10^3kg");
  testStr2.superscriptNumbersForPgplot();

  COUT("str = " << testStr2);

  return 0;

  String testStr1("3e6Msolar");

  testStr1.replace("Msolar", "");

  COUT("string is now " << testStr1);

  return 0;

  String str2("# beta_x beta_y beta_t0 beta_rc beta_eta beta_phi ptsrc0_x ptsrc0_y ptsrc0_flux _LOGP");

  COUT(str2.findNextInstanceOf("", false, " ", true, false));

  return 0;

  String val("0.01:-5:5'");
  val.advanceToNextNonWhitespaceChar();
  String tok1 = val.findNextInstanceOf(" ", false, ":", true, true);
  String mean, min, max, units;

  // If the string still contains another ':', then we assume we
  // are parsing a 'mean:min:max' specification, else just 'min:max'

  if(val.remainder().contains(":")) {
    mean = tok1;
    val.advanceToNextNonWhitespaceChar();
    min = val.findNextInstanceOf(" ", false, ":", true, true);
  } else {
    min = tok1;
  }

  val.advanceToNextNonWhitespaceChar();
  max   = val.findNextNumericString();
  units = val.remainder();
  units.strip(' ');

  COUT("mean = '" << mean << "' min = '" << min << "' max = '" << max << "' units = '" << units << "'");
  return 0;

  String ants("sza,bima");
  COUT(ants.remainder().contains(","));
  COUT(ants.findNextInstanceOf("", false, ",", true, false));
  COUT(ants.remainder().contains(","));
  COUT(ants.findNextInstanceOf(",", true, " ", false));
  COUT(ants.remainder().contains(","));

  return 0;

  String range("400:300,400:500");
  COUT(range.contains(","));
  COUT(range.findNextInstanceOf("", false, ",", true, false));
  COUT(range.contains(","));
  COUT(range.findNextInstanceOf(",", true, " ", false));
  COUT(range.contains(","));

  return 0;

  String uname("~eml/this/is/a/test/~eml");

  uname.expandTilde();

  COUT("uname is now: " << uname);

  return 0;

  String line(" adddataset name=test");

  // Get to the first non-whitespace char

  line.advanceToNextNonWhitespaceChar();

  // Skip the addmodel token

  COUT("Next string = '" << line.findNextString() << "'");
  COUT("Next string = '" << line.findNextString() << "'");

  return 0;

  String str1("23.5 +- 0.1mJy");

  COUT(str1.findNextInstanceOf("", false, "+-", true, true));
  str1.advanceToNextNonWhitespaceChar();
  COUT(str1.findNextNumericString());
  COUT(str1.remainder());

  str1 = "23.5 +- 0.1 mJy";

  COUT(str1.findNextInstanceOf("", false, "+-", true, true));
  str1.advanceToNextNonWhitespaceChar();
  COUT(str1.findNextNumericString());
  units = str1.remainder();
  units.strip(' ');
  COUT(units);
  COUT(units.isEmpty() << " " << units.size() << " '" << units[0] << "'");

  str1 = "23.5+-0.1";

  COUT(str1.findNextInstanceOf("", false, "+-", true, true));
  str1.advanceToNextNonWhitespaceChar();
  COUT(str1.findNextNumericString());
  units = str1.remainder();
  COUT(units.isEmpty() << " " << units.size() << " '" << units[0] << "'");
  units.strip(' ');
  COUT(units);
  COUT(units.isEmpty() << " " << units.size() << " '" << units[0] << "'");

  return 0;

  String str("-23:34:01.1");

  COUT(str.findNextInstanceOf("", false, ":", true, true));

  COUT(str.contains(":"));

  COUT(str.findNextInstanceOf("", false, ":", true, true));

  COUT(str.contains(":"));

  return 0;
}
