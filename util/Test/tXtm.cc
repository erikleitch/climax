#include <iostream>
#include <iomanip>
#include <termios.h>
#include <sys/ioctl.h>
#include <stdio.h>
#include <cmath>
#include <unistd.h>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Debug.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "rad",        "0", "d", "Radians"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  XtermManip xtm;

  std::string fg = xtm.getFg();
  std::string bg = xtm.getBg();

  //  xtm.setBg("teal");
  //  xtm.setFg("olive");

  std::cout << "\e[" << "38;2;80;80;80" << "m";
  //  std::cout << "\e[" << "95" << "m";

  std::cout << "this is a test" << std::endl;

  sleep(5);
  xtm.setHexBg(bg);
  xtm.setHexFg(fg);

  return 0;
}
