#include <iostream>
#include <fstream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Exception.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "file",        "logo.txt", "s", "Logo file"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

#define COLOR(color)\
  macroXtm.fg("black") << macroXtm.bg(color) << macroXtm.textMode("bold")

#define DEFAULT \
  macroXtm.bg("default") << macroXtm.fg("default") << macroXtm.textMode("normal")

std::string getLogo();

unsigned loadLogo(std::string str, unsigned iStart, unsigned iStop, 
		  std::string borderColor, std::string contentColor, std::string wipeColor)
{
  std::istringstream istr(str);

  std::ostringstream os;

  std::string s;

  gcp::util::XtermManip macroXtm;

  unsigned nline=0;
  while(getline(istr, s)) {
    nline++;
    os.str("");

    os << COLOR(wipeColor) << "  " << DEFAULT;

    for(unsigned i=0; i < s.size(); i++) {

      if(i < iStart)
	os << COLOR(wipeColor) << "  " << DEFAULT;
      else if(i > iStop || s[i] == ' ')
	os << COLOR(borderColor) << "  " << DEFAULT;
      else
	os << COLOR(contentColor) << "  " << DEFAULT;
    }

    COUT(os.str());
  }

  COUT(DEFAULT);

  return nline;
}

void getDimensions(std::string fileName, unsigned& nRow, unsigned& nCol)
{
  std::istringstream istr(fileName);
  std::string s;

  gcp::util::XtermManip macroXtm;

  unsigned nline=0;
  unsigned len=0;

  while(getline(istr, s)) {
    nline++;
    len = s.size() > len ? s.size() : len;
  }

  nCol = len;
  nRow = nline;

  return;
}

void flash(std::string file, unsigned ntime)
{
  gcp::util::XtermManip macroXtm;
  struct timespec ts;
  ts.tv_sec  = 0;
  ts.tv_nsec = 50000000;

  unsigned nline;

  unsigned nRowLogo, nColLogo;
  getDimensions(file, nRowLogo, nColLogo);

  for(unsigned i=0; i < ntime; i++) {
    nline = loadLogo(file, 0, nColLogo, "blue", "magenta", "default");

    macroXtm.moveCursorUp(nline+1);

    nanosleep(&ts, 0);
    loadLogo(file, 0, nColLogo, "magenta", "blue", "default");

    macroXtm.moveCursorUp(nline+1);
    nanosleep(&ts, 0);
  }
}

void stripe(std::string file, unsigned ntime)
{
  gcp::util::XtermManip macroXtm;
  struct timespec ts;
  ts.tv_sec  = 0;
  ts.tv_nsec = 5000000;

  unsigned nRowLogo, nColLogo;
  getDimensions(file, nRowLogo, nColLogo);

  unsigned nline=0;
  for(unsigned i=0; i <= nColLogo; i++) {
    nline = loadLogo(file, 0, i, "blue", "magenta", "default");

    macroXtm.moveCursorUp(nline+1);

    nanosleep(&ts, 0);
  }
}

void destripe(std::string file, unsigned ntime)
{
  gcp::util::XtermManip macroXtm;
  struct timespec ts;
  ts.tv_sec  = 0;
  ts.tv_nsec = 5000000;

  unsigned nRowLogo, nColLogo;
  getDimensions(file, nRowLogo, nColLogo);

  unsigned nline=0;
  for(unsigned i=0; i <= nColLogo; i++) {
    nline = loadLogo(file, i, nColLogo, "blue", "magenta", "default");

    if(i < nColLogo) {
      macroXtm.moveCursorUp(nline+1);
    }

    nanosleep(&ts, 0);
  }
}

int Program::main()
{
  gcp::util::XtermManip macroXtm;
  unsigned nCol = macroXtm.getNCol();
  unsigned nRow = macroXtm.getNRow();

  std::string file = getLogo();

  unsigned nRowLogo, nColLogo;
  getDimensions(file, nRowLogo, nColLogo);

  macroXtm.hideCursor();

  if(nCol >= 2*nColLogo && nRow >= nRowLogo+2) {
    stripe(file, 4);
    flash(file, 4);
    destripe(file, 4);
  }

  macroXtm.showCursor();

  return 0;
}

std::string getLogo()
{
  std::ostringstream os;
  
  os << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl
     << "x                                    x" << std::endl
     << "x               x     x              x" << std::endl
     << "x             x xx   xx   xx         x" << std::endl
     << "x       x     x x x x x  x  x        x" << std::endl
     << "x  xxx  x     x x  x  x x    x x   x x" << std::endl
     << "x x   x x     x x     x xxxxxx  x x  x" << std::endl
     << "x x     x     x x     x x    x   x   x" << std::endl
     << "x x   x x     x x     x x    x  x x  x" << std::endl
     << "x  xxx  xxxxx x x     x x    x x   x x" << std::endl
     << "x                                    x" << std::endl
     << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;

  return os.str();
}

