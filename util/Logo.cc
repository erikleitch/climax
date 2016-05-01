#include "gcp/util/Exception.h"
#include "gcp/util/Logo.h"
#include "gcp/util/XtermManip.h"

#include <sstream>
#include <stdlib.h>

#define LOGO_COL1 "magenta" 
#define LOGO_COL2 "blue"

using namespace std;

using namespace gcp::util;

#define COLOR(color)							\
  macroXtm_.fg("black") << macroXtm_.bg(color) << macroXtm_.textMode("bold")

#define DEFAULT								\
  macroXtm_.bg("default") << macroXtm_.fg("default") << macroXtm_.textMode("normal")

/**.......................................................................
 * Constructor.
 */
Logo::Logo() 
{
  setDone(true);
  nline_ = 0;
  nLinePrinted_ = 0;
}

/**.......................................................................
 * Destructor.
 */
Logo::~Logo() {}

unsigned Logo::loadLogo(std::string str, unsigned iStart, unsigned iStop, 
			std::string borderColor, std::string contentColor, std::string wipeColor, 
			unsigned nBuffer)
{
  std::istringstream istr(str);
  std::ostringstream os;
  std::string s;

  gcp::util::XtermManip macroXtm;

  unsigned nline=0;
  nLinePrinted_ = 0;

  while(getline(istr, s)) {

    os.str("");

    os << COLOR(wipeColor);
    for(unsigned i=0; i < nBuffer; i++)
      os << " ";
    os << DEFAULT;

    for(unsigned i=0; i < s.size(); i++) {

      if(i < iStart)
	os << COLOR(wipeColor) << "  " << DEFAULT;
      else if(i > iStop || s[i] == ' ')
	os << COLOR(borderColor) << "  " << DEFAULT;
      else
	os << COLOR(contentColor) << "  " << DEFAULT;
    }

    COUT(os.str());
    nline++;
    ++nLinePrinted_;
  }

  COUT(DEFAULT);

  return nline;
}

void Logo::getDimensions(std::string fileName, unsigned& nRow, unsigned& nCol)
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

void Logo::flash(std::string file, unsigned ntime, unsigned nBuffer)
{
  gcp::util::XtermManip macroXtm;
  struct timespec ts;
  ts.tv_sec  = 0;
  ts.tv_nsec = 50000000;

  unsigned nline;

  unsigned nRowLogo, nColLogo;
  getDimensions(file, nRowLogo, nColLogo);

  for(unsigned i=0; i < ntime; i++) {
    nline_ = loadLogo(file, 0, nColLogo, LOGO_COL2, LOGO_COL1, "default", nBuffer);

    macroXtm.moveCursorUp(nline_+1);

    nanosleep(&ts, 0);
    loadLogo(file, 0, nColLogo, LOGO_COL1, LOGO_COL2, "default", nBuffer);

    macroXtm.moveCursorUp(nline_+1);
    nanosleep(&ts, 0);
  }
}

void Logo::stripe(std::string file, unsigned nBuffer)
{
  gcp::util::XtermManip macroXtm;
  struct timespec ts;
  ts.tv_sec  = 0;
  ts.tv_nsec = 5000000;

  unsigned nRowLogo, nColLogo;
  getDimensions(file, nRowLogo, nColLogo);

  for(unsigned i=0; i <= nColLogo; i++) {
    nline_ = loadLogo(file, 0, i, LOGO_COL2, LOGO_COL1, "default", nBuffer);

    macroXtm.moveCursorUp(nline_+1);

    nanosleep(&ts, 0);
  }
}

void Logo::destripe(std::string file, unsigned nBuffer)
{
  gcp::util::XtermManip macroXtm;
  struct timespec ts;
  ts.tv_sec  = 0;
  ts.tv_nsec = 5000000;

  unsigned nRowLogo, nColLogo;
  getDimensions(file, nRowLogo, nColLogo);

  for(unsigned i=0; i <= nColLogo; i++) {
    nline_ = loadLogo(file, i, nColLogo, LOGO_COL2, LOGO_COL1, "default", nBuffer);

    if(i < nColLogo) {
      macroXtm.moveCursorUp(nline_+1);
    }
    
    nanosleep(&ts, 0);
  }
}

void Logo::display()
{
  gcp::util::XtermManip macroXtm;
  unsigned nCol = macroXtm.getNCol();
  unsigned nRow = macroXtm.getNRow();

  std::string file = getLogo();

  unsigned nRowLogo, nColLogo;
  getDimensions(file, nRowLogo, nColLogo);

  unsigned buffer = (nCol - 2*nColLogo)/2;

  macroXtm.hideCursor();

  if(nCol >= 2*nColLogo && nRow >= nRowLogo+2) {
    stripe(file, buffer);
    flash(file, 4, buffer);
    destripe(file, buffer);
  }

  macroXtm.showCursor();
}

void Logo::display2()
{
  macroXtm_.saveCursor();

  setDone(false);

  unsigned nCol = macroXtm_.getNCol();
  unsigned nRow = macroXtm_.getNRow();

  file_ = getLogo();
  std::string copy = file_;

  unsigned nRowLogo, nColLogo;
  getDimensions(file_, nRowLogo, nColLogo);

  buffer_ = (nCol - 2*nColLogo)/2;

  macroXtm_.hideCursor();

  struct timespec ts;

  unsigned ndiff=randomize(copy, file_, false);
  nline_ = loadLogo(copy, 0, nColLogo, LOGO_COL2, LOGO_COL1, "default", buffer_);

  ts.tv_sec = 0;

  ts.tv_nsec =  50000000;
  nanosleep(&ts, 0);
  ts.tv_nsec =  50000000;

  unsigned iter=0;
  while(ndiff > 0) {
    ++iter;
    macroXtm_.moveCursorUp(nline_+1);
    nanosleep(&ts, 0);
    loadLogo(copy, 0, nColLogo, LOGO_COL2, LOGO_COL1, "default", buffer_);
    ndiff=randomize(copy, file_, true, iter);
  }

  macroXtm_.moveCursorUp(nline_+1);

  flash(file_, 4, buffer_);
  destripe(file_, buffer_);

  macroXtm_.showCursor();

  setDone(true);
}

/**.......................................................................
 * Reset the window if we were in the middle of displaying a logo
 */
void Logo::reset()
{
  // Restore the cursor to the last saved position

  macroXtm_.restoreCursor();

  // Now blank up to the full size of the xterm

  unsigned nColWindow = macroXtm_.getNCol();
  unsigned nRowWindow = macroXtm_.getNRow();
  std::ostringstream os;
  
  for(unsigned i=0; i < nColWindow; i++)
    os << DEFAULT << " ";

  for(unsigned i=0; i < nRowWindow; i++)
    COUT(os.str());

  // Finally, reveal the cursor

  macroXtm_.showCursor();

  setDone(true);
}

bool Logo::coinFlip(unsigned iter)
{
  unsigned i = rand();
  if(iter < 10)
    return (i > RAND_MAX/2);
  else
    return (i > RAND_MAX/4);
}

char Logo::coinFlip(char c1, char c2)
{
  unsigned i = rand();
  if(i <= RAND_MAX/2)
    return c1;
  else
    return c2;
}

unsigned Logo::randomize(std::string& str, std::string& target, bool check, unsigned iter)
{
  unsigned ndiff=0;
  for(unsigned i=0; i < str.size(); i++) {
    char c = str[i];
    if(!check) {
      if(c != '\n') {
	str[i] = coinFlip('x', ' ');
	++ndiff;
      }
    } else {
      if(c != '\n') {
	if(str[i] != target[i]) {
	  ++ndiff;
	  if(coinFlip(iter)) {
	    str[i] = coinFlip('x', ' ');
	  }
	}
      }
    }
  }
  return ndiff;
}

std::string Logo::getLogo()
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

bool Logo::isDone()
{
  bool retVal = false;

  doneGuard_.lock();
  retVal = done_;
  doneGuard_.unlock();

  return retVal;
}

void Logo::setDone(bool done)
{
  doneGuard_.lock();
  done_ = done;
  doneGuard_.unlock();
}
