#include "gcp/util/Exception.h"
#include "gcp/util/XtermManip.h"

#include "gcp/util/CoProc.h"
#include "gcp/util/FdSet.h"
#include "gcp/util/String.h"

#include <iostream>
#include <unistd.h>

#include <stdlib.h>
#include <stdio.h>
#include <termios.h>
#include <sys/ioctl.h>

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
XtermManip::XtermManip() 
{
  initializeOptions();
}

void XtermManip::initializeOptions()
{
  //------------------------------------------------------------
  // Color foreground options
  //------------------------------------------------------------

  fg_[C_BLACK]   = "30";
  fg_[C_RED]     = "31";
  fg_[C_GREEN]   = "32";
  fg_[C_YELLOW]  = "33";
  fg_[C_BLUE]    = "34";
  fg_[C_MAGENTA] = "35";
  fg_[C_CYAN]    = "36";
  fg_[C_WHITE]   = "37";
  fg_[C_DEFAULT] = "39";
  fg_[C_ORANGE]  = "FE9A2E";
  fg_[C_OLIVE]   = "808000";
  fg_[C_TEAL]    = "008080";
  
  //------------------------------------------------------------
  // Color background options
  //------------------------------------------------------------

  bg_[C_BLACK]   = "40";
  bg_[C_RED]     = "41";
  bg_[C_GREEN]   = "42";
  bg_[C_YELLOW]  = "43";
  bg_[C_BLUE]    = "44";
  bg_[C_MAGENTA] = "45";
  bg_[C_CYAN]    = "46";
  bg_[C_WHITE]   = "47";
  bg_[C_DEFAULT] = "49";
  bg_[C_ORANGE]  = "FE9A2E";
  bg_[C_OLIVE]   = "808000";
  bg_[C_TEAL]    = "008080";

  //------------------------------------------------------------
  // Color key names
  //------------------------------------------------------------

  colorNames_["black"]   = C_BLACK;
  colorNames_["red"]     = C_RED;
  colorNames_["green"]   = C_GREEN;
  colorNames_["yellow"]  = C_YELLOW;
  colorNames_["blue"]    = C_BLUE;
  colorNames_["magenta"] = C_MAGENTA;
  colorNames_["cyan"]    = C_CYAN;
  colorNames_["white"]   = C_WHITE;
  colorNames_["default"] = C_DEFAULT;
  colorNames_["orange"]  = C_ORANGE;
  colorNames_["olive"]   = C_OLIVE;
  colorNames_["teal"]    = C_TEAL;

  //------------------------------------------------------------
  // Text modes
  //------------------------------------------------------------

  textModes_[TEXT_DEFAULT]        =  "0"; 
  textModes_[TEXT_BOLD]           =  "1"; 
  textModes_[TEXT_NORMAL]         = "22";
  textModes_[TEXT_UNDERLINED]     =  "4"; 
  textModes_[TEXT_NOT_UNDERLINED] = "24";
  textModes_[TEXT_BLINK]          =  "5"; 
  textModes_[TEXT_STEADY]         = "25";
  textModes_[TEXT_INVERSE]        =  "7"; 
  textModes_[TEXT_POSITIVE]       = "27";
  textModes_[TEXT_INVISIBLE]      =  "8"; 
  textModes_[TEXT_VISIBLE]        = "28";

  //------------------------------------------------------------
  // Text mode names
  //------------------------------------------------------------

  textModeNames_["default"]        = TEXT_DEFAULT;
  textModeNames_["bold"]           = TEXT_BOLD;          
  textModeNames_["normal"]         = TEXT_NORMAL;        
  textModeNames_["underlined"]     = TEXT_UNDERLINED;    
  textModeNames_["not underlined"] = TEXT_NOT_UNDERLINED;
  textModeNames_["blink"]          = TEXT_BLINK;         
  textModeNames_["steady"]         = TEXT_STEADY;        
  textModeNames_["inverse"]        = TEXT_INVERSE;       
  textModeNames_["positive"]       = TEXT_POSITIVE;      
  textModeNames_["invisible"]      = TEXT_INVISIBLE;     
  textModeNames_["visible"]        = TEXT_VISIBLE;       
}

/**.......................................................................
 * Destructor.
 */
XtermManip::~XtermManip() {}

std::string XtermManip::restoreDefaultColors()
{
  return "\033[0m";
}

void XtermManip::setFg(std::string name)
{
  setFg(getColorKey(name));
}

void XtermManip::setBg(std::string name)
{
  setBg(getColorKey(name));
}

void XtermManip::setHexBg(std::string hexCode)
{
  std::cout << "\033]11;#" << hexCode << "\007";
}

void XtermManip::setHexFg(std::string hexCode)
{
  std::cout << "\033]10;#" << hexCode << "\007";
}

void XtermManip::setCursorColor(std::string name)
{
  setCursorColor(getColorKey(name));
}

void XtermManip::setCursorColor(ColorKey key)
{
  std::map<ColorKey, std::string>::iterator slot;

  slot = fg_.find(key);

  if(slot != fg_.end()) {
    COUT("Issuign command with se = " << slot->second);
    std::cout << "\033]12;" << slot->second << "\007";
  }
}

void XtermManip::setFg(ColorKey key)
{
  std::map<ColorKey, std::string>::iterator slot;

  slot = fg_.find(key);

  if(slot != fg_.end()) {
    if(slot->second.size() > 2) {
      std::cout << "\033]10;#" << slot->second << "\007";
    } else {
      std::cout << "\e[" << slot->second << "m";
    }
  }
}

std::string XtermManip::fg(std::string name)
{
  return fg(getColorKey(name));
}

std::string XtermManip::fg(ColorKey key)
{
  std::map<ColorKey, std::string>::iterator slot;
  std::ostringstream os;

  slot = fg_.find(key);

  if(slot != fg_.end()) {
    if(slot->second.size() > 2) {
      os << "\033]10;#" << slot->second << "\007";
    } else {
      os << "\e[" << slot->second << "m";
    }
  }

  return os.str();
}

void XtermManip::setBg(ColorKey key)
{
  std::map<ColorKey, std::string>::iterator slot;

  slot = bg_.find(key);

  if(slot != bg_.end()) {

    if(slot->second.size() > 2) {
      std::cout << "\033]11;#" << slot->second << "\007";
    } else {
      std::cout << "\e[" << slot->second << "m";
    }

  }
}

std::string XtermManip::bg(std::string name)
{
  return bg(getColorKey(name));
}

std::string XtermManip::bg(ColorKey key)
{
  std::map<ColorKey, std::string>::iterator slot;
  std::ostringstream os;

  slot = bg_.find(key);

  if(slot != bg_.end()) {
    if(slot->second.size() > 2) {
      os << "\033]11;#" << slot->second << "\007";
    } else {
      os << "\e[" << slot->second << "m";
    }
  }

  return os.str();
}

XtermManip::ColorKey XtermManip::getColorKey(std::string name)
{
  std::map<std::string, ColorKey>::iterator slot;

  slot = colorNames_.find(name);

  if(slot != colorNames_.end()) {
    return slot->second;
  }

  ThrowError("No such key: " << name);
}

void XtermManip::setTextMode(std::string name)
{
  setTextMode(getTextModeKey(name));
}

void XtermManip::setTextMode(TextModeKey key)
{
  std::map<TextModeKey, std::string>::iterator slot;

  slot = textModes_.find(key);

  if(slot != textModes_.end()) {
    std::cout << "\e[" << slot->second << "m";
  }
}

std::string XtermManip::textMode(std::string name)
{
  return textMode(getTextModeKey(name));
}

std::string XtermManip::textMode(TextModeKey key)
{
  std::map<TextModeKey, std::string>::iterator slot;
  std::ostringstream os;

  slot = textModes_.find(key);

  if(slot != textModes_.end()) {
    os << "\e[" << slot->second << "m";
  }

  return os.str();
}

XtermManip::TextModeKey XtermManip::getTextModeKey(std::string name)
{
  std::map<std::string, TextModeKey>::iterator slot;

  slot = textModeNames_.find(name);

  if(slot != textModeNames_.end()) {
    return slot->second;
  }

  ThrowError("No such key: " << name);
}

void XtermManip::saveCursor()
{
  cout << "\e7";
}

void XtermManip::restoreCursor()
{
  cout << "\e8";
}

void XtermManip::getCursorPosition()
{
  setRawMode();

  FdSet fdSet;
  fdSet.registerReadFd(STDIN_FILENO);

  cout << "\e[6n";

  int nready = select(fdSet.size(), fdSet.readFdSet(), 0, 0, 0);

  COUT("Found: " << nready);
}

void XtermManip::clearAbove()
{
  cout << "\e[1J";
}

void XtermManip::moveCursorUp(unsigned nline)
{
  std::ostringstream os;

  os << nline;

  cout << "\e[" << os.str() << "A";
}

void XtermManip::hideCursor()
{
  cout << "\e[?25l";
}

void XtermManip::showCursor()
{
  cout << "\e[?25h";
}

void XtermManip::moveCursorDown(unsigned nline)
{
  std::ostringstream os;

  os << nline;

  cout << "\e[" << os.str() << "B";
}

void XtermManip::setRawMode()
{
  int fd = STDIN_FILENO;
  struct termios t;

  if (tcgetattr(fd, &t) < 0) 
    ThrowSysError("tcgetattr");
  
  t.c_lflag &= ~ICANON;
  
  if (tcsetattr(fd, TCSANOW, &t) < 0)
    ThrowSysError("tcsetattr");

  setbuf(stdin, NULL);
}

unsigned XtermManip::getNRow()
{
  struct winsize ws;
  ioctl(0, TIOCGWINSZ, &ws);
  return ws.ws_row;
}

unsigned XtermManip::getNCol()
{
  struct winsize ws;
  ioctl(0, TIOCGWINSZ, &ws);
  return ws.ws_col;
}

std::string XtermManip::exec(std::string cmd) 
{
  FILE* pipe = popen(cmd.c_str(), "r");
  if (!pipe) return "ERROR";
  char buffer[128];
  std::string result = "";
  while(!feof(pipe)) {
    if(fgets(buffer, 128, pipe) != NULL)
      result += buffer;
  }
  pclose(pipe);
  return result;
}

int getNbyte(int fd) 
{
  int request = FIONREAD;
  int nByte=0;
  
  if(ioctl(fd, request, &nByte) != 0) {
    ThrowSysError("ioctl()");
  }
  
  return nByte;
}

std::string XtermManip::getBg()
{
  return getColor(11);
}

std::string XtermManip::getFg()
{
  return getColor(10);
}

std::string XtermManip::getColor(unsigned int color)
{
  ostringstream os;
  os << SCRIPTDIR << "/getColor.csh " << color;
  CoProc proc(os.str());

  FILE* stdIn  = proc.stdIn()->writeFp();
  FILE* stdOut = proc.stdOut()->readFp();

  fclose(stdIn);

  char c;
  os.str("");

  while((c = (char)fgetc(stdOut)) != EOF)
    os << c;

  String str(os.str());

  String r = str.findNextInstanceOf("rgb:", true, "/", true, false);
  String g = str.findNextInstanceOf("/", true, "/", true, true);
  String b = str.remainder();

  ostringstream colorOs;
  colorOs << r[0] << r[1] << g[0] << g[1] << b[0] << b[1];

  return colorOs.str();
}

