#ifndef GCP_UTIL_XTERMMANIP_H
#define GCP_UTIL_XTERMMANIP_H

/**
 * @file XtermManip.h
 * 
 * Tagged: Sun Jan 20 13:34:20 NZDT 2008
 * 
 * @version: $Revision: 1.1 $, $Date: 2012/05/02 23:44:52 $
 * 
 * @author Erik Leitch
 */

#include <iostream>
#include <map>
#include <string>

namespace gcp {
  namespace util {

    class XtermManip {
    public:

      /**
       * Constructor.
       */
      XtermManip();

      /**
       * Destructor.
       */
      virtual ~XtermManip();

      struct Option {
	unsigned key_;
	std::string str_;

	Option(unsigned key, std::string str) {
	  key_ = key;
	  str_ = str;
	}
      };

      //------------------------------------------------------------
      // Color manipulation
      //------------------------------------------------------------

      enum ColorKey {
	C_BLACK,
	C_RED,
	C_GREEN,
	C_YELLOW,
	C_BLUE,
	C_MAGENTA,
	C_CYAN,
	C_WHITE,
	C_DEFAULT,
	C_ORANGE,
	C_OLIVE,
	C_TEAL
      };

      std::map<ColorKey, std::string> fg_;
      std::map<ColorKey, std::string> bg_;

      std::map<std::string, ColorKey> colorNames_;

      std::string restoreDefaultColors();

      void setFg(ColorKey key);
      void setFg(std::string key);

      void setBg(ColorKey key);
      void setBg(std::string key);
      
      void setHexFg(std::string hexCode);
      void setHexBg(std::string hexCode);

      void setCursorColor(std::string name);
      void setCursorColor(ColorKey key);

      std::string fg(ColorKey key);
      std::string fg(std::string key);

      std::string bg(ColorKey key);
      std::string bg(std::string key);

      std::string getColor(unsigned int color);
      std::string getFg();
      std::string getBg();

      ColorKey getColorKey(std::string name);

      std::string exec(std::string cmd);

      //------------------------------------------------------------
      // Text mode manipulation
      //------------------------------------------------------------

      enum TextModeKey {
	TEXT_DEFAULT,

	TEXT_BOLD,
	TEXT_NORMAL,

	TEXT_UNDERLINED,
	TEXT_NOT_UNDERLINED,

	TEXT_BLINK,
	TEXT_STEADY,

	TEXT_INVERSE,
	TEXT_POSITIVE,

	TEXT_INVISIBLE,
	TEXT_VISIBLE,
      };

      std::map<TextModeKey, std::string> textModes_;
      std::map<std::string, TextModeKey> textModeNames_;

      void setTextMode(TextModeKey key);
      void setTextMode(std::string name);

      std::string textMode(TextModeKey key);
      std::string textMode(std::string name);

      TextModeKey getTextModeKey(std::string name);

      //------------------------------------------------------------
      // Cursor modes
      //------------------------------------------------------------

      enum CursorCntlKey {
	CURS_GET_POSITION,
      };

      void getCursorPosition();

      void saveCursor();
      void restoreCursor();
      void moveCursorUp(unsigned nline);
      void moveCursorDown(unsigned nline);
      void showCursor();
      void hideCursor();

      void clearAbove();

      static unsigned getNRow();
      static unsigned getNCol();

    private:

      void initializeOptions();
      void setRawMode();

    }; // End class XtermManip

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_XTERMMANIP_H
