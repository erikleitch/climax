// $Id: $

#ifndef GCP_UTIL_LOGO_H
#define GCP_UTIL_LOGO_H

/**
 * @file Logo.h
 * 
 * Tagged: Wed Sep 18 14:37:41 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include <string>

namespace gcp {
  namespace util {

    class Logo {
    public:

      /**
       * Constructor.
       */
      Logo();

      /**
       * Destructor.
       */
      virtual ~Logo();

      std::string getLogo();
      unsigned loadLogo(std::string str, unsigned iStart, unsigned iStop, 
			std::string borderColor, std::string contentColor, std::string wipeColor, unsigned nBuffer);
      void getDimensions(std::string fileName, unsigned& nRow, unsigned& nCol);
      void flash(std::string file, unsigned ntime, unsigned nBuffer);
      void stripe(std::string file, unsigned nBuffer);
      void destripe(std::string file, unsigned nBuffer);

      unsigned randomize(std::string& str, std::string& target, bool check, unsigned iter=0);
      char coinFlip(char c1, char c2);
      bool coinFlip(unsigned iter=0);

      void display();
      void display2();

      Mutex doneGuard_;

      void setDone(bool done);
      bool isDone();
      void reset();

      bool done_;
      unsigned buffer_;
      std::string file_;
      XtermManip macroXtm_;
      unsigned nline_;
      unsigned nLinePrinted_;

    }; // End class Logo

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_LOGO_H
