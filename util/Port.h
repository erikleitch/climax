#ifndef GCP_UTIL_PORT_H
#define GCP_UTIL_PORT_H

/**
 * @file Port.h
 * 
 * Tagged: Mon May 10 16:41:20 PDT 2004
 * 
 * @author Erik Leitch
 */
#include <string>
#include <deque>
#include "gcp/util/String.h"
#include "gcp/util/Vector.h"

namespace gcp {
  namespace util {
    
    class Port {
    public:
      
      static const unsigned int MAX_RCV_BUFFER = 300;
      
      // The default telnet port number
      
      static const unsigned TELNET_PORT_NO = 23;
      
      /**
       * Constructor
       */
      Port(int fd = -1);
      
      /**
       * Destructor
       */
      virtual ~Port();
      
      /**
       * Return the file descriptor associated with our port connection
       */
      inline int getFd() {
	return fd_;
      }
      
      /**
       * Connect to the port
       */
      virtual int connect() {return 0;};
      
      /**
       * Write a message to the port
       */
      void writeString(std::string& message, int fd=-1);
      void writeBytes(Vector<unsigned char>& buffer);
      static void writeBytes(Vector<unsigned char>& buffer, int fd);
      
      /**  
       * Read a message from the port
       */
      unsigned int readBytes(unsigned char *message, int fd=-1);
      unsigned int readBytes(Vector<unsigned char>& buffer);
      
      // Read tne next byte from the serial port 
      
      unsigned char getNextByte();

      std::string readString(int fd=-1);
      
      bool concatenateString(std::ostringstream& os, int fd=-1, bool cont=true);
      void concatenateChar(std::ostringstream& os, int fd=-1);
      int getNbyte(int fd=-1);
      
      /**
       * terminate a read when any of the following characters are
       * read from the port
       */
      void terminateAt(std::string strip);
      
      /**
       * Strip any of the following characters from data read from the
       * port
       */
      void strip(std::string strip);
      
      /**
       * Characters we mustn't strip.  Note that this will override
       * duplicate characters implied by stripUnprintable()
       */
      void dontStrip(std::string strip);
      
      /**
       * Strip any unprintable characters
       */
      void stripUnprintable(bool strip);
      
      /**
       * Append the passed string to the end of each line
       */
      void append(std::string append);
      
      /**
       * Set no buffering for this stream
       */
      void setNoBuf();
      
      /**
       * Set line buffering for this stream
       */
      void setLineBuf();
      
    protected:
      
      int fd_;
      
    private:
      
      String termStr_;
      String stripStr_;
      String dontStripStr_;
      String appendStr_;
      bool stripUnprintable_;
      
    }; // End class Port
    
  } // End namespace util
} // End namespace gcp




#endif // End #ifndef GCP_UTIL_PORT_H
