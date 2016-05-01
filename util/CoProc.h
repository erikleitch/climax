// $Id: CoProc.h,v 1.1.1.1 2009/07/06 23:57:25 eml Exp $

#ifndef GCP_UTIL_COPROC_H
#define GCP_UTIL_COPROC_H

/**
 * @file CoProc.h
 * 
 * Tagged: Fri Feb  9 09:32:45 NZDT 2007
 * 
 * @version: $Revision: 1.1.1.1 $, $Date: 2009/07/06 23:57:25 $
 * 
 * @author Erik Leitch
 */

#include <iostream>
#include <vector>

namespace gcp {
  namespace util {

    class CoProc;

    // Define a class for managing a spawned process whose stdout and
    // stderr are captured by the parent

    class CoProc {
    public:

      class Pipe {

      public:

	~Pipe();

	FILE* readFp();
	FILE* writeFp();

	int readFd();
	int writeFd();

	void closeForWrite();
	void closeForRead();

	void close();

      private:

	friend class CoProc;

	int fd_[2];     // pair of file descriptors

	Pipe();
      };

      /**
       * Constructor.
       */
      CoProc(std::string exe);

      CoProc(char *exe);

      CoProc(char *exe, FILE** stdInFp, FILE** stdOutFp, FILE** stdErrFp);

      CoProc(char *exe, char** argv, 
	     FILE** stdInFp, FILE** stdOutFp, FILE** stdErrFp);

      void fork(char *exe, char** argv, 
		FILE** stdInFp, FILE** stdOutFp, FILE** stdErrFp);

      CoProc::Pipe* getPipe();

      CoProc::Pipe* stdIn();
      CoProc::Pipe* stdOut();
      CoProc::Pipe* stdErr();

      /**
       * Destructor.
       */
      virtual ~CoProc();
      
      // Split a command string into separate tokens

      static std::vector<std::string> split(std::string command);

      int readFd();
      int writeFd();

    private:

      int pid_;
      CoProc::Pipe* stdIn_;
      CoProc::Pipe* stdOut_;
      CoProc::Pipe* stdErr_;

    }; // End class CoProc
    
  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_COPROC_H
