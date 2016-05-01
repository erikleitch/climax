#ifndef GCP_PROGRAM_PROGRAM_H
#define GCP_PROGRAM_PROGRAM_H

/**
 * @file Program.h
 * 
 * Tagged: Mon Jun 14 13:22:06 CDT 2004
 * 
 * @author Erik Leitch
 */
#include <string>

#define END_OF_KEYWORDS  "end_of_keywords"

namespace gcp {
  namespace program {

    /**
     * @typedef KeyTabEntry Convenience structure to hold the
     * key=value/type/help-string information.
     */
    struct KeyTabEntry {
      const char *key;   // keyword (name used on commad line)
      const char *val;   // default value for keyword
      const char *type;  // type of value - could be 'f', 'd', 'b', 's'
      const char *help;  // help string - used to provide detailed 
    };
    
    /**
     * @typedef Keyword
     * private structure used by Program to hold all relevant information
     * per keyword.
     * @see Program
     */
    struct Keyword {
      std::string  keyValue;  // pointer to original R/O [not used]
      char    option;    // short option [not used]
      std::string  key;       // keyword (same as long option)
      std::string  val;       // value
      std::string  help;      // help
      int     count;     // update/reference count; 0=original default
      int     upd;       // 0=read 1=unread original 2=unread updated
      int     system;    // system keyword ? (0=no 1=yes)
    };

    class Program {
    public:
      
      /**
       * Static array of user-defined keywords
       */
      static KeyTabEntry keywords[];  // key,val,type,help string
				      // tuples for program keys

      /**
       * Destructor.
       */
      virtual ~Program();
      
      /**
       * debugging routine 
       */
      void show(void);     
      
      /**
       * @brief get string value of a parameter
       * @return string value of a parameter
       *
       */
      std::string getStringParameter(std::string key); 

      /**
       * @brief get double value of a parameter
       * @return double value of a parameter
       */
      double getDoubleParameter(std::string key);

      /**
       * @brief get int value of a parameter
       * @return int value of a parameter
       */
      int getIntegerParameter(std::string key);

      /**
       * @brief get boolean value of a parameter
       * @return boolean value of a parameter
       */
      bool getBooleanParameter(std::string key);

      /**
       * @brief does keyword have a value
       * @return does the keyword have a non-blank (assigned) value
       */
      bool hasValue(std::string key);

      /**
       * @brief 
       * @return has the keyword been used above their default definition
       * After keywords are instantiated from a KeyTabEntry (as defined
       * by the author of the program),
       */
      bool isDefault(std::string key);
      
      /**
       * @brief usage count of a keyword
       * @return number of times a keyword has been accessed
       */
      int count(std::string key);

      /**
       * Run a program
       */
      int run(int argc, char* argv[]);

      /**
       * @brief Class global method to get the process-wide singleton
       * instance of Program.
       *
       * Class global method that returns a reference to the process-wide
       * single instance of class Program. Basically a factory method that
       * is tailored to return the same instance for every call.
       *
       * @return Program& - reference to the singleton instance of class
       * 		      carma::util::Program.
       *
       * @see Singleton pattern in "Design Patterns", by 
       *      Gamma, Helm, Johnson, and Vlissides.
       */
      static Program &getProgram(); // get singleton instance of Program

    private:

      static KeyTabEntry system_[];   // key,val,type,help string
				      // tuples for system keys
      static std::string version_;      // N.M date author (and/or CVS id's)
      static std::string usage_;        // one line usage
      static std::string description_;  // a longer multi-line description

      std::string progName_;               // short program name without path

      int argc_;                      // (from main) argument count
      char **argv_;                   // (from main) array of CL strings
      int ddLoc_;                     // location of argv "--" (0 means
				      // absent)
      
      int nkeys_;		      // total number of keywords + 1
				      // (argv0, prog, sys)
      int npkeys_;                    // number of program keywords
      int nskeys_;                    // number of system keywords
      Keyword* keys_;                  // 0=not used 1=first keyword etc.

      // Private methods

      /**
       * @user main() program, returns the program exit status
       * @return exit status requested (normally 0)
       * Program::main is the routine written by the user. 
       */
      int main(void);
      void initializeUsage();

      /**
       * Making the constructor a private member function makes Program 
       * a singleton, The only way anybody can get at Program is via 
       * Progam::getProgram(), which returns a reference to the 
       * singleton instance.
       */
      Program();
      
      /**
       * initialize is called by run() to setup the Program environment:
       * initialize/setup, initial check if system-only information requested,
       * parse the command line, ...
       */
      int initialize();

      void printGreeting();
      std::string fillLine(unsigned char c, unsigned length);
      std::string formatLine(std::string line,
                             unsigned leftPad, unsigned rightPad, unsigned length);

      /**
       * called by run, this function makes sure Carma programs are
       * properly released from the carma environment
       */
      void terminate();

      /**
       * @return 
       */
      std::string getProgname(void);

      /**
       * @return 
       */
      int getArgc(void);

      /**
       * @return 
       */
      char **getArgv(void);

      /**
       * Add a key to the array of known keys
       */
      void addKey(int i, KeyTabEntry &kt, int system=0);

      /**
       * Return the index of the requested key in the keyword array
       */
      int findKey(const std::string);

    }; // End class Program
    
  } // End namespace program
} // End namespace gcp



#endif // End #ifndef GCP_PROGRAM_PROGRAM_H
