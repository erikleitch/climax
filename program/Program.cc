#include <iostream>
#include <cassert>

#include "gcp/program/Program.h"
#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/Logger.h"

#include <string.h>

#define PROGRAM_VERSION "Program:: $Date: 2010/07/13 17:56:49 $ $Revision: 1.1.1.1 $"
using namespace std;
using namespace gcp::program;
using namespace gcp::util;

// local functions, currently defined near the end of this file.
//      they are all old NEMO routines and exist because of a lack
//      of a real string class

static string parname(string arg);
static string parvalue(string arg);
static string parhelp(string arg);
static bool   getbvalue(string sp);

std::string Program::version_("none available");
std::string Program::usage_("none available");
std::string Program::description_("none available");

/**
 * the system keywords, known to all programs
 */
KeyTabEntry Program::system_[] = {
  { "logfile",    "(none)", "s", "Logfile prefix"},
  { "logdir",     ".",      "s", "Logfile directory"},
  { "greeting",   "t",      "b", "Print greeting on startup"},
  {END_OF_KEYWORDS},
};

/**.......................................................................
 * Constructor.
 */
Program::Program() {}

/**.......................................................................
 * Destructor.
 */
Program::~Program() {}

/**.......................................................................
 * Initialize the program
 */
int Program::initialize(void)
{
  bool posflag = false;  // force key=val, don't allow positional
			 // matching like NEMO
  string pn(argv_[0]), name;
  int i, islash;
  
  initializeUsage();

  // Store the program name
  
  islash = pn.find_last_of('/');
  if (islash >= 0)
    progName_ = pn.substr(islash+1);
  else
    progName_ = pn;
  
  // Count program keys
  
  for (npkeys_=0; ; npkeys_++)   
    if (strcmp(keywords[npkeys_].key,END_OF_KEYWORDS) == 0) break;
  
  // Count system keys

  for (nskeys_=0; ; nskeys_++)    
    if (strcmp(system_[nskeys_].key,END_OF_KEYWORDS) == 0) break;
  
  nkeys_ = npkeys_ + nskeys_ + 1;      // one extra for the "argv0"
  
  // Store information about each key in a Keyword struct
  
  keys_ = new Keyword[nkeys_];       
  keys_[0].key = "argv0";
  keys_[0].val = progName_;
  keys_[0].help = "the name of the Program";
  
  for (i=0; i<nkeys_; i++)
    keys_[i].count = -1;     // signal that keyword was not initialized
  
  for (i=0; i<npkeys_; i++)
    addKey(i+1,keywords[i],0);
  
  for (i=0; i<nskeys_; i++)
    addKey(npkeys_+i+1, system_[i], 1);
  
  // Parse allowed command line arguments
  
  if (argc_ > 1) {           // handle common system shortcuts
    int n=strlen(argv_[1]);
    
    if (strcmp (argv_[1], "--") == 0) {
      
      // rest of the commandline is supposed to be skipped
      
      if (argc_==2) return 1;
      
    } else if (strncmp (argv_[1], "--help",n) == 0) {
      cerr << PROGRAM_VERSION << endl;
      cerr << "Current system options:" << endl;
      cerr << " --help         this help" << endl;
      cerr << " --keywords     show program keys, values and help" << endl;
      cerr << " --system       show system keys, values and help" << endl;
      cerr << " --version      show program version" << endl;
      cerr << " --usage        show program usage line" << endl;
      cerr << " --description  show program description section" << endl;
      cerr << " --show         show some debugging info" << endl;
      
      if (argc_==2) return 1;
      
    } else if (strncmp (argv_[1], "--keywords",n) == 0) {
      
      for (int i=0; i < npkeys_; i++) {
	cout << keys_[1+i].key << "=" << keys_[1+i].val << endl;
	cout << "\t" << keys_[1+i].help << endl;
      }
      
      if (argc_==2) return 1;
      
    } else if (strncmp (argv_[1], "--system",n) == 0) {
      
      for (int i=0; i < nskeys_; i++) {
	cout << keys_[npkeys_+1+i].key << "=" << keys_[npkeys_+1+i].val << endl;
	cout << "\t" << keys_[1+npkeys_+i].help << endl;
      }
      
      if (argc_==2) return 1;
      
    } else if (strncmp (argv_[1], "--version",n) == 0) {
      cout << progName_ << " : " << version_ << endl;
      if (argc_==2) return 1;
      return 0;
    } else if (strncmp (argv_[1], "--usage",n) == 0) {
      cout <<  usage_ << endl;
      if (argc_==2) return 1;
    } else if (strncmp (argv_[1], "--description",n) == 0) {
      cout << description_ << endl;
      if (argc_==2) return 1;
    } else if (strncmp (argv_[1], "--show",n) == 0) {
      show();
      if (argc_==2) return 1;
    }
  }
  
  // check if any  '???' values remain
  // ??? is NEMO's way to enforce a value, CARMA doesn't need that ???
  
  int missing = 0;
  for (int i=1; i<nkeys_; i++)
    if (keys_[i].val == "???") {    
      missing++;
      break;
    }
  if (missing) {
    cerr << "Insufficient parameters" << endl;
    cerr << "Usage: " << progName_;
    bool otherargs = false;
    for (int i=1; i<nkeys_; i++)
      if (keys_[i].val == "???")
	cerr << keys_[i].key <<  "=???" ;
      else
	otherargs = true;
    cerr << (otherargs ? " ...\n" : "\n");
    cerr << "The above required keywords have no value" << endl;
    exit(1);
  }
  
  ddLoc_ = 0;  // location of "--", if present (0=absent)
  
  for (int i=1; i<argc_ ; i++) {        // parse command line
    if (strcmp(argv_[i],"--")==0) {
      // cerr << "Found -- : skipping rest of cmdline: " << i << endl;
      ddLoc_ = i;
      break;
    }
    name = parname(string(argv_[i]));
    posflag = posflag && (name.size() == 0);
    if (posflag) {
      if (i >= nkeys_) {
	cerr << "Too many un-named arguments" << endl;
	exit(1);
      }
      // should free old one here too
      keys_[i].val = argv_[i];
      keys_[i].count++;
    } else {
      if (name.size() == 0) {
	cerr << "Parameter " << name << "must be named" << endl;
	exit(1);
      }
      int j = findKey(name);
      if  (j >= 0) {
	if(keys_[j].count) {
	  cerr << "Parameter duplicated" << endl;
	  exit(1);
	}
	// free old value here too
	keys_[j].val = parvalue(argv_[i]);
	keys_[j].count++;
      } else {
	cerr << "Parameter \"" << name << "\" unknown" << endl;
	exit(1);
      } // j>0
    } // if(posflag)
  } // for (i)

  // Print a greeting message                                                                                                                                                          
  if(getBooleanParameter("greeting")) {
    printGreeting();
  }

  // See if a logfile was requested

  try {
    if(!isDefault("logfile")) {
      gcp::util::Logger::setLogFilePrefix(getStringParameter("logfile"));
      gcp::util::Logger::setLogFileDirectory(getStringParameter("logdir"));
      gcp::util::Logger::openLogFile();
    }
  } catch(...) {
    COUT("Here about to exit");
    exit(1);
  }

  return 0;
}

void Program::printGreeting()
{
  ostringstream welcome, type, usage, desc, listing, help;

  welcome << "Welcome to " << progName_ << ". Hit cntl-C to exit.";
  type    << "If you do not know how to use this program, type:";
  usage   << "'" << progName_ << " --usage'       for a synopsis";
  desc    << "'" << progName_ << " --description' for a detailed description";
  listing << "'" << progName_ << " --keywords'    for a listing of available keywords";
  help    << "'" << progName_ << " --help'        for this and other info";

  // Get the length of the longest line, plus padding                                                                                                                                  

  unsigned leftPad   = 1;
  unsigned rightPad  = 1;
  unsigned indentPad = 5;
  unsigned length    = listing.str().size() + leftPad + indentPad + rightPad;

  std::cout << std::endl
            << fillLine('-', length)                                << std::endl
            << formatLine(welcome.str(), leftPad, rightPad, length) << std::endl << std::endl
            << formatLine(type.str(),    leftPad, rightPad, length) << std::endl << std::endl
            << formatLine(usage.str(),   leftPad+indentPad, rightPad, length) << std::endl
            << formatLine(desc.str(),    leftPad+indentPad, rightPad, length) << std::endl
            << formatLine(listing.str(), leftPad+indentPad, rightPad, length) << std::endl
            << formatLine(help.str(),    leftPad+indentPad, rightPad, length) << std::endl
            << fillLine('-', length)                                << std::endl << std::endl;
}

std::string Program::fillLine(unsigned char c, unsigned length)
{
  ostringstream os;
  for(unsigned i=0; i < length; i++)
    os << c;
  return os.str();
}

std::string Program::formatLine(std::string line,
                                unsigned leftPad, unsigned rightPad, unsigned length)
{
  ostringstream os;
  for(unsigned i=0; i < leftPad; i++)
    os << " ";

  os << line;

  for(unsigned i=0; i < length-(leftPad+line.size()+rightPad); i++)
    os << " ";

  return os.str();
}

/**.......................................................................
 * Terminate the program
 */
void Program::terminate() {}

/**.......................................................................
 * Run a program
 */
int Program::run(int argc, char* argv[])
{
  LogStream errStr;
  int  returnValue = 0;  // success = 0, failure = -1
  
  argc_ = argc;
  argv_ = argv;
  
  try  {
    returnValue = initialize(); // initialize all
    if (returnValue) 
      return 0;  // nothing to do after all ... return right now
  } catch (...) {
    returnValue = -1;
  }
  
#if 0
  std::string defaultFg;
  std::string defaultBg;
  XtermManip xtm;

  defaultFg = xtm.getFg();
  defaultBg = xtm.getBg();

  xtm.setHexBg("000000");
#endif

  try  {
    returnValue = this->main();  // call Program::main() -- the real thing --
  } catch (gcp::util::Exception& err)  {

    fflush(stdout);

#if 0
    errStr.initMessage(true);
    errStr << err.what();
    errStr.report();
#endif

    COUT(err.what());

    returnValue = -1;

#if 0
    xtm.setHexBg(defaultBg);
    xtm.setHexFg(defaultFg);
#endif

  } catch(...) {
    errStr.appendMessage(true, "Caught an unknown exception");
    errStr.report();
    returnValue = -1;
  }

#if 0
  xtm.setHexBg(defaultBg);
  xtm.setHexFg(defaultFg);
#endif
  
  try  {
    Program::terminate();                  // all done
  } catch (...)  {
    returnValue = -1;
  }
  
  return returnValue;
}

/**.......................................................................
 * Show a listing of this program
 */
void Program::show(void) 
{ 
  string *sp;
  int i;
  
  cout << "SHOW-begin"     << endl;
  cout << " Program: "     << progName_ << endl;
  cout << " Usage: "       << usage_    << endl;
  cout << " Version: "     << version_  << endl;
  cout << " Description: " << endl << description_ << endl;
  
  for (i=1; i < argc_; i++)
    cout << " arg(" << i << ") = " << argv_[i] << endl;
  
  for (i=1; i < nkeys_; i++) {
    cout << " key(" << i << ") = " << keys_[i].key;
    if (keys_[i].system) cout << " (**system key**)";
    cout << endl;
    cout << " def(" << i << ") = " << keys_[i].val << endl;
    cout << " help(" << i << ") = " << keys_[i].help << endl;
  }
  
  cout << "SHOW-end" << endl;
}

/**.......................................................................
 * main() is defined here, application programmers will need to define
 * Program::main() as well as a set of (actually optional) keywords,
 * usage, version, description that the "keys" program parses into C++
 * code to be linked in during program compilation
 */
int main(int argc, char* argv[]) 
{ 
  int retval;
  Program &program = Program::getProgram();
  
  retval = program.run(argc, argv);
  
  // delete program.Keyword;
  
  delete &program;
  
  return retval;
}

/**.......................................................................
 * Add a key to the list of known keys
 */
void Program::addKey(int i, KeyTabEntry &kt, int system)
{
  int j;
  
  if (i >= nkeys_) {
    cerr << "addKey internal error i=" << i << endl;
    exit(1);
  }
  
  keys_[i].keyValue = "";   // currently not used
  keys_[i].option = 0;      // currently not used
  keys_[i].key = kt.key;
  keys_[i].val = kt.val;
  keys_[i].help = kt.help;
  keys_[i].count = 0;       // before addKey, count was -1
  keys_[i].upd = 1;
  keys_[i].system = system;
  
  // test internal consistencies, duplicate keys etc.
  
  for (j=0; j < nkeys_; j++) {
    if (i!=j && keys_[j].count >= 0) {
      if (keys_[j].key == keys_[i].key)
	cerr << "Keyword " << kt.key << " is duplicated" << endl;
    }
  }
}

/**.......................................................................
 * Return true if the named key has a value.
 */
bool Program::hasValue(string key) 
{
  int i = findKey(key);
  
  if (i<0) {
    cerr << key << " Unknown keyword" << endl;
    exit(1);
  }
  
  return keys_[i].val.size() > 0;
}

/**.......................................................................
 * Return true if the named key has its default value
 */
bool Program::isDefault(string key) 
{
  int i = findKey(key);
  
  if (i<0) {
    cerr << key << " Unknown keyword" << endl;
    exit(1);
  }
  
  return count(key)==0;
}

/**.......................................................................
 * Return the string version of a keyword
 */
string Program::getStringParameter (string key) 
{
  int i = findKey(key);
  
  if (i<0) {
    cerr << key << " Unknown keyword" << endl;
    exit(1);
  }
  
  keys_[i].upd = 0;        // mark keyword as read
  return keys_[i].val;
}

/**.......................................................................
 * Return the value of an integer keyword
 */
int Program::getIntegerParameter (string key) 
{
  string sp = getStringParameter(key);
  return atoi(sp.c_str());
}

/**.......................................................................
 * Return the value of a double keyword
 */
double Program::getDoubleParameter (string key) 
{
  string sp = getStringParameter(key);
  return atof(sp.c_str());
}

/**.......................................................................
 * Return the value of a boolean keyword
 */
bool Program::getBooleanParameter (string key) 
{
  static string yes = "tTyY1";
  static string no  = "fFnN0";
  string sp = getStringParameter(key);
  size_t nt = yes.find(sp[0]);
  size_t nf = no.find(sp[0]);
  
  if (nt >= 0 && nt < yes.length()) return true;
  if (nf >= 0 && nf < yes.length()) return false;
  cerr << "Syntax error for boolean flag " << key << "=" << sp << endl;
  cerr << "Expecting one of " << yes << " or " << no << endl;
  exit(1);
}

/**.......................................................................
 * Count the number of keywords
 */
int Program::count(string key) 
{
  int i = findKey(key);
  if (i<0) {
    cerr <<  key << " Illegal key" << endl;
    exit(1);
  }
  return keys_[i].count;
}

/**.......................................................................
 * findKey:  scan valid keywords and return index (>=0) for a keyname
 *           return -1 if none found
 *	     Optionally match can be done in minimum match mode
 */
int Program::findKey(const string name)
{
  int i, j, l, count, last;
  
  if (nkeys_ <= 0) return -1;                /* quick: no keywords at all */
  for (i = 0; i < nkeys_; i++)             /* try exact match */
    if (keys_[i].key == name) return i;
  
#if 0
  // lotsa compile booboos still....
#if defined(MINMATCH)
  l = strlen(name);                       /* try minimum match */
  count = 0;                              /* count # matches on */
  for (i=1; i < nkeys_; i++) {               /* any of the program keys */
    if (keys_[i].key.compare(name,l)) {
      last = i;
      count++;
    }
  }
  if (count==1) {
    cerr << "Resolving partially matched keyword " <<  name 
	 << " into ", keys_[last].key << "=" << endl;
    return last;
  } else if (count > 1) {
    cerr << "### Minimum match failed, found: " ;
    for (j=0; j < nkeys_; j++)
      if (strncmp(keys_[j].key,name,l)==0)
	cerr << keys_[j].key << " ";
    cerr << endl;
    cerr << "Ambiguous keyword %s=" << name << endl;
    exit(1);
  }
#endif
#endif
  
  return -1;          /* if all else fails: return -1 */
}

// local helper functions because we don't have a fancy string class (yet)
// these come pretty much straight from NEMO

static string parname(string arg)
{
  int ieq = arg.find_first_of('=');
  if (ieq < 0)
    cerr << "parname " << arg << " has no equals" << endl;
  return arg.substr(0,ieq);
}

/*
 * PARVALUE: extract value from name=value string, skipping initial whitespace
 *           if no value, pointer to 0 (last element) is returned!
 *           if no key, pointer to (local copy of ) first element is returned
 *  ???      BUG:  when HELPVEC is not defined, helpvec also returned
 *           Note: returns unsafe pointer into the input string
 *           
 */
static string parvalue(string arg)
{
  int ieq = arg.find_first_of('=');
  if (ieq < 0)
    cerr << "parvalue " << arg << " has no equals" << endl;
  return arg.substr(ieq+1);
}

/* 
 * PARHELP: extract help from a defv[] "keyword=value\n help" string
 *          If no help part (or \n), returns zero string.
 *          Note: returns pointer into original string, assumed
 *          to be R/O (might be in text space of code)
 */
static string parhelp(string arg)
{
  int ieq = arg.find_first_of('\n');
  if (ieq < 0)
    cerr << "parhelp " << arg << " has no newline" << endl;
  return arg.substr(ieq+1);
}

static bool getbvalue(string sp) {
  static string yes = "tTyY1";
  static string no  = "fFnN0";
  size_t nt = yes.find(sp[0]);
  size_t nf = no.find(sp[0]);
  
  if (nt >= 0 && nt < yes.length()) return true;
  if (nf >= 0 && nf < yes.length()) return false;
  cerr << "Syntax error for boolean flag " << sp << endl;
  cerr << "Expecting one of " << yes << " or " << no << endl;
  exit(1);
  
}

/**.......................................................................
 * Get a pointer to a static instance of the Program class
 */
Program &Program::getProgram()     //  get that only instance 
{
  static Program *program = 0;     //  the singleton instance

  if (program == 0)  {             //  first time around
    //cerr << "Program start " << endl;
    program = new Program();
  } 

  assert (program != 0);
  
  return *program;
}
