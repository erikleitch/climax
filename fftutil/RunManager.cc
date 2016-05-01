#include "gcp/datasets/DataSet1D.h"

#include "gcp/fftutil/RunManager.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/util/ChisqVariateGaussApprox.h"
#include "gcp/util/GaussianVariate.h"
#include "gcp/util/JointGaussianVariate.h"
#include "gcp/util/RangeParser.h"

#include "cpgplot.h"

#include <stack>

using namespace std;
using namespace gcp::datasets;
using namespace gcp::util;
using namespace gcp::models;

/**.......................................................................
 * Constructor.
 */
RunManager::RunManager()
{
  runMarkov_           = true;

  compChisq_           = false;
  displayDataSets_     = false;
  writeDataSets_       = false;
  datasetsInitialized_ = false;
  modelsInitialized_   = false;
  loadOutputFile_      = false;
  varPlot_             = "hist";
  debug_               = false;
  interactive_         = false;
  incBurnIn_           = false;
  printConvergence_    = false;
  runtoConvergence_    = false;
  targetVariance_ = 0.01;
  updateMethod_        = 1;

  genData_             = false;

  nBurn_               = 1000;
  nTry_                = 10000;

  nTimesAtThisPoint_   = 0;

  nModelThread_        = 0;
  modelPool_           = 0;
  modelCpus_.resize(0);

  nDataThread_         = 0;
  dataPool_            = 0;
  dataCpus_.resize(0);

  pgplotDev_           = "/xs";
  nBin_                = 30;
  runType_             = 1;

  stat_                = Model::STAT_MODE;
  nSigma_              = 1.0;

  //------------------------------------------------------------
  // Add 'parameters' that we know about
  //------------------------------------------------------------

  general_.addParameter("General notes on syntax", DataType::STRING, 
			"In the notes that follow, the characters {} are used to denote "
			"optional parameters, and are not to be interpreted literally. "
			"Likewise literal strings are encapsulated between single quotes '', "
			"and the quotation marks are not to be interpreted as part of the "
			"string.  Finally, all valid lines in CLiMax run files must be "
			"terminated with a semicolon, which is not included in the literal "
			"strings listed below.  Thus when an instruction says Use like: 'seed = "
			"val' it means that your run file should contain a line that looks "
			"like \n\n\t\t\tseed = 3;\n\n");

  general_.addParameter("General notes on runfiles", DataType::STRING, 
			"The code is intended to be largely self-documenting.  What this means "
			"in practice is that any erroneous statement should generate a help "
			"message showing you what the valid options are for that statement. "
			"For example, if you don't know what model types are available in the "
			"addmodel command, you can type a random string in the type field, and "
			"the code should generate an error message telling you what the valid "
			"options are.  Likewise for unit strings where appropriate.\n\n"

			"Lastly, some options are mutually exclusive.  For example, with no other arguments, CLiMax will try to run "
			"a Markov chain for the models and datasets that are currently defined, but if display=true "
			"is specified, CLiMax will attempt to display the model, data and residuals.  In this case, "
			"fixed values must be specified for all models, which obviously precludes running a Markov chain.\n");

  docs_.addParameter("adddataset",   DataType::STRING, "Directive to add a data set.  Use like: 'adddataset name=nameStr type=typeStr'");
  docs_.addParameter("addmodel",     DataType::STRING, "Directive to add a model.  Use like: 'addmodel name=nameStr type=typeStr'");
  docs_.addParameter("remmodel",     DataType::STRING, "Directive to remove a model from the data prior to fitting or displaying.  Use like: 'remmodel name=nameStr type=typeStr'");
  docs_.addParameter("addobs",       DataType::STRING, "Directive to add an observation object that can be used for simulation.  Use like 'addobs name=nameStr'");
  docs_.addParameter("comp",         DataType::BOOL,   "If set to 'true', chi-squared will be reported for any defined datasets and models");
  docs_.addParameter("debug",        DataType::BOOL,   "If set to 'true', print debugging information");
  docs_.addParameter("help",         DataType::STRING, "This output");
  docs_.addParameter("display",      DataType::BOOL,   "If set to 'true', the datasets will be displayed, along with models and residuals");
  docs_.addParameter("writedata",    DataType::BOOL,   "If set to 'true', any requested data, model or residual datasets will be written out");
  docs_.addParameter("gendata",      DataType::STRING, "Directive to generate simulated data for defined models.  Use like 'gendata file=fileName dataset=name {sigma=val} "
		     "where fileName = the output file name, dataset = the dataset for which to generate simulated data, and "
		     "sigma = (optional) value to assign fixed errors, where the dataset type is written to use them.  Note that "
		     "visdatasets use the observation object for simulating data, specified via the 'addobs' directive.");
  docs_.addParameter("include",      DataType::STRING, "A directive to include another run file in the current one.  Use like: 'include filename'");
  docs_.addParameter("d.exclude(m)", DataType::STRING, "A directive to exclude a model from a dataset. Use like 'd.exclude(m)' where d = dataset, m = model to exclude");
  docs_.addParameter("nbin",         DataType::UINT,   "The number of bins into which the data will be binned for Markov run histograms.  Use like 'nbin = 100'");
  docs_.addParameter("nburn",        DataType::UINT,   "The length of the burn-in sequence for the Markov chain.  Use like 'nburn = 3000'. The jumping distribution will be tuned during the "
		     "burn-in period, and these samples will be discarded from any output file.");
  docs_.addParameter("ntry",         DataType::UINT,   "The total length of the Markov chain to run.  Use like 'ntry = 10000'");
  docs_.addParameter("nmodelthread", DataType::UINT,   "The number of threads in the model pool.  Use like 'nmodelthread = 10'");
  docs_.addParameter("modelcpus",    DataType::STRING, "A list of cpus to which the model threads should be bound.  Use like 'modelcpus = 1,2,3'");
  docs_.addParameter("ndatathread",  DataType::UINT,   "The number of threads in the data pool.  Use like 'ndatathread = 10'");
  docs_.addParameter("datacpus",     DataType::STRING, "A list of cpus to which the data threads should be bound.  Use like 'datacpus = 1,2,3'");
  docs_.addParameter("output",       DataType::STRING, "If specified, the output file for Markov chain runs.  Use like 'output file=fileName'");
  docs_.addParameter("incburnin",    DataType::BOOL,   "If true, include burn-in samples in plots/output file (default is false)");
  docs_.addParameter("varplot",      DataType::STRING, "The type of variable plot to produce.  One of: 'hist' (default), 'line' or 'power'");
  docs_.addParameter("seed",         DataType::UINT,   "If specified, the random number generator will be explicitly seeded with this value.  Use like 'seed = value'");
  docs_.addParameter("nsigma",       DataType::DOUBLE, "The width of the statistic specified by the 'stat' variable, in units of gaussian sigma (default is 1.0)");
  docs_.addParameter("stat",         DataType::STRING, "The statistic to display for fit parameters.  Recognized values are:\n\n"
		     "     maxl  - Maximum likelihood value and nsigma confidence interval (the default)\n"
		     "     upper - Upper limit (nsigma)\n"
		     "     lower - Lower limit (nsigma)\n\n"
		     "where nsigma is specified by the the 'nsigma' parameter");
  docs_.addParameter("psfont",       DataType::STRING, "The Postscript font to use for postscript plots.  Use like 'psfont = Courier-Bold'");
  docs_.addParameter("dev",          DataType::STRING, "The PGPLOT output device for Markov run plots.  Use like 'dev = /xs'");
  docs_.addParameter("interactive",  DataType::BOOL,   "If set to 'true', summary plots will be interactive.");
  docs_.addParameter("store",        DataType::BOOL,   "If set to 'true', accepted values will be stored internally (so that histogram plots can be generated, etc).  Use like 'store = true' (the default). "
		     "For long runs, this can of course introduce a significant memory footprint, so you may wish to set this "
		     "to false if you are writing the output chain to disk anyway");
  docs_.addParameter("load",         DataType::STRING, "Directive to load a file containing an output Markov chain.  Use like: 'load file=fileName {discard=n} {name=str}', "
		     "where file = the name of the file containing the output chain, discard = the (optional) number of sample to discard on read-in, "
		     "and name = (optional) name to assign the model being read in (ignored if the chain file includes a model name)");
  docs_.addParameter("d.par = val",  DataType::STRING, "Assignment of a dataset parameter value, where d = dataset name, par = parameter name, val = value.  On error, each dataset will list "
		     "valid parameter names");
  docs_.addParameter("m.var = val {units}",  DataType::STRING, "Assignment of a model variable value, where m = model name, var = variable name, val = value. "
		     "On error, each model will list valid variable names.  If a variable is a physical quantity, the units are required. "
		     "On error, each variable will list recognized units");
  docs_.addParameter("m.var = min:max {unitstr}",             DataType::STRING, "Assignment of a uniform prior for a model variate, where m = model name, var = variable name, min:max = the range of the prior.");
  docs_.addParameter("m.var = start:min:max {unitstr}",       DataType::STRING, "Assignment of a uniform prior for a model variate, and a starting mean for the jumping distribution");
  docs_.addParameter("m.var = mean +- sigma {unitstr}",       DataType::STRING, "Assignment of a Gaussian prior for a model variate");
  docs_.addParameter("m.var = start:mean +- sigma {unitstr}", DataType::STRING, "Assignment of a Gaussian prior for a model variate, and a starting mean for the jumping distribution");
  docs_.addParameter("m.var = mean +- sigma x min:max {unitstr}", DataType::STRING, "Assignment of a truncated Gaussian prior for a model variate.  Infinities are allowed in the specification: i.e., if min is '-inf', no lower bound will be applied to the prior.  Likewise if max = 'inf', no upper bound will be applied to the prior.");
  docs_.addParameter("m.var[display] = true",                 DataType::STRING, "Set to 'true' or 'false' to control display of this variate in histogram plots");
  docs_.addParameter("m.var[order] = n",                      DataType::STRING, "Control the order of display of variates in histogram plots");
  docs_.addParameter("m.var[range] = min:max {unitstr}",      DataType::STRING, "Control the range in which variates in histogram plots will be displayed");
  docs_.addParameter("m.var[units] = unitstr",                DataType::STRING, "Control the units used for a variate");
  docs_.addParameter("printconvergence",                      DataType::BOOL,   "If true, estimate chain convergence");
  docs_.addParameter("runtoconvergence",                      DataType::BOOL,   "If true, run to convergence or 'ntry', whichever comes first");
  docs_.addParameter("update",                                DataType::UINT,   "Which update method to try?  Default is 0.  Mean update is 1");
  docs_.addParameter("targetvariance",                        DataType::DOUBLE, "Fractional variance to target (a number < 1).  Default is 0.01.  Used as a convergence criteria with 'convergence' and 'runtoconvergence' keywords");
}

/**.......................................................................
 * Destructor.
 */
RunManager::~RunManager() 
{
  if(modelPool_) {
    delete modelPool_;
    modelPool_ = 0;
  }

  if(dataPool_) {
    delete dataPool_;
    dataPool_ = 0;
  }

  //------------------------------------------------------------
  // Delete any persistent memory kept by FFTW under the hood
  //------------------------------------------------------------

  fftw_cleanup();
}

/**.......................................................................
 * Initialize convenience variables for running Markov chains
 */
void RunManager::initializeMarkovSpecificVariables()
{
  acceptFracLowLim_  = mm_.nVar() > 1 ? 0.15 : 0.4;
  acceptFracHighLim_ = mm_.nVar() > 1 ? 0.35 : 0.6;

  foundUpdate_ = false;
  nPerUpdate_  = 100;
  nUpdate_     = 0;
  iLastUpdate_ = 0;
  nAcceptedSinceLastUpdate_ = 0;
  nTrySinceLastUpdate_      = 0;

  nConverge_   = nTry_ > 1000 ? 1000 : nTry_ / 10;
  if(nConverge_ == 0)
    nConverge_ = 1;
  
  overallTime_ = 0.0;
  sampleTime_  = 0.0;
  likeTime_    = 0.0;
  tuneTime_    = 0.0;
}

/**.......................................................................
 * Main method of this object.  Execute a run file.
 */
void RunManager::run() 
{
  parseFile(runFile_);

  //------------------------------------------------------------
  // Runfiles can do several things:
  // 
  // If requested to generate fake data for a particular dataset type
  // and model, do it now
  //------------------------------------------------------------

  if(genDataLines_.size() > 0) {

    for(unsigned i=0; i < genDataLines_.size(); i++) {
      String str(genDataLines_[i]);
      getGenDataArgs(str);

      if(mm_.genData()) {
	mm_.generateFakeData();
      }
    }

    runMarkov_ = false;
  }

  //------------------------------------------------------------
  // If computing chisq, do that now
  //------------------------------------------------------------
  
  if(compChisq_) {
    computeChisq();
    runMarkov_ = false;
  }    

  //------------------------------------------------------------
  // If displaying datasets, do that now
  //------------------------------------------------------------
  
  if(displayDataSets_) {

    PgUtil::open(pgplotDev_);
    PgUtil::setOverplot(true);
    PgUtil::setWedge(true);
    PgUtil::setCharacterHeight(0.5);
    PgUtil::setInteractive(interactive_);

    dm_.initializeForDisplay();
    assertCurrentModel("display", true);

    populatePlotWithDataSets(0.1, 0.9, 0.1, 0.9);
    PgUtil::close();

    runMarkov_ = false;
  }

  //------------------------------------------------------------
  // If loading data, do that now
  //------------------------------------------------------------
  
  if(loadOutputFile_) {
    createMarkovDisplay();
    runMarkov_ = false;
  }

  //------------------------------------------------------------
  // Else by default run a Markov chain with the datasets and models
  // that were specified
  //------------------------------------------------------------

  if(runMarkov_) {

    try {
      if(nTry_ <= nBurn_)
	ThrowSimpleColorError("You have specified ntry = " << nTry_ 
			      << " and nburn = " << nBurn_ 
			      << ".  ntry should be > nburn", "red");

      mm_.setThreadPool(modelPool_);
      mm_.initializeForMarkovChain(nTry_, incBurnIn_ ? nTry_ : nTry_ - nBurn_, runFile_);

      //------------------------------------------------------------
      // Sanity check that any datasets were initialized
      //------------------------------------------------------------

      if(dm_.dataSetMap_.size() == 0)
	COUTCOLOR(std::endl << "Warning: You are running a Markov chain but no datasets have been loaded" << std::endl, "red");

      runMarkov();

      if(mm_.nAccepted_ > 0) {
	createMarkovDisplay();
      } else {
	ThrowSimpleColorError("No samples were accepted by the Markov chain", "red");
      }
    } catch(...) {
    }
  }

  //------------------------------------------------------------
  // If writing out data, do that now
  //------------------------------------------------------------

  if(writeDataSets_) {
    dm_.initializeForDisplay();

    if(runMarkov_)
      assertBestFitModel(true);
    else
      assertCurrentModel("writing", true);

    dm_.writeData();
  }

  return;
}

void RunManager::setRunFile(std::string runFile)
{
  String str(runFile);
  str.expandTilde();

  runFile_ = str.str();
}

void RunManager::setNtry(unsigned nTry)
{
  nTry_ = nTry;
}

void RunManager::setNburn(unsigned nBurn)
{
  nBurn_ = nBurn;
}

void RunManager::setNModelThread(unsigned nModelThread)
{
  nModelThread_ = nModelThread;
}

void RunManager::setNDataThread(unsigned nDataThread)
{
  nDataThread_ = nDataThread;
}

void RunManager::initializeThreadPools()
{
  if(nModelThread_ > 0) {
    if(modelCpus_.size() > 0)
      modelPool_ = new ThreadPool(nModelThread_, modelCpus_);
    else
      modelPool_ = new ThreadPool(nModelThread_);

    modelPool_->spawn();
  }
  
  if(nDataThread_ > 0) {
    if(dataCpus_.size() > 0)
      dataPool_ = new ThreadPool(nDataThread_, dataCpus_);
    else
      dataPool_ = new ThreadPool(nDataThread_);

    dataPool_->spawn();
  }
}
 
/**.......................................................................
 * Parse a config file
 */
void RunManager::parseFile(std::string fileName)
{
  //------------------------------------------------------------
  // Pre-process the parameter file.  This looks for directives that
  // affect subsequent parsing, like thread-pool allocation
  //------------------------------------------------------------

  processFile(fileName, true);

  //------------------------------------------------------------
  // Now initialize thread pools
  //------------------------------------------------------------

  initializeThreadPools();

  //------------------------------------------------------------
  // And process the rest of the file
  //------------------------------------------------------------

  processFile(fileName, false);

  //------------------------------------------------------------
  // Only once models have been defined can we process model-dependent
  // lines
  //------------------------------------------------------------

  processModelDependentLines();
  processModelLoadDependentLines();

  if(loadOutputFile_) {
    mm_.addDerivedVariates();
    mm_.fillDerivedVariates();
  }

  //------------------------------------------------------------
  // Now load any data sets that were specified
  //------------------------------------------------------------

  loadData();

  //------------------------------------------------------------
  // Check for position information.  This checks if ra/dec parameters
  // have been specified from an external source, and if not sets
  // those parameters to values that have been read in from the data,
  // so that they will be available to the calling interface
  //------------------------------------------------------------

  dm_.checkPosition();

  //------------------------------------------------------------
  // Only once data sets have been loaded, can we make variable
  // assignments that depend on them
  //------------------------------------------------------------

  processDatasetDependentLines();

  //------------------------------------------------------------
  // Check for position information one more time, in case a parameter
  // was specified externally as a symbolic name that wasn't defined
  // until processDatasetDependentLines() was called.  I.e., this
  // allows setting the ra/dec for a dataset to the ra/dec of another
  // dataset
  //------------------------------------------------------------

  dm_.checkPosition(true);

  //------------------------------------------------------------
  // Now call display method if any datasets should be displayed
  //------------------------------------------------------------

  try {
    dm_.displayIfRequested();
  } catch(Exception& err) {
    XtermManip xtm;
    COUT(COLORIZE(xtm, "red", "Error displaying datasets: ") << std::endl << err.what() <<
	 COLORIZE(xtm, "red", std::endl << "Attempting to continue" << std::endl));
  }

  //------------------------------------------------------------
  // Finally, check model setups for sense
  //------------------------------------------------------------

  if(!loadOutputFile_)
    mm_.checkSetup();

  //------------------------------------------------------------
  // Remove any models 
  //------------------------------------------------------------

  removeModels();
}

/**.......................................................................
 * Called to process a run file
 */
void RunManager::processFile(std::string fileName, bool preprocess)
{
  static String line;
  unsigned nLine=0;
  bool processNext = true;
  std::stack<bool> processNextTable;
  processNextTable.push(true);
  XtermManip xtm;

  std::ifstream fin;
  fin.open(fileName.c_str(), ios::in);
  
  if(!fin) {
    ThrowSimpleColorError(std::endl << "Unable to open file: " << fileName, "red");
  }
  
  try {

    //------------------------------------------------------------
    // Iterate through the file
    //------------------------------------------------------------
    
    while(!fin.eof()) {

      line.initialize();
      getline(fin, line.str());

      ++nLine;

      //------------------------------------------------------------
      // Handle if directives
      //------------------------------------------------------------

      if(line.contains("if") && line.contains("{")) {
	String ifClause = line.findNextInstanceOf("(", true, ")", true, true);
	ifClause.strip(" ");
	processNextTable.push(ifClause.toBool());

	String lineTest = line.findNextInstanceOf("{", true, " ", false, true);
	String rem = line.remainder();
	rem.strip(' ');
	rem.advanceToNextNonWhitespaceChar();

	if(!rem.isEmpty())
	  ThrowSimpleError(COLORIZE(xtm, "red", std::endl << "Found extraneous characters following a valid parenthetical, while processing line " << nLine << " of " << fileName << ":" 
				    << std::endl << std::endl << "'" << rem << "'") << 

			   COLORIZE(xtm, "green", std::endl << std::endl << "Each parenthetical clause should appear on a separate line"));

	line = lineTest;
      }

      //------------------------------------------------------------
      // Handle else directives
      //------------------------------------------------------------

      if(line.contains("}") && line.contains("else") && line.contains("{")) {

	if(processNextTable.size() > 1) {
	  String lineTest = line.findNextInstanceOf("{", true, " ", false, true);
	  String rem = line.remainder();
	  rem.strip(' ');
	  rem.advanceToNextNonWhitespaceChar();
	  
	  if(!rem.isEmpty())
	    ThrowSimpleError(COLORIZE(xtm, "red", std::endl << "Found extraneous characters following a valid else statement, while processing line " << nLine << " of " << fileName << ":" 
				      << std::endl << std::endl << "'" << rem << "'") << 
			     
			     COLORIZE(xtm, "green", std::endl << std::endl << "Each else statement should appear on a separate line"));
	  
	  line = lineTest;
	  
	  // If everything checks out, we pop the last condtional, and
	  // push a new one onto the stack:
	  
	  bool process = processNextTable.top();
	  processNextTable.pop();
	  processNextTable.push(!process);
	} else {
	  ThrowSimpleColorError("Encountered 'else' statement with no preceding 'if' clause", "red");
	}
      }

      //------------------------------------------------------------
      // Handle if/else terminations
      //------------------------------------------------------------

      if(processNextTable.size() > 1 && line.contains("}")) {
	String lineTest = line.findNextInstanceOf(" ", false, "}", true, true);
	processNextTable.pop();

	String rem = line.remainder();
	rem.strip(' ');
	rem.advanceToNextNonWhitespaceChar();

	if(!rem.isEmpty())
	  ThrowSimpleError(COLORIZE(xtm, "red", std::endl << "Found extraneous characters following a valid parenthetical, while processing line " << nLine << " of " << fileName << ":" 
				    << std::endl << std::endl << "'" << rem << "'") <<
			   COLORIZE(xtm, "green", std::endl << std::endl << "Each parenthetical clause should be terminated on a separate line"));

	line = lineTest;
      }

      processNext = processNextTable.top();

      //------------------------------------------------------------
      // Done handling if directives.  Check for correct termination
      // of lines
      //------------------------------------------------------------

      line.advanceToNextNonWhitespaceChar();

      String testStr(line);
      testStr.strip(' ');

      if(testStr.isEmpty())
	continue;

      if(!line.isEmpty() && !isComment(line) && !line.contains(";")) {
	ThrowSimpleError(COLORIZE(xtm, "red", std::endl << "Found an unterminated line while processing line " << nLine << " of " << fileName << ":" 
				  << std::endl << std::endl << "'" << line << "'") << 
			 COLORIZE(xtm, "green", std::endl << std::endl << "All valid directives should end in ';'"
				  << std::endl << std::endl << "(Comments can be inserted in a runfile by prepending '//' to the line)"));
      }
      
      line = line.findNextInstanceOf("", false, ";", true, true);
      line.advanceToNextNonWhitespaceChar();
      
      //------------------------------------------------------------
      // Don't process comment lines
      //------------------------------------------------------------
      
      if(line[0] == '/')
	continue;

      try {

	if(processNext) {
	  if(preprocess)
	    preProcessLine(line);
	  else
	    processLine(line);
	}

      } catch(Exception& err) {
	XtermManip xtm;
	ThrowSimpleError(COLORIZE(xtm, "red", std::endl << err.what() << std::endl) <<
			 COLORIZE(xtm, "red", std::endl << "while processing line '" << line << "'" << std::endl));
      }

    }
  } catch(Exception& err) {
    fin.close();
    ThrowSimpleError(err.what());
  }
  
  fin.close();
}

bool RunManager::isComment(String& line)
{
  return line[0] == '/' && line[1] == '/';
}

/**.......................................................................
 * Called to process a single line from a config file
 */
void RunManager::processLine(String& line)
{
  String copy(line);
  String firstToken = copy.findNextInstanceOf(" ", false, "=", true);
  firstToken.strip(' ');

  //------------------------------------------------------------
  // Add model line
  //------------------------------------------------------------

  if(line.contains("addmodel")) {

    addModel(line, false);

    //------------------------------------------------------------
    // Add dataset line
    //------------------------------------------------------------

  } else if(line.contains("adddataset")) {

    addDataset(line);

    //------------------------------------------------------------
    // Add obs line
    //------------------------------------------------------------

  } else if(line.contains("addobs")) {

    std::string name, type;

    getObsName(line, name);
    checkIfNameAlreadyExists(name);

    COUT("Adding obs by name '" << name << "'");
    ObsInfo* obs = om_.addObs(name);

    //------------------------------------------------------------
    // Remove model line
    //------------------------------------------------------------

  } else if(line.contains("remmodel")) {

    addModel(line, true);

    //------------------------------------------------------------
    // Output file line
    //------------------------------------------------------------
      
  } else if(line.contains("output ") && !line.contains("displayoutput")) {

    getOutputArgs(line);
	
    //------------------------------------------------------------
    // Generate fake data line
    //------------------------------------------------------------
      
  } else if(line.contains("gendata")) {
      
    genData_ = true;
    genDataLines_.push_back(line.str());
      
    //------------------------------------------------------------
    // Compute chisq line
    //------------------------------------------------------------
      
  } else if(line.contains("comp") && !line.contains("compton")) {
      
    getComputeChisqArgs(line);

    //------------------------------------------------------------
    // Get convergence parameter
    //------------------------------------------------------------
      
  } else if(firstToken == "printconvergence") {
      
    printConvergence_ = (getStrippedVal(line).toLower().str() == "true");

    //------------------------------------------------------------
    // Get runtoconvergence parameter
    //------------------------------------------------------------
      
  } else if(firstToken == "runtoconvergence") {
      
    runtoConvergence_ = (getStrippedVal(line).toLower().str() == "true");

    //------------------------------------------------------------
    // Get convergencevar parameter
    //------------------------------------------------------------
      
  } else if(firstToken == "targetvariance") {
      
    targetVariance_ = getStrippedVal(line).toDouble();

    if(targetVariance_ >= 1.0)
      ThrowSimpleColorError("Invalid fractional target variance: " << targetVariance_ << ".  Should be < 1", "red");

    //------------------------------------------------------------
    // Get incburnin parameter
    //------------------------------------------------------------
      
  } else if(line.contains("incburnin")) {
      
    incBurnIn_ = (getStrippedVal(line).toLower().str() == "true");

    //------------------------------------------------------------
    // Get debug parameter
    //------------------------------------------------------------
      
  } else if(line.contains("debug") && !line.contains(".debug")) {
      
    debug_ = (getStrippedVal(line).toLower().str() == "true");

    //------------------------------------------------------------
    // Get interactive parameter
    //------------------------------------------------------------

  } else if(line.contains("interactive") && !line.contains(".interactive")) {
      
    interactive_ = (getStrippedVal(line).toLower().str() == "true");

    //------------------------------------------------------------
    // Power spectrum plot
    //------------------------------------------------------------
      
  } else if(firstToken.contains("varplot")) {

    varPlot_ = getStrippedVal(line).toLower().str();
      
    String varStr(varPlot_);
    if(varStr != "hist" && varStr != "line" && varStr != "power")
      ThrowSimpleColorError("Unrecognized varplot type: '" << varPlot_ << "'. Should be one of 'hist', 'line' or 'power'", "red");

    //------------------------------------------------------------
    // Display datsets line
    //------------------------------------------------------------
      
  } else if(firstToken == "display") {
      
    getDisplayDataSetsArgs(line);
      
  } else if(firstToken == "writedata") {
      
    getWriteDataSetsArgs(line);
      
  } else if(line.contains("update")) {
      
    updateMethod_ = getStrippedVal(line).toInt();

    //------------------------------------------------------------
    // Include line
    //------------------------------------------------------------
    
  } else if(line.contains("include") && !line.contains("includeifs")) {
    
    // Skip the include token
    
    line.findNextString();

    String fileName = line.findNextString();
    fileName.strip(' ');
    fileName.expandTilde();

    processFile(fileName.str(), false);
      
    //------------------------------------------------------------
    // Exclude line
    //------------------------------------------------------------
    
  } else if(line.contains("exclude") && !line.contains("excludeants") && !line.contains("excludeifs")) {
      
    if(!datasetsInitialized_) {
      datasetDependentLines_.push_back(line.str());
      return;
    } else {
      getExcludeArgs(line);
    }

  } else if(line.contains("load ")) {
      
    getLoadDataArgs(line);

    //------------------------------------------------------------
    // Other lines we recognize
    //------------------------------------------------------------
    
  } else if(!line.isEmpty()) {
    
    //------------------------------------------------------------
    // Token = value statement
    //------------------------------------------------------------
    
    if(line.contains("=")) {

      String tok, val;
      getTokVal(line, tok, val);

      //------------------------------------------------------------
      // Intercept global variable assignments
      //------------------------------------------------------------

      if(!tok.contains(".")) {

	//------------------------------------------------------------
	// Ignore pre-process directives
	//------------------------------------------------------------

	if(tok.contains("nmodelthread")) {
	  return;
	} else if(tok.contains("ndatathread")) {
	  return;

	  //------------------------------------------------------------
	  // Explicitly set nburn
	  //------------------------------------------------------------

	} else if(tok.contains("nburn")) {
	  unsigned nBurnOld = nBurn_;
	  unsigned nBurnNew = val.toInt();

	  nBurn_ = nBurnNew;
	  return;

	  //------------------------------------------------------------
	  // Explicitly set ntry
	  //------------------------------------------------------------

	} else if(tok.contains("ntry")) {
	  unsigned nTryOld = nTry_;
	  unsigned nTryNew = val.toInt();

	  nTry_ = nTryNew;
	  return;

	  //------------------------------------------------------------
	  // Explicitly seed the sampler random number generator
	  //------------------------------------------------------------

	} else if(tok.contains("seed")) {
	  Sampler::seed(val.toInt());
	  return;

	  //------------------------------------------------------------
	  // Get the statistic to display
	  //------------------------------------------------------------

	} else if(tok.contains("stat")) {

	  String statVal = String::toLower(val.str());

	  if(statVal.contains("maxl"))
	    stat_ = Model::STAT_MODE;
	  else if(statVal.contains("upper"))
	    stat_ = Model::STAT_UPPER_LIMIT;
	  else if(statVal.contains("lower"))
	    stat_ = Model::STAT_LOWER_LIMIT;
	  else {
	    ThrowSimpleError("Unrecognized statistic: " << val 
			     << ".  Should be one of: 'maxl', 'upper' or 'lower'");
	  }

	  return;

	  //------------------------------------------------------------
	  // Get the width of the requested statistic, in sigma
	  //------------------------------------------------------------

	} else if(tok.contains("nsigma")) {

	  nSigma_ = val.toDouble();
	  return;

	  //------------------------------------------------------------
	  // Set the pgplot device for Markov chain plots
	  //------------------------------------------------------------

	} else if(tok.contains("dev")) {
	  val.strip(' ');
	  pgplotDev_ = val.str();
	  return;

	  //------------------------------------------------------------
	  // Set the postscript font for postscript output files
	  //------------------------------------------------------------

	} else if(tok.contains("psfont")) {
	  val.strip(' ');
	  PgUtil::setPostscriptFontName(val.str());
	  return;

	  //------------------------------------------------------------
	  // Set the number of bins for markov histograms
	  //------------------------------------------------------------

	} else if(tok.contains("nbin")) {
	  nBin_ = val.toInt();
	  return;

	} else if(tok.contains("runtype")) {
	  runType_ = val.toInt();
	  return;

	} else if(tok.contains("store")) {
	  val.strip(' ');
	  mm_.setStore(val.str() == "true");
	  return;

	} else if(tok.contains("modelcpus")) {
	  return;

	} else if(tok.contains("datacpus")) {
	  return;

	} else {
	  std::ostringstream os;
	  docs_.listParameters(os);
	  ThrowSimpleColorError(std::endl << "Unrecognized directive '" << line << "'" << std::endl << std::endl 
				<< "Recognized directives are: " << std::endl << std::endl << os.str(), "red");
	}
      }

      //------------------------------------------------------------
      // Don't proceed if this assignment depends on a dataset
      // parameter that has not yet been read in
      //------------------------------------------------------------

      if(assignmentDependsOnDataset(val) && !datasetsInitialized_) {
	datasetDependentLines_.push_back(line.str());
	return;
      }

      //------------------------------------------------------------
      // Don't proceed if this assignment depends on a model that has
      // (potentially) not yet been loaded
      //------------------------------------------------------------

      if(assignmentDependsOnModelLoad(tok) && !modelsInitialized_) {
	modelLoadDependentLines_.push_back(line.str());
	return;
      }

      //------------------------------------------------------------
      // Check for different types of tokens.
      //
      // Correlations are specified like:          var1 * var2 = val
      // Variable assignmemnts are specified like: var = val
      // 
      //------------------------------------------------------------

      if(tok.contains("*")) {

	std::string varName1, varName2;
	parseCorrelationVariables(tok, varName1, varName2);
	mm_.addComponentCorrelation(varName1, varName2, val.toDouble());

	// Else this is either an increment or an assignment

      } else if(line.contains("+=")) {
	parseVariableIncrement(tok, val);
      } else {
	parseVariableAssignmentNew(line, tok, val);
      }
	
    } else if(line.contains("help")) {
      std::ostringstream os1;
      std::ostringstream os2;
      general_.listParameters(os1, false);
      docs_.listParameters(os2, true);

      XtermManip xtm;
      ThrowSimpleColorError(std::endl 
			    << "A few general notes about running CLiMax: " << std::endl << std::endl << os1.str() << std::endl
			    << COLORIZE(xtm, "green", "Recognized directives are:") << std::endl << std::endl << os2.str(), "green");
    } else {
      std::ostringstream os;
      docs_.listParameters(os);
      ThrowSimpleColorError(std::endl 
			    << "Unrecognized directive '" << line << "'" << std::endl << std::endl 
			    << "Recognized directives are: " << std::endl << std::endl << os.str(), "red");
    }
      
  }
}

/**.......................................................................
 * Called to pre-process a single line from a config file
 */
void RunManager::preProcessLine(String& line)
{
  String copy(line);
  String firstToken = copy.findNextInstanceOf(" ", false, "=", true);
  firstToken.strip(' ');

  //------------------------------------------------------------
  // Include line
  //------------------------------------------------------------
    
  if(line.contains("include") && !line.contains("includeifs")) {
    
    // Skip the include token
    
    line.findNextString();
    
    String fileName = line.findNextString();
    fileName.strip(' ');
    fileName.expandTilde();
    
    processFile(fileName.str(), true);

    //------------------------------------------------------------
    // Token = value statement
    //------------------------------------------------------------

  } else if(line.contains("=")) {
    
    String tok, val;
    getTokVal(line, tok, val);
    
    //------------------------------------------------------------
    // Explicitly set nburn
    //------------------------------------------------------------
    
    if(tok.contains("nmodelthread")) {
      setNModelThread(val.toInt());

    } else if(tok.contains("ndatathread")) {
      setNDataThread(val.toInt());

    } else if(tok.contains("modelcpus")) {
      val.strip(' ');
      
      RangeParser parser;
      std::ostringstream os;
      String cpuStr(val);
      
      if(!cpuStr.contains("*")) {
	os << "[" << cpuStr.str() << "]";
	cpuStr = os.str();
      }
      
      modelCpus_ = parser.extractIndexRange(cpuStr);

    } else if(tok.contains("datacpus")) {
      val.strip(' ');
      
      RangeParser parser;
      std::ostringstream os;
      String cpuStr(val);
      
      if(!cpuStr.contains("*")) {
	os << "[" << cpuStr.str() << "]";
	cpuStr = os.str();
      }
      
      dataCpus_ = parser.extractIndexRange(cpuStr);
    }
  }
}

/**.......................................................................
 * Parse a model specification
 */
void RunManager::getModelArgs(String& line, String& name, String& type, String& file, String& discard)
{
  //------------------------------------------------------------
  // Get to the first non-whitespace char
  //------------------------------------------------------------

  line.advanceToNextNonWhitespaceChar();

  //------------------------------------------------------------
  // Skip the addmodel token
  //------------------------------------------------------------

  line.findNextString();

  String nt;
  do {

    line.advanceToNextNonWhitespaceChar();
    nt = line.findNextString();

    if(!nt.isEmpty()) {
      String tok,val;
      getTokVal(nt, tok, val);

      if(tok.str() == "name") {
	name = val;
      } else if(tok.str() == "type") {
	type = val;
      } else if(tok.str() == "file") {
	file = val;
      } else if(tok.str() == "discard") {
	discard = val;
      } else {
	ThrowError("Unrecognized token: " << tok);
      }
    }
  } while(!nt.isEmpty());
}

void RunManager::getDataSetNameAndType(String& line, std::string& name, std::string& type)
{
  // Get to the first non-whitespace char

  line.advanceToNextNonWhitespaceChar();

  // Skip the addmodel token

  line.findNextString();

  String nt;
  do {

    nt = line.findNextString();

    if(!nt.isEmpty()) {
      String tok,val;
      getTokVal(nt, tok, val);

      if(tok.str() == "name") {
	name = val.str();
      } else if(tok.str() == "type") {
	type = val.str();
      } else {
	ThrowColorError(std::endl << "Unrecognized token: " << tok << std::endl << std::endl
			<< "Use like: 'adddataset name=nameStr type=typeStr'" << std::endl << std::endl
			<< "While parsing line: '" << line << "'", "red");
      }
    }
  } while(!nt.isEmpty());
}

/**.......................................................................
 * Parse an expression of the form addobs name=whatever
 */
void RunManager::getObsName(String& line, std::string& name)
{
  // Get to the first non-whitespace char

  line.advanceToNextNonWhitespaceChar();

  // Skip the addobs token

  line.findNextString();

  String nt;
  do {

    nt = line.findNextString();

    if(!nt.isEmpty()) {
      String tok,val;
      getTokVal(nt, tok, val);

      if(tok.str() == "name") {
	val.strip(' ');
	name = val.str();
      } else {
	ThrowError("Unrecognized token: " << tok);
      }
    }
  } while(!nt.isEmpty());
}

void RunManager::getOutputArgs(String& line)
{
  line.advanceToNextNonWhitespaceChar();

  // Skip the output token
  
  line.findNextString();
  
  String nt;
  do {

    nt = line.findNextString();

    if(!nt.isEmpty()) {
      String tok,val;
      getTokVal(nt, tok, val);

      if(tok.str() == "file") {
	mm_.setOutputFileName(val.str());
      } else {
	ThrowError("Unrecognized token: " << tok);
      }
    }

  } while(!nt.isEmpty());

}

/**.......................................................................
 * Parse an expression of the form:
 *
 *   gendata file=whatever dataset=ds sigma=n
 *
 * Used for generating fake data from a model
 */
void RunManager::getGenDataArgs(String& line)
{
  line.advanceToNextNonWhitespaceChar();

  // Skip the output token
  
  line.findNextString();
  
  bool hasFile = false;
  bool hasDataSet = false;

  String nt;
  do {

    nt = line.findNextString();

    if(!nt.isEmpty()) {
      String tok,val;
      getTokVal(nt, tok, val);
      
      tok.strip(' ');
      val.strip(' ');

      if(tok.str() == "file") {
	mm_.setGenDataFile(val.str());
	mm_.setGenData(true);
	hasFile = true;
      } else if(tok.str() == "dataset") {
	DataSet* dataSet = dm_.getDataSet(val.str());
	mm_.setGenDataDataSet(dataSet);
	hasDataSet = true;
      } else if(tok.str() == "sigma") {
	mm_.setGenDataSigma(val.toDouble());
      } else {
	ThrowError("Unrecognized token: '" << tok << "'");
      }
    }

    line.advanceToNextNonWhitespaceChar();

  } while(!nt.isEmpty());

  if(!hasFile)
    ThrowSimpleColorError("No output file specified, while parsing line: '" << line << "'" << std::endl << std::endl
			  << "Usage: gendata file=fileName dataset=datasetName sigma=sigmaVal", "red");

  if(!hasDataSet)
    ThrowSimpleColorError("No dataset specified, while parsing line: '" << line << "'" << std::endl << std::endl
			  << "Usage: gendata file=fileName dataset=datasetName sigma=sigmaVal", "red");
}

void RunManager::getLoadDataArgs(String& line)
{
  line.advanceToNextNonWhitespaceChar();

  // Skip the output token
  
  line.findNextString();
  
  String file, name;
  unsigned discard=0;
  
  String nt;
  do {

    nt = line.findNextString();

    if(!nt.isEmpty()) {
      String tok,val;
      getTokVal(nt, tok, val);
      
      tok.strip(' ');
      val.strip(' ');
      
      if(tok.str() == "file") {
	file = val;
      } else if(tok.str() == "name") {
	name = val;
      } else if(tok.str() == "discard") {
	discard = val.toInt();
      } else {
	ThrowError("Unrecognized token: '" << tok << "'");
      }
    }
  } while(!nt.isEmpty());

  //------------------------------------------------------------
  // Perform any initialization for loading an output file
  //------------------------------------------------------------

  mm_.initializeForOutput(runFile_);

  mm_.loadOutputFile(file.str(), name.str(), discard);

  runMarkov_      = false;
  loadOutputFile_ = true;
}

void RunManager::getComputeChisqArgs(String& line)
{
  String tok,val;
  getTokVal(line, tok, val);
      
  tok.strip(' ');
  val.strip(' ');
  
  if(tok.str() == "comp") {
    
    if(val.str() == "true") {
      compChisq_ = true;
    } else {
      compChisq_ = false;
    }
    
  } else {
    ThrowError("Unrecognized token: '" << tok << "'");
  }

}

void RunManager::getDisplayDataSetsArgs(String& line)
{
  String tok,val;
  getTokVal(line, tok, val);
      
  tok.strip(' ');
  val.strip(' ');
  
  if(tok.str() == "display") {
    
    if(val.str() == "true") {
      displayDataSets_ = true;
    } else {
      displayDataSets_ = false;
    }
    
  } else {
    ThrowError("Unrecognized token: '" << tok << "'");
  }

}

void RunManager::getWriteDataSetsArgs(String& line)
{
  String tok,val;
  getTokVal(line, tok, val);
      
  tok.strip(' ');
  val.strip(' ');
  
  if(tok.str() == "writedata") {
    
    if(val.str() == "true") {
      writeDataSets_ = true;
    } else {
      writeDataSets_ = false;
    }
    
  } else {
    ThrowError("Unrecognized token: '" << tok << "'");
  }

}

/**.......................................................................
 * Get the value of the next token
 */
void RunManager::getTokVal(String& line, String& tok, String& val)
{
  line.advanceToNextNonWhitespaceChar();

  if(line.contains("+="))
    tok = line.findNextInstanceOf("", false, "+=", true, true);
  else
    tok = line.findNextInstanceOf("", false, "=", true, true);

  if(tok.isEmpty())
    ThrowError("Empty token name");

  val = line.remainder();

  if(val.isEmpty())
    ThrowError("Empty value");

  // Expand any shell ~ directives before returning

  val.expandTilde();

  // And strip out any whitespace

  tok.strip(' ');
  val.strip(' ');
}

void RunManager::parseCorrelationVariables(String& token, 
					   std::string& varName1, std::string& varName2)
{
  token.strip(' ');

  String varStr1, varStr2;
  varStr1 = token.findNextInstanceOf("", false, "*", true, false);
  varStr2 = token.findNextInstanceOf("*", true, " ", false, true);

  varName1 = varStr1.str();
  varName2 = varStr2.str();
}

/**.......................................................................
 * Parse the name of a variable, potentially of the form: 
 * 
 *    model.var[aspect]
 */
void RunManager::parseVarname(String& mVarName, String& modelName, String& varName, String& aspectName)
{
  if(!mVarName.contains(".")) {
    ThrowError("Invalid variable name: " << mVarName << " (should be owner.varname)");
  }

  modelName = mVarName.findNextInstanceOf(" ", false, ".", true,  true);

  if(mVarName.remainder().contains("[")) {
    varName   = mVarName.findNextInstanceOf(" ", false,  "[", true, false);
    aspectName = mVarName.findNextInstanceOf("[", true,  "]", true, true);
  } else {
    varName   = mVarName.remainder();
  }

  modelName.strip(' ');
  varName.strip(' ');
  aspectName.strip(' ');
}

/**.......................................................................
 * Parse the name of a model or dataset variable, potentially of the form: 
 * 
 *   model.var[aspect]
 *
 * or
 *
 *   dataset.parameter
 *
 * or
 *
 *   obs.parameter
 *
 * or
 * 
 *   dataset.dataset.parameter
 *
 * or 
 *
 *   dataset.obs.parameter
 */
void RunManager::parseVarname(String& mVarName, 
			      DataSet** dataSet, Model** model, ObsInfo** obs, 
			      String& varName,   String& aspectName)
{
  if(!mVarName.contains(".")) 
    ThrowError("Invalid variable name: " << mVarName << " (should be owner.varname)");

  *dataSet = 0;
  *model   = 0;
  *obs     = 0;

  //------------------------------------------------------------
  // Now check if this is a model, dataset, or obs object that we know
  // about
  //------------------------------------------------------------

  bool isObj = true;
  String name;
  while(mVarName.remainder().contains(".") && isObj) {
    name  = mVarName.findNextInstanceOf(" ", false, ".", true,  true);
    isObj = getObject(name, model, dataSet, obs);
  }

  if(mVarName.remainder().contains("[")) {
    varName    = mVarName.findNextInstanceOf(" ", false,  "[", true, false);
    aspectName = mVarName.findNextInstanceOf("[", true,  "]", true, true);
  } else {
    varName    = mVarName.remainder();
  }

  //------------------------------------------------------------
  // If the last string we parsed wasn't a valid object, concatenate
  // it with the varname
  //------------------------------------------------------------

  if(!isObj) {
    std::ostringstream os;
    os << name << "." << varName;
    varName = os.str();
  }

  // Check that the last item we read wasn't also an object (ie, obs)

  isObj = getObject(varName, model, dataSet, obs);

  if(isObj) {
    varName = "";
  }

  varName.strip(' ');
  aspectName.strip(' ');
}

/**.......................................................................
 * Return the object (if any) with the requested name
 */
bool RunManager::getObject(String& objName, 
			   Model** model, DataSet** dataSet, ObsInfo** obs)
{
  objName.strip(' ');

  //------------------------------------------------------------
  // If we already have a dataset, this object might be another
  // dataset managed by it, or an obs
  //------------------------------------------------------------

  if(*dataSet != 0) {

    // Check for obs

    if(objName == "obs") {
      *obs = &((*dataSet)->getObs());
      return true;
    }

    // Check for nested dataset

    DataSetManager* dm = dynamic_cast<DataSetManager*>(*dataSet);
    if(dm) {
      try {
	*dataSet = dm->getDataSet(objName.str());
	return true;
      } catch(Exception& err) {
      }
    }

    return false;
  }

  //-------------------------------------------------------------
  // If we got here, then we are at the root of the hierarchy.  Check
  // for model, dataset or obs
  //-------------------------------------------------------------

  try {
    *model = mm_.getModel(objName.str());
    return true;
  } catch(Exception& err) {
  }
  
  try {
    *dataSet = dm_.getDataSet(objName.str());
    return true;
  } catch(Exception& err) {
  }

  try {
    *obs = om_.getObs(objName.str());
    return true;
  } catch(Exception& err) {
  }

  return false;
}

/**.......................................................................
 * Parse a variable assignment of the form 'var = expr;'
 */
void RunManager::parseVariableAssignment(String& str, String& tok, String& valStr)
{
  //------------------------------------------------------------
  // First parse the variable name
  //------------------------------------------------------------

  tok.strip(' ');

  String ownerName, varName, aspectName;
  parseVarname(tok, ownerName, varName, aspectName);

  tok.resetToBeginning();

  //------------------------------------------------------------
  // Now check if the owner is a model or a dataset, or neither
  //------------------------------------------------------------

  Model*   model   = 0;
  DataSet* dataSet = 0;
  ObsInfo* obs     = 0;

  try {
    model = mm_.getModel(ownerName.str());
  } catch(Exception& err) {
  }

  if(model) {
    try {
      return parseModelVariableAssignment(tok, valStr);

      //------------------------------------------------------------
      // Catching an exception here means we have a valid model, but
      // invalid variable or parameter name.  In this case we print
      // valid names and rethrow.
      //------------------------------------------------------------

    } catch(Exception& err) {

      if(!err.printHelp()) {
	ThrowSimpleError(err.what());
      } else {
	std::ostringstream os1, os2, os3;
	XtermManip xtm;

	model->general_.listParameters(os1);
	model->listComponents(os2);
	model->listParameters(os3);

	if(model->parameterMap_.size() > 0) {

	  if(os1.str().size() > 0) {
	    ThrowSimpleError(COLORIZE(xtm, "red", std::endl << "Unrecognized model component: '" << tok << "'" << std::endl)  
			     << COLORIZE(xtm, "green", std::endl << "Model " << model->name_ << " has description: ")
			     << std::endl << std::endl << os1.str()
			     << COLORIZE(xtm, "green", std::endl << "components: ")
			     << std::endl << std::endl << os2.str()
			     << COLORIZE(xtm, "green", std::endl << "and parameters: ")
			     << std::endl << std::endl << os3.str() << std::endl);
	  } else {
	    ThrowSimpleError(COLORIZE(xtm, "red", std::endl << "Unrecognized model component: '" << tok << "'" << std::endl)  
			     << COLORIZE(xtm, "green", std::endl << "Model " << model->name_ << " has components: ")
			     << std::endl << std::endl << os2.str()
			     << COLORIZE(xtm, "green", std::endl << "and parameters: ")
			     << std::endl << std::endl << os3.str() << std::endl);
	  }

	} else {

	  ThrowSimpleError(COLORIZE(xtm, "red", std::endl << "Unrecognized model component: '" << str << "'" << std::endl) << 
			   COLORIZE(xtm, "green", std::endl << "Model " << model->name_ << " has components: " 
				    << std::endl << std::endl << os2.str() << std::endl));

	}
      }
    }
  }

  //------------------------------------------------------------
  // If we got here, then no model matches -- check datasets
  //------------------------------------------------------------

  try {
    dataSet = dm_.getDataSet(ownerName.str());
  } catch(Exception& err) {
  }

  if(dataSet) {
    return parseDataSetVariableAssignment(dataSet, varName, valStr);
  }

  try {
    obs = om_.getObs(ownerName.str());
  } catch(Exception& err) {
  }

  if(obs)
    return parseObsVariableAssignment(obs, varName, valStr);

  ThrowColorError(std::endl << "No model or dataset matches the name: '" << ownerName << "'", "red");
}

/**.......................................................................
 * Parse a variable assignment of the form 'var = expr;'
 */
void RunManager::parseVariableAssignmentNew(String& str, String& tok, String& valStr)
{
  //------------------------------------------------------------
  // First parse the variable name
  //------------------------------------------------------------

  tok.strip(' ');

  Model*   model   = 0;
  DataSet* dataSet = 0;
  ObsInfo* obs     = 0;

  String varName, aspectName;

  parseVarname(tok, &dataSet, &model, &obs, varName, aspectName);

  tok.resetToBeginning();

  //------------------------------------------------------------
  // Now check if the owner is a model or a dataset, or neither
  //------------------------------------------------------------

  if(model) {
    try {
      return parseModelVariableAssignment(tok, valStr);

      //------------------------------------------------------------
      // Catching an exception here means we have a valid model, but
      // invalid variable or parameter name.  In this case we print
      // valid names and rethrow.
      //------------------------------------------------------------

    } catch(Exception& err) {

      if(!err.printHelp()) {
	ThrowSimpleError(err.what());
      } else {
	std::ostringstream os1, os2, os3;
	XtermManip xtm;

	model->general_.listParameters(os1);
	model->listComponents(os2);
	model->listParameters(os3);

	if(model->parameterMap_.size() > 0) {

	  if(os1.str().size() > 0) {
	    ThrowSimpleError(COLORIZE(xtm, "red", std::endl << "Unrecognized model component: '" << tok << "'" << std::endl)  
			     << COLORIZE(xtm, "green", std::endl << "Model " << model->name_ << " has description: ")
			     << std::endl << std::endl << os1.str()
			     << COLORIZE(xtm, "green", std::endl << "components: ")
			     << std::endl << std::endl << os2.str()
			     << COLORIZE(xtm, "green", std::endl << "and parameters: ")
			     << std::endl << std::endl << os3.str() << std::endl);
	  } else {
	    ThrowSimpleError(COLORIZE(xtm, "red", std::endl << "Unrecognized model component: '" << tok << "'" << std::endl)  
			     << COLORIZE(xtm, "green", std::endl << "Model " << model->name_ << " has components: ")
			     << std::endl << std::endl << os2.str()
			     << COLORIZE(xtm, "green", std::endl << "and parameters: ")
			     << std::endl << std::endl << os3.str() << std::endl);
	  }

	} else {

	  ThrowSimpleError(COLORIZE(xtm, "red", std::endl << "Unrecognized model component: '" << str << "'" << std::endl) << 
			   COLORIZE(xtm, "green", std::endl << "Model " << model->name_ << " has components: " 
				    << std::endl << std::endl << os2.str() << std::endl));

	}
      }
    }
  }

  //------------------------------------------------------------
  // If we got here, then no model matches -- check datasets
  //------------------------------------------------------------

  if(obs) {

    //------------------------------------------------------------
    // If we were passed an empty variable name, we are setting the
    // obs itself
    //------------------------------------------------------------

    if(varName.isEmpty()) {
      return dataSet->setObs(om_.getObs(valStr.str()));

      //------------------------------------------------------------
      // Else just set a parameter of this obs
      //------------------------------------------------------------

    } else {
      return parseObsVariableAssignment(obs, varName, valStr);
    }
  }

  if(dataSet) {
    return parseDataSetVariableAssignment(dataSet, varName, valStr);
  }

  {
    String ownerName, varName, aspectName;
    parseVarname(tok, ownerName, varName, aspectName);
    ThrowSimpleColorError(std::endl << "No model, dataset or obs matches the name: '" << ownerName << "'", "red");
  }
}

/**.......................................................................
 * Parse an assignment to obs parameters
 */
void RunManager::parseObsVariableAssignment(ObsInfo* obs, String& varName, String& valStr)
{
  // Next parse the expression

  valStr.strip(' ');

  // If the beginning of this string is numeric, then also check for units

  if(valStr.isNumeric(valStr[0])) {

    String numStr  = valStr.findNextNumericString(":,.");
    String unitStr = valStr.remainder();

    obs->setParameter(varName.str(), numStr.str(), unitStr.str());

  } else {
    obs->setParameter(varName.str(), valStr.str());
  }

}

/**.......................................................................
 * Parse an assignment to cosmo parameters
 */
void RunManager::parseCosmoVariableAssignment(Cosmology* cosmo, String& varName, String& valStr)
{
  // Next parse the expression

  valStr.strip(' ');

  // If this is a numeric string, then also check for units

  if(valStr.isNumeric()) {

    String numStr  = valStr.findNextNumericString(":,.");
    String unitStr = valStr.remainder();

    cosmo->setParameter(varName.str(), numStr.str(), unitStr.str());

  } else {
    cosmo->setParameter(varName.str(), valStr.str());
  }

}

void RunManager::parseDataSetObsVariableAssignment(DataSet* dataSet, String& varName, String& valStr)
{
  String ownerName, parName, aspectName;
  parseVarname(varName, ownerName, parName, aspectName);

  ObsInfo& obs = dataSet->getObs();
  parseObsVariableAssignment(&obs, parName, valStr);
}

void RunManager::parseModelCosmoVariableAssignment(Model* model, String& varName, String& valStr)
{
  String ownerName, parName, aspectName;
  parseVarname(varName, ownerName, parName, aspectName);

  Cosmology& cosmo = model->getCosmo();
  parseCosmoVariableAssignment(&cosmo, parName, valStr);
}

/**.......................................................................
 * Parse an assignment to a dataset parameters
 */
void RunManager::parseDataSetVariableAssignment(DataSet* dataSet, String& varName, String& inputValStr)
{
  //------------------------------------------------------------
  // Strip spaces from the value
  //------------------------------------------------------------

  inputValStr.strip(' ');

  String valStr = getVariableValueString(inputValStr);

  //------------------------------------------------------------
  // Check for special variables first.  These are assignment-only,
  // like dataset parameters, but we have to check for them
  // separately, since we want to allow things like ra = 12:34:34.2
  // which would otherwise get interpreted as a prior specification
  //------------------------------------------------------------

  if(varName == "ra" && valStr.contains(":")) {
    dataSet->setParameter(varName.str(), valStr.str());
    return;
  } else if(varName == "dec" && valStr.contains(":")) {
    dataSet->setParameter(varName.str(), valStr.str());
    return;
  } else if(varName == "region") {
    dataSet->setParameter(varName.str(), valStr.str());
    return;
  }

  //------------------------------------------------------------
  // If this is a numeric string, then also check for units.  But if
  // it appears to be a list, just pass it unparsed to the dataset
  //------------------------------------------------------------

  if(!valStr.contains(",") && valStr.isNumeric()) {

    String numStr  = valStr.findNextNumericString();
    String unitStr = valStr.remainder();

    dataSet->setParameter(varName.str(), numStr.str(), unitStr.str());

  } else {

    dataSet->setParameter(varName.str(), valStr.str());

  }
}

/**.......................................................................
 * Parse an assignment to a model variable
 */
void RunManager::parseModelVariableAssignment(String& tokStr, String& inputValStr)
{
  Variate* var = 0;
  ParameterManager::Parameter* par = 0;

  String ownerName, varName, aspectName;
  parseVarname(tokStr, ownerName, varName, aspectName);

  std::ostringstream os;
  os << ownerName << "." << varName;

  String tok = os.str();

  try {
    var = mm_.getVar(tok.str());
  } catch(Exception& err) {
  }

  try {
    par = mm_.getParameter(tok.str(), false);
  } catch(Exception& err) {
  }

  if(var == 0 && par == 0) {
    ThrowColorErrorHelp(tok << " matches no known model or parameter", "red");
  }

  String valStr = getVariableValueString(inputValStr);

  if(var) {

    //------------------------------------------------------------
    // Check for special aspects first
    //------------------------------------------------------------

    if(!aspectName.isEmpty()) {

      //------------------------------------------------------------
      // Allow the user to control display of variates
      //------------------------------------------------------------

      if(aspectName == "display") {

	valStr.strip(' ');

	if(valStr == "false") {
	  var->setDisplay(false);
	} else if(valStr == "true") {
	  var->setDisplay(true);
	} else {
	  ThrowColorError("Invalid display state: '" << valStr << "'.  Should be 'true' or 'false'", "red");
	}

	return;
      }

      //------------------------------------------------------------
      // Allow the user to add derived variates
      //------------------------------------------------------------

      if(aspectName == "derive") {

	if(!var->isDerivable()) {
	  ThrowColorError(var->name_ << " isn't a derivable variate", "red");
	}

	valStr.strip(' ');

	if(valStr == "false") {
	  var->wasSpecified_ = false;
	} else if(valStr == "true") {
	  var->wasSpecified_ = true;
	  var->isDerived_    = true;

	  var->wasRequested_ = true;

	  // Default to displaying derived variates.  This will be
	  // used to distinguish internally between display of derived
	  // variates requested by the user and display of malleable
	  // derived variates required internally

	  var->setDisplay(true);
	  var->isVisible_ = true;

	} else {
	  ThrowColorError("Invalid derive state: '" << valStr << "'.  Should be 'true' or 'false'", "red");
	}

	return;
      }

      //------------------------------------------------------------
      // Allow the user to change the order in which this variate is
      // displayed
      //------------------------------------------------------------

      if(aspectName == "order") {
	var->setDisplayOrder(valStr.toUint());
      }

      //------------------------------------------------------------
      // Allow the user to change the display range.  We allow any
      // units to be used, even if the units are currently different
      // from the units in which the range is specified
      // ------------------------------------------------------------

      if(aspectName == "range") {
	String meanStr, minStr, maxStr, units;
	parseUniformPrior(valStr, meanStr, minStr, maxStr, units);

	double displayMin = var->getVal(minStr.toDouble(), units.str());
	double displayMax = var->getVal(maxStr.toDouble(), units.str());

	var->setDisplayRange(var->getUnitVal(displayMin), var->getUnitVal(displayMax));
      }

      if(aspectName == "prior") {
	String meanStr, minStr, maxStr, units;
	parseUniformPrior(valStr, meanStr, minStr, maxStr, units);

	double min = var->getVal(minStr.toDouble(), units.str());
	double max = var->getVal(maxStr.toDouble(), units.str());

	// Reverse min/max if specification was reversed.

	if(min > max) {
	  double tmp = min;
	  min = max;
	  max = tmp;
	}
	
	// If a mean was specified, set it explicitly, else, infer it
	// from the min/max
	
	double mean=0.0;
	if(!meanStr.isEmpty()) {
	  mean = var->getVal(meanStr.toDouble(), units.str());
	} else {
	  mean = (max + min)/2;
	}
	
	var->prior().setType(Distribution::DIST_UNIFORM);
	var->prior().setUniformXMin(min);
	var->prior().setUniformXMax(max);

	var->wasSpecified_ = true;

	// And set isDerived always so that malleable variates are
	// automatically set correctly
	
	var->isDerived_    = false;
	
	try {
	  var->setIsVariable(true);
	} catch(...) {
	  ThrowColorError("You must set a fixed value for '" << tok.str() << "'", "red");
	}
	
	var->samplingDistribution().setType(Distribution::DIST_GAUSS);
	var->samplingDistribution().setGaussMean(mean);
	var->samplingDistribution().setGaussSigma(0.1*(max - min));
      }

      //------------------------------------------------------------
      // Allow the user to change the units
      //------------------------------------------------------------

      if(aspectName == "units") {

	//------------------------------------------------------------
	// If we have already set a display range for this variate,
	// get the range first in native units, and restore it in the
	// new units
	//------------------------------------------------------------

	double displayMin;
	double displayMax;

	if(var->hasRange_) {
	  displayMin = var->getVal(var->displayMin_, var->units());
	  displayMax = var->getVal(var->displayMax_, var->units());
	}
	
	//------------------------------------------------------------
	// Now set the new units
	//------------------------------------------------------------

	valStr.strip(' ');
	var->setUnits(valStr.str());

	//------------------------------------------------------------
	// Now restore the range
	//------------------------------------------------------------

	if(var->hasRange_) 
	  var->setDisplayRange(var->getUnitVal(displayMin), var->getUnitVal(displayMax));
      }

      //------------------------------------------------------------
      // Any mention of a derived variate by name probably means that
      // the user wants to derive it
      //------------------------------------------------------------

      if(var->isDerived_) {
	var->wasSpecified_ = true;
      }

      return;
    }

    //------------------------------------------------------------
    // Check for special variables first.  These are assignment-only,
    // like dataset parameters, but we have to check for them
    // separately, since we want to allow things like ra = 12:34:34.2
    // which would otherwise get interpreted as a prior specification
    //------------------------------------------------------------

    String tokVal = tok.findNextInstanceOf(".", true, " ", false, true);

    if(tokVal == "ra" && valStr.contains(":")) {
      var->setVal(valStr.str());
      var->wasSpecified_ = true;
      return;
    } else if(tokVal == "dec" && valStr.contains(":")) {
      var->setVal(valStr.str());
      var->wasSpecified_ = true;
      return;
    } else if(tokVal == "region") {
      var->setVal(valStr.str());
      var->wasSpecified_ = true;
      return;
    }

    //------------------------------------------------------------
    // If the expression is of the form 'var = min:max;' or 
    // 'var = mean:min:max;', then we assume
    // a uniform prior is being specified
    //------------------------------------------------------------

    if(valStr.contains(":") && !valStr.contains("+-")) {

      if(!var->canBePrimary())
	ThrowColorError("Variate " << var->name_ << " is a derived variate only", "red");

      String startStr, minStr, maxStr, units;
      double start, min, max;

      parseUniformPrior(valStr, startStr, minStr, maxStr, units);

      var->setUnits(units.str());

      min = var->getVal(minStr.toDouble(), units.str());
      max = var->getVal(maxStr.toDouble(), units.str());

      // Reverse min/max if specification was reversed.

      if(min > max) {
	double tmp = min;
	min = max;
	max = tmp;
      }

      // If a start was specified, set it explicitly, else, infer it
      // from the min/max

      if(!startStr.isEmpty()) {
	start = var->getVal(startStr.toDouble(), units.str());

	if(start < min || start > max) {
	  ThrowSimpleColorError("Starting point: " << startStr 
				<< " lies outside the prior range " 
				<< minStr << ":" << maxStr, "red");
	}

      } else {
	start = (max + min)/2;
      }

      var->prior().setType(Distribution::DIST_UNIFORM);
      var->prior().setUniformXMin(min);
      var->prior().setUniformXMax(max);

      var->wasSpecified_ = true;

      // And set isDerived always so that malleable variates are
      // automatically set correctly

      var->isDerived_    = false;

      try {
	var->setIsVariable(true);
      } catch(...) {
	ThrowColorError("You must set a fixed value for '" << tok.str() << "'", "red");
      }

      var->samplingDistribution().setType(Distribution::DIST_GAUSS);
      var->samplingDistribution().setGaussMean(start);
      var->samplingDistribution().setGaussSigma(0.1*(max - min));

      //------------------------------------------------------------
      // Else if of the form 'var = mean +- sigma;' we assume a gaussian
      // prior is being specified
      //------------------------------------------------------------

    } else if(valStr.contains("+-")) {

      try {

	if(!var->canBePrimary())
	  ThrowColorError("Variate " << var->name_ << " is a derived variate only", "red");

	String startStr, meanStr, sigmaStr, units, minStr, maxStr;
	double start, mean, sigma;

	parseGaussianPrior(valStr, startStr, meanStr, sigmaStr, units, minStr, maxStr);

	var->setUnits(units.str());

	mean  = var->getVal(meanStr.toDouble(),  units.str());
	sigma = var->getVal(sigmaStr.toDouble(), units.str());

	if(minStr.isEmpty() && maxStr.isEmpty()) {
	  var->prior().setType(Distribution::DIST_GAUSS);
	  var->prior().setGaussMean(mean);
	  var->prior().setGaussSigma(sigma);
	} else {

	  var->prior().setType(Distribution::DIST_TRUNC_GAUSS);
	  var->prior().setGaussMean(mean);
	  var->prior().setGaussSigma(sigma);

	  double min = var->getVal(minStr.toDouble(),  units.str());
	  double max = var->getVal(maxStr.toDouble(),  units.str());

	  // Reverse min/max if specification was reversed.

	  if(min > max) {
	    double tmp = min;
	    min = max;
	    max = tmp;
	  }

	  var->prior().setUniformXMin(min);
	  var->prior().setUniformXMax(max);
	}

	var->wasSpecified_ = true;

	// And set isDerived always so that malleable variates are
	// automatically set correctly

	var->isDerived_    = false;

	try {
	  var->setIsVariable(true);
	} catch(...) {
	  ThrowColorError("You must set a fixed value for '" << tok.str() << "'", "red");
	}

	if(!startStr.isEmpty())
	  start = var->getVal(startStr.toDouble(),  units.str());
	else
	  start = mean;

	var->samplingDistribution().setType(Distribution::DIST_GAUSS);
	var->samplingDistribution().setGaussMean(start);
	var->samplingDistribution().setGaussSigma(sigma);

      } catch(Exception& err) {
	COUT("Caught an error: " << err.what());
	throw err;
      }

      //------------------------------------------------------------
      // Else this is just a simple 'var = val;' statement
      //------------------------------------------------------------

    } else {

      if(!var->canBePrimary())
	ThrowColorError("Variate " << var->name_ << " is a derived variate only", "red");

      String valueStr, units;
      String testStr = valStr;
      testStr.strip(' ');
      String testValStr  = testStr.findNextNumericString();
      String testUnitStr = testStr.remainder();

      //------------------------------------------------------------
      // If the value we are assigning is a symbolic name, this means
      // that we are being told to derive this variate from another.  
      //------------------------------------------------------------

      if(testValStr.isEmpty() && testStr.contains(".")) {
	Variate* srcVar  = mm_.getVar(testStr.str());
	var->deriveFrom(*srcVar);
	return;
      }

      //------------------------------------------------------------
      // Else we are assigning a fixed value to this component
      //------------------------------------------------------------

      if(!testValStr.isEmpty() && !testUnitStr.isEmpty()) {
	parseVal(valStr, valueStr, units);
	units.strip(' ');

	var->setUnits(units.str());
	var->setVal(valueStr.toDouble(), units.str());

      } else {
	valStr.advanceToNextNonWhitespaceChar();
	valStr.strip(' ');
	var->setVal(valStr.str());
      }

      var->wasSpecified_ = true;
      
      // And set isDerived always so that malleable variates are
      // automatically set correctly
      
      var->isDerived_    = false;

      // Likewise set isVariable to false so that malleable variates
      // are correctly designated

      var->setIsVariable(false);
    }

    //------------------------------------------------------------
    // Else this is a parameter assignment
    //------------------------------------------------------------

  } else {

    String valueStr, units;
    String testStr = valStr;

    testStr.strip(' ');

    if(varName.contains("cosmo.")) {
      Model* model = mm_.getModel(ownerName.str());
      parseModelCosmoVariableAssignment(model, varName, valStr);
      return;
    }

    if(!testStr.contains(",") && testStr.isNumeric()) {
      parseVal(valStr, valueStr, units);
      mm_.setParameter(tok.str(), valueStr.str(), units.str());
    } else {
      valStr.advanceToNextNonWhitespaceChar();
      valStr.strip(' ');
      mm_.setParameter(tok.str(), valStr.str());
    }
  }
}

/**.......................................................................
 * Parse a variable assignment of the form 'var += expr;'
 */
void RunManager::parseVariableIncrement(String& tok, String& valStr)
{
  // First parse the variable name

  tok.strip(' ');

  String ownerName, varName, aspectName;
  parseVarname(tok, ownerName, varName, aspectName);

  //------------------------------------------------------------
  // Now check if the owner is a model or a dataset, or neither
  //------------------------------------------------------------

  Model*   model   = 0;
  DataSet* dataSet = 0;

  try {
    model = mm_.getModel(ownerName.str());
  } catch(Exception& err) {
  }

  if(model) {
    ThrowColorError("Component increment is not currently supported for models", "red");
  }

  //------------------------------------------------------------
  // If we got here, then no model matches -- check datasets
  //------------------------------------------------------------

  try {
    dataSet = dm_.getDataSet(ownerName.str());
  } catch(Exception& err) {
  }

  if(dataSet)
    return parseDataSetVariableIncrement(dataSet, varName, valStr);

  ThrowColorError("No dataset matches the name: '" << ownerName << "'", "red");
}

/**.......................................................................
 * Parse an increment to a dataset parameter
 */
void RunManager::parseDataSetVariableIncrement(DataSet* dataSet, String& varName, String& valStr)
{
  // Next parse the expression

  valStr.strip(' ');

  //------------------------------------------------------------
  // Check for special variables first.  
  //------------------------------------------------------------

  if(varName == "ra" && valStr.contains(":")) {
    dataSet->incrementParameter(varName.str(), valStr.str());
    return;
  } else if(varName == "dec" && valStr.contains(":")) {
    dataSet->incrementParameter(varName.str(), valStr.str());
    return;
  }

  //------------------------------------------------------------
  // If this is a numeric string, then also check for units
  //------------------------------------------------------------

  if(valStr.isNumeric()) {

    String numStr  = valStr.findNextNumericString();
    String unitStr = valStr.remainder();

    dataSet->incrementParameter(varName.str(), numStr.str(), unitStr.str());

  } else {
    dataSet->incrementParameter(varName.str(), valStr.str());
  }

}

/**.......................................................................
 * Parse a statement of the form 'var = min:max;' specifying a uniform
 * prior
 */
void RunManager::parseUniformPrior(String& val, String& mean, String& min, String& max, String& units)
{
  val.advanceToNextNonWhitespaceChar();
  String tok1 = val.findNextInstanceOf(" ", false, ":", true, true);

  // If the string still contains another ':', then we assume we
  // are parsing a 'mean:min:max' specification, else just 'min:max'

  if(val.remainder().contains(":")) {
    mean = tok1;
    val.advanceToNextNonWhitespaceChar();
    min = val.findNextInstanceOf(" ", false, ":", true, true);
  } else {
    min = tok1;
  }

  val.advanceToNextNonWhitespaceChar();
  max   = val.findNextNumericString();
  units = val.remainder();
  units.strip(' ');
}

/**.......................................................................
 * Parse a statement of the form 'var = mean +- sigma;' specifying a
 * Gaussian prior
 */
void RunManager::parseGaussianPrior(String& val, String& start, String& mean, String& sigma, String& units, String& minStr, String& maxStr)
{
  val.advanceToNextNonWhitespaceChar();

  //------------------------------------------------------------
  // Check if the start of the string contains the ':' character.  If
  // so, assume we are specifying a starting position
  //------------------------------------------------------------

  String copy = val;
  String copyTest = val.findNextInstanceOf(" ", false, "+-", true, true);
  
  if(copyTest.contains(":")) {
    start = copyTest.findNextInstanceOf(" ", false, ":", true, true);
  }

  //------------------------------------------------------------
  // Now read the mean +- sigma part of the prior specification
  //------------------------------------------------------------

  copyTest.advanceToNextNonWhitespaceChar();
  mean  = copyTest.remainder();

  val.advanceToNextNonWhitespaceChar();
  sigma = val.findNextNumericString();

  //------------------------------------------------------------
  // If there is an 'x' character in the 
  //------------------------------------------------------------

  std::ostringstream os;
  if(val.remainder().contains("x")) {
    if(val.remainder().contains(":")) {
      minStr = val.findNextInstanceOf("x", true, ":", true, true);
      minStr.strip(' ');
      val.advanceToNextNonWhitespaceChar();

      //------------------------------------------------------------
      // Have to be careful here -- max can be a non-numeric string
      // (if inf was specified, say), so we explicitly check for inf
      // first
      //------------------------------------------------------------

      String infTest = val.toLower();
      if(infTest.contains("inf")) {
	maxStr = val.findNextInstanceOf(" ", false, "inf", true, true);
	os << maxStr.str() << "inf";
	maxStr = os.str();
      } else {
	maxStr = val.findNextNumericString();
      }

    } else {
      ThrowSimpleColorError("Invalid prior specification: " << val, "red");
    }
  }

  units = val.remainder();
  units.strip(' ');
}
 
/**.......................................................................
 * Parse a simple assignment of the form 'var = val;'
 */
void RunManager::parseVal(String& val, String& value, String& units)
{
  val.advanceToNextNonWhitespaceChar();

  value = val.findNextNumericString();
  units = val.remainder();
  units.strip(' ');
}

void RunManager::checkIfNameAlreadyExists(std::string name)
{
  // Check that no model or dataset already exists with that name

  bool exists=false;
  try {
    mm_.getModel(name);
    exists = true;
  } catch(...) {
  }

  if(exists) {
    ThrowColorError("A model with name: " << name << " has already been defined", "red");
  }

  try {
    dm_.getDataSet(name);
    exists = true;
  } catch(...) {
  }

  if(exists) {
    ThrowColorError("A dataset with name: " << name << " has already been defined", "red");
  }
}

/**.......................................................................
 * Run a Markov chain with the currently defined models and datasets
 */
void RunManager::runMarkov()
{
  //------------------------------------------------------------
  // Sanity check that nBurn < nTry
  //------------------------------------------------------------

  if(nTry_ <= nBurn_)
    ThrowSimpleColorError("You have specified ntry = " << nTry_ 
			  << " and nburn = " << nBurn_ 
			  << ".  ntry should be > nburn", "red");
  
  Probability likeCurr, likePrev, propDensCurr, propDensPrev;
  ChisqVariate chisq;
  bool converged = false;

  overallTimer_.start();

  //------------------------------------------------------------
  // Initialize variables we need for running Markov chains
  //------------------------------------------------------------

  initializeMarkovSpecificVariables();

  //------------------------------------------------------------
  // Main loop -- perform nTry_ iterations of the MH algorithm
  //------------------------------------------------------------

  COUT("");
  unsigned nTry=0;
  for(unsigned i=0; i < nTry_ && !converged; i++, nTry++) {
    
    //------------------------------------------------------------
    // If it's time to print our progress, do it now
    //------------------------------------------------------------

    printProgress(i);

    //------------------------------------------------------------
    // Generate a new sample
    //------------------------------------------------------------

    generateNewSample();

    //------------------------------------------------------------
    // See if this sample should be accepted
    //------------------------------------------------------------

    if(acceptMetropolisHastings(i, propDensCurr, propDensPrev, likeCurr, likePrev, chisq)) {

      ++nAcceptedSinceLastUpdate_;
      likePrev     = likeCurr;
      propDensPrev = propDensCurr;

      //------------------------------------------------------------
      // If this sample was accepted, and it's time to tune the
      // jumping distribution, do it now
      //------------------------------------------------------------

      if(timeToTune(i))
	tuneJumpingDistribution(i, propDensCurr, likeCurr);

      //------------------------------------------------------------
      // Store the latest accepted sample if we are not still in the
      // burn-in sequence.  We call storeMultiplicity() first so that
      // it gets installed for the previous accepted sample (not the
      // current one, for which we don't yet know the multiplicity)
      //------------------------------------------------------------

      if(i >= nBurn_ || incBurnIn_) {
	mm_.storeMultiplicity(nTimesAtThisPoint_);
	mm_.store(likeCurr, chisq);
      }

      nTimesAtThisPoint_ = 1;

    } else {
      ++nTimesAtThisPoint_;
      mm_.revert();
    }

    //------------------------------------------------------------
    // Check if the chain has converged
    //------------------------------------------------------------
    
    converged = checkConvergence(i);
  }

  //------------------------------------------------------------
  // Write the multiplicity of the last accepted sample (if any)
  //------------------------------------------------------------

  mm_.storeMultiplicity(nTimesAtThisPoint_);

  overallTimer_.stop();

  unsigned nTotal = incBurnIn_ ? nTry : (nTry - nBurn_);

  COUTCOLOR(std::endl, "yellow");
  COUTCOLOR("Elapsed time:                       "     << setprecision(1) << std::fixed << overallTimer_.deltaInSeconds() << "s", "yellow");
  COUTCOLOR("Time spent sampling:                  "   << std::setw(5) << std::right << setprecision(1) << std::fixed << sampleTime_                    << "s", "yellow");
  COUTCOLOR("Time spent tuning:                    "   << std::setw(5) << std::right << setprecision(1) << std::fixed << tuneTime_                      << "s", "yellow");
  COUTCOLOR("Time spent calculating likelihoods:   "   << std::setw(5) << std::right << setprecision(1) << std::fixed << likeTime_                      << "s", "yellow");
  COUTCOLOR("Time spent adding models:               " << std::setw(5) << std::right << setprecision(1) << std::fixed << dm_.addModelTime_              << "s", "yellow");
  COUTCOLOR("Time spent computing chisq:             " << std::setw(5) << std::right << setprecision(1) << std::fixed << dm_.computeChisqTime_          << "s", "yellow");
  COUTCOLOR(std::endl << "Fraction accepted:                  " <<(double)(mm_.nAccepted_)/(nTotal) << std::endl, "yellow");

#if 1
  dm_.debugPrint();
  mm_.debugPrint();
#endif
}

/**.......................................................................
 * Create the histogram + residual display
 */
void RunManager::createMarkovDisplay()
{
  PgUtil::open(pgplotDev_);
  PgUtil::setPrompt(false);
  PgUtil::setOverplot(true);
  PgUtil::setColormap("heat");
  PgUtil::setReverseX(false);

  //------------------------------------------------------------
  // Display the variate histograms
  //------------------------------------------------------------

  PgUtil::clearPgManager();
  mm_.histogramVariates(nBin_, true, stat_, nSigma_, varPlot_);

  //------------------------------------------------------------
  // Now add the best-fit model and display the datasets with
  // residuals
  //------------------------------------------------------------

  try {
    PgUtil::setInteractive(interactive_);
    displayDataSets();
  } catch(Exception& err) {
    XtermManip xtm;
    COUT(COLORIZE(xtm, "red", "Error displaying datasets: ") << err.what() <<
	 COLORIZE(xtm, "red", std::endl << "Attempting to continue" << std::endl));
  }

  PgUtil::close();

  //------------------------------------------------------------
  // Print the model
  //------------------------------------------------------------

  unsigned nVar = 0;
  mm_.printModel("Maximum Likelihood", nBin_, stat_, nSigma_, printConvergence_, targetVariance_, nVar);

  //------------------------------------------------------------
  // Finally, print the chisq and evidence
  //------------------------------------------------------------

  if(!loadOutputFile_) {
    ChisqVariate chisq = dm_.computeChisq();
    chisq.setChisq(chisq.chisq(), chisq.nDof() - nVar);
    COUTCOLOR(std::endl << "Maximum Likelihood model has " << chisq, "green");
  }

  COUTCOLOR(std::endl << "Ln Evidence: " << mm_.estimateLnEvidence() << std::endl, "green");
}

/**.......................................................................
 * Create the histogram only display
 */
void RunManager::createOutputDisplay()
{
  PgUtil::open(pgplotDev_);
  PgUtil::setOverplot(true);
  PgUtil::setColormap("heat");
  PgUtil::setReverseX(false);

  //------------------------------------------------------------
  // Display the variate histograms
  //------------------------------------------------------------

  mm_.histogramVariates(nBin_, false, stat_, nSigma_);
}

/**.......................................................................
 * Set the model to the best-fit values, and add it to all datasets
 */
void RunManager::assertBestFitModel(bool finalize)
{
  mm_.setValues(mm_.bestFitSample_, true);
  dm_.clearModel();
  dm_.addModel(mm_);

  if(finalize) {

    // Perform any final actions needed for display once the models have
    // been added
    
    dm_.finalizeForDisplay();
  }
}

void RunManager::assertCurrentModel(std::string reason, bool finalize)
{
  mm_.initializeForModelAssertion(reason);
  dm_.clearModel();
  dm_.addModel(mm_);

  if(finalize) {

    // Perform any final actions needed for display once the models have
    // been added
    
    dm_.finalizeForDisplay();
  }
}

/**.......................................................................
 * Display data, model and residuals
 */
void RunManager::displayDataSets()
{
  //------------------------------------------------------------
  // Calculate where the dataset plot should be located
  //------------------------------------------------------------

  unsigned nVar = mm_.nVar();

  if(nVar <= 1)
    nVar = 2;

  double xstart = 0.1;
  double ystart = 0.1;
  double xsep   = 0.01;
  double ysep   = 0.01;

  double xmin, xmax, ymin, ymax;

  xmax = 1.0;
  ymax = 0.98;

  double dx = (xmax - xstart - nVar * xsep)/nVar;
  double dy = (ymax - ystart - nVar * ysep)/nVar;

  if(nVar%2 == 0) {
    xmin = xstart + nVar/2 * (dx + xsep) + 3*xsep;
    //    xmax = xstart + nVar   * (dx + xsep) - xsep;
    ymin = ystart + nVar/2 * (dy + ysep) + ysep + 3*ysep;
  } else {
    xmin = xstart + (nVar/2 + 1) * (dx + xsep);
    ymin = ystart + (nVar/2 + 1) * (dy + ysep) + ysep;
  }
  
  //  xmax = xstart + nVar   * (dx + xsep) - xsep;
  //  ymax = 1.0;

  //------------------------------------------------------------
  // Now assert the best-fit model and populate the plot with datasets
  //------------------------------------------------------------

  dm_.initializeForDisplay();
  assertBestFitModel();

  populatePlotWithDataSets(xmin, xmax, ymin, ymax);
}

/**.......................................................................
 * Populate a triangle plot with data/model/res plots for all datasets
 * in the upper right quadrant
 */
void RunManager::populatePlotWithDataSets(double gxmin, double gxmax, double gymin, double gymax)
{
  PgUtil::setVp(false);

  cpgsvp(gxmin, gxmax, gymin, gymax);

  std::vector<unsigned> ixVec;
  std::vector<unsigned> iyVec;
  unsigned nside;

  getPlotIndices(nside, ixVec, iyVec);

  if((gymax - gymin) < 0.5)
    //    PgUtil::setCharacterHeight((0.5/nside) * (gymax-gymin)/0.25);
    PgUtil::setCharacterHeight((0.5/nside) * (gymax-gymin)/0.4);
  else
    PgUtil::setCharacterHeight(0.5);

  double dx = (gxmax - gxmin)/nside;
  double dy = (gymax - gymin)/nside;

  unsigned iPlot = 0;
  for(std::map<std::string, DataSet*>::iterator iter=dm_.dataSetMap_.begin(); iter != dm_.dataSetMap_.end(); iter++, iPlot++) {

    double xmin = gxmin + ixVec[iPlot] * dx;
    double xmax = xmin + dx;
    double ymax = gymax - iyVec[iPlot] * dy;
    double ymin = ymax - dy;
    
    cpgsvp(xmin, xmax, ymin, ymax);

    DataSet* ds = iter->second;

    if(ds->dataSetType_ & DataSetType::DATASET_1D) {
      cpgbox("BC", 0, 0, "BC", 0,0);
      displayDataSet1D(ds, xmin, xmax, ymin, ymax);
    } else {
      cpgwnad(0,1,0,1);
      float x1, x2, y1, y2;
      cpgqvp(0, &x1, &x2, &y1, &y2);
      displayDataSet2D(ds, x1, x2, y1, y2);
    }
  }
}

void RunManager::getPlotIndices(unsigned& nside, std::vector<unsigned>& ixVec, std::vector<unsigned>& iyVec)
{
  unsigned nDataSet = dm_.dataSetMap_.size();
  nside = (int)ceil(sqrt((double)nDataSet));

  if(nDataSet == 2) {
    ixVec.push_back(0);ixVec.push_back(1);
    iyVec.push_back(0);iyVec.push_back(1);
  } else if(nDataSet == 3) {
    ixVec.push_back(0);ixVec.push_back(0);ixVec.push_back(1);
    iyVec.push_back(0);iyVec.push_back(1);iyVec.push_back(1);
  } else if(nDataSet == 5) {
    ixVec.push_back(0);ixVec.push_back(0);ixVec.push_back(1);ixVec.push_back(2);ixVec.push_back(2);
    iyVec.push_back(0);iyVec.push_back(2);iyVec.push_back(1);iyVec.push_back(0);iyVec.push_back(2);
  } else {
    for(unsigned i=0; i < nDataSet; i++) {
      ixVec.push_back(i % nside);
      iyVec.push_back(i / nside);
    }
  }
}

/**.......................................................................
 * Display information for a 1D dataset
 */
void RunManager::displayDataSet1D(DataSet* ds, double xvp1, double xvp2, double yvp1, double yvp2)
{
  //------------------------------------------------------------
  // We will display the residuals as 20% the height of the data plot
  //------------------------------------------------------------

  double dyvp = (yvp2 - yvp1)/5;

  PgUtil::setXTick(true);
  PgUtil::setYTick(true);
  PgUtil::setXTickLabeling(false);
  PgUtil::setYTickLabeling(true);

  PgUtil::setWin(true);
  PgUtil::setBox(true);

  PgUtil::setXLabel(false);
  PgUtil::setYLabel(true);
  PgUtil::setLogPlot(false);

  //------------------------------------------------------------
  // Display the dataset
  //------------------------------------------------------------

  cpgsvp(xvp1, xvp2, yvp1+dyvp, yvp2);
  ds->display();

  //------------------------------------------------------------
  // Display the residuals
  //------------------------------------------------------------

  float xw1, xw2, yw1, yw2;
  cpgqwin(&xw1, &xw2, &yw1, &yw2);

  PgUtil::setTraceColor(6);
  PgUtil::setXTickLabeling(true);
  PgUtil::setXLabel(true);

  cpgsvp(xvp1, xvp2, yvp1, yvp1+dyvp);
  ds->displayResiduals();

  //------------------------------------------------------------
  // Finally, display the composite model
  //------------------------------------------------------------

  PgUtil::setWin(false);
  PgUtil::setBox(false);

  PgUtil::clearTraceColor();
  PgUtil::clearBoxColor();

  // For 1D datasets, resample the x-axis to the finest sampling found
  // in the data

  DataSet1D* ds1d = dynamic_cast<DataSet1D*>(ds);
  double xmin  = ds1d->x_[0];
  double xmax  = ds1d->x_[0];
  double dxmin = ds1d->x_[1] - ds1d->x_[0];
  double xcurr, xprev, dx;

  for(unsigned i=0; i < ds1d->x_.size(); i++) {

    xcurr = ds1d->x_[i];

    if(i > 0) {
      dx = xcurr - xprev;
      dxmin = dxmin < dx ? dxmin : dx;
    }

    xmin = xcurr < xmin ? xcurr : xmin;
    xmax = xcurr > xmax ? xcurr : xmax;

    xprev = xcurr;
  }

  dxmin = 0.1;

  unsigned n = (unsigned)ceil(((xmax - xmin) / dxmin));

  std::vector<double> xSave = ds1d->x_;
  std::vector<double> ySave = ds1d->compositeModel_;

  ds1d->x_.resize(n);
  ds1d->compositeModel_.resize(n);

  for(unsigned i=0; i < n; i++) {
    xcurr = xmin + dxmin * i;
    ds1d->x_[i] = xcurr;
  }

  mm_.setValues(mm_.bestFitSample_, true);
  dm_.clearModel();
  dm_.addModel(mm_);

  PgUtil::setTraceColor(6);
  cpgsvp(xvp1, xvp2, yvp1+dyvp, yvp2);
  cpgswin(xw1, xw2, yw1, yw2);

  ds->displayCompositeModel();
  PgUtil::clearTraceColor();

  // Now reset the internal data arrays

  ds1d->x_ = xSave;
  ds1d->compositeModel_ = ySave;
}

/**.......................................................................
 * Display function for 2D data sets
 */
void RunManager::displayDataSet2D(DataSet* ds, double xvp1, double xvp2, double yvp1, double yvp2)
{
  double dxvp = (xvp2 - xvp1)/(2.3);
  double dyvp = (yvp2 - yvp1)/(2.3);

  PgUtil::setWin(true);
  PgUtil::setBox(true);

  PgUtil::setTick(true);
  PgUtil::setLabel(true);
  PgUtil::setXTick(true);
  PgUtil::setYTick(true);
  PgUtil::setXTickLabeling(true);
  PgUtil::setYTickLabeling(true);
  
  PgUtil::setXLabel(false);
  PgUtil::setYLabel(true);
  PgUtil::setWedge(true);

  ds->insertDisplayModels();

  cpgsvp(xvp1, xvp1+dxvp, yvp2-dyvp, yvp2);
  ds->display();

  PgUtil::setXTickLabeling(true);
  PgUtil::setXLabel(true);
  PgUtil::setWedge(false);

  cpgsvp(xvp1, xvp1+dxvp, yvp1+0.1*dyvp, yvp1+dyvp+0.1*dyvp);
  ds->displayCompositeModel();

  PgUtil::setYTickLabeling(false);
  PgUtil::setYLabel(false);
  PgUtil::setWedge(false);

  cpgsvp(xvp2-dxvp, xvp2, yvp1+0.1*dyvp, yvp1+dyvp+0.1*dyvp);
  ds->displayResiduals();

  // Reset any legacy display settings

  PgUtil::setZmin(0.0);
  PgUtil::setZmax(0.0);
  PgUtil::clearTraceColor();
  PgUtil::clearBoxColor();

  ds->clearDisplayModels();
}

/**.......................................................................
 * Method to update the jumping distribution based on a quadratic
 * estimator for the posterior
 */
void RunManager::updateHessian(Probability startProb)
{
  // Store the current vector of values

  Vector<double> storedVal = mm_.currentSample_;
  Vector<double> mean      = mm_.mean_;
  Vector<double> sigmas    = mm_.sigma_;
  Matrix<double> invCov    = mm_.invCov_;
  double detC              = mm_.detC_;

  unsigned nVar = storedVal.size();
  Vector<double> x0(nVar), x1(nVar), x2(nVar), tmp(nVar);
  double lp0, lp1, lp2, dlp1, dlp2, d2;
  double val1, val2;

  Vector<double> newSigmas = sigmas;
  double newsigval;
  for(unsigned iVar=0; iVar < nVar; iVar++) {

    // Generate a new random point from this variate's conditional
    // distribution

    tmp = storedVal;
    Probability prob = startProb;

    x0 = tmp;
    lp0 = prob.lnValue();

    if(!getNextAcceptedSample(iVar, tmp, mean, invCov, prob)) {
      return;
    }

    x1 = tmp;
    lp1 = prob.lnValue();

    if(!getNextAcceptedSample(iVar, tmp, mean, invCov, prob)) {
      return;
    }

    x2 = tmp;
    lp2 = prob.lnValue();

    // Arrange the samples in increasing order

    val1 = (x1[iVar] + x0[iVar])/2;
    val2 = (x2[iVar] + x1[iVar])/2;

    // Calculate dlogp/dx at each point

    dlp1 = (lp1 - lp0) / (x1[iVar] - x0[iVar]);
    dlp2 = (lp2 - lp1) / (x2[iVar] - x1[iVar]);

    // Finally, calculate d2logp/dx2

    d2 = (dlp2 - dlp1) / (val2 - val1);

    newsigval = 1.0/sqrt(-d2);

    // Can't invert nans

    if(!finite(newsigval)) {
      return;
    }

    newSigmas[iVar] = newsigval;
  }

  mm_.setSamplingSigmas(newSigmas);
  mm_.updateSamplingSigmas();
}

/**.......................................................................
 * Update the jumping distribution by evaluating the Hessian
 */
bool RunManager::updateHessian2(Probability startProb)
{
  // Store the current vector of values

  Vector<double> storedVal = mm_.currentSample_;
  Vector<double> mean      = mm_.mean_;
  Vector<double> sigmas    = mm_.sigma_;
  Matrix<double> invCov    = mm_.invCov_;
  double detC              = mm_.detC_;

  unsigned nVar = storedVal.size();
  Vector<double> x0(nVar), x1(nVar), x2(nVar), tmp(nVar);
  double lp0, lp1, lp2, dlp1, dlp2, d2;
  double val1, val2;
  double xlo, xmid, xhi;
  double lplo, lpmid, lphi;

  Vector<double> newSigmas = sigmas;
  double newsigval;

  for(unsigned iVar=0; iVar < nVar; iVar++) {

    // Generate a new random point from this variate's conditional
    // distribution

    tmp = storedVal;
    Probability prob = startProb;

    x0 = tmp;
    lp0 = prob.lnValue();

    if(!getNextAcceptedSample(iVar, tmp, mean, invCov, prob)) {
      return false;
    }

    x1 = tmp;
    lp1 = prob.lnValue();

    if(x0[iVar] < x1[iVar]) {
      xlo  = x0[iVar];
      lplo = lp0;

      xhi  = x1[iVar];
      lphi = lp1;
    } else {
      xlo  = x1[iVar];
      lplo = lp1;

      xhi  = x0[iVar];
      lphi = lp0;
    }

    if(!getNextAcceptedSample(iVar, tmp, mean, invCov, prob)) {
      return false;
    }

    x2  = tmp;
    lp2 = prob.lnValue();

    if(x2[iVar] < xlo) {
      xmid  = xlo;
      lpmid = lplo;

      xlo  = x2[iVar];
      lplo = lp2;
    } else if(x2[iVar] > xhi) {
      xmid  = xhi;
      lpmid = lphi;

      xhi  = x2[iVar];
      lphi = lp2;
    } else {
      xmid  = x2[iVar];
      lpmid = lp2;
    }

    val1 = (xmid +  xlo)/2;
    val2 = (xhi  + xmid)/2;

    // Calculate dlogp/dx at each point

    dlp1 = (lpmid -  lplo) / (xmid -  xlo);
    dlp2 = (lphi  - lpmid) / (xhi  - xmid);

    // Finally, calculate d2logp/dx2

    d2 = (dlp2 - dlp1) / (val2 - val1);

    newsigval = 1.0/sqrt(-d2);

    // Can't invert nans

    if(!finite(newsigval)) {
      return false;
    }

    newSigmas[iVar] = newsigval;
  }

  mm_.setSamplingSigmas(newSigmas);
  mm_.updateSamplingSigmas();

  return true;
}

/**.......................................................................
 * Update the jumping distribution by evaluating the Hessian
 */
bool RunManager::updateHessian2Mean(Probability startProb)
{
  double gamma = 0.1;

  // Store the current vector of values

  Vector<double> storedVal = mm_.currentSample_;
  Vector<double> mean      = mm_.mean_;
  Vector<double> sigmas    = mm_.sigma_;
  Matrix<double> invCov    = mm_.invCov_;
  double detC              = mm_.detC_;

  unsigned nVar = storedVal.size();
  Vector<double> x0(nVar), x1(nVar), x2(nVar), tmp(nVar);
  double lp0, lp1, lp2, dlp1, dlp2, d2;
  double val1, val2;
  double xlo, xmid, xhi;
  double lplo, lpmid, lphi;

  Vector<double> newSigmas = sigmas;
  Vector<double> newMeans  = mean;

  double newsigval;

  for(unsigned iVar=0; iVar < nVar; iVar++) {

    // Generate a new random point from this variate's conditional
    // distribution

    tmp = storedVal;
    Probability prob = startProb;

    x0 = tmp;
    lp0 = prob.lnValue();

    if(!getNextAcceptedSample(iVar, tmp, mean, invCov, prob)) {
      return false;
    }

    x1 = tmp;
    lp1 = prob.lnValue();

    if(x0[iVar] < x1[iVar]) {
      xlo  = x0[iVar];
      lplo = lp0;

      xhi  = x1[iVar];
      lphi = lp1;
    } else {
      xlo  = x1[iVar];
      lplo = lp1;

      xhi  = x0[iVar];
      lphi = lp0;
    }

    if(!getNextAcceptedSample(iVar, tmp, mean, invCov, prob)) {
      return false;
    }

    x2  = tmp;
    lp2 = prob.lnValue();

    if(x2[iVar] < xlo) {
      xmid  = xlo;
      lpmid = lplo;

      xlo  = x2[iVar];
      lplo = lp2;
    } else if(x2[iVar] > xhi) {
      xmid  = xhi;
      lpmid = lphi;

      xhi  = x2[iVar];
      lphi = lp2;
    } else {
      xmid  = x2[iVar];
      lpmid = lp2;
    }

    val1 = (xmid +  xlo)/2;
    val2 = (xhi  + xmid)/2;

    // Calculate dlogp/dx at each point

    dlp1 = (lpmid -  lplo) / (xmid -  xlo);
    dlp2 = (lphi  - lpmid) / (xhi  - xmid);

    // Finally, calculate d2logp/dx2

    d2 = (dlp2 - dlp1) / (val2 - val1);

    newsigval = 1.0/sqrt(-d2);

    // Can't invert nans

    if(!finite(newsigval)) {
      return false;
    }

    newSigmas[iVar] = newsigval;

    // And the estimated mean of the distribution is offset from the
    // current sample position by the first derivative

    // Estimate dlogp/dx at the midpoint by interpolating between
    // these two

    // X-values at the points where the first derivatives were
    // evaluated

    double xval1 = (xmid +  xlo)/2;
    double xval2 = (xhi  + xmid)/2;

    double dlp0 = dlp1 + (dlp2 - dlp1) / (xval2 - xval1) * (xmid - xval1);
    double s2  = -1.0/d2;
    double newmeanval = xmid + gamma * dlp0 * s2;
    newMeans[iVar]  = newmeanval;
  }

  mm_.setSamplingSigmas(newSigmas);
  mm_.updateSamplingSigmas();

  COUT("Setting sampling means from " << mm_.currentSample_ << " to " << newMeans);

  mm_.previousSample_ = newMeans;
  mm_.setSamplingMeans(newMeans);
  mm_.updateSamplingMeans();

  mm_.currentSample_ = newMeans;
  mm_.setValues(tmp, true);

  return true;
}

/**.......................................................................
 * Version of updateHessian() that samples the surface locally
 * (regardless of probability of sampled points) and updates the mean
 * of the jumping distribution accordingly.
 */
bool RunManager::updateHessian3(Probability startProb, double frac, double gamma)
{
  // Store the current vector of values

  Vector<double> storedValues = mm_.currentSample_;
  Vector<double> means        = mm_.mean_;
  Vector<double> sigmas       = mm_.sigma_;

  unsigned nVar = storedValues.size();
  double lp0, lp1, lp2, dlp0, dlp1, dlp2, d2, s2;
  double xval1, xval2;
  double xlo, xmid, xhi;
  double lplo, lpmid, lphi;

  Vector<double> newSigmas = sigmas;
  Vector<double> newMeans  = means;
  Vector<double> sample;

  double newsigval;
  double newmeanval;

  //------------------------------------------------------------
  // Iterate over all variates, sampling the likelihood surface at two
  // points about the current value to determine the curvature of the
  // likelihood surface
  //------------------------------------------------------------

  Probability prob;
  double fracIter;

  for(unsigned iVar=0; iVar < nVar; iVar++) {

    double mean  = means[iVar];
    double sigma = sigmas[iVar];

    // For each variate, we start at the last accepted point and its
    // corresponding probability

    sample = storedValues;
    prob   = startProb;

    // Midpoint of the sample will be the starting point

    xmid  = sample[iVar];
    lpmid = prob.lnValue();
    
    // Low sample will be at xmid - frac*sigma.  We iterate in case
    // xmid - frac*sigma happens to walk off the edge of our prior, in
    // which case prob = 0.0.  In this case, we shrink the step size
    // until we get a valid sample with non-zero probability

    fracIter = frac;
    do {
      xlo = xmid - fracIter*sigma;
      sample[iVar] = xlo;
      getNextSample(sample, prob);
      lplo = prob.lnValue();
      fracIter /= 2;
    } while(!(prob > 0.0));

    // High sample will be at xmid + frac*sigma.  We iterate in case
    // xmid + frac*sigma happens to walk off the edge of our prior, in
    // which case prob = 0.0.  In this case, we shrink the step size
    // until we get a valid sample with non-zero probability

    fracIter = frac;
    do {
      xhi = xmid + fracIter*sigma;
      sample[iVar] = xhi;
      getNextSample(sample, prob);
      lphi = prob.lnValue();
      fracIter /= 2;
    } while(!(prob > 0.0));

    // Now we have all three samples -- calculate derivatives

    // X-values at the points where the first derivatives were
    // evaluated

    xval1 = (xmid +  xlo)/2;
    xval2 = (xhi  + xmid)/2;

    // Calculate dlogp/dx at each of these points

    dlp1 = (lpmid -  lplo) / (xmid -  xlo);
    dlp2 = (lphi  - lpmid) / (xhi  - xmid);

    // Estimate dlogp/dx at the midpoint by interpolating between
    // these two

    dlp0 = dlp1 + (dlp2 - dlp1) / (xval2 - xval1) * (xmid - xval1);

    // Finally, calculate d2logp/dx2

    d2 = (dlp2 - dlp1) / (xval2 - xval1);

    // Estimate of sigma^2 is 1.0/(-d2logp/dx2);

    s2  = -1.0/d2;

    // sigma is sqrt() of that

    newsigval  = sqrt(s2);

    // And the estimated mean of the distribution is offset from the
    // current sample position by the first derivative

    newmeanval = xmid + gamma * dlp0 * s2;

    // Set the mew sampling means and sigmas, but only if we have
    // detected positive curvature (s2 > 0.0)

    if(s2 > 0.0) {
      newSigmas[iVar] = newsigval;
      newMeans[iVar]  = newmeanval;
    }
  }

  mm_.previousSample_ = newMeans;
  mm_.setSamplingMeans(newMeans);

  if(lpmid > 2 * nVar * -80.0) {
    mm_.setSamplingSigmas(newSigmas);
    mm_.updateSamplingSigmas();

    return true;
  }

  return false;
}

/**.......................................................................
 * Sample from the conditional distribution of the specified variate
 * (iVar) until an accepted sample is found
 */
bool RunManager::getNextAcceptedSample(unsigned iVar, Vector<double>& tmp, Vector<double>& mean, 
				       Matrix<double>& invCov, Probability& prob)
{
  Probability rat, propDensCurr, likeCurr;

  // Store the values on entry -- we will reset at the start of each
  // iteration until an accepted draw is found

  Probability storedProb    = prob;
  Vector<double>& storedVal = tmp;

  // Iterate drawing samples from the conditional distribution of the
  // requested variate until one is accepted

  bool accepted = false;
  unsigned ntry=1000, itry=0;

  do {

    tmp = storedVal;
    storedProb = prob;

    double alpha = Sampler::generateUniformSample(0.0, 1.0);
    
    dm_.clearModel();

    mm_.previousSample_ = tmp;
    Sampler::multiVariateSampleIterator(iVar, tmp, mean, invCov);
    mm_.currentSample_ = tmp;
    mm_.setValues(tmp, true);
    
    propDensCurr = mm_.priorPdf();
    rat = propDensCurr/prob;
    
    if(rat > alpha) {

      likeCurr = dm_.likelihood(mm_);
      rat *= likeCurr;
      
      if(rat > alpha) {
	accepted = true;
	prob = likeCurr*propDensCurr;
	return true;
      } else {
	accepted = false;
	tmp  = storedVal;
	prob = storedProb;
	mm_.revert();
      }
    }

    ++itry;
    
  } while(!accepted && itry < ntry);

  return false;
}

/**.......................................................................
 * Calculate the probability of the current sample
 */
void RunManager::getNextSample(Vector<double>& sample, Probability& prob)
{
  dm_.clearModel();

  mm_.currentSample_  = sample;
  mm_.setValues(sample, true);
    
  prob = mm_.priorPdf() * dm_.likelihood(mm_);    
}

/**.......................................................................
 * Compute chisq for the current model
 */
void RunManager::computeChisq()
{
  assertCurrentModel("Chisq computation", false);
  mm_.printModel("Current");
  ChisqVariate chisq = dm_.computeChisq();
  COUTCOLOR(std::endl << chisq, "green");
}

void RunManager::displayModel()
{
  dm_.clearModel();
  dm_.addModel(mm_);

  dm_.display();

  PgUtil::setZmin(0.0);
  PgUtil::setZmax(0.0);

  dm_.displayCompositeModel();
}

/**.......................................................................
 * Return the minimum chisq
 */
ChisqVariate RunManager::getMinimumChisq()
{
  assertBestFitModel();
  return dm_.computeChisq();
}

/**.......................................................................
 * Return the natural log of the Bayesian evidence
 */
double RunManager::getLnEvidence()
{
  return mm_.estimateLnEvidence();
}

/**.......................................................................
 * Take an input value string for a variable assignment, and see if
 * the value is a symbolic name.  If so, get the value of the variable
 * that that symbolic name corresponds to and return that value
 * instead.
 */
String RunManager::getVariableValueString(String& inputValStr)
{
  String valStr = inputValStr;

  if(valStr.contains(".")) {

    String ownerName, varName, aspectName;
    parseVarname(valStr, ownerName, varName, aspectName);

    // Now check if the owner is a model or a dataset, or neither
    
    Model*   model   = 0;
    DataSet* dataSet = 0;
    
    //------------------------------------------------------------
    // Check for models first
    //------------------------------------------------------------

    try {
      model = mm_.getModel(ownerName.str());
    } catch(Exception& err) {
    }

    // If this was a valid model.variate specification, return the
    // value of that variable as a (possibly) unit'd string
    
    if(model) {

      // Check for model variates first

      try {
	std::ostringstream os;
	os << ownerName << "." << varName;

	Variate* var  = mm_.getVar(os.str());

	//------------------------------------------------------------
	// If the specified model component is variable , this has a
	// special meaning (derive one from the other), that will not
	// be handled here.  In this case, just return the name of
	// the variate as specified
	//
	// Else if the component is fixed, return the value instead
	//------------------------------------------------------------

	if(var->isVariable()) {
	  return inputValStr;
	} else {
	  os.str("");
	  os << var->getUnitVal() << " " << var->units();
	  valStr = os.str();
	  return valStr;
	}
      } catch(Exception& err) {
      }

      // Check for model parameters next

      try {
	std::ostringstream os;
	os << ownerName << "." << varName;
	ParameterManager::Parameter* par  = mm_.getParameter(os.str(), true);
	os.str("");
	os << par->data_ << " " << par->units_;
	valStr = os.str();
	return valStr;
      } catch(Exception& err) {
      }

    }
    
    //------------------------------------------------------------
    // If we got here, then no model matches -- check datasets
    //------------------------------------------------------------
    
    try {
      dataSet = dm_.getDataSet(ownerName.str());
      DataSet::Parameter* par = dataSet->getParameter(varName.str(), true);
    } catch(Exception& err) {
    }
    
    //------------------------------------------------------------
    // If this was a valid dataset.parameter specification, return the
    // value of that variable as a (possibly) unit'd string
    //------------------------------------------------------------

    if(dataSet) {
      DataSet::Parameter* par = dataSet->getParameter(varName.str(), true);
      std::ostringstream os;
      os << par->data_ << " " << par->units_;
      valStr = os.str();

      return valStr;
    }
  }

  //------------------------------------------------------------
  // If we got here, then either the value string doesn't look like a
  // symbolic name, or doesn't resolve to any symbols we know about,
  // so pass it back unmodified
  //------------------------------------------------------------

  return inputValStr;
}

/**.......................................................................
 * Return true if this assignment depends on a dataset, in which case
 * we must process it after datasets have been read in.
 */
bool RunManager::assignmentDependsOnDataset(String& inputVal)
{
  String valStr = inputVal;

  if(valStr.contains(".")) {
    String ownerName, varName, aspectName;
    parseVarname(valStr, ownerName, varName, aspectName);

    try {
      DataSet* dataSet = dm_.getDataSet(ownerName.str());
      return true;
    } catch(Exception& err) {
    }
  }

  return false;
}

/**.......................................................................
 * Return true if this assignment depends on a model, in which case
 * we must process it after models have been loaded.
 */
bool RunManager::assignmentDependsOnModelLoad(String& inputVal)
{
  String valStr = inputVal;

  if(valStr.contains("[order]")) {
    return true;
  }

  if(valStr.contains("[range]")) {
    return true;
  }

  if(valStr.contains("[prior]")) {
    return true;
  }

  return false;
}

/**.......................................................................
 * Load any data sets that were specified
 */
void RunManager::loadData()
{
  dm_.loadData(genData_);
  datasetsInitialized_ = true;
}

/**.......................................................................
 * Remove any fixed models from the datasets
 */
void RunManager::removeModels()
{
  dm_.clearModel();
  dm_.remModel(mm_);
  dm_.clearModel();
}

/**.......................................................................
 * Process any lines that depend on datasets having been read in
 */
void RunManager::processDatasetDependentLines()
{
  for(unsigned iLine=0; iLine < datasetDependentLines_.size(); iLine++) {
    String line(datasetDependentLines_[iLine]);
    processLine(line);
  }
}

/**.......................................................................
 * Process any lines that depend on models being defined
 */
void RunManager::processModelDependentLines()
{
  modelsInitialized_ = true;
  for(unsigned iLine=0; iLine < modelDependentLines_.size(); iLine++) {
    String line(modelDependentLines_[iLine]);
    processLine(line);
  }
}

/**.......................................................................
 * Process any lines that depend on models being defined
 */
void RunManager::processModelLoadDependentLines()
{
  for(unsigned iLine=0; iLine < modelLoadDependentLines_.size(); iLine++) {
    String line(modelLoadDependentLines_[iLine]);
    processLine(line);
  }
}

/**.......................................................................
 * Return true if we count this character as 'numeric'
 */
bool RunManager::isNumeric(char c)
{
  return isdigit(c) || c=='-' || c=='+';
}

/**.......................................................................
 * Return true if our first data set is a 1D dataset
 */
bool RunManager::is1DDataSet()
{
  DataSet* ds = dm_.dataSetMap_.begin()->second;
  return ds->dataSetType_ & DataSetType::DATASET_1D;
}

/**.......................................................................
 * Return a vector of residuals for 1D data
 */
std::vector<double> RunManager::get1DResiduals()
{
  std::vector<double> res;

  if(dm_.dataSetMap_.empty()) {
    ThrowError("No data sets have been loaded");
    return res;
  }

  DataSet* ds = dm_.dataSetMap_.begin()->second;

  if(ds->dataSetType_ & DataSetType::DATASET_1D) {
    DataSet1D* ds1d = (DataSet1D*)ds;
    return ds1d->getResiduals();
  } else {
    ThrowError("Data set is not 1D");
    return res;
  }
}

/**.......................................................................
 * Return a vector of x values for 1D data
 */
std::vector<double> RunManager::get1DXData()
{
  std::vector<double> x;

  if(dm_.dataSetMap_.empty()) {
    ThrowError("No data sets have been loaded");
    return x;
  }

  DataSet* ds = dm_.dataSetMap_.begin()->second;

  if(ds->dataSetType_ & DataSetType::DATASET_1D) {
    DataSet1D* ds1d = (DataSet1D*)ds;
    return ds1d->getXData();
  } else {
    ThrowError("Data set is not 1D");
    return x;
  }
}

std::vector<std::string> RunManager::listModelNames()
{
  return mm_.listModelNames("Maximum Likelihood", nBin_, Model::STAT_MODE);
}

std::vector<std::string> RunManager::listModelUnits()
{
  return mm_.listModelUnits("Maximum Likelihood", nBin_, Model::STAT_MODE);
}

std::vector<std::string> RunManager::listModelStrVals()
{
  return mm_.listModelStrVals("Maximum Likelihood", nBin_, Model::STAT_MODE);
}

std::vector<Model::VarVal> RunManager::listModelVals()
{
  return mm_.listModelVals("Maximum Likelihood", nBin_, stat_, nSigma_);
}

std::map<std::string, Model::VarVal> RunManager::listModel()
{
  return mm_.listModel("Maximum Likelihood", nBin_, stat_, nSigma_);
}

unsigned RunManager::getChainLength()
{
  return mm_.getChainLength();
}

/**.......................................................................
 * Process a line of the form:
 *
 * dataset.exclude(model);
 */
void RunManager::getExcludeArgs(String& line)
{
  line.strip(' ');
  String dataSetStr = line.findNextInstanceOf("", false, ".exclude", true, true);
  String modelStr   = line.findNextInstanceOf("(", true, ")", true, true);

  DataSet* dataSet = 0;
  Model* model = 0;

  dataSet = dm_.getDataSet(dataSetStr.str());
  model   = mm_.getModel(modelStr.str());

  dataSet->exclude(model);
}

/**.......................................................................
 * Process an addmodel line
 */
void RunManager::addModel(String& line, bool remove)
{
  String name, type, file, discard;
  getModelArgs(line, name, type, file, discard);

  checkIfNameAlreadyExists(name.str());

  //------------------------------------------------------------
  // Now check how the model was specified.  If a name and type were
  // given, add the model normally.  If a file was specified, assume
  // that we are loading a model from an output file.
  //------------------------------------------------------------

  if(!type.isEmpty()) {
    Model* model = mm_.addModel(type.str(), name.str(), remove);
    model->setThreadPool(modelPool_);
  } else {
    mm_.loadOutputFile(file.str(), name.str(), discard.isEmpty() ? 0.0 : discard.toInt());
    runMarkov_ = false;
    loadOutputFile_ = true;
  }
}

/**.......................................................................
 * Process an adddataset line
 */
void RunManager::addDataset(String& line)
{
  std::string name, type;

  getDataSetNameAndType(line, name, type);
  checkIfNameAlreadyExists(name);

  DataSet* dataSet = dm_.addDataSet(type, name);
  dataSet->setThreadPool(dataPool_);
}

/**.......................................................................
 * Print a message about our progress
 */
void RunManager::printProgress(unsigned i)
{
  static bool printNewline = true;
  static unsigned ndigit = ceil(log10((double)nTry_)) + 1;
  static ostringstream filler;

  if(i==0) {
    unsigned nSpace = 31 + (ndigit+1)*2;
    filler << '\r';
    for(unsigned iSpace=0; iSpace < nSpace; iSpace++) 
      filler << ' ';
    filler << std::ends;
  }

  //------------------------------------------------------------
  // We will print progress at 1% intervals
  //------------------------------------------------------------

  static unsigned nPrint = nTry_ > 100 ? (nTry_ / 100) : 1;

  if((i+1) % nPrint == 0 || i==nBurn_-1 || i==nTry_-1) {
    
    if(i > nBurn_-1 && printNewline) {
      COUT("");
      printNewline = false;
    }
    
    if(i > nBurn_-1) {
      COUTCOLORNNL(std::cout << filler.str()
		   << "\rIteration: " << std::setw(ndigit) << std::right << (i+1) << "/" << nTry_ << " (" << std::right << setprecision(0) << std::fixed << (100*(double)(i+1)/nTry_) << "%)", "green");
    } else {
      COUTCOLORNNL(std::cout << filler.str()
		   << "\rBurn-in:   " << std::setw(ndigit) << std::right << (i+1) << "/" << nTry_ << " (" << std::right << setprecision(0) << std::fixed << (100*(double)(i+1)/nTry_) << "%)", "green");
    }
    
    fflush(stdout);
  }
}

/**.......................................................................
 * Generate the next sample in the Markov chain
 */
void RunManager::generateNewSample()
{
  //------------------------------------------------------------
  // Clear the models for all datasets
  //------------------------------------------------------------
  
  dm_.clearModel();
  
  //------------------------------------------------------------
  // Generate a new sample for all model parameters
  //------------------------------------------------------------
  
  sampleTimer_.start();

  mm_.sample();

  sampleTimer_.stop();
  sampleTime_ += sampleTimer_.deltaInSeconds();

  ++nTrySinceLastUpdate_;
}

/**.......................................................................
 * Return true if the current sample would be accepted by the MH algorithm
 */
bool RunManager::acceptMetropolisHastings(unsigned i,
					  Probability& propDensCurr, Probability& propDensPrev, 
					  Probability& likeCurr,     Probability& likePrev, 
					  ChisqVariate& chisq)
{
  //------------------------------------------------------------
  // Generate the uniform variate we will use to accept or reject
  // this sample
  //------------------------------------------------------------
  
  double alpha = Sampler::generateUniformSample(0.0, 1.0);
  
  //-----------------------------------------------------------
  // Calculate the pdf of the proposed sample
  //------------------------------------------------------------

  propDensCurr = mm_.priorPdf();

#if PRIOR_DEBUG
  COUT("propDensCurr = " << propDensCurr);
#endif

  //------------------------------------------------------------
  // We don't accept this sample if it is prohibited by the prior
  //------------------------------------------------------------

  if(!(propDensCurr > 0)) { 
#if PRIOR_DEBUG
    COUT("Returning because propdens = 0: sample = " << mm_.currentSample_);
#endif
    return false;
  }

  //------------------------------------------------------------
  // The ratio 
  //
  //        P(D|H1) * P(H1)
  //  r  =  ---------------
  //        P(D|H2) * P(H2)
  // 
  // is used to determine whether the new proposed sample (H2) will
  // be accepted.  In the Metropolis-Hastings algorithm, this ratio
  // is compared to a uniform deviate alpha and the sample is
  // accepted if rat > alpha.  
  //
  //------------------------------------------------------------

  likeTimer_.start();

  dm_.likelihood(mm_, likeCurr, chisq);

  likeTimer_.stop();
  likeTime_ += likeTimer_.deltaInSeconds();

  if(i > 0) {

    Probability rat = (propDensCurr/propDensPrev) * (likeCurr/likePrev);

    if(rat > alpha) {
#if PRIOR_DEBUG
      COUT("Accepted because rat = " << rat << " alpha = " << alpha << " prop = " << propDensCurr << " like - " << likeCurr);
#endif
      return true;
    } else {
#if PRIOR_DEBUG
      COUT("Returning because rat = " << rat << " alpha = " << alpha << " prop = " << propDensCurr << " prev = " << propDensPrev << " like - " << likeCurr << " prev = "  << likePrev);
#endif
      return false;
    }
  } else {
    return true;
  }
}

/**.......................................................................
 * Return true if we should attempt to tune the jumping distribution
 * on this iteration
 */
bool RunManager::timeToTune(unsigned i)
{
  return (i > 0) && (i < nBurn_) && (i - iLastUpdate_ >= nPerUpdate_);
}

/**.......................................................................
 * Tune the jumping distribution
 */
void RunManager::tuneJumpingDistribution(unsigned i, Probability& propDensCurr, Probability& likeCurr)
{
  double acceptFrac = (double)(nAcceptedSinceLastUpdate_)/nTrySinceLastUpdate_;

  //------------------------------------------------------------
  // Accepted fraction is lower than we want -- use quadratic
  // estimator to tune the jumping distribution, or narrow the
  // jumping distribution if we are already in a
  // high-probability part of the posterior
  //------------------------------------------------------------

  if(acceptFrac < acceptFracLowLim_) {

    //------------------------------------------------------------
    // If we have already tuned the distribution several times,
    // just narrow the width to increase the number of samples
    // accepted
    //------------------------------------------------------------

    if(nUpdate_ > 3) {
	    
      Vector<double> currentSigmas = mm_.sigma_;
      Vector<double> newSigmas = currentSigmas / 2;
      mm_.setSamplingSigmas(newSigmas);
      mm_.updateSamplingSigmas();

      //------------------------------------------------------------
      // Else tune the distribution using a quadratic estimator
      // of the posterior
      //------------------------------------------------------------

    } else {

      tuneTimer_.start();

#if 1
      if(updateHessian2(likeCurr*propDensCurr)) {
	++nUpdate_;
	foundUpdate_ = true;
      }
#else
      if(updateMethod_ == 1) {
	if(updateHessian2(likeCurr*propDensCurr)) {
	  ++nUpdate_;
	  foundUpdate_ = true;
	}
      } else if(updateMethod_ == 2) {
	COUT("Using mean update method: nacc = " << mm_.nAccepted_);
	if(updateHessian2Mean(likeCurr*propDensCurr)) {
	  ++nUpdate_;
	  foundUpdate_ = true;
	}
      } else {
	COUT("Using no update method: nacc = " << mm_.nAccepted_);
	nUpdate_ = 3;
	foundUpdate_ = true;
      }
#endif

      tuneTimer_.stop();
      tuneTime_ += tuneTimer_.deltaInSeconds();
    }

    //------------------------------------------------------------
    // Else accepted fraction is higher than we want.  Widen the
    // jumping distribution so that more samples are rejected
    //------------------------------------------------------------

  } else if(acceptFrac > acceptFracHighLim_) {
    Vector<double> currentSigmas = mm_.sigma_;
    Vector<double> newSigmas = currentSigmas * 2;
    mm_.setSamplingSigmas(newSigmas);
    mm_.updateSamplingSigmas();
  }

  nAcceptedSinceLastUpdate_ = 0;
  nTrySinceLastUpdate_      = 0;
  iLastUpdate_ = i;
}

double RunManager::lnLikelihood()
{
  dm_.likelihood(mm_, convProb_, convChisq_);
  return convProb_.lnValue();
}

/**.......................................................................
 * Calculate the probability of the passed sample
 */
double RunManager::lnLikelihood(Vector<double>& sample)
{
  //------------------------------------------------------------
  // Clear all models that are currently set for datasets
  //------------------------------------------------------------

  dm_.clearModel();

  //------------------------------------------------------------
  // Set the external sample
  //------------------------------------------------------------

  mm_.externalSample(sample);

  //------------------------------------------------------------
  // And compute the likelihood of the current sample
  //------------------------------------------------------------

  dm_.likelihood(mm_, convProb_, convChisq_);

  return convProb_.lnValue();
}

/**.......................................................................
 * Calculate the probability of the passed sample, in unit'd values
 */
double RunManager::lnLikelihoodUnits(Vector<double>& sample)
{
  //------------------------------------------------------------
  // Clear all models that are currently set for datasets
  //------------------------------------------------------------

  dm_.clearModel();

  //------------------------------------------------------------
  // Set the external sample
  //------------------------------------------------------------

  mm_.externalSampleUnits(sample);

  //------------------------------------------------------------
  // And compute the likelihood of the current sample
  //------------------------------------------------------------

  dm_.likelihood(mm_, convProb_, convChisq_);

  return convProb_.lnValue();
}

void RunManager::initializeForMarkovChain()
{
  mm_.store_ = false;
  mm_.initializeForMarkovChain(1e6, 1e6, "");
}

std::vector<Variate*>& RunManager::getVariableComponents()
{
  return mm_.variableComponents_;
}

String RunManager::getStrippedVal(String& line)
{
  String tok,val;
  getTokVal(line, tok, val);
      
  tok.strip(' ');
  val.strip(' ');

  return val;
}

void RunManager::printModel()
{
  mm_.printModel("Current");
}

void RunManager::fillChain(std::string varName, double** dPtr)
{
  mm_.fillChain(varName, dPtr);
}

void RunManager::fillMultiplicity(unsigned** dPtr)
{
  mm_.fillMultiplicity(dPtr);
}

void RunManager::fillLnLikelihood(double** dPtr)
{
  mm_.fillLnLikelihood(dPtr);
}

bool RunManager::checkConvergence(unsigned i)
{
  if(!runtoConvergence_)
    return false;

  if(!(i > nBurn_ && (i % nConverge_ == 0)))
    return false;

  return mm_.converged(targetVariance_);
}
