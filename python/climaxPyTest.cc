/**.......................................................................
 * Python file for reading from GCP archive data files.
 */
#include "gcp/fftutil/RunManager.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/util/DataType.h"
#include "gcp/util/Directives.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Logo.h"
#include "gcp/util/Timer.h"

#include "gcp/models/ArnaudModel.h"
#include "gcp/models/CosmologyModel.h"

#include "gcp/python/PyParser.h"
#include "gcp/python/PyHandler.h"

#include "cpgplot.h"

#include <Python.h>
#include <iostream>  
#include <vector>  
#include <map>  

#if DIR_HAVE_NUMPY
#include "arrayobject.h"
#endif

using namespace std;
using namespace gcp::python;
using namespace gcp::util;
using namespace gcp::models;

PyObject* buildVarDict(Model::VarVal& val);
PyObject* buildVarDict(Model::VarVal& val, unsigned len, char** cPtr);

std::map<std::string, Model::VarVal> getModelMap();
std::map<std::string, Model::VarVal> listModel(std::string prefix, unsigned nBin, Model::Stat stat);

#if DIR_HAVE_NUMPY
int numpyTypeOf(gcp::util::DataType::Type dataType);
std::vector<npy_intp> getNumPyInds(std::vector<int>& dims);
#endif

PyObject* getArray(std::vector<int> dims, gcp::util::DataType::Type dataType, char** data);
PyObject* getArray(unsigned len, gcp::util::DataType::Type dataType, char** data);

/**.......................................................................
 * Entry point from the Python environment
 *
 * SPYDOC
 *
 * Run climax with the specified parameter file
 * 
 * Use like: 
 * 
 *    dict = climaxPyTest.run(runFile)
 *
 * Where:
 *
 *   dict    - A dictionary of returned parameter estimates and chi-squared
 *   runFile - The text file controlling the Markov chain run
 *
 * EPYDOC
 */
static PyObject* run(PyObject* self, PyObject* args)
{
  Timer timer;
  timer.start();

  PyObject* obj = 0;

  std::ostringstream usage;

  usage << "Usage: dict = " << "climaxPyTest."
	<< "run(runfile)";

  try {

    PyParser modelFileParser(PyParser::getArrayItem(args, 0));
    std::string runFile = modelFileParser.getString();
  
    gcp::util::RunManager rm;

    rm.setRunFile(runFile);
    rm.run();

    obj = PyDict_New();

    try {
      //------------------------------------------------------------
      // Add the chisq
      //------------------------------------------------------------
      
      ChisqVariate chisq = rm.getMinimumChisq();

      PyObject* chisqObj = Py_BuildValue("f", chisq.reducedChisq());
      PyDict_SetItemString(obj, "chisq", chisqObj);
      Py_DECREF(chisqObj);

      PyObject* nDofObj = Py_BuildValue("f", (float)chisq.nDof());
      PyDict_SetItemString(obj, "nDof", nDofObj);
      Py_DECREF(nDofObj);

      //------------------------------------------------------------
      // Add the ln evidence
      //------------------------------------------------------------
      
      double lnEvidence = rm.getLnEvidence();
      PyObject* evidenceObj = Py_BuildValue("f", lnEvidence);
      PyDict_SetItemString(obj, "lnEvidence", evidenceObj);
      Py_DECREF(evidenceObj);

    } catch(Exception& err) {
      COUT("Caught an error: " << err.what());
    }

    //------------------------------------------------------------
    // Add the residuals
    //------------------------------------------------------------

    try {
      std::vector<double> x   = rm.get1DXData();
      std::vector<double> res = rm.get1DResiduals();

      std::vector<unsigned> dims;
      dims.push_back(res.size());

      PyObject* xObj   = PyHandler::getFlatList(dims);
      PyObject* resObj = PyHandler::getFlatList(dims);

      for(unsigned i=0; i < res.size(); i++) {
	PyObject* currXObj = PyFloat_FromDouble(x[i]);
	PyList_SetItem(xObj, i, currXObj);

	PyObject* currResObj = PyFloat_FromDouble(res[i]);
	PyList_SetItem(resObj, i, currResObj);
      }

      PyDict_SetItemString(obj, "x",         xObj);
      PyDict_SetItemString(obj, "residuals", resObj);

      Py_DECREF(xObj);
      Py_DECREF(resObj);

    } catch(...) {
    }

    //------------------------------------------------------------
    // Add entries for the best-fit model parameter values
    //------------------------------------------------------------

    std::map<std::string, Model::VarVal> valueMap = rm.listModel();

    for(std::map<std::string, Model::VarVal>::iterator iter=valueMap.begin(); iter != valueMap.end(); iter++) {
      Model::VarVal val = iter->second;

      PyObject* varObj = buildVarDict(val);

      std::string name = iter->first;
      PyDict_SetItemString(obj, name.c_str(), varObj);
    }

  } catch(Exception& err) {
    COUTCOLOR(err.what(), "red");
    COUTCOLOR(std::endl << usage.str() << std::endl, "green");
    obj = Py_BuildValue("c", ' ');
  }

  return obj;
}

static PyObject* getChain(PyObject* self, PyObject* args)
{
  Timer timer;
  timer.start();

  PyObject* obj = 0;

  std::ostringstream usage;

  usage << "Usage: dict = " << "climaxPyTest."
	<< "run(runfile)";

  try {

    if(PyParser::getSize(args) > 1) {
      PyParser printParser(PyParser::getArrayItem(args, 1));
      if(printParser.getBoolVal() == false) {
	Logger::installStdoutPrintFn(Logger::dontPrint);
      }
    }

    PyParser modelFileParser(PyParser::getArrayItem(args, 0));
    std::string runFile = modelFileParser.getString();
  
    gcp::util::RunManager rm;

    rm.setRunFile(runFile);
    rm.run();

    obj = PyDict_New();

    try {
      //------------------------------------------------------------
      // Add the chisq
      //------------------------------------------------------------
      
      ChisqVariate chisq = rm.getMinimumChisq();

      PyObject* chisqObj = Py_BuildValue("f", chisq.reducedChisq());
      PyDict_SetItemString(obj, "chisq", chisqObj);
      Py_DECREF(chisqObj);

      PyObject* nDofObj = Py_BuildValue("f", (float)chisq.nDof());
      PyDict_SetItemString(obj, "nDof", nDofObj);
      Py_DECREF(nDofObj);

      //------------------------------------------------------------
      // Add the ln evidence
      //------------------------------------------------------------
      
      double lnEvidence = rm.getLnEvidence();
      PyObject* evidenceObj = Py_BuildValue("f", lnEvidence);
      PyDict_SetItemString(obj, "lnEvidence", evidenceObj);
      Py_DECREF(evidenceObj);

    } catch(Exception& err) {
      COUT("Caught an error: " << err.what());
    }

    //------------------------------------------------------------
    // Add entries for the best-fit model parameter values
    //------------------------------------------------------------

    std::map<std::string, Model::VarVal> valueMap = rm.listModel();

    unsigned chainLength = rm.getChainLength();

    for(std::map<std::string, Model::VarVal>::iterator iter=valueMap.begin(); iter != valueMap.end(); iter++) {
      Model::VarVal val = iter->second;
      char* cPtr=0;
      PyObject* varObj  = buildVarDict(val, chainLength, &cPtr);
      std::string name  = iter->first;

      if(!val.isFixed_) {
	double* dPtr = (double*) cPtr;
	rm.fillChain(name, &dPtr);
      }

      PyDict_SetItemString(obj, name.c_str(), varObj);
    }

    //------------------------------------------------------------
    // Lastly, add entries for the lnLikelihood and multiplicity
    //------------------------------------------------------------

    char* cPtr=0;
    double* dPtr=0;
    PyObject* tmpObj = 0;

    tmpObj = getArray(chainLength, DataType::DOUBLE, &cPtr);
    dPtr = (double*) cPtr;
    rm.fillLnLikelihood(&dPtr);
    PyDict_SetItemString(obj, "lnLikelihood", tmpObj);

    tmpObj = getArray(chainLength, DataType::UINT, &cPtr);
    unsigned* uPtr = (unsigned*) cPtr;
    rm.fillMultiplicity(&uPtr);
    PyDict_SetItemString(obj, "multiplicity", tmpObj);
    
  } catch(Exception& err) {
    COUTCOLOR(err.what(), "red");
    COUTCOLOR(std::endl << usage.str() << std::endl, "green");
    obj = Py_BuildValue("c", ' ');
  }

  return obj;
}

static PyObject* runTest(PyObject* self, PyObject* args)
{
  Timer timer;
  timer.start();

  PyObject* obj = 0;

  std::ostringstream usage;

  usage << "Usage: dict = " << "climaxPyTest."
	<< "run(runfile)";

  try {

    PyParser modelFileParser(PyParser::getArrayItem(args, 0));
    std::string runFile = modelFileParser.getString();
  
    gcp::util::RunManager rm;

    rm.setRunFile(runFile);
    rm.run();

    obj = PyDict_New();
    Py_DECREF(obj);

    try {
      //------------------------------------------------------------
      // Add the chisq
      //------------------------------------------------------------
      
      ChisqVariate chisq = rm.getMinimumChisq();

      PyObject* chisqObj = Py_BuildValue("f", chisq.reducedChisq());
      PyDict_SetItemString(obj, "chisq", chisqObj);
      Py_DECREF(chisqObj);

      PyObject* nDofObj = Py_BuildValue("f", (float)chisq.nDof());
      PyDict_SetItemString(obj, "nDof", nDofObj);
      Py_DECREF(nDofObj);

      //------------------------------------------------------------
      // Add the ln evidence
      //------------------------------------------------------------
      
      double lnEvidence = rm.getLnEvidence();
      PyObject* evidenceObj = Py_BuildValue("f", lnEvidence);
      PyDict_SetItemString(obj, "lnEvidence", evidenceObj);
      Py_DECREF(evidenceObj);

    } catch(Exception& err) {
      COUT("Caught an error: " << err.what());
    }

    //------------------------------------------------------------
    // Add the residuals
    //------------------------------------------------------------

    try {
      std::vector<double> x   = rm.get1DXData();
      std::vector<double> res = rm.get1DResiduals();

      std::vector<unsigned> dims;
      dims.push_back(res.size());

      PyObject* xObj   = PyHandler::getFlatList(dims);
      PyObject* resObj = PyHandler::getFlatList(dims);

      for(unsigned i=0; i < res.size(); i++) {
	PyObject* currXObj = PyFloat_FromDouble(x[i]);
	PyList_SetItem(xObj, i, currXObj);

	PyObject* currResObj = PyFloat_FromDouble(res[i]);
	PyList_SetItem(resObj, i, currResObj);
      }

      PyDict_SetItemString(obj, "x",         xObj);
      PyDict_SetItemString(obj, "residuals", resObj);

      Py_DECREF(xObj);
      Py_DECREF(resObj);

    } catch(...) {
    }

    //------------------------------------------------------------
    // Add entries for the best-fit model parameter values
    //------------------------------------------------------------

    std::vector<Model::VarVal> vals    = rm.listModelVals();
    std::vector<std::string>   units   = rm.listModelUnits();
    std::vector<std::string>   names   = rm.listModelNames();
    std::vector<std::string>   strVals = rm.listModelStrVals();

    for(unsigned iVar=0; iVar < vals.size(); iVar++) {
      Model::VarVal& val = vals[iVar];
      COUT("Val " << iVar << " has name " << names[iVar]);
    }

    return obj;

  } catch(Exception& err) {
    COUTCOLOR(err.what(), "red");
    COUTCOLOR(std::endl << usage.str() << std::endl, "green");
    obj = Py_BuildValue("c", ' ');
  }

  return obj;
}

std::map<std::string, Model::VarVal> listModel(std::string prefix, unsigned nBin, Model::Stat stat)
{
  std::map<std::string, Model::VarVal> modelMap;

  for(unsigned iVar=0; iVar < 4; iVar++) {
    Model::VarVal val;

    val.stat_   = Model::STAT_MODE;
    val.refVal_ = 0.0;
    val.units_  = "test";
    
    std::ostringstream os;
    os << "name" << iVar;

    modelMap[os.str()] = val;
  }

  return modelMap;
}

std::map<std::string, Model::VarVal> getModelMap()
{
  std::map<std::string, Model::VarVal> modelMap;

  for(unsigned iVar=0; iVar < 4; iVar++) {
    Model::VarVal val;

    val.stat_   = Model::STAT_MODE;
    val.refVal_ = 0.0;
    val.units_  = "test";
    
    std::ostringstream os;
    os << "name" << iVar;

    modelMap[os.str()] = val;
  }

  return modelMap;
}

/**.......................................................................
 * Entry point from the Python environment
 *
 * SPYDOC
 *
 * Use like: 
 * 
 *    dict = climaxPyTest.run(runFile)
 *
 * Where:
 *
 *   dict    - A dictionary of returned parameter estimates and chi-squared
 *   runFile - The text file controlling the Markov chain run
 *
 * EPYDOC
 */
static PyObject* dAMpc(PyObject* self, PyObject* args)
{
  PyObject* obj = 0;

  PyParser zParser(PyParser::getArrayItem(args, 0));
  double z = zParser.getFloat();

  HubbleConstant H0;
  H0.setKmPerSecPerMpc(70);

  Cosmology cosmo;
  cosmo.setH0(H0);

  cosmo.setOmegaM(0.3);
  cosmo.setOmegaL(0.7);
  cosmo.setRedshift(z);

  try {
    Length DA = cosmo.angularDiameterDistance();

    obj = Py_BuildValue("f", DA.Mpc());
  } catch(Exception& err) {
    COUT(err.what());
    obj = Py_BuildValue("f", 0.0);
  }

  return obj;
}

/**.......................................................................
 * Entry point from the Python environment
 *
 * SPYDOC
 *
 * Use like: 
 * 
 *    rescale = climaxPyTest.arnaudRescale(m500, z)
 *
 * Where:
 *
 *   m500    -  The M500 for which to determine the rescale factor
 *   rescale -  The scale factor
 *
 * EPYDOC
 */
static PyObject* arnaudRescale(PyObject* self, PyObject* args)
{
  PyObject* obj = 0;

  PyParser mParser(PyParser::getArrayItem(args, 0));
  double m500 = mParser.getFloat();

  PyParser zParser(PyParser::getArrayItem(args, 1));
  double z = zParser.getFloat();

  ArnaudModel   arnaud;

  arnaud.specifyValue("m500", m500, "Msolar");
  arnaud.specifyValue("normalizationFrequency", 30.0, "GHz");
  arnaud.specifyValue("spectralType", "sz");

  arnaud.specifyDerivedVariate("Mtot500");

  CosmologyModel cosmo;

  cosmo.getVar("H0")->setVal(70, "km/s/Mpc");
  cosmo.getVar("H0")->wasSpecified_ = true;

  cosmo.getVar("z")->setVal(z, "");
  cosmo.getVar("z")->wasSpecified_ = true;

  cosmo.getVar("omegaM")->setVal(0.3, "");
  cosmo.getVar("omegaM")->wasSpecified_ = true;

  cosmo.getVar("omegaL")->setVal(0.7, "");
  cosmo.getVar("omegaL")->wasSpecified_ = true;

  cosmo.update();

  arnaud.initializeCosmology(&cosmo);

  arnaud.checkSetup();
  arnaud.updateVariableMap();

  arnaud.deriveVariates();

  arnaud.interpolateForRescale();

  try {
    double rescale = arnaud.rescale_.val_;
    obj = Py_BuildValue("f", rescale);
  } catch(Exception& err) {
    COUT(err.what());
    obj = Py_BuildValue("f", 0.0);
  }

  return obj;
}

static PyObject* arnaudY500(PyObject* self, PyObject* args)
{
  PyObject* obj = 0;

  PyParser mParser(PyParser::getArrayItem(args, 0));
  double m500 = mParser.getFloat();

  PyParser zParser(PyParser::getArrayItem(args, 1));
  double z = zParser.getFloat();

  PyParser rescaleParser(PyParser::getArrayItem(args, 2));
  double rescale = rescaleParser.getFloat();

  ArnaudModel   arnaud;

  arnaud.specifyValue("m500", m500, "Msolar");
  arnaud.specifyValue("normalizationFrequency", 30.0, "GHz");
  arnaud.specifyValue("spectralType", "sz");
  arnaud.specifyValue("rescale", rescale, "");

  arnaud.specifyDerivedVariate("Ysph500");

  CosmologyModel cosmo;

  cosmo.getVar("H0")->setVal(70, "km/s/Mpc");
  cosmo.getVar("H0")->wasSpecified_ = true;

  cosmo.getVar("z")->setVal(z, "");
  cosmo.getVar("z")->wasSpecified_ = true;

  cosmo.getVar("omegaM")->setVal(0.3, "");
  cosmo.getVar("omegaM")->wasSpecified_ = true;

  cosmo.getVar("omegaL")->setVal(0.7, "");
  cosmo.getVar("omegaL")->wasSpecified_ = true;

  cosmo.update();

  arnaud.initializeCosmology(&cosmo);

  arnaud.checkSetup();

  arnaud.updateVariableMap();
  arnaud.deriveVariates();

  try {
    double y500 = arnaud.ySph500_.squaredMpc();
    obj = Py_BuildValue("f", y500);
  } catch(Exception& err) {
    COUT(err.what());
    obj = Py_BuildValue("f", 0.0);
  }

  return obj;
}

static PyObject* arnaudYsph(PyObject* self, PyObject* args)
{
  PyObject* obj = 0;

  PyParser mParser(PyParser::getArrayItem(args, 0));
  double m500 = mParser.getFloat();

  PyParser zParser(PyParser::getArrayItem(args, 1));
  double z = zParser.getFloat();

  PyParser rescaleParser(PyParser::getArrayItem(args, 2));
  double rescale = rescaleParser.getFloat();

  PyParser outerParser(PyParser::getArrayItem(args, 3));
  double outerRadMpc = outerParser.getFloat();

  ArnaudModel   arnaud;

  arnaud.specifyValue("m500", m500, "Msolar");
  arnaud.specifyValue("normalizationFrequency", 30.0, "GHz");
  arnaud.specifyValue("spectralType", "sz");
  arnaud.specifyValue("rescale", rescale, "");

  arnaud.specifyParameter("innerRadius", 0.1, "Mpc");
  arnaud.specifyParameter("outerRadius", outerRadMpc, "Mpc");

  arnaud.specifyDerivedVariate("Ysph");

  CosmologyModel cosmo;

  cosmo.getVar("H0")->setVal(70, "km/s/Mpc");
  cosmo.getVar("H0")->wasSpecified_ = true;

  cosmo.getVar("z")->setVal(z, "");
  cosmo.getVar("z")->wasSpecified_ = true;

  cosmo.getVar("omegaM")->setVal(0.3, "");
  cosmo.getVar("omegaM")->wasSpecified_ = true;

  cosmo.getVar("omegaL")->setVal(0.7, "");
  cosmo.getVar("omegaL")->wasSpecified_ = true;

  cosmo.update();

  arnaud.initializeCosmology(&cosmo);

  arnaud.checkSetup();
  arnaud.updateVariableMap();

  arnaud.deriveVariates();

  try {
    double ySph = arnaud.ySph_.squaredMpc();
    obj = Py_BuildValue("f", ySph);
  } catch(Exception& err) {
    COUT(err.what());
    obj = Py_BuildValue("f", 0.0);
  }

  return obj;
}

static PyObject* arnaudR500(PyObject* self, PyObject* args)
{
  PyObject* obj = 0;

  PyParser mParser(PyParser::getArrayItem(args, 0));
  double m500 = mParser.getFloat();

  PyParser zParser(PyParser::getArrayItem(args, 1));
  double z = zParser.getFloat();

  PyParser rescaleParser(PyParser::getArrayItem(args, 2));
  double rescale = rescaleParser.getFloat();

  ArnaudModel   arnaud;

  arnaud.specifyValue("m500", m500, "Msolar");
  arnaud.specifyValue("normalizationFrequency", 30.0, "GHz");
  arnaud.specifyValue("spectralType", "sz");
  arnaud.specifyValue("rescale", rescale, "");


  arnaud.specifyDerivedVariate("R500");
  CosmologyModel cosmo;

  cosmo.getVar("H0")->setVal(70, "km/s/Mpc");
  cosmo.getVar("H0")->wasSpecified_ = true;

  cosmo.getVar("z")->setVal(z, "");
  cosmo.getVar("z")->wasSpecified_ = true;

  cosmo.getVar("omegaM")->setVal(0.3, "");
  cosmo.getVar("omegaM")->wasSpecified_ = true;

  cosmo.getVar("omegaL")->setVal(0.7, "");
  cosmo.getVar("omegaL")->wasSpecified_ = true;

  cosmo.update();

  arnaud.initializeCosmology(&cosmo);

  arnaud.checkSetup();
  arnaud.updateVariableMap();

  arnaud.deriveVariates();

  try {
    double ySph = arnaud.r500_.Mpc();
    obj = Py_BuildValue("f", ySph);
  } catch(Exception& err) {
    COUT(err.what());
    obj = Py_BuildValue("f", 0.0);
  }

  return obj;
}

/**.......................................................................
 * Build a variate dictionary
 */
PyObject* buildVarDict(Model::VarVal& val)
{
  PyObject* obj = PyDict_New();
  PyObject* tmpObj = 0;

  //------------------------------------------------------------
  // Build the val member of this dictionary, as a string or double,
  // depending on the variate type
  //------------------------------------------------------------

  if(val.isStr_) {
    tmpObj = Py_BuildValue("s", val.strVal_.c_str());
    PyDict_SetItemString(obj, "val", tmpObj);
  } else {
    tmpObj = Py_BuildValue("f", val.refVal_);
    PyDict_SetItemString(obj, "val", tmpObj);
  }

  //------------------------------------------------------------
  // If the component is fixed, just set the value, else the 1-sigma
  // condifidence interval
  //------------------------------------------------------------

  tmpObj = Py_BuildValue("b", val.isFixed_);
  PyDict_SetItemString(obj, "fixed", tmpObj);
    
  if(!val.isFixed_) {
    tmpObj = Py_BuildValue("f", val.lowVal_);
    PyDict_SetItemString(obj, "low", tmpObj);

    tmpObj = Py_BuildValue("f", val.highVal_);
    PyDict_SetItemString(obj, "high", tmpObj);

    std::string stat = Model::statString(val.stat_);

    tmpObj = Py_BuildValue("s", stat.c_str());
    PyDict_SetItemString(obj, "stat", tmpObj);
  }

  tmpObj = Py_BuildValue("s", val.units_.c_str());
  PyDict_SetItemString(obj, "units", tmpObj);
 
  return obj;
}

/**.......................................................................
 * Build a variate dictionary
 */
PyObject* buildVarDict(Model::VarVal& val, unsigned len, char** cPtr)
{
  PyObject* obj = PyDict_New();
  PyObject* tmpObj = 0;

  //------------------------------------------------------------
  // Build the val member of this dictionary, as a string or double,
  // depending on the variate type
  //------------------------------------------------------------

  if(val.isStr_) {
    tmpObj = Py_BuildValue("s", val.strVal_.c_str());
    PyDict_SetItemString(obj, "val", tmpObj);
  } else {
    tmpObj = Py_BuildValue("f", val.refVal_);
    PyDict_SetItemString(obj, "val", tmpObj);
  }

  //------------------------------------------------------------
  // If the component is fixed, just set the value, else the 1-sigma
  // confidence interval
  //------------------------------------------------------------

  tmpObj = Py_BuildValue("b", val.isFixed_);
  PyDict_SetItemString(obj, "fixed", tmpObj);
    
  if(!val.isFixed_) {
    tmpObj = Py_BuildValue("f", val.lowVal_);
    PyDict_SetItemString(obj, "low", tmpObj);

    tmpObj = Py_BuildValue("f", val.highVal_);
    PyDict_SetItemString(obj, "high", tmpObj);

    std::string stat = Model::statString(val.stat_);

    tmpObj = Py_BuildValue("s", stat.c_str());
    PyDict_SetItemString(obj, "stat", tmpObj);

    //------------------------------------------------------------
    // Now build a numpy array to contain the values in the chain
    //------------------------------------------------------------

#if 1
    tmpObj = getArray(len, DataType::DOUBLE, cPtr);
    PyDict_SetItemString(obj, "vals", tmpObj);
#endif
  }

  tmpObj = Py_BuildValue("s", val.units_.c_str());
  PyDict_SetItemString(obj, "units", tmpObj);
 
  return obj;
}

/**.......................................................................
 * SPYDOC
 *
 * Return a Climax run manager to the python environment.  
 *
 * Use like: 
 * 
 *    rm = climaxPyTest.newRunManager()
 *
 * Context:
 * 
 *   Part of the Climax external calling interface.  This function returns the
 *   run manager object that must be passed into all other external calling 
 *   interface functions;
 *
 *     parseFile()
 *     lnLikelihood()
 *     delRunManager()
 * 
 *   NB: when you are finished with this object, you can delete it by calling
 *
 *     delRunManager()
 *
 * EPYDOC
 */
static PyObject* newRunManager(PyObject* self, PyObject* args)
{
  gcp::util::RunManager* rm = new RunManager();
  PyObject* obj = PyCObject_FromVoidPtr(rm, 0);

  return obj;
}

/**.......................................................................
 * SPYDOC
 *
 * Delete a Climax run manager previously allocated by newRunManager()
 * 
 * Use like: 
 * 
 *    climaxPyTest.newRunManager(rm)
 *
 * Where:
 *
 *    rm  -  is the run manager previously allocated by newRunManager()
 *
 * Context:
 * 
 *    Part of the Climax external calling interface.  This function deletes the
 *    run manager object that was allocated by newRunManager()
 *
 * EPYDOC
 */
static PyObject* delRunManager(PyObject* self, PyObject* args)
{
  gcp::util::RunManager* rm = new RunManager();
  delete rm;
  rm = 0;
  return args;
}

/**.......................................................................
 * SPYDOC
 *
 * Load a Climax parameter file into a previously allocated run manager
 * 
 * Use like: 
 * 
 *    climaxPyTest.parseFile(rm, file)
 *
 * Where:
 *
 *    rm   - is the run manager previously allocated by newRunManager()
 *
 *    file - is the parameter file to load
 *
 * Context:
 * 
 *    Part of the Climax external calling interface.  This function loads a parameter
 *    file defining datasets and models into the run manager object that was previously 
 *    allocated by newRunManager()
 *
 * EPYDOC
 */
static PyObject* parseFile(PyObject* self, PyObject* args)
{
  try {
    PyObject* obj    = PyParser::getArrayItem(args, 0);
    RunManager* rm = (RunManager*)PyCObject_AsVoidPtr(obj);
    
    bool doLogo=true;
    if(PyParser::getSize(args) > 2) {
      PyParser logoParser(PyParser::getArrayItem(args, 2));
      doLogo = logoParser.getBoolVal();
    }

    if(doLogo) {
      Logo logo;
      logo.display();
    }
      
    PyParser fileParser(PyParser::getArrayItem(args, 1));
    rm->parseFile(fileParser.getString());
    rm->initializeForMarkovChain();

    return args;

  } catch(Exception& err) {
    COUT(err.what());
    PyObject* obj = Py_BuildValue("f", 0.0);
    return obj;
  }
}

/**.......................................................................
 * SPYDOC
 *
 * Return a dictionary of variable components and expected units
 * 
 * Use like: 
 * 
 *    d = climaxPyTest.getVariableComponents(rm)
 *
 * Where:
 *
 *    rm  -  is the run manager previously allocated by newRunManager()
 *
 * Context:
 *
 *    After a parameter file has been loaded into the run manager, this 
 *    function will return the length, order and units in which variable components
 *    should be specified when calling lnLikelihood() (as a numpy double array)
 *
 * EPYDOC
 */
static PyObject* getVariableComponents(PyObject* self, PyObject* args)
{
  PyObject* ret=0;

  try {
    PyObject* obj  = PyParser::getArrayItem(args, 0);
    RunManager* rm = (RunManager*)PyCObject_AsVoidPtr(obj);
    
    std::vector<Variate*>& vc = rm->getVariableComponents();
    ret = PyList_New(vc.size());
    
    for(unsigned i=0; i < vc.size(); i++) {
      PyObject* tmpObj = Py_BuildValue("ssd", vc[i]->name_.c_str(), vc[i]->units_.c_str(), vc[i]->getUnitVal());
      PyList_SET_ITEM(ret, i, tmpObj);
    }
  } catch(Exception& err) {
    ret = PyList_New(0);
  }

  return ret;
}

/**.......................................................................
 * SPYDOC
 *
 * Return the log likelihood of the passed array of parameter values
 * 
 * Use like: 
 * 
 *    lnlike = climaxPyTest.lnLikelihood(rm, arr)
 *
 * Where:
 *
 *    rm  -  is the run manager previously allocated by newRunManager()
 * 
 *    arr -  is a 1D numpy double array of current variable parameter values
 *           specified in the order and units returned by getVariableComponents()
 *
 * Context:
 *
 *    After a parameter file has been loaded into the run manager, this 
 *    function will return the log likelihood of the passed array of parameter values
 *
 * EPYDOC
 */
static PyObject* lnLikelihood(PyObject* self, PyObject* args)
{
  try {
    PyObject* obj  = PyParser::getArrayItem(args, 0);
    RunManager* rm = (RunManager*)PyCObject_AsVoidPtr(obj);

    PyObject* arr  = PyParser::getArrayItem(args, 1);

    unsigned n = PyParser::getNumpyLength(arr, 0);
    double* dptr = PyParser::getNumpyDoublePtr(arr, true);

    gcp::util::Vector<double> vec(n);

    for(unsigned i=0; i < n; i++)
      vec[i] = dptr[i];
    
    PyObject* ret = Py_BuildValue("f", rm->lnLikelihoodUnits(vec));

#if 0
    rm->printModel();
#endif

    return ret;

  } catch(Exception& err) {
    COUT(err.what());
    PyObject* obj = Py_BuildValue("f", 0.0);
    return obj;
  }
}

static PyObject* pgplotTest(PyObject* self, PyObject* args)
{
  cpgslct(1);
  PgUtil::advance();

  PyObject* obj = Py_BuildValue("f", 0.0);
  return obj;
}

#if DIR_HAVE_NUMPY
int numpyTypeOf(gcp::util::DataType::Type dataType)
{
  switch (dataType) {
  case DataType::BOOL:
    return PyArray_BOOL;
    break;
  case DataType::UCHAR:
    return PyArray_UBYTE;
    break;
  case DataType::CHAR:
    return PyArray_BYTE;
    break;
  case DataType::USHORT:
    return PyArray_USHORT;
    break;
  case DataType::SHORT:
    return PyArray_SHORT;
    break;
  case DataType::UINT:
    return PyArray_UINT;
    break;
  case DataType::INT:
    return PyArray_INT;
    break;
  case DataType::ULONG:
    return PyArray_ULONG;
    break;
  case DataType::LONG:
    return PyArray_LONG;
    break;
  case DataType::FLOAT:
    return PyArray_FLOAT;
    break;
  case DataType::DOUBLE:
    return PyArray_DOUBLE;
    break;
  case DataType::STRING:
    return PyArray_STRING;
    break;
  case DataType::DATE:
    return PyArray_DOUBLE;
    break;
  case DataType::COMPLEX_FLOAT:
    return PyArray_CFLOAT;
    break;
  case DataType::COMPLEX_DOUBLE:
    return PyArray_CDOUBLE;
    break;
  default:
    return PyArray_UBYTE;
    break;
  }
}

/**.......................................................................
 * Return the dimensions of this register block
 */
std::vector<npy_intp> getNumPyInds(std::vector<int>& dims)
{
  std::vector<npy_intp> npyInds(dims.size());

  for(unsigned iDim=0; iDim < dims.size(); iDim++)
    npyInds[iDim]  = 0;

  return npyInds;
}
#endif

PyObject* getArray(unsigned len, gcp::util::DataType::Type dataType, char** data)
{
  std::vector<int> dims(1);
  dims[0] = len;
  return getArray(dims, dataType, data);
}

PyObject* getArray(std::vector<int> dims, gcp::util::DataType::Type dataType, char** data)
{
  PyObject* obj = 0;

#if DIR_HAVE_NUMPY

  int type = numpyTypeOf(dataType);
  obj = PyArray_FromDims(dims.size(), &dims[0], type);

  std::vector<npy_intp> npyInds = getNumPyInds(dims);

  if(data != 0) 
    *data = (char*)PyArray_GetPtr((PyArrayObject*)obj, &npyInds[0]);

#else

  ThrowError("Numpy environment is not defined");

#endif

  return obj;
}
