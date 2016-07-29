#include "gcp/util/Exception.h"
#include "gcp/util/ParameterManager.h"
#include "gcp/util/String.h"

#include <iomanip>

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ParameterManager::ParameterManager() {}

/**.......................................................................
 * Destructor.
 */
ParameterManager::~ParameterManager() 
{
  for(unsigned i=0; i < parameterVec_.size(); i++) {
    delete parameterVec_[i];
    parameterVec_[i] = 0;
  }
}

/**.......................................................................
 * Add a parameter to our map.  If resizable=true, we create a new
 * vector of the specified type.
 */
void ParameterManager::remParameter(std::string name)
{
  if(parameterMap_.find(name) == parameterMap_.end()) {
    COUTCOLOR("No parameter found with name: " << name, "red");
    return;
  }

  std::map<std::string, Parameter*>::iterator pIter = parameterMap_.find(name);
  Parameter* param = pIter->second;

  parameterMap_.erase(pIter);

  std::map<Parameter*, std::string>::iterator rIter = reverseParameterMap_.find(param);
  reverseParameterMap_.erase(rIter);

  std::vector<Parameter*>::iterator vIter = parameterVec_.begin();
  for(; vIter != parameterVec_.end(); vIter++) {
    Parameter* curr = *vIter;
    if(curr == param)
      break;
  }

  if(vIter != parameterVec_.end()) {
    parameterVec_.erase(vIter);
    delete param;
  }
}

/**.......................................................................
 * Add a parameter to our map.  If resizable=true, we create a new
 * vector of the specified type.
 */
void ParameterManager::addParameter(std::string name, DataType::Type type, std::string comment, bool resizable)
{
  if(parameterMap_.find(name) != parameterMap_.end()) {
    ThrowColorError("Another parameter with name: " << name << " already exists for this dataset", "red");
  }

  Parameter* param = new Parameter();

  param->data_.type_ = type;
  param->data_.allowAssignmentFromDifferentType(false);

  param->setResizable(resizable);
  param->setComment(comment);
  param->owner_ = this;

  parameterMap_[name] = param;
  reverseParameterMap_[param] = name;

  parameterVec_.push_back(param);
}

/**.......................................................................
 * Add another parameter manager's parameters to our hierarchy of
 * parameters (with a name prefix)
 */
void ParameterManager::addParameter(std::string name, ParameterManager& pm)
{
  std::ostringstream os;
  for(std::map<std::string, Parameter*>::iterator iter=pm.parameterMap_.begin(); 
      iter != pm.parameterMap_.end(); iter++) {
    os.str("");
    os << name << "." << iter->first;

    std::string parName = os.str();
    if(parameterMap_.find(parName) != parameterMap_.end())
      ThrowColorError("Another parameter with name: " << parName << " already exists for this object ", "red");

    parameterMap_[os.str()] = iter->second;
    reverseParameterMap_[iter->second] = os.str();
  }
}

/**.......................................................................
 * Add another parameter manager's parameters to our hierarchy of
 * parameters (without a prefix)
 */
void ParameterManager::addParameter(ParameterManager& pm)
{
  std::ostringstream os;
  for(std::map<std::string, Parameter*>::iterator iter=pm.parameterMap_.begin(); 
      iter != pm.parameterMap_.end(); iter++) {
    os.str("");
    os << iter->first;

    std::string parName = os.str();
    if(parameterMap_.find(parName) != parameterMap_.end())
      ThrowColorError("Another parameter with name: " << parName << " already exists for object " << name_, "red");

    parameterMap_[os.str()] = iter->second;
    reverseParameterMap_[iter->second] = os.str();
  }
}

/**.......................................................................
 * Update another parameter manager's parameters to our hierarchy of
 * parameters
 */
void ParameterManager::updateParameter(std::string name, ParameterManager& pm)
{
  std::ostringstream os;
  for(std::map<std::string, Parameter*>::iterator iter=pm.parameterMap_.begin(); 
      iter != pm.parameterMap_.end(); iter++) {
    os.str("");
    os << name << "." << iter->first;

    if(parameterMap_.find(os.str()) == parameterMap_.end()) {
      *parameterMap_[os.str()] = *iter->second;
    }
  }
}

//------------------------------------------------------------
// Accessor methods for parameters
//------------------------------------------------------------

std::string ParameterManager::getStringVal(std::string name)
{
  return getParameter(name, true)->data_.getStringVal();
}

bool ParameterManager::getBoolVal(std::string name)
{
  return getParameter(name, true)->data_.getBoolVal();
}

unsigned char ParameterManager::getUcharVal(std::string name)
{
  return getParameter(name, true)->data_.getUcharVal();
}

char ParameterManager::getCharVal(std::string name)
{
  return getParameter(name, true)->data_.getCharVal();
}

unsigned short ParameterManager::getUshortVal(std::string name)
{
  return getParameter(name, true)->data_.getUshortVal();
}

short ParameterManager::getShortVal(std::string name)
{
  return getParameter(name, true)->data_.getShortVal();
}

unsigned int ParameterManager::getUintVal(std::string name)
{
  return getParameter(name, true)->data_.getUintVal();
}

int ParameterManager::getIntVal(std::string name)
{
  return getParameter(name, true)->data_.getIntVal();
}

unsigned long ParameterManager::getUlongVal(std::string name)
{
  return getParameter(name, true)->data_.getUlongVal();
}

long ParameterManager::getLongVal(std::string name)
{
  return getParameter(name, true)->data_.getLongVal();
}

float ParameterManager::getFloatVal(std::string name)
{
  return getParameter(name, true)->data_.getFloatVal();
}

double ParameterManager::getDoubleVal(std::string name)
{
  return getParameter(name, true)->data_.getDoubleVal();
}

ParameterManager::Parameter* ParameterManager::getParameter(std::string name, bool checkVal)
{
  if(parameterMap_.find(name) == parameterMap_.end()) {

    std::ostringstream os;
    listParameters(os);

    XtermManip xtm;

    std::ostringstream os2;
    os2.str("");
    os2 << COLORIZE(xtm, "red", std::endl << "No parameter '" << name << "' found for object '" << name_ << "'" << std::endl << std::endl);
    os2 << COLORIZE(xtm, "green", "Recognized parameters for " << name_ << " are: " << std::endl << std::endl << os.str());

    ThrowSimpleError(os2.str());
  }

  Parameter* param = parameterMap_.find(name)->second;

  if(checkVal && !param->data_.hasValue()) {
    ThrowSimpleColorError(std::endl << "No value has been specified for parameter: " << name_ << "." << name, "red");
  }

  return param;
}

/**.......................................................................
 * Set the value of a parameter
 */
void ParameterManager::setParameter(std::string name, std::string val, std::string units, bool external)
{
  Parameter* param = getParameter(name, false);
  DataType& dat = param->data_;

  String valStr(val);

  switch(dat.type_) {
  case DataType::INT:
    dat = valStr.toInt();
    break;
  case DataType::UINT:
    dat = (unsigned int)valStr.toInt();
    break;
  case DataType::FLOAT:
    dat = valStr.toFloat();
    break;
  case DataType::DOUBLE:
    dat = valStr.toDouble();
    break;
  case DataType::STRING:
    dat = valStr.str();
    break;
  case DataType::BOOL:
    dat = valStr.toBool();
    break;
  default:
    ThrowColorError("Unsupported conversion from string to type: " << dat.type_, "red");
    break;
  }

  // Set the units too

  if(units == "")
    param->units_ = " ";
  else
    param->units_ = units;
}

/**.......................................................................
 * Increment the value of a parameter
 */
void ParameterManager::incrementParameter(std::string name, std::string val, std::string units, bool external)
{
  Parameter* param = getParameter(name, false);
  DataType& dat = param->data_;

  String valStr(val);

  switch(dat.type_) {
  case DataType::STRING:
    {
      if(dat.hasValue_) {
	std::ostringstream os;
	os << dat.data_.str << " " << valStr.str();
	dat = os.str();
      } else {
	dat = valStr.str();
      }
    }
    break;
  default:
    ThrowColorError("Unsupported increment operator on parameter of type: " << dat.type_, "red");
    break;
  }

  // Set the units too

  param->units_ = units;
}

/**.......................................................................
 * List parameters managed by this object
 */
void ParameterManager::listParameters(std::ostringstream& os, bool sort)
{
  int maxLen = 0, size;

  for(std::map<std::string, Parameter*>::iterator iter = parameterMap_.begin(); iter != parameterMap_.end(); iter++) {
    size = iter->first.size();
    maxLen = (maxLen > size) ? maxLen : size;
  }

  unsigned lineWidth = XtermManip::getNCol() - (maxLen+2);
  
  if(sort) {
    for(std::map<std::string, Parameter*>::iterator iter = parameterMap_.begin(); iter != parameterMap_.end(); iter++)
      formatParameter(os, iter->first, iter->second->comment_, maxLen, lineWidth);
  } else {
    for(unsigned i=0; i < parameterVec_.size(); i++)
      formatParameter(os, reverseParameterMap_[parameterVec_[i]], parameterVec_[i]->comment_, maxLen, lineWidth);
  }
}

/**.......................................................................
 * Reformat the comment string to line up with the end of the label
 * column
 */
void ParameterManager::formatParameter(std::ostringstream& os, std::string name, std::string comment,
				       unsigned maxLen, unsigned lineWidth)
{
  String commentStr(comment);
  XtermManip xtm;
  std::vector<std::string> lines;

  //------------------------------------------------------------
  // Construct a string we will insert every time we encounter an
  // explicit newline
  //------------------------------------------------------------

  std::ostringstream newlineOs, osLine, osWord;

  for(unsigned i=0; i < maxLen+2; i++)
    newlineOs << " ";

  unsigned nInCurrentLine=0, nInCurrentWord=0;

  osLine.str("");
  osWord.str("");

  //------------------------------------------------------------
  // Iterate through this comment string, accumulating characters into
  // words, adding new lines as we exceed the allowed line width
  //------------------------------------------------------------

  for(unsigned i=0; i < commentStr.size(); i++) {

    char c = commentStr[i];
    
    //------------------------------------------------------------
    // We will start a new word whever we hit a space
    //------------------------------------------------------------

    if(c == ' ') {
      appendLineAndWord(osLine, nInCurrentLine, osWord, nInCurrentWord, newlineOs, lines, lineWidth, false);

      //------------------------------------------------------------
      // We will start a new word and a new line when we hit an explicit
      // line break
      //------------------------------------------------------------
    
    } else if(c == '\n') {
      appendLineAndWord(osLine, nInCurrentLine, osWord, nInCurrentWord, newlineOs, lines, lineWidth, true);

      //------------------------------------------------------------
      // Else we are just appending to the current word
      //------------------------------------------------------------
      
    } else {
      osWord << c;
      ++nInCurrentWord;
    }
  }

  //------------------------------------------------------------
  // Finish the last line and word we may have been accumulating
  //------------------------------------------------------------

  appendLineAndWord(osLine, nInCurrentLine, osWord, nInCurrentWord, newlineOs, lines, lineWidth, true);

  os << COLORIZE(xtm, "yellow", std::left << std::setw(maxLen+2) << name);
  
  for(unsigned iLine = 0; iLine < lines.size(); iLine++)
    os << COLORIZE(xtm, "green", lines[iLine]) << std::endl;
}

/**.......................................................................
 * Append a word and line to the passed array
 */
void ParameterManager::appendLineAndWord(std::ostringstream& osLine, unsigned& nInLine, 
					 std::ostringstream& osWord, unsigned& nInWord, 
					 std::ostringstream& filler, std::vector<std::string>& lines,
					 unsigned lineWidth, bool atLineBreak)
{
  //------------------------------------------------------------ 
  // If there's room in the line for this word, just append it to the
  // current line
  //------------------------------------------------------------
  
  if(nInLine + nInWord < lineWidth) {
    osLine << osWord.str() << " ";
    nInLine += nInWord + 1;

    //------------------------------------------------------------
    // If not, append the line, and add the word to the new line, with
    // a space
    //------------------------------------------------------------
    
  } else {
    lines.push_back(osLine.str());

    osLine.str("");
    osLine << filler.str() << osWord.str() << " ";
    nInLine = nInWord + 1;
  }

  //------------------------------------------------------------
  // If we were at a line break, add a new line now
  //------------------------------------------------------------
  
  if(atLineBreak) {
    lines.push_back(osLine.str());
    osLine.str("");
    osLine << filler.str();
    nInLine = 0;
  }
  
  //------------------------------------------------------------
  // Initialize the word count
  //------------------------------------------------------------
  
  osWord.str("");
  nInWord = 0;
}

void ParameterManager::Parameter::setComment(std::string comment)
{
  comment_ = comment;
}

void ParameterManager::Parameter::setResizable(bool resizable)
{
  if(!resizable)
    return;

  resizable_      = true;
  data_.isArray_  = true;
  data_.hasValue_ = false;

  switch (data_.type_) {
  case DataType::BOOL:
    data_.ptr_ = new std::vector<bool>(1);
    break;
  case DataType::UCHAR:
    data_.ptr_ = new std::vector<unsigned char>(1);
    break;
  case DataType::CHAR:
    data_.ptr_ = new std::vector<char>(1);
    break;
  case DataType::USHORT:
    data_.ptr_ = new std::vector<unsigned short>(1);
    break;
  case DataType::SHORT:
    data_.ptr_ = new std::vector<short>(1);
    break;
  case DataType::INT:
    data_.ptr_ = new std::vector<int>(1);
    break;
  case DataType::UINT:
    data_.ptr_ = new std::vector<unsigned int>(1);
    break;
  case DataType::LONG:
    data_.ptr_ = new std::vector<long>(1);
    break;
  case DataType::ULONG:
    data_.ptr_ = new std::vector<unsigned long>(1);
    break;
  case DataType::FLOAT:
    data_.ptr_ = new std::vector<float>(1);
    break;
  case DataType::DOUBLE:
    data_.ptr_ = new std::vector<double>(1);
    break;
  default:
    ThrowError("Unsupported array type: " << data_.type_);
    break;
  }
}

void ParameterManager::Parameter::deleteResizable()
{
  if(!resizable_)
    return;

  switch (data_.type_) {
  case DataType::BOOL:
    delete (std::vector<bool>*) data_.ptr_;
    data_.ptr_ = 0;
    break;
  case DataType::UCHAR:
    delete (std::vector<unsigned char>*) data_.ptr_;
    data_.ptr_ = 0;
    break;
  case DataType::CHAR:
    delete (std::vector<char>*) data_.ptr_;
    data_.ptr_ = 0;
    break;
  case DataType::USHORT:
    delete (std::vector<unsigned short>*) data_.ptr_;
    data_.ptr_ = 0;
    break;
  case DataType::SHORT:
    delete (std::vector<short>*) data_.ptr_;
    data_.ptr_ = 0;
    break;
  case DataType::INT:
    delete (std::vector<int>*) data_.ptr_;
    data_.ptr_ = 0;
    break;
  case DataType::UINT:
    delete (std::vector<unsigned int>*) data_.ptr_;
    data_.ptr_ = 0;
    break;
  case DataType::LONG:
    delete (std::vector<long>*) data_.ptr_;
    data_.ptr_ = 0;
    break;
  case DataType::ULONG:
    delete (std::vector<unsigned long>*) data_.ptr_;
    data_.ptr_ = 0;
    break;
  case DataType::FLOAT:
    delete (std::vector<float>*) data_.ptr_;
    data_.ptr_ = 0;
    break;
  case DataType::DOUBLE:
    delete (std::vector<double>*) data_.ptr_;
    data_.ptr_ = 0;
    break;
  default:
    ThrowError("Unsupported array type: " << data_.type_);
    break;
  }
}

/**.......................................................................
 * Parse a string that represents a range, list or single value
 */
std::vector<double> ParameterManager::parseRange(String& val)
{
  std::vector<double> vals;

  //------------------------------------------------------------
  // Check for range specification
  //------------------------------------------------------------
  
  if(val.contains(":")) {
    String startStr = val.findNextInstanceOf("", false, ":", true, true);
    String stopStr  = val.findNextInstanceOf("", false, ":", true, true);
    String deltaStr = val.remainder();

    double startVal = startStr.toDouble();
    double stopVal  = stopStr.toDouble();
    double deltaVal = deltaStr.toDouble();
    
    unsigned nVal = (unsigned)((stopVal - startVal)/deltaVal + 1);
    vals.resize(nVal);
    
    for(unsigned iVal=0; iVal < nVal; iVal++) {
      vals[iVal] = startVal + iVal*deltaVal;
    }
    
    //------------------------------------------------------------
    // Else a list
    //------------------------------------------------------------
    
  } else if(val.contains(",")) {
    
    String nextVal;
    
    do {
      nextVal = val.findNextStringSeparatedByChars(",", true);
      
      if(!nextVal.isEmpty())
	vals.push_back(nextVal.toDouble());
      
    } while(!nextVal.isEmpty());
    
    // Else a single value
    
  } else {
    vals.resize(1);
    vals[0] = val.toDouble();
  }
  
  return vals;
}

/**.......................................................................
 * Copy common parameters from another parameter manager
 */ 
void ParameterManager::copyParameters(ParameterManager* pm, std::map<std::string, std::string>& excludedParameters, bool all)
{
  //------------------------------------------------------------
  // Iterate over the other manager's list of parameters
  //------------------------------------------------------------

  for(std::map<std::string, Parameter*>::iterator iter=pm->parameterMap_.begin(); iter != pm->parameterMap_.end(); ++iter) {
    std::string name = iter->first;
    Parameter* par   = iter->second;

    //------------------------------------------------------------
    // If this parameter is found in our map, copy its value
    //------------------------------------------------------------

    if((parameterMap_.find(name) != parameterMap_.end()) && (excludedParameters.find(name) == excludedParameters.end())) {
      Parameter* ourPar = getParameter(name);

      if(!ourPar->data_.hasValue() || all) {
	*ourPar = *par;
      }
    }
  }
}

/**.......................................................................
 * Set parameter values from another parameter manager
 */ 
void ParameterManager::setParameters(ParameterManager* pm, std::map<std::string, std::string>& excludedParameters, bool all)
{
  //------------------------------------------------------------
  // Iterate over the other manager's list of parameters
  //------------------------------------------------------------

  for(unsigned i=0; i < pm->parameterVec_.size(); i++) {
    Parameter* par = pm->parameterVec_[i];
    std::string name = pm->reverseParameterMap_[par];
    
    COUT("PM: setParameters pm = " << pm << " this = " << this << " name = " << name);
    //------------------------------------------------------------
    // If this parameter is found in our map, set its value
    //------------------------------------------------------------

    std::ostringstream os;
    if((parameterMap_.find(name) != parameterMap_.end()) && (excludedParameters.find(name) == excludedParameters.end())) {
        if(par->data_.hasValue()) {
            Parameter* ourPar = getParameter(name);
            os.str("");
            os << par->data_;
            setParameter(name, os.str(), par->units_);
        }
    }
  }
}

/**.......................................................................
 * Return true if this name matches any parameter (minimum matching allowed)
 */
bool ParameterManager::matches(std::string name)
{
  String nameStr(name);
  nameStr = nameStr.toLower();

  for(std::map<std::string, Parameter*>::iterator iter=parameterMap_.begin(); iter != parameterMap_.end(); ++iter) {
    String parStr(iter->first);
    parStr = parStr.toLower();

    if(parStr.contains(nameStr.str()))
      return true;
  }

  return false;
}

/**.......................................................................
 * Return true if this name unique matches any parameter (minimum matching allowed)
 */
bool ParameterManager::unique(std::string name)
{
  String nameStr(name);
  nameStr = nameStr.toLower();
  unsigned nMatch=0;
  nameStr.strip(' ');

  for(std::map<std::string, Parameter*>::iterator iter=parameterMap_.begin(); iter != parameterMap_.end(); ++iter) {
    String parStr(iter->first);
    parStr = parStr.toLower();
    parStr.strip(' ');

#if 0
    if(parStr.contains(nameStr.str()))
      ++nMatch;
#else
    if(parStr.str() == nameStr.str())
      ++nMatch;
#endif
  }

  if(nMatch == 0)
    ThrowSimpleError("No parameter matches: " << name);

  return nMatch==1;
}

/**.......................................................................
 * Return the parameter that matches this name
 */
std::string ParameterManager::getMatch(std::string name)
{
  if(!unique(name))
    ThrowSimpleError("'" << name << "' matches more than one parameter");

  String nameStr(name);
  nameStr = nameStr.toLower();
  nameStr.strip(' ');

  for(std::map<std::string, Parameter*>::iterator iter=parameterMap_.begin(); iter != parameterMap_.end(); ++iter) {
    String parStr(iter->first);
    parStr = parStr.toLower();
    parStr.strip(' ');

#if 0
    if(parStr.contains(nameStr.str()))
      return name;
#else
    if(parStr.str() == nameStr.str())
      return name;
#endif
  }

  return "";
}
