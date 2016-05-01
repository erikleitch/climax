#include "gcp/util/ConformableQuantity.h"
#include "gcp/util/Exception.h"
#include "gcp/util/String.h"

#include <string.h>

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ConformableQuantity::ConformableQuantity() 
{
  typeName_ = "ConformableQuantity";
}

/**.......................................................................
 * Destructor.
 */
ConformableQuantity::~ConformableQuantity() {}

void ConformableQuantity::setVal(double val, std::string units)
{
  val_ = getVal(val, units);
  hasValue_ = true;
}

/**.......................................................................
 * setVal(string) defaults to the base-class stub unless overloaded by
 * inheritors
 */
void ConformableQuantity::setVal(std::string val)
{
  COUT("Inside CQ::setVal with val = '" << val << "'");
  Variate::setVal(val);
}

/**.......................................................................
 * Return the native value for the corresponding unit'ed type
 */
double ConformableQuantity::getVal(double val, std::string units)
{
  unsigned nmin;

  for(std::map<std::string, Conversion>::iterator iter=unitConversions_.begin();
      iter != unitConversions_.end(); iter++) {

    std::string str  = iter->first;
    Conversion  conv = iter->second;

    if(strcasecmp(str.c_str(), units.c_str())==0) 
      return val * conv.scaleFactor_ + conv.offset_;
 
  }
      
  if(units.length() == 1 && !isprint(units[0]))
    units = " ";

  {
    XtermManip xtm;
    ThrowSimpleError(COLORIZE(xtm, "red", std::endl << "Unrecognized unit: '" << units << "' for variate " << name_ 
			      << std::endl << std::endl)
		     << COLORIZE(xtm, "green", "Recognized units are: " << std::endl << std::endl
				 << getUnitsString()));
  }

  return 0.0;
}

/**.......................................................................
 * Set the internal multiplier to convert from native units to
 * user-specified units
 */
void ConformableQuantity::setUnits(std::string units)
{
  unsigned nmin;

  //------------------------------------------------------------
  // See if this units string matches any known units
  //------------------------------------------------------------

  for(std::map<std::string, Conversion>::iterator iter=unitConversions_.begin();
      iter != unitConversions_.end(); iter++) {

    std::string str  = iter->first;
    Conversion  conv = iter->second;

    if(strcasecmp(str.c_str(), units.c_str())==0) {

      unitsToNative_ = conv;
      units_         = str;

      return;
    }
  }

  //------------------------------------------------------------
  // If not, see if this is a scaling derivable from units that we
  // know about
  //------------------------------------------------------------

  String unitsStr(units);
  for(std::map<std::string, Conversion>::iterator iter=unitConversions_.begin();
      iter != unitConversions_.end(); iter++) {

    if(unitsStr.contains(iter->first)) {

      String unitsCopy(units);
      unitsCopy.replace(iter->first, "");
      unitsCopy.strip(' ');

      double scale = 1.0;

      try {

	scale = unitsCopy.toDouble();
	
	// This means that the remainder of the string is a valid
	// scale factor.  Insert a new conversion and allow the units

	Conversion  conv = iter->second;
	
	unitsToNative_ = addConversion(units, conv.scaleFactor_ * scale, conv.offset_);
	units_         = units;

	return;
      } catch(Exception& err) {
      } catch(...) {
      }

    }
  }

  if(units.size() == 0 || units[0] == '\0') {
    ThrowColorError(std::endl << "No unit specified." << std::endl
		    << "Recognized units are: " << std::endl << std::endl
		    << getUnitsString(), "red");
  } else {
    ThrowColorError(std::endl << "Unrecognized unit: '" << units << "'" << std::endl << std::endl
		    << "Recognized units are: " << std::endl << std::endl
		    << getUnitsString(), "red");
  }

  return;
}

/**.......................................................................
 * Add a unit conversion to the list of conversions we know about
 */
Variate::Conversion ConformableQuantity::addConversion(std::string units, double scaleFactor, double offset)
{
  unitConversions_[units].scaleFactor_ = scaleFactor;
  unitConversions_[units].offset_      = offset;

  return unitConversions_[units];
}

std::string ConformableQuantity::getUnitsString()
{
  std::ostringstream os;

  for(std::map<std::string, Conversion>::iterator iter=unitConversions_.begin();
      iter != unitConversions_.end(); iter++) {

    std::string str = iter->first;
    os << "    " << str << std::endl;
  }

  return os.str();
}

void ConformableQuantity::printUnits()
{
  COUT("Recognized units are: ");

  for(std::map<std::string, Conversion>::iterator iter=unitConversions_.begin();
      iter != unitConversions_.end(); iter++) {

    std::string str = iter->first;
    COUT(str);
  }
}

bool ConformableQuantity::isValidUnit(std::string units)
{
  //------------------------------------------------------------
  // See if this units string matches any known units
  //------------------------------------------------------------

  for(std::map<std::string, Conversion>::iterator iter=unitConversions_.begin();
      iter != unitConversions_.end(); iter++) {

    std::string str  = iter->first;
    Conversion  conv = iter->second;

    if(strcasecmp(str.c_str(), units.c_str())==0) {
      return true;
    }
  }

  //------------------------------------------------------------
  // If not, see if this is a scaling derivable from units that we
  // know about
  //------------------------------------------------------------

  String unitsStr(units);
  for(std::map<std::string, Conversion>::iterator iter=unitConversions_.begin();
      iter != unitConversions_.end(); iter++) {

    if(unitsStr.contains(iter->first)) {

      String unitsCopy(units);
      unitsCopy.replace(iter->first, "");
      unitsCopy.strip(' ');

      double scale = 1.0;

      try {
	scale = unitsCopy.toDouble();
	return true;
      } catch(Exception& err) {
	return false;
      } catch(...) {
	return false;
      }
    }
  }

  return false;
}


double ConformableQuantity::getUnitVal(std::string units)
{
  Conversion* conv = findConversion(units);

  if(conv == 0)
    ThrowColorError("No valid conversion exists for units '" << units << "'", "red");

  return (val_ - conv->offset_) / conv->scaleFactor_;
}

Variate::Conversion* ConformableQuantity::findConversion(std::string units)
{
  //------------------------------------------------------------
  // See if this units string matches any known units
  //------------------------------------------------------------

  for(std::map<std::string, Conversion>::iterator iter=unitConversions_.begin();
      iter != unitConversions_.end(); iter++) {

    std::string str  = iter->first;
    Conversion*  conv = &iter->second;

    if(strcasecmp(str.c_str(), units.c_str())==0) {
      return conv;
    }
  }

  //------------------------------------------------------------
  // If not, see if this is a scaling derivable from units that we
  // know about
  //------------------------------------------------------------

  String unitsStr(units);
  for(std::map<std::string, Conversion>::iterator iter=unitConversions_.begin();
      iter != unitConversions_.end(); iter++) {

    if(unitsStr.contains(iter->first)) {

      String unitsCopy(units);
      unitsCopy.replace(iter->first, "");
      unitsCopy.strip(' ');

      double scale = 1.0;

      try {
	scale = unitsCopy.toDouble();
	
	// This means that the remainder of the string is a valid
	// scale factor.  Insert a new conversion and allow the units

	Conversion  conv = iter->second;
	addConversion(units, conv.scaleFactor_ * scale, conv.offset_);
	return &unitConversions_[units];

      } catch(Exception& err) {
	COUT(err.what());
      } catch(...) {
      }

    }
  }

  return 0;
}

void ConformableQuantity::printConversions()
{
  for(std::map<std::string, Conversion>::iterator iter=unitConversions_.begin();
      iter != unitConversions_.end(); iter++) {
    std::string str  = iter->first;
    COUT("Found conversion " << str.c_str());
  }
}

#if 0
ConformableQuantity::ConformableQuantity(const ConformableQuantity& cq)
{
  *this = cq;
}

ConformableQuantity::ConformableQuantity(ConformableQuantity& cq)
{
  *this = cq;
}

void ConformableQuantity::operator=(const ConformableQuantity& cq)
{
  *this = (ConformableQuantity&) cq;
}

void ConformableQuantity::operator=(ConformableQuantity& cq)
{
  finite_ = cq.finite_;
  typeName_ = cq.typeName_;
  unitConversions_ = cq.unitConversions_;
}
#endif
