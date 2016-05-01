#include "gcp/util/VariableUnitQuantity.h"
#include "gcp/util/String.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
VariableUnitQuantity::VariableUnitQuantity() 
{
  unitlessAllowed_ = false;
}

/**.......................................................................
 * Destructor.
 */
VariableUnitQuantity::~VariableUnitQuantity() {}

void VariableUnitQuantity::allowUnitless(bool allow)
{
  unitlessAllowed_ = allow;
}

bool VariableUnitQuantity::unitlessAllowed()
{
  return unitlessAllowed_;
}

void VariableUnitQuantity::setVal(double val, std::string units)
{
  if(stringIsEmpty(units) && unitlessAllowed_) {
    val_          = val;
    hasValue_     = true;
    wasSpecified_ = true;
  } else {
    ConformableQuantity::setVal(val, units);
  }
}

void VariableUnitQuantity::setVal(std::string val) 
{
  if(unitlessAllowed_) {
    String valStr(val);
    val_ = valStr.toDouble();
    hasValue_     = true;
    wasSpecified_ = true;
  } else {
    ConformableQuantity::setVal(val);
  }
}

double VariableUnitQuantity::getVal(double val, std::string units)
{
  if(stringIsEmpty(units) && unitlessAllowed_) {
    return val;
  } else {
    return ConformableQuantity::getVal(val, units);
  }
}

void VariableUnitQuantity::setUnits(std::string units)
{
  if(!(stringIsEmpty(units) && unitlessAllowed_)) {
    ConformableQuantity::setUnits(units);
  }
}

bool VariableUnitQuantity::stringIsEmpty(std::string& str)
{
  return (str.size()==0 || (str.size()==1 && str[0]=='\0') || str == " ");
}

void VariableUnitQuantity::operator=(const VariableUnitQuantity& vuq)
{
  operator=((VariableUnitQuantity&) vuq);
}

void VariableUnitQuantity::operator=(VariableUnitQuantity& vuq)
{
  val_             = vuq.val_;
  unitlessAllowed_ = vuq.unitlessAllowed_;
}
