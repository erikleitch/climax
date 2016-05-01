#include "gcp/util/GenericNormalization.h"
#include "gcp/util/Flux.h"
#include "gcp/util/String.h"
#include "gcp/util/Temperature.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
GenericNormalization::GenericNormalization() 
{
  allowUnitless(true);

  knownQuantities_.push_back(new Temperature());
  knownQuantities_.push_back(new Flux());
}

/**.......................................................................
 * Destructor.
 */
GenericNormalization::~GenericNormalization() 
{
  for(unsigned i=0; i < knownQuantities_.size(); i++) {
    delete knownQuantities_[i];
  }
}

/**.......................................................................
 * Set val
 */
void GenericNormalization::setVal(double val, std::string units)
{
  String unitStr(units);

  if(unitStr.isEmpty()) {
    VariableUnitQuantity::setVal(val, units);
  } else {
    for(unsigned i=0; i < knownQuantities_.size(); i++) {
      ConformableQuantity* cq = knownQuantities_[i];
      if(cq->isValidUnit(units)) {
	currentQuantity_ = cq;
	cq->setVal(val, units);
	break;
      }
    }

    if(currentQuantity_ == 0)
      ThrowColorError("Unrecognized unit '" << units << "'" << " for any quantity", "red");
  }

  val_   = val;
  units_ = units;

  return;
}

/**.......................................................................
 * Get val
 */
double GenericNormalization::getVal(std::string units)
{
  //------------------------------------------------------------
  // If the units are empty, check that our current units are also
  // empty.  
  //------------------------------------------------------------

  String unitStr(units);

  if(unitStr.isEmpty()) {

    if(currentQuantity_ != 0)
      ThrowColorError("Attempt to retrieve a unitless value from a unit'd quantity", "red");

    return VariableUnitQuantity::getUnitVal(units);

  } else {

    if(currentQuantity_ == 0)
      ThrowSimpleColorError(std::endl << "Attempt to retrieve a unit'd quantity from a unitless quantity", "red");

    if(currentQuantity_->isValidUnit(units)) 
      return currentQuantity_->getUnitVal(units);

    ThrowSimpleColorError(std::endl << "Unit '" << units << "'" << " is not a valid unit for a " << currentQuantity_->typeName_ << " object", "red");
  }

  return 0.0;
}


GenericNormalization::GenericNormalization(const GenericNormalization& gn)
{
  *this = (GenericNormalization&) gn;
}

GenericNormalization::GenericNormalization(GenericNormalization& gn)
{
  *this = gn;
}

void GenericNormalization::operator=(const GenericNormalization& gn)
{
  *this = (GenericNormalization&) gn;
}

void GenericNormalization::operator=(GenericNormalization& gn)
{
  val_   = gn.val_;
  units_ = gn.units_;

  currentQuantity_ = 0;

  setVal(val_, units_);
}
