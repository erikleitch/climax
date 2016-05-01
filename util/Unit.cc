#include "gcp/util/Unit.h"
#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/String.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Unit::Unit() {}

/**.......................................................................
 * Destructor.
 */
Unit::~Unit() {}

/**.......................................................................
 * Return true if the passed name is a recognized name for this
 * unit
 */
bool Unit::isThisUnit(std::string unitName)
{
  addNames();

  for(unsigned iName=0; iName < names_.size(); iName++) {
    if(unitName == names_[iName])
      return true;
  }
  return false;
}

/**.......................................................................
 * Add a name to this list
 */
void Unit::addName(std::string name)
{
  names_.push_back(name);
}

void Unit::addNames()
{
}

Unit::Units Unit::stringToUnits(std::string units)
{
  String unitStr(units);
  unitStr.strip(" ");
  unitStr = unitStr.toLower();

  if(unitStr == "jy") {
    return Unit::UNITS_JY;
  } else if(unitStr == "none") {
    return Unit::UNITS_NONE;
  } else if(unitStr == "mjy") {
    return Unit::UNITS_MILLIJY;
  } else if(unitStr == "jy/beam" || unitStr == "jy/bm") {
    return Unit::UNITS_JYBEAM;
  } else if(unitStr == "jy/sr") {
    return Unit::UNITS_JYSR;
  } else if(unitStr == "mjy/sr") {
    String origStr(units);
    if(origStr.contains("mJy"))
      return Unit::UNITS_MILLIJYSR;
    else if(origStr.contains("MJy"))
      return Unit::UNITS_MEGAJYSR;
    else {
      // Assume MJy/sr for now
      return Unit::UNITS_MEGAJYSR;
    //      ThrowSimpleColorError("Ambiguous unit: " << origStr, "red");
    }
  } else if(unitStr == "micro-kelvin") {
    return Unit::UNITS_UK;
  }  else if(unitStr == "muk") {
    return Unit::UNITS_UK;
  } else if(unitStr == "mk") {
    return Unit::UNITS_MILLIK;
  } else if(unitStr == "kelvin" || unitStr == "k") {
    return Unit::UNITS_K;
  } else if(unitStr == "comptony") {
    return Unit::UNITS_Y;
  } else if(unitStr == "snr") {
    return Unit::UNITS_SNR;
  } else if(unitStr == "counts") {
    return Unit::UNITS_COUNTS;
  } else if(unitStr == "meters") {
    return Unit::UNITS_METERS;
  } else if(unitStr == "feet") {
    return Unit::UNITS_FEET;
  } else if(unitStr == "kilo-feet") {
    return Unit::UNITS_KFEET;
  } else if(unitStr == "/cm**2/s") {
    return Unit::UNITS_COUNT_FLUX;
  } else {
    return Unit::UNITS_UNKNOWN;
  }
}

std::string Unit::unitsToString(Unit::Units units)
{
  if(units == Unit::UNITS_JY) {
    return "jy";
  } else if(units == Unit::UNITS_NONE) {
    return "none";
  } else if(units == Unit::UNITS_MILLIJY) {
    return "mjy";
  } else if(units == Unit::UNITS_JYBEAM) {
    return "jy/beam";
  } else if(units == Unit::UNITS_JYSR) {
    return "jy/sr";
  } else if(units == Unit::UNITS_MILLIJYSR) {
    return "mJy/sr";
  } else if(units == Unit::UNITS_MEGAJYSR) {
    return "MJy/sr";
  } else if(units == Unit::UNITS_UK) {
    return "micro-kelvin";
  } else if(units == Unit::UNITS_UK) {
    return "muk";
  } else if(units == Unit::UNITS_MILLIK) {
    return "mk";
  } else if(units == Unit::UNITS_K) {
    return "kelvin";
  } else if(units == Unit::UNITS_Y) {
    return "comptony";
  } else if(units == Unit::UNITS_SNR) {
    return "snr";
  } else if(units == Unit::UNITS_COUNTS) {
    return "counts";
  } else if(units == Unit::UNITS_METERS) {
    return "meters";
  } else if(units == Unit::UNITS_FEET) {
    return "feet";
  } else if(units == Unit::UNITS_KFEET) {
    return "kilo-feet";
  } else if(units == Unit::UNITS_COUNT_FLUX) {
    return "counts/cm**2/s";
  } else {
    return "unknown";
  }
}
