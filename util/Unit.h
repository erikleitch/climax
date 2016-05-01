// $Id: Unit.h,v 1.3 2012/05/08 21:58:00 eml Exp $

#ifndef GCP_UTIL_UNIT_H
#define GCP_UTIL_UNIT_H

/**
 * @file Unit.h
 * 
 * Tagged: Sun Oct 22 19:22:27 PDT 2006
 * 
 * @version: $Revision: 1.3 $, $Date: 2012/05/08 21:58:00 $
 * 
 * @author Erik Leitch.
 */
#include <string>
#include <vector>

namespace gcp {
  namespace util {

    // This class is intended as a base-class for all unit classes
    // defined in inheritors of ConformableQuantity class

    class Unit {
    public:

      //------------------------------------------------------------
      // Possible input units that we know about
      //------------------------------------------------------------
      
      enum Units {
	UNITS_UNKNOWN,
	UNITS_NONE,
	UNITS_JY,
	UNITS_MILLIJY,
	UNITS_JYSR,
	UNITS_MILLIJYSR,
	UNITS_MEGAJYSR,
	UNITS_JYBEAM,
	UNITS_UK,
	UNITS_MILLIK,
	UNITS_K,
	UNITS_Y,
	UNITS_SNR,
	UNITS_COUNTS,
	UNITS_METERS,
	UNITS_FEET,
	UNITS_KFEET,
	UNITS_COUNT_FLUX
      };

      /**
       * Constructor.
       */
      Unit();

      /**
       * Destructor.
       */
      virtual ~Unit();

      // Return true if the passed name is a recognized name for this
      // unit

      bool isThisUnit(std::string unitName);

      static Unit::Units stringToUnits(std::string);
      static std::string unitsToString(Unit::Units units);

    protected:

      // A collection of valid name for this unit

      std::vector<std::string> names_;
	
      // Add a name to this list

      void addName(std::string name);

      virtual void addNames();

    }; // End class Unit

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_UNIT_H
