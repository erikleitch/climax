// $Id: Variate.h,v 1.4 2012/05/31 16:22:56 eml Exp $

#ifndef GCP_UTIL_VARIATE_H
#define GCP_UTIL_VARIATE_H

/**
 * @file Variate.h
 * 
 * Tagged: Fri Nov 19 13:57:31 PST 2010
 * 
 * @version: $Revision: 1.4 $, $Date: 2012/05/31 16:22:56 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Distribution.h"

#include <list>
#include <map>

#define VAR_DERIVE_FN(fn) void (fn)(void* args)

namespace gcp {
  namespace util {

    // Class intended as a base-class for variates of any type

    class Variate {
    public:

      // Struct to handle simple unit conversions of the form:
      //
      //  native = units * scaleFactor_ + offset_

      struct Conversion {
	double scaleFactor_;
	double offset_;

	Conversion() {
	  scaleFactor_ = 1.0;
	  offset_      = 0.0;
	}

	Conversion(const Conversion& conv) {
	  *this = conv;
	}

	Conversion(Conversion& conv) {
	  *this = conv;
	}

	void operator=(const Conversion& conv) {
	  operator=((Conversion&)conv);
	}

	void operator=(Conversion& conv) {
	  scaleFactor_ = conv.scaleFactor_;
	  offset_      = conv.offset_;
	}

      };

      /**
       * Constructor.
       */
      Variate();

      /**
       * Destructor.
       */
      virtual ~Variate();

      //------------------------------------------------------------
      // Statistics of the current value of this variate, for the
      // sampling distribution
      //------------------------------------------------------------

      virtual Probability samplingPdf();
      virtual Probability samplingCdf();
      Probability samplingPte();

      virtual Probability pdf();
      virtual Probability cdf();
      Probability pte();

      //------------------------------------------------------------
      // Statistics of the current value of this variate, for the
      // prior distribution
      //------------------------------------------------------------

      Probability priorPdf();
      Probability priorCdf();
      Probability priorPte();

      //------------------------------------------------------------
      // Statistics of the joint distribution
      //------------------------------------------------------------

      Probability jointPdf();

      void operator=(double val);

      void operator=(const Variate& var);
      void operator=(Variate& var);

      // Return the sampling distribution

      Distribution& samplingDistribution() {
	return samplingDistribution_;
      }

      // Return a handle to the prior distribution

      Distribution& prior() {
	return prior_;
      }

      void setIsVariable(bool variable);
      void allowToVary(bool allow);

      // Return true if this variate is allowed to vary

      bool& isVariable() {
	return variable_;
      }

      bool isEnumerated() {
	return isEnumerated_;
      }

      // Sample this variate from its sampling distribution

      void sample();

      // Return the value of the current sample

      double& value();

      void addCorrelationCoefficient(Variate* var, double coeff);
      double getCorrelationCoefficient(Variate* var);

      virtual void setVal(double val, std::string units) {
	val_          = val;
	hasValue_     = true;
	wasSpecified_ = true;
	units_        = units;
      }

      virtual double getVal(double val, std::string units) {
	return val;
      }

      virtual void setVal(std::string val) {
	ThrowColorError("You cannot specify a value (" << val << ") for this variate as a unitless string", "red");
      }

      virtual void setUnits(std::string units) {
	return;
      }

      virtual double getUnitVal() {
	return getUnitVal(val_);
      }

      virtual double getUnitVal(double val) {
	return (val - unitsToNative_.offset_) / unitsToNative_.scaleFactor_;
      }

      virtual double geUnitVal(std::string units) {
	ThrowError("Inheritor has not defined any getUnitVal()");
	return 0.0;
      }

      double getRawVal(double unitVal) {
	return unitVal * unitsToNative_.scaleFactor_ + unitsToNative_.offset_;
      }

      virtual std::string getStringVal() {
	return "";
      };

      std::string units() {
	return units_;
      }

      void setDisplayRange(double min, double max) {
	hasRange_ = true;
	displayMin_ = min;
	displayMax_ = max;
      }

      void setDisplay(bool display);

      void setDisplayOrder(unsigned order) {
	displayOrder_ = order;
      }

      void setName(std::string owner, std::string name);

      // Return true if this variate can be specified without units

      virtual bool unitlessAllowed();

      virtual void printUnits() {};

      bool isDerivable();
      bool canBePrimary();

      bool hasValue_;

      // True if a use requested that this variate be derived

      bool wasRequested_;

      // True if a value for this variate was specified, or if it was
      // explicitly requested to derive this variate

      bool wasSpecified_;

      // True if this variate is derived from another
      // variate.  False if it is a primary variate

      bool isDerived_;

      // True if this variate is malleable (can be primary or derived)

      bool isMalleable_;

      // True if this variate should be made visible to user space
      // (set to false for variates that may be required
      // intermediates, but which we don't want to expose to users)

      bool isVisible_;

      // True if this variate is required by another derived variate

      bool isRequired_;

      // True if this variate is a prerequisite (to be calculated only
      // once).  False if it must be calculated at each step (in a
      // Markov chain, for example)

      bool isPrerequisite_;

      // True if this variable was loaded from a file

      bool loadedFromFile_;

      // True if this variate is used by any model

      bool isUsed_;

      // Utility variables for storing a name for this variate and an
      // explanatory comment

      std::string name_;
      std::string comment_;

      void throwIfNotSpecified();

      // Convenience variables for displaying this variate

      double displayMin_;
      double displayMax_;
      bool hasRange_;
      unsigned displayOrder_;
      bool display_;
      bool isDefaultDisplay_;

    public:
      
      friend class Distribution;

      // Every variate will have a value

      double val_;

      // Every variate can have a sampling distribution (that it is
      // drawn from)

      Distribution samplingDistribution_;

      // Every variate can have a prior distribution

      Distribution prior_;

      // Every variate can be variable or fixed (for use with Models)

      bool variable_;
      bool allowedToVary_;

      // Some variates can be enumerated

      bool isEnumerated_;

      // Every variate can have an arbitrary map of correlations with other variates

      std::map<Variate*, double> correlationCoefficientMap_;

      // Every variate can have an arbitrary multiplier and offset to
      // convert specified units to native units

      Conversion unitsToNative_;
      std::string units_;

    public:

      //------------------------------------------------------------
      // Infrastructure for use in calculating derived variates
      //------------------------------------------------------------

      std::list<Variate*> dependsOn_;
      std::list<Variate*> dependedOnBy_;

      void dependsOn(Variate& var);
      void doesntDependOn(Variate& var);

      void dependedOnBy(Variate& var);
      void notDependedOnBy(Variate& var);

      void removeDependency(Variate* var, std::list<Variate*>& vLlist);

      void listVarsDependedOn();
      bool doesntDependOnAnyone();
      bool doesntDependOnPrimaries();

      void clearDependencies();
      void checkForCircularDependencies(Variate* varTest=0, Variate* varStart=0);

      VAR_DERIVE_FN(*deriveFn_);
      void* deriveArgs_;

      void printDependencies();
      void deriveWith(VAR_DERIVE_FN(deriveFn), void* args);
      void derive();
      bool onlyDependsOnVars(std::map<Variate*,Variate*>& vars);
      bool isRequired();
      void removeIrrelevantDependencies();

      void deriveFrom(Variate& var);

      Variate* variateToDeriveFrom_;

      static VAR_DERIVE_FN(copyVal);

      virtual void plotPdf(double min, double max, unsigned npt);

    }; // End class Variate

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_VARIATE_H
