#include "gcp/util/Exception.h"
#include "gcp/util/Variate.h"

#include "gcp/pgutil/PgUtil.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Variate::Variate() 
{
  unitsToNative_.scaleFactor_ = 1.0;
  unitsToNative_.offset_      = 0.0;

  units_            = " ";

  allowedToVary_    = true;
  hasValue_         = false;
  isEnumerated_     = false;
  hasRange_         = false;

  wasRequested_     = false;
  wasSpecified_     = false;
  isDerived_        = false;
  isMalleable_      = false;
  isVisible_        = true;
  isRequired_       = false;
  isPrerequisite_   = false;
  loadedFromFile_   = false;

  isUsed_           = true;

  displayMin_       = 0.0;
  displayMax_       = 0.0;
  displayOrder_     = 0;

  display_          = true;
  isDefaultDisplay_ = true;

  variateToDeriveFrom_ = 0;
}

/**.......................................................................
 * Destructor.
 */
Variate::~Variate() {}

Probability Variate::pdf()
{
  return samplingDistribution_.pdf(val_);
}

Probability Variate::cdf()
{
  return samplingDistribution_.cdf(val_);
}

Probability Variate::pte()
{
  return samplingDistribution_.pte(val_);
}

Probability Variate::samplingPdf()
{
  return samplingDistribution_.pdf(val_);
}

Probability Variate::samplingCdf()
{
  return samplingDistribution_.cdf(val_);
}

Probability Variate::samplingPte()
{
  return samplingDistribution_.pte(val_);
}

Probability Variate::priorPdf()
{
  return prior_.pdf(val_);
}

Probability Variate::priorCdf()
{
  return prior_.cdf(val_);
}

Probability Variate::priorPte()
{
  return prior_.pte(val_);
}

Probability Variate::jointPdf()
{
  return samplingDistribution_.pdf(val_) * prior_.pdf(val_);
}

void Variate::sample()
{
  val_ = samplingDistribution_.sample();
}

double& Variate::value()
{
  return val_;
}

void Variate::operator=(double val)
{
  val_ = val;
}

void Variate::addCorrelationCoefficient(Variate* var, double coeff)
{
  correlationCoefficientMap_[var] = coeff;
}

double Variate::getCorrelationCoefficient(Variate* var)
{
  std::map<Variate*, double>::iterator iter = correlationCoefficientMap_.find(var);

  if(iter == correlationCoefficientMap_.end()) {
    ThrowError("No correlation with var " << var << " specified");
  }

  return iter->second;
}

/**.......................................................................
 * This is a check for unit'ed quantities.  Only bare variates are
 * allowed to be specified as unitless.  Inheritors that are
 * ConformableQuantities should not allow values to be specified
 * without units
 */
bool Variate::unitlessAllowed()
{
  return true;
}

void Variate::setIsVariable(bool variable)
{
  if(variable) {

    if(allowedToVary_) {
      variable_ = variable;
    } else {
      ThrowColorError("You must specify a fixed value for this variable", "red");
    }

  } else {
    variable_ = variable;
  }

}

void Variate::allowToVary(bool allow)
{
  allowedToVary_ = allow;
}

void Variate::throwIfNotSpecified()
{
  if(!wasSpecified_) {
    ThrowSimpleColorError("Variate " << name_ << " wasn't specified", "red")
  }
}

bool Variate::isDerivable()
{
  return isDerived_ || isMalleable_;
}

bool Variate::canBePrimary()
{
  return !isDerived_ || isMalleable_;
}

void Variate::setDisplay(bool display)
{
  display_ = display;
  isDefaultDisplay_ = false;
}

//=======================================================================
// Infrastructure for handling derived variates
//=======================================================================

void Variate::dependsOn(Variate& var)
{
  dependsOn_.push_back(&var);
  var.dependedOnBy(*this);
}

void Variate::doesntDependOn(Variate& var)
{
  removeDependency(&var, dependsOn_);
  var.notDependedOnBy(*this);
}

void Variate::notDependedOnBy(Variate& var)
{
  removeDependency(&var, dependedOnBy_);
}

void Variate::removeDependency(Variate* var, std::list<Variate*>& vList)
{
  std::list<Variate*>::iterator rmIter;
  bool found = false;

  for(std::list<Variate*>::iterator iter=vList.begin(); iter != vList.end(); iter++) {
    if(*iter == var) {
      rmIter = iter;
      found = true;
      break;
    }
  }

  if(found)
    vList.erase(rmIter);
}

void Variate::dependedOnBy(Variate& var)
{
  dependedOnBy_.push_back(&var);
}

void Variate::listVarsDependedOn()
{
  COUT("Var " << name_ << " depends on");
  for(std::list<Variate*>::iterator iter=dependsOn_.begin(); iter != dependsOn_.end(); iter++) {
    COUT((*iter)->name_);
  }
}

bool Variate::doesntDependOnPrimaries()
{
  for(std::list<Variate*>::iterator iter = dependsOn_.begin(); iter != dependsOn_.end(); iter++) {
    Variate* var = *iter;

    COUT("Checking dependencies: derived = " << var->isDerived_);

    if(!var->isDerived_)
      return false;
  }

  return true;
}

bool Variate::doesntDependOnAnyone()
{
  return dependsOn_.size() == 0;
}

void Variate::clearDependencies()
{
  dependsOn_.clear();
  dependedOnBy_.clear();
}

void Variate::checkForCircularDependencies(Variate* varTest, Variate* varStart)
{
  COUT("this = " << this << " Inside cfcD woth varTest = " << varTest 
       << " varStart = " << varStart << " size = " << dependsOn_.size());

  for(std::list<Variate*>::iterator iter=dependsOn_.begin(); iter != dependsOn_.end(); iter++) {
    Variate* var = *iter;

    if(var == varTest || this == varTest) {
      ThrowError("Var " << varTest->name_ 
		 << " has a circular dependency (ultimately depends on " << name_
		 << " but " << name_ << " depends on " << varTest->name_ << ")");
    }

    //------------------------------------------------------------
    // If we hit the same var that we started on, stop traversing the
    // dependency tree.  This means that there is a circular
    // dependency, but not for us
    //------------------------------------------------------------

    if(var == varStart)
      return;

    if(varTest == 0)
      var->checkForCircularDependencies(this, var);
    else
      var->checkForCircularDependencies(varTest, varStart);
  }

}

void Variate::deriveWith(VAR_DERIVE_FN(deriveFn), void* args)
{
  deriveFn_   = deriveFn;
  deriveArgs_ = args;
}

void Variate::derive()
{
  if(deriveFn_ == 0) {
    ThrowColorError("No function has been specified to derive variate " << name_, "red");
  } else {
    deriveFn_(deriveArgs_);
  }
}

/**.......................................................................
 * Return true if this variate was either specified, or is depended on
 * by a variate that was specified
 */
bool Variate::isRequired()
{
  if(wasSpecified_) {
    return true;
  }

  for(std::list<Variate*>::iterator iter=dependedOnBy_.begin(); iter != dependedOnBy_.end(); iter++) {
    Variate* var = *iter;

    if(var->isRequired()) {
      return true;
    }
  }

  return false;
}

bool Variate::onlyDependsOnVars(std::map<Variate*,Variate*>& vars)
{
  for(std::list<Variate*>::iterator iter=dependsOn_.begin(); iter != dependsOn_.end(); iter++) {
    Variate* varTest = *iter;

    if(vars.find(varTest) == vars.end()) {
      return false;
    }
  }

  return true;
}

/**.......................................................................
 * If we depend on a variate that is malleable, remove it from our
 * list of dependencies if it is not currently being derived.
 */
void Variate::removeIrrelevantDependencies()
{
  std::list<Variate*>::iterator rmIter;
  bool found;

  do {
    found = false;
    rmIter = dependsOn_.end();

    for(std::list<Variate*>::iterator iter=dependsOn_.begin(); iter != dependsOn_.end(); iter++) {
      Variate* varTest = *iter;
      if(!varTest->isDerived_) {
	rmIter = iter;
	break;
      }
    }

    if(rmIter != dependsOn_.end()) {
      dependsOn_.erase(rmIter);
      found = true;
    }

  } while(found == true);

  do {
    found = false;
    rmIter = dependedOnBy_.end();
    for(std::list<Variate*>::iterator iter=dependedOnBy_.begin(); iter != dependedOnBy_.end(); iter++) {
      Variate* varTest = *iter;
      if(!varTest->isDerived_) {
	rmIter = iter;
	break;
      }
    }

    if(rmIter != dependedOnBy_.end()) {
      dependedOnBy_.erase(rmIter);
      found = true;
    }

  } while(found == true);

}

void Variate::printDependencies()
{
  COUTCOLOR("Var " << name_ << " depends on: ", "green");
  for(std::list<Variate*>::iterator iter=dependsOn_.begin(); iter != dependsOn_.end(); iter++) {
    Variate* var = *iter;
    COUTCOLOR(var->name_ << " isRequired: " << var->isRequired(), "yellow");
  }

  COUTCOLOR("Var " << name_ << " is depended on by: ", "green");
  for(std::list<Variate*>::iterator iter=dependedOnBy_.begin(); iter != dependedOnBy_.end(); iter++) {
    Variate* var = *iter;
    COUTCOLOR(var->name_ << " isRequired: " << var->isRequired(), "yellow");
  }

  COUTCOLOR("Var " << name_ << " is required: " << isRequired(), "green");
}

void Variate::operator=(const Variate& var)
{
  operator=((Variate&) var);
}

void Variate::operator=(Variate& var)
{
  COUTCOLOR("Inside Variate operator=", "red");
  ThrowError("Using undefined base-class Variate operator=()");
}

/**.......................................................................
 * A static method for deriving this variate from another variate
 */
VAR_DERIVE_FN(Variate::copyVal)
{
  Variate* var = (Variate*)args;

  if(var->variateToDeriveFrom_)
    var->val_ = var->variateToDeriveFrom_->val_;
}

/**.......................................................................
 * Set up this variate to be derived from another variate
 */
void Variate::deriveFrom(Variate& var)
{
  isDerived_ = true;
  wasSpecified_ = var.wasSpecified_;

  if(var.isDerived_ && var.doesntDependOnPrimaries()) {
    XtermManip xtm;
    ThrowSimpleError(COLORIZE(xtm, "red",
			      "You can't derive " << name_ << " from " << var.name_ << " because " << var.name_ << " is already derived variate with no primary dependencies"));
  }

  setIsVariable(true);
  samplingDistribution().setType(Distribution::DIST_GAUSS);
  setUnits(var.units());

  dependsOn(var);

  variateToDeriveFrom_ = &var;
  deriveWith(Variate::copyVal, this);
}

void Variate::plotPdf(double min, double max, unsigned npt)
{
  double dx = (max-min)/(npt-1);

  std::vector<double> x(npt);
  std::vector<double> y(npt);

  for(unsigned i=0; i < npt; i++) {
    x[i] = min + dx*i;
    y[i] = samplingDistribution_.pdf(x[i]).value();
  }

  PgUtil::linePlot(x, y);
}

void Variate::setName(std::string owner, std::string name)
{
  std::ostringstream os;
  os << owner << "." << name;
  name_ = os.str();
}
