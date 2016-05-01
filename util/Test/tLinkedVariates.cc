#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/UniformVariate.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Sampler.h"

#include "gcp/pgutil/PgUtil.h"

#include <vector>
#include <list>

using namespace std;
using namespace gcp::util;
using namespace gcp::program;


KeyTabEntry Program::keywords[] = {
  { "val",   "0.0", "d", "Val"},
  { "mean",  "0.0", "d", "Mean"},
  { "sigma", "0.0", "d", "Sigma"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

class LinkedVariate {
public:

  LinkedVariate() {
  }

  LinkedVariate(std::string name) {
    name_ = name;
  }

  LinkedVariate(const LinkedVariate& var) {
    *this = var;
  }

  void operator=(const LinkedVariate& var) {
    operator=((LinkedVariate&) var);
  }

  void operator=(LinkedVariate& var) {
    name_ = var.name_;
    dependsOn_ = var.dependsOn_;
    dependedOnBy_ = var.dependedOnBy_;
  }

  void dependsOn(LinkedVariate& var) {
    dependsOn_.push_back(&var);
    var.dependedOnBy(*this);
  }

  void dependedOnBy(LinkedVariate& var) {
    dependedOnBy_.push_back(&var);
  }

  void clearDependencies() {
    dependsOn_.clear();
    dependedOnBy_.clear();
  }

  bool doesntDependOnAnyone()
  {
    return dependsOn_.size() == 0;
  }

  void listVarsDependedOnBy()
  {
    COUT("Var " << name_ << " is depended on by");
    for(std::list<LinkedVariate*>::iterator iter=dependedOnBy_.begin(); iter != dependedOnBy_.end(); iter++) {
      COUT((*iter)->name_);
    }
  }

  void listVarsDependedOn()
  {
    COUT("Var " << name_ << " depends on");
    for(std::list<LinkedVariate*>::iterator iter=dependsOn_.begin(); iter != dependsOn_.end(); iter++) {
      COUT((*iter)->name_);
    }
  }

  void checkForCircularDependencies(LinkedVariate* varTest=0)
  {
    for(std::list<LinkedVariate*>::iterator iter=dependsOn_.begin(); iter != dependsOn_.end(); iter++) {
      LinkedVariate* var = *iter;

      if(var == varTest) {
	ThrowError("Var " << varTest->name_ << " has a circular dependency (ultimately depends on " << name_ 
		   << " but " << name_ << " depends on " << varTest->name_ << ")");
      }

      if(varTest == 0)
	var->checkForCircularDependencies(this);
      else
	var->checkForCircularDependencies(varTest);
    }
  }

  std::string name_;
  std::list<LinkedVariate*> dependsOn_;
  std::list<LinkedVariate*> dependedOnBy_;
};

bool onlyDependsOnVarsAlreadyComputed(LinkedVariate* var, std::map<LinkedVariate*,LinkedVariate*>& alreadyComputedVars);
std::vector<LinkedVariate*> organizeVars(std::vector<LinkedVariate*>& vars);

int Program::main()
{
  std::vector<LinkedVariate*> vars;

  LinkedVariate var1("var1"), var2("var2"), var3("var3"), var4("var4");

  var3.dependsOn(var2);
  var2.dependsOn(var1);
  var2.dependsOn(var4);
  var1.dependsOn(var4);
  var4.dependsOn(var3);
  
  var2.listVarsDependedOn();
  var1.listVarsDependedOnBy();

  COUT("Here 1");
  var1.checkForCircularDependencies();
  COUT("Here 2");
  var2.checkForCircularDependencies();
  COUT("Here 3");
  var3.checkForCircularDependencies();

  vars.push_back(&var1);
  vars.push_back(&var2);
  vars.push_back(&var3);
  vars.push_back(&var4);

  std::vector<LinkedVariate*> orderedVars = organizeVars(vars);

  COUT("Vars can be computed in the following order: " );
  for(unsigned iVar=0; iVar < orderedVars.size(); iVar++) {
    COUT(orderedVars[iVar]->name_);
  }

  return 0;
}

std::vector<LinkedVariate*> organizeVars(std::vector<LinkedVariate*>& vars)
{
  std::vector<LinkedVariate*> orderedVars;
  std::map<LinkedVariate*, LinkedVariate*> alreadyComputedVars;

  unsigned nVarsAdded;

  do {
    nVarsAdded = 0;
    for(unsigned iVar=0; iVar < vars.size(); iVar++) {
      LinkedVariate* var = vars[iVar];
      
      if(alreadyComputedVars.find(var) == alreadyComputedVars.end()) {
	if(var->doesntDependOnAnyone() || onlyDependsOnVarsAlreadyComputed(var, alreadyComputedVars)) {
	  orderedVars.push_back(var);
	  alreadyComputedVars[var] = var;
	  ++nVarsAdded;
	}
      }

    }
  } while(nVarsAdded > 0);

  return orderedVars;
}

bool onlyDependsOnVarsAlreadyComputed(LinkedVariate* var, std::map<LinkedVariate*,LinkedVariate*>& alreadyComputedVars)
{
  for(std::list<LinkedVariate*>::iterator iter=var->dependsOn_.begin(); iter != var->dependsOn_.end(); iter++) {
    LinkedVariate* varTest = *iter;

    if(alreadyComputedVars.find(varTest) == alreadyComputedVars.end())
      return false;
  }

  return true;
}
