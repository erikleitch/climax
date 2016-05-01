#include "gcp/fftutil/ObsManager.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ObsManager::ObsManager() {}

/**.......................................................................
 * Destructor.
 */
ObsManager::~ObsManager() 
{
  for(std::map<std::string, ObsInfo*>::iterator obs = obsMap_.begin(); obs != obsMap_.end(); obs++) {
    delete obs->second;
    obs->second = 0;
  }
}

/**.......................................................................
 * Add a model to the map of models we are managing
 */
ObsInfo* ObsManager::addObs(std::string obsName)
{
  ObsInfo* obs = 0;

  // Make sure the obs hasn't already been declared:

  bool exists = false;

  try {
    obs = getObs(obsName);
    exists = true;
  } catch(...) {
    exists = false;
  }

  if(exists)
    ThrowColorError(std::endl << "An observation by the name of " << obsName << " has already been initialized", "red");

  obs = new ObsInfo();

  // Default this to 1 for now.

  obs->setNumberOfStokesParameters(1);

  // Set this obs' handle 

  obs->name() = obsName;

  // Insert a new node in our list of models

  obsMap_[obsName] = obs;

  // Add parameters needed for the parsing interface

  obs->addParameters();

  // And return it

  return obs;
}

/**.......................................................................
 * Return an observation by name
 */
ObsInfo* ObsManager::getObs(std::string name)
{
  std::map<std::string, ObsInfo*>::iterator obs = obsMap_.find(name);

  if(obs == obsMap_.end())
    ThrowSimpleColorError(std::endl << "No obs named: " << name << " has been initialized", "red");

  return obs->second;
}
