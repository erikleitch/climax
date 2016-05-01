#include "gcp/models/CosmologyModel.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
CosmologyModel::CosmologyModel() 
{
  initialize();
}

/**.......................................................................
 * Destructor.
 */
CosmologyModel::~CosmologyModel() {}

void CosmologyModel::initialize()
{
  addComponent(H0_);
  addComponent(z_);
  addComponent(omegaM_);
  addComponent(omegaL_);

  z_.allowUnitless(true);
  omegaM_.allowUnitless(true);
  omegaL_.allowUnitless(true);

  addComponentName(H0_,     "H0",     "The Hubble constant");
  addComponentName(z_,      "z",      "The redshift");
  addComponentName(omegaM_, "omegaM", "Omega matter");
  addComponentName(omegaL_, "omegaL", "Omega lambda");

  wasComputed_ = false;
  isVariable_  = false;
  isUsed_      = false;
  wasPruned_   = false;

  //------------------------------------------------------------
  // And initialize all components to fixed
  //------------------------------------------------------------

  initializeComponentsToFixed();
}

/**.......................................................................
 * Called whenever any variable components of this model are updated
 */
void CosmologyModel::update()
{
  //  COUT("this = " << this << " Inside update() : needs = " << needsUpdating());

  if(needsUpdating()) {

    //    COUT("Setting z to " << z_.value());

    cosmology_.setRedshift(z_.value());
    cosmology_.setOmegaM(omegaM_.value());
    cosmology_.setOmegaL(omegaL_.value());
    cosmology_.setH0(H0_);

    dA_ = cosmology_.angularDiameterDistance();
    H_  = cosmology_.H();

    wasComputed_ = true;
  }
}

bool CosmologyModel::needsUpdating()
{
  //------------------------------------------------------------
  // If we have not yet ever computed the cosmology, check the inputs
  // that we need
  //------------------------------------------------------------

  if(!wasComputed_) {

    checkVar("z");
    checkVar("omegaM");
    checkVar("omegaL");
    checkVar("H0");

    isVariable_ = 
      (z_.isVariable()) || 
      (omegaM_.isVariable()) || 
      (omegaL_.isVariable()) || 
      (H0_.isVariable());

    return true;

  } else {

    pruneVariates();
 
    return isVariable_ && isUsed_;
  }
}

/**.......................................................................
 * Compute the possible range of dA represented by this cosmology
 */
void CosmologyModel::angularDiameterDistanceRange(Length& daMin, Length& daMax)
{
  if(!isVariable_) {
    daMin = cosmology_.angularDiameterDistance();
    daMax = cosmology_.angularDiameterDistance();
  } else {
    daMin = cosmology_.angularDiameterDistance() / 2;
    daMax = cosmology_.angularDiameterDistance() * 2;
  }
}

/**.......................................................................
 * Mark our variates as unused is nobody is using this model
 */
void CosmologyModel::pruneVariates()
{
  if(wasPruned_)
    return;

  //------------------------------------------------------------
  // If this model is not used by anybody, mark our variates as unused
  //------------------------------------------------------------

  if(!isUsed_) {
    for(unsigned iVar=0; iVar < componentVec_.size(); iVar++) {
      Variate* var = componentVec_[iVar];
      var->isUsed_ = false;
    }
  }

  wasPruned_ = true;
}
