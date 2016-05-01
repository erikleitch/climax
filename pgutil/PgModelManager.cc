#include "gcp/pgutil/PgModelManager.h"
#include "cpgplot.h"

using namespace std;

using namespace gcp::util;

enum {
  B_NORM=0,
  B_LINE=1,
  B_RECT=2,
  B_YRNG=3,
  B_XRNG=4,
  B_YVAL=5,
  B_XVAL=6,
  B_CROSS=7
};

/**.......................................................................
 * Constructor.
 */
PgModelManager::PgModelManager() 
{
  display_ = true;
}

/**.......................................................................
 * Destructor.
 */
PgModelManager::~PgModelManager() {}

/**.......................................................................
 * Display all models managed by this class
 */
void PgModelManager::display() {
  for(unsigned iMod=0; iMod < models_.size(); iMod++) {
    models_[iMod].precalculateShape();
    models_[iMod].draw();
  }
};

/**.......................................................................
 * Remove a model from the list of models maintained by this class
 */
void PgModelManager::removeModel(float x, float y) {
  float sep,minSep;
  int minInd=0;

  std::vector<PgModel>::iterator minIter=models_.end();

  for(std::vector<PgModel>::iterator iter=models_.begin(); iter != models_.end(); iter++) {
    PgModel& model = *iter;
    sep = sqrt((model.xMid_ - x)*(model.xMid_ - x) + (model.yMid_ - y)*(model.yMid_ - y));
    if(iter == models_.begin()) {
      minSep = sep;
      minIter = iter;
    } else {
      if(sep < minSep) {
	minSep = sep;
	minIter = iter;
      }
    }
  }
    
  if(minIter != models_.end())
    models_.erase(minIter);
}

void PgModelManager::getModel(float xstart, float ystart, bool& read, char& key, std::string unit, Trans& trans)
{
  bool done   = false;
  bool cancel = false;
  float xMid, yMid;
  float xRad1, yRad1;
  float xRad2, yRad2;
  float xtemp, ytemp;

  float xIn=xstart;
  float yIn=ystart;
  xtemp = xIn;
  ytemp = yIn;

  std::ostringstream removeHelp;
  removeHelp << "Use one of the following: " << std::endl
	     << "  A (left mouse) -- Select a model to remove" << std::endl
	     << "Or use 'C' to cancel model selection";

  std::ostringstream addHelp;
  addHelp << "Use one of the following: " << std::endl
	 << "  A -- Arnaud model" << std::endl
	 << "  B -- Symmetric Beta model" << std::endl
	 << "  E -- Elliptical Beta model" << std::endl
	 << "  G -- Gaussian model" << std::endl
	 << "  P -- Point source model" << std::endl
	 << "Or use 'C' to cancel model selection";

  std::ostringstream actionHelp;
  actionHelp << "Use one of the following: " << std::endl
	     << "  A -- Add a new model component" << std::endl
	     << "  R -- Remove the model component nearest the cursor" << std::endl
	     << "  D -- Toggle displaying model components" << std::endl
	     << "  P -- Print model components" << std::endl
	     << "Or use 'C' to cancel model selection";

  COUT(actionHelp.str());

  bool haveModelType = false;
  bool getCenter  = false;
  bool getRadius1 = false;
  bool getRadius2 = false;
  int  cursorMode = B_NORM;
  bool accepted   = false;
  bool getAction  = true;
  bool add = false;
  bool remove = false;
  PgModel model;
  
  int ci;
  cpgqci(&ci);
  cpgsci(5);

  do {
    accepted = 0;

    cpgband(cursorMode, 0, xIn, yIn, &xtemp, &ytemp, &key);

    if(islower((int) key))
      key = (char) toupper((int) key);

    //------------------------------------------------------------
    // If determining what the user wants to do
    //------------------------------------------------------------

    if(getAction) {
      switch(key) {
      case 'A':
	getAction = false;
	add = true;
	COUT(addHelp.str());
	break;
      case 'R':
	cursorMode = B_CROSS;
	getAction = false;
	remove = true;
	COUT(removeHelp.str());
	break;
      case 'P':
	printModels(unit, trans);
	accepted = true;
	read = false;
	break;
      case 'D':
	display_ = !display_;
	accepted = true;
	read = false;
	key = 'L';
	break;
      case 'C':
	accepted = true;
	cancel   = true;
	break;
      default:
	COUT(actionHelp.str());
	break;
      }

      //------------------------------------------------------------
      // Else if removing a model component
      //------------------------------------------------------------

    } else if(remove) {

      switch(key) {
      case 'A':
	COUT("Removing model");
	removeModel(xtemp, ytemp);
	accepted = true;
	read = false;
	key = 'L';
	COUT("Removing model done");
	break;
      case 'C':
	accepted = true;
	cancel   = true;
	break;
      default:
	COUT(removeHelp.str());
	break;
      }

      //------------------------------------------------------------
      // Else if defining a new model component
      //------------------------------------------------------------

    } else if(add) {
      switch(key) {
	// Arnaud model
      case 'A':

	if(getCenter) {
	  getCenter  = false;
	  getRadius1 = true;

	  // Store the center of the model

	  model.xMid_ = xtemp;
	  model.yMid_ = ytemp;
	  model.peak_ = trans.valNearestToPoint(xtemp, ytemp);

	  COUT("Setting peak to: " << model.peak_);

	  // And set the anchor point to the midpoint

	  xIn = xMid=xtemp;
	  yIn = yMid=ytemp;

	  if(model.type_ == PgModel::TYPE_DELTA) {
	    cursorMode = B_NORM;
	    accepted = true;
	    add = true;
	  } else {
	    cursorMode = B_LINE;
	    COUT("Click to define the core radius");
	  }

	  break;

	} else if(getRadius1) {
	  model.xRad1_ = xtemp;
	  model.yRad1_ = ytemp;
	  getRadius1 = false;

	  // If we just got a radius, we are done for symmetric models

	  if(model.type_ == PgModel::TYPE_ARNAUD || model.type_ == PgModel::TYPE_BETA) {
	    accepted = true;
	    add = true;
	  } else {
	    cpgmove(model.xMid_, model.yMid_);
	    cpgdraw(model.xRad1_, model.yRad1_);
	    getRadius2 = true;
	    COUT("Click to define the other axis");
	  }

	  break;

	  // Else we need to get a second radius for assymetric models

	} else if(getRadius2) {

	  model.xRad2_ = xtemp;
	  model.yRad2_ = ytemp;
	  getRadius2 = false;
	  accepted = true;
	  add = true;

	  break;
	}

	// Point-source model

      case 'P':

	// Symmetric beta model

      case 'B':

	// Elliptical beta model

      case 'E':

	// Gaussian model

      case 'G':
	model.type_ = PgModel::keyToType(key);
	COUT("Model typew is now " << model.type_);
	haveModelType = true;
	getCenter = true;
	COUT("Click to locate the model center");
	cursorMode = B_CROSS;
	break;
      case 'C':
	accepted = true;
	cancel = true;
	break;
      default:
	COUT(addHelp.str());
	break;
      }
    }
  } while(!accepted);

  if(add && !cancel) {
    model.rectify();
    model.draw();
    models_.push_back(model);
  }
  
  cpgsci(ci);
}

void PgModelManager::printModels(string& unit, Trans& trans)
{
  for(unsigned iMod=0; iMod < models_.size(); iMod++)
    models_[iMod].print(unit, iMod, trans);
}

void PgModelManager::addModel(PgModel& model)
{
  models_.push_back(model);
}
