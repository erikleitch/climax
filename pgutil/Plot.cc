#include "gcp/pgutil/PgUtil.h"
#include "gcp/pgutil/Plot.h"

#include "gcp/util/Stats.h"
#include "gcp/util/String.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Plot::Plot() {}

Plot::Plot(PgUtil* parent, unsigned ndata, std::string xlab, std::string ylab, std::string title)
{
  parent_ = parent;

  xlab_   = xlab;
  ylab_   = ylab;
  title_  = title;

  nx_ = ndata;
  ny_ = ndata;

  ndata_ = ndata;

  initialize();
}

void Plot::initialize()
{
  //------------------------------------------------------------
  // Query if we have a cursor or not, and whether the plot was
  // already open
  //------------------------------------------------------------

  anchor_.x_ = 0.0;
  anchor_.y_ = 0.0;

  parent_->queryDevice(wasopen_);
  docurs_ = parent_->haveCursor();

  cpgsch(parent_->ch_);

  dofull_    = true;
  cancel_    = true;
  key_       = G_DIS;
  autoscale_ = true;

  // Store the current viewport and window

  cpgqvp(0, &xvp1_, &xvp2_, &yvp1_, &yvp2_);
  cpgqwin(&xwin1_, &xwin2_, &ywin1_, &ywin2_);
}

/**.......................................................................
 * Destructor.
 */
Plot::~Plot() {}

/**.......................................................................
 * Get new plot boundaries if requested, and recompute displayed ranges
 */
void Plot::zoom(char keyFull, ZoomType type)
{
  PlotBound bound;
  Point pt1, pt2;
  getRange(keyFull, type, bound, pt1, pt2);

  if(!cancel_) {

    //------------------------------------------------------------
    // First erase the current display
    //------------------------------------------------------------
    
    COUT("Erasing display...");
    display(true);
    COUT("Erasing display...done");

    //------------------------------------------------------------
    // Now set the new boundary
    //------------------------------------------------------------

    boundCurr_ = bound;
	
    //------------------------------------------------------------
    // Compute the new greyscale boundaries, if this is a 2D plot
    //------------------------------------------------------------

    if(plotType_ == PLOT_PLOT2D)
      computeGreyscaleRange();

    //------------------------------------------------------------
    // And redisplay
    //------------------------------------------------------------

    COUT("Redrawing display...");
    display();
    COUT("Redrawing display... done");
  }
}

void Plot::computeGreyscaleRange()
{
  if(autoscale_ || !dofull_){
    float x,y;
    int first=1,ind;
	
    double min, max;
    for(unsigned i=0; i < nx_; i++)
      for(unsigned j=0; j < ny_; j++) {
	x = bound_.xrng_.absMin() + dx_*i;
	y = bound_.yrng_.absMin() + dy_*j;
	ind = i+j*nx_;
	    
	//------------------------------------------------------------
	// Only iterate over currently displayed points
	//------------------------------------------------------------

	if(boundCurr_.contains(x, y)) {

	  if(isnan(data_[ind]))
	    continue;

	  if(first) {
	    min = max = data_[ind];
	    first = 0;
	  }

	  min = MIN(min, data_[ind]);
	  max = MAX(max, data_[ind]);
	}
      }

    zrngCurr_.setTo(min, max);

  } else {
    zrngCurr_ = zrng_;
  }
}

std::vector<float> Plot::getCurrentlyDisplayedData(unsigned& nxSub, unsigned& nySub)
{
  unsigned iXMin = (boundCurr_.xrng_.absMin() - bound_.xrng_.absMin())/fabs(dx_);
  unsigned iXMax = (boundCurr_.xrng_.absMax() - bound_.xrng_.absMin())/fabs(dx_);

  unsigned iYMin = (boundCurr_.yrng_.absMin() - bound_.yrng_.absMin())/fabs(dy_);
  unsigned iYMax = (boundCurr_.yrng_.absMax() - bound_.yrng_.absMin())/fabs(dy_);

  nxSub = (iXMax - iXMin) + 1;
  nySub = (iYMax - iYMin) + 1;

  std::vector<float> subData(nxSub*nySub);

  for(unsigned i=iXMin; i <= iXMax; i++) {
    for(unsigned j=iYMin; j <= iYMax; j++) {
      unsigned ind = i+j*nx_;
      unsigned indSub = (i-iXMin) + (j-iYMin)*nxSub;
      subData[indSub] = data_[ind];
    }
  }

  return subData;
}

/**.......................................................................
 * Get new plot boundaries if requested, and recompute displayed ranges
 */
void Plot::getRange(char keyFull, ZoomType type, PlotBound& bound, Point& pt1, Point& pt2)
{
  static string mess1("\n %c - Select start and end vertices using this key\n");
  static string mess2(" %c - Select the whole plot\n");
  static string mess3(" %c - Abort selection\n");

  unsigned uType = (unsigned)type;
  float xtemp, ytemp;

  int oldcol;
  cpgqci(&oldcol);
  cpgsci(5);

  dofull_ = false;
  cancel_ = false;

  bool accepted;

  int cursorMode;

  switch(type) {
  case ZOOM_BOTH:
    cursorMode = B_RECT;
    break;
  case ZOOM_X:
    cursorMode = B_XRNG;
    break;
  case ZOOM_Y:
    cursorMode = B_YRNG;
    break;
  case ZOOM_LINE:
    cursorMode = B_LINE;
    break;
  default:
    break;
  }

  xtemp = xpos_[0];
  ytemp = ypos_[0];

  //------------------------------------------------------------
  // Starting point is the current bound.  Note that this assignment
  // also preserves the sense of the current boundary, so we don't
  // have to worry about reversing axes if plotting with -ve x values
  // to the right, for example
  //------------------------------------------------------------

  bound = boundCurr_;

  for(unsigned iter = 0; iter < 2 && !dofull_ && !cancel_; iter++) {

    do {

      accepted = false;
      cpgband(iter==0 ? B_NORM : cursorMode, 0, xtemp, ytemp, &xtemp, &ytemp, &key_);
      if(islower((int) key_))
	key_ = (char) toupper((int) key_);

      xpos_[iter] = xtemp;
      ypos_[iter] = ytemp;

      if(iter == 0) {
	pt1.x_ = xtemp;
	pt1.y_ = ytemp;
      } else {
	pt2.x_ = xtemp;
	pt2.y_ = ytemp;
      }

      switch(key_) {
      case G_CAN:      // Abort box selection
	accepted = cancel_ = true;
	break;
      case G_QUIT:     // Quit now 
	if(!wasopen_)
	  cpgend();
	return;
	break;
      case G_CUR:             // Accept the selected start vertex
	accepted=true;
	break;
      default:            // Unexpected cursor input key - show usage 
	if(key_ == keyFull) {
	  accepted = dofull_ = true;
	} else {
	  fprintf(stdout,mess1.c_str(), G_CUR);
	  fprintf(stdout,mess2.c_str(), keyFull);
	  fprintf(stdout,mess3.c_str(), G_CAN);
	}
	break;
      };
    } while(!accepted);
  };

  if(!cancel_) {

    //------------------------------------------------------------
    // If the full range was requested, revert to stored ranges
    //------------------------------------------------------------

    if(dofull_) {
      bound = bound_;

      // Otherwise set the new boundaries.  Don't worry about the
      // sense of the plot axes, since this is preserved in the range
      // objects themselves

    } else {
	
      if(uType & ZOOM_X)
	bound.xrng_.setTo(xpos_[0], xpos_[1]);
      
      if(uType & ZOOM_Y)
	bound.yrng_.setTo(ypos_[0], ypos_[1]);

      if(uType & ZOOM_LINE)
	bound.setTo(xpos_[0], xpos_[1], ypos_[0], ypos_[1]);
    }
  }

  pt1.x_ = xpos_[0];
  pt1.y_ = ypos_[0];
  pt2.x_ = xpos_[1];
  pt2.y_ = ypos_[1];
}

void Plot::getPoint(ZoomType type, Point& pt)
{
  static string mess1("\n %c - Select start and end vertices using this key\n");
  static string mess3(" %c - End selection\n");

  unsigned uType = (unsigned)type;
  float xtemp, ytemp;

  int oldcol;
  cpgqci(&oldcol);
  cpgsci(5);

  cancel_ = false;

  bool accepted;

  int cursorMode;

  switch(type) {
  case ZOOM_BOTH:
    cursorMode = B_RECT;
    break;
  case ZOOM_X:
    cursorMode = B_XRNG;
    break;
  case ZOOM_Y:
    cursorMode = B_YRNG;
    break;
  case ZOOM_LINE:
    cursorMode = B_LINE;
    break;
  default:
    break;
  }

  // Start at the current buffered cursor position

  xtemp = anchor_.x_;
  ytemp = anchor_.y_;

  for(unsigned iter = 1; iter < 2 && !dofull_ && !cancel_; iter++) {

    do {

      accepted = false;
      cpgband(iter==0 ? B_NORM : cursorMode, 0, xtemp, ytemp, &xtemp, &ytemp, &key_);
      if(islower((int) key_))
	key_ = (char) toupper((int) key_);

      anchor_.x_ = xtemp;
      anchor_.y_ = ytemp;

      pt = anchor_;

      switch(key_) {
      case G_CAN:      // Abort box selection
	accepted = cancel_ = true;
	break;
      case G_QUIT:     // Quit now 
	cancel_ = true;
	return;
	break;
      case G_CUR:             // Accept the selected start vertex
	accepted=true;
	break;
      default:            // Unexpected cursor input key - show usage 
	fprintf(stdout,mess1.c_str(), G_CUR);
	fprintf(stdout,mess3.c_str(), G_CAN);
	break;
      };
    } while(!accepted);
  };
}

/**.......................................................................
 * Get new plot boundaries if requested, and recompute displayed ranges
 */
void Plot::getLine(ZoomType type, Point& pt1, Point& pt2)
{
  static string mess1("\n %c - Select start and end vertices using this key\n");
  static string mess3(" %c - End selection\n");

  unsigned uType = (unsigned)type;
  float xtemp, ytemp;

  int oldcol;
  cpgqci(&oldcol);
  cpgsci(5);

  dofull_ = false;
  cancel_ = false;

  bool accepted;

  int cursorMode;

  switch(type) {
  case ZOOM_BOTH:
    cursorMode = B_RECT;
    break;
  case ZOOM_X:
    cursorMode = B_XRNG;
    break;
  case ZOOM_Y:
    cursorMode = B_YRNG;
    break;
  case ZOOM_LINE:
    cursorMode = B_LINE;
    break;
  default:
    break;
  }

  // Start at the current buffered cursor position

  xtemp = anchor_.x_;
  ytemp = anchor_.y_;

  for(unsigned iter = 0; iter < 2 && !dofull_ && !cancel_; iter++) {

    do {

      accepted = false;
      cpgband(iter==0 ? B_NORM : cursorMode, 0, xtemp, ytemp, &xtemp, &ytemp, &key_);
      if(islower((int) key_))
	key_ = (char) toupper((int) key_);

      anchor_.x_ = xtemp;
      anchor_.y_ = ytemp;

      if(iter == 0) {
	pt1 = anchor_;
      } else {
	pt2 = anchor_;
      }

      switch(key_) {
      case G_CAN:      // Abort box selection
	accepted = cancel_ = true;
	break;
      case G_QUIT:     // Quit now 
	cancel_ = true;
	return;
	break;
      case G_CUR:             // Accept the selected start vertex
	accepted=true;
	break;
      default:            // Unexpected cursor input key - show usage 
	fprintf(stdout,mess1.c_str(), G_CUR);
	fprintf(stdout,mess3.c_str(), G_CAN);
	break;
      };
    } while(!accepted);
  };
}

/**.......................................................................
 * Get new plot boundaries if requested, and recompute displayed ranges
 */
void Plot::getEndOfRange(char keyFull, ZoomType type, PlotBound& bound, Point& pt1, Point& pt2)
{
  static string mess1("\n %c - Select start and end vertices using this key\n");
  static string mess2(" %c - Select the whole plot\n");
  static string mess3(" %c - Abort selection\n");

  unsigned uType = (unsigned)type;
  float xtemp, ytemp;

  int oldcol;
  cpgqci(&oldcol);
  cpgsci(5);

  dofull_ = false;
  cancel_ = false;

  bool accepted;

  int cursorMode;

  switch(type) {
  case ZOOM_BOTH:
    cursorMode = B_RECT;
    break;
  case ZOOM_X:
    cursorMode = B_XRNG;
    break;
  case ZOOM_Y:
    cursorMode = B_YRNG;
    break;
  case ZOOM_LINE:
    cursorMode = B_LINE;
    break;
  default:
    break;
  }

  xtemp = pt1.x_;
  ytemp = pt1.y_;
  xpos_[0] = xtemp;
  ypos_[0] = ytemp;

  //------------------------------------------------------------
  // Starting point is the current bound.  Note that this assignment
  // also preserves the sense of the current boundary, so we don't
  // have to worry about reversing axes if plotting with -ve x values
  // to the right, for example
  //------------------------------------------------------------

  bound = boundCurr_;

  for(unsigned iter = 1; iter < 2 && !dofull_ && !cancel_; iter++) {

    do {

      accepted = false;
      cpgband(iter==0 ? B_NORM : cursorMode, 0, xtemp, ytemp, &xtemp, &ytemp, &key_);
      if(islower((int) key_))
	key_ = (char) toupper((int) key_);

      xpos_[iter] = xtemp;
      ypos_[iter] = ytemp;

      if(iter == 0) {
	pt1.x_ = xtemp;
	pt1.y_ = ytemp;
      } else {
	pt2.x_ = xtemp;
	pt2.y_ = ytemp;
      }

      switch(key_) {
      case G_CAN:      // Abort box selection
	accepted = cancel_ = true;
	break;
      case G_QUIT:     // Quit now 
	if(!wasopen_)
	  cpgend();
	return;
	break;
      case G_CUR:             // Accept the selected start vertex
	accepted=true;
	break;
      default:            // Unexpected cursor input key - show usage 
	if(key_ == keyFull) {
	  accepted = dofull_ = true;
	} else {
	  fprintf(stdout,mess1.c_str(), G_CUR);
	  fprintf(stdout,mess2.c_str(), keyFull);
	  fprintf(stdout,mess3.c_str(), G_CAN);
	}
	break;
      };
    } while(!accepted);
  };

  if(!cancel_) {

    //------------------------------------------------------------
    // If the full range was requested, revert to stored ranges
    //------------------------------------------------------------

    if(dofull_) {
      bound = bound_;

      // Otherwise set the new boundaries.  Don't worry about the
      // sense of the plot axes, since this is preserved in the range
      // objects themselves

    } else {
	
      if(uType & ZOOM_X)
	bound.xrng_.setTo(xpos_[0], xpos_[1]);
      
      if(uType & ZOOM_Y)
	bound.yrng_.setTo(ypos_[0], ypos_[1]);

      if(uType & ZOOM_LINE)
	bound.setTo(xpos_[0], xpos_[1], ypos_[0], ypos_[1]);
    }
  }

  pt1.x_ = xpos_[0];
  pt1.y_ = ypos_[0];
  pt2.x_ = xpos_[1];
  pt2.y_ = ypos_[1];
}

/**.......................................................................
 * Get a numerical string from the cursor
 */
void Plot::getNumber(int& num, bool& wasSigned)
{
  int keymax=5;
  ostringstream os;
  int i=0;

  cancel_ = false;

  fprintf(stdout,"\nEnter a number, followed by a carriage return.\n");

  do {

    cpgband(B_NORM, 0, xpos_[0], ypos_[0], &xpos_[0], &ypos_[0], &key_);

    if(islower((int) key_))
      key_ = (char) toupper((int) key_);

    switch(key_) {

    case G_QUIT: 
    case G_CAN: 
      cancel_ = true;
      break;
    case G_HELP:
      fprintf(stdout,"\nEnter a number, followed by a carriage return.\n");
      break;
    case '\n': case '\r':
      os << '\0';
      break;
    default:
      os << key_;
      ++i;
      break;
    }
  } while(!(key_ == '\n' || key_ == '\r') && i < keymax);

  // If a number
  
  String str(os.str());
  wasSigned = str.contains("-") || str.contains("+");

  char* ptr=0;
  num = (int)strtod(os.str().c_str(), &ptr);
  
  if(*ptr != '\0') {
    fprintf(stderr,"Not a number.\n");
    cancel_ = true;
  }
}

void Plot::setupPlotBoundaries()
{
  //------------------------------------------------------------
  // Advance to the next page if not overplotting
  //------------------------------------------------------------

  if(!parent_->overplot_)
    cpgpage();
  
  //------------------------------------------------------------
  // Set the viewport if requested, else revert to the viewport that
  // was set externally
  //------------------------------------------------------------

  if(parent_->vp_) {
    cpgvstd();
  } else {
    cpgsvp(xvp1_, xvp2_, yvp1_, yvp2_);
  }

  //------------------------------------------------------------
  // Set the window if requested
  //------------------------------------------------------------

  if(parent_->win_) {
    cpgswin(boundCurr_.plotXmin(), boundCurr_.plotXmax(), boundCurr_.plotYmin(), boundCurr_.plotYmax());
  }

  //------------------------------------------------------------
  // Adjust the viewport if requested
  //------------------------------------------------------------
  
  if(parent_->wnad_) {
    float xw1,xw2,yw1,yw2;
    cpgqwin(&xw1, &xw2, &yw1, &yw2);
    cpgwnad(xw1, xw2, yw1, yw2);
  }
}
