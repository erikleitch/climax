#include "gcp/pgutil/Plot2D.h"
#include "gcp/pgutil/ContourPlot.h"
#include "gcp/pgutil/LinePlot.h"

#include "gcp/util/Stats.h"

#include "cpgplot.h"

#include <string.h>

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

Plot2D::Plot2D() {}

/**.......................................................................
 * Constructor.
 */
Plot2D::Plot2D(PgUtil* parent, 
	       int ndata, float *zdata, int nx,int ny, 
	       float xmin, float xmax, float ymin, float ymax, 
	       float *flag,float z1, float z2, 
	       std::string xlab, std::string ylab, std::string title, std::string unit, bool expand) :
  Plot(parent, ndata, xlab, ylab, title)
{
  data_     = zdata;

  nx_       = nx;
  ny_       = ny;

  bound_.setTo(xmin, xmax, ymin, ymax, parent->reverseX_, parent->reverseY_);
  zrng_.setTo(z1, z2, ndata, zdata);

  autoscale_ = (z1==z2);

  unit_    = unit;

  plotType_ = PLOT_PLOT2D;

  initialize(expand);
}

Plot2D::Plot2D(Plot2D& plot2d)
{
  *this = plot2d;
}

/**.......................................................................
 * Destructor.
 */
Plot2D::~Plot2D() {}

/**.......................................................................
 * Plot this image
 */
void Plot2D::plot()
{
  key_ = G_DIS;

  if(docurs_) {
    
    fprintf(stdout,"For HELP, hit the \'%c\' key on your keyboard\n", G_HELP);
    
    do {
      cancel_ = false;
      
      switch(key_) {

	//------------------------------------------------------------
	// Zoom the plot
	//------------------------------------------------------------
	
      case G_ZOOM:
	zoom(G_ZOOM, ZOOM_BOTH);
	break;

      case G_HORI:
	zoom(G_HORI, ZOOM_X);
	break;

      case G_VERT:
	zoom(G_VERT, ZOOM_Y);
	break;

      case G_RAD:
	//linePlot();
	
	multiLinePlot();
	break;

	//------------------------------------------------------------
	// Re-display the plot
	//------------------------------------------------------------

      case G_DIS:
	display();
	break;
	
	//------------------------------------------------------------
	// Compute statistics on a selected region of the plot.
	//------------------------------------------------------------

      case G_STAT:
	printStats();
	break;

	//------------------------------------------------------------
	// Define a new model component
	//------------------------------------------------------------

      case G_MODEL:
	parent_->pgManager_.getModel(xpos_[0], ypos_[0], read_, key_, unit_, trans_);
	break;
	
      case G_CONT:
	{
	  int num=0;
	  bool wasSigned=false;
	  getNumber(num, wasSigned);

	  if(!cancel_) {
	    
	    //------------------------------------------------------------
	    // Extract only the data within the current bound
	    //------------------------------------------------------------

	    unsigned nxSub, nySub;
	    std::vector<float> subData = getCurrentlyDisplayedData(nxSub, nySub);

	    ContourPlot cp(*this, num, wasSigned);

	    // Clear this plot

	    display(true);

	    // Now plot the new object
	    
	    cp.plot();
	    
	    // Now clear the new object

	    cp.display(true);

	    // Now redisplay us

	    display();
	  }
	}
	break;

	//------------------------------------------------------------
	// Print information about the nearest point
	//------------------------------------------------------------

      case G_INS:
	getInfo();
	break;
      case G_HELP:     // Print usage info 
	fprintf(stdout,"\nYou requested help by pressing \'%c\'.\n", G_HELP);
	fprintf(stdout,"All cursor positions are entered with \'%c\' key (Left mouse button)\n", G_CUR);
	fprintf(stdout,"\n %c - Select a sub-image to be displayed\n", G_ZOOM);
	fprintf(stdout," %c - Calculate statistics on a subset of data\n", G_STAT);
	fprintf(stdout," %c - Redisplay current plot\n", G_DIS);
	fprintf(stdout," %c - Add/remove/display models\n", G_MODEL);
	fprintf(stdout," %c - Fiddle contrast & brightness\n", G_FID);
	fprintf(stdout," %c - Use greyscale\n", G_GREY);
	fprintf(stdout," %c - Use rainbow colormap\n", G_RAIN);
	fprintf(stdout," %c - Use heat colormap\n", G_HEAT);
	fprintf(stdout," %c - Make a contour plot\n", G_CONT);
	fprintf(stdout," %c - Toggle crosshair cursor\n", G_CROS);
	fprintf(stdout,"\nTo end this session hit the \'%c\' key (Right mouse button)\n", G_QUIT);
	fprintf(stdout,"\n");
	break;
      default :
	break;
      case G_GREY:
	cmap_ = grey_;
	display();
	break;
      case G_RAIN:
	cmap_ = rain_;
	display();
	break;
      case G_HEAT:
	cmap_ = heat_; 
	display();
	break;
      case G_CROS:
	{
	  if(cursor_ == B_CROSS)
	    cursor_ = B_NORM;
	  else
	    cursor_ = B_CROSS;
	}
	break;
      case G_FID:
	resetContrast();
	break;
      }

      cpgband(cursor_, 0, xpos_[0], ypos_[0], &xpos_[0], &ypos_[0], &key_);

      if(islower((int) key_))
	key_ = (char) toupper((int) key_);
      
    } while(key_ != G_QUIT);
    
    //------------------------------------------------------------
    // Else if no cursor, just plot and exit.
    //------------------------------------------------------------
    
    
  } else {
    display();
  }
  
  //------------------------------------------------------------
  // Close PGPLOT only if it was opened in this function
  //------------------------------------------------------------
  
  if(!wasopen_)
    cpgend(); 
}

void Plot2D::initialize(bool expand)
{
  //------------------------------------------------------------
  // Setup internal variables
  //------------------------------------------------------------

  for(unsigned i=0; i < n_std_cmap; i++) {
    if(strcmp(std_cmaps[i].name.c_str(),"grey")==0) 
      grey_ = &std_cmaps[i];
    if(strcmp(std_cmaps[i].name.c_str(),"rainbow")==0) 
      rain_ = &std_cmaps[i];
    if(strcmp(std_cmaps[i].name.c_str(),"heat")==0) 
      heat_ = &std_cmaps[i];
  }

  cmap_ = parent_->cmap_;

  bright_   =  0.5;
  contrast_ = -1.0;
  cursor_   =  B_NORM;

  //------------------------------------------------------------
  // Now store information about the data we are displaying
  //------------------------------------------------------------

  dx_ = bound_.xrng_.range()/(nx_);
  dy_ = bound_.yrng_.range()/(ny_);

  // Set the transformation matrix for the data array.

  if(bound_.dataXmin() < bound_.dataXmax()) {
    tr_[0] = bound_.absXmin() - 0.5*dx_;
    tr_[1] = dx_;
    tr_[2] = 0.0;
  } else {
    tr_[0] = bound_.absXmax() + 0.5*dx_;
    tr_[1] = -dx_;
    tr_[2] = 0.0;
  }

  if(bound_.dataYmin() < bound_.dataYmax()) {
    tr_[3] = bound_.absYmin() - 0.5*dy_;
    tr_[4] = 0.0;
    tr_[5] = dy_;
  } else {
    tr_[3] = bound_.absYmax() + 0.5*dy_;
    tr_[4] = 0.0;
    tr_[5] = -dy_;
  }

  i1_ = j1_ = 1;
  i2_ = nx_;
  j2_ = ny_;

  //------------------------------------------------------------
  // Now that the transformation matrix is set up, recompute display
  // limits if requested
  //------------------------------------------------------------

  if(parent_->usedefs_) {
    double xmin, xmax, ymin, ymax;
    xmin = (parent_->xmin_ == parent_->xmax_) ? bound_.plotXmin() : parent_->xmin_;
    xmax = (parent_->xmin_ == parent_->xmax_) ? bound_.plotXmax() : parent_->xmax_;
    ymin = (parent_->ymin_ == parent_->ymax_) ? bound_.plotYmin() : parent_->ymin_;
    ymax = (parent_->ymin_ == parent_->ymax_) ? bound_.plotYmax() : parent_->ymax_;

    bound_.setTo(xmin, xmax, ymin, ymax);
  }

  //------------------------------------------------------------
  // Store transform info
  //------------------------------------------------------------

  trans_.dx_    = fabs(dx_);
  trans_.dy_    = fabs(dy_);
  trans_.nx_    = nx_;
  trans_.ny_    = ny_;
  trans_.ndata_ = ndata_;
  trans_.xmins_ = bound_.absXmin();
  trans_.ymins_ = bound_.absYmin();
  trans_.xmaxs_ = bound_.absXmax();
  trans_.ymaxs_ = bound_.absYmax();
  trans_.zdata_ = data_;
  trans_.reverseX_ = parent_->reverseX_;
  trans_.reverseY_ = parent_->reverseY_;

  //------------------------------------------------------------
  // Expand the plot boundaries.  One more dx adds an extra pixel-wide
  // space around the image
  //------------------------------------------------------------

  if(expand) {
    bound_.xrng_.expandAbs(1.0*dx_);
    bound_.yrng_.expandAbs(1.0*dy_);
  }

  boundCurr_ = bound_;
  zrngCurr_  = zrng_;

  //------------------------------------------------------------
  // Keep a copy of the array with all zeros
  //------------------------------------------------------------

  zeros_.resize(ndata_);
  for(unsigned i=0; i < ndata_; i++)
    zeros_[i] = zrng_.absMin();
}


/**.......................................................................
 * Display the greyscale plot
 */
void Plot2D::display(bool erase)
{
  //------------------------------------------------------------
  // Set the viewport and window
  //------------------------------------------------------------

  if(erase) {

    cpgsci(0);

    parent_->drawBox(erase ? 0 : 1);
    parent_->drawLabels(xlab_, ylab_, title_);
    parent_->drawWedge(zrngCurr_.absMin(), zrngCurr_.absMax(), unit_);

    float x1,x2,y1,y2;
    cpgqvp(0, &x1, &x2, &y1, &y2);

    cpgsvp(0,1,0,1);
    float xp1,xp2,yp1,yp2;
    cpgqvp(3, &xp1, &xp2, &yp1, &yp2);

    float pixPerNdcy = yp2-yp1;
    float pixPerNdcx = xp2-xp1;

    float ndcxToAdd = (pixPerNdcy / 40 * 4) / pixPerNdcx;

    cpgsvp(x1, x2+ndcxToAdd, y1, y2);

  } else {
    setupPlotBoundaries();
  }

  //------------------------------------------------------------
  // Set up colormap
  //------------------------------------------------------------

  if(erase) {
    double bright   = bright_;
    double contrast = contrast_;

    bright_   =  0.5;
    contrast_ = -1.0;

    setupColormap();
    bright_   = bright;
    contrast_ = contrast;
  } else {
    setupColormap();
  }

  //------------------------------------------------------------
  // If requested, override min/max greyscale with user-specified values
  //------------------------------------------------------------

  if(parent_->zmin_ != parent_->zmax_) {
    zrngCurr_.setTo(parent_->zmin_, parent_->zmax_);
  }

  //------------------------------------------------------------
  // Display the image
  //------------------------------------------------------------

  cpgimag(erase ? &zeros_[0] : data_, nx_, ny_, i1_, i2_, j1_, j2_, zrngCurr_.absMax(), zrngCurr_.absMin(), tr_);
  cpgsci(1);

  //------------------------------------------------------------
  // Draw models on top, if displaying models
  //------------------------------------------------------------

  float x1,x2,y1,y2;
  cpgqwin(&x1,&x2,&y1,&y2);

  if(!erase && parent_->pgManager_.display_) {
    parent_->pgManager_.display();
  }
  
  //------------------------------------------------------------
  // Redraw the plot box
  //------------------------------------------------------------

  parent_->drawBox(erase ? 0 : 1);

  //------------------------------------------------------------
  // Redraw labels
  //------------------------------------------------------------

  if(erase)
    cpgsci(0);

  parent_->drawLabels(xlab_, ylab_, title_);

  //------------------------------------------------------------
  // Redraw the wedge
  //------------------------------------------------------------

  if(!erase)
    parent_->drawWedge(zrngCurr_.absMin(), zrngCurr_.absMax(), unit_);
  
  //------------------------------------------------------------
  // Draw a header if using a header
  //------------------------------------------------------------

  parent_->drawHeader();
}
	
/**.......................................................................
 * Print statistics about the plot
 */
void Plot2D::printStats()
{
  float xtemp,ytemp;

  static string mess1("\n %c - Select start and end vertices using this key\n");
  static string mess2(" %c - Select the whole plot\n");
  static string mess3(" %c - Abort selection\n");

  int oldcol;
  cpgqci(&oldcol);
  cpgsci(5);

  dofull_ = 0;
  cancel_ = 0;

  for(unsigned iter = 0;iter<2 && !dofull_ && !cancel_;iter++) {
    bool accepted;
    do {
      accepted = 0;

      cpgband((iter==0) ? B_NORM : B_RECT, 0, xtemp, ytemp, &xtemp,
	      &ytemp, &key_);

      if(islower((int) key_))
	key_ = (char) toupper((int) key_);

      xpos_[iter] = xtemp;
      ypos_[iter] = ytemp;

      switch(key_) {
      case G_STAT:
	accepted = dofull_ = true;
	break;
      case G_CAN:      // Abort box selection
	accepted = cancel_ = 1;
	break;
      case G_QUIT:     // Quit now
	
	// Close PGPLOT only if it was opened in this function

	if(!wasopen_)
	  cpgend(); 
	return;
	break;
      case G_CUR:             // Accept the selected start vertex
	accepted=1;
	break;
      default:            // Unexpected cursor input key - show usage 
	fprintf(stdout,mess1.c_str(), G_CUR);
	fprintf(stdout,mess2.c_str(), G_ZOOM);
	fprintf(stdout,mess3.c_str(), G_CAN);
	break;
      };
    } while(!accepted);
  };
  
  if(dofull_) {
    xpos_[0] = boundCurr_.plotXmin();
    ypos_[0] = boundCurr_.plotYmin();
    xpos_[1] = boundCurr_.plotXmax();
    ypos_[1] = boundCurr_.plotYmax();
  }

  trans_.printStats(xpos_[0], xpos_[1], ypos_[0], ypos_[1]);

  //------------------------------------------------------------
  // And reset the old color
  //------------------------------------------------------------

  cpgsci(oldcol);
}

/**.......................................................................
 * Print info about the nearest point
 */
void Plot2D::getInfo()
{
  // Find the nearest grid point

  float val = trans_.valNearestToPoint(xpos_[0], ypos_[0]);
  std::string xstr, ystr;

  if(parent_->coordCallback_) {
    std::string xstr, ystr;
    (*(parent_->coordCallback_))(xpos_[0], ypos_[0], xstr, ystr, parent_->coordCallbackArgs_);
    fprintf(stdout,"Pixel value at (%s, %s) is %g %s\n", xstr.c_str(), ystr.c_str(), val, unit_.c_str());
    fprintf(stdout,"Pixel value at (%lf, %lf) is %g %s\n", xpos_[0], ypos_[0], val, unit_.c_str());
  } else {
    fprintf(stdout,"Pixel value at (%lf, %lf) is %g %s\n", xpos_[0], ypos_[0], val, unit_.c_str());
  }
}

/**.......................................................................
 * Install a colormap
 */
void Plot2D::setupColormap()
{
  cpgctab(cmap_->l,cmap_->r,cmap_->g,cmap_->b,cmap_->n,contrast_,bright_);
}

/**.......................................................................
 * Fiddle the contrast/brightness
 */
void Plot2D::resetContrast()
{
  double xmid = boundCurr_.xrng_.midPoint();
  double ymid = boundCurr_.yrng_.midPoint();
  
  contrast_ = 5.0 * (ypos_[0]-ymid)/(ypos_[0] < ymid ? (boundCurr_.plotYmin()-ymid) : -(boundCurr_.plotYmax()-ymid));
  bright_   = 0.5 + 1.0 * (fabs(contrast_)+1.0)*((xpos_[0] - boundCurr_.plotXmax())/(boundCurr_.plotXmin() - boundCurr_.plotXmax()) - 0.5);
  
  setupColormap();
  
  //cpgband(cursor_, 0, xpos_[0], ypos_[0], &xpos_[0], &ypos_[0], &key_);
  read_ = true;
  display();
}

void Plot2D::multiLinePlot()
{
  PlotBound boundSave     = bound_;
  PlotBound boundCurrSave = boundCurr_;
  Point pt1, pt2;
  std::vector<float> xdata;
  std::vector<float> ydata;

  bool first = true;
  do {
    getLine(xdata, ydata, pt1, pt2, first);

    cpgmove(pt1.x_, pt1.y_);
    cpgdraw(pt2.x_, pt2.y_);

    pt1 = pt2;
  } while(!cancel_);

  LinePlot lp(parent_, xdata.size(), &xdata[0], &ydata[0], 0, xlab_, unit_, title_, true, false);
  display(true);
  
  PgUtil::setWnad(false);
  lp.plot();
  lp.display(true);
  PgUtil::setWnad(true);

  display();

  bound_ = boundSave;
  boundCurr_ = boundCurrSave;
}

void Plot2D::linePlot()
{
  PlotBound bound;
  Point pt1, pt2;
  getRange(G_RAD, ZOOM_LINE, bound, pt1, pt2);

  std::vector<float> xdata;
  std::vector<float> ydata;

  //------------------------------------------------------------ 
  // First get the angle the user-defined line makes with the x-axis
  //------------------------------------------------------------

  double dx = parent_->reverseX_ ? (pt1.x_ - pt2.x_) : (pt2.x_ - pt1.x_);
  double dy = (pt2.y_ - pt1.y_);
  Angle theta;
  theta.setRadians(atan2(dy,dx));

  double ct = cos(theta.radians());
  double st = sin(theta.radians());

  // Now compute n points along this line

  double len = bound.range();
  unsigned n = len/fabs(dx_) + 1;

  for(unsigned i=0; i < n; i++) {
    float xrot = dx_ * i;
    float yrot = 0.0;

    // Rotate back into the unrotated, shifted coordinate system

    float xshift = xrot * ct - yrot * st;
    float yshift = xrot * st + yrot * ct;

    // And unshift

    float x = xshift + (parent_->reverseX_ ? pt2.x_ : pt1.x_);
    float y = yshift + (parent_->reverseY_ ? pt2.x_ : pt1.y_);

    xdata.push_back(xrot);
    //ydata.push_back(y);
    ydata.push_back(trans_.convolveAroundPoint(x, y, data_));
  }

  LinePlot lp(parent_, xdata.size(), &xdata[0], &ydata[0], 0, xlab_, ylab_, title_, true, false);
  display(true);

  PgUtil::setWnad(false);
  lp.plot();
  lp.display(true);
  PgUtil::setWnad(true);

  display();
}

void Plot2D::getLine(std::vector<float>& xdata, std::vector<float>& ydata, Point& pt1, Point& pt2, bool& first)
{
  if(first) {
    Plot::getLine(ZOOM_LINE, pt1, pt2);
    first = false;
  } else {
    xpos_[0] = pt1.x_;
    ypos_[0] = pt1.y_;
    getPoint(ZOOM_LINE, pt2);
  }

  PlotBound bound;
  bound.setTo(pt1.x_, pt2.x_, pt1.y_, pt2.y_);

  //------------------------------------------------------------ 
  // First get the angle the user-defined line makes with the x-axis
  //------------------------------------------------------------

  double dx = parent_->reverseX_ ? (pt1.x_ - pt2.x_) : (pt2.x_ - pt1.x_);
  double dy = (pt2.y_ - pt1.y_);
  Angle theta;
  theta.setRadians(atan2(dy,dx));

  double ct = cos(theta.radians());
  double st = sin(theta.radians());

  // Now compute n points along this line

  double len = bound.range();
  unsigned n = len/fabs(dx_) + 1;

  double xStart = xdata.size() > 0 ? xdata[xdata.size()-1] : 0.0;

  for(unsigned i=0; i < n; i++) {
    float xrot = dx_ * i;
    float yrot = 0.0;

    // Rotate back into the unrotated, shifted coordinate system

    float xshift = xrot * ct - yrot * st;
    float yshift = xrot * st + yrot * ct;

    // And unshift

    float x = xshift + (parent_->reverseX_ ? pt2.x_ : pt1.x_);
    float y = yshift + (parent_->reverseY_ ? pt2.x_ : pt1.y_);

    xdata.push_back(xrot + xStart);
    ydata.push_back(trans_.convolveAroundPoint(x, y, data_));
  }
}
