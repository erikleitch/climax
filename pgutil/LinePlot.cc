#include "gcp/pgutil/LinePlot.h"
#include "gcp/pgutil/PgUtil.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
LinePlot::LinePlot(PgUtil* parent,
		   int narr, float* xarr, float* yarr, float* earr,
		   std::string xlab, std::string ylab, std::string title, 
		   bool doLine, bool doErr) :
  Plot(parent, narr, xlab, ylab, title)
{
  if(doErr && earr == 0)
    ThrowSimpleColorError("Error plot was requested, but no errors were specified", "red");

  xdata_    = xarr;
  data_     = yarr;
  edata_    = earr;

  doLine_   = doLine;
  doErr_    = doErr;

  plotType_ = PLOT_PLOT1D;

  initialize();
}

LinePlot::~LinePlot() {}

/**.......................................................................
 * Destructor.
 */

void LinePlot::initialize()
{
  //------------------------------------------------------------
  // We store a copy of the data, in case log plotting was requested
  //------------------------------------------------------------

  xs_.resize(ndata_);
  ys_.resize(ndata_);

  for(unsigned i=0; i < ndata_; i++) {
    if(parent_->logPlot_) {
      xs_[i] = log10(xdata_[i]);
      ys_[i] = log10(data_[i]);
    } else {
      xs_[i] = xdata_[i];
      ys_[i] = data_[i];
    }
  }

  //------------------------------------------------------------
  // If an error array was specified, store the 
  //------------------------------------------------------------

  if(doErr_) {
    for(unsigned i=0; i < ndata_; i++) {
      yDataErrHi_.resize(ndata_);
      yDataErrLo_.resize(ndata_);
      
      yDataErrHi_[i] = data_[i] + edata_[i];
      yDataErrLo_[i] = data_[i] - edata_[i];
    }
  }

  //------------------------------------------------------------
  // Calculate the min/max of the data
  //------------------------------------------------------------

  double xmin, xmax, ymin, ymax;

  xmin = xmax = xdata_[0];
  ymin = (edata_ ? data_[0] - edata_[0] : data_[0]);
  ymax = (edata_ ? data_[0] - edata_[0] : data_[0]);

  for(unsigned i=0; i < ndata_;i++) {
    xmin = (xmin < xdata_[i]) ? xmin : xdata_[i];
    xmax = (xmax > xdata_[i]) ? xmax : xdata_[i];
    float ylo = edata_ ? data_[i]-edata_[i] : data_[i];
    float yhi = edata_ ? data_[i]+edata_[i] : data_[i];
    ymin = (ymin < ylo) ? ymin : ylo;
    ymax = (ymax > yhi) ? ymax : yhi;
  }

  if(parent_->logPlot_) {
    xmin = log10(xmin);
    xmax = log10(xmax);
    ymin = log10(ymin);
    ymax = log10(ymax);
  }

  xmin -= 0.1*(xmax-xmin);
  xmax += 0.1*(xmax-xmin);

  ymin -= 0.1*(ymax-ymin);
  ymax += 0.1*(ymax-ymin);

  if(ymin == ymax) {
    ymin -= 0.1*ymin;
    ymax += 0.1*ymax;
  }

  if(parent_->usedefs_) {

    if(parent_->xmin_ != parent_->xmax_) {
      xmin = parent_->logPlot_ ? log10(parent_->xmin_) : parent_->xmin_;
      xmax = parent_->logPlot_ ? log10(parent_->xmax_) : parent_->xmax_;
    }

    if(parent_->ymin_ != parent_->ymax_) {
      ymin = parent_->logPlot_ ? log10(parent_->ymin_) : parent_->ymin_;
      ymax = parent_->logPlot_ ? log10(parent_->ymax_) : parent_->ymax_;
    }
  }

  bound_.setTo(xmin, xmax, ymin, ymax, parent_->reverseX_, parent_->reverseY_);
  boundCurr_ = bound_;
}

void LinePlot::plot()
{
  if(docurs_) {

    int cursor = B_NORM;
    printf("For HELP, hit the \'%c\' key on your keyboard\n", G_HELP);

    do {

      switch(key_) {
      case G_ZOOM:
	zoom(G_ZOOM, ZOOM_BOTH);
	break;

      case G_HORI:
	zoom(G_HORI, ZOOM_X);
	break;

      case G_VERT:
	zoom(G_VERT, ZOOM_Y);
	break;

      case G_DIS:
	display();
	break;

      case G_INS:
	printValNearestToPoint(xpos_[0], ypos_[0]);
	break;

      case G_FLG:
	cacheHighestPoint();
	labelPoints(false);
	break;

      case G_BEAM:
	cacheNHighestPeaks();
	labelPoints(false);
	break;

      case G_CUT:
	display(true);
	clearLabeledPoints();
	display();
	break;

	//------------------------------------------------------------
	// Toggle crosshair cursor
	//------------------------------------------------------------

      case G_CROS:

	if(cursor == B_CROSS)
	  cursor = B_NORM;
	else
	  cursor = B_CROSS;
	break;

      case G_HELP:     /* Print usage info */
	printf("\nYou requested help by pressing \'%c\'.\n", G_HELP);
	printf("All cursor positions are entered with \'%c\' key (Left mouse button)\n", G_CUR);
	printf(" %c - Clear labeled peaks\n", G_CUT);
	printf(" %c - Label highest peak\n", G_FLG);
	printf(" %c - Label highest peak with just a number\n", G_NXT);
	printf(" %c - Label highest peaks over a range\n", G_BEAM);
	printf(" %c - Select X range to be displayed (hit %c twice for full range)\n", G_HORI, G_HORI);
	printf(" %c - Select Y range to be displayed (hit %c twice for full range)\n", G_VERT, G_VERT);
	printf(" %c - Select a sub-image to be displayed.\n", G_ZOOM);
	printf(" %c - Redisplay current plot.\n", G_DIS);
	printf("\nTo end this session hit the \'%c\' key (Right mouse button)\n", G_QUIT);
	printf("\n");
	break;
#if 0
      case KEY_FLG:
	v_lwrite(narr, xarr, yarr, xmin, xmax, ymin,ymax);
	redisp = 0;
	break;
      case KEY_NXT:
	v_lnum(narr, xarr, yarr, xmin, xmax, ymin,ymax);
	redisp = 0;
	break;
      case KEY_BEAM: /* Label the ten highest peaks */
	v_lten(narr, &xs[0], &ys[0], xmin, xmax, ymin,ymax);
	redisp = 0;
	break;
#endif
      default :
	break;
      }

      cpgband(cursor, 0, xpos_[0], ypos_[0], &xpos_[0], &ypos_[0], &key_);

      if(islower((int) key_))
	key_ = (char) toupper((int) key_);

    } while(key_ != G_QUIT);
  } else {
    display();
  }

  //------------------------------------------------------------
  // Close PGPLOT only if it was opened in this function
  //------------------------------------------------------------
  
  if(!wasopen_)
    cpgend(); 

  return;
}

/**.......................................................................
 * Display the line plot
 */
void LinePlot::display(bool erase)
{
  //------------------------------------------------------------
  // Set the viewport and window
  //------------------------------------------------------------

  if(erase)
    cpgsci(0);
  else
    cpgsci(1);

  setupPlotBoundaries();

  labelPoints(erase);

  //------------------------------------------------------------
  // Redraw the plot box
  //------------------------------------------------------------

  parent_->drawBox(erase ? 0 : 1);

  //------------------------------------------------------------
  // Label the plot
  //------------------------------------------------------------

  parent_->drawLabels(xlab_, ylab_, title_);

  //------------------------------------------------------------
  // Now draw the data
  //------------------------------------------------------------

  if(parent_->useTraceCi_)
    cpgsci(erase ? 0 : parent_->traceCi_);
  else
    cpgsci(erase ? 0 : 10);

  if(doLine_) {
    if(parent_->logPlot_) {
      cpgline(ndata_, &xs_[0], &ys_[0]);
    } else {
      cpgline(ndata_, xdata_, data_);
    }
  } else {

    if(doErr_) {
      cpgerry(ndata_, xdata_, &yDataErrHi_[0], &yDataErrLo_[0], 1.0);
    } else {
      cpgpt(ndata_, xdata_, data_, -1);
    }
  }

  // And reset the color to the box color

  if(parent_->useBoxCi_)
    cpgsci(parent_->boxCi_);
  else
    cpgsci(1);
}

void LinePlot::printValNearestToPoint(float x, float y)
{
  bool first=true;
  float dist,distmin,xtmp,ytmp,val;
  float xmin,ymin;

  for(unsigned i=0; i < ndata_; i++) {

    xtmp = xdata_[i];
    ytmp = data_[i];

    dist = sqrt((xtmp-x)*(xtmp-x) + (ytmp-y)*(ytmp-y));

    if(dist < distmin || first) {
      xmin = xtmp;
      ymin = ytmp;
      distmin = dist;
      first = false;
    }
  }
}

void LinePlot::cacheHighestPoint()
{
  PlotBound bound;
  Point pt1, pt2;
  getRange(G_HORI, ZOOM_X, bound, pt1, pt2);

  if(!cancel_) {

    //------------------------------------------------------------
    // Now find the highest point over this x range
    //------------------------------------------------------------
    
    bool first=true;
    unsigned imax;
    float maxVal;
    
    for(unsigned i=0; i < ndata_; i++) {
      
      float x = xdata_[i];
      float y = data_[i];
      
      if(bound.containsX(x)) {
	
	if(first) {
	  imax = i;
	  maxVal = y;
	  first = false;
	}
	
	if(y > maxVal) {
	  imax = i;
	  maxVal = y;
	}
      }
    }

    float y1 = maxVal +   boundCurr_.yrng_.range()/32;
    float y2 = maxVal + 5*boundCurr_.yrng_.range()/32;
    float xVal = xdata_[imax];

    cachePoint(xVal, xVal, y1, y2, true);
  }
}

void LinePlot::cachePoint(float x1, float x2, float y1, float y2, bool doLine)
{
  std::ostringstream os;
  os << setprecision(3) << std::fixed << x1;
      
  LabeledPoint pt;
  pt.bound_.setTo(x1, x2, y1, y2);
  pt.label_ = os.str();
  pt.doLine_ = doLine;

  labeledPoints_.push_back(pt);
}

/*.......................................................................
 * Module for v_lplot.
 */
void LinePlot::cacheNHighestPeaks()
{
  int n;
  bool wasSigned;

  getNumber(n, wasSigned);

  if(cancel_)
    return;

  PlotBound bound;
  Point pt1, pt2;
  getRange(G_HORI, ZOOM_X, bound, pt1, pt2);

  if(!cancel_) {

    int *index=NULL;
    
    // Sort the arrays

    PgUtil::indexx(ndata_, data_, &index);

    int i=0;
    for(unsigned j=0,i=ndata_-1;i > 0;i--) {
      float xval = xdata_[index[i]];
      float yval = data_[index[i]];

      if(j < n && bound.containsX(xval)) {
	float y = yval + boundCurr_.yrng_.range()/32;
	cachePoint(xval, xval, y, y, false);
	j++;
      }
    }

    if(index)
      free(index);
  }
}

void LinePlot::clearLabeledPoints()
{
  labeledPoints_.resize(0);
}

void LinePlot::labelPoints(bool erase)
{
  int oldcol;
  cpgqci(&oldcol);
  cpgsci(erase ? 0 : 5);

  for(unsigned i=0; i < labeledPoints_.size(); i++) {
    LabeledPoint& pt = labeledPoints_[i];

    if(boundCurr_.containsX(pt.bound_.absXmin())) {

      if(pt.doLine_) {
	cpgmove(pt.bound_.absXmin(), pt.bound_.absYmin());
	cpgdraw(pt.bound_.absXmax(), pt.bound_.absYmax());
      }

      cpgptxt(pt.bound_.absXmax(), pt.bound_.absYmax() + boundCurr_.yRange()/64,0,0.5, pt.label_.c_str());
    }

  }

  cpgsci(oldcol);
}
