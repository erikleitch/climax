#include "gcp/pgutil/GreyscalePlot.h"

#include "cpgplot.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
GreyscalePlot::GreyscalePlot(PgUtil* parent, 
			     int ndata, float *zdata, int nx,int ny, 
			     float xmin, float xmax, float ymin, float ymax, 
			     float *flag,float z1, float z2,
 			     std::string xlab, std::string ylab, std::string title, std::string unit) :
  Plot2D(parent, ndata, zdata, nx, ny, xmin, xmax, ymin, ymax, flag, z1, z2, xlab, ylab, title, unit)
{
}

/**.......................................................................
 * Destructor.
 */
GreyscalePlot::~GreyscalePlot() {}

/**.......................................................................
 * Display the greyscale plot
 */
void GreyscalePlot::display(bool erase)
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

  float x1,x2,y1,y2;
  cpgqwin(&x1,&x2,&y1,&y2);
  float xvp1,xvp2,yvp1,yvp2;
  cpgqwin(&xvp1,&xvp2,&yvp1,&yvp2);

#if 0
  COUT("erase = " << erase << " Drawing with tr_ = " << tr_);
  COUT(tr_[0]);
  COUT(tr_[1]);
  COUT(tr_[2]);
  COUT(tr_[3]);
  COUT(tr_[4]);
  COUT(tr_[5]);
#endif

  cpgimag(erase ? &zeros_[0] : data_, nx_, ny_, i1_, i2_, j1_, j2_, zrngCurr_.absMax(), zrngCurr_.absMin(), tr_);

  cpgsci(1);

  //------------------------------------------------------------
  // Draw models on top, if displaying models
  //------------------------------------------------------------
  
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
