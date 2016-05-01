#include "gcp/pgutil/ContourPlot.h"

#include "cpgplot.h"

using namespace std;

using namespace gcp::util;

ContourPlot::ContourPlot()
{
}

/**.......................................................................
 * Constructor.
 */
ContourPlot::ContourPlot(PgUtil* parent, 
			 int ndata, float *zdata, int nx,int ny, 
			 float xmin, float xmax, float ymin, float ymax, 
			 float *flag,float z1, float z2, int ncontour, bool wasSigned,
			 std::string xlab, std::string ylab, std::string title, std::string unit, bool expand) :
  Plot2D(parent, ndata, zdata, nx, ny, xmin, xmax, ymin, ymax, flag, z1, z2, xlab, ylab, title, unit, expand)
{
  contours_.resize(abs(ncontour));
  nContour_ = ncontour;
  wasSigned_ = wasSigned;
  initialize(expand);
  setupContours();
}

ContourPlot::ContourPlot(Plot2D& plot2d, int ncontour, bool wasSigned) :
  Plot2D(plot2d)
{
  contours_.resize(abs(ncontour));
  nContour_ = ncontour;
  wasSigned_ = wasSigned;
  setupContours();
}

/**.......................................................................
 * Destructor.
 */
ContourPlot::~ContourPlot() {}

void ContourPlot::setupContours()
{
  //------------------------------------------------------------
  // If plotting positive-only contours, range from 0 - max
  //------------------------------------------------------------

  if(wasSigned_ && nContour_ > 0) {

    float dz = (zrng_.absMax() - 0.0)/(contours_.size()-1);
    for(unsigned icont=0; icont < contours_.size(); icont++)
      contours_[icont] = zrng_.absMax() - nContour_*dz + icont * dz;

    //------------------------------------------------------------
    // If plotting negative-only contours, range from min - 0
    //------------------------------------------------------------

  } else if(wasSigned_ && nContour_ < 0) {

    float dz = (0.0 - zrng_.absMin())/(contours_.size()-1);
    for(unsigned icont=0; icont < contours_.size(); icont++)
      contours_[icont] = zrng_.absMin() + icont * dz;

    //------------------------------------------------------------
    // Else range from min - max
    //------------------------------------------------------------

  } else {

    float dz = (zrng_.absMax() - zrng_.absMin())/(contours_.size()-1);
    for(unsigned icont=0; icont < contours_.size(); icont++)
      contours_[icont] = zrng_.absMin() + icont * dz;

  }
}

/**.......................................................................
 * Display the greyscale plot
 */
void ContourPlot::display(bool erase)
{
  //------------------------------------------------------------
  // Set the viewport and window
  //------------------------------------------------------------

  if(erase)
    cpgsci(0);
  else
    cpgsci(1);

  setupPlotBoundaries();

  //------------------------------------------------------------
  // Set up colormap
  //------------------------------------------------------------

  setupColormap();

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
  COUT("Contour display X1 = " << x1 << " x2 = " << x2  << " y1 = " << y1 << " y2 = " << y2);
  COUT("Contour display tr_ = " << tr_[0] << " " << tr_[1] << " " << tr_[2] << " " << tr_[3] << " " << tr_[4] << " " << tr_[5]);

  if(erase)
    cpgimag(erase ? &zeros_[0] : data_, nx_, ny_, i1_, i2_, j1_, j2_, zrngCurr_.absMax(), zrngCurr_.absMin(), tr_);
  else {
    cpgcons(data_, nx_, ny_, i1_, i2_, j1_, j2_, &contours_[0], contours_.size(), tr_);

    cpgsci(5);
    float conlev = 10.0;
    cpgcont(data_, nx_, ny_, i1_, i2_, j1_, j2_, &conlev, 1, tr_);
  }

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
  // Draw a header if using a header
  //------------------------------------------------------------

  parent_->drawHeader();
}

