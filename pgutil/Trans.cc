#include "gcp/pgutil/Trans.h"
#include "gcp/util/Exception.h"

#include <cmath>
#include <iostream>

#include <stdio.h>

using namespace std;
using namespace gcp::util;

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

/**.......................................................................
 * Constructor.
 */
Trans::Trans() 
{
  reverseX_ = false;
  reverseY_ = false;
}

/**.......................................................................
 * Destructor.
 */
Trans::~Trans() {}

void Trans::rectify(float& min, float& max)
{
  if(max < min) {
    float tmp = min;
    min = max;
    max = tmp;
  }
}

void Trans::printStats(float x1, float x2, float y1, float y2)
{
  rectify(xmins_, xmaxs_);
  rectify(ymins_, ymaxs_);

  rectify(x1, x2);
  rectify(y1, y2);

  bool first = true;
  bool firstGZ = true;
  double min = 0.0;
  double minGZ = 0.0;
  double max = 0.0;
  double mean = 0.0;
  double sd = 0.0;
  unsigned n = 0;

  for(unsigned i=0; i < ndata_; i++) {
    unsigned yind = i/nx_;
    unsigned xind = i - yind*nx_;

    double xtemp = xmins_ + dx_/2 + xind*dx_;
    double ytemp = ymins_ + dy_/2 + yind*dy_;

    if(xtemp >= x1 && xtemp <= x2 && ytemp >= y1 && ytemp <= y2) {

      if(first) {
	min = max = zdata_[i];
	first = false;
      }

      min = MIN(zdata_[i],min);
      max = MAX(zdata_[i],max);
      mean += (zdata_[i] - mean)/(n+1);
      ++n;

      if(zdata_[i] > 0.0) {
	if(firstGZ) {
	  minGZ = zdata_[i];
	  firstGZ = false;
	} else {
	  minGZ = MIN(zdata_[i], minGZ);
	}
      }
    }
  }

  sd = 0.0;
  n = 0;

  for(unsigned i=0; i < ndata_; i++) {
    unsigned yind = i/nx_;
    unsigned xind = i - yind*nx_;

    double xtemp = xmins_ + dx_/2 + xind*dx_;
    double ytemp = ymins_ + dy_/2 + yind*dy_;

    if(xtemp >= x1 && xtemp <= x2 && ytemp >= y1 && ytemp <= y2) {
      sd += ((zdata_[i] - mean)*(zdata_[i]-mean) - sd)/(n+1);
      ++n;
    }
  }

  if(n > 1)
    sd = sqrt(sd*n/(n-1));
  else 
    sd = 0.0f;

  fprintf(stdout, "\n\n\t\tmean\t=\t%g\n\t\tsd\t=\t%g\n\t\tmin\t=\t%g\n\t\tmin > 0\t=\t%g\n\t\tmax\t=\t%g\n\t\tnpts\t=\t%d\n", mean, sd, min, minGZ, max, n);
}

float Trans::valNearestToPoint(float x, float y)
{
  bool first=true;
  float dist,distmin,xtmp,ytmp,val;
  int xmin,ymin;

  rectify(xmins_, xmaxs_);
  rectify(ymins_, ymaxs_);

  for(unsigned i=0; i < ndata_; i++) {
    unsigned yind = i/nx_;
    unsigned xind = i - yind*nx_;

    xtmp = xmins_ + dx_/2 + xind*dx_;
    ytmp = ymins_ + dy_/2 + yind*dy_;

    dist = sqrt((xtmp-x)*(xtmp-x) + (ytmp-y)*(ytmp-y));

    if(dist < distmin || first) {
      xmin = xind;
      ymin = yind;
      distmin = dist;
      first = false;
      val = zdata_[i];
    }
  }

  return val;
}

float Trans::convolveAroundPoint(float x, float y, float* data)
{
  int convMaskInPixels = 2;
  double convSigInPixels = convMaskInPixels/2.35;

  int iXNear;
  int iYNear;

  //  COUT("Convol reversx = " << reverseX_ << " y = " << reverseY_);

  iXNear = ceil((x - xmins_) / dx_);
  iYNear = ceil((y - ymins_) / dy_);

  if(iXNear < 0)
    iXNear = 0;
  if(iXNear > nx_-1)
    iXNear = nx_+1;

  if(iYNear < 0)
    iYNear = 0;
  if(iYNear > ny_-1)
    iYNear = ny_+1;

  int iXMin = iXNear > convMaskInPixels ? iXNear-convMaskInPixels : 0;
  int iXMax = iXNear < nx_-1-convMaskInPixels ? iXNear+convMaskInPixels : nx_-1;

  int iYMin = iYNear > convMaskInPixels ? iYNear-convMaskInPixels : 0;
  int iYMax = iYNear < ny_-1-convMaskInPixels ? iYNear+convMaskInPixels : ny_-1;

  double val=0.0;
  double wt=0.0;
  double wtSum=0.0;

  double sigx = convSigInPixels * fabs(dx_);
  double sigy = convSigInPixels * fabs(dy_);

  for(unsigned iX=iXMin; iX < iXMax; iX++) {
    for(unsigned iY=iYMin; iY < iYMax; iY++) {

      unsigned ind = iY * nx_ + iX;

      double xCurr = (xmins_ + dx_/2 + iX*dx_);
      double yCurr = (ymins_ + dy_/2 + iY*dy_);

      double dist2 = (x - xCurr)*(x - xCurr) + (y - yCurr)*(y - yCurr);
      
      wt = exp(-dist2/(2*sigx*sigy));
      val += (data[ind] - val)*wt / (wtSum + wt);
      wtSum += wt;
    }
  }

  return val;
}


