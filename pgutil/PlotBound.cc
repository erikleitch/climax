#include "gcp/pgutil/PlotBound.h"
#include <cmath>

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
PlotBound::PlotBound() {}

/**.......................................................................
 * Destructor.
 */
PlotBound::~PlotBound() {}

void PlotBound::expandPerc(double perc)
{
  xrng_.expandPerc(perc);
  yrng_.expandPerc(perc);
}

void PlotBound::expandAbs(double val)
{
  xrng_.expandAbs(val);
  yrng_.expandAbs(val);
}

void PlotBound::setTo(double xmin, double xmax, double ymin, double ymax,
		      bool reverseX, bool reverseY)
{
  xrng_.setTo(xmin, xmax, reverseX);
  yrng_.setTo(ymin, ymax, reverseY);
}

double PlotBound::maxRadialExtent()
{
  double dx = xrng_.absMax() - xrng_.midPoint();
  double dy = yrng_.absMax() - yrng_.midPoint();

  return sqrt(dx*dx + dy*dy);
}

bool PlotBound::contains(double x, double y)
{
  return containsX(x) && containsY(y);
}

bool PlotBound::containsX(double x)
{
  return xrng_.contains(x);
}

bool PlotBound::containsY(double y)
{
  return yrng_.contains(y);
}

double PlotBound::absXmin()
{
  return xrng_.absMin();
}

double PlotBound::absXmax()
{
  return xrng_.absMax();
}

double PlotBound::absYmin()
{
  return yrng_.absMin();
}

double PlotBound::absYmax()
{
  return yrng_.absMax();
}

double PlotBound::plotXmin()
{
  return xrng_.plotMin();
}

double PlotBound::plotXmax()
{
  return xrng_.plotMax();
}

double PlotBound::plotYmin()
{
  return yrng_.plotMin();
}

double PlotBound::plotYmax()
{
  return yrng_.plotMax();
}

double PlotBound::dataXmin()
{
  return xrng_.dataMin();
}

double PlotBound::dataXmax()
{
  return xrng_.dataMax();
}

double PlotBound::dataYmin()
{
  return yrng_.dataMin();
}

double PlotBound::dataYmax()
{
  return yrng_.dataMax();
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, PlotBound& b)
{
  os << "(" << b.absXmin() << ", " << b.absXmax() << ") (" << b.absYmin() << ", " << b.absYmax() << ")";
  return os;
}

double PlotBound::xRange()
{
  return xrng_.range();
}

double PlotBound::yRange()
{
  return yrng_.range();
}

double PlotBound::range()
{
  return sqrt(xrng_.range() * xrng_.range() + yrng_.range() * yrng_.range());
}
