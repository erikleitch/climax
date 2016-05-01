#include "gcp/pgutil/PlotAxisRange.h"
#include "gcp/util/Exception.h"

#include <cmath>

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
PlotAxisRange::PlotAxisRange() 
{
  reverse_ = false;
  isSpecified_ = false;
}

/**.......................................................................
 * Destructor.
 */
PlotAxisRange::~PlotAxisRange() {}

void PlotAxisRange::setMinMax(double min, double max)
{
  //------------------------------------------------------------
  // Always store the min max in true (absolute) order
  //------------------------------------------------------------

  dataMin_ = min;
  dataMax_ = max;

  min_ = (min < max ? min : max);
  max_ = (min > max ? min : max);

  isSpecified_ = false;
}

void PlotAxisRange::setTo(double min, double max)
{
  setMinMax(min, max);
}

/**.......................................................................
 * Set the range of this axis.  The min/max specifies the true data
 * range, in increasing order of data index.  
 *
 * The parameter 'reverse' specifies whether we should display the
 * data with the reverse sense to the natural order.  Thus if min <
 * max, the data will be displayed in increasing order if reverse = false.
 *
 * But if min > max, the data will be displayed in decreasing order if
 * reverse = false;
 */
void PlotAxisRange::setTo(double min, double max, bool reverse)
{
  setMinMax(min, max);
  reverse_ = reverse ^ (min > max);
}

void PlotAxisRange::setTo(double min, double max, unsigned ndata, float* data)
{
  setMinMax(min, max);

  bool autoscale = (min_ == max_);

  if(autoscale) {
    min_ = max_ = data[0];

    for(unsigned i=0; i < ndata; i++) {
      if(!isnan(data[i])) {
	min_ = min_ < data[i] ? min_ : data[i];
	max_ = max_ > data[i] ? max_ : data[i];
      }
    }
  
    if(min_ == max_) {
      min_ -= 0.1*min_;
      max_ += 0.1*max_;
    }
  }
}

void PlotAxisRange::setTo(double min, double max, unsigned ndata, float* data, bool reverse)
{
  setMinMax(min, max);
  reverse_ = reverse;

  bool autoscale = (min_ == max_);

  if(autoscale) {
    min_ = max_ = data[0];

    for(unsigned i=0; i < ndata; i++) {
      if(!isnan(data[i])) {
	min_ = min_ < data[i] ? min_ : data[i];
	max_ = max_ > data[i] ? max_ : data[i];
      }
    }
  
    if(min_ == max_) {
      min_ -= 0.1*min_;
      max_ += 0.1*max_;
    }
  }
}

void PlotAxisRange::expandPerc(double perc)
{
  double rng = range();

  max_ += rng*perc;
  min_ -= rng*perc;
}

void PlotAxisRange::expandAbs(double val)
{
  max_ += val;
  min_ -= val;
}


double PlotAxisRange::midPoint()
{
  return  (min_ + max_)/2;
}

double PlotAxisRange::range()
{
  return  (max_ - min_);
}


double PlotAxisRange::absMin()
{
  return min_;
}

double PlotAxisRange::absMax()
{
  return max_;
}

double PlotAxisRange::dataMin()
{
  return dataMin_;
}

double PlotAxisRange::dataMax()
{
  return dataMax_;
}

double PlotAxisRange::plotMin()
{
  return reverse_ ? max_ : min_;
}

double PlotAxisRange::plotMax()
{
  return reverse_ ? min_ : max_;
}

bool PlotAxisRange::contains(double val)
{
  return val >= absMin() && val <= absMax();
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, PlotAxisRange& par)
{
  os << "(" << par.absMin() << ", " << par.absMax() << ")";
  return os;
}
