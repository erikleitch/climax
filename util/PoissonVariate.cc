#include "gcp/util/PoissonVariate.h"

#include "gcp/pgutil/PgUtil.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
PoissonVariate::PoissonVariate() 
{
  initialize();
}

void PoissonVariate::initialize()
{
  samplingDistribution().setType(Distribution::DIST_POISS);
  val_ = 0.0;
}

/**.......................................................................
 * Destructor.
 */
PoissonVariate::~PoissonVariate() {}

void PoissonVariate::setValue(unsigned k, double mean)
{
  val_ = (double)k;
  setMean(mean);
}

void PoissonVariate::setMean(double mean)
{
  samplingDistribution().setPoissMean(mean);
}

double PoissonVariate::mean()
{
  return samplingDistribution().getPoissMean();
}

void PoissonVariate::plotPdf(double min, double max, unsigned n)
{
  unsigned imax = (unsigned)max;
  unsigned imin = (unsigned)min;

  n = imax-imin+1;

  std::vector<double> x(n);
  std::vector<double> y(n);

  unsigned ind=0;
  for(unsigned i=imin; i <= imax; i++, ind++) {
    x[ind] = (double)i;
    y[ind] = samplingDistribution_.pdf((double)x[ind]).value();
  }

  PgUtil::linePlot(x, y);
}
