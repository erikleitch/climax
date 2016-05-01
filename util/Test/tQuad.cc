#include <iostream>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/Exception.h"
#include "gcp/util/QuadraticInterpolatorPositiveAngle.h"
#include "gcp/util/QuadraticInterpolatorSignedAngle.h"
#include "gcp/util/Stats.h"
#include "gcp/util/String.h"

#include "gcp/pgutil/PgUtil.h"


using namespace std;
using namespace gcp::util;using namespace gcp::program;void Program::initializeUsage() {};
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "angle",     "0.0",                        "s", "angle to convert"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

int Program::main()
{
  QuadraticInterpolator* ra_   = new QuadraticInterpolatorPositiveAngle(0.0, true);
  QuadraticInterpolator*  dec_ = new QuadraticInterpolatorSignedAngle(0.0, true);

  QuadraticInterpolator* ra1_   = new QuadraticInterpolatorPositiveAngle(0.0, false);
  QuadraticInterpolator*  dec1_ = new QuadraticInterpolatorSignedAngle(0.0, false);

  double m1 = 56881.9871723476;
  double m2 = 56881.9878669310;
  double m3 = 56881.9885615143;

  double m4 = 56881.9892611423;
  double m5 = 56881.9899563297;

#if 0
  m1 = 0.9871723476;
  m2 = 0.9878669310;
  m3 = 0.9885615143;

  m4 = 0.9892611423;
  m5 = 0.9899563297;
#endif

  ra_->extend( m1, 2.35426474069928);
  ra_->extend( m2, 2.35426474069928);
  ra_->extend( m3, 2.35426474107692);

  dec_->extend(m1, -0.830222279849653);
  dec_->extend(m2, -0.830222279849653);
  dec_->extend(m3, -0.83022227790325);

  ra1_->extend( m1, 2.35426474069928);
  ra1_->extend( m2, 2.35426474069928);
  ra1_->extend( m3, 2.35426474107692);

  dec1_->extend(m1, -0.830222279849653);
  dec1_->extend(m2, -0.830222279849653);
  dec1_->extend(m3, -0.83022227790325);

#if 1
  double xmax = m3;
  double xmin = m1;
#else
  ra_->extend(m4, 2.35426474128128);
  ra_->extend(m5, 2.35426474148438);
  
  dec_->extend(m4, -0.830222276850304);
  dec_->extend(m5, -0.830222275804074);
  
  double xmin = m3;
  double xmax = m5;
#endif

  double xw = (xmax - xmin);
  unsigned npt=1000;

  double dx = xw/(npt-1);

  std::vector<double> xs(npt);
  std::vector<double> ys1(npt);
  std::vector<double> ys2(npt);

  std::vector<double> ys11(npt);
  std::vector<double> ys21(npt);

  //------------------------------------------------------------
  // First plot the quadratic approximation
  //------------------------------------------------------------

  for(unsigned i=0; i < npt; i++) {
    xs[i] = xmin + dx*i;
    ys1[i] = ra_->evaluate(xs[i]);
    ys2[i] = dec_->evaluate(xs[i]);

    ys11[i] = ra1_->evaluate(xs[i]);
    ys21[i] = dec1_->evaluate(xs[i]);

    //    COUT(std::setprecision(15) << "xs[i] = " << xs[i] << " ys1[i] = " << ys1[i]);

    xs[i] = i;
  }

  //  COUT("Min of dec = " << std::setprecision(15) << Stats::min(ys2));
  //  COUT("Max of dec = " << std::setprecision(15) << Stats::max(ys2));

  PgUtil::open("/xs");
  PgUtil::setTraceColor(10);
  PgUtil::linePlot(xs, ys1);
  PgUtil::setOverplot(true);
  PgUtil::setTraceColor(2);
  PgUtil::linePlot(xs, ys11);

  double ymean = Stats::mean(ys2);

  PgUtil::setYmin(ymean - 1e-6*ymean);
  PgUtil::setYmax(ymean + 1e-6*ymean);
  PgUtil::setUsedefs(true);

  PgUtil::setOverplot(false);
  PgUtil::setTraceColor(10);
  PgUtil::linePlot(xs, ys2);
  PgUtil::setOverplot(true);
  PgUtil::setTraceColor(2);
  PgUtil::linePlot(xs, ys21);

  delete ra_;
  delete dec_;

  return 0;
}
