#include "gcp/models/Generic1DGaussian.h"

#include "gcp/fftutil/DataSetType.h"

using namespace std;

using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Generic1DGaussian::Generic1DGaussian() 
{
  dataSetType_ = DataSetType::DATASET_1D;

  addComponent(norm_);
  addComponent(mean_);
  addComponent(sigma_);

  addComponentName(norm_, "norm",   "The amplitude of the gaussian (normalization)");
  addComponentName(mean_, "mean",   "The mean of the gaussian");
  addComponentName(sigma_, "sigma", "Gaussian sigma");

  initializeComponentsToFixed();
}

/**.......................................................................
 * Destructor.
 */
Generic1DGaussian::~Generic1DGaussian() {}

void Generic1DGaussian::setNorm(double norm)
{
  norm_ = norm;
}

void Generic1DGaussian::setMean(double mean)
{
  mean_ = mean;
}

void Generic1DGaussian::setSigma(double sigma)
{
  sigma_ = sigma;
}

double Generic1DGaussian::eval(double x)
{
  double dx = x - mean_.value();
  double sig = sigma_.value();

  return norm_.value() * exp(-dx*dx/(2*sig*sig));
}
