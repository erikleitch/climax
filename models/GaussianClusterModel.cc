#include "gcp/models/GaussianClusterModel.h"

using namespace std;
using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
GaussianClusterModel::GaussianClusterModel()
{
  dataSetType_ |= DataSetType::DATASET_RADIO;

  addComponent(minSigma_);
  addComponent(majSigma_);

  addComponentName(minSigma_, "minsigma");
  addComponentName(majSigma_, "majsigma");

  initializeComponentsToFixed();
}

/**.......................................................................
 * Destructor.
 */
GaussianClusterModel::~GaussianClusterModel() {}

/**.......................................................................
 * Set the major axis sigma
 */
void GaussianClusterModel::setMajSigma(gcp::util::Angle majSigma)
{
  majSigma_ = majSigma;
}

/**.......................................................................
 * Specify the major axis as a FWHM
 */
void GaussianClusterModel::setMajFwhm(gcp::util::Angle majFwhm)
{
  setMajSigma(majFwhm/sqrt(8*log(2.0)));
}

/**.......................................................................
 * Set the minor axis sigma
 */
void GaussianClusterModel::setMinSigma(gcp::util::Angle minSigma)
{
  minSigma_ = minSigma;
}

/**.......................................................................
 * Specify the minor axis as a FWHM
 */
void GaussianClusterModel::setMinFwhm(gcp::util::Angle minFwhm)
{
  setMinSigma(minFwhm/sqrt(8*log(2.0)));
}

/**.......................................................................
 * Set both sigmas
 */
void GaussianClusterModel::setSigma(gcp::util::Angle sigma)
{
  majSigma_ = sigma;
  minSigma_ = sigma;
}

/**.......................................................................
 * Set both FWHMs
 */
void GaussianClusterModel::setFwhm(gcp::util::Angle fwhm)
{
  setMinSigma(fwhm/sqrt(8*log(2.0)));
  setMajSigma(fwhm/sqrt(8*log(2.0)));
}

/**.......................................................................
 * Evaluate the dimensionless shape of this profile
 */
double GaussianClusterModel::radioEnvelope(double xRad, double yRad)
{
  double minSigRad = minSigma_.radians();
  double majSigRad = majSigma_.radians();
  double arg = 0.5 * ((xRad * xRad)/(majSigRad * majSigRad) + (yRad * yRad)/(minSigRad * minSigRad));

  return exp(-arg);
}

void GaussianClusterModel::fillSzImage(Image& image, Frequency& frequency)
{
  return Generic2DAngularModel::fillImage(DataSetType::DATASET_RADIO, image, &frequency);
}
