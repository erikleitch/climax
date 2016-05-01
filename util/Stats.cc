#include "gcp/util/Exception.h"
#include "gcp/util/Sampler.h"
#include "gcp/util/Stats.h"

#include "gcp/pgutil/PgUtil.h"

#include <cmath>

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Stats::Stats() {}

/**.......................................................................
 * Destructor.
 */
Stats::~Stats() {}

void Stats::checkIndices(unsigned n, int& iStart, int& iStop)
{
  if(iStart == 0 && iStop == 0)
    iStop = n-1;

  if(iStart > n-1 || iStop > n-1)
    ThrowError("Invalid start/stop indices: " << iStart << "/" << iStop << " for vector of length " << n);
}

/**.......................................................................
 * Calculate the min
 */
double Stats::min(std::vector<double>& vec, int iStart, int iStop)
{
  checkIndices(vec.size(), iStart, iStop);

  double minVal = vec[iStart];

  for(unsigned i=iStart; i <= iStop; i++) 
    minVal = vec[i] < minVal ? vec[i] : minVal;

  return minVal;
}

/**.......................................................................
 * Calculate the max
 */
double Stats::max(std::vector<double>& vec, int iStart, int iStop)
{
  checkIndices(vec.size(), iStart, iStop);

  double maxVal = vec[iStart];

  for(unsigned i=iStart; i <= iStop; i++)
    maxVal = vec[i] > maxVal ? vec[i] : maxVal;

  return maxVal;
}

/**.......................................................................
 * Calculate the index of the max
 */
double Stats::max(std::vector<double>& vec, unsigned& ind, int iStart, int iStop)
{
  checkIndices(vec.size(), iStart, iStop);

  double maxVal = vec[iStart];
  ind = iStart;

  for(unsigned i=iStart; i <= iStop; i++) {
    if(vec[i] > maxVal) {
      maxVal = vec[i];
      ind = i;
    }
  }

  return maxVal;
}

void Stats::minmax(std::vector<double>& vec, double& minVal, unsigned& iMin, double& maxVal, unsigned& iMax, int iStart, int iStop)
{
  checkIndices(vec.size(), iStart, iStop);

  minVal = maxVal = vec[iStart];
  iMin = iMax = iStart;

  for(unsigned i=iStart; i <= iStop; i++) {

    if(vec[i] > maxVal) {
      maxVal = vec[i];
      iMax = i;
    }

    if(vec[i] < minVal) {
      minVal = vec[i];
      iMin = i;
    }

  }
}

/**.......................................................................
 * Calculate the mean
 */
double Stats::mean(std::vector<unsigned>& vec, std::vector<unsigned>* multiplicity, int iStart, int iStop)
{
  std::vector<double> dvec(vec.size());

  for(unsigned i=0; i < vec.size(); i++)
    dvec[i] = (double)vec[i];

  return mean(dvec, multiplicity, iStart, iStop);
}

/**.......................................................................
 * Calculate the mean
 */
double Stats::mean(std::vector<double>& vec, std::vector<unsigned>* multiplicity, int iStart, int iStop)
{
  checkIndices(vec.size(), iStart, iStop);

  double val = 0.0;
  double wt = 1.0, wtSum = 0.0;
  for(unsigned i=iStart; i <= iStop; i++) {
    wt = (multiplicity ? multiplicity->at(i) : 1.0);
    wtSum += wt;
    val += wt*(vec[i] - val)/wtSum;
  }
  return val;
}

/**.......................................................................
 * Calculate the sum
 */
double Stats::sum(std::vector<unsigned>& vec, std::vector<unsigned>* multiplicity, int iStart, int iStop)
{
  checkIndices(vec.size(), iStart, iStop);

  double val = 0.0;
  for(unsigned i=iStart; i <= iStop; i++)
    val += vec[i];

  return val;
}

/**.......................................................................
 * Calculate the sum
 */
double Stats::sum(std::vector<double>& vec, std::vector<unsigned>* multiplicity, int iStart, int iStop)
{
  checkIndices(vec.size(), iStart, iStop);

  double val = 0.0;
  double wt = 1.0, wtSum = 0.0;
  for(unsigned i=iStart; i <= iStop; i++) {
    wt = (multiplicity ? multiplicity->at(i) : 1.0);
    wtSum += wt;
    val += wt*(vec[i] - val)/wtSum;
  }

  return val*wtSum;
}

/**.......................................................................
 * Calculate the rms
 */
double Stats::rms(std::vector<double>& vec, std::vector<unsigned>* multiplicity, int iStart, int iStop)
{
  return sqrt(variance(vec, multiplicity, iStart, iStop));
}

double Stats::rms(std::vector<unsigned>& vec, std::vector<unsigned>* multiplicity, int iStart, int iStop)
{
  return sqrt(variance(vec, multiplicity, iStart, iStop));
}

double Stats::variance(std::vector<double>& vec, std::vector<unsigned>* multiplicity, int iStart, int iStop)
{
  checkIndices(vec.size(), iStart, iStop);

  if((iStop - iStart + 1) < 2)
    return 0.0;

  // Get the weighted mean

  double mn = mean(vec, multiplicity, iStart, iStop);

  // Now accumulate the weighted average of the variance

  double val;
  double wt = 1.0;
  double sum = 0.0, wtSum = 0.0, wtSum2 = 0.0;

  for(unsigned i=iStart; i <= iStop; i++) {
    val = vec[i] - mn;
    val *= val;

    wt = (multiplicity ? multiplicity->at(i) : 1.0);

    // Accumulate the sum

    sum += wt*val;

    // And increment the sums we need to convert the sum to a variance

    wtSum  += wt;
    wtSum2 += wt*wt;
  }

  return sum * wtSum2 / (wtSum*wtSum - wtSum2);
}

/**.......................................................................
 * Calculate the variance 
 */
double Stats::variance(std::vector<unsigned>& vec, std::vector<unsigned>* multiplicity, int iStart, int iStop)
{
  checkIndices(vec.size(), iStart, iStop);

  if((iStop - iStart + 1) < 2)
    return 0.0;

  // Get the weighted mean

  double mn = mean(vec, multiplicity, iStart, iStop);

  // Now accumulate the weighted average of the variance

  double val;
  double wt = 1.0;
  double sum = 0.0, wtSum = 0.0, wtSum2 = 0.0;

  for(unsigned i=iStart; i <= iStop; i++) {
    val = vec[i] - mn;
    val *= val;

    wt = (multiplicity ? multiplicity->at(i) : 1.0);

    // Accumulate the sum

    sum += wt*val;

    // And increment the sums we need to convert the sum to a variance

    wtSum  += wt;
    wtSum2 += wt*wt;
  }

  return sum * wtSum2 / (wtSum*wtSum - wtSum2);
}

/**.......................................................................
 * Calculate the mode (most likely value)
 */
double Stats::mode(std::vector<double>& vec, unsigned nbin, std::vector<unsigned>* multiplicity, int iStart, int iStop)
{
  std::vector<double> hist(nbin);
  std::vector<double> x(nbin);

  histogram(vec, nbin, x, hist, multiplicity, iStart, iStop);

  unsigned indMax;
  (void)max(hist, indMax, iStart, iStop);

  return x[indMax];
}

/**.......................................................................
 * Return a histogram of the input data vector
 */
void Stats::histogram(std::vector<double>& vec, unsigned nbin, std::vector<double>& x, std::vector<double>& y, 
		      std::vector<unsigned>* multiplicity, int iStart, int iStop)
{
  checkIndices(vec.size(), iStart, iStop);

  double minVal = min(vec, iStart, iStop);
  double maxVal = max(vec, iStart, iStop);
  double midVal = (minVal + maxVal)/2;

  //------------------------------------------------------------
  // Check for finite data
  //------------------------------------------------------------

  if(!isfinite(minVal) || !isfinite(maxVal)) {
    ThrowError("Data not finite");
  }

  //------------------------------------------------------------
  // And construct the histogram.
  //------------------------------------------------------------

  double dx = (maxVal-minVal)/(nbin-1);

  if(fabs(dx/midVal) < 1e-12) {
    minVal = 0.9*minVal;
    maxVal = 1.1*maxVal;
    dx = (maxVal-minVal)/(nbin-1);
  }

  x.resize(nbin);
  y.resize(nbin);

  for(unsigned i=0; i < nbin; i++) {
    y[i] = 0.0;
    x[i] = minVal + i*dx;
  }
  
  //------------------------------------------------------------
  // Construct the histogram, optionally multiplying by a multiplicity
  // if provided
  //------------------------------------------------------------

  for(unsigned i=iStart; i <= iStop; i++) {
    int ind = (int)floor((vec[i]-minVal)/dx + 0.5);

    if(ind < 0 || ind > nbin-1) {
      y[0]   += 1.0 * (multiplicity ? multiplicity->at(i) : 1);
    } else {
      y[ind] += 1.0 * (multiplicity ? multiplicity->at(i) : 1);
    }
  }
}

/**.......................................................................
 * Calculate the N-sigma confidence interval about the specified point
 */
void Stats::confidenceIntervalN(std::vector<double>& vec, unsigned nbin, double nSigma, double refVal, 
				double& lowVal, double& highVal, std::vector<unsigned>& multiplicity, 
				int iStart, int iStop)
{
  return confidenceIntervalN(vec, nbin, nSigma, refVal, lowVal, highVal, &multiplicity, iStart, iStop);
}

void Stats::confidenceIntervalN(std::vector<double>& vec, unsigned nbin, double nSigma, double refVal, 
				double& lowVal, double& highVal, std::vector<unsigned>* multiplicity, 
				int iStart, int iStop)
{
  checkIndices(vec.size(), iStart, iStop);

  double perc = 1.0 - 2*Sampler::gaussCdf(-nSigma, 0.0, 1.0);

  std::vector<double> x, y;

  histogram(vec, nbin, x, y, multiplicity, iStart, iStop);

  double hSum = multiplicity ? sum(*multiplicity) : (double)vec.size();

  double minVal = x[iStart];
  double dx = x[iStart+1] - x[iStart];

  //------------------------------------------------------------
  // Get the index of the reference value and add the contribution of
  // this bin to the total
  //------------------------------------------------------------

  int indRef = (int)floor((refVal - minVal)/dx + 0.5);

  if(indRef < 0 || indRef > nbin-1) {
    ThrowError("Invalid index");
  }

  double frac = y[indRef] / hSum;

  //------------------------------------------------------------
  // Now iterate symmetrically from the reference value until perc is
  // enclosed
  //------------------------------------------------------------

  bool stop = false, lowStop=false, highStop=false;
  int iHigh=indRef, iLow=indRef;

  for(int i=1; frac < perc && !stop; i++) {

    iHigh = indRef + i;

    if(iHigh < nbin) {
      frac += y[iHigh]/hSum;
    } else {
      iHigh = nbin - 1;
      highStop = true;
    }

    iLow = indRef - i;

    if(iLow >= 0) {
      frac += y[iLow]/hSum;
    } else {
      iLow = 0;
      lowStop = true;
    }
    
    stop = lowStop && highStop;
  }

  lowVal  = x[iLow];
  highVal = x[iHigh];
}
