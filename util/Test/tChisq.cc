#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Chisq.h"
#include "gcp/util/ChisqVariate.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Sampler.h"

#include "gcp/pgutil/PgUtil.h"

#include <vector>

using namespace std;
using namespace gcp::util;
using namespace gcp::program;


KeyTabEntry Program::keywords[] = {
  { "ntrial",     "1000", "i", "Number of trials"},
  { "nsamp",      "10", "i", "Number of samples per trial"},
  { "chisq",      "42", "d", "Chisquared value to test"},
  { "ndof",       "51", "i", "Number of dof"},
  { "seed",       "1",  "i", "Seed"},
  { "nhist",      "10", "i", "number of histogram bins"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

void pte(double chisq, unsigned ndof)
{
  ChisqVariate cv;
  cv.setChisq(chisq, ndof);
  COUT("PTE = " << cv.samplingPte() << " 1.0-PTE = " << 1.0-cv.samplingPte().value() << " pdf = " << cv.pdf() << " ln = " << cv.pdf().lnValue() << " " 
       << ((double)(ndof)/2-1)*log(chisq) - chisq/2);

  
  double val1 = (((double)(53312)/2-1)*log(54541.0) - 54541.0/2);
  double val2 = (((double)(53312)/2-1)*log(54540.0) - 54540.0/2);

  COUT("diff = " << val1 - val2);

  unsigned npt=1000;
  double xmin=0;
  double xmax=3;
  double dx = (xmax-xmin)/(npt-1);

  std::vector<double> x(npt);
  std::vector<double> y(npt);
  
  for(unsigned i=0; i < x.size(); i++) {
    x[i] = xmin + dx*i;
    y[i] = cv.samplingDistribution().pdf(x[i]).value();
  }

  PgUtil::linePlot(x, y);
}


int Program::main()
{
#if 1
  pte(Program::getDoubleParameter("chisq"), Program::getIntegerParameter("ndof"));
  return 0;
#endif

  unsigned ntrial = Program::getIntegerParameter("ntrial");
  unsigned nsamp  = Program::getIntegerParameter("nsamp");

  Sampler::seed(Program::getIntegerParameter("seed"));

  ChisqVariate chisq;
  double chisqSum = 0.0;

  std::vector<double> samples;
  std::vector<double> samples2;

  // In each trial, we will average together nSamp samples

  double meanVal=0.0, rmsVal=0.0;
  double chisqVal = 0.0;

  for(unsigned i=0; i < ntrial; i++) {

    //    nsamp = (int)Sampler::generateUniformSample(5,20);
    samples = Sampler::generateGaussianSamples(0.01, nsamp);

#if 0
#if 0
    samples2 = Sampler::generateGaussianSamples(0.01, 2*nsamp);
    
    for(unsigned iSamp=0; iSamp < samples.size(); iSamp++) {
      samples[iSamp] = (samples2[2*iSamp] + samples2[2*iSamp+1])/2;
    }
#endif
    std::vector<float>x;
    std::vector<float>y;

    Sampler::histogram(samples, Program::getIntegerParameter("nhist"), x, y);
    PgUtil::linePlot(x.size(), &x[0], &y[0]);
    std::vector<float> ygpdf(x.size());

    for(unsigned i=0; i < x.size(); i++) {
      ygpdf[i] = Sampler::gaussPdf(x[i], 0.0, 0.01/sqrt(2));
    }

    PgUtil::setOverplot(true);
    PgUtil::linePlot(x.size(), &x[0], &ygpdf[0]);
#endif

    for(unsigned iSamp=0; iSamp < samples.size(); iSamp++) {
      meanVal += samples[iSamp];
    }
    meanVal /= samples.size();

    for(unsigned iSamp=0; iSamp < samples.size(); iSamp++) {
      rmsVal += (samples[iSamp] - meanVal) * (samples[iSamp] - meanVal);
    }
    rmsVal = sqrt(rmsVal / (samples.size() - 1));

    //    COUT("mean = " << meanVal << " rms = " << rmsVal << " err = " << rmsVal / sqrt((double)samples.size()));

    double m1=0.0, m2=0.0, mean=0.0, rms=0.0, err=0.0, wtSum=0.0, wt2Sum=0.0;
    //    double wt = 1.0/(0.01*0.01);
    double wt = 1.0;

    for(unsigned iSamp=0; iSamp < nsamp; iSamp++) {
      double val = samples[iSamp];
      m1    += (val - m1) * wt / (wtSum + wt);
      wtSum += wt;
    }

    mean  = m1;
    wtSum = 0.0;

    for(unsigned iSamp=0; iSamp < nsamp; iSamp++) {
      double val = (samples[iSamp] - mean) * (samples[iSamp] - mean);
      m2     += (val - m2) * wt / (wtSum + wt);
      wtSum  += wt;
      wt2Sum += wt*wt;
    }

    if(nsamp > 1) {
      double rmsPrefac = (wtSum * wtSum) / (wtSum * wtSum - wt2Sum);
      double errPrefac = wt2Sum / (wtSum * wtSum - wt2Sum);
      rms = sqrt(rmsPrefac * m2);
      err = sqrt(errPrefac * m2);
    } else {
      err = 0.01;
    }

    err =0.01/sqrt(nsamp);

    //    COUT("mean = " << mean << " rms = " << rms << " err = " << err);

    double cont = (mean/err) * (mean/err);

#if 0
    chisq += (cont * 1.0/(nsamp * nsamp));
    chisqSum += 1.0/(nsamp*nsamp);
#else
    chisq += cont;
#endif
  }

  COUT("Chisq is now: " << chisq);
  COUT("Probability to exceed: " << chisq.samplingPte());

#if 1

#if 0
  unsigned n = 100;
#else
  unsigned n = 2*ntrial;
#endif

  std::vector<double> x(n);
  std::vector<double> y1(n);
  std::vector<double> y2(n);

#if 0
  float dx = (float)(ntrial/100);
#else
  float dx = (float)(1.0);
#endif

  for(unsigned i=0; i < n; i++) {

#if 0
    x[i] = (float)ntrial - n/2 * dx + i * dx;
#else
    x[i] = (float)i;
#endif

    COUT("x = " << x[i]);
    y1[i] = Sampler::chisqPdf(x[i], ntrial+2);
#if 0
    y2[i] = Sampler::gaussPdf(x[i], (double)ntrial, sqrt(2.0*ntrial));
#else
    y2[i] = Sampler::gaussPdf(x[i], (double)ntrial, sqrt((double)ntrial));
#endif
  }

  PgUtil::linePlot(x, y1);
  PgUtil::setOverplot(true);
  PgUtil::setTraceColor(2);
  PgUtil::linePlot(x, y2);
#endif
  return 0;
}
