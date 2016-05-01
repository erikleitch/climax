#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Distribution.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Sampler.h"
#include "gcp/pgutil/PgUtil.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "rad",        "0", "d", "Radians"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

int Program::main()
{
  std::vector<double> noise       = Sampler::generateGaussianSamples(0.1, 100000);
  std::vector<double> signal      = Sampler::generateGaussianSamples(0.1, 100000);
  std::vector<double> signalnoise = Sampler::generateGaussianSamples(0.1, 100000);

  std::vector<double> samples(signal.size() + noise.size());

  for(unsigned i=0; i < signal.size(); i++)
    signal[i] += 1.0;

  for(unsigned i=0; i < noise.size(); i++)
    samples[i] = noise[i];

  for(unsigned i=0; i < signal.size(); i++)
    samples[i + noise.size()-1] = (signal[i] + signalnoise[i] - 0.7*signal[i]);

  PgUtil::histogram(&samples, 100);


  return 0;
}
