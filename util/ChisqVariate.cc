#include "gcp/util/ChisqVariate.h"
#include "gcp/util/Sampler.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ChisqVariate::ChisqVariate() 
{
  initialize();
}

/**.......................................................................
 * Destructor.
 */
ChisqVariate::~ChisqVariate() {}

void ChisqVariate::initialize()
{
  samplingDistribution().setType(Distribution::DIST_CHISQ);
  val_ = 0.0;
  setNdof(0);
}

void ChisqVariate::setChisq(double chisq, unsigned nDof)
{
  val_ = chisq / nDof;
  setNdof(nDof);
}

void ChisqVariate::setNdof(unsigned nDof)
{
  samplingDistribution().setChisqNdof(nDof);
}

unsigned ChisqVariate::nDof()
{
  return samplingDistribution().chisqNdof();
}

double ChisqVariate::chisq()
{
  return val_ * nDof();
}

double ChisqVariate::reducedChisq()
{
  return val_;
}

/**.......................................................................
 * Coadd two chi-squares
 */
void ChisqVariate::operator+=(const ChisqVariate& chisq)
{
  operator+=((ChisqVariate&) chisq);
}

void ChisqVariate::operator+=(ChisqVariate& chisq)
{
  double n1 = (double)nDof();
  double n2 = (double)chisq.nDof();

  // Work with reduced chi-squares to avoid roundoff errors

  double r1 = val_;
  double r2 = chisq.val_;

  // Compute the new reduced chisq

  val_ = (n1*r1 + n2*r2)/(n1 + n2);

  setNdof((unsigned)(n1+n2));
}

bool ChisqVariate::operator<(const ChisqVariate& chisq)
{
  return operator<((ChisqVariate&) chisq);
}

bool ChisqVariate::operator<(ChisqVariate& chisq)
{
  return val_ < chisq.val_;
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, const ChisqVariate& chisq)
{
  return operator<<(os, (ChisqVariate&)chisq);
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, ChisqVariate& chisq)
{
  os << "reduced chi-squared = " << chisq.reducedChisq() << " nDof = " << chisq.nDof() << " PTE = " << setprecision(8) << chisq.pte();
  return os;
}

/**.......................................................................
 * Return the pdf of this value of chi-squared
 */
Probability ChisqVariate::pdf()
{
  return Variate::pdf();
}

/**.......................................................................
 * Return the equivalent Gaussian likelihood of this value of chi-squared
 */
Probability ChisqVariate::likelihood()
{
  Probability prob;
  prob.setLnValue(-chisq()/2);
  return prob;
}

void ChisqVariate::operator=(const ChisqVariate& var)
{
  operator=((ChisqVariate&) var);
}

void ChisqVariate::operator=(ChisqVariate& var)
{
  val_ = var.val_;
  setNdof(var.nDof());
}
