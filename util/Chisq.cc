#include "gcp/util/Chisq.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Sampler.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Chisq::Chisq() 
{
  initialize();
}

/**.......................................................................
 * Destructor.
 */
Chisq::~Chisq() {}

void Chisq::initialize()
{
  reducedChisq_ = 0.0;
  nChisq_       = 0;
}

void Chisq::operator+=(double chisqVal)
{
  reducedChisq_ += (chisqVal - reducedChisq_) / (nChisq_ + 1);
  nChisq_       += 1;
}

/**.......................................................................
 * Coadd two chi-squares
 */
void Chisq::operator+=(const Chisq& chisq)
{
  operator+=((Chisq&) chisq);
}

void Chisq::operator+=(Chisq& chisq)
{
  double n1 = (double)nChisq_;
  double n2 = (double)chisq.nChisq_;

  // Work with reduced chi-squares to avoid roundoff errors

  double r1 = reducedChisq_;
  double r2 = chisq.reducedChisq_;

  reducedChisq_ = n1/(n1 + n2) * r1 + n2/(n1+n2) * r2;
  nChisq_       = nChisq_ + chisq.nChisq_; 
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, const Chisq& chisq)
{
  return operator<<(os, (Chisq&)chisq);
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, Chisq& chisq)
{
  os << "Reduced chi-square = " << chisq.reducedChisq() << " nChisq = " << chisq.nChisq();
  return os;
}

unsigned Chisq::nChisq()
{
  return nChisq_;
}

double Chisq::chisq()
{
  return (nChisq_ * reducedChisq_);
}

double Chisq::reducedChisq()
{
  return reducedChisq_;
}

double Chisq::probabilityToExceed()
{
  double a = (double)(nChisq())/2;
  double x = chisq()/2;

  double val;
  Sampler::gammp(a, x, &val);

  return 1.0 - val;
}

void Chisq::setChisq(double chisq, unsigned nDof)
{
  reducedChisq_ = chisq / nDof;
  nChisq_       = nDof;
}
