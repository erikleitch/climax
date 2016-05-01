#include "gcp/util/Percent.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Percent::Percent() 
{
  val_ = 0.0;
}

Percent::Percent(const Max1& unit, double percent) 
{
  setPercentMax1(percent);
}

Percent::Percent(const Max100& unit, double percent) 
{
  setPercentMax100(percent);
}

/**.......................................................................
 * Destructor.
 */
Percent::~Percent() {}

void Percent::setPercentMax1(double percent)
{
  val_ = percent;
}

void Percent::setPercentMax100(double percent)
{
  val_ = percent / 100;
}

double Percent::percentMax1()
{
  return val_;
}

double Percent::percentMax100()
{
  return val_ * 100;
}

void Percent::operator=(const Percent& var)
{
  operator=((Percent&) var);
}

void Percent::operator=(Percent& var)
{
  val_ = var.val_;
}

