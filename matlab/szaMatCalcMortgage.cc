/**.......................................................................
 * MATLAB Mex file for calculating UVW, given an array of antenna
 * locations, a source declination and a list of HA
 *
 * Use like:
 *
 * d=gcpMatCalcUvw(lat, antLoc(ENU), refLoc(ENU), dec(rad), HA(hours))
 *
 */
#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"
#include "gcp/array/code/share/slalib/slalib.h"
#include "gcp/util/Coordinates.h"

#include "mex.h"
#include "matrix.h"

#include <iostream.h>
#include <math.h>

using namespace std;
using namespace gcp::util;
using namespace gcp::matlab;

/**.......................................................................
 * Entry point from the matlab environment
 */
void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  gcp::util::Logger::installStdoutPrintFn(MexHandler::stdoutPrintFn);
  gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);
  
  gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
  gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
  gcp::util::ErrHandler::installLogFn(MexHandler::logFn);
  
  MexParser parser(prhs[0]);

  if(!parser.isStruct())
    ThrowError("Input should be a struct");

  double price               = *parser.getFieldAsDouble("price");
  double down                = *parser.getFieldAsDouble("down");
  double principal           = price - down;
  double annualInterestRate  = *parser.getFieldAsDouble("mortrate");
  double monthlyInterestRate = annualInterestRate / 12;
  double npt                 = *parser.getFieldAsDouble("npt");
  double nMonth              = *parser.getFieldAsDouble("nyear") * 12;
  double nMonthEval          = *parser.getFieldAsDouble("nyeareval") * 12;
  double extra1              = *parser.getFieldAsDouble("extra1");
  double nmonthextra1        = *parser.getFieldAsDouble("nmonthextra1");
  double extra2              = *parser.getFieldAsDouble("extra2");
  double nmonthextra2        = *parser.getFieldAsDouble("nmonthextra2");
  
  //------------------------------------------------------------
  // Mortgage calculation
  //------------------------------------------------------------

  double fiducialMonthlyMortgagePayment = 
    principal * monthlyInterestRate / (1.0 - pow((1 + monthlyInterestRate), -nMonth));

  double interest;
  double fiducialMonthlyPayment;
  double currentMonthlyPayment;
  double paid=0.0;

  double annualPropertyTaxRate  = *parser.getFieldAsDouble("ptax");
  double monthlyPropertyTaxRate = annualPropertyTaxRate / 12;
  double monthlyPropertyTax = monthlyPropertyTaxRate * price;

  double taxReturn, effectiveMonthlyPayment;

  double totalRent = 0.0;
  double effectiveTotalValue = 0.0;

  double monthlySalary = *parser.getFieldAsDouble("salary");

  // The principal in our savings account

  double savingsPrincipal = *parser.getFieldAsDouble("savings");

  // Decrement by the down payment, plus money we pay in points

  savingsPrincipal -= down;

  // Each point is 1% of the principal amount

  savingsPrincipal -= npt * principal * 0.01;

  COUT("Savings principal is: " << savingsPrincipal);

  double annualSavingsInterestRate = *parser.getFieldAsDouble("saverate");
  double monthlySavingsInterestRate = annualSavingsInterestRate / 12;

  // The savings interest we made this month

  double monthlyEarnedSavingsInterest = 0.0;
  double monthlySavings;
  double annualFedTaxRate    = *parser.getFieldAsDouble("fedtax");
  double annualStateTaxRate  = *parser.getFieldAsDouble("statetax");
  double sixMonthPropertyTax = 0.0;

  //------------------------------------------------------------
  // Now calculate the 'amortization table'
  //------------------------------------------------------------

  if(nMonthEval > nMonth)
    nMonth = nMonthEval;

  int dims[2] = {1,1};
  plhs[0] = mxCreateStructArray(2,dims,0,NULL);

  double* monthlyPaymentPtr     = MexHandler::addNamedDoubleStructField(plhs[0], "monthly", 1);
  double* monthlyPropertyTaxPtr = MexHandler::addNamedDoubleStructField(plhs[0], "proptax", 1);
  double* extraPtr              = MexHandler::addNamedDoubleStructField(plhs[0], "extra",   1);
  double* totalPtr              = MexHandler::addNamedDoubleStructField(plhs[0], "total",   1);

  double* paidPtr       = MexHandler::addNamedDoubleStructField(plhs[0], "paid",       (unsigned)nMonth);
  double* principalPtr  = MexHandler::addNamedDoubleStructField(plhs[0], "principal",  (unsigned)nMonth);
  double* interestPtr   = MexHandler::addNamedDoubleStructField(plhs[0], "interest",   (unsigned)nMonth);
  double* savingsPtr    = MexHandler::addNamedDoubleStructField(plhs[0], "savings",    (unsigned)nMonth);
  double* effPaymentPtr = MexHandler::addNamedDoubleStructField(plhs[0], "effpayment", (unsigned)nMonth);

  *monthlyPaymentPtr     = fiducialMonthlyPayment;
  *monthlyPropertyTaxPtr = monthlyPropertyTaxRate * price;
  *extraPtr              = extra1 + extra2;
  *totalPtr              = fiducialMonthlyPayment + (monthlyPropertyTaxRate * price) + extra1 + extra2;

  double currentExtra;

  for(unsigned iMonth=0; iMonth < nMonth; iMonth++) {

    if(iMonth < nmonthextra1)
      currentExtra = extra1;
    else
      currentExtra = 0.0;

    if(iMonth < nmonthextra2)
      currentExtra += extra2;
    else
      currentExtra += 0.0;

    // This month's interest on the principal
 
    interest  = principal * monthlyInterestRate;

    // Calculate what our fiducial mortgage-only (+extra) payment is
    // for this month

    fiducialMonthlyPayment = fiducialMonthlyMortgagePayment + currentExtra;

    // If the total amount we still need to pay is less than our
    // fiducial monthly payment, then reduce the monthly payment
    // accordingly (ie, this is our last payment)
 
    if(principal + interest < fiducialMonthlyPayment) {
      fiducialMonthlyPayment = principal + interest;
    }

    if(iMonth == 0) {
      *monthlyPaymentPtr     = fiducialMonthlyMortgagePayment;
      *monthlyPropertyTaxPtr = monthlyPropertyTaxRate * price;
      *extraPtr              = currentExtra;
      *totalPtr              = fiducialMonthlyMortgagePayment + (monthlyPropertyTaxRate * price) + currentExtra;
    }

    // The amount by which we can reduce our monthly withholding,
    // because we can claim it back at the end of the year

    taxReturn  =   annualFedTaxRate * interest;
    taxReturn += annualStateTaxRate * interest;

    // Accumulate the property tax, to be paid in 6-month increments

    sixMonthPropertyTax += monthlyPropertyTax;

    // Our effective monthly payment for this month

    effectiveMonthlyPayment = fiducialMonthlyPayment + monthlyPropertyTax - taxReturn;

    // The remaining principal, after this month's payment
 
    principal -= (fiducialMonthlyPayment - interest);

    if(principal < 1)
      principal = 0.0;

    // The total amount we've paid, after this month's payment

    paid += effectiveMonthlyPayment;

    // The effective total value of the house = the money we put down,
    // plus what we've already paid into it, plus the remaining
    // principal

    effectiveTotalValue = down + principal + paid;

    //------------------------------------------------------------
    // Savings calculation
    //------------------------------------------------------------

    // Difference between our salary and money we have to pay is our
    // savings.  Here we increment our salary by the amount of money
    // by which we can reduce our monthly withholding, since we can
    // claim back interest.

    monthlySavings = (monthlySalary + taxReturn) - fiducialMonthlyPayment;
    monthlyEarnedSavingsInterest = savingsPrincipal * monthlySavingsInterestRate;

    savingsPrincipal += monthlySavings;
    savingsPrincipal += monthlyEarnedSavingsInterest;

    // Every six months, we have to pay property tax, which comes out
    // of our savings

    if((iMonth + 1) % 6 == 0) {
      savingsPrincipal -= sixMonthPropertyTax;
      sixMonthPropertyTax = 0.0;
    }

    paidPtr[iMonth]       = paid + down;
    principalPtr[iMonth]  = principal;
    interestPtr[iMonth]   = interest;
    savingsPtr[iMonth]    = savingsPrincipal;
    effPaymentPtr[iMonth] = effectiveMonthlyPayment;
  }

  return;
}
