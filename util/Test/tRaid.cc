#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Mass.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "startyear",            "2014",  "d", "Start year"},
  { "endyear",              "2018",  "d", "Start year"},
  { "gbperday",             "120",   "d", "GB per day"},
  { "tbperdisk",            "6",     "d", "TB per disk"},
  { "pricepertb",           "80.0",   "d", "$ per TB"},
  { "raidlevel",            "10",    "i", "RAID level"},
  { "diskfailrate",         "0.025", "d", "Rate of disk failures per month (estimated from CARMA)"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

double predictedCostPerGb(double year)
{
  return pow(10, -0.2505 * (year - 1980) + 6.304);
}

int nDiskRaid(int raidLevel, int tbPerDisk, unsigned tbTotal)
{
  unsigned nDisk = ceil(tbTotal/tbPerDisk);

  switch(raidLevel) {
  case 0:
    return nDisk;
  case 5:
    if(nDisk < 3)
      return 3;
    else
      return nDisk + 1; // RAID5 requires data striped across N disks, plus 1 for parity
    break;
  case 10:
    return 2*nDisk;     // RAID10 requires at least one mirror of the complete data set
    break;
  default:
    ThrowError("I don't curently support other RAID options");
    break;
  }
}

double calculateAnnualCostForStorage(double gbPerDay, double tbPerDisk, double pricePerTb, unsigned raidLevel, double failratepermonth)
{
  double costRaw, costAncillary;

  double tbTotalRaw = gbPerDay * 365.0 / 1000 * 5/3.5; // multiply by ratio of gzip to xz compression
  unsigned nDiskRaw = nDiskRaid(raidLevel, tbPerDisk, tbTotalRaw);

  double tbTotalAncillary = tbTotalRaw * 2.0;
  unsigned nDiskAncillary = nDiskRaid(0, tbPerDisk, tbTotalRaw);

  unsigned nDisk = nDiskRaw + nDiskAncillary;
  nDisk += ceil(failratepermonth * nDisk * 12);

  double tbTotal = nDisk * tbPerDisk;

  COUT("tbRaw = " << tbTotalRaw << " (gz) " << tbTotalRaw * 3.5/5 << " (xz) " << " tbTotalAncillary = " << tbTotalAncillary << " tbTotal = " << tbTotal);

  return tbTotal * pricePerTb;
}

int Program::main()
{
  double gbPerDay    = Program::getDoubleParameter("gbperday");
  double tbPerDisk   = Program::getDoubleParameter("tbperdisk");
  double pricePerTb  = Program::getDoubleParameter("pricepertb");
  unsigned raidLevel = Program::getIntegerParameter("raidlevel");
  unsigned startYear = Program::getIntegerParameter("startyear");
  unsigned endYear   = Program::getIntegerParameter("endyear");
  double failRate    = Program::getDoubleParameter("diskfailrate");

  for(unsigned iYear = startYear; iYear < endYear; iYear++) {
    COUT("Cost for year: " << iYear << " is " << calculateAnnualCostForStorage(gbPerDay, tbPerDisk, pricePerTb, raidLevel, failRate));
  }

  return 0;
}
