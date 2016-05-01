#include <iostream>
#include <iomanip>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/CurlUtils.h"
#include "gcp/util/Exception.h"
#include "gcp/util/String.h"
#include "gcp/util/Sort.h"
#include "gcp/util/Date.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "nmess",        "100", "s", "Source to fetch"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

double mjd(std::string scriptDate);

int Program::main()
{
  std::ostringstream os;
  os << "http://logstash.ovro.pvt:9200/logstash-2014.08.28/_search?q=host:boot.carma.pvt&size=" << Program::getIntegerParameter("nmess") << "&q=message:SCRIPT";

  CurlUtils cu;
  String str = cu.getUrl(os.str());

  std::stringstream dateOs;
  std::vector<string> dateStr;
  std::vector<string> dateStrSorted;
  std::map<std::string, std::string> dateMap;

  String mess;
  String date;
  do {
    mess = str.findNextInstanceOf("_source\":{\"message\":\"", true, "\",\"@", true, true);

    if(!mess.contains("SCRIPT"))
      continue;

    if(!mess.isEmpty()) {
      //      COUT(mess);
      date = str.findNextInstanceOf("carma_utc_datetime\":\"", true, "\"", true, true);

      double mjdD = mjd(date.str());
      dateOs.str("");
      dateOs << std::setw(15) << std::setprecision(13) << std::left << std::setfill('0') << mjdD;

      while(dateMap.find(dateOs.str()) != dateMap.end()) {
	mjdD += 1e-8;
	dateOs.str("");
	dateOs << std::setw(15) << std::setprecision(12) << mjdD;
      }

      dateMap[dateOs.str()] = mess.str();
      dateStr.push_back(dateOs.str());
    }
  } while(!mess.isEmpty());

  dateStrSorted = Sort::sort(dateStr);

  for(unsigned i=0; i < dateStr.size(); i++) {
    COUT(dateMap[dateStrSorted[i]]);
  }

  return 0;
}

double mjd(std::string scriptDate)
{
  String dateStr = scriptDate;

  std::ostringstream dateOs;

  String mon  = dateStr.findNextString();
  String day  = dateStr.findNextString();
  String time = dateStr.findNextString();
  String year = dateStr.findNextString();

  dateOs << day << " " << mon << " " << year << " " << time;

  Date date;
  date.setToDateAndTime(dateOs.str());

  return date.mjd();
}
