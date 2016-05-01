#include <iostream>
#include <iomanip>

#include <iostream>
#include <fstream>

#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/CurlUtils.h"
#include "gcp/util/Exception.h"
#include "gcp/util/String.h"
#include "gcp/util/Sort.h"
#include "gcp/util/Date.h"
#include "gcp/util/DirList.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "path",      "100", "s", "Source to fetch"},
  { "start",     "28-Aug-2014:07:00:00", "s", "Start UTC"},
  { "stop",      "28-Aug-2014:08:00:00", "s", "Stop  UTC"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

double mjd(std::string scriptDate);

void parseFile(std::string fileName, double mjdStart, double mjdStop, 
	       std::vector<string>& dateVec, std::map<std::string, std::string>& dateMap);

int Program::main()
{
  DirList dirList(Program::getStringParameter("path"), false);
  Date start, stop;

  start.setToDateAndTime(Program::getStringParameter("start"));
  stop.setToDateAndTime( Program::getStringParameter("stop"));

  std::list<DirList::DirEnt> entries = dirList.getFiles();

  std::stringstream dateOs;
  std::vector<string> dateVec;
  std::vector<string> dateVecSorted;
  std::map<std::string, std::string> dateMap;

  std::ostringstream os;
  for(std::list<DirList::DirEnt>::iterator iter=entries.begin(); iter != entries.end(); iter++) {
    os.str("");
    os << iter->path_ << "/" << iter->name_;
    //    COUT("path = " << os.str());
    parseFile(os.str(), start.mjd(), stop.mjd(), dateVec, dateMap);
  }

  //------------------------------------------------------------
  // Done parsing files.  Sort the lines we found
  //------------------------------------------------------------

  dateVecSorted = Sort::sort(dateVec);

  for(unsigned i=0; i < dateVec.size(); i++) {
    COUT(dateMap[dateVecSorted[i]]);
  }

}

void parseFile(std::string fileName, double mjdStart, double mjdStop, 
	       std::vector<string>& dateVec, std::map<std::string, std::string>& dateMap)
{
  std::ifstream inputFile(fileName.c_str());
  std::ostringstream dateOs;

  std::string s;
  String str, mon, day, time, year;
  while(getline(inputFile, s)) {

    //    COUT("Read line: " << s);
    String str(s);

    //------------------------------------------------------------
    // Parse the date
    //------------------------------------------------------------

    mon  = str.findNextString();
    str.advanceToNextNonWhitespaceChar();

    day  = str.findNextString();
    str.advanceToNextNonWhitespaceChar();

    time = str.findNextString();
    str.advanceToNextNonWhitespaceChar();

    year = str.findNextInstanceOf("{", true, "}", true, true);

    //------------------------------------------------------------
    // Reformat it to something that we can parse
    //------------------------------------------------------------

    dateOs.str("");
    dateOs << day << " " << mon << " " << year << " " << time;

    //    COUT("date = " << dateOs.str());
    Date date;
    date.setToDateAndTime(dateOs.str());

    //------------------------------------------------------------
    // Check if the mjd is between the dates we want
    //------------------------------------------------------------

    double mjd = date.mjd();

    if(mjd >= mjdStart && mjd <= mjdStop) {

      double mjdD = date.mjd();

      dateOs.str("");
      dateOs << std::setw(15) << std::setprecision(13) << std::left << std::setfill('0') << mjdD;
      
      while(dateMap.find(dateOs.str()) != dateMap.end()) {
	mjdD += 1e-8;
	dateOs.str("");
	dateOs << std::setw(15) << std::setprecision(12) << mjdD;
      }
      
      dateMap[dateOs.str()] = str.str();
      dateVec.push_back(dateOs.str());
    }
  }

  inputFile.close();
}

