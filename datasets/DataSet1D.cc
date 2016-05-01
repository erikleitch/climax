#include "gcp/datasets/DataSet1D.h"

#include "gcp/fftutil/DataSetType.h"
#include "gcp/fftutil/Generic1DModel.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/util/Stats.h"
#include "gcp/util/String.h"

#include <fstream>

using namespace std;

using namespace gcp::datasets;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
DataSet1D::DataSet1D() 
{
  addParameter("n",      DataType::UINT, "If generating data, the length of the array to generate");

  setParameter("relativex", "false");
  setParameter("n", "100");

  dataSetType_ = DataSetType::DATASET_1D;
  displayInitialized_ = false;
}

/**.......................................................................
 * Destructor.
 */
DataSet1D::~DataSet1D() {}

/**.......................................................................
 * Initialize this data set from a file
 */
void DataSet1D::loadData(bool simulate)
{
  if(simulate && !getParameter("file", false)->data_.hasValue()) {
    setupForSim();
  } else {
    loadData(getStringVal("file"));
  }
}

void DataSet1D::setupForSim()
{
  unsigned n = getUintVal("n");

  x_.resize(n);
  y_.resize(n);
  error_.resize(n);

  // Resize internal arrays used for managing model components

  compositeModel_.resize(n);

  // Initialize the x array to integer 1-n

  for(unsigned i=0; i < n; i++)
    x_[i] = (double)i;
}

void DataSet1D::loadData(std::string fileName)
{
  bool relative = getBoolVal("relativex");

  static String line;

  std::ifstream fin;
  fin.open(fileName.c_str(), ios::in);

  if(!fin) {
    ThrowColorError(std::endl << "Unable to open file: " << fileName, "red");
  }

  unsigned nLine=0;
  String xStr, yStr, eStr;
  double firstXVal, xVal, yVal;
  bool first = true;

  while(!fin.eof()) {

    ++nLine;
    line.initialize();
    getline(fin, line.str());

    line.advanceToNextNonWhitespaceChar();

    // Ignore empty lines

    if(line.isEmpty()) 
      continue;      

    // Ignore empty lines

    if(line[0] == '/')
      continue;

    // We expect three columns of numbers, separated by any of a
    // number of allowable separators

    xStr = line.findNextStringSeparatedByChars(" \t,", true);

    if(xStr.isEmpty()) 
      ThrowError("Unable to read x-value on line: " << nLine << ". Perhaps you are using an invalid separator?");

    yStr = line.findNextStringSeparatedByChars(" \t,", true);

    if(yStr.isEmpty()) 
      ThrowError("Unable to read y-value on line: " << nLine << ". Perhaps you are using an invalid separator?");

    eStr = line.findNextStringSeparatedByChars(" \t,", true);

    if(eStr.isEmpty()) 
      ThrowError("Unable to read error value on line: " << nLine << ". Perhaps you are using an invalid separator?");

    xVal = xStr.toDouble();
    yVal = yStr.toDouble();

    if(first) {
      firstXVal = xVal;
      first = false;
    }

    if(relative) {
      xVal -= firstXVal;
    }

    x_.push_back(xVal);
    y_.push_back(yVal);
    error_.push_back(eStr.toDouble());
  }

  fin.close();

  // Resize internal arrays used for managing model components

  compositeModel_.resize(x_.size());
}

/**.......................................................................
 * Clear any composite model that is currently set
 */
void DataSet1D::clearModel()
{
  for(unsigned iDat=0; iDat < x_.size(); iDat++)
    compositeModel_[iDat] = 0.0;

  if(displayInitialized_)
    for(unsigned iDat=0; iDat < xDisplay_.size(); iDat++)
      compositeModelDisplay_[iDat] = 0.0;
}

/**.......................................................................
 * Define what it means to add a model to a 1D data set
 */
void DataSet1D::addModel(gcp::util::Model& model)
{
  Generic1DModel& model1D = dynamic_cast<Generic1DModel&>(model);

  if(applies(model1D)) {

    for(unsigned iDat=0; iDat < x_.size(); iDat++) 
      compositeModel_[iDat] += model1D.eval(x_[iDat]);

#if 0

    static bool first=true;
    if(first) {
      initializeForDisplay();
      PgUtil::setInteractive(false);
      PgUtil::open("/xs");
      first = false;
    }

    PgUtil::linePlot(xDisplay_, compositeModel_, "", "", "", true);
#endif

    //------------------------------------------------------------
    // If we are displaying, compute the model on the display grid too
    //------------------------------------------------------------

    if(displayInitialized_) {
      for(unsigned iDat=0; iDat < xDisplay_.size(); iDat++) 
	compositeModelDisplay_[iDat] += model1D.eval(xDisplay_[iDat]);
    }
  }
}

/**.......................................................................
 * Compute chisq for a 1D model
 */
gcp::util::ChisqVariate DataSet1D::computeChisq()
{
  ChisqVariate chisq;

  double val;
  for(unsigned iDat=0; iDat < x_.size(); iDat++) {
    val = (y_[iDat] - compositeModel_[iDat]) / error_[iDat];
    chisq += val*val;
  }

  return chisq;
}

void DataSet1D::display()
{
  std::ostringstream osx, osy, ost;
  osx << name_ << ".x";
  osy << name_ << ".y";
  ost << "Dataset " << name_;
  PgUtil::errPlot(x_, y_, error_, osx.str(), osy.str(), ost.str());
}

void DataSet1D::displayCompositeModel()
{
  PgUtil::linePlot(xDisplay_, compositeModelDisplay_, "", "", "", true);
}

void DataSet1D::displayResiduals()
{
  std::vector<double> res(x_.size());

  for(unsigned i=0; i < x_.size(); i++)
    res[i] = y_[i] - compositeModel_[i];

  std::ostringstream osx, osy, ost;
  osx << name_ << ".x";
  osy << name_ << ".y";
  ost << "Dataset " << name_;

  PgUtil::errPlot(x_, res, error_, osx.str(), "Res", "");
}

std::vector<double> DataSet1D::getXData()
{
  return x_;
}

std::vector<double> DataSet1D::getResiduals()
{
  std::vector<double> res(x_.size());

  for(unsigned i=0; i < x_.size(); i++)
    res[i] = y_[i] - compositeModel_[i];

  return res;
}

/**.......................................................................
 * Define what it means to generate noise for a 1D data set
 */
void DataSet1D::simulateData(double sigma)
{
  vector<double> errs(x_.size()), errs2(x_.size());

  PgUtil::linePlot(compositeModel_);

  for(unsigned i=0; i < x_.size(); i++) {
    errs[i] = Sampler::generateGaussianSample(sigma);
    errs2[i] = errs[i]*errs[i]/(sigma*sigma);
    compositeModel_[i] += errs[i];
  }

  COUT("rms  = " << Stats::rms(errs));
  COUT("chi2 = " << Stats::sum(errs2));

  PgUtil::linePlot(compositeModel_);
}

/**.......................................................................
 * Define what it means to write a 1D composite model to a file
 */
void DataSet1D::writeCompositeModelToFile(std::string fileName, double sigma)
{
  std::ofstream fout;
  fout.open(fileName.c_str(), ios::out);

  if(!fout) {
    ThrowError("Unable to open file: " << fileName);
  }

  for(unsigned i=0; i < y_.size(); i++)
    fout << x_[i] << " " << compositeModel_[i] << " " << sigma << std::endl;

  fout.close();
}

void DataSet1D::initializeForDisplay()
{
  unsigned n = getUintVal("n");
  xDisplay_.resize(n);
  compositeModelDisplay_.resize(n);

  double xMin = Stats::min(x_);
  double xMax = Stats::max(x_);

  double dx = (xMax - xMin)/(n-1);
  for(unsigned i=0; i < n; i++) {
    xDisplay_[i] = xMin + dx*i;
    compositeModelDisplay_[i] = 0.0;
  }

  displayInitialized_ = true;
}
