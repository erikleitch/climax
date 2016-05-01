#include "gcp/datasets/Array2DDataSet.h"

#include "gcp/fftutil/Generic2DAngularModel.h"

using namespace std;

using namespace gcp::datasets;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Array2DDataSet::Array2DDataSet() 
{
  hasErrorVal_         = false;
  hasData_             = false;

  addParameter("frequency", DataType::DOUBLE, "When simulating, the frequency at which to calculate the model");
  addParameter("error",     DataType::DOUBLE, "Noise sigma to use");
  addParameter("n",         DataType::UINT,   "The number of samples to use when simulating data for this dataset");
  addParameter("size",      DataType::DOUBLE, "The size of the image to use when simulating data for this dataset"); 
  addParameter("units",     DataType::STRING, "The units of the x and y data axes");
  addParameter("thetaExc",  DataType::DOUBLE, "If specified, data within thetaExc of the current model center will be excluded from the fit");
}

/**.......................................................................
 * Destructor.
 */
Array2DDataSet::~Array2DDataSet() {}

void Array2DDataSet::initializeErrors()
{
  if(getParameter("error", false)->data_.hasValue_) {
    hasErrorVal_ = true;
    errorVal_    = getDoubleVal("error");  
  }
}

void Array2DDataSet::loadData(bool simulate)
{
  if(simulate) {
    initializeForSim();
  } else {
    loadDataFromFile();
  }
}

void Array2DDataSet::initializeForSim()
{
}

/**.......................................................................
 * Load data from a specified file
 */
void Array2DDataSet::loadDataFromFile()
{
  //------------------------------------------------------------
  // Initialize errors for this data set
  //------------------------------------------------------------

  initializeErrors();

  //------------------------------------------------------------
  // Initialize our internal image from a file
  //------------------------------------------------------------

  initializeFromTextFile(getStringVal("file"));

  initializeImage();

  if(getParameter("frequency", false)->data_.hasValue_)
    frequency_.setVal(getDoubleVal("frequency"), getParameter("frequency", true)->units_);
}

void Array2DDataSet::initializeImage()
{
  double xmin = x_[0];
  double xmax = x_[0];
  double ymin = y_[0];
  double ymax = y_[0];

  for(unsigned i=0; i < x_.size(); i++) {
    xmin = (xmin < x_[i] ? xmin : x_[i]);
    xmax = (xmax > x_[i] ? xmax : x_[i]);
    ymin = (ymin < y_[i] ? ymin : y_[i]);
    ymax = (ymax > y_[i] ? ymax : y_[i]);
  }

  // Approximate resolution of the data:

  unsigned npix = (unsigned)sqrt((double)x_.size());

  // Store the axis units

  axisUnits_.setUnits(getStringVal("units"));
  axisUnits_.setVal(1.0, getStringVal("units"));

  // And allocate an image we can use to interpolate the data onto for
  // display, etc.

  Angle xw, yw;
  xw.setVal((xmax - xmin), getStringVal("units"));
  yw.setVal((ymax - ymin), getStringVal("units"));

  image_.xAxis().setAngularSize(xw);
  image_.xAxis().setNpix(npix);

  image_.yAxis().setAngularSize(yw);
  image_.yAxis().setNpix(npix);

  image_.raRefPix_  = npix/2;
  image_.decRefPix_ = npix/2;
}

/**.......................................................................
 * Initialize this data set from a fits file
 */
void Array2DDataSet::initializeFromTextFile(std::string fileName)
{
  unsigned nLine = countLines(fileName);

  static String line;

  std::ifstream fin;
  fin.open(fileName.c_str(), ios::in);

  if(!fin) {
    ThrowColorError(std::endl << "Unable to open file: " << fileName, "red");
  }

  String xStr, yStr, dStr, eStr;
  double firstXVal, xVal, yVal;
  bool first = true;

  x_.resize(nLine);
  y_.resize(nLine);
  data_.resize(nLine);
  error_.resize(nLine);

  nLine = 0;
  unsigned iLine=0;
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

    dStr = line.findNextStringSeparatedByChars(" \t,", true);

    if(dStr.isEmpty()) 
      ThrowError("Unable to read data value on line: " << nLine << ". Perhaps you are using an invalid separator?");

    eStr = line.findNextStringSeparatedByChars(" \t,", true);

    if(eStr.isEmpty() && !hasErrorVal_) {
      ThrowError("Unable to read error value on line: " << nLine << ". Perhaps you are using an invalid separator?" 
		 << "(You must either specify a fixed error value with the parameter 'error' or the file must contain a per-point error");
    }

    x_[iLine] = xStr.toDouble();
    y_[iLine] = yStr.toDouble();
    data_[iLine] = dStr.toDouble();

    if(!eStr.isEmpty())
      error_[iLine] = eStr.toDouble();

    iLine += 1;
  }

  fin.close();

  // Resize internal arrays used for managing model components

  compositeModel_.resize(x_.size());
  modelComponent_.resize(x_.size());
}

/**.......................................................................
 * Initialize this data set from a fits file
 */
unsigned Array2DDataSet::countLines(std::string fileName)
{
  static String line;

  std::ifstream fin;
  fin.open(fileName.c_str(), ios::in);

  if(!fin) {
    ThrowColorError(std::endl << "Unable to open file: " << fileName, "red");
  }

  unsigned nLine=0;

  while(!fin.eof()) {

    line.initialize();
    getline(fin, line.str());

    line.advanceToNextNonWhitespaceChar();

    // Ignore empty lines

    if(line.isEmpty()) 
      continue;      

    // Ignore empty lines

    if(line[0] == '/')
      continue;

    ++nLine;
  }

  fin.close();

  return nLine;}

void Array2DDataSet::clearModel()
{
  compositeModel_ = 0.0;
  hasData_ = false;
}

/**.......................................................................
 * Generic method to compute chisq for a 2D data set (defaults to
 * image data set -- inheritors can override)
 *
 * We allow either a fixed error (hasErrorVal_ == true) to be
 * specified, or a per-pixel image (hasErrorVal_ == false).
 * 
 * Additionally, if an exclusion radius was specified, we exclude data
 * within this radius of the current model center.
 */
gcp::util::ChisqVariate Array2DDataSet::computeChisq()
{
  ChisqVariate chisq;
  double val;

  unsigned n = x_.size();
  double errorVal;
  
  for(unsigned i=0; i < n; i++) {
    
    errorVal = hasErrorVal_ ? errorVal_ : error_[i];
    
    val = (data_[i] - compositeModel_[i]) / errorVal;

    //    COUT("d[" << i << "] = " << data_[i] << " m[" << i << "] = " << compositeModel_[i] 
    //	 << " delta = " << data_[i] - compositeModel_[i]);

    chisq += val*val;
  }
  
  return chisq;
}

void Array2DDataSet::display()
{
  PgUtil::setWnad(true);

  if(getParameter("cmap", false)->data_.hasValue()) {
    PgUtil::setColormap(getStringVal("cmap"));
  } else {
    PgUtil::setColormap("grey");
  }

  image_.fillFrom(axisUnits_, x_, y_, data_);
  image_.display();

  if(getParameter("dataimage", false)->data_.hasValue()) {
    writeToTextFile(data_, getStringVal("dataimage"));
  }
}

void Array2DDataSet::displayCompositeModel()
{
  PgUtil::setInteractive(true);


  image_.fillFrom(axisUnits_, x_, y_, compositeModel_);
  image_.display();

  if(getParameter("modelimage", false)->data_.hasValue()) {
    writeToTextFile(compositeModel_, getStringVal("modelimage"));
  }
}

void Array2DDataSet::displayResiduals()
{
  std::valarray<double> res = data_ - compositeModel_;
  image_.fillFrom(axisUnits_, x_, y_, res);

  //  PgUtil::setWnad(false);
  //  PgUtil::errPlot(x_, res, error_);

  PgUtil::setInteractive(true);
  PgUtil::setZmin(res.min());
  PgUtil::setZmax(res.max());

  image_.display();

  PgUtil::setZmin(0.0);
  PgUtil::setZmax(0.0);

  if(getParameter("resimage", false)->data_.hasValue()) {
    writeToTextFile(res, getStringVal("resimage"));
  }
}

/**.......................................................................
 * Add a model to this dataset
 */
void Array2DDataSet::addModel(gcp::util::Model& model)
{
  //------------------------------------------------------------
  // Only proceed if this model applies to this data set
  //------------------------------------------------------------

  if(applies(model)) {

    Generic2DAngularModel& model2D = dynamic_cast<Generic2DAngularModel&>(model);
    
    //------------------------------------------------------------
    // Now proceed with the calculation
    //------------------------------------------------------------
    
    model2D.fillArray(dataSetType_, axisUnits_, x_, y_, modelComponent_, &frequency_);
    
    if(!hasData_) {
      if(model2D.modelType_ & ModelType::MODEL_ADDITIVE) {
	compositeModel_ = modelComponent_;
	hasData_ = true;
      } else {
	ThrowSimpleColorError("Can't assign a multiplicative model without any additive models defined", "red");
      }
    } else {
      if(model2D.modelType_ & ModelType::MODEL_ADDITIVE)
	compositeModel_ += modelComponent_;
      else {
	compositeModel_ *= modelComponent_;
      }
    }

#if 0
    displayCompositeModel();
#endif
  }
}

/**.......................................................................
 * Base-class method to initialize data that depend on absolute position
 */
void Array2DDataSet::initializePositionDependentData()
{
}

void Array2DDataSet::simulateData(double sigma)
{
  vector<double> errs(x_.size()), errs2(x_.size());

  displayCompositeModel();

  for(unsigned i=0; i < x_.size(); i++) {
    errs[i] = Sampler::generateGaussianSample(sigma);
    errs2[i] = errs[i]*errs[i]/(sigma*sigma);
    compositeModel_[i] += errs[i];
  }

  displayCompositeModel();
}

void Array2DDataSet::writeToTextFile(std::valarray<double>& d, std::string fileName)
{
  std::ofstream fout;
  fout.open(fileName.c_str(), ios::out);

  if(!fout) {
    ThrowError("Unable to open file: " << fileName);
  }

  for(unsigned i=0; i < x_.size(); i++)
    fout << x_[i] << " " << y_[i] << " " << d[i] << std::endl;

  fout.close();
}
