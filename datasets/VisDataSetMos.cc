#include "gcp/datasets/VisDataSet.h"
#include "gcp/datasets/VisDataSetMos.h"

#include "gcp/models/PtSrcModel.h"

#include "gcp/util/DataType.h"
#include "gcp/util/Geoid.h"

#include "cpgplot.h"

using namespace gcp::datasets;
using namespace gcp::util;
using namespace std;

/**.......................................................................
 * Constructor.
 */
VisDataSetMos::VisDataSetMos() 
{
  //-----------------------------------------------------------------------
  // VisDataSetMos inherits the same interface as VisDataSet
  //-----------------------------------------------------------------------

  addParameter("perc",           DataType::DOUBLE, "If assigned, the data will be gridded to achieve this percent correlation");
  addParameter("npix",           DataType::UINT,   "The number of pixels to use when gridding the data");		  
  addParameter("size",           DataType::DOUBLE, "The size of the image to which the data will be gridded");
  
  addParameter("wtscale",        DataType::STRING, "A factor by which to scale the weights in the UVF file (or 'auto' to estimate from the data)");
  addParameter("uvmin",          DataType::DOUBLE, "The minimum UV radius to use in chi-squared computation");	
  addParameter("uvmax",          DataType::DOUBLE, "The maximum UV radius to use in chi-squared computation");
  addParameter("store",          DataType::BOOL,    "True to store data internally on readin");

  addParameter("useanttypes",    DataType::STRING, "List of antenna types to include on read-in");
  addParameter("uvtaper",        DataType::STRING, "UV-taper to apply to the weights");

  //------------------------------------------------------------
  // And has some of its own parameters
  //------------------------------------------------------------

  remParameter("file");
  addParameter("file",           DataType::STRING, "A list of data files to include in this dataset.  To add files to the mosaic, use 'file += filename'.  Note that parameters can be adjusted for individual datasets in the filelist by name. i.e., if 'mosname' is the name of your mosaicked dataset, then parameters of the first dataset can be accessed as 'mosname.dataset1.parname'");

  addParameter("power",          DataType::DOUBLE, "Point on the beam to which images should be displayed (defaults to 0.5)");
  addParameter("wtmin",          DataType::DOUBLE,  "Minimum valid weight for the combined map (used for mosaicking only)");

  //------------------------------------------------------------
  // Default to displaying to the half-power point
  //------------------------------------------------------------

  setParameter("power", "0.5");

  addParameter("forcewt",        DataType::BOOL,   "If true, force weights to match the data variance (wtscale must be set to 'auto')");
  addParameter("displaybeam",    DataType::BOOL,   "If true, display the synthesized beam on read-in");

  //------------------------------------------------------------
  // Add position information for this object too
  //------------------------------------------------------------

  addParameter("reversedelays",  DataType::BOOL,   "Flip the sense of the delays on read-in");
  addParameter("reversedisplay", DataType::BOOL,   "Flip the sense of the x-display");

  //------------------------------------------------------------
  // Default to displaying like difmap
  //------------------------------------------------------------

  setParameter("reversedisplay", "false");

  //------------------------------------------------------------
  // Add names for image output files
  //------------------------------------------------------------

  remParameter("dataimage");
  addParameter("dataimage",      DataType::STRING, "If specified, a FITS image of the data (SNR) will be written to this file");

  remParameter("modelimage");
  addParameter("modelimage",     DataType::STRING, "If specified, a FITS image of the model (SNR) will be written to this file");

  remParameter("resimage");
  addParameter("resimage",       DataType::STRING, "If specified, a FITS image of the residuals (SNR) will be written to this file");

  addParameter("noiseimage",     DataType::STRING, "If specified, a FITS image of the noise (Jy/beam) will be written to this file");

  addParameter("datauvf",        DataType::STRING,  "UVF file to output data visibilities");  
  addParameter("resuvf",         DataType::STRING,  "UVF file to output residual visibilities");  
  addParameter("modeluvf",       DataType::STRING,  "UVF file to output model visibilities");  

  addParameter("clean",          DataType::BOOL,    "If true, a clean model image will be displayed.  Note: if clean = true, you must force the data to a fixed resolution using parameters 'size' and 'npix' (i.e., you cannot use the 'perc' option, which is the default)");
  addParameter("cleantype",      DataType::STRING,  "If clean=true, the type of clean image to display, one of: 'model' or 'delta'.  If cleantype = model, then any currently defined models will be used to construct the clean image.  If cleantype = delta, then a delta-function model will be iteratively constructed using the (Hogbom) CLEAN algorithm.");
  addParameter("cleaniter",      DataType::UINT,    "If clean=true and cleantype=delta, the number of clean iterations to perform (default is 100)");
  addParameter("cleangain",      DataType::DOUBLE,  "If clean=true and cleantype=delta, the clean gain to use when subtracting model components (default is 0.05)");
  addParameter("cleancutoff",    DataType::DOUBLE,  "If clean=true and cleantype=delta, then cleaning will stop when the maximum absolute residual reaches this value (SNR).  Default is no cutoff (0.0)");
  addParameter("cleanwindow",    DataType::STRING,  "If clean=true and cleantype=delta, and no model is specified, then this sets a clean window to be searched to iteratively build up a model.  Format is: 'cleanwindow = [xmin:xmax, ymin:ymax, units]', or 'cleanwindow = [xpos +- xdelta, ypos +- ydelta, units]', where 'units' are the units in which the window is specified.  Alternately, 'cleanwindow = [tag +- delta, units]' can be used, where tag is one of: 'min', 'max' or 'abs' to set up a clean window about the minimum, maximum or absolute maximum in the image. Use 'cleanwindow += [xmin:xmax, ymin:ymax, units]' to specify more than one clean window.");

  //------------------------------------------------------------
  // Set defaults for some parameters
  //------------------------------------------------------------

  setParameter("cleanwindow", "[abs +- 0.025, abs +- 0.025, deg]");
  setParameter("cleaniter", "100");
  setParameter("cleangain", "0.05");
  setParameter("cleantype", "model");
  setParameter("cleancutoff", "0.0");
}

/**.......................................................................
 * Destructor.
 */
VisDataSetMos::~VisDataSetMos() {}

/**.......................................................................
 * Overload the base-class loadData method to define what happens when
 * loadData() is called on this object.  
 * 
 * The intention is that multiple vis datasets can be defined as part
 * of a VisDataSetMos, via the "fileList" parameter, which can be a
 * list of filenames.  
 *
 * This class instantiates the datasets as the fileList parameter is
 * specified, so at this point, all of our datasets should exist, but
 * need to be loaded().  Prior to loading, we will also force them to
 * inherit the values of relevant parameters
 */
void VisDataSetMos::loadData(bool simulate)
{
  if(!getParameter("file", false)->data_.hasValue())
    ThrowSimpleColorError("You must specify a file list via the file parameter", "red");

  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;

    //------------------------------------------------------------
    // Initialize this dataset's parameters to ours (but don't
    // override the position, which must be filled in from the data
    // load or by explicit specification from the parsing interface,
    // or the file name, which will have been specified when the
    // dataset was instantiated).
    // 
    // We also copy only parameters that haven't been otherwise
    // specified (false argument to copyParameters() below)
    //------------------------------------------------------------
    
    std::map<std::string, std::string> exc;
    
    exc["ra"]   = "ra";
    exc["dec"]  = "dec";
    exc["file"] = "file";
    
    dataSet->copyParameters(this, exc, false);
  }

  //------------------------------------------------------------
  // Now load the data
  //------------------------------------------------------------

  DataSetManager::loadData(simulate);

  if(dataSetMap_.size() == 1) {
    VisDataSet* dataSet = (VisDataSet*)dataSetMap_.begin()->second;
    dataSet->insertSynthesizedBeamModelForPlots(pgManager_);
  }
}

/**.......................................................................
 * Define what it means to display data for a mosaicked data set
 */
void VisDataSetMos::display()
{
  if(getParameter("cmap", false)->data_.hasValue()) {
    PgUtil::setColormap(getStringVal("cmap"));
  }

  display(VisDataSet::ACC_DATA);
}

/**.......................................................................
 * Define what it means to display the beam for a mosaicked data set
 */
void VisDataSetMos::displayBeam()
{
  PgUtil::setInteractive(interactive_);
  display(VisDataSet::ACC_BEAM);
}

/**.......................................................................
 * Define what it means to display the model for a mosaicked data set
 */
void VisDataSetMos::displayCompositeModel()
{  
  PgUtil::setInteractive(interactive_);

  bool doClean = false;

  if(getParameter("clean", false)->data_.hasValue())
    doClean = getBoolVal("clean");

  if(doClean) {
    if(getStringVal("cleantype") == "model")
      display(VisDataSet::ACC_CLEAN);
    else if(getStringVal("cleantype") == "delta")
      clean();
    else
      ThrowSimpleColorError("Unrecognized cleantype: " << getStringVal("cleantype"), "red");

  } else {
    display(VisDataSet::ACC_MODEL);
  }
}

/**.......................................................................
 * Define what it means to display residuals for a mosaicked data set
 */
void VisDataSetMos::displayResiduals()
{
  if(getParameter("zmin", false)->data_.hasValue() && getParameter("zmax", false)->data_.hasValue()) {
    PgUtil::setZmin(getDoubleVal("zmin"));
    PgUtil::setZmax(getDoubleVal("zmax"));
  }

  display(VisDataSet::ACC_RES);
}

/**.......................................................................
 * Get the requested mosaicked image
 */
void VisDataSetMos::getImage(VisDataSet::AccumulatorType type, Image& image, Image& noise)
{
  image   = getImage();
  Image wtimage = image;

  PgUtil::setInteractive(interactive_);

  image.setUnits(Unit::UNITS_SNR);

  image.zero();
  image.invalidate();
  wtimage.zero();
  wtimage.invalidate();

  //------------------------------------------------------------
  // Accumulate images from our datasets now
  //------------------------------------------------------------

  for(std::map<std::string, gcp::util::DataSet*>::iterator iter=dataSetMap_.begin(); iter != dataSetMap_.end(); iter++) {
    gcp::datasets::VisDataSet* visDataSet = (gcp::datasets::VisDataSet*)iter->second;
    visDataSet->getImage(image, wtimage, type);
  }

  //------------------------------------------------------------
  // Take the sqrt() of the variance to get the noise
  //------------------------------------------------------------

  Image wtsqrt = wtimage.getSqrt();

  if(getParameter("wtmin", false)->data_.hasValue()) {
    wtsqrt.invalidateDataLessThan(getDoubleVal("wtmin"));
    wtsqrt.setInvalidPixelsTo(0.0);
  }

  //------------------------------------------------------------
  // Divide the image by the noise to get SNR
  //------------------------------------------------------------


  image /= wtsqrt;
  image.setInvalidPixelsTo(0.0);

  if(type == VisDataSet::ACC_BEAM) { 
    image /= image.max();
  }

  noise = 1.0/wtsqrt;
}

/**.......................................................................
 * Define what it means to display data for a mosaicked data set
 */
void VisDataSetMos::display(VisDataSet::AccumulatorType type)
{
  Image image, noise;

  if(type == VisDataSet::ACC_CLEAN) {
    getImage(type, image, noise);
    Image res;
    getImage(VisDataSet::ACC_RES, res, noise);
    image += res;
  } else {
    getImage(type, image, noise);
  }

  //------------------------------------------------------------
  // Adjust the window to display the real data aspect ratio
  //------------------------------------------------------------

  PgUtil::setInteractive(interactive_);
  PgUtil::setWnad(true);

  //------------------------------------------------------------
  // If a colormap was specified, set it now
  //------------------------------------------------------------

  if(getParameter("cmap", false)->data_.hasValue()) {
    PgUtil::setColormap(getStringVal("cmap"));
  } else {
    PgUtil::setColormap("grey");
  }

  //------------------------------------------------------------
  // Default to displaying the full data range, unless a range was
  // specified.  If displaying the beam, ignore user-specified range
  //------------------------------------------------------------

  if(type == VisDataSet::ACC_DATA) {
    zmin_ = image.min();
    zmax_ = image.max();
  }

  if(type != VisDataSet::ACC_BEAM) {
    if(getParameter("zmin", false)->data_.hasValue() && getParameter("zmax", false)->data_.hasValue()) {
      zmin_ = getDoubleVal("zmin");
      zmax_ = getDoubleVal("zmax");
    }
  } else {
    zmin_ = 0.0;
    zmax_ = 0.0;
    image.units_    = Unit::UNITS_NONE;
    image.hasUnits_ = false;
  }

  VisDataSet::setWedge(type, zmin_, zmax_, getParameter("zmin", false)->data_.hasValue() && getParameter("zmax", false)->data_.hasValue());

  PgUtil::setZmin(zmin_);
  PgUtil::setZmax(zmax_);

  //------------------------------------------------------------
  // And display it
  //------------------------------------------------------------

  bool reverse = getBoolVal("reversedisplay");

  PgUtil::useHeader(true);
  PgUtil::setHeader(VisDataSet::typeString(type), PgUtil::JUST_LEFT);

  if(reverse)
    image.display();
  else
    image.difmapDisplay();

  PgUtil::useHeader(false);

  //------------------------------------------------------------
  // If requested, write out the image too
  //------------------------------------------------------------

  switch(type) {
  case VisDataSet::ACC_DATA:

    if(getParameter("dataimage", false)->data_.hasValue()) {
      image.writeToFitsFile(getStringVal("dataimage"));
    }

    break;
  case VisDataSet::ACC_MODEL:
  case VisDataSet::ACC_CLEAN:

    if(getParameter("modelimage", false)->data_.hasValue()) {
      image.writeToFitsFile(getStringVal("modelimage"));
    }

    break;
  case VisDataSet::ACC_RES:
    if(getParameter("resimage", false)->data_.hasValue()) {

      image.writeToFitsFile(getStringVal("resimage"));
    }

    if(getParameter("noiseimage", false)->data_.hasValue()) {
      noise.setUnits(Unit::UNITS_JYBEAM);
      noise.writeToFitsFile(getStringVal("noiseimage"));
    }

    break;
  default:
    break;
  }
}

void VisDataSetMos::writeData()
{
  for(std::map<std::string, gcp::util::DataSet*>::iterator iter=dataSetMap_.begin(); iter != dataSetMap_.end(); iter++) {
    gcp::datasets::VisDataSet* visDataSet = (gcp::datasets::VisDataSet*)iter->second;
    visDataSet->writeData();
  }
}

/**.......................................................................
 * Return an image with appropriate position and size to display the
 * mosaicked data
 */
void VisDataSetMos::getCenterPosition(HourAngle& raCenter, Declination& decCenter, Angle& xSize, Angle& ySize)
{
  double decDegMean  = 0.0;
  double raHoursMean = 0.0;

  bool first=true;

  //------------------------------------------------------------
  // Iterate over data sets
  //------------------------------------------------------------

  double decMinDeg, decMaxDeg;
  double raMinDeg, raMaxDeg;

  Angle size;
  size.setName(name_, "size");
  size.setVal(getDoubleVal("size"), getParameter("size", true)->units_);

  for(std::map<std::string, gcp::util::DataSet*>::iterator iter=dataSetMap_.begin(); 
      iter != dataSetMap_.end(); iter++) {

    gcp::datasets::VisDataSet* visDataSet = (gcp::datasets::VisDataSet*)iter->second;

    //------------------------------------------------------------
    // Store the largest primary beam size found in this dataset
    //------------------------------------------------------------

    Angle fwhm = visDataSet->estimateLargestPrimaryBeamFwhm();

    Declination decCurr = visDataSet->dec_;
    HourAngle   raCurr  = visDataSet->ra_;
    
    // Now compute the (rough) RA and Dec that correspond to the
    // extrema of the primary beam half-width

    double decCurrMinDeg = decCurr.degrees() - size.degrees()/2;
    double decCurrMaxDeg = decCurr.degrees() + size.degrees()/2;

    double raCurrMinDeg = raCurr.degrees() - (size.degrees()/2)/cos(decCurr.radians());
    double raCurrMaxDeg = raCurr.degrees() + (size.degrees()/2)/cos(decCurr.radians());
    
    // Now store the min/max RA DEC over all

    if(first) {

      decMinDeg = decCurrMinDeg;
      decMaxDeg = decCurrMaxDeg;
      raMinDeg  = raCurrMinDeg;
      raMaxDeg  = raCurrMaxDeg;
      first = false;

    } else {

      decMinDeg = decMinDeg < decCurrMinDeg ? decMinDeg : decCurrMinDeg;
      decMaxDeg = decMaxDeg > decCurrMaxDeg ? decMaxDeg : decCurrMaxDeg;

      raMinDeg  = raMinDeg < raCurrMinDeg ? raMinDeg : raCurrMinDeg;
      raMaxDeg  = raMaxDeg > raCurrMaxDeg ? raMaxDeg : raCurrMaxDeg;

    }
  }

  //------------------------------------------------------------
  // Now we have the extrema of all mosaicked data sets.  Compute the
  // center and size of an image large enough to encompass them all
  //------------------------------------------------------------

  double decCenterDeg = (decMinDeg + decMaxDeg)/2;
  double raCenterDeg  = (raMinDeg  + raMaxDeg)/2;

  decCenter.setDegrees(decCenterDeg);
  raCenter.setDegrees(raCenterDeg);

  //------------------------------------------------------------
  // Finally, get the total size of the image in each dimension
  //------------------------------------------------------------

  xSize.setDegrees((raMaxDeg - raMinDeg) * cos(decCenter.radians()));
  ySize.setDegrees(decMaxDeg - decMinDeg);

  //------------------------------------------------------------
  // Reset the center position to be the mean of the pointings
  //------------------------------------------------------------

  getMeanPosition(raCenter, decCenter);
}

/**.......................................................................
 * Compute the mean position of this mosaic
 */
void VisDataSetMos::getMeanPosition(HourAngle& raMean, Declination& decMean)
{
  //------------------------------------------------------------
  // Iterate over data sets, computing the mean RA, DEC position
  //------------------------------------------------------------

  double decMeanDeg = 0.0;
  double raMeanDeg  = 0.0;

  Lla lla;
  lla.setDatum(DATUM_SPHERE);
  lla.setCoordSystem(COORD_GEOCENTRIC);

  Geoid geoid(DATUM_SPHERE);
  lla.altitude_ = geoid.majorAxis();
  LengthTriplet xyzMean, xyz;
  xyzMean.setCoordSystem(COORD_ECF);

  double X,Y,Z;
  double XMean, YMean, ZMean;

  unsigned iDataSet=0;
  for(std::map<std::string, gcp::util::DataSet*>::iterator iter=dataSetMap_.begin(); 
      iter != dataSetMap_.end(); iter++) {

    gcp::datasets::VisDataSet* visDataSet = (gcp::datasets::VisDataSet*)iter->second;

    double raCurrDeg  =  visDataSet->ra_.degrees();
    double decCurrDeg = visDataSet->dec_.degrees();

    lla.longitude_.setDegrees(raCurrDeg);
    lla.latitude_.setDegrees(decCurrDeg);

    XMean = xyzMean.X_.meters();
    YMean = xyzMean.Y_.meters();
    ZMean = xyzMean.Z_.meters();

    //    COUTCOLOR("LLA = " << std::endl << lla, "green");
    xyz = geoid.geocentricLlaAndEnhToEcf(lla);
    //    COUTCOLOR("XYZ = " << std::endl << xyz, "magenta");

    X = xyz.X_.meters();
    Y = xyz.Y_.meters();
    Z = xyz.Z_.meters();

    XMean += (X - XMean) / (iDataSet + 1);
    YMean += (Y - YMean) / (iDataSet + 1);
    ZMean += (Z - ZMean) / (iDataSet + 1);

    xyzMean.X_.setMeters(XMean);
    xyzMean.Y_.setMeters(YMean);
    xyzMean.Z_.setMeters(ZMean);

    ++iDataSet;
  }

  //  raMeanDeg  += ( raCurrDeg -  raMeanDeg) / (iDataSet + 1);
  //  decMeanDeg += (decCurrDeg - decMeanDeg) / (iDataSet + 1);

  //  COUTCOLOR("XYZ Mean = " << std::endl << xyzMean, "magenta");
  lla = geoid.ecfToGeocentricLla(xyzMean);
  //  COUTCOLOR("LLA Mean = " << std::endl << lla, "green");

  raMean.setDegrees(lla.longitude_.degrees());
  decMean.setDegrees(lla.latitude_.degrees());

  //  COUTCOLOR("RA  Mean = " << raMean, "magenta");
  //  COUTCOLOR("DEC Mean = " << decMean, "magenta");
}

/**.......................................................................
 * Return an image with appropriate position and size to display the
 * mosaicked data
 */
Image VisDataSetMos::getImage()
{
  Angle size;
  size.setName(name_, "size");
  size.setVal(getDoubleVal("size"), getParameter("size", true)->units_);

  HourAngle raCenter;
  Declination decCenter;
  Angle xSize, ySize;

  getCenterPosition(raCenter, decCenter, xSize, ySize);

  //  COUT("RA center = " << raCenter);
  //  COUT("DEC center = " << decCenter);

  //------------------------------------------------------------
  // If a position was explicitly specified, override the position we
  // just calculated
  //------------------------------------------------------------

  if(getParameter("dec", false)->data_.hasValue()) {
    decCenter = dec_;
  }

  if(getParameter("ra", false)->data_.hasValue()) {
    raCenter = ra_;
  }

  unsigned npix;
  npix = getUintVal("npix");

  //------------------------------------------------------------
  // And get the nearest number of pixels at the current resolution 
  //------------------------------------------------------------

  double res = size.degrees()/npix;

  unsigned xNpix = (unsigned)ceil(xSize.degrees() / res);
  unsigned yNpix = (unsigned)ceil(ySize.degrees() / res);

  if(xNpix % 2 != 0) 
    xNpix += 1;

  if(yNpix % 2 != 0) 
    yNpix += 1;

  xSize.setDegrees(xNpix * res);
  ySize.setDegrees(yNpix * res);

  //------------------------------------------------------------
  // Now construct the image
  //------------------------------------------------------------

  Image image;
  image.xAxis().setNpix(xNpix);
  image.xAxis().setAngularSize(xSize);

  image.yAxis().setNpix(yNpix);
  image.yAxis().setAngularSize(ySize);

  image.setRaDecFft(raCenter, decCenter);
  image.zero();
  image.hasData_ = true;

  return image;
}

void VisDataSetMos::displayIfRequested()
{
  double vpmin = 0.1;
  double vpmax = 0.9;
  double vsep  = 0.15;
  double dvp   = (vpmax-vpmin-vsep)/2;
  double vp11  = vpmin;
  double vp12  = vpmin+dvp;
  double vp21  = vpmin+dvp+vsep;
  double vp22  = vpmax;

  PgUtil::initialize();

  if(getParameter("display", false)->data_.hasValue()) {
    if(getBoolVal("display")) {
      
      if(getParameter("dev", false)->data_.hasValue())
	PgUtil::open(getStringVal("dev"));
      else
	PgUtil::open("/xs");

      PgUtil::setVp(false);
      PgUtil::setOverplot(true);

      cpgsvp(vp11, vp12, vp21, vp22);
      display();

      cpgsvp(vp21, vp22, vp21, vp22);
      displayBeam();


      PgUtil::close();

      // Reset any legacy display settings

      PgUtil::setZmin(0.0);
      PgUtil::setZmax(0.0);
    }
  }
}

gcp::util::DataSet* VisDataSetMos::getDataSet(std::string name)
{
  for(std::map<std::string, gcp::util::DataSet*>::iterator iter=dataSetMap_.begin(); iter != dataSetMap_.end(); iter++) {
    gcp::util::DataSet* dataSet = iter->second;

    if(dataSet->name_ == name)
      return dataSet;
  }

  ThrowSimpleColorError("Invalid dataset name: " << name << " for dataset manager " << name_, "red");
  return 0;
}

/**.......................................................................
 * Initialize new datasets that we will manage from a list of files
 */
void VisDataSetMos::initializeDataSets(std::string fileList)
{
  //------------------------------------------------------------
  // Parse the list of filenames
  //------------------------------------------------------------

  std::vector<std::string> fileVec = getFileList(fileList);

  //------------------------------------------------------------
  // Now iterate over the list, instantiating datasets
  //------------------------------------------------------------

  for(unsigned i=0; i < fileVec.size(); i++) {

    String fileStr(fileVec[i]);

    //------------------------------------------------------------
    // Now iterate over the current file list, adding datasets
    //------------------------------------------------------------
    
    std::ostringstream name;
    unsigned iDataSet = dataSetMap_.size()+1;
    name << "dataset" << iDataSet;
    
    gcp::util::DataSet* dataSet=0;
    
    if(fileStr.contains("uvf")) {
      dataSet = DataSetManager::addDataSet("uvf", name.str());
    } else {
      ThrowColorError("Unable to determine the type of this dataset: " << fileStr, "red");
    }

    // Add a 'parameter' corresponding to this dataset, so it appears
    // in the parameter listing

    ostringstream os;
    os << "The object corresponding to dataset " << iDataSet;
    addParameter(name.str(),           DataType::STRING, os.str());

    //------------------------------------------------------------
    // Set the file name for this dataset -- will be needed when we
    // actually load the data
    //------------------------------------------------------------
    
    dataSet->setParameter("file", fileStr.str());
  }
}

/**.......................................................................
 * Parse a list of files into individual file names
 */
std::vector<std::string> VisDataSetMos::getFileList(std::string fileList)
{
  std::vector<std::string> fileVec;

  //------------------------------------------------------------
  // Iterate through the file parameter, looking for discrete file
  // names that were specified
  //------------------------------------------------------------

  String fileListStr = fileList;
  String nextFile;

  fileListStr.advanceToNextNonWhitespaceChar();

  do {

    nextFile = fileListStr.findNextInstanceOf("", false, " ", true, true);
    fileListStr.advanceToNextNonWhitespaceChar();

    if(!nextFile.isEmpty()) {
      fileVec.push_back(nextFile.str());
    } else {
      String remain = fileListStr.remainder();
      if(!remain.isEmpty()) {
	fileVec.push_back(remain.str());
      }
    }

  } while(!nextFile.isEmpty());

  return fileVec;
}

//=======================================================================
// Overloaded methods of ParameterManager to intercept parameter
// declarations that require internal action before loading data
//=======================================================================

/**.......................................................................
 * Set a parameter
 */
void VisDataSetMos::setParameter(std::string name, std::string val, std::string units)
{
  String nameStr(name);

  if(nameStr.contains("dataset")) {
    String dataSetName = nameStr.findNextInstanceOf("", false, ".", true, true);
    String parName     = nameStr.remainder();

    DataSet* dataSet = getDataSet(dataSetName.str());
    dataSet->setParameter(parName.str(), val, units);

  } else if(nameStr.contains("file")) {
    initializeDataSets(val);
    gcp::util::ParameterManager::setParameter(name, val, units);
  } else {
    gcp::util::ParameterManager::setParameter(name, val, units);
  }
}

/**.......................................................................
 * Increment a parameter
 */
void VisDataSetMos::incrementParameter(std::string name, std::string val, std::string units)
{
  String nameStr(name);

  if(nameStr.contains("dataset")) {
    String dataSetName = nameStr.findNextInstanceOf("", false, ".", true, true);
    String parName     = nameStr.remainder();

    DataSet* dataSet = getDataSet(dataSetName.str());
    dataSet->incrementParameter(parName.str(), val, units);

    //------------------------------------------------------------
    // If a file list is being incremented, add another dataset to
    // correspond to it
    //------------------------------------------------------------

  } else if(nameStr.contains("file")) {
    initializeDataSets(val);
    gcp::util::ParameterManager::incrementParameter(name, val, units);

    //------------------------------------------------------------ 
    // Else
    // just call down to the base-class method to increment the
    // parameter
    //------------------------------------------------------------

  } else {
    gcp::util::ParameterManager::incrementParameter(name, val, units);
  }
}

/**.......................................................................
 * Perform an iterative image-plane CLEAN of the image
 */
void VisDataSetMos::clean()
{
  //------------------------------------------------------------
  // First, clear the model to get rid of any legacy models that might
  // have been installed before us
  //------------------------------------------------------------

  clearModel();

  //------------------------------------------------------------
  // Iterate, removing delta functions
  //------------------------------------------------------------

  double val, flux, fluxTotal=0.0;
  Angle xOff, yOff;
  unsigned iMax;
  Image image, noise;
  std::vector<Model*> deltaComponents;

  getImage(VisDataSet::ACC_RES, image, noise);

  //------------------------------------------------------------
  // Parse window specifications now.  Pass in the image so that if
  // the user specified 'peak' or '-peak', we can replace with the
  // offsets corresponding to those
  //------------------------------------------------------------

  std::vector<gcp::util::Image::Window> windows = VisDataSet::parseWindows(getStringVal("cleanwindow"), image);

  double cleanGain   = getDoubleVal("cleangain");
  double cleanCutoff = getDoubleVal("cleancutoff");
  unsigned cleanIter = getUintVal("cleaniter");

  if(cleanGain > 1.0)
    ThrowSimpleColorError("CLEAN gain should be <= 1.0", "red");

#if 1
  COUT("Cleaning with niter = " << cleanIter << " gain = " << cleanGain << " and windows: ");
  for(unsigned iWin=0; iWin < windows.size(); iWin++) 
    COUT(windows[iWin] << std::endl);
#endif

  //------------------------------------------------------------
  // Iterate, removing delta functions
  //------------------------------------------------------------

  unsigned ndigit = ceil(log10((double)cleanIter));

  for(unsigned i=0; i < cleanIter; i++) {

    COUTCOLORNNL("\rClean iteration: " << std::setw(ndigit) << std::right << (i+1) << "/" 
    		 << cleanIter << " (" << std::right << setprecision(0) << std::fixed << (100*(double)(i+1)/cleanIter) << "%)", "green");

    //------------------------------------------------------------
    // Get the image and display prior to removal of the delta component
    //------------------------------------------------------------

    getImage(VisDataSet::ACC_RES, image, noise);
    image.getMax(val, xOff, yOff, iMax, windows, Image::TYPE_ABS);

    if(fabs(image.data_[iMax]) < cleanCutoff)
      break;

    flux = cleanGain*image.data_[iMax]*noise.data_[iMax];
    Model* model = VisDataSet::addDeltaFunctionModel(deltaComponents, flux, xOff, yOff);
    fluxTotal += flux;

    //------------------------------------------------------------
    // Now remove the model from the dataset
    //------------------------------------------------------------
    
    addModel(*model);
  }

  COUT("");

  //------------------------------------------------------------
  // Done with clean -- display restored image
  //------------------------------------------------------------

  COUT("Added models with combined flux of: " << fluxTotal);
  display(VisDataSet::ACC_CLEAN);

  for(unsigned iModel=0; iModel < deltaComponents.size(); iModel++) {
    delete deltaComponents[iModel];
    deltaComponents[iModel] = 0;
  }
}

void VisDataSetMos::addModel(gcp::util::Model& model)
{
  //------------------------------------------------------------
  // Check if this model has an absolute position.  If not, set it up
  // to be the same as our center position
  //------------------------------------------------------------

  for(std::map<std::string, gcp::util::DataSet*>::iterator diter = dataSetMap_.begin();
      diter != dataSetMap_.end(); diter++) {
    DataSet* dataSet = diter->second;

    if(dataSet->applies(model)) {

      Generic2DAngularModel* model2d = (Generic2DAngularModel*) &model;

      //------------------------------------------------------------
      // If this model has no absolute position, set its position to
      // be our mean position
      //------------------------------------------------------------

      model2d->checkPosition();

      if(!model2d->hasAbsolutePosition_) {
	model2d->setRa(ra_);
	model2d->setDec(dec_);
      }

      //------------------------------------------------------------
      // Now add the model with modified position
      //------------------------------------------------------------

      dataSet->addModel(model);
    }
  }
}

void VisDataSetMos::checkPosition(bool override)
{
  //------------------------------------------------------------
  // Call the underlying DSM method first so that all our managed
  // datasets are set up properly
  //------------------------------------------------------------

  DataSetManager::checkPosition(override);

  //------------------------------------------------------------
  // Calculate a mean position from the datasets we are managing
  //------------------------------------------------------------

  HourAngle raMean;
  Declination decMean;

  getMeanPosition(raMean, decMean);

  //------------------------------------------------------------
  // But don't override the position if it was explicitly set
  //------------------------------------------------------------
  
  if(getParameter("dec", false)->data_.hasValue())
    decMean = dec_;

  if(getParameter("ra", false)->data_.hasValue())
    raMean = ra_;

  //------------------------------------------------------------
  // Finally, set the result as our default position.  For mosaics,
  // the position of our maps will be recalculated to the center of a
  // map that is large enough to encompass the extrema of our
  // datasets, but this mean position will be used as the absolute
  // position for any model component that has not explicitly
  // specified its position, since relative model positions are not
  // defined for a mosaic
  //------------------------------------------------------------

  setRa(raMean);
  setDec(decMean);
}
