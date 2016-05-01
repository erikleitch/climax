#include "fftutil/ImageManager.h"

#include "util/FitsReader.h"
#include "util/FitsBinTableReader.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ImageManager::ImageManager() 
{
  addParameter("scale",       DataType::STRING, "Scale factor to apply to the image: image = (image + offset) * scale");
  addParameter("offset",      DataType::STRING, "Offset to apply to the image: image = (image + offset) * scale");
  addParameter("npix",        DataType::STRING, "The number of pixels to use when simulating data for this dataset");			   
  addParameter("size",        DataType::STRING, "The size of the image to use when simulating data for this dataset");		   
  addParameter("trans",       DataType::STRING, "A transformation to apply to the image");		   
  addParameter("sigmasmooth", DataType::DOUBLE, "Gaussian sigma of a smoothing kernel to apply to the image");		   
  addParameter("sigmaapod",   DataType::DOUBLE, "Gaussian sigma of a smoothing kernel to apply to the image");		   
  addParameter("units",       DataType::STRING, "Units for the image data");		   
  addParameter("region",      DataType::STRING, 
	       "Region of the input image to use.  Format is 'ixmin:ixmax, iymin:iymax;' or 'nx, ny'");

  addParameter("extname",     DataType::STRING, "If no image is contained in the primary header unit, this specifies FITS extension from which to read the image (defaults to 'EVENTS')");
  addParameter("xcol",        DataType::STRING, "If initialized from FITS binary table, the x-column to use (defaults to 'x')");
  addParameter("ycol",        DataType::STRING, "If initialized from FITS binary table, the x-column to use (defaults to 'y')");
  addParameter("dcol",        DataType::STRING, "If initialized from FITS binary table, the data column to use (defaults to 'energy')");

  setParameter("scale",  "1.0");
  setParameter("offset", "0.0");

  setParameter("extname", "EVENTS");
  setParameter("xcol",    "x");
  setParameter("ycol",    "y");
  setParameter("dcol",    "energy");
}

/**.......................................................................
 * Destructor.
 */
ImageManager::~ImageManager() {}

Image ImageManager::initializeImage(std::string fileName, Image& image, ObsInfo* obs)
{
  bool rebin = true;

  if(isImage(fileName)) {
    image.initializeFromFitsFile(fileName, obs);
  } else {
    image.initializeFromFitsTable(fileName, getStringVal("extname"), 
				  getStringVal("xcol"), getStringVal("ycol"), getStringVal("dcol"), obs);
  }

  //------------------------------------------------------------
  // Whatever the reference pixel in the FITS file, set the center of
  // the returned image to be the new reference pixel
  //------------------------------------------------------------

  image.setRaDec(image.getRa(), image.getDec());

  return initializeImageParameters(image, rebin);
}

Image ImageManager::initializeImage(std::string fileName, ObsInfo* obs)
{
  Image image;
  bool rebin = true;

  if(isImage(fileName)) {
    image.initializeFromFitsFile(fileName, obs);
  } else {
    image.initializeFromFitsTable(fileName, getStringVal("extname"), 
				  getStringVal("xcol"), getStringVal("ycol"), getStringVal("dcol"), obs);
  }

  //------------------------------------------------------------
  // Whatever the reference pixel in the FITS file, set the center of
  // the returned image to be the new reference pixel
  //------------------------------------------------------------

  image.setRaDec(image.getRa(), image.getDec());

  return initializeImageParameters(image, rebin);
}

/**.......................................................................
 * Initialize the image, based on parameter values
 */
Image ImageManager::initializeImageParameters(Image& image, bool rebinIfRequested)
{
  Image imageCut, imageRet;

  //------------------------------------------------------------
  // Apply scaling/offset, if any
  //------------------------------------------------------------

  double offset = getVal(image, "offset");
  double scale  = getVal(image, "scale");

  image += offset;
  image *= scale;

  COUTCOLOR(std::endl << "Native image size is:  " << image.xAxis().getNpix() << " x " << image.yAxis().getNpix() << " pixels ("
	    << std::setprecision(2) << image.xAxis().getAngularSize().degrees() << " x " << image.yAxis().getAngularSize().degrees() << " degrees)", "magenta");
  imageCut = image;

  //------------------------------------------------------------
  // See if a region was specified
  //------------------------------------------------------------

  if(getParameter("region", false)->data_.hasValue_) {
    String region = getStringVal("region");
    imageCut = extractRegion(region, image);
    COUTCOLOR("Extracting image size: " << imageCut.xAxis().getNpix() << " x " << imageCut.yAxis().getNpix() << " pixels ("
	      << std::setprecision(2) << imageCut.xAxis().getAngularSize().degrees() << " x " << imageCut.yAxis().getAngularSize().degrees() << " degrees)", "magenta");
  }

  //------------------------------------------------------------
  // If a size was specified, set (or re-set) the size of the image
  //------------------------------------------------------------
  
  if(getParameter("size", false)->data_.hasValue_) {

    std::ostringstream os;
    os << getStringVal("size") << getParameter("size", true)->units_;
    String val(os.str());
    val.strip(' ');

    Angle xSize, ySize;

    if(val.contains("x")) {
      String xsize = val.findNextNumericString();
      val.advanceToNextNumericChar();
      String ysize = val.findNextNumericString();
      String units = val.remainder();

      xSize.setVal(xsize.toDouble(), units.str());
      ySize.setVal(ysize.toDouble(), units.str());

    } else {
      String size  = val.findNextNumericString();
      String units = val.remainder();

      xSize.setVal(size.toDouble(), units.str());
      ySize.setVal(size.toDouble(), units.str());

    }

    imageCut.xAxis().setAngularSize(xSize);
    imageCut.yAxis().setAngularSize(ySize);
  }

  //------------------------------------------------------------
  // If npix was specified, rebin the image now
  //------------------------------------------------------------
  
  imageRet = imageCut;
  
  if(rebinIfRequested && getParameter("npix", false)->data_.hasValue_) {
    String val(getStringVal("npix"));
    
    unsigned xNpix, yNpix;
    if(val.contains("x")) {
      String xnpix = val.findNextNumericString();
      val.advanceToNextNumericChar();
      String ynpix = val.findNextNumericString();

      xNpix = xnpix.toInt();
      yNpix = ynpix.toInt();

    } else {
      xNpix = val.toInt();
      yNpix = xNpix;
    }

    imageRet.xAxis().setNpix(xNpix);
    imageRet.yAxis().setNpix(yNpix);

    
    if(imageCut.hasAbsolutePosition_) {
      imageRet.setRaDec(imageCut.getRa(), imageCut.getDec());
    }

    COUTCOLOR("Binning image to " << xNpix << " x " << yNpix << " pixels", "magenta");
    imageRet.fillFrom(imageCut, Image::OPER_BIN);
  }

  //------------------------------------------------------------
  // If units were specified, override whatever we read from the FITS
  // file now
  //------------------------------------------------------------

  if(getParameter("units", false)->data_.hasValue_)
    imageRet.setUnits(getStringVal("units"));

  //------------------------------------------------------------  
  // If a smoothing kernel was specified, apply it now
  //------------------------------------------------------------  

  if(getParameter("sigmaapod", false)->data_.hasValue_) {
    double sigma = getDoubleVal("sigmaapod");
    std::string units = getParameter("sigmaapod", true)->units_;
    Angle sig;
    sig.setVal(sigma, units);
    Image gauss = imageRet;
    gauss.createGaussianImage(1.0, sig);
    imageRet *= gauss;

    COUTCOLOR("Apodizing with sigma = " << sigma << " " << units, "magenta");
  }

  if(getParameter("sigmasmooth", false)->data_.hasValue_) {
    double sigma = getDoubleVal("sigmasmooth");
    std::string units = getParameter("sigmasmooth", true)->units_;
    Angle sig;
    sig.setVal(sigma, units);
    Image gauss = imageRet;
    gauss.createGaussianImage(1.0, sig);
    imageRet.convolve(gauss);

    COUTCOLOR("Convolving with sigma = " << sigma << " " << units, "magenta");
  }

  //------------------------------------------------------------
  // If a transformation was specified, apply it now
  //------------------------------------------------------------

  if(getParameter("trans", false)->data_.hasValue_) {
    String trans = getStringVal("trans");
    trans.toLower();
    trans.strip(' ');
    
    if(trans == "sqrt") {
      imageRet.sqrt();
    } else if(trans == "ln") {
      imageRet.ln();
    } else if(trans == "log10") {
      imageRet.log10();
    } else if(trans.contains("^")) {
      double ex = trans.findNextInstanceOf("^", true, " ", false).toDouble();
      imageRet.power(ex);
    } else {
      XtermManip xtm;
      ThrowSimpleError(COLORIZE(xtm, "red", "Unsupported transform type: " << trans << std::endl << std::endl)
		       << COLORIZE(xtm, "green", "Supported transforms are: " << std::endl << std::endl
				   << "  sqrt  (Take the square root)" << std::endl
				   << "  ^ex   (Raise to the power of ex)" << std::endl
				   << "  log10 (Base-10 logarithm)" << std::endl
				   << "  ln    (Natural log)" << std::endl));
    }
  }

  return imageRet;
}

bool ImageManager::isImage(std::string fileName)
{
  FitsReader reader;
  reader.open(fileName);

  int naxis = reader.getLongKey("NAXIS");

  if(naxis == 2)
    return true;

  return false;
}

bool ImageManager::isEvent(std::string fileName)
{
  FitsBinTableReader reader(fileName);

  int naxis = reader.getLongKey("NAXIS");

  if(naxis == 0)
    return reader.hasTable("EVENTS");

  return false;
}

void ImageManager::getNpix(unsigned& nx, unsigned& ny)
{
  if(getParameter("npix", false)->data_.hasValue_) {
    String val(getStringVal("npix"));
    
    unsigned xNpix, yNpix;
    if(val.contains("x")) {
      String xnpix = val.findNextNumericString();
      val.advanceToNextNumericChar();
      String ynpix = val.findNextNumericString();
      
      nx = xnpix.toInt();
      ny = ynpix.toInt();
      
    } else {
      nx = val.toInt();
      ny = nx;
    }
  }
}

double ImageManager::getVal(Image& image, std::string par)
{
  double val;
  String str = getStringVal(par);
  str.strip(' ');

  std::ostringstream usage;
  usage << "Scale/offset can be a number, or one of the following symbolic names: " << std::endl << std::endl
	<< "  min     (the minimum of the image)" << std::endl
	<< "  mingz   (the minimum of the image above zero)" << std::endl
	<< "  max     (the maximum of the image)" << std::endl
	<< "  mean    (the mean of the image)" << std::endl
	<< "  rms     (the rms of the image)" << std::endl
	<< "  absmin  (the absolute values of the minimum)" << std::endl
	<< "  absmax  (the absolute values of the maximum)" << std::endl
	<< "  absmean (the absolute values of the mean)" << std::endl << std::endl
	<< "Any of the above can also be prefixed with '-' sign" << std::endl;

  if(str == "min")
    val = image.min();
  else if(str == "-min")
    val = -image.min();
  else if(str == "mingz")
    val = image.minGreaterThanZero();
  else if(str == "-mingz")
    val = -image.minGreaterThanZero();
  else if(str == "absmin")
    val = fabs(image.min());
  else if(str == "-absmin")
    val = -fabs(image.min());

  else if(str == "max")
    val = image.max();
  else if(str == "-max")
    val = -image.max();
  else if(str == "absmax")
    val = fabs(image.max());
  else if(str == "-absmax")
    val = -fabs(image.max());

  else if(str == "mean")
    val = image.mean();
  else if(str == "-mean")
    val = -image.mean();
  else if(str == "absmean")
    val = fabs(image.mean());
  else if(str == "-absmean")
    val = -fabs(image.mean());

  else if(str == "rms")
    val = image.rms();
  else if(str == "-rms")
    val = -image.rms();

  else if(str.isNumeric())
    val = str.toDouble();
  else {
    XtermManip xtm;
    ThrowSimpleError(COLORIZE(xtm, "red", "Invalid scale/offset specification: '" << str << "'" << std::endl << std::endl)
		     << COLORIZE(xtm, "green", usage.str()));
  }

  return val;
}

/**.......................................................................
 * Parse an expression of the form: 
 *
 *   x,y [units]
 */
void ImageManager::parseCoordinate(String& valStr, double& x, double& y, String& unit)
{
  String xStr = valStr.findNextInstanceOf("", false, ",", true, true);
  String yStr = valStr.findNextNumericString(":,.");

  x = xStr.toDouble();
  y = yStr.toDouble();

  unit = valStr.remainder();
}

/**.......................................................................
 * Parse an expression of the form
 *
 *   val [units]
 */
void ImageManager::parseVal(String& valStr, double& val, String& unit)
{
  val  = valStr.findNextNumericString(":,.").toDouble();
  unit = valStr.remainder();
}

/**.......................................................................
 * Parse a range expression of the form
 *
 *       val [units]
 *
 *  Or
 * 
 *   min:max [units]
 */
void ImageManager::parseRange(String& valStr, double& valMin, double& valMax, String& unit)
{
  //------------------------------------------------------------
  // See if this is of minVal:maxVal [units] form
  //------------------------------------------------------------

  if(valStr.contains(":")) {
    valMin = valStr.findNextInstanceOf("", false, ":", true, true).toDouble();
    valMax = valStr.findNextNumericString(":,.").toDouble();
    unit   = valStr.remainder();

    //------------------------------------------------------------
    // Or just 'val [units]'
    //------------------------------------------------------------

  } else {
    double val;
    parseVal(valStr, val, unit);
    valMin = val;
    valMax = val;
  }
}

/**.......................................................................
 * Parse an expression of the form
 *
 *   xrange, yrange [units]
 *
 * where xrange and yrange can be any valid 'range' specification
 */
void ImageManager::parseRanges(String& valStr, 
			       double& xMin, double& xMax, 
			       double& yMin, double& yMax, String& unit)
{
  String xStr = valStr.findNextInstanceOf("", false, ",", true, true);
  String yStr = valStr.remainder();

  parseRange(xStr, xMin, xMax, unit);
  parseRange(yStr, yMin, yMax, unit);
}

/**.......................................................................
 * Parse a region specification
 *
 * Valid region specifications are:
 *
 *   ixmin:ixmax, iymin:iymax;    // An explicit pixel region
 *   xmin:xmax, ymin:ymax units;  // A region specified in physical units
 *
 * Or:
 * 
 *   center +- delta; 
 *
 * Where center can be specified in pixels or physical units, i.e.:
 *
 *   512,512 +- delta;
 *   0.1, 0.2 degrees +- delta;
 *
 * Are both valid specifications, and delta can be specified in pixels
 * or physical units, as either a single delta, or separate deltas for
 * x and y, i.e.:
 *
 *   center +- 256;
 *   center +- 0.1 deg;
 *   center +- 256, 512;
 *   center +- 0.1, 0.2 deg;
 *
 */
Image ImageManager::extractRegion(String& region, Image& image)
{
  int ixmin, ixmax, iymin, iymax;

  std::ostringstream usage;
  usage << "Region should be specified as: " << std::endl << std::endl
	<< "  xmin:xmax,ymin:ymin [unit]" << std::endl << std::endl
	<< "or " << std::endl << std::endl
	<< "  xcenter,ycenter [unit] +- range [unit]" << std::endl << std::endl
	<< "where [] signifies an optional unit specification (e.g., 'degrees'), and range can be " << std::endl << std::endl
	<< "  val            (symmetric range for both x and y, about xcenter, ycenter)" << std::endl
	<< "  minval:maxval  (asymmetric range for both x and y, about xcenter, ycenter)" << std::endl
	<< "  xrange, yrange (separate ranges for x and y, specified as either of the above)" << std::endl << std::endl
	<< "If no units are specified, they will be assumed to be pixels.  "
	<< "Note that a requested range that exceeds the native image boundaries will be clipped to the actual data range.";

  try {

    //------------------------------------------------------------
    // Specified as x,y +- dx,dy ?
    //------------------------------------------------------------

    if(region.contains("+-")) {
      String center = region.findNextInstanceOf("", false, "+-", true, true);
      String range  = region.remainder();
      
      double xCenter, yCenter;
      String centerUnit;
      parseCoordinate(center, xCenter, yCenter, centerUnit);
      
      convertRange(range, image, ixmin, ixmax, iymin, iymax, true, xCenter, yCenter, centerUnit);
    } else {
      convertRange(region, image, ixmin, ixmax, iymin, iymax, false);
    }

  } catch(Exception& err) {
    XtermManip xtm;
    ThrowSimpleError(COLORIZE(xtm, "red", "Invalid range specification: '" << region << "'" << std::endl << std::endl)
		     << COLORIZE(xtm, "green", usage.str()));
  }

  return image.extract(ixmin, ixmax, iymin, iymax, true);
}

/**.......................................................................
 * Parse a region specified as: 
 *
 *   val [units]      In which case xmin = xmax = ymin = ymax = val
 *
 * Or 
 * 
 *  xrange, yrange [units]
 *
 * Where xrange or yrange can be:
 *
 *             val    In which case min = max = val
 *   minVal:maxVal    In which case min = minVal, max = maxVal
 *
 * And return the absolute pixel boundary corresponding to this region
 */
void ImageManager::convertRange(String& range, Image& image, 
				int& ixmin, int& ixmax, int& iymin, int& iymax,
				bool isRel, double xCenterVal, double yCenterVal, String centerUnit)
{
  String unit;
  double xMinVal, xMaxVal, yMinVal, yMaxVal;

  //------------------------------------------------------------
  // If there is no comma at all, then this is a single range
  //------------------------------------------------------------

  if(!range.contains(",")) {

    if(!isRel)
      ThrowSimpleColorError("Invalid absolute range: '" << range << "'", "red");

    parseRange(range, xMinVal, xMaxVal, unit);

    yMinVal = xMinVal;
    yMaxVal = xMaxVal;

    //------------------------------------------------------------
    // If there is a comma, then this is two ranges
    //------------------------------------------------------------

  } else {

    String xRange = range.findNextInstanceOf("", false, ",", true, false);
    String yRange = range.findNextInstanceOf(",", true, " ", false);

    if(!isRel && !xRange.contains(":"))
      ThrowSimpleColorError("Invalid absolute range: '" << range << "'", "red");

    if(!isRel && !yRange.contains(":"))
      ThrowSimpleColorError("Invalid absolute range: '" << range << "'", "red");

    parseRange(xRange, xMinVal, xMaxVal, unit);
    parseRange(yRange, yMinVal, yMaxVal, unit);
  }

  //------------------------------------------------------------
  // If no unit was specified, treat the range as pixels
  //------------------------------------------------------------
  
  if(unit.isEmpty()) {

    ixmin = (int) xMinVal;
    ixmax = (int) xMaxVal;
    iymin = (int) yMinVal;
    iymax = (int) yMaxVal;
      
    if(!isRel) {
      if(ixmin > ixmax || iymin > iymax)
	ThrowSimpleColorError("Invalid range specification: '" << range << "'", "red");
    }

    //------------------------------------------------------------
    // If this is relative to a center position, add the center now
    //------------------------------------------------------------

    if(isRel) {
      int ixCenter, iyCenter;

      //------------------------------------------------------------
      // Was center already specified in pixels?
      //------------------------------------------------------------
      
      if(centerUnit.isEmpty()) {

	ixCenter = (int) xCenterVal;
	iyCenter = (int) yCenterVal;
	
	//------------------------------------------------------------
	// Else center was in angular units -- convert to pixels
	//------------------------------------------------------------
	
      } else {

	Angle xCenter, yCenter;

	xCenter.setVal(xCenterVal, centerUnit.str());
	yCenter.setVal(yCenterVal, centerUnit.str());
	image.getPixelRef(xCenter, yCenter, ixCenter, iyCenter, false, false);
      }
      
      // Shift the boundaries relative to the center pixel

      ixmin = ixCenter - ixmin;
      ixmax = ixCenter + ixmax;

      iymin = iyCenter - iymin;
      iymax = iyCenter + iymax;
    }

    //------------------------------------------------------------
    // Else the range is in physical units -- convert to pixels
    //------------------------------------------------------------
    
  } else {

    Angle xMin, xMax, yMin, yMax;
    Angle xCenter, yCenter;
    
    //------------------------------------------------------------
    // Set the angular values of the range
    //------------------------------------------------------------

    xMin.setVal(xMinVal, unit.str());
    xMax.setVal(xMaxVal, unit.str());
    yMin.setVal(yMinVal, unit.str());
    yMax.setVal(yMaxVal, unit.str());

    if(isRel) {
      
      //------------------------------------------------------------
      // Was center specified in pixels?  If so, convert to angular
      // position first:
      //------------------------------------------------------------
      
      if(centerUnit.isEmpty()) {
	image.getPosition(xCenterVal, yCenterVal, xCenter, yCenter);

	//------------------------------------------------------------
	// Else center was specified in units
	//------------------------------------------------------------
	
      } else {
	xCenter.setVal(xCenterVal, centerUnit.str());
	yCenter.setVal(yCenterVal, centerUnit.str());
      }

      image.getPixelObj(xCenter - xMin, yCenter - yMin, ixmin, iymin, false, false);
      image.getPixelObj(xCenter + xMax, yCenter + yMax, ixmax, iymax, false, false);
      
      //------------------------------------------------------------
      // Else this is an absolute range
      //------------------------------------------------------------

    } else {
      image.getPixelRef(xMin, yMin, ixmin, iymin, false, false);
      image.getPixelRef(xMax, yMax, ixmax, iymax, false, false);
    }
  }
}
