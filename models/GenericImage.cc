#include "gcp/models/GenericImage.h"
#include "gcp/util/DataType.h"

using namespace std;

using namespace gcp::models;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
GenericImage::GenericImage() 
{
  addParameter("file",        DataType::STRING, "The file containing the image");
  addParameter("background",  DataType::DOUBLE, "Background count rate");
  addParameter("thetaMinErr", DataType::DOUBLE, "If specified, data outside of thetaMinErr will be used to estimate the background.  Note that this can conflict with 'region'");
  addParameter("errRegion",   DataType::STRING, "A region from which to estimate background.  Note that this can conflict with 'region'");
  addParameter("errImage",    DataType::STRING, "The image from which the background will be estimated, if thetaMinErr or errRegion are set.  Set to 'pre' to estimate from the original image, 'post' to estimate from the final image");
  addParameter("oper",        DataType::STRING, "The operation used to compare this model to datasets at a given offset.  Use 'nearest' to extract the nearest pixel value, or 'interpolate' to interpolate at the requested position (the default)");

  addParameter("display",     DataType::BOOL,   "Set to 'true' to display this model on readin");
  addParameter("interactive", DataType::BOOL,   "If true, make plots interactive");
  addParameter("dev",         DataType::STRING, "Pgplot device to use when displaying this model");
  addParameter("cmap",        DataType::STRING, "Colormap to use for displaying images: 'grey', 'heat' or 'rainbow'");

  addParameter(imageManager_);

  setParameter("errImage", "post");
  setParameter("oper",     "interpolate");

  normalization_.allowUnitless(true);

  initializeComponentsToFixed();

  imageInitialized_ = false;

  setupChecked_ = false;
}

/**.......................................................................
 * Destructor.
 */
GenericImage::~GenericImage() {}

/**.......................................................................
 * Unitless envelope function for this model -- just interpolate the
 * image at the requested offset
 */
double GenericImage::genericEnvelope(double xRad, double yRad)
{
  return envelope(xRad, yRad);
}

/**.......................................................................
 * Unitless envelope function for this model -- just interpolate the
 * image at the requested offset
 */
double GenericImage::radioEnvelope(double xRad, double yRad)
{
  return envelope(xRad, yRad);
}

/**.......................................................................
 * Unitless envelope function for this model -- just interpolate the
 * image at the requested offset
 */
double GenericImage::xrayImageEnvelope(double xRad, double yRad)
{
  return envelope(xRad, yRad);
}

/**.......................................................................
 * Unitless envelope function for this model -- just interpolate the
 * image at the requested offset
 */
double GenericImage::envelope(double xRad, double yRad)
{
  double val=0.0;
  bool valid = false;
  xOffCalc_.setRadians(xRad);
  yOffCalc_.setRadians(yRad);

  if(near_) {
    image_.getNearestData(xOffCalc_, yOffCalc_, val, valid, false);
  } else
    image_.interpolateData(xOffCalc_, yOffCalc_, val, valid, false);

  return val;
}

/**.......................................................................
 * Initialize the image
 */
void GenericImage::initializeImage()
{
  if(imageInitialized_)
    return;

  //------------------------------------------------------------
  // First read in the image from the named file
  //------------------------------------------------------------

  Image preImage;
  image_ = imageManager_.initializeImage(getStringVal("file"), preImage);

  double background;
  //------------------------------------------------------------
  // Throw if we can't determine the background
  //------------------------------------------------------------

  if(!getParameter("thetaMinErr", false)->data_.hasValue() && !getParameter("errRegion", false)->data_.hasValue()) {
    if(getParameter("background", false)->data_.hasValue())
      ThrowSimpleColorError("You must specify either a background via the 'background' parameter, or use 'thetaMinErr' or 'errRegion' to estimate the background from the data", "red");
  }

  //------------------------------------------------------------
  // Check if a fixed error value was specified
  //------------------------------------------------------------

  if(getParameter("background", false)->data_.hasValue()) {
    if(getParameter("thetaMinErr", false)->data_.hasValue() || getParameter("errRegion", false)->data_.hasValue())
      ThrowSimpleColorError("You can specify a fixed background value via the 'background' parameter, or use 'thetaMinErr' or 'errRegion' to estimate the background from the data, but not both", "red");

    background = getDoubleVal("background");  
  }

  //------------------------------------------------------------
  // Check if a minimum radius was specified for estimating the error.
  // If it was, estimate the error now.  The image we use is set by
  // the "errImage" parameter.
  //------------------------------------------------------------

  if(getParameter("thetaMinErr", false)->data_.hasValue_) {

    if(getParameter("errRegion", false)->data_.hasValue_)
      ThrowSimpleColorError("You can use 'thetaMinErr' or 'errRegion' to estimate the background from the data, but not both", "red");

    COUT("");

    Angle thetaMinErr;
    thetaMinErr.setVal(getDoubleVal("thetaMinErr"), getParameter("thetaMinErr", true)->units_);

    Image errImage;
    if(getStringVal("errImage") == "pre") {
      COUTCOLOR("Estimating background from pixels > " << thetaMinErr << " from the center of the original image", "magenta");
      errImage = preImage;
    } else {
      COUTCOLOR("Estimating background from pixels > " << thetaMinErr << " from the center of the extracted image", "magenta");
      errImage = image_;
    }

    Angle rmin, rmax;
    rmin.setVal(getDoubleVal("thetaMinErr"), getParameter("thetaMinErr", true)->units_);
    double xRad = errImage.xAxis().getAngularSize().radians();
    double yRad = errImage.yAxis().getAngularSize().radians();
    rmax.setRadians(2*sqrt(xRad*xRad + yRad*yRad));

    background = errImage.mean(rmin, rmax);
    COUTCOLOR("Estimating background to be         " << background, "magenta");

    COUT("");
  }

  //------------------------------------------------------------
  // Check if a region was specified for estimating the error.
  //------------------------------------------------------------

  if(getParameter("errRegion", false)->data_.hasValue_) {

    if(getParameter("thetaMinErr", false)->data_.hasValue_)
      ThrowSimpleColorError("You can use 'thetaMinErr' or 'errRegion' to estimate the background from the data, but not both", "red");

    String region = getStringVal("errRegion");

    COUT("");
    
    Image errImage;
    if(getStringVal("errImage") == "pre") {
      COUTCOLOR("Estimating background from region " << region << " of the original image", "magenta");
      errImage = imageManager_.extractRegion(region, preImage);
    } else {
      COUTCOLOR("Estimating background from region " << region << " of the extracted image", "magenta");
      errImage = imageManager_.extractRegion(region, image_);
    }

    background = errImage.mean();
    COUTCOLOR("Estimating background to be         " << background, "magenta");

    COUT("");
  }

  //------------------------------------------------------------
  // Now that we have determined the background, subtract it, and
  // normalize the result to unity at the peak
  //------------------------------------------------------------

  image_ -= background;
  image_ /= image_.max();

  image_.hasAbsolutePosition_ = false;

  imageInitialized_ = true;

  //------------------------------------------------------------
  // Finally, determine the operation we will use to extract data from
  // our image
  //------------------------------------------------------------

  String oper = getStringVal("oper");
  if(oper.contains("near"))
    near_ = true;
  else if(oper.contains("inter"))
    near_ = false;
  else {
    ThrowSimpleColorError("Unsupported operation type: '" << oper << "'.  Use either 'nearest' or 'interpolate'", "red");
  }
}

void GenericImage::difmapDisplay()
{
  image_.difmapDisplay();
}

void GenericImage::display()
{
  image_.display();
}

void GenericImage::checkSetup()
{
  if(setupChecked_)
    return;

  Generic2DAngularModel::checkSetup();

  initializeImage();

  displayIfRequested();

  setupChecked_ = true;
}

void GenericImage::displayIfRequested()
{  
  if(getParameter("display", false)->data_.hasValue()) {

    if(getBoolVal("display")) {

      if(getParameter("dev", false)->data_.hasValue())
	PgUtil::open(getStringVal("dev"));
      else
	PgUtil::open("/xs");

      PgUtil::setInteractive(getBoolVal("interactive"));

      image_.display();

      PgUtil::close();

      // Reset any legacy display settings

      PgUtil::setZmin(0.0);
      PgUtil::setZmax(0.0);
    }
  }
}
