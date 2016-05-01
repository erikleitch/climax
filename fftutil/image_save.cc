#include "gcp/fftutil/Dft2d.h"
#include "gcp/fftutil/Image.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/util/Astrometry.h"
#include "gcp/util/Constants.h"
#include "gcp/util/FitsImageReader.h"
#include "gcp/util/FitsImageWriter.h"
#include "gcp/util/String.h"
#include "gcp/util/SzCalculator.h"
#include "gcp/util/Planck.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace std;

using namespace gcp::util;

const Angle Image::zero_ = Angle(Angle::Degrees(), 0.0);

/**.......................................................................
 * Constructors
 */
Image::Image() 
{
  initialize();
}

Image::Image(unsigned npix, gcp::util::Angle size)
{
  initialize();

  xAxis().setNpix(npix);
  yAxis().setNpix(npix);

  xAxis().setAngularSize(size);
  yAxis().setAngularSize(size);
}

Image::Image(const Image& image) 
{
  *this = (Image&) image;
}

Image::Image(Image& image) 
{
  *this = image;
}

Image::Image(std::vector<float>& data, Unit::Units units)
{
  initializeFromArray(data, units);
}

void Image::assignDataFrom(const Image& image)
{
  assignDataFrom((Image&)image);
}

void Image::assignDataFrom(Image& image)
{
  if(image.hasData_) {
    if(data_.size() != image.data_.size()) {
      ThrowColorError("Cannot assign data from images of different sizes", "red");
    }
    
    data_ = image.data_;
    hasData_ = true;
  }
}


void Image::operator=(const Image& image) 
{
  *this = (Image&)image;
}

void Image::operator=(Image& image) 
{
  initialize();

  data_.resize(image.data_.size());

  for(unsigned i=0; i < data_.size(); i++)
    data_[i] = image.data_[i];

  xAxis_ = image.xAxis_;
  yAxis_ = image.yAxis_;
  
  units_         = image.units_;
  nativeToJy_    = image.nativeToJy_;

  hasData_       = image.hasData_;
  hasUnits_      = image.hasUnits_;
  hasNativeToJy_ = image.hasNativeToJy_;

  dataMin_       = image.dataMin_;
  dataMax_       = image.dataMax_;
}

void Image::initialize()
{
  data_.resize(0);
  units_               = Unit::UNITS_UNKNOWN;
  hasData_             = false;
  hasUnits_            = false;
  hasNativeToJy_       = false;
  hasAbsolutePosition_ = false;

  dataMin_       = 0.0;
  dataMax_       = 0.0;

  axes_ = ImageAxis::AXIS_NONE;

  xAxis_.setAxisType(Axis::AXIS_X);
  yAxis_.setAxisType(Axis::AXIS_Y);

  xAxis_.parent_ = this;
  yAxis_.parent_ = this;

  // Explicitly zero the data array

  zero();
}

/**.......................................................................
 * Destructor.
 */
Image::~Image() {}


/**.......................................................................
 * Fill this image with a uniform value
 */
void Image::createUniformImage(unsigned nx, unsigned ny, float value)
{
  resize(nx, ny);

  for(unsigned i=0; i < nx*ny; i++) {
    data_[i] = value;
  }

  hasData_ = true;
}

/**.......................................................................
 * Fill this image with a border
 */
void Image::createEdgeImage(unsigned nx, unsigned ny)
{
  unsigned ind;

  resize(nx, ny);

  for(unsigned i=0; i < nx*ny; i++) {
    data_[i] = 0.0;
  }

  for(unsigned i=0; i < nx; i++) {
    ind = (0) * nx + i;
    data_[ind] = 1.0;
 
    ind = (ny-1) * nx + i;
    data_[ind] = 1.0;
  }

  for(unsigned j=0; j < ny; j++) {
    ind = j * nx + (0);
    data_[ind] = 1.0;
 
    ind = j * nx + (nx-1);
    data_[ind] = 1.0;
  }

  hasData_ = true;
}

/**.......................................................................
 * Fill this image with a delta function
 */
void Image::createDeltaFunctionImage(unsigned nx, unsigned ny, float val, int xval, int yval)
{
  resize(nx, ny);

  for(unsigned i=0; i < nx*ny; i++) {
    data_[i] = 0.0;
  }
 
  // Default to center, if no location was given

  if(xval < 0)
    xval = nx/2;

  if(yval < 0)
    yval = ny/2;

  // Install a point source at the center of the image

  unsigned ind = (yval) * nx + xval;
  data_[ind] = val;

  hasData_ = true;
}

/**.......................................................................
 * Fill this image with a delta function
 */
void Image::createDeltaFunctionImage(float val, Angle xAngle, Angle yAngle)
{
  double dxRad    = xAxis().getAngularResolution().radians();
  double dyRad    = yAxis().getAngularResolution().radians();
  unsigned nx     = xAxis().getNpix();
  unsigned ny     = yAxis().getNpix();

  for(unsigned i=0; i < nx*ny; i++) {
    data_[i] = 0.0;
  }
 
  // Find the nearest pixel to the requested position

  int iX = (int)round(xAngle.radians() / dxRad) + nx/2;
  int iY = (int)round(yAngle.radians() / dxRad) + ny/2;

  if(iX < 0 || iX > nx-1) {
    ThrowError("Requested position lies outside the range covered by this image");
  }

  if(iY < 0 || iY > ny-1) {
    ThrowError("Requested position lies outside the range covered by this image");
  }

  // Install a point source at the requested location

  unsigned ind = (iY) * nx + iX;

  data_[ind] = val;

  hasData_ = true;
}
/**.......................................................................
 * Fill this image with a gaussian (sigma in pixels)
 */
void Image::createGaussianImage(unsigned nx, unsigned ny, double sigma)
{
  resize(nx, ny);
  
  double xval, yval, arg;
  for(int iy=0; iy < ny; iy++) {
    for(int ix=0; ix < nx; ix++) {

      xval = (double)(ix) - (double)(nx/2);
      yval = (double)(iy) - (double)(ny/2);
      arg = (xval * xval + yval * yval)/(2*sigma*sigma);
      data_[iy * nx + ix] = exp(-arg);
    }
  }

  hasData_ = true;
}

/**.......................................................................
 * Fill this image with a gaussian (sigma in pixels)
 */
void Image::createGaussianImage(unsigned nx, unsigned ny, double sigmax, double sigmay)
{
  resize(nx, ny);
  
  double xval, yval, arg;
  for(int iy=0; iy < ny; iy++) {
    for(int ix=0; ix < nx; ix++) {

      xval = (double)(ix) - (double)(nx/2);
      yval = (double)(iy) - (double)(ny/2);
      arg = (xval * xval)/(2*sigmax*sigmax) + (yval * yval)/(2*sigmay*sigmay);
      data_[iy * nx + ix] = exp(-arg);
    }
  }

  hasData_ = true;
}

/**.......................................................................
 * Fill this image with a gaussian (sigma in Angle)
 */
void Image::createGaussianImage(double amp, Angle sigma, Angle xOffset, Angle yOffset)
{
  double sigmaRad = sigma.radians();
  double dxRad    = xAxis().getAngularResolution().radians();
  double dyRad    = yAxis().getAngularResolution().radians();
  unsigned nx     = xAxis().getNpix();
  unsigned ny     = yAxis().getNpix();

  double xval, yval, arg;
  for(int iy=0; iy < ny; iy++) {
    for(int ix=0; ix < nx; ix++) {

      xval = (double)(ix) - (double)(nx/2);
      yval = (double)(iy) - (double)(ny/2);

      xval *= dxRad;
      yval *= dyRad;

      xval -= xOffset.radians();
      yval -= yOffset.radians();

      arg = (xval * xval + yval * yval)/(2*sigmaRad*sigmaRad);
      data_[iy * nx + ix] = amp * exp(-arg);
    }
  }

  hasData_ = true;
}

Angle Image::gaussianSigma(Length& diameter, Frequency& freq)
{
  double radians = (1.02 * freq.centimeters() / diameter.centimeters()) / sqrt(8*log(2.0));
  Angle sigma(Angle::Radians(), radians);
  return sigma;
}

/**.......................................................................
 * Fill this image with an approximate primary beam for an antenna of
 * given diameter, at the requested frequency
 */
void Image::createGaussianPrimaryBeam(Length& diameter, Frequency& freq)
{
  double sigma = gaussianSigma(diameter, freq).radians();

  double dx    = xAxis_.getAngularResolution().radians();
  double dy    = yAxis_.getAngularResolution().radians();

  unsigned nx = xAxis_.getNpix();
  unsigned ny = yAxis_.getNpix();

  double xval, yval, arg;

  for(int iy=0; iy < ny; iy++) {
    for(int ix=0; ix < nx; ix++) {

      xval = ((double)(ix) - (double)(nx/2)) * dx;
      yval = ((double)(iy) - (double)(ny/2)) * dy;
      arg = (xval * xval + yval * yval)/(2*sigma*sigma);
      data_[iy * nx + ix] = exp(-arg);
    }
  }

  hasData_ = true;
}

/**.......................................................................
 * Add a 2-D cosine function (period in pixels) to this image
 */
void Image::addCosine(double period, Angle phase)
{
  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();

  double xval, yval, arg;
  double xpval, ypval;
  double cp = cos(phase.radians());
  double sp = sin(phase.radians());

  for(int iy=0; iy < ny; iy++) {
    for(int ix=0; ix < nx; ix++) {

      xval = (double)(ix) - (double)(nx/2);
      yval = (double)(iy) - (double)(ny/2);

      xpval =   xval * cp + yval * sp;
      ypval = - xval * sp + yval * cp;

      data_[iy * nx + ix] += cos(2*M_PI*xpval/period);
    }
  }
}

/**.......................................................................
 * Fill this image with a 2-D cosine function (period in pixels)
 */
void Image::createCosineImage(unsigned nx, unsigned ny, double period, Angle phase)
{
  createUniformImage(nx, ny, 0.0);
  addCosine(period, phase);

  hasData_ = true;
}

/**.......................................................................
 * Fill this image with a 2-D sine function (period in pixels)
 */
void Image::createSineImage(unsigned nx, unsigned ny, double period, Angle phase)
{
  createUniformImage(nx, ny, 0.0);
  addSine(period, phase);

  hasData_ = true;
}

/**.......................................................................
 * Add a 2D sine function to this image
 */
void Image::addSine(double period, Angle phase)
{
  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();

  double xval, yval, arg;
  double xpval, ypval;
  double cp = cos(phase.radians());
  double sp = sin(phase.radians());

  for(int iy=0; iy < ny; iy++) {
    for(int ix=0; ix < nx; ix++) {

      xval = (double)(ix) - (double)(nx/2);
      yval = (double)(iy) - (double)(ny/2);

      xpval =   xval * cp + yval * sp;
      ypval = - xval * sp + yval * cp;

      data_[iy * nx + ix] += sin(2*M_PI*xpval/period);
    }
  }
}

void Image::initialize(Dft2d& dft)
{
  *this = dft.getImage();
}

/**.......................................................................
 * Initialize the data in this image from an external data array
 */
void Image::initializeFromArray(std::vector<float>& data, Unit::Units units)
{
  data_.resize(data.size());

  for(unsigned i=0; i < data.size(); i++) {
    data_[i] = data[i];
  }

  setUnits(units);
  hasData_ = true;
}

/**.......................................................................
 * Initialize the data in this image from an external data array
 */
void Image::initializeFromArray(double* data, unsigned ndata, Unit::Units units)
{
  data_.resize(ndata);

  for(unsigned i=0; i < ndata; i++) {
    data_[i] = data[i];
  }

  setUnits(units);
  hasData_ = true;
}

/**.......................................................................
 * Initialize the data in this image from a FITS image file
 */
void Image::initializeFromFitsFile(std::string fileName)
{
  gcp::util::FitsImageReader reader(fileName);
  FitsImageReader::Image image;
  reader.readData(image);

  FitsImageReader::Axis& xFitsAxis = reader.axes_[0];
  FitsImageReader::Axis& yFitsAxis = reader.axes_[1];

  resize(xFitsAxis.n_, yFitsAxis.n_);

  String xAxisType(xFitsAxis.type_);
  String xAxisTypeComment(xFitsAxis.typeComment_);

  if(xAxisType.contains("RA---SIN") || xAxisTypeComment.contains("degrees")) {
    xAxis().setAngularSize(Angle(Angle::Degrees(), xFitsAxis.n_ * xFitsAxis.delta_));
    ra_.setDegrees(xFitsAxis.refVal_);
    hasAbsolutePosition_ = true;
    COUT("RA is now: " << ra_);
  }

  String yAxisType(yFitsAxis.type_);
  String yAxisTypeComment(yFitsAxis.typeComment_);

  if(yAxisType.contains("DEC--SIN") || yAxisTypeComment.contains("degrees")) {
    yAxis().setAngularSize(Angle(Angle::Degrees(), yFitsAxis.n_ * yFitsAxis.delta_));
    dec_.setDegrees(yFitsAxis.refVal_);
    hasAbsolutePosition_ = true;
    COUT("DEC is now: " << dec_);
  }

  reader.printAxes();
  reader.printHeaderCards();

  for(unsigned i=0; i < data_.size(); i++) {
    data_[i] = image.data_[i];
  }

  hasData_ = true;

  units_ = Unit::stringToUnits(reader.units());

  if(units_ != Unit::UNITS_UNKNOWN) {
    hasUnits_ = true;
  }
}

/**.......................................................................
 * Write this image to a FITS file
 */
void Image::writeToFitsFile(std::string fileName)
{
  gcp::util::FitsImageWriter writer(fileName);

  std::vector<long> axisLengths(2);

  axisLengths[0] = xAxis().getNpix();
  axisLengths[1] = yAxis().getNpix();

  //------------------------------------------------------------
  // Write BITPIX and NAXIS, etc.
  //------------------------------------------------------------

  writer.writeStandardKeys(-32, axisLengths);

  //------------------------------------------------------------
  // Write random FITS keywords that seem standard
  //------------------------------------------------------------

  int pcount = 0;
  writer.putKey("PCOUNT", pcount, "Parameter count");

  int gcount = 1;
  writer.putKey("GCOUNT", gcount, "Group count");

  //------------------------------------------------------------
  // Now write keywords for the X-axis
  //------------------------------------------------------------

  if(xAxis().hasAngularSize()) {
    std::string xtype("X");
    writer.putKey("CTYPE1", xtype, "Axis name (degrees)");
    float crpix = (float)(xAxis().getNpix())/2;
    writer.putKey("CRPIX1", crpix, "Reference pixel");
    float crval = 0.0;
    writer.putKey("CRVAL1", crval, "Reference value (degrees)");
    float crdelt = xAxis().getAngularResolution().degrees();
    writer.putKey("CDELT1", crdelt, "Pixel increment (degrees)");
    float crota = 0.0;
    writer.putKey("CROTA1", crota, "Axis rotation");

  } else {
    std::string xtype("X");
    writer.putKey("CTYPE1", xtype, "Axis name (degrees)");
    float crpix = (float)(xAxis().getNpix())/2;
    writer.putKey("CRPIX1", crpix, "Reference pixel");
    float crval = (float)(xAxis().getNpix())/2;
    writer.putKey("CRVAL1", crval, "Reference value (pixels)");
    float crdelt = 1.0;
    writer.putKey("CDELT1", crdelt, "Pixel increment");
    float crota = 0.0;
    writer.putKey("CROTA1", crota, "Axis rotation");
  }

  // Now the Y-axis

  if(yAxis().hasAngularSize()) {
    std::string ytype("Y");
    writer.putKey("CTYPE2", ytype, "Axis name (degrees)");
    float crpix = (float)(yAxis().getNpix())/2;
    writer.putKey("CRPIX2", crpix, "Reference pixel");
    float crval = 0.0;
    writer.putKey("CRVAL2", crval, "Reference value (degrees)");
    float crdelt = yAxis().getAngularResolution().degrees();
    writer.putKey("CDELT2", crdelt, "Pixel increment (degrees)");
    float crota = 0.0;
    writer.putKey("CROTA2", crota, "Axis rotation");

  } else {
    std::string ytype("Y");
    writer.putKey("CTYPE2", ytype, "Axis name (degrees)");
    float crpix = (float)(yAxis().getNpix())/2;
    writer.putKey("CRPIX2", crpix, "Reference pixel");
    float crval = (float)(yAxis().getNpix())/2;
    writer.putKey("CRVAL2", crval, "Reference value (pixels)");
    float crdelt = 1.0;
    writer.putKey("CDELT2", crdelt, "Pixel increment");
    float crota = 0.0;
    writer.putKey("CROTA2", crota, "Axis rotation");
  }

  // Write the units, if they are known

  writer.putKey("BUNIT", unitToString(), "Unit of measurement");

  // Now write the data

  writer.writeData(data_);
}

/**.......................................................................
 * Initialize the data in this image from a binary image file
 */
void Image::initializeFromBinaryImageFile(std::string fileName)
{
  int fd = -1;

  fd = open(fileName.c_str(), O_RDONLY);

  if(fd < 0) {
    COUT("Unable to open file: " << fileName);
  }

  struct stat buf;
  fstat(fd, &buf);

  unsigned npix = (unsigned)sqrt((double)(buf.st_size/sizeof(float)));

  resize(npix, npix);

  read(fd, (void*)&data_[0], npix*npix*sizeof(float));

  hasData_ = true;

  close(fd);
}

void Image::resize()
{
  data_.resize(xAxis().getNpix() * yAxis().getNpix());
}

void Image::resize(unsigned nx, unsigned ny)
{
  xAxis().setNpix(nx);
  yAxis().setNpix(ny);
}

/**.......................................................................
 * Initialize this object from another Image object.  This method initializes the size and angular resolution to match.
 */
void Image::initialize(Image& image)
{
  initialize();
  xAxis().setNpix(image.xAxis().getNpix());
  yAxis().setNpix(image.yAxis().getNpix());
  xAxis().setAngularSize(image.xAxis().getAngularSize());
  yAxis().setAngularSize(image.yAxis().getAngularSize());
}

/**.......................................................................
 * Determine the conversion from native image units to Jy
 */
void Image::convertToJy(Frequency& nu)
{
  (*this) *= nativeToJy(nu);
  setUnits(Unit::UNITS_JY);
}

/**.......................................................................
 * Determine the conversion from native image units to Jy
 */
double Image::nativeToJy(Frequency& nu)
{
  double nativeToJy;
  double dxr    = xAxis_.getAngularResolution().radians();
  double dyr    = yAxis_.getAngularResolution().radians();

  switch (getUnits()) {

    // If the image is already in Jy, we don't have to do anything

  case Unit::UNITS_JY:
    nativeToJy = 1.0;
    break;

  case Unit::UNITS_MJYSR:
    nativeToJy = 1e6 * dxr * dyr;
    break;

    // If the image is in microK, we have to convert to Planck
    // intensity, and scale by the pixel size

  case Unit::UNITS_UK:

    nativeToJy = Planck::JyPerSrPerKPlanck(nu, Constants::Tcmb_) 
      * dxr * dyr * 1e-6;
    break;

  case Unit::UNITS_K:

    nativeToJy = Planck::JyPerSrPerKPlanck(nu, Constants::Tcmb_) 
      * dxr * dyr;
    break;

  case Unit::UNITS_Y:
    {
      Intensity intensity = SzCalculator::comptonYToIntensity(nu);
      nativeToJy = intensity.JyPerSr() * dxr * dyr;
    }
    break;

  default:
    ThrowError("Unrecognized units: " << units_);
    break;
  }

  return nativeToJy;
}

/**.......................................................................
 * Determine the conversion from native image units to Jy
 */
void Image::setNativeToJy(Frequency& nu)
{
  setNativeToJy(nativeToJy(nu));
}

/**.......................................................................
 * Set a conversion factor from native units to Jansky
 */
void Image::setNativeToJy(double jyconv)
{
  nativeToJy_ = jyconv;
  hasNativeToJy_ = true;
}

/**.......................................................................
 * Specifiy the units of this image's data
 */
void Image::setUnits(Unit::Units units)
{
  units_    = units;
  hasUnits_ = true;
}

/**.......................................................................
 * Return the units of this image's data
 */
Unit::Units Image::getUnits()
{
  if(hasUnits_) {
    return units_;
  } else {
    ThrowError("No units have been specified for this image");
  }
}

/**.......................................................................
 * Display this image
 */
void Image::display()
{
  if(!xAxis().hasNpix()) {
    ThrowError("x-axis length hasn't been specified");
  }

  if(!yAxis().hasNpix()) {
    ThrowError("y-axis length hasn't been specified");
  }

  std::string unit = unitToString();
  std::string xlab = "Pixel Index";
  std::string ylab = "Pixel Index";

  double xmin = 0;
  double xmax = xAxis().getNpix()-1;
  
  if(xAxis().hasAngularSize()) {
    xmin = -xAxis().getAngularSize().degrees()/2;
    xmax =  xAxis().getAngularSize().degrees()/2;
    xlab = "Degrees";
  }

  double ymin = 0;
  double ymax = yAxis().getNpix()-1;

  if(yAxis().hasAngularSize()) {
    ymin = -yAxis().getAngularSize().degrees()/2;
    ymax =  yAxis().getAngularSize().degrees()/2;
    ylab = "Degrees";
  }
  

  PgUtil::greyScale(data_.size(), &data_[0], xAxis().getNpix(), yAxis().getNpix(),
		    xmin, xmax, ymin, ymax,
		    0, 0, 0,
		    (char*)xlab.c_str(), (char*)ylab.c_str(), "", (char*)unit.c_str());
  //  PgUtil::wnad(0, 1, 0, 1);
}

/**.......................................................................
 * Return the rms of the image data
 */
double Image::rms()
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  double imMean = mean();
  double imSdev = 0.0;
  double val;
  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();
  unsigned n = nx * ny;

  for(unsigned i=0; i < nx * ny; i++) {
    val = data_[i] - imMean;
    imSdev += (val * val - imSdev) / (i + 1);
  }

  imSdev = (n > 1) ? sqrt(n * imSdev / (n - 1)) : 0.0;

  return imSdev;
}

/**.......................................................................
 * Return the max of the image data
 */
double Image::max()
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  double imMax = 0.0;
  double val;
  for(unsigned i=0; i < data_.size(); i++) {
    val = data_[i];
    imMax = (val > imMax) ? val : imMax;
  }

  return imMax;
}

/**.......................................................................
 * Return the min of the image data
 */
double Image::min()
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  double imMin = 0.0;
  double val;
  for(unsigned i=0; i < data_.size(); i++) {
    val = data_[i];
    imMin = (val < imMin) ? val : imMin;
  }

  return imMin;
}

/**.......................................................................
 * Return the mean of the image data
 */
double Image::mean()
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  double imMean = 0.0;
  for(unsigned i=0; i < data_.size(); i++) {
    imMean += (data_[i] - imMean) / (i + 1);
  }

  return imMean;
}

/**.......................................................................
 * Return the sum of the image data
 */
double Image::sum()
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  double imSum = 0.0;
  for(unsigned i=0; i < data_.size(); i++) {
    imSum += data_[i];
  }

  return imSum;
}

//-----------------------------------------------------------------------
// Arithmetic operators on image data -- operators with constants
//-----------------------------------------------------------------------

void Image::operator+=(double incr)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  data_ += incr;
}

void Image::operator-=(double decr)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  data_ -= decr;
}

void Image::operator*=(double mult)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  data_ *= mult;
}

void Image::operator/=(double div)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  data_ /= div;
}

void Image::operator=(double val)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  data_ = val;
}

Image Image::operator+(double incr)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  Image ret = *this;
  ret += incr;
  
  return ret;
}

Image Image::operator-(double decr)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  Image ret = *this;
  ret -= decr;
  
  return ret;
}

Image Image::operator*(double mult)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  Image ret = *this;
  ret *= mult;
  
  return ret;
}

Image Image::operator/(double div)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  Image ret = *this;
  ret /= div;
  
  return ret;
}


//-----------------------------------------------------------------------
// Arithmetic operators on image data -- operators with Images
//-----------------------------------------------------------------------

void Image::operator+=(Image& image)
{
  if(!hasData_ || !image.hasData_) {
    ThrowError("Image contains no data");
  }

  if(image.hasData_) {
    if(hasAbsolutePosition_ && image.hasAbsolutePosition_) {
      addAbsolute(image);
      COUT("Adding absolute");
    } else {
      data_ += image.data_;
    }

    hasData_ = true;
  }
}

void Image::operator-=(Image& image)
{
  if(!hasData_ || !image.hasData_) {
    ThrowError("Image contains no data");
  }

  data_ -= image.data_;
}

void Image::operator*=(const Image& image) {
  return operator*=((Image&)image);
}

void Image::operator*=(Image& image) {

  if(!hasData_) {
    ThrowError("This image contains no data");
  }

  if(!image.hasData_) {
    ThrowError("That image contains no data");
  }

  data_ *= image.data_;
}

void Image::operator/=(Image& image)
{
  if(!hasData_ || !image.hasData_) {
    ThrowError("Image contains no data");
  }

  data_ /= image.data_;
}

Image Image::operator+(Image& image)
{
  if(!hasData_ || !image.hasData_) {
    ThrowError("Image contains no data");
  }

  Image ret = *this;
  ret += image;

  return ret;
}

Image Image::operator-(Image& image)
{
  if(!hasData_ || !image.hasData_) {
    ThrowError("Image contains no data");
  }

  Image ret = *this;
  ret -= image;
  
  return ret;
}

Image Image::operator*(Image& image)
{
  if(!hasData_ || !image.hasData_) {
    ThrowError("Image contains no data");
  }

  Image ret = *this;
  ret *= image;

  return ret;
}

Image Image::operator/(Image& image)
{
  if(!hasData_ || !image.hasData_) {
    ThrowError("Image contains no data");
  }

  Image ret = *this;
  ret /= image;

  return ret;
}

/**.......................................................................
 * Convolve this image with another one
 */
void Image::convolve(Image& image, bool zeropad)
{
  if(!hasData_ || !image.hasData_) {
    ThrowError("Image contains no data");
  }

  Dft2d dft1, dft2;

  if(zeropad) {
    dft1.xAxis().setZeropadFactor(2);
    dft1.yAxis().setZeropadFactor(2);
    dft2.xAxis().setZeropadFactor(2);
    dft2.yAxis().setZeropadFactor(2);
  }

  dft1.initialize(*this);
  dft2.initialize(image);

  dft1.normalize(true);

  dft1.computeForwardTransform();
  dft2.computeForwardTransform();

  dft1.complexMultiply(dft2, false);
  dft1.shift();
  dft1.computeInverseTransform();

  *this = dft1.getImage(zeropad);
}

/**.......................................................................
 * Plot the specified row of this image
 */
void Image::plotRow(unsigned iRow)
{
  unsigned ind;
  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();

  std::vector<float> xarr(nx);
  std::vector<float> yarr(ny);

  for(unsigned iCol=0; iCol < ny; iCol++) {
    ind  = iCol * ny + iRow;
    xarr[iCol] = (float)iCol;
    yarr[iCol] = data_[ind];
  }

  PgUtil::linePlot(xarr.size(), &xarr[0], &yarr[0], 0, "Column index", "Data", "");
}

/**.......................................................................
 * Plot the specified column of this image
 */
void Image::plotColumn(unsigned iCol)
{
  unsigned ind;
  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();

  std::vector<float> xarr(nx);
  std::vector<float> yarr(ny);

  for(unsigned iRow=0; iRow < ny; iRow++) {
    ind  = iCol * ny + iRow;
    xarr[iRow] = (float)iRow;
    yarr[iRow] = data_[ind];
  }

  PgUtil::linePlot(xarr.size(), &xarr[0], &yarr[0], 0, "Row index", "Data", "");
}

/**.......................................................................
 * Return the x-axis descriptor of this image
 */
Image::Axis& Image::xAxis()
{
  return xAxis_;
}

/**.......................................................................
 * Return the y-axis descriptor of this image
 */
Image::Axis& Image::yAxis()
{
  return yAxis_;
}

std::string Image::unitToString()
{
  switch (units_) {
  case Unit::UNITS_JY:
    return "Jy";
    break;
  case Unit::UNITS_MJYSR:
    return "MJy/sr";
    break;
  case Unit::UNITS_UK:
    return "Micro-Kelvin";
    break;
  case Unit::UNITS_K:
    return "Kelvin";
    break;
  case Unit::UNITS_Y:
    return "Compton Y";
    break;
  default:
    return "(Unknown units)";
    break;
  }
}

/**.......................................................................
 * Set the number of pixels in this axis
 */
void Image::Axis::setNpix(unsigned n)
{
  ImageAxis::setNpix(n);
  parent_->updateAxis(type_);
}

void Image::updateAxis(unsigned axis)
{
  unsigned axesMask = (unsigned)axes_;
  unsigned axisMask = (unsigned)axis;
  axes_ = (axesMask | axisMask);

  checkAxes();
}

void Image::checkAxes()
{
  if(axes_ == ImageAxis::AXIS_BOTH) {
    resize();
  }
}

void Image::Axis::operator=(const ImageAxis& axis)
{
  return operator=((ImageAxis&) axis);
}

void Image::Axis::operator=(ImageAxis& axis)
{
  if(axis.hasAngularSize())
    setAngularSize(axis.getAngularSize());

  if(axis.hasNpix())
    setNpix(axis.getNpix());
}

void Image::setNpix(unsigned npix)
{
  xAxis().setNpix(npix);
  yAxis().setNpix(npix);
}

void Image::setAngularSize(Angle size)
{
  xAxis().setAngularSize(size);
  yAxis().setAngularSize(size);
}

void Image::zero()
{
  data_ = 0.0;
}

void Image::setValue(float val)
{
  data_ = val;
}

float& Image::val(unsigned ix, unsigned iy)
{
  if(ix > xAxis_.getNpix()-1 || iy > yAxis_.getNpix()-1) {
    ThrowError("Invalid pixel index: " << ix << "," << iy);
  }

  unsigned imInd = iy * xAxis_.getNpix() + ix;

  return data_[imInd];
}

bool Image::hasData()
{
  return hasData_;
}

bool Image::axesAreEquivalent(Image& image)
{
  return (xAxis_ == image.xAxis_) && (yAxis_ == image.yAxis_);
}

void Image::createGaussianImageFullSpecification(double amp, 
						 Angle majSigma, Angle minSigma, 
						 Angle rotAngle,
						 Angle xOffset, Angle yOffset)
{
  unsigned nx     = xAxis().getNpix();
  unsigned ny     = yAxis().getNpix();

  double dxRad    = xAxis().getAngularResolution().radians();
  double dyRad    = yAxis().getAngularResolution().radians();

  double cRotAng  = cos(rotAngle.radians());
  double sRotAng  = sin(rotAngle.radians());

  double xOffRad  = xOffset.radians();
  double yOffRad  = yOffset.radians();

  double minSigRad2 = minSigma.radians() * minSigma.radians();
  double majSigRad2 = majSigma.radians() * majSigma.radians();

  double x, y, xp, yp, xpp, ypp, arg;
  for(int iy=0; iy < ny; iy++) {
    for(int ix=0; ix < nx; ix++) {

      x = (double)(ix) - (double)(nx/2);
      y = (double)(iy) - (double)(ny/2);

      x *= dxRad;
      y *= dyRad;

      // Coordinate in the untranslated frame

      xp = x - xOffRad;
      yp = y - yOffRad;

      // Coordinate in the unrotated frame

      xpp =   xp * cRotAng + yp * sRotAng;
      ypp = - xp * sRotAng + yp * cRotAng;

      arg = 0.5 * ((xpp * xpp)/majSigRad2 + (ypp * ypp)/minSigRad2);

      data_[iy * nx + ix] = amp * exp(-arg);
    }
  }

  hasData_ = true;
  
}

void Image::setHasData(bool hasData)
{
  hasData_ = hasData;
}

/**.......................................................................
 * Interpolate this image at the requested point
 */
void Image::interpolateData(Angle& xOff, Angle& yOff, double& val, bool& valid)
{
  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();

  // First get the nearest pixel to the requested point

  int ixNear = (int)(xOff / xAxis().getAngularResolution()) + nx/2;
  int iyNear = (int)(yOff / yAxis().getAngularResolution()) + ny/2;

  unsigned convMaskInPixels = Dft2d::convMaskInPixels_;
  double   convSigInPixels  = Dft2d::convSigInPixels_;

  val = 0.0;

  if(ixNear < convMaskInPixels || (xAxis().getNpix()-1 - ixNear) < convMaskInPixels) {
    valid = false;
    return;
  }

  if(iyNear < convMaskInPixels || (yAxis().getNpix()-1 - iyNear) < convMaskInPixels) {
    valid = false;
    return;
  }

  unsigned ixStart = ixNear - convMaskInPixels;
  unsigned ixStop  = ixNear + convMaskInPixels;

  unsigned iyStart = iyNear - convMaskInPixels;
  unsigned iyStop  = iyNear + convMaskInPixels;

  double dxRad = xAxis().getAngularResolution().radians();
  double dyRad = yAxis().getAngularResolution().radians();

  // Get the requested position in pixel coordinates

  double xReqPix = xOff.radians() / dxRad;
  double yReqPix = yOff.radians() / dyRad;

  // Finally, initialize sums

  double wtSum = 0.0;
  double s2    = 2 * convSigInPixels * convSigInPixels;

  // Finally, accumulate the running mean over the pixel mask

  for(unsigned ix = ixStart; ix <= ixStop; ix++) {
    for(unsigned iy = iyStart; iy <= iyStop; iy++) {

      unsigned imInd = iy * nx + ix;

      // Calculate the coordinate of the center of this pixel

      double xPix = (double)(ix) - double(nx/2);
      double yPix = (double)(iy) - double(ny/2);
      
      // Calculate the delta of the requested point from the pixel
      // center

      double dxPix = xReqPix - xPix;
      double dyPix = yReqPix - yPix;

      double arg = (dxPix * dxPix + dyPix * dyPix)/s2;

      double wt = exp(-arg);

      // Accumulate running mean
	
      double pixVal = data_[imInd];

      val += ((pixVal - val) * wt) / (wtSum + wt);

      wtSum += wt;
    }
  }

  valid = true;
  return;
}

/**.......................................................................
 * Method to set an absolute position for this image
 */
void Image::setRaDec(HourAngle& ra, Declination& dec)
{
  ra_  = ra;
  dec_ = dec;
  hasAbsolutePosition_ = true;
}

/**.......................................................................
 * Find the x-y shift, in integer pixels, between these two images
 */
void Image::addAbsolute(Image& image)
{
  if(image.hasData_) {
    ThrowColorError("Image contains no data");
  }

  if(xAxis_.getAngularResolution() != image.xAxis_.getAngularResolution()) {
    ThrowColorError("X angular resolutions of the images don't match", "red");
  }

  if(yAxis_.getAngularResolution() != image.yAxis_.getAngularResolution()) {
    ThrowColorError("Y angular resolutions of the images don't match", "red");
  }

  //------------------------------------------------------------
  // Get the separation, in pixels, between the centers of the two
  // images
  //------------------------------------------------------------

  Angle xSep, ySep;
  gcp::util::Astrometry::flatSkyApproximationSeparations(xSep, ySep, image.ra_, image.dec_, ra_, dec_);

  int xDeltaPix = (int)rint(xSep / xAxis_.getAngularResolution());
  int yDeltaPix = (int)rint(ySep / yAxis_.getAngularResolution());

  //------------------------------------------------------------
  // Now iterate over all pixels in the image being added, adding them
  // into the correct pixels in the target image
  //------------------------------------------------------------

  int targetXNpix = xAxis_.getNpix();
  int targetYNpix = yAxis_.getNpix();

  int xNpix = image.xAxis_.getNpix();
  int yNpix = image.yAxis_.getNpix();

  int iXTarget, iYTarget;

  // We should determine up-front the overlap of the two images and
  // only iterate over those pixels, but I don't have the energy to
  // figure this out right now, so just iterating over all pixels and
  // checking if they lie within the target image

  for(int iX=0; iX < xNpix; iX++) {

    iXTarget = targetXNpix/2 + xDeltaPix + (iX - xNpix/2);

    if(iXTarget >= 0 && iXTarget < targetXNpix) {

      for(int iY=0; iY < yNpix; iY++) {
	iYTarget = targetYNpix/2 + yDeltaPix + (iY - yNpix/2);

	if(iYTarget >= 0 && iYTarget < targetYNpix) {
	  unsigned ind       = iY * xNpix + iX;
	  unsigned indTarget = iYTarget * targetXNpix + iXTarget;

	  data_[indTarget] += image.data_[ind];
	}

      }

    }
  }

  hasData_ = true;
}

void Image::zeroBeyondRadius(Angle radius)
{
  unsigned nx = xAxis().getNpix();
  Angle xRes = xAxis().getAngularResolution();

  unsigned ny = yAxis().getNpix();
  Angle yRes = yAxis().getAngularResolution();

  for(unsigned ix=0; ix < nx; ix++) {
    double dx = ((double)(ix) - (double)(nx)/2) * xRes.radians();
    for(unsigned iy=0; iy < ny; iy++) {
      double dy = ((double)(iy) - (double)(ny)/2) * yRes.radians();

      double rad = sqrt(dx*dx + dy*dy);

      if(rad > radius.radians()) {
	unsigned ind = iy * nx + ix;
	data_[ind] = 0.0;
      }
    }
  }
}
