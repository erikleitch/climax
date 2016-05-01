#include "gcp/fftutil/Dft2d.h"
#include "gcp/fftutil/Image.h"
#include "gcp/fftutil/ObsInfo.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/util/Astrometry.h"
#include "gcp/util/Constants.h"
#include "gcp/util/FitsImageReader.h"
#include "gcp/util/FitsBinTableReader.h"
#include "gcp/util/FitsImageWriter.h"
#include "gcp/util/Stats.h"
#include "gcp/util/String.h"
#include "gcp/util/SzCalculator.h"
#include "gcp/util/Planck.h"

#include "gcp/util/PtSrcFinder.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace std;

using namespace gcp::util;

const Angle Image::zero_ = Angle(Angle::Degrees(), 0.0);

//=======================================================================
// An explanatory note about images: pixels are numbered from 0 -->
// N-1.  If dx represents the width of a pixel, then the axis value at
// the center of pixel ipix is given by (ipix + 0.5) * dx.  For an
// even number of pixels, the center of the image is at pixel value
// N/2 (2 in the example below).  For an odd number of pixels, the
// center of the image corresponds to N/2 + 0.5 (2.5 if N = 5, for
// example).  Ie, in pixel units:
// 
//                  4 _______________________
//                   |     |     |     |     |
//                   |     |     |     |     |
//                  3 -----------------------
//                   |     |     |     |     |
//                   |     |     |     |     |
//                  2 -----------------------
//                   |     |     |     |     |
//                   |     |     |     |     |
//                  1 -----------------------
//                   |     |     |     |     |
//                   | 0.5 | 1.5 | 2.5 | 3.5 |
//                    -----------------------
//                   0     1     2     3     4
//                   
//=======================================================================

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

void Image::assignDataFrom(const Image& image, unsigned iXStart, unsigned iYStart)
{
  assignDataFrom((Image&)image, iXStart, iYStart);
}

void Image::assignDataFrom(Image& image, unsigned iXStart, unsigned iYStart)
{
  unsigned thatNx = image.xAxis().getNpix();
  unsigned thatNy = image.yAxis().getNpix();

  unsigned thisNx = xAxis().getNpix();
  unsigned thisNy = yAxis().getNpix();

  unsigned lastXIndex = iXStart + thatNx - 1;
  unsigned lastYIndex = iYStart + thatNy - 1;

  if(lastXIndex > thisNx-1 || lastYIndex > thisNy-1) {
    ThrowColorError("Cannot copy data from image of size: " 
		    << thatNx << "x" << thatNy
		    << " into image of size "
		    << thisNx << "x" << thisNy
		    << " at starting index (" << iXStart << ", " << iYStart << ")",
		    "red");
  }
    
  for(unsigned ix=0; ix < thatNx; ix++) {
    for(unsigned iy=0; iy < thatNy; iy++) {

      unsigned iX = iXStart + ix;
      unsigned iY = iYStart + iy;
      unsigned thisInd = iY * thisNx + iX;
      unsigned thatInd = iy * thatNx + ix;

      data_[thisInd] = image.data_[thatInd];
    }
  }

  // Copy some of the ancillary data from the assigned image (but we
  // don't want to call copyAncillaryData() here because it will
  // overwrite, for example, the reference pixel, which can be
  // different between the two images

  hasData_  = true;

  units_    = image.units_;
  hasUnits_ = image.hasUnits_;
    
  frequency_    = image.frequency_;
  hasFrequency_ = image.hasFrequency_;
}

void Image::copyAncillaryData(Image& image)
{
  units_               = image.units_;
  hasUnits_            = image.hasUnits_;
  nativeToJy_          = image.nativeToJy_;

  hasData_             = image.hasData_;
  hasUnits_            = image.hasUnits_;
  hasNativeToJy_       = image.hasNativeToJy_;

  ra_                  = image.ra_;
  raRefPix_            = image.raRefPix_;
  dec_                 = image.dec_;
  decRefPix_           = image.decRefPix_;
  hasAbsolutePosition_ = image.hasAbsolutePosition_;

  dataMin_             = image.dataMin_;
  dataMax_             = image.dataMax_;

  frequency_           = image.frequency_;
  hasFrequency_        = image.hasFrequency_;
}

void Image::operator=(const Image& image) 
{
  *this = (Image&)image;
}

void Image::operator=(Image& image) 
{
  initialize();

  xAxis_.operator=((Axis&)image.xAxis_);
  yAxis_.operator=((Axis&)image.yAxis_);

  data_.resize(image.data_.size());
  wtSum_.resize(image.wtSum_.size());
  n_.resize(image.n_.size());
  valid_.resize(image.valid_.size());

  for(unsigned i=0; i < data_.size(); i++) {
    data_[i]  = image.data_[i];
    n_[i]     = image.n_[i];
    wtSum_[i] = image.wtSum_[i];
    valid_[i] = image.valid_[i];
  }
  
  copyAncillaryData(image);
}

void Image::initialize()
{
  isLatLng_ = false;
  data_.resize(0);
  n_.resize(0);
  wtSum_.resize(0);
  hasData_             = false;
  hasNativeToJy_       = false;
  raRefPix_            = 0.0;
  decRefPix_           = 0.0;

  dataMin_             = 0.0;
  dataMax_             = 0.0;

  axes_ = ImageAxis::AXIS_NONE;

  xAxis_.setAxisType(Axis::AXIS_X);
  yAxis_.setAxisType(Axis::AXIS_Y);

  xAxis_.setProjection("RA---SIN");
  yAxis_.setProjection("DEC--SIN");

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
  double radians = (1.02 * freq.centimeters() / diameter.centimeters()) / ::sqrt(8*log(2.0));
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
  n_.resize(data.size());
  wtSum_.resize(data.size());
  valid_.resize(data.size());

  for(unsigned i=0; i < data.size(); i++) {
    data_[i]  = data[i];
    n_[i]     = 1;
    wtSum_[i] = 1;
    valid_[i] = 1;
  }

  setUnits(units);
  hasData_ = true;
}

/**.......................................................................
 * Initialize the data in this image from an external data array
 */
void Image::initializeFromArray(std::valarray<double>& data, Unit::Units units)
{
  data_.resize(data.size());
  n_.resize(data.size());
  wtSum_.resize(data.size());
  valid_.resize(data.size());

  for(unsigned i=0; i < data.size(); i++) {
    data_[i]  = data[i];
    n_[i]     = 1;
    wtSum_[i] = 1;
    valid_[i] = 1;
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
  n_.resize(ndata);
  wtSum_.resize(ndata);

  for(unsigned i=0; i < ndata; i++) {
    data_[i]  = data[i];
    n_[i]     = 1;
    wtSum_[i] = 1;
    valid_[i] = 1;
  }

  setUnits(units);
  hasData_ = true;
}

/**.......................................................................
 * Initialize the data in this image from a FITS image file
 */
void Image::initializeFromFitsFile(std::string fileName, ObsInfo* obs)
{
  //------------------------------------------------------------
  // Determine some things about this file before trying to read it
  //------------------------------------------------------------

  FitsReader initReader;
  initReader.open(fileName);
  unsigned nAxis = initReader.getLongKey("NAXIS");
  initReader.close();

  if(nAxis == 0) {

    if(initReader.nHdu() > 1) {
      double raDegMean=0.0;
      double decDegMean=0.0;
      unsigned nImage=0;
      Image image;

      double gxmin, gxmax, gymin, gymax;

      for(unsigned iHdu=1; iHdu < initReader.nHdu(); iHdu++) {
	double xmin, xmax, ymin, ymax;

	image.initializeFromFitsFile(fileName, iHdu, obs);

	raDegMean  += (image.ra_.degrees()  -  raDegMean) / (nImage + 1);
	decDegMean += (image.dec_.degrees() - decDegMean) / (nImage + 1);
	
	xmin = image.ra_.degrees() - image.xAxis().getAngularSize().degrees() /2 ;
	xmax = image.ra_.degrees() + image.xAxis().getAngularSize().degrees() /2 ;
	swapIfNeeded(xmin, xmax);
	  
	ymin = image.dec_.degrees() - image.yAxis().getAngularSize().degrees() /2 ;
	ymax = image.dec_.degrees() + image.yAxis().getAngularSize().degrees() /2 ;
	swapIfNeeded(ymin, ymax);

	if(iHdu == 0) {
	  gxmin = xmin;
	  gxmax = xmax;
	  gymin = ymin;
	  gymax = ymax;
	}

	gxmin = gxmin < xmin ? gxmin : xmin;
	gxmax = gxmax > xmax ? gxmax : xmax;
	gymin = gymin < ymin ? gymin : ymin;
	gymax = gymax > ymax ? gymax : ymax;

	nImage++;

      }

      HourAngle    raMean;
      raMean.setDegrees(raDegMean);
      Declination decMean;
      decMean.setDegrees(decDegMean);

      Angle xSize(Angle::Degrees(), gxmax-gxmin);
      Angle ySize(Angle::Degrees(), gymax-gymin);

      unsigned nx = (ceil)(xSize/image.xAxis().getAngularResolution());
      unsigned ny = (ceil)(ySize/image.yAxis().getAngularResolution());

      xSize.setDegrees(nx * image.xAxis().getAngularResolution().degrees());
      ySize.setDegrees(ny * image.yAxis().getAngularResolution().degrees());

      Image composite;

      composite.xAxis().setAngularSize(xSize);
      composite.yAxis().setAngularSize(ySize);

      composite.xAxis().setNpix(nx);
      composite.yAxis().setNpix(ny);

      composite.setRaDec(raMean, decMean);

      for(unsigned iHdu=1; iHdu < initReader.nHdu(); iHdu++) {
	Image image;
	image.initializeFromFitsFile(fileName, iHdu);
	composite.addImage(image, OPER_AVERAGE);
	composite.display();
      }

      *this = composite;
    }

  } else {
    initializeFromFitsFile(fileName, 0, obs);
  }
}

void Image::initializeFromFitsFile(std::string fileName, unsigned iHdu, ObsInfo* obs)
{
  FitsImageReader reader(fileName, iHdu);

  FitsImageReader::Image image;

  reader.readData(image);

  FitsImageReader::Axis& xFitsAxis = reader.axes_[0];
  FitsImageReader::Axis& yFitsAxis = reader.axes_[1];

  resize(xFitsAxis.n_, yFitsAxis.n_);

  String xAxisType(xFitsAxis.type_);
  String xAxisTypeComment(xFitsAxis.typeComment_);

  if(xAxisType.contains("RA---TAN") || xAxisType.contains("RA---SIN") || xAxisType.contains("RA---ZPN") ||
     xAxisTypeComment.contains("degrees")) {

    xAxis().setAngularSize(Angle(Angle::Degrees(), xFitsAxis.n_ * xFitsAxis.delta_));
    xAxis().setProjection(xAxisType.str());
    ra_.setDegrees(xFitsAxis.refVal_);
    raRefPix_ = xFitsAxis.refPix_;
    hasAbsolutePosition_ = true;
  }

  String yAxisType(yFitsAxis.type_);
  String yAxisTypeComment(yFitsAxis.typeComment_);

  if(yAxisType.contains("DEC--TAN") || yAxisType.contains("DEC--SIN") || yAxisType.contains("DEC--ZPN") ||
     yAxisTypeComment.contains("degrees")) {

    yAxis().setAngularSize(Angle(Angle::Degrees(), yFitsAxis.n_ * yFitsAxis.delta_));
    yAxis().setProjection(yAxisType.str());
    dec_.setDegrees(yFitsAxis.refVal_);
    decRefPix_ = yFitsAxis.refPix_;

    hasAbsolutePosition_ = true;
  }

  for(unsigned i=0; i < data_.size(); i++) {
    data_[i] = image.data_[i];
  }

  hasData_ = true;

  units_ = Unit::stringToUnits(reader.units());

  if(units_ != Unit::UNITS_UNKNOWN) {
    hasUnits_ = true;
  }

  //------------------------------------------------------------
  // Initialize information about the observation
  //------------------------------------------------------------

  if(obs) {

    try {
      obs->setSourceName(    reader.object());
    } catch(Exception& err) {
    };

    try {
      obs->setObsRa(         reader.obsra());
    } catch(...) {};

    try {
      obs->setObsDec(        reader.obsdec());
    } catch(...) {};

    try {
      obs->setObsEquinox(    reader.equinox());
    } catch(...) {};

    try {
      obs->setTelescopeName( reader.telescope());
    } catch(...) {};

    try {
      obs->setInstrumentName(reader.instrument());
    } catch(...) {};

  }
}

void Image::initializeFromFitsTable(std::string fileName, std::string colName)
{
  gcp::util::FitsBinTableReader reader(fileName);

  reader.getNextTableInfo();

  std::vector<long> naxes = reader.getColNaxes(colName);
  
  if(naxes.size() != 2) 
    ThrowError("Column " << colName << " of file " << fileName << " does not appear to be image data");

  resize(naxes[0], naxes[1]);

  std::vector<double> data = reader.getDoubleData(colName);

  for(unsigned i=0; i < data.size(); i++)
    data_[i] = data[i];

  hasData_  = true;
}

void Image::initializeFromFitsTable(std::string fileName, std::string extName, std::string xColName, std::string yColName, std::string dataColName, ObsInfo* obs)
{
  gcp::util::FitsBinTableReader reader(fileName, extName);

  FitsReader::Axis xFitsAxis = reader.getAxis("RA--");
  FitsReader::Axis yFitsAxis = reader.getAxis("DEC--");

  std::vector<double> xVals = reader.getData(xColName);
  std::vector<double> yVals = reader.getData(yColName);
  std::vector<double> dVals = reader.getData(dataColName);

  // Convert axes to offsets relative to the fiducial position

  for(unsigned i=0; i < xVals.size(); i++) {
    xVals[i] = xFitsAxis.refVal_ + (xVals[i] - xFitsAxis.refPix_) * xFitsAxis.delta_;
    yVals[i] = yFitsAxis.refVal_ + (yVals[i] - yFitsAxis.refPix_) * yFitsAxis.delta_;
  }

  double xMin = Stats::min(xVals);
  double xMax = Stats::max(xVals);

  double yMin = Stats::min(yVals);
  double yMax = Stats::max(yVals);

  double xWidth = xMax - xMin;
  double yWidth = yMax - yMin;

  // Make the image large enough that the minimum value corresponds at
  // least to the center of the first pixel, and the max corresponds
  // at least to the center of the last pixel (i.e., +1)
  //
  // But only resize if this image doesn't already have a size

  unsigned nx, ny;

  try {
    nx = xAxis().getNpix();
    ny = yAxis().getNpix();
  } catch(...) {
    nx = ceil(::sqrt((double)xVals.size()) + 1);
    ny = ceil(::sqrt((double)yVals.size()) + 1);
    resize(nx, ny);
  }

  double dx = xWidth / nx;
  double dy = yWidth / nx;

  //------------------------------------------------------------
  // Now calculate where the reference value lies in our pixel space
  //
  //  i.e., xRef = (0.5 + iRefPix)*dx;
  //------------------------------------------------------------

  String xAxisType(xFitsAxis.type_);
  String xAxisTypeComment(xFitsAxis.typeComment_);

  if(xAxisType.contains("RA---TAN") || xAxisType.contains("RA---SIN") || xAxisType.contains("RA---ZPN")) {
    xAxis().setAngularSize(Angle(Angle::Degrees(), nx*dx));
    xAxis().setProjection(xAxisType.str());

    // We will set the center pixel as the reference

    double raRefPix  = fabs(xFitsAxis.refVal_ - xMin) / dx;
    raRefPix_ = nx/2;
    ra_.setDegrees(xFitsAxis.refVal_ + (raRefPix_ - raRefPix) * dx);

    hasAbsolutePosition_ = true;

    if(obs)
      obs->setObsRa(ra_);
  }

  String yAxisType(yFitsAxis.type_);
  String yAxisTypeComment(yFitsAxis.typeComment_);

  if(yAxisType.contains("DEC--TAN") || yAxisType.contains("DEC--SIN") || yAxisType.contains("DEC--ZPN")) {
    yAxis().setAngularSize(Angle(Angle::Degrees(), ny*dy));
    yAxis().setProjection(yAxisType.str());

    double decRefPix = fabs(yFitsAxis.refVal_ - yMin) / dy;
    decRefPix_ = ny/2;
    dec_.setDegrees(yFitsAxis.refVal_ + (decRefPix_ - decRefPix) * dy);

    hasAbsolutePosition_ = true;

    if(obs)
      obs->setObsDec(dec_);
  }

  // Finally, load the data from the vector of (possibly) irregularly sampled data

  Angle axisUnits;
  axisUnits.setUnits("degrees");
  axisUnits.setVal(1.0, "degrees");

  fillFrom(axisUnits, xVals, yVals, dVals);

  hasData_  = true;

  //------------------------------------------------------------
  // Initialize information about the observation
  //------------------------------------------------------------

  if(obs) {
    obs->setSourceName(    reader.getKeyword(extName, "OBJECT"));
    obs->setTelescopeName( reader.getKeyword("TELESCOP"));
    obs->setInstrumentName(reader.getKeyword("INSTRUME"));
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
    float crval = 0.0;

    if(hasAbsolutePosition_) {
      xtype = xAxis().getProjection();
      crval = ra_.degrees();
    }

    writer.putKey("CTYPE1", xtype, "Axis name (degrees)");
    float crpix = (float)(xAxis().getNpix())/2;
    writer.putKey("CRPIX1", crpix, "Reference pixel");

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

  //------------------------------------------------------------
  // Now the Y-axis
  //------------------------------------------------------------

  if(yAxis().hasAngularSize()) {

    float crval = 0.0;
    std::string ytype("Y");
    if(hasAbsolutePosition_) {
      ytype = yAxis().getProjection();
      crval = dec_.degrees();
    }

    writer.putKey("CTYPE2", ytype, "Axis name (degrees)");
    float crpix = (float)(yAxis().getNpix())/2;
    writer.putKey("CRPIX2", crpix, "Reference pixel");

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

  //------------------------------------------------------------
  // Write the units, if they are known
  //------------------------------------------------------------

  writer.putKey("BUNIT", unitToString(), "Unit of measurement");

  //------------------------------------------------------------
  // Write the units, if they are known
  //------------------------------------------------------------

  if(hasFrequency_)
    writer.putKey("RESTFREQ", frequency_.Hz(), "Frequency in Hz");

  //------------------------------------------------------------
  // Now write the data
  //------------------------------------------------------------

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

  unsigned npix = (unsigned)::sqrt((double)(buf.st_size/sizeof(float)));

  resize(npix, npix);

  read(fd, (void*)&data_[0], npix*npix*sizeof(float));

  hasData_ = true;

  close(fd);
}

void Image::resize()
{
  data_.resize(xAxis().getNpix() * yAxis().getNpix());
  n_.resize(xAxis().getNpix() * yAxis().getNpix());
  wtSum_.resize(xAxis().getNpix() * yAxis().getNpix());
  valid_.resize(xAxis().getNpix() * yAxis().getNpix());

  valid_ = 1;

  // Initialize these to the midpoint of the array

  raRefPix_  = double(xAxis().getNpix())/2 - 0.5;
  decRefPix_ = double(yAxis().getNpix())/2 - 0.5;
}

void Image::initializeRefPixFftConvention()
{
  raRefPix_  = double(xAxis().getNpix())/2;
  decRefPix_ = double(yAxis().getNpix())/2;
}

void Image::initializeRefPixNonFftConvention()
{
  raRefPix_  = double(xAxis().getNpix())/2 - 0.5;
  decRefPix_ = double(yAxis().getNpix())/2 - 0.5;
}

void Image::resize(unsigned nx, unsigned ny)
{
  axes_ = ImageAxis::AXIS_NONE;
  xAxis().setNpix(nx);
  yAxis().setNpix(ny);

  // Initialize these to the (non FFT convention) midpoint of the array

  raRefPix_  = double(nx)/2 - 0.5;
  decRefPix_ = double(ny)/2 - 0.5;
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
 * Initialize this object
 */
void Image::initialize(Angle& xSize, Angle& ySize, unsigned nXPix, unsigned nYPix)
{
  initialize();
  xAxis().setNpix(nXPix);
  yAxis().setNpix(nYPix);
  xAxis().setAngularSize(xSize);
  yAxis().setAngularSize(ySize);
}

/**.......................................................................
 * Determine the conversion from native image units to Jy
 */
void Image::setNativeToJy(Frequency& nu, SolidAngle& beam)
{
  setNativeToJy(nativeToJy(nu, beam));
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
 * Display this image
 */
void Image::display(bool setWin)
{
  PgUtil::setWin(setWin);
  display(true, true);
}

void Image::display(bool normalDisplay, bool greyScale)
{
  if(!hasData_)
    ThrowError("This image contains no data to display");

  PgUtil::setPrompt(false);
  PgUtil::setWnad(true);

  if(!xAxis().hasNpix()) 
    ThrowError("x-axis length hasn't been specified");

  if(!yAxis().hasNpix())
    ThrowError("y-axis length hasn't been specified");

  std::string unit = unitToString();
  std::string xlab = xAxis().name();
  std::string ylab = yAxis().name();

  double xmin, xmax, ymin, ymax;

  xmin =  xMin();
  xmax =  xMax();
  ymin =  yMin();
  ymax =  yMax();

  if(!normalDisplay) {
    PgUtil::setReverseX(true);
  }

  xmin *= xAxis().scale();
  xmax *= xAxis().scale();
  ymin *= yAxis().scale();
  ymax *= yAxis().scale();

  //------------------------------------------------------------
  // Set appropriate callbacks
  //------------------------------------------------------------

  if(hasAbsolutePosition_) {
    if(isLatLng_)
      PgUtil::setCoordCallback(Image::latLongCallback, this);
    else
      PgUtil::setCoordCallback(Image::angularCallback, this);
  }

  PgUtil::setUnitCallback(Image::unitCallback, this);

  if(normalDisplay)
    PgUtil::setReverseX(false);
  else
    PgUtil::setReverseX(true);

  if(greyScale) {
    PgUtil::greyScale(data_.size(), &data_[0], xAxis().getNpix(), yAxis().getNpix(),
		      xmin, xmax, ymin, ymax,
		      0, 0, 0,
		      (char*)xlab.c_str(), (char*)ylab.c_str(), "", (char*)unit.c_str());
  } else {
    PgUtil::contour(data_.size(), &data_[0], xAxis().getNpix(), yAxis().getNpix(),
		    xmin, xmax, ymin, ymax,
		    0, 0, 0, 10, 
		    (char*)xlab.c_str(), (char*)ylab.c_str(), "", (char*)unit.c_str());
  }

#if 0
  PtSrcFinder finder;
  std::vector<HourAngle> ras;
  std::vector<Declination> decs;
  COUT("Looking for NVSS sources");
  Angle rad;
  rad.setDegrees("00:10:00");
  finder.findSources("/Users/eml/projects/climax/climaxSrc/", "nvss", ra_, dec_, rad, ras, decs);

  for(unsigned iSrc=0; iSrc < ras.size(); iSrc++) {
    double xoff = ras[iSrc].degrees()  - ra_.degrees();
    double yoff = decs[iSrc].degrees() - dec_.degrees();
    PgUtil::plotPoint(xoff, yoff, 1);
  }

  COUT("Looking for FIRST sources");
  finder.findSources("/Users/eml/projects/climax/climaxSrc/", "first", ra_, dec_, rad, ras, decs);

  for(unsigned iSrc=0; iSrc < ras.size(); iSrc++) {
    double xoff = ras[iSrc].degrees()  - ra_.degrees();
    double yoff = decs[iSrc].degrees() - dec_.degrees();
    PgUtil::plotPoint(xoff, yoff, 2);
  }
  COUT("Done");
#endif

  if(hasAbsolutePosition_)
    PgUtil::setCoordCallback(0, 0);

  PgUtil::setUnitCallback(0, 0);
}

void Image::logDisplay()
{
  logDisplay(true,true);
}

void Image::logDisplay(bool normalDisplay, bool greyScale)
{
  if(!hasData_)
    ThrowError("This image contains no data to display");

  PgUtil::setPrompt(false);

  if(!xAxis().hasNpix()) 
    ThrowError("x-axis length hasn't been specified");

  if(!yAxis().hasNpix())
    ThrowError("y-axis length hasn't been specified");

  std::string unit = unitToString();
  std::string xlab = xLab();
  std::string ylab = yLab();

  double xmin, xmax, ymin, ymax;

  if(normalDisplay) {
    xmin =  xMin();
    xmax =  xMax();
    ymin =  yMin();
    ymax =  yMax();
  } else {
    xmin = -xMin();
    xmax = -xMax();
    ymin =  yMin();
    ymax =  yMax();
  }

  if(hasAbsolutePosition_)
    PgUtil::setCoordCallback(Image::angularCallback, this);

  if(normalDisplay)
    PgUtil::setReverseX(false);
  else
    PgUtil::setReverseX(true);

  std::valarray<float> pldata = data_;
  for(unsigned i=0; i < data_.size(); i++) {
    if(pldata[i] > 0.0) {
      pldata[i] = pow((double)pldata[i], (double)0.3);
    } else {
      pldata[i] = 0.0;
    }
  }

  if(greyScale) {
    PgUtil::greyScale(data_.size(), &pldata[0], xAxis().getNpix(), yAxis().getNpix(),
		      xmin, xmax, ymin, ymax,
		      0, 0, 0,
		      (char*)xlab.c_str(), (char*)ylab.c_str(), "", (char*)unit.c_str());
  } else {
    PgUtil::contour(data_.size(), &pldata[0], xAxis().getNpix(), yAxis().getNpix(),
		    xmin, xmax, ymin, ymax,
		    0, 0, 0, 10, 
		    (char*)xlab.c_str(), (char*)ylab.c_str(), "", (char*)unit.c_str());
  }

  if(hasAbsolutePosition_)
    PgUtil::setCoordCallback(0, 0);
}

/**.......................................................................
 * Display this image the way difmap displays things
 */
void Image::difmapDisplay()
{
  display(false, true);
}

void Image::contour()
{
  PgUtil::setPrompt(false);

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

  //------------------------------------------------------------
  // Set appropriate callbacks
  //------------------------------------------------------------

  if(hasAbsolutePosition_)
    PgUtil::setCoordCallback(Image::angularCallback, this);

  PgUtil::setUnitCallback(Image::unitCallback, this);
    
  //------------------------------------------------------------
  // Make contour plot
  //------------------------------------------------------------

  PgUtil::contour(data_.size(), &data_[0], xAxis().getNpix(), yAxis().getNpix(),
		  xmin, xmax, ymin, ymax,
		  0, 0, 0, 10, 
		  xlab, ylab, "", unit);

  //------------------------------------------------------------
  // Unset appropriate callbacks
  //------------------------------------------------------------

  if(hasAbsolutePosition_)
    PgUtil::setCoordCallback(0, 0);

  PgUtil::setUnitCallback(0, 0);
}

/**.......................................................................
 * Return the rms of the image data within an annulus
 */
double Image::rms(Angle& rmin, Angle& rmax)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  double imMean = mean(rmin, rmax);
  double imSdev = 0.0;
  double val;

  unsigned nx  = xAxis().getNpix();
  unsigned ny  = yAxis().getNpix();

  double dxRad = xAxis().getAngularResolution().radians();
  double dyRad = yAxis().getAngularResolution().radians();

  double rminRad = rmin.radians();
  double rmaxRad = rmax.radians();
  
  double x, y, rad;
  unsigned npt=0;
  for(int iy=0; iy < ny; iy++) {
    for(int ix=0; ix < nx; ix++) {

      x = (double)(ix) - (double)(nx/2);
      y = (double)(iy) - (double)(ny/2);

      x *= dxRad;
      y *= dyRad;

      rad = ::sqrt(x*x + y*y);

      if(rad >= rminRad && rad < rmaxRad) {
	val = data_[iy * nx + ix] - imMean;
	imSdev += (val * val - imSdev) / (npt + 1);
	npt++;
      }
    }
  }

  imSdev = (npt > 1) ? ::sqrt(npt * imSdev / (npt - 1)) : 0.0;

  if(npt == 0) {
    ThrowError("No pixels lie within " << rmin << " <= r < " << rmax);
  }

  return imSdev;
}

/**.......................................................................
 * Return the mean of the image data within an annulus
 */
double Image::mean(Angle& rmin, Angle& rmax)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  double imMean = 0.0;
  double val;

  unsigned nx  = xAxis().getNpix();
  unsigned ny  = yAxis().getNpix();

  double dxRad = xAxis().getAngularResolution().radians();
  double dyRad = yAxis().getAngularResolution().radians();

  double rminRad = rmin.radians();
  double rmaxRad = rmax.radians();

  double x, y, rad;
  unsigned npt=0;
  for(int iy=0; iy < ny; iy++) {
    for(int ix=0; ix < nx; ix++) {

      x = (double)(ix) - (double)(nx/2);
      y = (double)(iy) - (double)(ny/2);

      x *= dxRad;
      y *= dyRad;

      rad = ::sqrt(x*x + y*y);

      if(rad >= rminRad && rad < rmaxRad) {
	imMean += (data_[iy * nx + ix] - imMean) / (npt + 1);
	npt++;
      }
    }
  }

  if(npt == 0) {
    ThrowError("No pixels lie within " << rmin << " <= r < " << rmax);
  }

  return imMean;
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

  imSdev = (n > 1) ? ::sqrt(n * imSdev / (n - 1)) : 0.0;

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

  double imMax = data_[0];
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

  double imMin = data_[0];
  double val;
  for(unsigned i=0; i < data_.size(); i++) {
    val = data_[i];
    imMin = (val < imMin) ? val : imMin;
  }

  return imMin;
}

/**.......................................................................
 * Return the min of the image data
 */
double Image::minGreaterThanZero()
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  bool first = true;
  double imMin;
  double val;
  for(unsigned i=0; i < data_.size(); i++) {
    val = data_[i];
    if(val > 0.0) {
      if(first) {
	imMin = val;
	first = false;
      } else {
	imMin = (val < imMin) ? val : imMin;
      }
    }
  }

  if(first)
    ThrowError("No pixel greater than zero found");

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
      addImage(image, OPER_ADD);
    } else {
      addWithValidation(image);
    }

    hasData_ = true;
  }
}

void Image::operator-=(Image& image)
{
  if(!hasData_ || !image.hasData_) {
    ThrowError("Image contains no data");
  }

  for(unsigned i=0; i < data_.size(); i++) {
    
    // Only divide this point if the pixel being added is valid
    
    if(image.valid_[i]) {
      
      // If the pixel we are adding to is also valid, then coadd them
      
      if(valid_[i]) {
	data_[i] -= image.data_[i];
      }

    } else {
      valid_[i] = 0;
    }
  }
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

  for(unsigned i=0; i < data_.size(); i++) {
    
    // Only divide this point if the pixel being added is valid
    
    if(image.valid_[i]) {
      
      // If the pixel we are adding to is also valid, then coadd them
      
      if(valid_[i]) {
	data_[i] /= image.data_[i];
      }

    } else {
      valid_[i] = 0;
    }
  }

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
void Image::convolve(Image& image, bool zeropad, Dft2d* dft)
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

  dft1.setPlan(dft);
  dft2.setPlan(dft);
  
  dft1.initialize(*this);
  dft2.initialize(image);

  dft1.normalize(true);

  dft1.computeForwardTransform();
  dft2.computeForwardTransform();

  dft1.complexMultiply(dft2, false);
  dft1.shift();

  dft1.computeInverseTransform();

  // Just copy the data, rather than the whole returned image,
  // otherwise we'll lose ancillary data in this image

  Image ret = dft1.getImage(false);
  data_ = ret.data_;
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
  case Unit::UNITS_JYSR:
    return "Jy/sr";
    break;
  case Unit::UNITS_JYBEAM:
    return "Jy/beam";
    break;
  case Unit::UNITS_MEGAJYSR:
    return "MJy/sr";
    break;
  case Unit::UNITS_MILLIJYSR:
    return "mJy/sr";
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
  case Unit::UNITS_SNR:
    return "SNR";
    break;
  case Unit::UNITS_COUNTS:
    return "COUNTS";
    break;
  case Unit::UNITS_COUNT_FLUX:
    return "COUNTS/cm\\u2\\d/s";
    break;
  case Unit::UNITS_METERS:
    return "METERS";
    break;
  case Unit::UNITS_FEET:
    return "FEET";
    break;
  case Unit::UNITS_KFEET:
    return "1000 FEET";
    break;
  case Unit::UNITS_NONE:
    return " ";
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

void Image::Axis::operator=(const Axis& axis)
{
  return operator=((ImageAxis&) axis);
}

void Image::Axis::operator=(Axis& axis)
{
  return operator=((ImageAxis&) axis);
}

void Image::Axis::operator=(const ImageAxis& axis)
{
  return operator=((ImageAxis&) axis);
}

void Image::Axis::operator=(ImageAxis& axis)
{
  if(axis.hasAngularSize()) {
    setAngularSize(axis.getAngularSize());
    sense_ = axis.getSense();
  }

  if(axis.hasNpix())
    setNpix(axis.getNpix());

  projection_ = axis.projection_;
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
  data_  = 0.0;
  wtSum_ = 0.0;
  n_     = 0;
}

void Image::setValue(float val)
{
  data_  = val;
  wtSum_ = 1;
  n_     = 1;
}

float& Image::val(unsigned ix, unsigned iy)
{
  if(ix > xAxis_.getNpix()-1 || iy > yAxis_.getNpix()-1) {
    ThrowError("Invalid pixel index: " << ix << "," << iy);
  }

  unsigned imInd = iy * xAxis_.getNpix() + ix;

  return data_[imInd];
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

  n_ = 1;
  hasData_ = true;
}

/**.......................................................................
 * Interpolate this image at the requested offset from the reference pixel
 */
void Image::interpolateData(Angle& xOff, Angle& yOff, double& val, bool& valid, bool doThrow)
{
  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();

  int xSense = xAxis().getSense();
  int ySense = yAxis().getSense();

  // First get the nearest pixel to the requested point

  int ixNear, iyNear;
  getPixelRef(xOff, yOff, ixNear, iyNear, false, doThrow);

  int      convMaskInPixels = Dft2d::convMaskInPixels_;
  double   convSigInPixels  = Dft2d::convSigInPixels_;

  val = 0.0;

  int xdelta = (int)nx - 1 - ixNear;
  int ydelta = (int)ny - 1 - iyNear;

  if(ixNear < convMaskInPixels || xdelta < convMaskInPixels) {
    val = 0.0;
    valid = false;
    return;
  }

  if(iyNear < convMaskInPixels || ydelta < convMaskInPixels) {
    val = 0.0;
    valid = false;
    return;
  }

  int ixStart = ixNear - convMaskInPixels;
  int ixStop  = ixNear + convMaskInPixels;

  int iyStart = iyNear - convMaskInPixels;
  int iyStop  = iyNear + convMaskInPixels;

  double dxRad = xAxis().getAngularResolution().radians();
  double dyRad = yAxis().getAngularResolution().radians();

  // Get the requested position in pixel coordinates

  double xReqPix = raRefPix_  + xSense * xOff.radians() / dxRad;
  double yReqPix = decRefPix_ + ySense * yOff.radians() / dyRad;

  // Finally, initialize sums

  double wtSum = 0.0;
  double s2    = 2 * convSigInPixels * convSigInPixels;

  // Finally, accumulate the running mean over the pixel mask

  for(int ix = ixStart; ix <= ixStop; ix++) {

    for(int iy = iyStart; iy <= iyStop; iy++) {

      unsigned imInd = iy * nx + ix;

      // Calculate the coordinate of the center of this pixel

      double xPix = (double)(ix);
      double yPix = (double)(iy);
      
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
 * Interpolate this image at the requested offset from the reference pixel
 */
void Image::getNearestData(Angle& xOff, Angle& yOff, double& val, bool& valid, bool doThrow)
{
  unsigned nx = xAxis_.n_;
  unsigned ny = yAxis_.n_;

  int xSense = xAxis_.sense_;
  int ySense = yAxis_.sense_;

  // First get the nearest pixel to the requested point

  int ixNear, iyNear;
  getPixelRef(xOff, yOff, ixNear, iyNear, false, doThrow);

  if(ixNear < 0 || ixNear > nx-1 || iyNear < 0 || iyNear > ny-1) {
    val = 0.0;
    valid = false;
  } else {
    val = data_[iyNear * nx + ixNear];
    valid = true;
  }

  return;
}

/**.......................................................................
 * Interpolate this image at the requested off from the reference pixel
 */
void Image::fillFrom(Angle& axisUnits, std::vector<double>& x, std::vector<double>& y, std::vector<double>& d)
{
  std::valarray<double> xvals(x.size());
  std::valarray<double> yvals(y.size());
  std::valarray<double> dvals(d.size());

  for(unsigned i=0; i < x.size(); i++) {
    xvals[i] = x[i];
    yvals[i] = y[i];
    dvals[i] = d[i];
  }
    
  fillFrom(axisUnits, xvals, yvals, dvals);
}

/**.......................................................................
 * Fill this image from an array of irregularly sampled data
 */
void Image::fillFrom(Angle& axisUnits, std::valarray<double>& x, std::valarray<double>& y, std::valarray<double>& d)
{
  // Explicitly zero the data array

  zero();

  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();

  int xSense = xAxis().getSense();
  int ySense = yAxis().getSense();

  //------------------------------------------------------------
  // Now iterate over all elements in the array
  //------------------------------------------------------------

  Angle xOff, yOff;
  Angle xRef, yRef;

  xRef.setDegrees(ra_.degrees());
  yRef.setDegrees(dec_.degrees());

  double dxRad = xAxis().getAngularResolution().radians();
  double dyRad = yAxis().getAngularResolution().radians();

  for(unsigned i=0; i < x.size(); i++) {

    // First get the nearest pixel to the requested point

    xOff.setVal(x[i], axisUnits.units());
    yOff.setVal(y[i], axisUnits.units());

    int ixNear = raRefPix_  + xSense * (floor)((xOff.radians() - xRef.radians()) / dxRad);
    int iyNear = decRefPix_ + ySense * (floor)((yOff.radians() - yRef.radians()) / dyRad);

    int      convMaskInPixels = Dft2d::convMaskInPixels_;
    double   convSigInPixels  = Dft2d::convSigInPixels_;

    int xdelta = (int)nx - 1 - ixNear;
    int ydelta = (int)ny - 1 - iyNear;

    int ixStart = ixNear - convMaskInPixels;
    int ixStop  = ixNear + convMaskInPixels;
    
    int iyStart = iyNear - convMaskInPixels;
    int iyStop  = iyNear + convMaskInPixels;
    
    // Get the requested position in pixel coordinates
    
    double xReqPix = raRefPix_  + xSense * (xOff.radians() - xRef.radians()) / dxRad;
    double yReqPix = decRefPix_ + ySense * (yOff.radians() - yRef.radians()) / dyRad;
    
    // Finally, initialize sums
    
    double wtSum = 0.0;
    double s2    = 2 * convSigInPixels * convSigInPixels;
    
    // Finally, accumulate the running mean over the pixel mask
    
    for(int ix = ixStart; ix <= ixStop; ix++) {

      if(ix < 0 || ix > nx-1)
	continue;

      for(int iy = iyStart; iy <= iyStop; iy++) {
	
	if(iy < 0 || iy > ny-1)
	  continue;
	
	unsigned imInd = iy * nx + ix;
	
	// Calculate the coordinate of the center of this pixel
	
	double xPix = (double)(ix);
	double yPix = (double)(iy);
	
	// Calculate the delta of the requested point from the pixel
	// center
	
	double dxPix = xReqPix - xPix;
	double dyPix = yReqPix - yPix;
	
	double arg = (dxPix * dxPix + dyPix * dyPix)/s2;
	
	double wt = exp(-arg);
	
	// Accumulate running mean
	
	data_[imInd] += ((d[i] - data_[imInd]) * wt) / (wtSum_[imInd] + wt);
	wtSum_[imInd] += wt;
	n_[imInd] += 1;
      }
    }
  }

  hasData_ = true;

  return;
}

/**.......................................................................
 * Bin this image at the requested point.  We are summing all pixels
 * in the image whose centers lie within xRes/yRes of the requested
 * offset position.
 *
 *
 */
void Image::binData(Angle& xOff, Angle& xRes, Angle& yOff, Angle& yRes, 
		    double& val, bool& valid)
{
  //------------------------------------------------------------
  // Store the number of pixels in each axis
  //------------------------------------------------------------

  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();

  //------------------------------------------------------------
  // Now get the nearest pixel to the requested point
  //------------------------------------------------------------

  double xResRad = xAxis().getAngularResolution().radians();
  double yResRad = yAxis().getAngularResolution().radians();

  int xSense = xAxis().getSense();
  int ySense = yAxis().getSense();

  int ixNear = raRefPix_  + xSense * (floor)(xOff / xAxis().getAngularResolution());
  int iyNear = decRefPix_ + ySense * (floor)(yOff / yAxis().getAngularResolution());

  //------------------------------------------------------------
  // Now compute the range of pixels over which we will iterate.  
  // This will be the number of pixels that fall within the boundary:
  // 
  // (xOff +- xRes/2, yOff +- yRes/2), plus 1
  //
  //------------------------------------------------------------

  int xMaskInPixels = (xRes.radians() / xResRad)/2 + 1;
  int yMaskInPixels = (yRes.radians() / yResRad)/2 + 1;

  // Adjust indices so we don't fall off the ends of the array

  unsigned ixStart = ixNear < xMaskInPixels ? 0 : ixNear - xMaskInPixels;
  unsigned ixStop  = (ixNear + xMaskInPixels) > (nx-1) ? nx-1 : ixNear + xMaskInPixels;

  unsigned iyStart = iyNear < yMaskInPixels ? 0 : iyNear - yMaskInPixels;
  unsigned iyStop  = (iyNear + yMaskInPixels) > (ny-1) ? ny-1 : iyNear + yMaskInPixels;

  // Now accumulate the sum over all pixels in the pixel mask that
  // fall within the bounds

  val = 0.0;

  // Find the boundaries of the pixel for which we want to compute the
  // binned value

  double xOffMin = xOff.radians() - 0.5 * xRes.radians();
  double xOffMax = xOff.radians() + 0.5 * xRes.radians();

  double yOffMin = yOff.radians() - 0.5 * yRes.radians();
  double yOffMax = yOff.radians() + 0.5 * yRes.radians();

  swapIfNeeded(xOffMin, xOffMax);
  swapIfNeeded(yOffMin, yOffMax);

  for(unsigned ix = ixStart; ix <= ixStop; ix++) {
    for(unsigned iy = iyStart; iy <= iyStop; iy++) {

      unsigned imInd = iy * nx + ix;

      // Calculate the offset of the center of this pixel from the
      // reference pixel

      double xOffCenter = xSense * ((double)(ix) - raRefPix_ ) * xResRad;
      double yOffCenter = ySense * ((double)(iy) - decRefPix_) * yResRad;

      // If this pixel's center lies within the boundary of the region
      // over which we're summing, include it in the sum

      if(xOffCenter >= xOffMin && xOffCenter < xOffMax &&
	 yOffCenter >= yOffMin && yOffCenter < yOffMax) {

	// Accumulate sum

	double pixVal = data_[imInd];
	val += pixVal;
	valid = true;
      }
    }
  }

  return;
}

/**.......................................................................
 * Method to set an absolute position for this image
 */
void Image::setLatLng(Angle lat, Angle lng)
{
  HourAngle ra;
  Declination dec;

  ra.setDegrees(lng.degrees());
  dec.setDegrees(lat.degrees());

  isLatLng_ = true;
  setRaDec(ra, dec);
}

/**.......................................................................
 * Method to set an absolute position for this image
 */
void Image::setRaDec(HourAngle ra, Declination dec)
{
  int nx = xAxis().getNpix();
  int ny = yAxis().getNpix();

  double raRefPix, decRefPix;

  // We default to non-FFT convention here.  If an image is n pixels on a
  // side, and pixel indices run from 0 to n-1, then the center of the
  // image is the n/2 - 0.5'th pixel.
  //
  // FFT convention is that the center of an n x n image is the center
  // of the n/2'th pixel, or n/2 in pixel units

  raRefPix  = (double)(nx)/2 - 0.5;
  decRefPix = (double)(ny)/2 - 0.5;

  setRaDec(ra, dec, (double)raRefPix, (double)decRefPix);
}

/**.......................................................................
 * Method to set an absolute position for this image, FFT convention
 */
void Image::setRaDecFft(HourAngle ra, Declination dec)
{
  int nx = xAxis().getNpix();
  int ny = yAxis().getNpix();

  double raRefPix, decRefPix;

  // FFT convention is that the center of an n x n image is the center
  // of the n/2'th pixel, or n/2 in pixel units

  raRefPix  = (double)(nx)/2;
  decRefPix = (double)(ny)/2;

  setRaDec(ra, dec, (double)raRefPix, (double)decRefPix);
}

void Image::setRaDec(HourAngle ra, Declination dec, double raRefPix, double decRefPix)
{
  ra_        = ra;
  dec_       = dec;
  raRefPix_  = raRefPix;
  decRefPix_ = decRefPix;

  hasAbsolutePosition_ = true;
}

/**.......................................................................
 * Use the validity flags when coadding pixel values
 */
void Image::addWithValidation(Image& image)
{
  for(unsigned i=0; i < data_.size(); i++)  {

    //------------------------------------------------------------   
    // Only add this point if the pixel being added is valid
    //------------------------------------------------------------   
    
    if(image.valid_[i]) {
      
      //------------------------------------------------------------   
      // If the pixel we are adding to is also valid, then coadd them
      //------------------------------------------------------------   

      if(valid_[i])
	data_[i] += image.data_[i];

      //------------------------------------------------------------   
      // Else assign this point if it's the first valid data in this pixel
      //------------------------------------------------------------   

      else {
	data_[i] = image.data_[i];
	valid_[i] = 1;
      }
    }
  }
}

/**.......................................................................
 * Use validity flags when assigning pixel values
 */
void Image::assignWithValidation(Image& image)
{
  for(unsigned i=0; i < data_.size(); i++) 
    if(image.valid_[i])
      data_[i] = image.data_[i];
}

/**.......................................................................
 * Add an image to our image, using the absolute positions to register
 * where in our image the passed image belongs
 */
void Image::addImage(Image& image, Operation oper)
{
  if(!image.hasData_) {
    ThrowColorError("Image contains no data", "red");
  }

  if(xAxis_.getAngularResolution() != image.xAxis_.getAngularResolution()) {
    ThrowColorError("X angular resolutions of the images don't match (this = " 
		    << xAxis_.getAngularResolution() << " that = " 
		    << image.xAxis_.getAngularResolution() << ")", "red");
  }

  if(yAxis_.getAngularResolution() != image.yAxis_.getAngularResolution()) {
    ThrowColorError("Y angular resolutions of the images don't match (this = " 
		    << yAxis_.getAngularResolution() << " that = " 
		    << image.yAxis_.getAngularResolution() << ")", "red");
  }

  //------------------------------------------------------------
  // Get the separation, in pixels, between the reference positions of
  // the two images
  //------------------------------------------------------------

  Angle xSep, ySep;
  getSeparationOfImageFromUs(xSep, ySep, image);

  //------------------------------------------------------------
  // Convert this absolute separation to a pixel delta, relative to
  // our reference pixel
  //------------------------------------------------------------

  int xDeltaRefPix = (int)rint(xSep / xAxis_.getAngularResolution() * xAxis_.getSense());
  int yDeltaRefPix = (int)rint(ySep / yAxis_.getAngularResolution() * yAxis_.getSense());

  //------------------------------------------------------------
  // A pixel in the passed image that is a positive offset relative to
  // its reference pixel will be a positive offset in the target image
  // if both images have the same sense, but a negative offset if the
  // images have opposite sense.
  //------------------------------------------------------------

  int xSense = xAxis_.getSense() * image.xAxis().getSense();
  int ySense = yAxis_.getSense() * image.yAxis().getSense();

  //------------------------------------------------------------
  // Now iterate over all pixels in the image being added, adding them
  // into the correct pixels in the target image
  //------------------------------------------------------------

  int targetXNpix   = xAxis_.getNpix();
  int targetYNpix   = yAxis_.getNpix();

  int srcXNpix      = image.xAxis_.getNpix();
  int srcYNpix      = image.yAxis_.getNpix();

  int targetXRefPix = raRefPix_;
  int targetYRefPix = decRefPix_;

  int srcXRefPix    = image.raRefPix_;
  int srcYRefPix    = image.decRefPix_;

  //------------------------------------------------------------
  // We should determine up-front the overlap of the two images and
  // only iterate over those pixels, but I don't have the energy to
  // figure this out right now, so just iterating over all pixels and
  // checking if they lie within the target image
  //------------------------------------------------------------

  int iXTarget, iYTarget;

  for(int iXSrc=0; iXSrc < srcXNpix; iXSrc++) {

    //------------------------------------------------------------
    // Get the x-pixel in the target image corresponding to the pixel in
    // the image to be added
    //------------------------------------------------------------

    iXTarget = targetXRefPix + (iXSrc - srcXRefPix) * xSense + xDeltaRefPix;

    if(iXTarget >= 0 && iXTarget < targetXNpix) {

      for(int iYSrc=0; iYSrc < srcYNpix; iYSrc++) {

	//------------------------------------------------------------
	// Get the y-pixel in the target image corresponding to the pixel in
	// the image to be added
	//------------------------------------------------------------

	iYTarget = targetYRefPix + (iYSrc - srcYRefPix) * ySense + yDeltaRefPix;

	if(iYTarget >= 0 && iYTarget < targetYNpix) {

	  unsigned indSrc    =    iYSrc *    srcXNpix +    iXSrc;
	  unsigned indTarget = iYTarget * targetXNpix + iXTarget;

	  switch (oper) {
	  case OPER_ADD:
	    if(image.valid_[indSrc]) {
	      if(valid_[indTarget]) {
		data_[indTarget] += image.data_[indSrc];
	      } else {
		data_[indTarget] = image.data_[indSrc];
		valid_[indTarget] = 1;
	      }
	    }
	    break;
	  case OPER_AVERAGE:
	    if(image.valid_[indSrc]) {
	      data_[indTarget] += (image.data_[indSrc] - data_[indTarget])/(n_[indTarget]+1);
	      n_[indTarget] += 1;
	      valid_[indTarget] = 1;
	    }
	    break;
	  case OPER_ASSIGN:
	    if(image.valid_[indSrc]) {
	      data_[indTarget] = image.data_[indSrc];
	      valid_[indTarget] = 1;
	    }
	    break;
	  default:
	    ThrowError("Unrecognized operation: " << oper);
	    break;
	  }
	    
	}

      }

    }
  }

  hasData_ = true;
}

/**.......................................................................
 * Fill this image with zeros beyond the specified radius
 */
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

      double rad = ::sqrt(dx*dx + dy*dy);

      if(rad > radius.radians()) {
	unsigned ind = iy * nx + ix;
	data_[ind] = 0.0;
      }
    }
  }
}

/**.......................................................................
 * Fill this image with NANs beyond the specified radius
 */
void Image::invalidate()
{
  for(unsigned i=0; i < valid_.size(); i++)
    valid_[i] = 0;
}

/**.......................................................................
 * Fill this image with NANs beyond the specified radius
 */
void Image::invalidateBeyondRadius(Angle radius)
{
  unsigned nx = xAxis().getNpix();
  Angle xRes = xAxis().getAngularResolution();

  unsigned ny = yAxis().getNpix();
  Angle yRes = yAxis().getAngularResolution();

  for(unsigned ix=0; ix < nx; ix++) {
    double dx = ((double)(ix) - (double)(nx)/2) * xRes.radians();
    for(unsigned iy=0; iy < ny; iy++) {
      double dy = ((double)(iy) - (double)(ny)/2) * yRes.radians();

      double rad = ::sqrt(dx*dx + dy*dy);

      if(rad > radius.radians()) {
	unsigned ind = iy * nx + ix;
	valid_[ind] = 0;
      }
    }
  }
}

/**.......................................................................
 * Extend this image by nxPix and nyPix on a side
 */
void Image::extendBy(unsigned nxPix, unsigned nyPix)
{
  std::valarray<float> data = data_;
  std::valarray<unsigned> n = n_;

  unsigned nxOld = xAxis().getNpix();
  unsigned nyOld = yAxis().getNpix();

  unsigned nxNew = nxOld + nxPix;
  unsigned nyNew = nyOld + nyPix;

  resize(nxNew, nyNew);

  for(unsigned ix=0; ix < nxOld; ix++) {
    for(unsigned iy=0; iy < nyOld; iy++) {
      
      unsigned indOld = iy * nxOld + ix;
      unsigned indNew = iy * nxNew + ix;
      
      data_[indNew] = data[indOld];
      n_[indNew] = n[indOld];
    }
  }

}

/**.......................................................................
 * Blank the outer nxPix and nyPix of an image
 */
void Image::blankOuterEdge(unsigned nxPix, unsigned nyPix)
{
  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();

  for(unsigned ix=0; ix < nx; ix++) {

    for(unsigned iy=0; iy < ny; iy++) {
	
      if(ix < nxPix || ix >= nx-nxPix) {
	unsigned ind = iy*nx + ix;
	data_[ind] = 0.0;
      }

	
      if(iy < nyPix || iy >= ny-nyPix) {
	unsigned ind = iy*nx + ix;
	data_[ind] = 0.0;
      }
	
    }
  }
}

PGUTIL_COORD_CALLBACK(Image::angularCallback)
{
  Image* image = (Image*)args;

  // Convert from angular offset to absolute position

  HourAngle& ra  = image->getRa();
  Angle&     dec = image->getDec();

  dec.addDegrees(y);

  if(fabs(dec.degrees()) < 90)
    ra.addDegrees(x / cos(dec.radians()));

  std::ostringstream os;
  os << ra;
  xstr = os.str();

  os.str("");
  os << dec;
  ystr = os.str();

  return;
}

PGUTIL_COORD_CALLBACK(Image::latLongCallback)
{
  Image* image = (Image*)args;

  // Convert from angular offset to absolute position

  HourAngle& ra  = image->getRa();
  Angle&     dec = image->getDec();

  ra.addDegrees(x/image->xAxis().scale());
  dec.addDegrees(y/image->yAxis().scale());

  std::ostringstream os;
  os << (Angle&)ra;
  xstr = os.str();

  os.str("");
  os << dec;
  ystr = os.str();

  return;
}

PGUTIL_UNIT_CALLBACK(Image::unitCallback)
{
  Image* image = (Image*)args;

  // Convert from angular offset to absolute position

  HourAngle& ra  = image->getRa();
  Angle&     dec = image->getDec();

  std::ostringstream os;
  os << x << " deg";
  xstr = os.str();

  os.str("");
  os << y << " deg";
  ystr = os.str();

  return;
}

/**.......................................................................
 * Fill this image by coadding pixels from the passed image
 */
void Image::fillFrom(Image& image, Operation oper)
{
  if(!image.hasData_) {
    ThrowError("Passed image contains no data");
  }

  // Get the flat-sky approximation for the separation of the two
  // image reference pixels

  unsigned xnpix = xAxis().getNpix();
  unsigned ynpix = yAxis().getNpix();

  Angle xres = xAxis().getAngularResolution();
  Angle yres = yAxis().getAngularResolution();

  int xsense = xAxis().getSense();
  int ysense = yAxis().getSense();

  Angle xSep, ySep;
  getSeparationOfUsFromImage(xSep, ySep, image);

  //------------------------------------------------------------
  // Now iterate over all pixels of the new image, interpolating from
  // the passed image
  //------------------------------------------------------------

  Angle xOff, yOff;
  bool anyValidData = false;

  for(unsigned iYpix=0; iYpix < ynpix; iYpix++) {
    for(unsigned iXpix=0; iXpix < xnpix; iXpix++) {

      //------------------------------------------------------------
      // Calculate the offset of this pixel in the reference frame of
      // the passed image.  This will be the separation of the reference
      // pixel of the two images, plus the separation of this pixel
      // from our reference pixel.
      //------------------------------------------------------------

      xOff.setDegrees(xSep.degrees() + xsense * xres.degrees() * ((double)iXpix - (double)raRefPix_));
      yOff.setDegrees(ySep.degrees() + ysense * yres.degrees() * ((double)iYpix - (double)decRefPix_));
      
      unsigned ind = iYpix * xnpix + iXpix;

      double datum=0.0;
      bool valid;

      switch (oper) {
      case OPER_BIN:
	image.binData(xOff, xres, yOff, yres, datum, valid);
	break;
      case OPER_INTERPOLATE:
	image.interpolateData(xOff, yOff, datum, valid);
	break;
      default:
	ThrowError("Unrecognized operation: " << oper);
	break;
      }

      data_[ind] = datum;
      anyValidData = anyValidData || valid;
    }
  }

  units_    = image.units_;
  hasUnits_ = image.hasUnits_;

  if(anyValidData)
    hasData_  = true;
}

/**.......................................................................
 * Take the sqrt of an image
 */
Image Image::getSqrt()
{
  Image ret;
  ret = *this;

  for(unsigned i=0; i < data_.size(); i++) {
    ret.data_[i]  = data_[i] >= 0.0 ? ::sqrt(data_[i]) : 0.0;
    ret.valid_[i] = valid_[i] && data_[i] >= 0.0;
  }

  return ret;
}

/**.......................................................................
 * Take the sqrt of an image
 */
void Image::sqrt()
{
  for(unsigned i=0; i < data_.size(); i++) {
    data_[i]  = data_[i] >= 0.0 ? ::sqrt(data_[i]) : 0.0;
    valid_[i] = valid_[i] && data_[i] >= 0.0;
  }
}

/**.......................................................................
 * Take the natural log of an image
 */
void Image::ln()
{
  for(unsigned i=0; i < data_.size(); i++) {
    data_[i]  = data_[i] > 0.0 ? ::log(data_[i]) : 0.0;
    valid_[i] = valid_[i] && data_[i] > 0.0;
  }
}

/**.......................................................................
 * Take the log base 10 of an image
 */
void Image::log10()
{
  for(unsigned i=0; i < data_.size(); i++) {
    data_[i]  = data_[i] > 0.0 ? ::log10(data_[i]) : 0.0;
    valid_[i] = valid_[i] && data_[i] > 0.0;
  }
}

/**.......................................................................
 * Raise this image to a power
 */
void Image::power(double ex)
{
  // Is this a fractional power?

  double res  = abs(ex) - abs((int)ex);
  bool isFrac = res > 0.0;

  for(unsigned i=0; i < data_.size(); i++) {
    data_[i]  = (!isFrac || (isFrac && data_[i] >= 0.0)) ? ::pow(data_[i], ex) : 0.0;
    valid_[i] = valid_[i] && (!isFrac || (isFrac && data_[i] >= 0.0));
  }
}

void Image::replaceLessThanLimWithVal(double lim, double val)
{
  for(unsigned i=0; i < data_.size(); i++) {
    if(data_[i] < lim) {
      data_[i] = val;
    }
  }
}

Image Image::replaceZerosWithMin()
{
  Image ret;
  ret = *this;

  bool first=true;
  double nonZeroMin;

  for(unsigned i=0; i < data_.size(); i++) {
    if(data_[i] > 0.0) {
      if(first) {
	nonZeroMin = data_[i];
	first = false;
      } else {
	nonZeroMin = data_[i] < nonZeroMin ? data_[i] : nonZeroMin;
      }
    }
  }

  for(unsigned i=0; i < data_.size(); i++) {
    if(!(data_[i] > 0.0)) {
      data_[i] = nonZeroMin;
    }
  }

  return ret;
}

SolidAngle Image::getAngularResolution()
{
  return xAxis().getAngularResolution() * yAxis().getAngularResolution();
}

/**.......................................................................
 * Return an image comprised of the requested portion of this image
 */
Image Image::extract(unsigned nxpix, unsigned nypix)
{
  int ixmin, ixmax, iymin, iymax;
  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();

  if(nxpix % 2 == 0) {
    ixmin = (int)nx/2 - nxpix/2;
    ixmax = (int)nx/2 + nxpix/2 - 1;
  } else {
    ixmin = (int)nx/2 - nxpix/2;
    ixmax = (int)nx/2 + nxpix/2;
  }

  if(nypix % 2 == 0) {
    iymin = (int)ny/2 - nypix/2;
    iymax = (int)ny/2 + nypix/2 - 1;
  } else {
    iymin = (int)ny/2 - nypix/2;
    iymax = (int)ny/2 + nypix/2;
  }

  return extract(ixmin, ixmax, iymin, iymax);
}

/**.......................................................................
 * Return an image comprised of the requested portion of this image
 */
Image Image::extract(Angle xPos, Angle yPos, Angle dx, Angle dy, bool clip)
{
  int ixmin, ixmax, iymin, iymax;
  getPixelObj(xPos - dx, yPos - dy, ixmin, iymin, false, false);
  getPixelObj(xPos + dx, yPos + dx, ixmax, iymax, false, false);

  return extract(ixmin, ixmax, iymin, iymax, clip);
}

int Image::clipXAxis(int iPix)
{
  if(iPix < 0)
    return 0;

  if(iPix > xAxis().getNpix()-1)
    return xAxis().getNpix()-1;

  return iPix;
}

int Image::clipYAxis(int iPix)
{
  if(iPix < 0)
    return 0;

  if(iPix > yAxis().getNpix()-1)
    return yAxis().getNpix()-1;

  return iPix;
}

Image Image::extract(int ixmin, int ixmax, int iymin, int iymax, bool clip)
{
  Image image;
  image = *this;

  swapIfNeeded(ixmin, ixmax);
  swapIfNeeded(iymin, iymax);

  if(!clip) {
    if(ixmin < 0 || ixmax > xAxis().getNpix()-1 || iymin < 0 || iymax > yAxis().getNpix()-1)
      ThrowSimpleError("Invalid pixel range for this image: must be "  
		       << "(0 - " << xAxis().getNpix()-1 
		       << ", 0 - " << yAxis().getNpix()-1 << ")");
  } else {
    ixmin = clipXAxis(ixmin);
    ixmax = clipXAxis(ixmax);

    iymin = clipYAxis(iymin);
    iymax = clipYAxis(iymax);
  }

  unsigned nxDest = ixmax-ixmin+1;
  unsigned nyDest = iymax-iymin+1;

  image.xAxis().setNpix(nxDest);
  image.yAxis().setNpix(nyDest);

  if(xAxis().hasAngularSize() && yAxis().hasAngularSize()) {
    image.xAxis().setAngularSize(nxDest * xAxis().getAngularResolution());
    image.yAxis().setAngularSize(nyDest * yAxis().getAngularResolution());
  }

  unsigned nxSrc = xAxis().getNpix();
  unsigned nySrc = yAxis().getNpix();

  for(unsigned iy=iymin; iy <= iymax; iy++) {
    for(unsigned ix=ixmin; ix <= ixmax; ix++) {
      unsigned indSrc  = iy * nxSrc + ix;
      unsigned indDest = (iy-iymin) * nxDest + (ix-ixmin);
      image.data_[indDest] = data_[indSrc];
    }
  }

  // We are free to choose the reference pixel for the extracted
  // image.  For convenience, we choose to set the new reference
  // position to be the center of the extracted image

  image.raRefPix_  = (double)(ixmax - ixmin)/2;
  image.decRefPix_ = (double)(iymax - iymin)/2;
  image.setRaDec(getRa((double)(ixmin + ixmax)/2), getDec((double)(iymin + iymax)/2));

  return image;
}

/**.......................................................................
 * Return the angular separation of the ref pixel in the passed image
 * from our ref pixel
 */
void Image::getSeparationOfImageFromUs(Angle& xSep, Angle& ySep, Image& image)
{
  if(hasAbsolutePosition_ && image.hasAbsolutePosition_)
    gcp::util::Astrometry::flatSkyApproximationSeparations(xSep, ySep,
							   image.ra_, image.dec_,
							   ra_, dec_);
}

/**.......................................................................
 * Return the angular separation of our ref pixel from the passed
 * image's ref pixel
 */
void Image::getSeparationOfUsFromImage(Angle& xSep, Angle& ySep, Image& image)
{
  if(hasAbsolutePosition_ && image.hasAbsolutePosition_)
    gcp::util::Astrometry::flatSkyApproximationSeparations(xSep, ySep,
							   ra_, dec_,
							   image.ra_, image.dec_);
}

void Image::swapIfNeeded(double& min, double& max)
{
  if(min > max) {
    double tmp = min;
    min = max;
    max = tmp;
  }
}

void Image::swapIfNeeded(int& min, int& max)
{
  if(min > max) {
    int tmp = min;
    min = max;
    max = tmp;
  }
}

/**.......................................................................
 * Return the RA of the center pixel of the image
 */
HourAngle& Image::getRa()
{
  unsigned nx = xAxis().getNpix();
  return getRa((double)(nx)/2 - 0.5);
}

/**.......................................................................
 * Get the RA corresponding to the requested pixel value
 */
HourAngle& Image::getRa(double ixPix)
{
  if(!hasAbsolutePosition_)
    ThrowError("Image has no position specified");

  raCurr_.setDegrees(ra_.degrees() + (ixPix - raRefPix_) * xAxis().getAngularResolution().degrees() * xAxis().getSense());
  return raCurr_;
}

/**.......................................................................
 * Return the Dec of the center pixel of the image
 */
Declination& Image::getDec()
{
  unsigned ny = yAxis().getNpix();
  return getDec((double)(ny)/2 - 0.5);
}

/**.......................................................................
 * Return the Dec of the requested pixel
 */
Declination& Image::getDec(double iyPix)
{
  if(!hasAbsolutePosition_)
    ThrowError("Image has no position specified");

  Declination dec;
  decCurr_.setDegrees(dec_.degrees() + (iyPix - decRefPix_) * yAxis().getAngularResolution().degrees() * yAxis().getSense());
  return decCurr_;
}

void Image::flipY()
{
  unsigned ind, indflip;
  double dat, datflip;

  unsigned nx     = xAxis().getNpix();
  unsigned ny     = yAxis().getNpix();

  for(int iy=0; iy < ny/2; iy++) {
    for(int ix=0; ix < nx; ix++) {

      ind     =       iy  * nx + ix;
      indflip = (ny-1-iy) * nx + ix;

      dat     = data_[ind];
      datflip = data_[indflip];

      data_[indflip] = dat;
      data_[ind]     = datflip;
    }
  }
}

ImageAxis& Image::xImageAxis()
{
  return xAxis_;
}

ImageAxis& Image::yImageAxis()
{
  return yAxis_;
}

void Image::setInvalidPixelsTo(double val)
{
  for(unsigned i=0; i < data_.size(); i++) 
    if(!valid_[i])
      data_[i]  = val;
}

void Image::invalidateDataLessThan(double val)
{
  for(unsigned i=0; i < data_.size(); i++) 
    if(data_[i] < val)
      valid_[i] = 0;
}

void Image::invalidateDataGreaterThan(double val)
{
  for(unsigned i=0; i < data_.size(); i++) 
    if(data_[i] > val)
      valid_[i] = 0;
}

std::string Image::xLab()
{
  if(xAxis().hasAngularSize())
    return "Degrees";
  else
    return "Pixel Index";
}

std::string Image::yLab()
{
  if(yAxis().hasAngularSize())
    return "Degrees";
  else
    return "Pixel Index";
}

/**.......................................................................
 * For display purposes, the minimum value is the value required to
 * display the left edge of the first pixel.  If the raRefPix_ is the
 * pixel whose center is zero, then this is given by:
 *
 * - (raRefPix_ + 0.5) * xRes, i.e.:
 *
 *                    _______________________
 *                   |     |     |     |     |
 *                   |     |     |     |     |
 *                   |-----------------------|
 *                   |     |     |     |     |
 *                   |     |     |     |     |
 *                   |-----------------------|
 *                   |  0  |  1  | ref | n-1 |  n
 *                   |     |     |     |     |
 *                   |-----------------------|
 *                   |     |     |     |     |
 *                   |     |     |     |     |
 *                    -----------------------
 *                                           
 *                -->|              |<-- 
 *           
 *                   - (ref + 0.5)*dx          
 */
double Image::xMin()
{
  if(xAxis().hasAngularSize()) {

    double xRes = xAxis().getAngularResolution().degrees();
    return -xAxis().getSense() * (raRefPix_ + 0.5) * xRes;

  } else {
    return -0.5;
  }
}

/**.......................................................................
 * For display purposes, the maximum value is the value required to
 * display the right edge of the last pixel.  If the raRefPix_ is the
 * pixel whose center is zero, then this is given by:
 *
 *  (npix - 0.5 - raRefPix_) * xRes, i.e.:
 *
 *                    _______________________
 *                   |     |     |     |     |
 *                   |     |     |     |     |
 *                   |-----------------------|
 *                   |     |     |     |     |
 *                   |     |     |     |     |
 *                   |-----------------------|
 *                   |  0  |  1  | ref | n-1 |  n
 *                   |     |     |     |     |
 *                   |-----------------------|
 *                   |     |     |     |     |
 *                   |     |     |     |     |
 *                    -----------------------
 *                                           
 *                               -->|        |<-- 
 *           
 *                               (n - 0.5 - ref)*dx
 */
double Image::xMax()
{
  if(xAxis().hasAngularSize()) {

    double xRes = xAxis().getAngularResolution().degrees();
    return +xAxis().getSense() * (xAxis().getNpix() - raRefPix_ - 0.5) * xRes;

  } else {
    return xAxis().getNpix()-0.5;
  }
}
/**.......................................................................
 * For display purposes, the minimum value is the value required to
 * display the bottom edge of the first pixel.  If the decRefPix_ is the
 * pixel whose center is zero, then this is given by:
 *
 *  - (decRefPix_ + 0.5) * yRes, i.e.:
 *
 *
 *                                   n
 *                     _______________________
 *                    |     |     | n-1 |     |
 *                    |     |     |     |     |
 *                    |-----------------------|
 *                    |     |     | n-2 |     |
 *               |    |     |     |     |     |
 *              \|/   |-----------------------|
 *               -    |     |     | ref |     |
 *                    |     |     |     |     |
 *  -(ref + 0.5)*dy   |-----------------------|
 *                    |     |     |  0  |     |
 *                    |     |     |     |     |
 *               _     -----------------------
 *	        /|\
 *	         | 
 */
double Image::yMin()
{
  if(yAxis().hasAngularSize()) {
    double yRes =  yAxis().getAngularResolution().degrees();
    return -yAxis().getSense() * (decRefPix_ + 0.5) * yRes;
  } else {
    return -0.5;
  }
}

/**.......................................................................
 * For display purposes, the maximum value is the value required to
 * display the tio edge of the last pixel.  If the dexRefPix_ is the
 * pixel whose center is zero, then this is given by:
 *
 * (npix - 0.5 - decRefPix_) * yRes, i.e.:
 *
 *               |    
 *              \|/                  n
 *               -     _______________________
 *                    |     |     | n-1 |     |
 *                    |     |     |     |     |
 *                    |-----------------------|
 * (n - 0.5 - ref)*dy |     |     | n-2 |     |
 *                    |     |     |     |     |
 *                    |-----------------------|
 *               _    |     |     | ref |     |
 *              /|\   |     |     |     |     |
 *               |    |-----------------------|
 *                    |     |     |  0  |     |
 *                    |     |     |     |     |
 *                     -----------------------
 */
double Image::yMax()
{
  if(yAxis().hasAngularSize()) {

    double yRes = yAxis().getAngularResolution().degrees();
    return +yAxis().getSense() * (yAxis().getNpix() - decRefPix_ - 0.5) * yRes;

  } else {
    return yAxis().getNpix()-0.5;
  }
}

/**.......................................................................
 * Divide this image by a constant factor
 */
Image gcp::util::operator/(double fac, Image& image)
{
  Image res;
  res = image;
  for(unsigned i=0; i < res.data_.size(); i++) {
    res.data_[i]  = res.valid_[i] ? fac / res.data_[i] : 0.0;
  }
  return res;
}

void Image::getMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax, MaxType type)
{
  Angle xOffMin, xOffMax, yOffMin, yOffMax;

  xOffMin.setDegrees(xMin());
  xOffMax.setDegrees(xMax());

  yOffMin.setDegrees(yMin());
  yOffMax.setDegrees(yMax());

  getMax(val, xoff, yoff, iMax, xOffMin, xOffMax, yOffMin, yOffMax, type);
}

void Image::getMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax, Angle& dBound, MaxType type)
{
  if(dBound.degrees() > 0.0) {
    Angle xOffMin, xOffMax, yOffMin, yOffMax;
    
    xOffMin.setDegrees(-dBound.degrees());
    xOffMax.setDegrees(+dBound.degrees());
    
    yOffMin.setDegrees(-dBound.degrees());
    yOffMax.setDegrees(+dBound.degrees());
    
    getMax(val, xoff, yoff, iMax, xOffMin, xOffMax, yOffMin, yOffMax, type);
  } else {
    getMax(val, xoff, yoff, iMax, type);
  }
}

void Image::getMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax, 
		   Window& win, MaxType type)
{
  return getMax(val, xoff, yoff, iMax, win.xMin_, win.xMax_, win.yMin_, win.yMax_, type);
}

/**.......................................................................
 * Iterate over windows, finding the max
 */
void Image::getMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax, 
		   std::vector<Window>& windows, MaxType type)
{
  double currMax;
  Angle currXoff, currYoff;
  unsigned currImax;
  
  for(unsigned iWin=0; iWin < windows.size(); iWin++) {
    Window& win = windows[iWin];

    getMax(currMax, currXoff, currYoff, currImax, win, type);

    if(iWin == 0 || fabs(currMax) > fabs(val)) {
      val  = currMax;
      xoff = currXoff;
      yoff = currYoff;
      iMax = currImax;
    }
  }
}

void Image::getMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax, Angle& xOffMin, Angle& xOffMax, Angle& yOffMin, Angle& yOffMax, MaxType type)
{
  switch(type) {
  case TYPE_ABS:
    getAbsMax(val, xoff, yoff, iMax, xOffMin, xOffMax, yOffMin, yOffMax);
    break;
  case TYPE_POS:
    getMax(val, xoff, yoff, iMax, xOffMin, xOffMax, yOffMin, yOffMax);
    break;
  case TYPE_NEG:
    getMin(val, xoff, yoff, iMax, xOffMin, xOffMax, yOffMin, yOffMax);
    break;
  default:
    ThrowError("Undefined max type");
    break;
  }
}

void Image::getAbsMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax, Angle& xOffMin, Angle& xOffMax, Angle& yOffMin, Angle& yOffMax)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();

  double xmin = xMin();
  double ymin = yMin();

  double xRes =  xAxis().getAngularResolution().degrees();
  double yRes =  yAxis().getAngularResolution().degrees();

  int xSense = xAxis().getSense();
  int ySense = yAxis().getSense();

  if(fabs(xOffMax.degrees() - xOffMin.degrees()) < xRes || fabs(yOffMax.degrees() - yOffMin.degrees()) < yRes)
    ThrowSimpleColorError("Box is smaller than the resolution of this image", "red");

  double xMinBound = xOffMin.degrees();
  double xMaxBound = xOffMax.degrees();
  double yMinBound = yOffMin.degrees();
  double yMaxBound = yOffMax.degrees();
  double xOffCurr, yOffCurr;

  bool first = true;

  for(unsigned iy=0; iy < ny; iy++) {
    for(unsigned ix=0; ix < nx; ix++) {

      unsigned ind  = iy * nx + ix;
      xOffCurr = xmin + ((double)(ix) + 0.5) * xRes;
      yOffCurr = ymin + ((double)(iy) + 0.5) * yRes;

      if(xOffCurr >= xMinBound && xOffCurr <= xMaxBound &&
	 yOffCurr >= yMinBound && yOffCurr <= yMaxBound) {

	if(fabs(data_[ind]) > fabs(val) || first) {
	  val  = data_[ind];
	  xoff.setDegrees(xOffCurr);
	  yoff.setDegrees(yOffCurr);
	  iMax = ind;
	  first = false;
	}

      }
    }
  }
}

void Image::getMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax, Angle& xOffMin, Angle& xOffMax, Angle& yOffMin, Angle& yOffMax)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();

  double xmin = xMin();
  double ymin = yMin();

  double xRes =  xAxis().getAngularResolution().degrees();
  double yRes =  yAxis().getAngularResolution().degrees();

  if(fabs(xOffMax.degrees() - xOffMin.degrees()) < xRes || fabs(yOffMax.degrees() - yOffMin.degrees()) < yRes)
    ThrowSimpleColorError("Box is smaller than the resolution of this image", "red");

  double xMinBound = xOffMin.degrees();
  double xMaxBound = xOffMax.degrees();
  double yMinBound = yOffMin.degrees();
  double yMaxBound = yOffMax.degrees();
  double xOffCurr, yOffCurr;

  bool first = true;

  for(unsigned iy=0; iy < ny; iy++) {
    for(unsigned ix=0; ix < nx; ix++) {

      unsigned ind  = iy * nx + ix;
      xOffCurr = xmin + ((double)(ix) + 0.5) * xRes;
      yOffCurr = ymin + ((double)(iy) + 0.5) * yRes;

      if(xOffCurr >= xMinBound && xOffCurr <= xMaxBound &&
	 yOffCurr >= yMinBound && yOffCurr <= yMaxBound) {

	if(data_[ind] > val || first) {
	  val  = data_[ind];
	  xoff.setDegrees(xOffCurr);
	  yoff.setDegrees(yOffCurr);
	  iMax = ind;
	  first = false;
	}
	
      }
    }
  }
}

void Image::getMin(double& val, Angle& xoff, Angle& yoff, unsigned& iMax, Angle& xOffMin, Angle& xOffMax, Angle& yOffMin, Angle& yOffMax)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();

  double xmin = xMin();
  double ymin = yMin();

  double xRes =  xAxis().getAngularResolution().degrees();
  double yRes =  yAxis().getAngularResolution().degrees();

  if(fabs(xOffMax.degrees() - xOffMin.degrees()) < xRes || fabs(yOffMax.degrees() - yOffMin.degrees()) < yRes)
    ThrowSimpleColorError("Box is smaller than the resolution of this image", "red");

  double xMinBound = xOffMin.degrees();
  double xMaxBound = xOffMax.degrees();
  double yMinBound = yOffMin.degrees();
  double yMaxBound = yOffMax.degrees();
  double xOffCurr, yOffCurr;

  bool first = true;

  for(unsigned iy=0; iy < ny; iy++) {
    for(unsigned ix=0; ix < nx; ix++) {

      unsigned ind  = iy * nx + ix;
      xOffCurr = xmin + ((double)(ix) + 0.5) * xRes;
      yOffCurr = ymin + ((double)(iy) + 0.5) * yRes;

      if(xOffCurr >= xMinBound && xOffCurr <= xMaxBound &&
	 yOffCurr >= yMinBound && yOffCurr <= yMaxBound) {

	if(data_[ind] < val || first) {
	  val  = data_[ind];
	  xoff.setDegrees(xOffCurr);
	  yoff.setDegrees(yOffCurr);
	  iMax = ind;
	  first = false;
	}

      }
    }
  }
}

void Image::getAbsMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();

  double xmin = xMin();
  double ymin = yMin();

  double xRes =  xAxis().getAngularResolution().degrees();
  double yRes =  yAxis().getAngularResolution().degrees();

  double xOffCurr, yOffCurr;

  bool first = true;

  for(unsigned iy=0; iy < ny; iy++) {
    for(unsigned ix=0; ix < nx; ix++) {

      unsigned ind  = iy * nx + ix;
      xOffCurr = xmin + ((double)(ix) + 0.5) * xRes;
      yOffCurr = ymin + ((double)(iy) + 0.5) * yRes;

      if(fabs(data_[ind]) > fabs(val) || first) {
	val  = data_[ind];
	xoff.setDegrees(xOffCurr);
	yoff.setDegrees(yOffCurr);
	iMax = ind;
	first = false;
      }
    }
  }
}

void Image::getMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();

  double xmin = xMin();
  double ymin = yMin();

  double xRes =  xAxis().getAngularResolution().degrees();
  double yRes =  yAxis().getAngularResolution().degrees();

  double xOffCurr, yOffCurr;

  bool first = true;

  for(unsigned iy=0; iy < ny; iy++) {
    for(unsigned ix=0; ix < nx; ix++) {

      unsigned ind  = iy * nx + ix;
      xOffCurr = xmin + ((double)(ix) + 0.5) * xRes;
      yOffCurr = ymin + ((double)(iy) + 0.5) * yRes;

      if(data_[ind] > val || first) {
	val  = data_[ind];
	xoff.setDegrees(xOffCurr);
	yoff.setDegrees(yOffCurr);
	iMax = ind;
	first = false;
      }
      
    }
  }
}

void Image::getMin(double& val, Angle& xoff, Angle& yoff, unsigned& iMax)
{
  if(!hasData_) {
    ThrowError("Image contains no data");
  }

  unsigned nx = xAxis().getNpix();
  unsigned ny = yAxis().getNpix();

  double xmin = xMin();
  double ymin = yMin();

  double xRes =  xAxis().getAngularResolution().degrees();
  double yRes =  yAxis().getAngularResolution().degrees();

  double xOffCurr, yOffCurr;

  bool first = true;

  for(unsigned iy=0; iy < ny; iy++) {
    for(unsigned ix=0; ix < nx; ix++) {

      unsigned ind  = iy * nx + ix;
      xOffCurr = xmin + ((double)(ix) + 0.5) * xRes;
      yOffCurr = ymin + ((double)(iy) + 0.5) * yRes;

      if(data_[ind] < val || first) {
	val  = data_[ind];
	xoff.setDegrees(xOffCurr);
	yoff.setDegrees(yOffCurr);
	iMax = ind;
	first = false;
      }
      
    }
  }
}

std::ostream& gcp::util::operator<<(std::ostream& os, const Image::Window& win)
{
  os << "xmin = " << win.xMin_ << std::endl;
  os << "xmax = " << win.xMax_ << std::endl;
  os << "ymin = " << win.yMin_ << std::endl;
  os << "ymax = " << win.yMax_ << std::endl;

  return os;
}

/**.......................................................................
 * Return the window corresponding to the range specifications
 */
Image::Window Image::getWindow(String& xRange, String& yRange, String& units)
{
  Angle xMin, xMax, yMin, yMax;
  parseRange(xRange, xMin, xMax, units, true);
  parseRange(yRange, yMin, yMax, units, false);

  Window win;

  win.xMin_.setDegrees(xMin.degrees());
  win.xMax_.setDegrees(xMax.degrees());
  win.yMin_.setDegrees(yMin.degrees());
  win.yMax_.setDegrees(yMax.degrees());

  return win;
}  
 
/**.......................................................................
 * Parse a range specification
 */
void Image::parseRange(String& range, Angle& rangeMin, Angle& rangeMax, String& units, bool isX)
{
  if(range.contains("+-")) 
    parseSymmRange(range, rangeMin, rangeMax, units, isX);
  else if(range.contains(":")) 
    parseFullRange(range, rangeMin, rangeMax, units);
  else
    ThrowError("Unrecognized range specification: " << range);
}

/**.......................................................................
 * Parse a range specification of the form val +- delt, allowing for
 * val to be a symbolic name like:
 *
 *   min -- The minimum
 *   max -- The maximum
 *   abs -- The absolute maximum
 */
void Image::parseSymmRange(String& range, Angle& rangeMin, Angle& rangeMax, String& units, bool isX)
{
  std::ostringstream os;
  String centerStr = range.findNextInstanceOf(" ", false, "+-", true, false);
  String deltaStr  = range.findNextInstanceOf("+-", true, " ", false, false);

  double val;
  Angle xPeak, yPeak, center, delta;
  unsigned iMax;

  if(centerStr.contains("min")) {
    getMin(val, xPeak, yPeak, iMax);
    delta.setVal(deltaStr.toDouble(), units.str());
  } else if(centerStr.contains("max")) {
    getMax(val, xPeak, yPeak, iMax);
    delta.setVal(deltaStr.toDouble(), units.str());
  } else if(centerStr.contains("abs")) {
    getAbsMax(val, xPeak, yPeak, iMax);
    delta.setVal(deltaStr.toDouble(), units.str());
  } else {
    center.setVal(centerStr.toDouble(), units.str());
    delta.setVal(deltaStr.toDouble(), units.str());
    rangeMin.setDegrees(center.degrees() - delta.degrees());
    rangeMax.setDegrees(center.degrees() + delta.degrees());
    return;
  }

  if(isX) {
    rangeMin.setDegrees(xPeak.degrees() - delta.degrees());
    rangeMax.setDegrees(xPeak.degrees() + delta.degrees());
  } else {
    rangeMin.setDegrees(yPeak.degrees() - delta.degrees());
    rangeMax.setDegrees(yPeak.degrees() + delta.degrees());
  }

}

/**.......................................................................
 * Parse a range specification of the form min:max
 */
void Image::parseFullRange(String& range, Angle& rangeMin, Angle& rangeMax, String& units)
{
  String minStr  = range.findNextInstanceOf(" ", false, ":", true, false);
  String maxStr  = range.findNextInstanceOf(":", true, " ", false, false);

  rangeMin.setVal(minStr.toDouble(), units.str());
  rangeMax.setVal(maxStr.toDouble(), units.str());
}

void Image::getPosition(int ix, int iy, Angle& xPos, Angle& yPos)
{
  // The reference pixel is the zero-point of the array, thus the
  // position of pixel ix, iy is returned as the center of pixel
  // ix, iy

  xPos.setRadians((((double)(ix) -  raRefPix_) * xAxis().getSense()) * xAxis().getAngularResolution().radians());
  yPos.setRadians((((double)(iy) - decRefPix_) * yAxis().getSense()) * yAxis().getAngularResolution().radians());
}

/**.......................................................................
 * Return the pixel corresponding to the specified offset.  Note that
 * the pixel 'value' is treated as the center of the pixel, thus the
 * reference pixel spans refPix - 0.5 <= val < refPix + 0.5.  We
 * therefore add 0.5 to the offset from the reference pixel and take
 * the floor of the resulting value to get the corresponding pixel
 * offset
 *
 *                    -----------------------
 *                   |  0  |  1  | ref | n-1 |  n
 *                   |     |     |     |     |
 *                    -----------------------
 */
void Image::getPixelObj(Angle xPos, Angle yPos, int& ix, int& iy, bool clip, bool doThrow)
{
  return getPixelRef(xPos, yPos, ix, iy, clip, doThrow);
}

void Image::getPixelRef(Angle& xPos, Angle& yPos, int& ix, int& iy, bool clip, bool doThrow)
{
  // The reference pixel is always the zero-point of the array

  double xVal = floor((xPos.radians() / xAxis().getAngularResolution().radians() * xAxis().getSense()) +  raRefPix_ + 0.5);
  double yVal = floor((yPos.radians() / yAxis().getAngularResolution().radians() * yAxis().getSense()) + decRefPix_ + 0.5);

  // If requested to clip, clip the image at the pixel boundaries

  if(clip) {

    if(xVal > xAxis().getNpix()-1)
      xVal = xAxis().getNpix()-1;
    if(yVal > yAxis().getNpix()-1)
      yVal = yAxis().getNpix()-1;
    if(xVal < 0)
      xVal = 0;
    if(yVal < 0)
      yVal = 0;

    // Else throw if the requested position doesn't lie on this image
    // and we were requested to throw.  Else we return the effective
    // pixel value, even if it doesn't lie on this image.

  } else if(doThrow) {
    if(xVal > xAxis().getNpix()-1 ||  yVal > yAxis().getNpix()-1 ||
       xVal < 0 || yVal < 0)
      ThrowSimpleColorError("Position: " << xPos << ", " << yPos << " lies outside of this image", "red");
  }

  ix = xVal;
  iy = yVal;
}

/**.......................................................................
 * Return the pixel deltas corresponding to the specified angular
 * deltas.  An offset of 0.5 in pixel units puts us at the lower
 * boundary of the next pixel, hence the + 0.5
 */
void Image::getPixelDelta(Angle xDelta, Angle yDelta, int& xPixDelta, int& yPixDelta)
{
  xPixDelta = floor((xDelta / xAxis().getAngularResolution()) + 0.5);
  yPixDelta = floor((yDelta / yAxis().getAngularResolution()) + 0.5);
}
