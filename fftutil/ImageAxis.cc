#include "gcp/fftutil/ImageAxis.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
ImageAxis::ImageAxis() 
{
  n_               = 0;
  nIsSpecified_    = false;
  sizeIsSpecified_ = false;
  sense_           = 1;
  scale_           = 1.0;
  name_            = "";
}

/**.......................................................................
 * Destructor.
 */
ImageAxis::~ImageAxis() {}


/**.......................................................................
 * Set the type of this axis
 */
void ImageAxis::setAxisType(unsigned type)
{
  type_ = type;
}

//-----------------------------------------------------------------------
// Npix methods
//-----------------------------------------------------------------------

void ImageAxis::setNpix(unsigned n)
{
  n_            = n;
  nIsSpecified_ = true;
}

bool ImageAxis::hasNpix()
{
  return nIsSpecified_;
}

unsigned ImageAxis::getNpix()
{
  if(nIsSpecified_) {
    return n_;
  } else {
    ThrowError("Number of pixels isn't specified for this axis: use setNpix()");
  }
}

//-----------------------------------------------------------------------
// Angular size methods
//-----------------------------------------------------------------------

void ImageAxis::setAngularSize(Angle size)
{
  // If the size was passed as a negative angle, set the sense to -1
  // However, don't change the sense if it was specified as a positive
  // angle.  This is because we often resize images (always with a
  // positive size) but want to preserve the sense of the axis when we
  // do.

  if(size.degrees() < 0.0) {
    size_.setDegrees(-1 * size.degrees());
    setSense(-1);
  } else {
    size_  =  size;
  }

  sizeIsSpecified_ = true;
}

void ImageAxis::setSense(int sense)
{
  sense_ = sense;
}

bool ImageAxis::hasAngularSize()
{
  return sizeIsSpecified_;
}

Angle& ImageAxis::getAngularSize()
{
  if(sizeIsSpecified_) {
    return size_;
  } else {
    ThrowError("No size is specified for this axis: use setAngularSize()");
  }
}

int ImageAxis::getSense()
{
  return sense_;
}

std::string ImageAxis::getProjection()
{
  return projection_;
}

void ImageAxis::setProjection(std::string projection)
{
  projection_ = projection;
}

/**.......................................................................
 * Set the spatial frequency resolution of this image
 */
void ImageAxis::setSpatialFrequencyResolution(double inverseRadians)
{
  Angle size;
  size.setRadians(1.0/inverseRadians);

  // Don't call inheritors methods -- we only want to couple
  // base-class methods here

  ImageAxis::setAngularSize(size);
}

/**.......................................................................
 * Return the angular resolution of this image
 */
Angle& ImageAxis::getAngularResolution()
{
  if(sizeIsSpecified_ && nIsSpecified_) {
    resolution_.setRadians(size_.radians() / n_);
    return resolution_;
  } else {
    ThrowError("Not enough information to calculate the angular resolution of this axis");
  }
}

/**.......................................................................
 * Return the spatial frequency resolution of this image
 */
double ImageAxis::getSpatialFrequencyResolution()
{
  if(!sizeIsSpecified_)
    ThrowError("No size is specified for this axis: use setAngularSize()");

  return 1.0/size_.radians();
}

/**.......................................................................
 * Return the minimum spatial frequency in this image
 */
double ImageAxis::getMinimumSpatialFrequency()
{
  return getSpatialFrequencyResolution();
}

/**.......................................................................
 * Return the maximum spatial frequency in this image
 */
double ImageAxis::getMaximumSpatialFrequency()
{
  if(!sizeIsSpecified_)
    ThrowError("No size is specified for this axis: use setAngularSize()");

  return getNpix() * 1.0/(2*size_.radians());
}

/**.......................................................................
 * Return true if the axes are equivalent in size and resolution
 */
bool ImageAxis::operator==(const ImageAxis& axis)
{
  return operator==((ImageAxis&) axis);
}

bool ImageAxis::operator==(ImageAxis& axis)
{
  bool sameNpix  = getNpix() == axis.getNpix();
  bool sameSize  = fabs((getAngularSize().radians() - axis.getAngularSize().radians()) / getAngularSize().radians()) < 1e-12;
  bool sameSense = (sense_ == axis.sense_);

  return sameNpix && sameSize && sameSense;
}

ImageAxis::ImageAxis(const ImageAxis& axis)
{
  *this = axis;
}

ImageAxis::ImageAxis(ImageAxis& axis)
{
  *this = axis;
}

void ImageAxis::operator=(const ImageAxis& axis)
{
  return operator=((ImageAxis&)axis);
}

void ImageAxis::operator=(ImageAxis& axis)
{
  size_            = axis.size_;
  resolution_      = axis.resolution_;
  n_               = axis.n_;
  sizeIsSpecified_ = axis.sizeIsSpecified_;
  nIsSpecified_    = axis.nIsSpecified_;
  type_            = axis.type_;
  sense_           = axis.sense_;
}

void ImageAxis::setScale(double scale, std::string name)
{
  scale_ = scale;
  name_  = name;
}

double ImageAxis::scale()
{
  return scale_;
}

std::string ImageAxis::name()
{
  if(scale_ == 1.0)
    return hasAngularSize() ? "Degrees" : "Pixel Index";
  else
    return name_;
}
