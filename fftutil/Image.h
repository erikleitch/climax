// $Id: Image.h,v 1.8 2012/05/31 16:22:56 eml Exp $

#ifndef GCP_UTIL_IMAGE_H
#define GCP_UTIL_IMAGE_H

/**
 * @file Image.h
 * 
 * Tagged: Tue Jun  1 22:32:16 PDT 2010
 * 
 * @version: $Revision: 1.8 $, $Date: 2012/05/31 16:22:56 $
 * 
 * @author tcsh: Erik Leitch.
 */
#include <string>
#include <valarray>
#include <vector>

#include "gcp/util/Angle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Length.h"
#include "gcp/util/SolidAngle.h"
#include "gcp/util/Unit.h"

#include "gcp/fftutil/Generic2DObject.h"
#include "gcp/fftutil/ImageAxis.h"
#include "gcp/fftutil/Model.h"

#include "gcp/pgutil/PgUtil.h"

namespace gcp {
  namespace util {

    class Dft2d;
    class ObsInfo;

    class Image : public Generic2DObject {
    public:

      enum MaxType {
	TYPE_NEG,
	TYPE_POS,
	TYPE_ABS
      };

      enum Param {
	PARAM_NONE   = 0x0,
	PARAM_UNITS  = 0x1,
	PARAM_JYCONV = 0x2,
      };

      enum Operation {
	OPER_NONE        = 0x0,
	OPER_ADD         = 0x1,
	OPER_AVERAGE     = 0x2,
	OPER_ASSIGN      = 0x4,
	OPER_BIN         = 0x8,
	OPER_INTERPOLATE = 0x10,
      };

      struct Window {
	Angle xMin_;
	Angle xMax_;
	Angle yMin_;
	Angle yMax_;

	friend std::ostream& operator<<(std::ostream& os, const Window& win);
      };

      struct Axis : public ImageAxis {

	Image* parent_;

	void setNpix(unsigned n);
	void operator=(ImageAxis& axis);
	void operator=(const ImageAxis& axis);
	void operator=(Axis& axis);
	void operator=(const Axis& axis);
      };

      // Constructors

      Image();
      Image(const Image& image);
      Image(Image& image);
      Image(std::vector<float>&, Unit::Units units);
      Image(unsigned npix, Angle size);

      void operator=(const Image& image);
      void operator=(Image& image);
      void initialize();

      // Destructor.

      virtual ~Image();

      // Methods of loading data into this image
      
      void initializeFromFitsFile(std::string fileName, ObsInfo* obs=0);
      void initializeFromFitsFile(std::string fileName, unsigned iHdu, ObsInfo* obs=0);
      void initializeFromFitsTable(std::string fileName, std::string colName);
      void initializeFromFitsTable(std::string fileName, std::string extName, std::string xColName, std::string yColName, std::string dataColName, ObsInfo* obs=0);
      void initializeFromBinaryImageFile(std::string fileName);
      void initializeFromArray(double* data, unsigned ndata, Unit::Units units=Unit::UNITS_UNKNOWN);
      void initializeFromArray(std::vector<float>& data, Unit::Units units=Unit::UNITS_UNKNOWN);
      void initializeFromArray(std::valarray<double>& data, Unit::Units units=Unit::UNITS_UNKNOWN);
      void initialize(Dft2d& dft);
      
      void writeToFitsFile(std::string fileName);

      // Methods of internally filling the image

      // Create a delta function image

      void createDeltaFunctionImage(unsigned nx, unsigned ny, float val=1.0, int xval=-1, int yval=-1);
      void createDeltaFunctionImage(float val=1.0, 
				    Angle xval=Angle(Angle::Degrees(), 0.0),
				    Angle yval=Angle(Angle::Degrees(), 0.0));

      void createUniformImage(unsigned nx, unsigned ny, float value=1.0);
      void createEdgeImage(unsigned nx, unsigned ny);

      void blankOuterEdge(unsigned nxPix, unsigned nyPix);
      void extendBy(unsigned nxPix, unsigned nyPix);

      void createGaussianImage(unsigned nx, unsigned ny, double sigma);
      void createGaussianImage(unsigned nx, unsigned ny, double sigmax, double sigmay);
      void createGaussianImage(double amp, Angle sigma, Angle xOffset=zero_, Angle yOffset=zero_);

      void createGaussianImageFullSpecification(double amp, 
						Angle majSigma,      Angle minSigma, 
						Angle rotAngle,
						Angle xOffset=zero_, Angle yOffset=zero_);

      void createSineImage(  unsigned nx, unsigned ny, double period, Angle phase);
      void createCosineImage(unsigned nx, unsigned ny, double period, Angle phase);

      void addSine(double period, Angle phase);      
      void addCosine(double period, Angle phase);      

      // These routines require that the user have already specified
      // an image size and npix

      void createGaussianPrimaryBeam(Length& diameter, Frequency& freq);
      void createUniformAperturePrimaryBeam(Length& diameter, Frequency& freq);
      static Angle gaussianSigma(Length& length, Frequency& freq);

      // Convolve this image with another

      void convolve(Image& im, bool zeropad=false, Dft2d* dft=0);

      // Resize this image

      void resize(unsigned nx, unsigned ny);

      // Initialize this object from another Image object

      void initialize(Image& image);

      // Or explicitly

      void initialize(Angle& xSize, Angle& ySize, unsigned nXPix, unsigned nYPix);

      // Set the parameters of this image

      void setNativeToJy(double jyconv);
      void setNativeToJy(Frequency& nu, SolidAngle& beam);

      void setAngularSize(Angle size);
      void setNpix(unsigned nPix);
      SolidAngle getAngularResolution();
      HourAngle&   getRa(double ixPix);
      Declination& getDec(double iyPix);
      HourAngle&   getRa();
      Declination& getDec();

      // Return an axis descriptor

      Axis& xAxis();
      Axis& yAxis();

      //------------------------------------------------------------
      // Main method to display this image.  
      // 
      // If normalDisplay = true,  display axes increasing to the right (top).  
      // If normalDisplay = false, display x-axis increasing to the left, per difmap.
      // If greyScale = false, make a contour plot instead of greyscale
      //------------------------------------------------------------

      void display(bool normalDisplay, bool greyScale);

      // Normal greyscale display

      void display(bool setWin=true);

      void logDisplay();
      void logDisplay(bool normalDisplay, bool greyScale);

      // Difmap-style greyscale display

      void difmapDisplay();

      // Normal-style contour display

      void contour();

      void plotRow(unsigned iRow);
      void plotColumn(unsigned iCol);

      double xMin();
      double xMax();
      double yMin();
      double yMax();

      std::string xLab();
      std::string yLab();

      // Return statistics of the pixels in this image

      double sum();
      double mean();
      double min();
      double minGreaterThanZero();
      double max();
      double rms();

      //------------------------------------------------------------
      // Functions to find the min/max/absmax for various window
      // specifications
      //------------------------------------------------------------

      void getMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax, MaxType type);
      void getMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax, Angle& dBound, MaxType type);
      void getMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax, 
	          Angle& xOffMin, Angle& xOffMax, Angle& yOffMin, Angle& yOffMax, MaxType type);
      void getMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax, 
		  Window& window, MaxType type);

      void getMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax, 
		  std::vector<Window>& windows, MaxType type);

      void getAbsMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax, 
	          Angle& xOffMin, Angle& xOffMax, Angle& yOffMin, Angle& yOffMax);
      void getMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax, 
	          Angle& xOffMin, Angle& xOffMax, Angle& yOffMin, Angle& yOffMax);
      void getMin(double& val, Angle& xoff, Angle& yoff, unsigned& iMax, 
	          Angle& xOffMin, Angle& xOffMax, Angle& yOffMin, Angle& yOffMax);

      void getAbsMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax);
      void getMax(double& val, Angle& xoff, Angle& yoff, unsigned& iMax);
      void getMin(double& val, Angle& xoff, Angle& yoff, unsigned& iMax);

      // Utility functions for parsing windows

      void parseFullRange(String& range, Angle& rangeMin, Angle& rangeMax, String& units);
      void parseSymmRange(String& range, Angle& rangeMin, Angle& rangeMax, String& units, bool isX);
      void parseRange(String& range, Angle& rangeMin, Angle& rangeMax, String& units, bool isX);
      Window getWindow(String& xRange, String& yRange, String& units);

      // Statistics within a specified annulus

      double mean(Angle& rmin, Angle& rmax);
      double rms(Angle& rmin, Angle& rmax);

      // Scale operators for Image class

      void operator+=(double incr);
      void operator-=(double decr);
      void operator*=(double mult);
      void operator/=(double div);
      void operator=(double val);

      Image operator+(double incr);
      Image operator-(double decr);
      Image operator*(double mult);
      Image operator/(double div);

      // Image operators for Image class

      void operator-=(Image& image);
      void operator+=(Image& image);
      void operator*=(const Image& image);
      void operator*=(Image& image);
      void operator/=(Image& image);

      Image operator+(Image& image);
      Image operator-(Image& image);
      Image operator*(Image& image);
      Image operator/(Image& image);

      Image getSqrt();
      void sqrt();
      void ln();
      void log10();
      void power(double ex);

      Image replaceZerosWithMin();
      void replaceLessThanLimWithVal(double lim, double val);

      bool axesAreEquivalent(Image& image);

      std::string unitToString();

      // Zero the image (set all pixel values to zero)

      void zero();
      void setValue(float val);

      void zeroBeyondRadius(Angle radius);
      void invalidateBeyondRadius(Angle radius);
      void invalidateDataLessThan(double val);
      void invalidateDataGreaterThan(double val);

      void invalidate();

      // Return a reference to the specified pixel

      float& val(unsigned ix, unsigned iy);

      // Fill this image from the passed image

      void fillFrom(Image& image, Operation oper);
      void fillFrom(Angle& axisUnits, std::valarray<double>& x, std::valarray<double>& y, std::valarray<double>& data);
      void fillFrom(Angle& axisUnits, std::vector<double>& x, std::vector<double>& y, std::vector<double>& data);

      // Interpolate data from this image at the requested offset from
      // the center

      void interpolateData(Angle& xOff, Angle& yOff, double& val, bool& valid, bool doThrow=true);
      void getNearestData(Angle& xOff, Angle& yOff, double& val, bool& valid, bool doThrow=true);
      void binData(Angle& xOff, Angle& xRes, Angle& yOff, Angle& yRes, double& val, bool& valid);

      // Assignment operators that only assign data

      void assignDataFrom(const Image& image, unsigned iXStart=0, unsigned iYStart=0);
      void assignDataFrom(Image& image, unsigned iXStart=0, unsigned iYStart=0);

      void copyAncillaryData(Image& image);


      void setRaDec(HourAngle ra, Declination dec);
      void setRaDecFft(HourAngle ra, Declination dec);
      void setLatLng(Angle lat, Angle lon);
      void setRaDec(HourAngle ra, Declination dec, double raRefPix, double decRefPix);

      static PGUTIL_COORD_CALLBACK(angularCallback);
      static PGUTIL_COORD_CALLBACK(latLongCallback);
      static PGUTIL_UNIT_CALLBACK(unitCallback);

      int clipXAxis(int iPix);
      int clipYAxis(int iPix);

      Image extract(int ixmin, int ixmax, int iymin, int iymax, bool clip=false);
      Image extract(Angle xPos, Angle yPos, Angle dx, Angle dy, bool clip=false);
      Image extract(unsigned nx, unsigned ny);
      void getSeparationOfImageFromUs(Angle& xSep, Angle& ySep, Image& image);
      void getSeparationOfUsFromImage(Angle& xSep, Angle& ySep, Image& image);
      void flipY();
      void assignWithValidation(Image& image);
      void addWithValidation(Image& image);
      void setInvalidPixelsTo(double val);

    public:

      friend class Dft2d;

      void addParameter(Param param);
      void checkParameter(Param param);

      void addImage(Image& image, Operation oper);

      std::valarray<float>    data_;  // The actual data
      std::valarray<float>    wtSum_; // Utility for storing weights
      std::valarray<unsigned> n_;     // number of pixels in the mean
      std::valarray<unsigned> valid_; // True if this pixel contains valid data (only used by xxxWithValidation() methods)

      static const Angle zero_;

      Axis xAxis_;                // Parameters of the x-axis
      Axis yAxis_;                // Parameters of the y-axis

      double nativeToJy_;         // A conversion factor from native
				  // units to Janskys
      bool hasNativeToJy_;

      float dataMin_;
      float dataMax_;

      unsigned axes_;

      //------------------------------------------------------------
      // Members for images with absolute position information
      //------------------------------------------------------------

      // By default, I associate this position with the literal
      // center of the image.  If the axes are even numbers of pixels, that
      // means the position midway between the center pixels.  If an
      // odd number of pixels, this is the position of the center pixel

      double raRefPix_;          // The reference pixel for this position
      double decRefPix_;         // The reference pixel for this position

      HourAngle raCurr_;         // A handle used to pass a requested RA to external callers, for optimization
      Declination decCurr_;      // A handle used to pass a requested DEC to external callers, for optimization

      void resize();

      void updateAxis(unsigned axis);
      void checkAxes();
      void swapIfNeeded(double& min, double& max);
      void swapIfNeeded(int& min, int& max);

      bool isLatLng_;

      // Get the angular coordinate of the specified pixel

      void getPosition(int ix, int iy, Angle& xPos, Angle& yPos);
      void getPixelRef(Angle& xPos, Angle& yPos, int& ix, int& iy, bool clip=false, bool doThrow=true);
      void getPixelObj(Angle xPos, Angle yPos, int& ix, int& iy, bool clip=false, bool doThrow=true);
      void getPixelDelta(Angle xDelta, Angle yDelta, int& xpixDelta, int& yPixDelta);

      void initializeRefPixFftConvention();
      void initializeRefPixNonFftConvention();

    private:

      ImageAxis& xImageAxis();
      ImageAxis& yImageAxis();

    }; // End class Image

    Image operator/(double fac, Image& image);

    std::ostream& operator<<(std::ostream& os, const Image::Window& win);    

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_IMAGE_H
