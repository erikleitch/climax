// $Id: Dft2d.h,v 1.8 2012/05/29 16:20:28 eml Exp $

#ifndef GCP_UTIL_DFT2D_H
#define GCP_UTIL_DFT2D_H

/**
 * @file Dft2d.h
 * 
 * Tagged: Thu Jun  3 11:41:14 PDT 2010
 * 
 * @version: $Revision: 1.8 $, $Date: 2012/05/29 16:20:28 $
 * 
 * @author Erik Leitch
 */
#include <fftw3.h>

#include "gcp/fftutil/Image.h"
#include "gcp/fftutil/ImageAxis.h"

#include "gcp/util/Mutex.h"
#include "gcp/util/Wavelength.h"

namespace gcp {
  namespace util {

    class Dft2d : public Generic2DObject {
    public:

      enum DataType {
	DATA_REAL,
	DATA_IMAG,
	DATA_ABS,
	DATA_ABS_LOG,
	DATA_UV
      };

      struct Axis : public ImageAxis {

	ImageAxis logicalAxis_;
	Dft2d* parent_;
	unsigned zeropadFactor_;

	Axis();
	void setZeropadFactor(unsigned zeropadFactor);
	void setNpix(unsigned n);
	void setAngularSize(Angle size);
	void setSpatialFrequencyResolution(double inverseRadians);
	double getSpatialFrequencyResolution();

	void operator=(const Axis& axis);
	void operator=(Axis& axis);
	void operator=(ImageAxis& axis);
      };

      // Constructor.

      Dft2d(bool optimize=true);
      Dft2d(int nx, int ny, bool optimize=true);
      Dft2d(Image& image, bool optimize=true);

      Dft2d(const Dft2d& dft);
      void operator=(const Dft2d& dft);
      void operator=(Dft2d& dft);

      virtual void operator*=(double mult);


      // Destructor.

      virtual ~Dft2d();

      void resize(unsigned int nx, unsigned int ny);
      void initialize(Image& image);
      void initializeForVis(Image& image);
      void initializeForVis(Angle& xSize, Angle& ySize, unsigned nXPix, unsigned nYPix);

      void setInput(Image& image);

      // True to normalize the DFT

      void normalize(bool norm);

      // Explicitly zero the data;

      void zero();

      // True to zeropad the DFT

      void zeropad(bool pad, unsigned fac=2);

      // Compute the transforms

      void computeForwardTransform(fftw_plan* fwdPlan=0);
      virtual void computeInverseTransform(fftw_plan* invPlan=0);

      // Return a pointer to the input data

      double* getInputData();

      // Return a pointer to the transformed data

      fftw_complex* getTransformDataPtr();
      unsigned getTransformLength();
      double*       getImageDataPtr();
      Image getImage(bool includeZeropaddedRegion=false);

      void removeMean();

      void complexMultiply(Dft2d& dft, bool conjugate=false);
      void checkConsistency(Dft2d& dft);

      void shift();
      void shiftBy(Angle& xoff, Angle& yoff);

      void plotInput(bool includeZeroPad=true);

      void plotData2D(DataType type);
      void plotReal();
      void plotImag();
      void plotAbs();
      void plotAbsLog();

      void plotRealU(unsigned iV);
      void plotRealV(unsigned iU);
      void plotImagU(unsigned iV);
      void plotImagV(unsigned iU);
      void plotAbsU(unsigned iV);
      void plotAbsV(unsigned iU);

      static double inverseGaussianWidth(double gaussianWidth);

      // Return an axis descriptor

      Axis& xAxis();
      Axis& yAxis();

      // Fill the transform array with a uniform plane

      void createUniformDft(unsigned nx, unsigned ny, double value);

      // Fill the transform array with a uniform disk

      void createUniformDiskDft(unsigned nx, unsigned ny, double radius);
      void createUniformDiskDft(double radius);

      // Fill the transform array with a uniform disk with a central obscuration

      void createBlockedApertureUniformDiskDft(unsigned nx, unsigned ny, double innerRadius, double outerRadius);
      void createBlockedApertureUniformDiskDft(double innerRadius, double outerRadius);

      // Fill the transform array with a 2-parameter tapered function, per Rohlfs & Wilson, 6.34

      void createTaperedApertureDft(unsigned nx, unsigned ny, Wavelength& lambda, Length& diameter, double K, unsigned p);
      void createTaperedApertureDft(Wavelength& lambda, Length& diameter, double K, unsigned p);

      // Fill the transform array with a tapered Bessel function

      void createJ0BesselFnDft(unsigned nx, unsigned ny, Wavelength& lambda, Length& diameter, double fractionalPower);
      void createJ0BesselFnDft(Wavelength& lambda, Length& diameter, double fractionalPower);

      // Fill the transform array with a blocked-aperture tapered Bessel function

      void createBlockedApertureJ0BesselFnDft(unsigned nx, unsigned ny, Wavelength& lambda, 
					      Length& innerDiameter, Length& outerDiameter,
					      double fractionalPower);
      void createBlockedApertureJ0BesselFnDft(Wavelength& lambda, Length& innerDiameter, Length& outerDiameter, double fractionalPower);

      // Calculate the J0 Bessel function

      static double j0BesselFn(double x);
      static double j0BesselFnRadius(double val, double eps = 1e-4);

      static double correlationLength(Length& diameter, Frequency& freq, double correlationCoeff);
      static double correlationLength(Length& diameter1, Length& diameter2, Frequency& freq, double correlationCoeff);


      void interpolateReImData(double u, double v, double& re, double& im, bool& valid);

      // Methods for manipulating an existing FT

      // Apply a gaussian taper to the low-frequency side of the
      // spectrum

      void highPass(double uvrPeak, double uvrSigma);

      // Apply a gaussian taper to the high-frequency side of the
      // spectrum

      void lowPass(double uvrPeak, double uvrSigma);

      // Apply a notch filter at the specified frequency

      void notch(double uvrCenter, double uvrSigma);

    public:

      static Mutex planGuard_;

      friend class Axis;

      // Parameters for convolving Fourier-plane data

      static const double convSigInPixels_;
      static const int    convMaskInPixels_;

      unsigned axes_;
      double* in_;            // The input array to be transformed
      fftw_complex* out_;     // The output of the transform

      fftw_plan forwardPlan_; // Instructions to fftw for the best forward method
      bool fwdPlanComputed_;

      fftw_plan inversePlan_; // Instructions to fftw for the best inverse method
			      // to FFT
      bool invPlanComputed_;

      bool precomputedPlan_;

      fftw_plan* fwdPlanTmp_;
      fftw_plan* invPlanTmp_;

      bool optimize_;         // True if we want fftw to perform expensive
			      // tests to determine the optimal FFT
			      // method.  Note that these tests will only
			      // be done on calls to resize()

      bool normalize_;
      unsigned zeropad_;

      unsigned nIn_;
      unsigned nOut_;
      unsigned nInZeroPad_;
      unsigned nOutZeroPad_;

      unsigned xOffset_;
      unsigned yOffset_;

      unsigned nx_;
      unsigned ny_;
      unsigned nxZeroPad_;
      unsigned nyZeroPad_;
      unsigned nxZeroPadPrev_;
      unsigned nyZeroPadPrev_;

      Axis xAxis_;
      Axis yAxis_;

      // True when this dft has been transformed

      bool isTransformed_;

      // Compute a plan for this fft

      void computePlan(fftw_plan& fwdPlan, fftw_plan& invPlan);
      void computePlan(unsigned flag); 
      void setPlan(fftw_plan forwardPlan, fftw_plan inversePlan);
      void normalizeTransform();
      void initialize();

      void getUData(unsigned iV, DataType type, std::vector<double>& xarr, std::vector<double>& yarr);
      void getVData(unsigned iU, DataType type, std::vector<double>& xarr, std::vector<double>& yarr);

      std::vector<double> getUVData(DataType type);

      void getUVData(unsigned dftInd, DataType type, double& u, double& v, double& val);
      void getUVData(unsigned iU, unsigned iV, DataType type, double& u, double& v, double& val);
      void getUVData(double u, double v, DataType type, double& uNearest, double& vNearest, double& val);
      void getUVData(int uInd, int iU, int vInd, int iV, double& re, double& im);

      // Return the uv radius of the index into the transform
      // array

      double uvRadius(unsigned dftIndex);
      double uvRadius(unsigned iU, unsigned iV);

      // Return the uv coordinate of the index into the transform
      // array

      void uvCoord(unsigned dftIndex, double& u, double& v);
      void uvCoord(unsigned iU, unsigned iV, double& u, double& v);

      // Return the u and v indices in the transform array
      // corresponding to this UV point

      void dftIndex(double u, double v, unsigned& dftInd, bool& conj);
      void uvIndex(double u, double v, unsigned& iU, unsigned& iV, bool& conj);
      void uvIndex(double u, double v, unsigned& iU, unsigned& iV, bool& conj, bool& valid);

      bool isZeroPad(unsigned iU, unsigned iV);
      virtual void resize();
      void updateAxis(unsigned axis);
      void checkAxes();
      void setPlan(Dft2d* dft=0);
      
      void setAngularSize(Angle size);
      void setSpatialFrequencyResolution(double inverseRadians);
      void setNpix(unsigned nPix);

      static unsigned nearestPowerOf2NotLessThan(double val);

    private:

      ImageAxis& xImageAxis();
      ImageAxis& yImageAxis();

    public:
      bool debugPrint_;

    }; // End class Dft2d

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_DFT2D_H
