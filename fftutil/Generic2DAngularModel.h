// $Id: Generic2DAngularModel.h,v 1.4 2012/05/16 18:00:50 eml Exp $

#ifndef GCP_UTIL_GENERIC2DANGULARMODEL_H
#define GCP_UTIL_GENERIC2DANGULARMODEL_H

/**
 * @file Generic2DAngularModel.h
 * 
 * Tagged: Fri May  4 10:35:22 PDT 2012
 * 
 * @version: $Revision: 1.4 $, $Date: 2012/05/16 18:00:50 $
 * 
 * @author Erik Leitch
 */
#include "gcp/fftutil/Dft2d.h"
#include "gcp/fftutil/Model.h"

#include "gcp/util/Area.h"
#include "gcp/util/Cosmology.h"
#include "gcp/util/Declination.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Length.h"
#include "gcp/util/SpectralType.h"
#include "gcp/util/SzCalculator.h"
#include "gcp/util/VariableUnitQuantity.h"

namespace gcp {
  namespace util {

    class UvDataGridder;

    //-----------------------------------------------------------------------
    // A generic base-class for 2D image util with angular coordinates
    //-----------------------------------------------------------------------

    class Generic2DAngularModel : public gcp::util::Model {
    public:
      
      //------------------------------------------------------------
      // A base-class for inheritors, to define data needed for model
      // evaluation in a multi-threaded context
      //------------------------------------------------------------

      class EvalData {
      public:

	EvalData() {
	};

	virtual ~EvalData() {
	};

      };

      //------------------------------------------------------------
      // A struct to define the per-thread data needed when executing
      // multi-threaded code within a thread pool
      //------------------------------------------------------------

      class ExecData {
      public: 

	Generic2DAngularModel* model_;
	Image*   image_;

	std::valarray<double>* xArr_;
	std::valarray<double>* yArr_;
	std::valarray<double>* dArr_;

	double angleConv_;
	unsigned iSegment_;
	unsigned nSegment_;
	unsigned iYStart_;
	unsigned iYStop_;
	double   prefactor_;
	unsigned type_;
	void*    params_;

	void*    evalData_;

	// Use for calculation

	unsigned nx_, ny_;
	double dxRad_, dyRad_;
	int xSense_, ySense_;
		    
	double cRotAng_, sRotAng_;
	double xOffRad_, yOffRad_;

	double x_, y_, xp_, yp_, xpp_, ypp_, arg_;
	int ix_, iy_;

	ExecData(Generic2DAngularModel* model) {
	  initialize();
	  model_    = model;

	  // And allocate eval data now

	  evalData_ = model->allocateEvalData();
	};

	void setTo(Image* image, 
		   unsigned iSegment, unsigned nSegment,
		   unsigned iYStart, unsigned iYStop) {
	  image_    = image;
	  iSegment_ = iSegment;
	  nSegment_ = nSegment;
	  iYStart_  = iYStart;
	  iYStop_   = iYStop;
	};

	void initialize()
	{
	  model_    = 0;
	  image_    = 0;
	  iSegment_ = 0;
	  nSegment_ = 0;
	  iYStart_  = 0;
	  iYStop_   = 0;
	  evalData_ = 0;
	  type_     = DataSetType::DATASET_UNKNOWN;
	  params_   = 0;
	};

	~ExecData() {
	  if(evalData_) {
	    //	    delete evalData_;
	    //	    evalData_ = 0;
	  }
	};

      };

      // Constructor.

      Generic2DAngularModel(ThreadPool* pool=0);
      
      // Destructor.

      virtual ~Generic2DAngularModel();

      // Initialize this object to sensible defaults

      void initialize();

      //------------------------------------------------------------
      // Properties generic to most model components
      //------------------------------------------------------------

      // Every 2D model can have an absolute position associated with
      // the center

      void setRa(gcp::util::HourAngle ra);
      void setDec(gcp::util::Declination dec);

      // Every 2D model can have an offset from center

      void setXOffset(gcp::util::Angle xOffset);
      void setYOffset(gcp::util::Angle yOffset);

      // Every 2D model can have a rotation angle (for non-symmetric
      // models)

      void setRotationAngle(gcp::util::Angle rotationAngle);
      
      // Every 2D model can have a normalization, possibly with units

      void setNormalization(double norm);

      // Every model can have a frequency at which that normalization
      // is defined

      void setNormalizationFrequency(gcp::util::Frequency normFreq);

      // Every 2D model can have a spectral index.  Interpretation of the
      // spectral index will depend on the type of normalization, ie.,
      // if normalization is in Temperature units, spectral index will
      // apply to temperature, if normalization is in Flux units,
      // spectral index will apply to flux).

      void setSpectralIndex(double alpha);

      //------------------------------------------------------------
      // Methods which should be instantiated by inheritors, where
      // appropriate
      //------------------------------------------------------------

      // Return the value to multiply the dimensionless model envelope
      // by to get physical units

      virtual double getEnvelopePrefactor(unsigned type, void* params);

      // Specific types of prefactors we know about

      double getRadioEnvelopePrefactor(Frequency* freq);
      double getSpectralIndexEnvelopePrefactor(Frequency* freq);
      double getSzEnvelopePrefactor(Frequency* freq);
      double getItohEnvelopePrefactor(Frequency* freq);

      virtual double pressureToComptonY() {
	return 1.0;
      }

      virtual double getXrayImageEnvelopePrefactor(Frequency* freq);

      // Return the value of the unity-normalized envelope at the
      // specified coordinate.  Ie, this is a shape function for this
      // model, modulo a normalization
      
      virtual double envelope(unsigned type, double xRad, double yRad);
      virtual double radioEnvelope(double xRad, double yRad);
      virtual double xrayImageEnvelope(double xRad, double yRad);
      virtual double genericEnvelope(double xRad, double yRad);

      virtual double envelope(void* execData, unsigned type, double xRad, double yRad);
      virtual double radioEnvelope(void* execData, double xRad, double yRad);
      virtual double xrayImageEnvelope(void* execData, double xRad, double yRad);
      virtual double genericEnvelope(void* execData, double xRad, double yRad);

      // Fill an image with externally specified parameters

      virtual void fillImage( unsigned type,  gcp::util::Image& image,          void* params=0);
      virtual void fillUvData(unsigned type, gcp::util::UvDataGridder& gridder, void* params=0);
      virtual void fillArray(unsigned type, Angle& axisUnits, 
			     std::valarray<double>& x, std::valarray<double>& y, std::valarray<double>& d, 
			     void* params=0);

      void fillImageSingleThread(unsigned type, gcp::util::Image& image, void* params);
      void fillImageMultiThread(ExecData* ed);
      unsigned initializeExecData(unsigned type, Image& image, void* params);
      static EXECUTE_FN(execFillImageMultiThread);

      void fillArraySingleThread(unsigned type, Angle& axisUnits, 
				 std::valarray<double>& x, std::valarray<double>& y, std::valarray<double>& d, 
				 void* params=0);
      void fillArrayMultiThread(ExecData* ed);
      unsigned initializeExecData(unsigned type, Angle& axisUnits, 
				  std::valarray<double>& x, std::valarray<double>& y, std::valarray<double>& d, 
				  void* params=0);
      static EXECUTE_FN(execFillArrayMultiThread);

      // Method to write fake 2D data from this model to a file

      void generateFake2DData(std::string fileName, Image& image, double sigma, std::string units);

      // Inherited method to initialize our thread pool

      void setThreadPool(ThreadPool* pool);

      // If inheritors support multi-threaded execution, they should
      // define this method to return an object containing the stack needed to
      // evaluate the envelope function

      virtual void* allocateEvalData();
      virtual void initializeEvalData(void* evalData);

      //------------------------------------------------------------
      // Store the absolute separation between an image and this model
      //------------------------------------------------------------

      void getAbsoluteSeparation(gcp::util::Image& image);
      void getAbsoluteSeparation(gcp::util::Dft2d& dft);

      //------------------------------------------------------------
      // Check our setup for sense
      //------------------------------------------------------------

      virtual void checkSetup();
      void checkSpectralSetup();
      void checkNormalizationSetup();
      bool normalizationImpliesSz();

      //------------------------------------------------------------
      // Overloaded sample() function from the Model base-class
      //------------------------------------------------------------

      virtual void sample();
      virtual void debugPrint();

      //------------------------------------------------------------
      // Integrate this model out to the specified radius
      //------------------------------------------------------------

      virtual gcp::util::SolidAngle solidAngleIntegral(unsigned type, Angle radius);
      virtual void solidAngleIntegral(unsigned type, Angle& radius, SolidAngle& ret);

      //------------------------------------------------------------
      // Calculate integrated Y out to the specified (physical) radius
      //------------------------------------------------------------

      gcp::util::Area integratedY(gcp::util::Cosmology& cosmo, gcp::util::Length& radius);
      void integratedY(double scaleFactor, Length& dA, Angle& theta, Area& ret);

      //------------------------------------------------------------
      // Get the conversion between native radio normalization units
      // and Compton-y
      //------------------------------------------------------------

      double getComptonYScaleFactor(Frequency& freq);
      double getComptonYScaleFactor();

      virtual void fillDerivedVariates() {};

      virtual PgModel pgModel();

      virtual void checkPosition();

    public:

      //------------------------------------------------------------
      // Put common model components here
      //------------------------------------------------------------

      // Any model component can have an RA/DEC

      gcp::util::HourAngle   ra_;
      bool hasRa_;

      gcp::util::Declination dec_;
      bool hasDec_;

      bool positionChecked_;
      bool hasAbsolutePosition_;

      // Used for calculation

      Angle xSep_;
      Angle ySep_;

      // Any model component can have an offset

      gcp::util::Angle xOffset_;
      gcp::util::Angle yOffset_;

      // Any model component can have a rotation angle

      gcp::util::Angle rotationAngle_;

      // Any model component can have an amplitude and associated
      // units.  We store this as a bare Variate, since we don't
      // automatically want to associate a physical quantity with it

      gcp::util::VariableUnitQuantity normalization_;
      gcp::util::VariableUnitQuantity radioNormalization_;
      gcp::util::VariableUnitQuantity xrayNormalization_;
      gcp::util::Unit::Units units_;

      // And a frequency at which the normalization is defined

      gcp::util::Frequency normalizationFrequency_;

      // And a spectral index

      gcp::util::VariableUnitQuantity spectralIndex_;

      // And a spectral type

      gcp::util::SpectralType spectralType_;

      // A gas-mass fraction

      gcp::util::VariableUnitQuantity fGas_;

      // The mean molecular weight

      gcp::util::VariableUnitQuantity mu_;

      // The mean molecular weight per free electron

      gcp::util::VariableUnitQuantity mue_;

      // The electronTemperature

      gcp::util::Temperature electronTemperature_;

    private:

      std::vector<ExecData*> execData_;

      SzCalculator szCalculator_;
      Temperature normTemperatureConv_;
      Intensity normIntensityConv_;

    }; // End class Generic2DAngularModel

  } // End namespace util
} // End namespace gcp

#endif // End #ifndef GCP_UTIL_GENERIC2DANGULARMODEL_H
