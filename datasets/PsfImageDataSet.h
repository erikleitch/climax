// $Id: $

#ifndef GCP_DATASETS_PSFIMAGEDATASET_H
#define GCP_DATASETS_PSFIMAGEDATASET_H

/**
 * @file PsfImageDataSet.h
 * 
 * Tagged: Wed Oct  2 13:54:33 PDT 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "datasets/DataSet2D.h"

#include "fftutil/Dft2d.h"
#include "fftutil/Image.h"
#include "fftutil/ImageManager.h"

#include "util/Angle.h"
#include "util/Distribution.h"

#define PSF_LK_FN(fn) double (fn)(void* obj, unsigned ind, unsigned& nDof)

namespace gcp {
  namespace datasets {

    class PsfImageDataSet : public DataSet2D {
    public:

      enum PsfType {
	PSF_UNKNOWN,
	PSF_NONE,
	PSF_REAL,
	PSF_GAUSS
      };

      //------------------------------------------------------------
      // Struct representing a single image from this dataset
      //------------------------------------------------------------

      struct PsfImageData {

	PsfImageData();
	virtual ~PsfImageData();

	//------------------------------------------------------------
	// Fitting interface
	//------------------------------------------------------------

	void addModel(gcp::util::Model& model, unsigned type);
	void clearModel();
	gcp::util::ChisqVariate computeChisq();
	gcp::util::ChisqVariate computeChisq2();

	//------------------------------------------------------------
	// Initialize this object from an image
	//------------------------------------------------------------

	void initialize(PsfImageDataSet* parent, gcp::util::Image& image, gcp::util::Image& preImage, gcp::util::Antenna& ant, gcp::util::Frequency& freq);
	void initializeFromImage(gcp::util::Image& imageToInit, gcp::util::Image& image);
	void initializeErrors(gcp::util::Image& image, gcp::util::Image& preImage);
	
	void assignPaddedImageIfNecessary(gcp::util::Image& image);

	void display();
	void displayCompositeModel();
	void displayResiduals();

	void simulateData(double sigma);
	void convolveModel();

	static PSF_LK_FN(fixedErrorGauss);
	static PSF_LK_FN(fixedErrorGaussApprox);
	static PSF_LK_FN(fixedErrorPoiss);
	static PSF_LK_FN(pixelErrorGauss);
	static PSF_LK_FN(pixelErrorPoiss);

	static PSF_LK_FN(fixedErrorPoissTest);

	double fixedErrorGauss(unsigned ind, unsigned& nDof);
	double fixedErrorGaussApprox(unsigned ind, unsigned& nDof);
	double fixedErrorPoiss(unsigned ind, unsigned& nDof);
	double pixelErrorGauss(unsigned ind, unsigned& nDof);
	double pixelErrorPoiss(unsigned ind, unsigned& nDof);

	double fixedErrorPoissTest(unsigned ind, unsigned& nDof);

	//------------------------------------------------------------
	// A function for evaluating the likelihood
	//------------------------------------------------------------

	PSF_LK_FN(*lkFn_);
	PSF_LK_FN(*lkFn2_);

	//------------------------------------------------------------
	// The parent of this object
	//------------------------------------------------------------

	PsfImageDataSet* parent_;

	//------------------------------------------------------------
	// Frequency and antenna descriptors
	//------------------------------------------------------------

	gcp::util::Frequency frequency_; // The frequency
	gcp::util::Antenna   antenna_;   // The antenna descriptor 

	//------------------------------------------------------------
	// A composite image-plane model for this frequency and Stokes
	// parameter
	//------------------------------------------------------------

	gcp::util::Image compositeModel_;

	//------------------------------------------------------------
	// A work array to save object creation and destruction during
	// construction of the composite Image-plane model
	//------------------------------------------------------------

	gcp::util::Image modelComponent_;

	//------------------------------------------------------------
	// A container for holding the transform of the composite
	// Image-plane model
	//------------------------------------------------------------

	gcp::util::Dft2d compositeModelDft_;

	//------------------------------------------------------------
	// A primary beam for this image, and its integral (pbSum_)
	//------------------------------------------------------------

	gcp::util::Image primaryBeam_;
	double pbSum_; 

	//------------------------------------------------------------
	// A container for holding the transform of the primary beam
	//------------------------------------------------------------

	gcp::util::Dft2d primaryBeamDft_;

	//------------------------------------------------------------
	// The data for this image
	//------------------------------------------------------------

	gcp::util::Image data_;

	//------------------------------------------------------------
	// The portion of the image that contains data
	//------------------------------------------------------------

	unsigned iXStart_;
	unsigned iYStart_;
	unsigned nX_;
	unsigned nY_;

	//------------------------------------------------------------
	// We can potentially have a separate error for each pixel
	//------------------------------------------------------------
	
	gcp::util::Image error_;
	
	//------------------------------------------------------------
	// Or a single error, either specified or estimated from the
	// data using thetaMinErr_
	//------------------------------------------------------------
	
	bool hasErrorVal_;
	bool hasBackground_;
	bool hasNoiseRms_;
	double background_;
	double noiseRms_;
	gcp::util::Angle thetaMinErr_;

	//------------------------------------------------------------
	// Members for excluding data from a certain region
	//------------------------------------------------------------
	
	gcp::util::Angle thetaMin_;
	gcp::util::Angle thetaMax_;
	bool excMin_;
	bool excMax_;

	//------------------------------------------------------------
	// PSF and Distribution types
	//------------------------------------------------------------

	PsfType psfType_;
	gcp::util::Distribution::Type distType_;

	//------------------------------------------------------------
	// The current model offset, in radians
	//------------------------------------------------------------

	double currentXoffRad_;
	double currentYoffRad_;
      };

      /**
       * Constructor.
       */
      PsfImageDataSet();

      /**
       * Destructor.
       */
      virtual ~PsfImageDataSet();

      //------------------------------------------------------------
      // Inherited dataset interface for Markov chains
      //------------------------------------------------------------

      void addModel(gcp::util::Model& model);
      void clearModel();
      gcp::util::ChisqVariate computeChisq();
      gcp::util::ChisqVariate computeChisq2();
      virtual void finalizeForDisplay();

      //------------------------------------------------------------
      // Utility methods for computeChisq()
      //------------------------------------------------------------

      void convolveModels();
      gcp::util::ChisqVariate accumulateChisq();
      gcp::util::ChisqVariate accumulateChisq2();

      //------------------------------------------------------------
      // Display methods
      //------------------------------------------------------------

      void setupForDisplay();
      virtual void display();
      virtual void displayCompositeModel();
      virtual void displayResiduals();

      //------------------------------------------------------------
      // Initialization methods
      //------------------------------------------------------------

      virtual void loadData(bool simulate);
      void loadDataFromFile();
      virtual void initializeForSim();
      virtual void initializeFromFitsFile(std::string fileName);
      void printFileStats();

      void initializePositionDependentData();
      void initFromImage(gcp::util::Image& image, gcp::util::Image& preImage);

      virtual void initImage(std::string fileName, gcp::util::Image& image) {};

      //------------------------------------------------------------
      // Simulate methods
      //------------------------------------------------------------

      virtual void simulateData(double sigma);
      void writeCompositeModelToFile(std::string fileName, double sigma);

      void setName(std::string name);

      //------------------------------------------------------------
      // Utility methods
      //------------------------------------------------------------
      
      static PsfType psfType(std::string type);
      static gcp::util::Distribution::Type distType(std::string type);
      static std::string distToString(gcp::util::Distribution::Type type);

      std::string displayHeaderString(std::string type);

    public:

      //------------------------------------------------------------
      // Data will be stored as an array of images
      //------------------------------------------------------------

      std::vector<PsfImageData> data_;
      PsfType psfType_;

      gcp::util::ImageManager imp_;

    }; // End class PsfImageDataSet

  } // End namespace datasets
} // End namespace gcp



#endif // End #ifndef GCP_DATASETS_PSFIMAGEDATASET_H
