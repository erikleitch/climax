// $Id: VisDataSet.h,v 1.11 2012/05/31 16:22:56 eml Exp $

#ifndef GCP_DATASETS_VISDATASET_H
#define GCP_DATASETS_VISDATASET_H

/**
 * @file VisDataSet.h
 * 
 * Tagged: Wed Jun 16 10:24:00 PDT 2010
 * 
 * @version: $Revision: 1.11 $, $Date: 2012/05/31 16:22:56 $
 * 
 * @author Erik Leitch
 */
#include "gcp/datasets/DataSetManager.h"
#include "gcp/datasets/DataSet2D.h"

#include "gcp/fftutil/Antenna.h"
#include "gcp/fftutil/FitsIoHandler.h"
#include "gcp/fftutil/Generic2DAngularModel.h"
#include "gcp/fftutil/Image.h"
#include "gcp/fftutil/ObsInfo.h"
#include "gcp/fftutil/Stokes.h"
#include "gcp/fftutil/UvDataGridder.h"

#include "gcp/util/Angle.h"
#include "gcp/util/BitMask.h"
#include "gcp/util/CondVar.h"
#include "gcp/util/ChisqVariate.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Declination.h"
#include "gcp/util/FitsUvfReader.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/Geoid.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Length.h"
#include "gcp/util/Mutex.h"
#include "gcp/util/Percent.h"
#include "gcp/util/SolidAngle.h"
#include "gcp/util/ThreadPool.h"

#include <map>
#include <vector>

#define ACC_DISPATCH_FN(fn) void (fn)(VisDataSet::VisData& data, bool& first, VisDataSet& src, VisDataSet& dest, \
				      VisDataSet::VisBaselineGroup& group, VisDataSet::VisStokesData& stokes, VisDataSet::VisFreqData& freq, gcp::util::ObsInfo::Vis& vis)

namespace gcp {
  namespace datasets {

    //------------------------------------------------------------
    // We inherit from DataSetManager (instead of DataSet2D) so that
    // subdatasets (ie, stacked) can be acccessed from the calling
    // interface
    //------------------------------------------------------------

    class VisDataSet : public gcp::datasets::DataSetManager {

    public:

      class VisBaselineGroup;
      class VisStokesData;
      class VisFreqData;

      //------------------------------------------------------------
      // To make maps of what data are we accumulating visibility data?
      //------------------------------------------------------------

      enum AccumulatorType {
	ACC_DATA,       // Image of the data
	ACC_MODEL,      // Image of the synthesized-beam filtered composite model
	ACC_RES,        // Image of the data residuals after the model has been subtracted
	ACC_BEAM,       // Image of the synthesized beam
	ACC_CLEAN,      // Image of the model, convolved with a Gaussian
	                // approximation to the synthesized beam, and with noise residuals added
	ACC_CLEANMODEL, // Accumulate the model only, filled aperture (no synthesized beam)
	ACC_CLEANBEAM,  // Accumulate the approximate synthesized beam
      };

      //------------------------------------------------------------
      // What type of image are we making?
      //------------------------------------------------------------

      enum ImageType {
	IMG_SKY,  // Beam-corrected mosaicked map (sky signal) over frequencies + antenna pairs
	IMG_SNR,  // (Beam-corrected) Significance map over frequencies + antenna pairs
	IMG_DIRTY // Combined dirty map -- no mosaicking of frequencies + antenna pairs
      };

      //=======================================================================
      // Struct defining an antenna associated with this visibility
      // data set
      //=======================================================================

      struct VisAntenna {
	gcp::util::Antenna* antenna_; // The antenna descriptor for this antenna
	bool isUnique_;               // True if this is a unique type of antenna
      };

      //=======================================================================
      // Struct for defining a single baseline
      //=======================================================================

      struct VisBaseline {
	int aipsBaselineIndex_;
	gcp::util::Antenna* ant1_;
	gcp::util::Antenna* ant2_;

	// For simulation only

	VisBaseline() {
	  aipsBaselineIndex_ = -1;
	  ant1_  = 0;
	  ant2_  = 0;
	};

	VisBaseline(unsigned aipsBaselineIndex, gcp::util::Antenna& ant1, gcp::util::Antenna& ant2) {
	  aipsBaselineIndex_ = aipsBaselineIndex;
	  ant1_  = &ant1;
	  ant2_  = &ant2;
	}
      };

      //=======================================================================
      // A struct for encapsulating data from a single baseline, at a
      // single stokes and frequency, for a single timestamp
      //=======================================================================

      struct VisData {
	double u_;
	double v_;
	double w_;

	double r_;

	double re_;
	double im_;
	double wt_;

	unsigned int baseline_;
	double jd_;
      };

      //=======================================================================
      // Convenience struct for calling execute methods
      //=======================================================================

      struct VisExecData {
	static gcp::util::Mutex dataAccessGuard_;
	VisDataSet* vds_;
	VisFreqData* vfd_;
	gcp::util::Antenna* ant1_;
	gcp::util::Antenna* ant2_;
	unsigned iGroup_;
	unsigned iStokes_;
	unsigned iFreq_;
	gcp::util::ChisqVariate* chisq_;
	gcp::util::Generic2DAngularModel* model_;

	VisExecData(VisDataSet* vds, VisFreqData* vfd, unsigned iGroup, unsigned iStokes, unsigned iFreq) {
	  initialize(vds, vfd, iGroup, iStokes, iFreq);
	};

	void initialize(VisDataSet* vds, VisFreqData* vfd, unsigned iGroup, unsigned iStokes, unsigned iFreq) {
	  vds_     = vds;
	  vfd_     = vfd;
	  iGroup_  = iGroup;
	  iStokes_ = iStokes;
	  iFreq_   = iFreq;
	  ant1_    = 0;
	  ant2_    = 0;
	  chisq_   = 0;
	  model_   = 0;
	}


	VisExecData(VisDataSet* vds, VisFreqData* vfd, gcp::util::Antenna* ant1, gcp::util::Antenna* ant2, 
		    unsigned iGroup, unsigned iStokes, unsigned iFreq) 
	{
	  initialize(vds, vfd, iGroup, iStokes, iFreq);
	  initialize(ant1, ant2);
	}

	void initialize(gcp::util::Antenna* ant1, gcp::util::Antenna* ant2) 
	{
	  ant1_ = ant1;
	  ant2_ = ant2;
	}

	VisExecData(VisDataSet* vds, VisFreqData* vfd, gcp::util::ChisqVariate* chisq, 
		    unsigned iGroup, unsigned iStokes, unsigned iFreq) {
	  initialize(vds, vfd, iGroup, iStokes, iFreq);
	  initialize(chisq);
	};

	void initialize(gcp::util::ChisqVariate* chisq) {
	  chisq_ = chisq;
	}

	VisExecData(VisDataSet* vds, VisFreqData* vfd, gcp::util::Generic2DAngularModel* model,
		    unsigned iGroup, unsigned iStokes, unsigned iFreq) {
	  initialize(vds, vfd, iGroup, iStokes, iFreq);
	  initialize(model);
	}

	void initialize(gcp::util::Generic2DAngularModel* model) {
	  model_ = model;
	}

      };

      //=======================================================================
      // For a given type of baseline at a single Stokes parameter and
      // frequency, this struct encapsulates all timestamps and
      // baselines.
      //=======================================================================

      struct VisFreqData {

	//------------------------------------------------------------
	// Image-plane model handling
	//------------------------------------------------------------

	// A composite image-plane model for this frequency and Stokes parameter

	gcp::util::Image compositeImageModel_;

	// A work array to save object creation and destruction during
	// construction of the composite Image-plane model

	gcp::util::Image imageModelComponent_;

	// A container for holding the transform of the composite
	// Image-plane model

	gcp::util::Dft2d compositeImageModelDft_;

	//------------------------------------------------------------
	// Fourier-plane model handling
	//------------------------------------------------------------

	// A work array for use during construction of the composite
	// Fourier-plane model

	gcp::util::UvDataGridder fourierModelComponent_;

	// A container for holding the composite Fourier-plane model

	gcp::util::UvDataGridder compositeFourierModelDft_;

	//------------------------------------------------------------
	// Data handling
	//------------------------------------------------------------

	// A container for holding the gridded visibility data

	gcp::util::UvDataGridder griddedData_;

	// A container for performing temporary operations, like
	// gridding residuals for a single frequency

	gcp::util::UvDataGridder utilityGridder_;

	// A primary beam for this frequency and Stokes parameter

	gcp::util::Image primaryBeam_;

	// The frequency of this data set

	gcp::util::Frequency frequency_;

	// The IF number of this frequency (will be used to optionally
	// include/exclude this IF)

	unsigned ifNo_;

	// An index that will be used to label this frequency across
	// all groups and stokes parameters

	unsigned globalIfIndex_;

	// The bandwidth

	gcp::util::Frequency bandwidth_;

	// Estimated primary beam halfwidth

	gcp::util::Angle primaryBeamHalfWidth_;

	// The maximum UV radius of any visibility in this group, in
	// inverse radians
	
	double uvrMax_;
	double uAbsMax_;
	double vAbsMax_;
	double reMean_;
	double imMean_;
	double estChisq_;
	double wtScale_;
	double wtSumTotal_;
	
	std::vector<double> wtSums_;
	std::vector<gcp::util::Angle> xShifts_;
	std::vector<gcp::util::Angle> yShifts_;

	// The number of visibilities in this group

	unsigned nVis_;

	// Temporary variable used by external classes

	unsigned iVis_;

	// The number of unflagged visibilities in this group

	unsigned nVisUsed_;

	// The Stokes container to which this object belongs

	VisStokesData* stokes_;

	// The group to which this object belongs

	VisBaselineGroup* group_;

	//------------------------------------------------------------
	// Simulation only
	//------------------------------------------------------------

	bool hasImage_;
	bool generatingFakeData_;
	
	//------------------------------------------------------------
	// For multithreading
	//------------------------------------------------------------

	VisExecData* execData_;

	//------------------------------------------------------------
	// For convenience, a copy of the combined synthesized beam of
	// this dataset
	//------------------------------------------------------------

	gcp::util::SolidAngle estimatedGlobalSynthesizedBeam_;

	//------------------------------------------------------------
	// And estimated synthesized beam widths for this VisFreqData
	// subset
	//------------------------------------------------------------

	gcp::util::Angle synthBeamMajSig_;
	gcp::util::Angle synthBeamMinSig_;
	gcp::util::Angle synthBeamRotAngle_;

	bool debug_;

	//------------------------------------------------------------
	// General methods
	//------------------------------------------------------------

	bool hasData() {
	  return generatingFakeData_ || nVisUsed_ > 0;
	}
	
	void operator=(const VisFreqData& data) {
	  *this = (VisFreqData&) data;
	}

	void operator=(VisFreqData& data);

	void mergeData(VisDataSet::VisFreqData& freq, 
		       gcp::util::Angle& xShift, gcp::util::Angle& yShift);

	void duplicate(VisFreqData& data);
	
	VisFreqData(const VisFreqData& data) {
	  *this = data;
	}

	VisFreqData(VisFreqData& data) {
	  *this = data;
	}

	VisFreqData() {
	  iVis_       = 0;
	  nVis_       = 0;
	  nVisUsed_   = 0;
	  uvrMax_     = 0.0;
	  uAbsMax_    = 0.0;
	  vAbsMax_    = 0.0;
	  reMean_     = 0.0;
	  imMean_     = 0.0;
	  estChisq_   = 0.0;
	  wtScale_    = 1.0;
	  wtSumTotal_ = 0.0;
	  debug_      = false;
	  group_      = 0;
	  stokes_     = 0;

	  // Simulation only

	  hasImage_           = false;
	  generatingFakeData_ = false;

	  execData_ = 0;
	}

	virtual ~VisFreqData() {
	  if(execData_) {
	    delete execData_;
	    execData_ = 0;
	  }
	}

	bool operator==(VisFreqData& freq);

	std::string formatString();

	// Add a model component to this data set

	void addModel(gcp::util::Generic2DAngularModel& model);
	void remModel();
	void addImagePlaneModel(gcp::util::Generic2DAngularModel& model);
	void addFourierPlaneModel(gcp::util::Generic2DAngularModel& model);

	// Clear all model components

	void clearModel();

	// Resize for a new correlation length

	void resize(double percentCorrelation, double correlationLength, bool isSim);

	// Resize to match an image

	void resize(gcp::util::Image& image, bool isSim, gcp::util::UvDataGridder** planGridder);

	bool isImagePlaneModel(gcp::util::Generic2DAngularModel& model);

	// Take a composite image-plane model and transform, prior to
	// computing chisq

	void transformModel();

	// Compute the chisq

	virtual gcp::util::ChisqVariate computeChisq();

	void accumulateMoments(bool first, VisDataSet::VisData& data, gcp::util::Angle& xShift, gcp::util::Angle& yShift);
	void accumulateVarianceStats(VisDataSet::VisData& data);

	void storeWtSum(gcp::util::Angle& xShift, gcp::util::Angle& yShift);

	//------------------------------------------------------------
	// Accumulate data, model, or residuals
	//------------------------------------------------------------

	void accumulate(gcp::util::UvDataGridder& gridder, VisDataSet::AccumulatorType type);
	void accumulateDirty(gcp::util::UvDataGridder& gridder, VisDataSet::AccumulatorType type);
	void accumulateClean(gcp::util::UvDataGridder& gridder, VisDataSet::AccumulatorType type);
	gcp::util::Image getCleanImage(VisDataSet::AccumulatorType type);

	//------------------------------------------------------------
	// Simulation only
	//------------------------------------------------------------

	void addImage(gcp::util::Image& image, gcp::util::UvDataGridder** planGridder);

	// Return true if an image has been installed for this
	// Frequency

	bool hasImage();

	// Take a composite image-plane model, multiply by the primary
	// beam, and transform

	void transformImage();

	// Take the gridded Fourier-plane data, and transform it to
	// the image plane

	void inverseTransformData();

	friend std::ostream& operator<<(std::ostream& os, VisFreqData& freq);
	friend std::ostream& operator<<(std::ostream& os, const VisFreqData& freq);
      };

      //=======================================================================
      // For a given type of baseline at a single Stokes parameter,
      // this struct encapsulates all frequencies, timestamps and
      // baselines.
      //=======================================================================

      struct VisStokesData {
       
	// Which Stokes parameter does this data represent?
       
	gcp::util::Stokes::Param stokes_;
       
	// A vector of data, sorted by Frequency, for this Stokes parameter
       
	std::vector<VisFreqData> freqData_;

	// An index that will be used to label this Stokes parameter across
	// all groups

	unsigned globalStokesIndex_;

	double uvrMax_;

	void checkFrequencyIndex(unsigned iFreq);

	VisStokesData() {
	  uvrMax_ = 0.0;
	}

	bool operator==(VisStokesData& stokes);

	friend std::ostream& operator<<(std::ostream& os, const VisStokesData& stokes);

	// Simulation only: Return true if an image has been installed
	// for all frequencies of this Stokes parameter

	bool hasImage();
      };

      //=======================================================================
      // For a given type of baseline, this struct encapsulates all
      // stokes parameters, frequencies, timestamps and baselines.
      //=======================================================================
     
      struct VisBaselineGroup {

	// How many baselines (per timestamp) are represented in this
	// group?

	unsigned nBaseline_;
	unsigned iBaseline_;

	// The two types of antennas that comprise this baseline group

	std::pair<gcp::util::Antenna, gcp::util::Antenna> antennaPair_;

	// The vector of Stokes parameters present

	std::vector<VisStokesData> stokesData_;

	friend std::ostream& operator<<(std::ostream& os, const VisBaselineGroup& group);

	double uvrMax_;

	VisDataSet* dataset_;

	//------------------------------------------------------------
	// Simulation only
	//------------------------------------------------------------

	std::vector<VisBaseline> baselines_;

	//------------------------------------------------------------
	// General methods
	//------------------------------------------------------------

	VisBaselineGroup() {
	  nBaseline_ = 0;
	  uvrMax_    = 0.0;
	  dataset_   = 0;
	}

	~VisBaselineGroup() {
	}

	// Initialize pertinent members of this struct 

	void initialize(unsigned nStokes, std::vector<gcp::util::Frequency>& frequencies, 
			std::vector<gcp::util::Frequency>& bandwidths, bool debug);

	void initialize(unsigned nStokes, unsigned nFrequency, unsigned nVis, bool debug);

	void initialize(std::map<gcp::util::Stokes::Param, std::map<double, VisDataSet::VisFreqData*> >& stokesMap);

	// Install a primary beam for this group

	void installPrimaryBeam(gcp::util::Image& image, gcp::util::Stokes stokes, unsigned iFreq);

	void checkStokesIndex(unsigned iStokes);

	void addBaseline(unsigned iBase, gcp::util::Antenna& ant1, gcp::util::Antenna& ant2);

	gcp::util::Angle estimatePrimaryBeamHalfWidth(gcp::util::Frequency& freq);

	//------------------------------------------------------------
	// Simulation only
	//------------------------------------------------------------

	void calculateVisibilities(unsigned& iVisGroup,
				   gcp::util::Declination& dec, 
				   gcp::util::HourAngle& ha); 

	//------------------------------------------------------------
	// An equivalency operator
	//------------------------------------------------------------

	bool operator==(VisBaselineGroup& group);
      };

      //=======================================================================
      // Methods of VisDataSet
      //=======================================================================

      //------------------------------------------------------------
      // Used for loading data into this object from another
      // VisDataSet. This map maintains the relationship between VisFreqData in
      // the external object and the corresponding VisFreqData object
      // in this one
      //------------------------------------------------------------

      std::map<VisFreqData*, VisFreqData*> freqMap_;

      static const double   fourierPlaneConvSig_;
      static const unsigned fourierPlaneNConv_;

      /**
       * Constructor.
       */
      VisDataSet(gcp::util::ThreadPool* pool=0);

      /**
       * Destructor.
       */
      virtual ~VisDataSet();

      // Get the observation info container for this data set

      gcp::util::ObsInfo& obs();

      // Initialize this data set from a file.

      void initializeFromFile(std::string fileName);
      void countData(std::string fileName);

      void initializeAndCountData(bool simulate);
      void initializeAndCountDataSingle(std::string fileName);
      void initializeAndCountDataMultiple();

      void initializeStoreParameters();
      void storeDataIfRequested(VisData& data, VisFreqData* destFreq, VisDataSet& src, gcp::util::ObsInfo::Vis& vis);
      void buildAntennaMap(std::vector<VisDataSet*>& datasets);
      void remapStokesAndFreqs();
      void initializeStoreArrays(std::vector<VisDataSet*>& datasets);
      void updateAntennaMap(VisDataSet* dataset);

      void getStoreIndices(VisData& data, VisFreqData* destFreq, VisDataSet& src, unsigned& iGroup, unsigned& iVis, unsigned& iDate);
      unsigned getAipsBaselineIndex(VisDataSet& dataset, unsigned aipsIndex);

      void determineUvMax(gcp::util::ObsInfo& obs);

      //------------------------------------------------------------
      // Data loading
      //------------------------------------------------------------

      // Inherited interface to load data into this dataset

      void loadData(bool simulate);

      // Wrapper around single/multiple data loading

      void loadData();

      // Load data from a single file

      void loadDataSingle(std::string fileName);

      // Load data from multiple datasets into this object

      void loadDataMultiple();

      // Load data, checking that the chisq matches the weights

      void loadDataWithChecks();

      // Return true if data should be reloaded, based on chisq
      // estimation

      bool reload();
      bool reloadSingle(std::string namePrefix, std::string file);
      bool reloadMultiple();
      
      void parseParameters();

      void printReImRms();
      void printReIm();

      double estimateWtScale();
      void printOccupiedIndices();

      void storeWtSums(gcp::util::Angle& xShift, gcp::util::Angle& yShift);

      void loadDataFromObs();

      void setupForSimulation(bool sim);
      void storeDataInternallyOnReadin(bool store);
      void releaseDataAfterReadin(bool release);

      // Compute the chi-squared of this dataset with the composite image-plane model

      virtual gcp::util::ChisqVariate computeChisq();

      // Install frequencies

      void installFrequencies(std::vector<gcp::util::Frequency>& freqs);
      void installBandwidth(gcp::util::Frequency bw);

      // Compute primary beams

      void computePrimaryBeams();

      // Shift data if requested

      void shiftIfRequested();

      // Compute the combined estimated synthesized beam

      void computeGlobalSynthesizedBeam();
      void insertSynthesizedBeamModelForPlots(gcp::util::PgModelManager& pgManager);
      void storeSynthesizedBeamModelForPlots(gcp::util::PgModelManager& pgManager);

      static gcp::util::Model* addDeltaFunctionModel(std::vector<gcp::util::Model*>& models, 
						     double flux, gcp::util::Angle& xOff, gcp::util::Angle& yOff);

      static std::vector<gcp::util::Image::Window> parseWindows(std::string windowSpec, gcp::util::Image& image);

      static std::string parseFileName(std::string fileStr, bool& shiftRequested, gcp::util::Angle& xoff, gcp::util::Angle& yoff);

      void addSynthesizedBeamDisplayModel(gcp::util::Image& image);

      static void addDisplayWindow(gcp::util::Image::Window& win);

      // Compute the estimated synthesized beam for each VisFreqData
      // subset

      void estimateSynthesizedBeams();

      // Plot methods

      void plotUv(  int iGroup=-1, int iStokes=-1, int iFreq=-1);
      void plotReal(int iGroup=-1, int iStokes=-1, int iFreq=-1);
      void plotImag(int iGroup=-1, int iStokes=-1, int iFreq=-1);
      void plotAbs( int iGroup=-1, int iStokes=-1, int iFreq=-1);
      void plotSimVis();

      void addModel(gcp::util::Model& model);
      void remModel();
      void clearModel();

      void displayPrimaryBeams();

      //------------------------------------------------------------
      // Methods used for simulation
      //------------------------------------------------------------

      // Add an image to the image-plane model for this frequency

      void addImage(gcp::util::Image& image, int iFrequency=-1, int iStokes=-1);

      void calculateSimulatedUvw();

      // Observe an image, overwriting internally stored visibilities

      void observe();

      // Observe an image, using the supplied observation object as
      // the

      void observe(gcp::util::ObsInfo& obs);

      void transformImages();

      void replaceVisibilities();
      void replaceVisibilities(VisDataSet::AccumulatorType type);

      void calculateVisibilities(gcp::util::ObsInfo& obs);
      void fillSimulationVisibilityArray(gcp::util::ObsInfo& obs);

      void estimateErrorInMeanFromData(bool estimate);

      void writeUvfFile(std::string fileName);
      void writeDataToFile(std::string fileName, VisDataSet::AccumulatorType type);

      void inverseTransformData();

      // Get the requested image of any type

      gcp::util::Image getImage(AccumulatorType accType, ImageType imgType);

      // Get the requested dirty image (average over frequency/antenna pairs)

      gcp::util::Image getImage(AccumulatorType type);

      // Get the requested mosaicked image

      void getImage(gcp::util::Image& image, gcp::util::Image& wtimage, AccumulatorType type);

      gcp::util::Image getCleanImage(VisDataSet::AccumulatorType type);

      // Accumulate a mosaicked image

      void accumulateBeamCorrectedImage(gcp::util::Image& image, gcp::util::Image& wtimage);

      virtual void displayIfRequested();
      virtual void initializeForDisplay();
      void uvPlot();
      void radPlot();

      static void setWedge(VisDataSet::AccumulatorType type, 
			   double zmin, double zmax, bool specified);

      void display();
      void displayBeam();
      void displayCompositeModel();
      void clean();

      void writeData();

      void displayResiduals();
      void display(VisDataSet::AccumulatorType type);

      // Inherited interface from DataSet

      virtual void simulateData(double sigma);
      virtual void writeCompositeModelToFile(std::string fileName, double sigma);

      virtual void setParameter(std::string name, std::string val, std::string units=" ");
      virtual void incrementParameter(std::string name, std::string val, std::string units=" ");
      
    protected:

      // A bit mask to keep track of which groups and frequencies are
      // done for a given operation

      gcp::util::BitMask doneMask_;

      // A condition variable we will block on to wait for doneness,
      // when executing in multi-thread mode

      gcp::util::CondVar doneVar_;

      // If true, we will estimte the error in the mean of each
      // gridded Fourier cell from the data themselves

      bool estimateErrInMeanFromData_;

      // If true, data will be stored internally during read-in

      bool storeDataInternally_;
      bool releaseDataAfterReadin_;

      // The vector of baseline groups encountered in this data set

      std::vector<VisBaselineGroup> baselineGroups_;
      
      // Maps used to convert between AIPS-style baseline indices, and
      // internal baseline group indices

      std::map<unsigned, int> baselineTagToGroupIndexMap_;
      std::map<unsigned, int> aipsBaselineIndexToGroupIndexMap_;
      std::map<unsigned, unsigned> aipsBaselineIndexCount_;

      // True if a weight is valid for a visibility

      bool goodWt(double wt);

      void guessAtAntennaType();
      void printFileStats(std::string fileName);

      double taper(double u, double v);

      double wtScale_;
      double intScale_;
      bool   autoScale_;

      //------------------------------------------------------------
      // Parameters for gridding the data to a fixed correlation
      // percentage
      //------------------------------------------------------------

      bool usePerc_;
      double percentCorrelation_;

      //------------------------------------------------------------
      // If usePerc_ == false, we use this to
      // store the image to which we will grid the data
      //------------------------------------------------------------

      gcp::util::Image image_; 

      //-----------------------------------------------------------------------
      // A gridder for storing the data from frequencies, all baseline
      // groups, for display.  if we ever actually deal with real
      // Stokes data, this will have to be per-Stokes parameter, but
      // for now we just lump everything together
      //-----------------------------------------------------------------------

      gcp::util::UvDataGridder utilityGridder_;

      int lastGroup_;
      int lastDate_;
      VisDataSet* lastDataSet_;
      double lastJd_;
      unsigned lastBaseline_;
      double writeJd_;

      // A map of antenna number in each dataset to antenna number in
      // our array

      std::map<VisDataSet*, std::map<unsigned, unsigned> > datasetAntMap_;

      //------------------------------------------------------------
      // The min/max UV data to grid on read-in
      //------------------------------------------------------------

      double uvMin_;
      double uvMax_;

      //------------------------------------------------------------
      // The absolute min/max of all good data
      //------------------------------------------------------------

      double uAbsMax_;
      double vAbsMax_;

      std::map<gcp::util::Antenna::AntennaType, gcp::util::Antenna::AntennaType> includedAntTypes_;
      std::map<unsigned, unsigned> excludedAntNos_;
      std::map<unsigned, unsigned> excludedIfNos_;
      std::map<unsigned, unsigned> includedIfNos_;

      bool   taper_;
      bool   taperInvert_;
      double taperSigma_;

      bool reverseDelays_;
      bool reverseDisplay_;
      bool forceWt_;

      double wtSumTotal_;

      double wtMin_;
      double wtMax_;

      bool shiftRequested_;

      //------------------------------------------------------------
      // If multiple files were specified, we are stacking data into
      // this dataset.  In this case, we will instantiate a vector of
      // datasets used for loading in data
      //------------------------------------------------------------

      std::vector<std::string> fileList_;
      std::vector<VisDataSet*> datasets_;

      // A map of all unique Stokes parameters we encountered

      std::map<gcp::util::Stokes::Param, unsigned> allStokes_;

      // A map of all unique Frequencies we encountered
      
      std::map<double, unsigned> allFreqs_;
      std::map<double, double>   allBws_;

      //------------------------------------------------------------
      // The maximum PB halfwidth of any data maintained by this
      // object
      //------------------------------------------------------------

      gcp::util::Angle maxPrimaryBeamHalfwidth_;

      gcp::util::Angle synthBeamMajSig_;
      gcp::util::Angle synthBeamMinSig_;
      gcp::util::Angle synthBeamRotAngle_;

      bool visibilitiesInitialized_;
      
    public:

      gcp::util::Angle xShift_;
      gcp::util::Angle yShift_;

      void determineUniqueBaselineGroupings(gcp::util::ObsInfo& obs);

      //------------------------------------------------------------
      // Return an estimate of the largest primary beam in this data
      // set
      //------------------------------------------------------------

      gcp::util::Angle estimateLargestPrimaryBeamFwhm();

      static std::string typeString(VisDataSet::AccumulatorType type);
      std::string displayHeaderString(VisDataSet::AccumulatorType type);

      void operator+=(VisDataSet& vds);

    public:

      void checkAntenna(gcp::util::Antenna& ant, std::vector<gcp::util::Antenna>& uniqueAnts);
      void checkBaselineGroupIndex(unsigned iGroup);

      void getVis(gcp::util::Dft2d::DataType type,
		  std::vector<float>& x, std::vector<float>& y, 
		  int iGroup, int iStokes, int iFreq);

      // Initialize dfts, if gridding the visibility data

      void initializeVisibilityArrays(double percentCorrelation);
      void initializeVisibilityArrays2(double percentCorrelation);
      void initializeVisibilityArrays(gcp::util::Image& image);
      void initializeEstChisq();
      void storeEstChisqAsWtScale();

      void purgeZeroSizedGroups();
      void purgeZeroSizedGroupsTest();

      // Transform image-plane model for computing chisq

      void transformModels();
      gcp::util::ChisqVariate accumulateChisq();

      void accumulateMoments(std::string fileName, bool first);
      void accumulateMoments(std::string fileName, bool first, bool init, VisDataSet& vds, ACC_DISPATCH_FN(*dispatchFn));

      void initializeForMomentAccumulation(bool first);
      void calculateErrorInMean();

      static ACC_DISPATCH_FN(dispatchInternal);
      static ACC_DISPATCH_FN(dispatchExternal);

      VisBaselineGroup& findMatch(VisBaselineGroup& group);

      //------------------------------------------------------------
      // Virtual functions that must be defined by inheritors
      //------------------------------------------------------------

      virtual void openFileReader(std::string fileName)                    = 0;
      virtual void closeFileReader()                                       = 0;
      virtual void updateFrequencyInformation()                            = 0;
      virtual void updateObservationInformation()                          = 0;
      virtual void updateObservationInformation(VisDataSet* dataset)       = 0;
      virtual void updateVisibilityInformation()                           = 0;
      virtual void getGroup(unsigned iGroup, gcp::util::ObsInfo::Vis& vis) = 0;
      virtual void initializeAntennaInformation(std::string fileName)      = 0;
      virtual void initializeFrequencyInformation(std::string fileName)    = 0;

      bool haveImages();
      bool canSimulate();

      // Accumulate the requested data into internal gridders

      void accumulate(VisDataSet::AccumulatorType type);

      // Accumulate the requested data into the passed gridder

      void accumulate(VisDataSet::AccumulatorType type, gcp::util::UvDataGridder& gridder);

      // Initialize global gridders to match the passed gridder

      void initializeGlobalGridders(gcp::util::UvDataGridder* gridder);

      gcp::util::Image getRequiredImage(double percentCorrelation);

      //-----------------------------------------------------------------------
      // Inherited method to initialize position dependent data for
      // this dataset
      //-----------------------------------------------------------------------

      void initializePositionDependentData();

      //-----------------------------------------------------------------------
      // Method to set a list of antenna types to include on read-in,
      // and to set a list of antenna numbers excluded
      //-----------------------------------------------------------------------

      void initializeIncludedAntennaTypes(std::string includedTypes);
      void initializeExcludedAntennaNumbers(std::string excludedNos);

      void initializeIncludedIfNumbers(std::string includedNos);
      void initializeExcludedIfNumbers(std::string excludedNos);

      void initializeTaper(std::string taperSpec);

      //-----------------------------------------------------------------------
      // True if this group involves included antennas
      //-----------------------------------------------------------------------

      bool antennaIsIncluded(gcp::util::Antenna& ant);
      bool antennaTypeIsIncluded(gcp::util::Antenna& ant);
      bool antennaNumberIsIncluded(gcp::util::Antenna& ant);
      bool isIncluded(VisBaselineGroup& group);

      bool ifNumberIsIncluded(unsigned ifNo);
      bool ifNumberIsExcluded(unsigned ifNo);
      bool useIfNumber(unsigned ifNo);

      void debugPrint();

      void getAntOs(std::ostream& os, std::vector<gcp::util::Antenna>& uniqueAnts);
      void getBaseOs(std::ostream& os);
      void getIfOs(std::ostream& os);

      //-----------------------------------------------------------------------
      // Private methods to do with multi-threading
      //-----------------------------------------------------------------------

    private:

      // Multi-thread-aware version of addModel

      void addModelMultiThread(VisFreqData& vfd, gcp::util::Generic2DAngularModel& model, 
			       unsigned iGroup, unsigned iStokes, unsigned iFreq);
      static EXECUTE_FN(execAddModel);

      // Multi-threaded version of computeChisq

      void computeChisqMultiThread(VisFreqData& vfd, gcp::util::ChisqVariate& chisq,
				   unsigned iGroup, unsigned iStokes, unsigned iFreq);
      static EXECUTE_FN(execComputeChisq);

      // Multi-threaded version of transformImage

      void transformImageMultiThread(VisFreqData& vfd, unsigned iGroup, unsigned iStokes, unsigned iFreq);
      static EXECUTE_FN(execTransformImage);

      // Multi-threaded version of transformModel

      void transformModelMultiThread(VisFreqData& vfd, unsigned iGroup, unsigned iStokes, unsigned iFreq);
      static EXECUTE_FN(execTransformModel);

      // Multi-threaded version of computePrimaryBeam

      void computePrimaryBeamMultiThread(VisFreqData& vfd, gcp::util::Antenna& ant1, gcp::util::Antenna& ant2, 
					 unsigned iGroup, unsigned iStokes, unsigned iFreq);
      static EXECUTE_FN(execComputePrimaryBeam);

      void registerDone(unsigned iGroup, unsigned iStokes, unsigned iFreq);
      void registerPending(unsigned iGroup, unsigned iStokes, unsigned iFreq);
      void waitUntilDone();
      void initWait();

      void mergeData(VisDataSet& vds, VisBaselineGroup& findGroup, VisStokesData& findStokes, VisFreqData& findFreq);

    protected:

      //------------------------------------------------------------
      // Methods used to load data from an external VisDataSet
      //------------------------------------------------------------

      void buildInternalMap(std::vector<VisDataSet*>& datasets);
      void initializeBaselineGroups(std::vector<VisDataSet*>& datasets);
      void mapFrequencies(std::vector<VisDataSet*>& datasets);
      VisFreqData& findMatch(VisBaselineGroup& group, VisStokesData& stokes, VisFreqData& freq);

      virtual void addDataSet(std::string file);
      DataSet* getDataSet(std::string name);


    }; // End class VisDataSet

    std::ostream& operator<<(std::ostream& os, VisDataSet::VisFreqData& freq);
    std::ostream& operator<<(std::ostream& os, const VisDataSet::VisFreqData& freq);
    std::ostream& operator<<(std::ostream& os, const VisDataSet::VisStokesData& stokes);
    std::ostream& operator<<(std::ostream& os, const VisDataSet::VisBaselineGroup& group);

  } // End namespace datasets
} // End namespace gcp



#endif // End #ifndef GCP_DATASETS_VISDATASET_H
