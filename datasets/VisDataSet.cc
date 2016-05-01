#include "gcp/datasets/VisDataSet.h"

#include "gcp/fftutil/Dft2d.h"
#include "gcp/fftutil/FitsIoHandler.h"

#include "gcp/models/PtSrcModel.h"
#include "gcp/models/Generic2DGaussian.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/util/Date.h"
#include "gcp/util/Delay.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Exception.h"
#include "gcp/util/FitsBinTableReader.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/RangeParser.h"
#include "gcp/util/Sampler.h"
#include "gcp/util/String.h"
#include "gcp/util/Timer.h"
#include "gcp/util/Wavelength.h"

#include "cpgplot.h"

using namespace gcp::datasets;
using namespace gcp::util;
using namespace std;

#if 1
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif


//#define TIMER_TEST
#define SIM_TIMER_TEST

#ifdef TIMER_TEST
Timer addmodeltimer1, addmodeltimer2, addmodeltimer3;
double amt1 =0.0, amt2= 0.0, amt3 = 0.0;
Timer cc1, cc2, cc3;
double cct1=0.0, cct2=0.0, cct3=0.0;
#endif

Mutex VisDataSet::VisExecData::dataAccessGuard_;

//=======================================================================
// Methods of VisDataSet
//=======================================================================

/**.......................................................................
 * Constructor.
 */
VisDataSet::VisDataSet(gcp::util::ThreadPool* pool) 
{
  pool_                      = pool;
  storeDataInternally_       = false;
  releaseDataAfterReadin_    = true;
  estimateErrInMeanFromData_ = false;
  uvMin_                     = -1;
  uvMax_                     = -1;
  uAbsMax_                   = 0.0;
  vAbsMax_                   = 0.0;
  obsWasSet_                 = false;

  wtScale_                   = 1.0;
  intScale_                  = 1.0;
  autoScale_                 = true;

  taper_                     = false;
  taperSigma_                = 0.0;
  taperInvert_               = false;

  reverseDelays_             = false;
  reverseDisplay_            = false;

  usePerc_                   = true;
  percentCorrelation_        = 0.95;
  forceWt_                   = false;

  wtSumTotal_                = 0.0;

  maxPrimaryBeamHalfwidth_.setRadians(0.0);

  shiftRequested_            = false;
  xShift_.setRadians(0.0);
  yShift_.setRadians(0.0);

  datasets_.resize(0);

  addParameter("perc",           DataType::DOUBLE,  "If assigned, only grid together data that are at least this correlated"); 
  addParameter("npix",           DataType::UINT,    "The number of pixels to use when gridding the data");			   
  addParameter("size",           DataType::DOUBLE,  "The size of the image to which the data will be gridded");		   
  addParameter("wtscale",        DataType::STRING,  "A factor by which to scale the weights in the UVF file (or 'auto' to estimate from the data)");			   
  addParameter("intscale",       DataType::DOUBLE,  "A factor by which to scale the integration time when generating simulated data");			   
  addParameter("uvmin",          DataType::DOUBLE,  "The minimum UV radius to use");		   
  addParameter("uvmax",          DataType::DOUBLE,  "The maximum UV radius to use");                   
  addParameter("store",          DataType::BOOL,    "True to store data internally on readin");
  addParameter("useanttypes",    DataType::STRING,  "List of antenna types to include on read-in.  Use like: 'useanttypes = sza,bima'");
  addParameter("excludeants",    DataType::STRING,  "List of antenna numbers to exclude.  Use like: 'excludeants = 1,3,5'");
  addParameter("excludeifs",     DataType::STRING,  "List of IF numbers to exclude.  Use like: 'excludeifs = 1-3,5'");
  addParameter("includeifs",     DataType::STRING,  "List of IF numbers to include.  Use like: 'includeifs = 1-3,5'");
  addParameter("uvtaper",        DataType::STRING,  "UV-taper to apply to the weights.  Use like: 'uvtaper val, uvrad', where uvrad is in lambda");
  addParameter("power",          DataType::DOUBLE,  "Point on the beam to which images should be displayed (defaults to 0.5)");
  addParameter("wtmin",          DataType::DOUBLE,  "Minimum valid weight for the combined map (used for mosaicking only)");
  addParameter("antplot",        DataType::BOOL,    "True to plot antennas on read-in");
  addParameter("reversedelays",  DataType::BOOL,    "Flip the sense of the delays on read-in");
  addParameter("reversedisplay", DataType::BOOL,    "Flip the sense of the x-display");
  addParameter("forcewt",        DataType::BOOL,    "If true, force weights to match the data variance (wtscale must be set to 'auto')");
  addParameter("displaybeam",    DataType::BOOL,    "If true, display the synthesized beam on read-in");
  addParameter("clean",          DataType::BOOL,    "If true, a clean model image will be displayed.  Note: if clean = true, you must force the data to a fixed resolution using parameters 'size' and 'npix' (i.e., you cannot use the 'perc' option, which is the default)");
  addParameter("cleantype",      DataType::STRING,  "If clean=true, the type of clean image to display, one of: 'model' or 'delta'.  If cleantype = model, then any currently defined models will be used to construct the clean image.  If cleantype = delta, then a delta-function model will be iteratively constructed using the (Hogbom) CLEAN algorithm.");
  addParameter("cleaniter",      DataType::UINT,    "If clean=true and cleantype=delta, the number of clean iterations to perform (default is 100)");
  addParameter("cleangain",      DataType::DOUBLE,  "If clean=true and cleantype=delta, the clean gain to use when subtracting model components (default is 0.05)");
  addParameter("cleancutoff",    DataType::DOUBLE,  "If clean=true and cleantype=delta, then cleaning will stop when the maximum absolute residual reaches this value (Jy).  Default is no cutoff (0.0)");
  addParameter("cleanwindow",    DataType::STRING,  "If clean=true and cleantype=delta, and no model is specified, then this sets a clean window to be searched to iteratively build up a model.  Format is: 'cleanwindow = [xmin:xmax, ymin:ymax, units]', or 'cleanwindow = [xpos +- xdelta, ypos +- ydelta, units]', where 'units' are the units in which the window is specified.  Alternately, 'cleanwindow = [tag +- delta, units]' can be used, where tag is one of: 'min', 'max' or 'abs' to set up a clean window about the minimum, maximum or absolute maximum in the image. Use 'cleanwindow += [xmin:xmax, ymin:ymax, units]' to specify more than one clean window.");
  addParameter("datauvf",        DataType::STRING,  "UVF file to output data visibilities");  
  addParameter("resuvf",         DataType::STRING,  "UVF file to output residual visibilities");  
  addParameter("modeluvf",       DataType::STRING,  "UVF file to output model visibilities");  

  remParameter("file");
  addParameter("file",        DataType::STRING, "The input file for this dataset.  An optional shift can also be specified.  I.e., 'file=name, shift=0.01,0.01 deg;' would cause the dataset to be shifted by 0.01 degree in x and y.");

  //------------------------------------------------------------
  // Set defaults for some parameters
  //------------------------------------------------------------

  setParameter("cleanwindow", "[abs +- 0.025, abs +- 0.025, deg]");
  setParameter("cleaniter",   "100");
  setParameter("cleangain",   "0.05");
  setParameter("cleantype",   "model");
  setParameter("cleancutoff", "0.0");

  dataSetType_ = DataSetType::DATASET_2D | DataSetType::DATASET_RADIO;
}

/**.......................................................................
 * Destructor.
 */
VisDataSet::~VisDataSet() 
{
  for(unsigned i=0; i < datasets_.size(); i++) {
    delete datasets_[i];
    datasets_[i] = 0;
  }
}

/**.......................................................................
 * Configure for use as a simulator.
 */
void VisDataSet::setupForSimulation(bool sim)
{
  if(sim) {
    storeDataInternallyOnReadin(true);
    releaseDataAfterReadin(false);
  } else {
    releaseDataAfterReadin(true);
  }
}

/**.......................................................................
 * If true, store data internally on readin.  We may either want to do
 * this because it is much faster when calculating moments to iterate
 * over an internal copy of the data than to re-read the data from a
 * file, or if simulating, so we have the original data to overwrite
 * with simulated visibilities.
 */
void VisDataSet::storeDataInternallyOnReadin(bool store)
{
  storeDataInternally_ = store;
}

/**.......................................................................
 * If true, release the internal copy of the data after moment
 * calculation is complete.  If using for simulation, we don't want to
 * do this!
 */
void VisDataSet::releaseDataAfterReadin(bool release)
{
  releaseDataAfterReadin_ = release;
}

/**.......................................................................
 * Initialize this object from a file
 */
void VisDataSet::initializeFromFile(std::string fileName)
{
  //------------------------------------------------------------
  // Initialize antenna information
  //------------------------------------------------------------

  initializeAntennaInformation(fileName);

  //------------------------------------------------------------
  // Initialize frequency information
  //------------------------------------------------------------

  initializeFrequencyInformation(fileName);

  //------------------------------------------------------------
  // Initialize relevant header information from the file
  //------------------------------------------------------------

  openFileReader(fileName);

  //------------------------------------------------------------
  // Update internal frequency members from what we just determined
  //------------------------------------------------------------

  updateFrequencyInformation();

  //------------------------------------------------------------
  // Update internal observation members from what we just determined
  //------------------------------------------------------------

  updateObservationInformation();

  //------------------------------------------------------------
  // Update internal information about the visibility data
  // from what we just determined
  //------------------------------------------------------------

  updateVisibilityInformation();

  //------------------------------------------------------------
  // Close the file reader
  //------------------------------------------------------------

  closeFileReader();

  //------------------------------------------------------------
  // Now print a summary of what we know, before reading the
  // visibility data
  //------------------------------------------------------------

  printFileStats(fileName);

  //------------------------------------------------------------
  // And guess at the type of antennas
  //------------------------------------------------------------

  guessAtAntennaType();
}

/**.......................................................................
 * Count the data by baseline
 */
void VisDataSet::countData(std::string fileName)
{
  unsigned ngoodvis = 0;
  double lastMjd;

  //------------------------------------------------------------
  // If we were told to store data internally, resize our internal
  // array of groups to accommodate the data.
  //------------------------------------------------------------

  if(storeDataInternally_) {
    obs_.visibilities_.resize(obs_.nGroup_);
  }

  //------------------------------------------------------------
  // Now read through the file
  //------------------------------------------------------------

  unsigned badBaseline = 0;
  unsigned badVis      = 0;
  unsigned nBadGroup   = 0;

  //------------------------------------------------------------
  // First examine what we know about the antennas to determine how
  // many unique baseline groupings there are.  
  //------------------------------------------------------------

  determineUniqueBaselineGroupings(obs_);

  //------------------------------------------------------------
  // Now read through the UVF file, counting how many baselines in
  // each group were actually encountered, and storing the maximum UV
  // radius, which will be used to determine how to compress the data,
  // if compression is requested
  //------------------------------------------------------------

  openFileReader(fileName);

  ObsInfo::Vis vis;
  Wavelength wavelength;
  TimeVal lastTs, currTs, delta;
  double deltaSec, avIntTime=0.0;
  unsigned nInt=0;

  for(unsigned iGroup=0; iGroup < obs_.nGroup_; iGroup++) {

    getGroup(iGroup, vis);
    
    //------------------------------------------------------------
    // Estimate the integration time from the timestamps
    //------------------------------------------------------------

    currTs.setMjd(vis.jd_ - 2400000.5);

    if(iGroup > 0) {
      try {
	delta = currTs - lastTs;
	deltaSec = delta.getTimeInSeconds();

	if(deltaSec > 0.0 && deltaSec < 100.0)
	  avIntTime += (deltaSec - avIntTime) / (nInt + 1);
	  
      } catch(...) {
      }
    }

    lastTs = currTs;

    //------------------------------------------------------------
    // If we are storing data internally, save a copy in the
    // visibility array
    //------------------------------------------------------------

    if(storeDataInternally_) {
      obs_.visibilities_[iGroup] = vis;
    }

    //------------------------------------------------------------
    // Find the baseline grouping to which this baseline belongs
    //------------------------------------------------------------

    int baselineGroupIndex = aipsBaselineIndexToGroupIndexMap_[vis.baseline_];

    //------------------------------------------------------------
    // Determine whether this visibility is valid
    //------------------------------------------------------------

    bool goodBaseline = false;
    if(baselineGroupIndex >= 0)
      goodBaseline = isfinite(vis.u_) && isfinite(vis.v_) && isfinite(vis.w_);
    else {
      goodBaseline = false;
      ++nBadGroup;
    }

    if(goodBaseline) {

      double mjd = vis.jd_ - 2400000.5;

      if(ngoodvis == 0) {
	obs_.mjdMin_ = mjd;
	obs_.mjdMax_ = mjd;
      } else {
	obs_.mjdMin_ = mjd < obs_.mjdMin_ ? mjd : obs_.mjdMin_;
	obs_.mjdMax_ = mjd > obs_.mjdMax_ ? mjd : obs_.mjdMax_;
      }

      VisBaselineGroup& group = baselineGroups_[baselineGroupIndex];

      //------------------------------------------------------------
      // Iterate over Stokes parameters and IFs to determine the maximum
      // UV radius
      //------------------------------------------------------------
      
      Delay delay;
      VisData data;
      unsigned iVis=0;
      unsigned nVis=0;

      for(unsigned iIf=0; iIf < obs_.nFreq_; iIf++) {
	for(unsigned iStokes=0; iStokes < obs_.nStokes_; iStokes++) {
	  
	  VisStokesData& stokes = group.stokesData_[iStokes];
	  VisFreqData& freq     = stokes.freqData_[iIf];
	  wavelength.setFrequency(obs_.frequencies_[iIf]);

	  delay.setDelayInSeconds(vis.u_);
	  data.u_  = delay.cm() / wavelength.cm();
	  delay.setDelayInSeconds(vis.v_);
	  data.v_  = delay.cm() / wavelength.cm();
	  delay.setDelayInSeconds(vis.w_);
	  data.w_  = delay.cm() / wavelength.cm();

	  data.r_ = sqrt(data.u_ * data.u_ + data.v_ * data.v_);
	  
	  data.re_ = vis.re_[iVis];
	  data.im_ = vis.im_[iVis];

	  //------------------------------------------------------------
	  // By convention, we assume that the weights represent the
	  // Re/Im variance, or half the variance on the total
	  // intensity
	  //
	  // We also apply any user-supplied scaling and taper at this
	  // point
	  //------------------------------------------------------------

	  data.wt_ = vis.wt_[iVis] * wtScale_ * taper(data.u_, data.v_);

	  iVis++;

	  //------------------------------------------------------------
	  // Only count data that are finite, have finite and non-zero
	  // weights, and fall within the specified UV min/max
	  //------------------------------------------------------------

	  bool goodVis = isfinite(data.re_) && isfinite(data.im_) && isfinite(data.wt_) && data.wt_ > 0.0;

	  goodVis &= (uvMin_ < 0.0 || data.r_ > uvMin_);
	  goodVis &= (uvMax_ < 0.0 || data.r_ < uvMax_);

	  if(debug_)
	    COUT("Counting " << baselineGroupIndex << " iif = " << iIf << " re = " << data.re_ << " im = " << data.im_ << " wt = " << data.wt_ << " u = " << data.u_ << " v = " << data.v_ << " goodVis = " << goodVis << " baseline = " << vis.baseline_);


	  if(goodVis) {

	    //------------------------------------------------------------
	    // Store the max UV radius of the visibilities for this
	    // frequency
	    //------------------------------------------------------------
	    
	    double uAbs = fabs(data.u_);
	    double vAbs = fabs(data.v_);

	    freq.uAbsMax_  = freq.uAbsMax_  > uAbs ? freq.uAbsMax_     : uAbs;
	    freq.vAbsMax_  = freq.vAbsMax_  > vAbs ? freq.vAbsMax_     : vAbs;

	    uAbsMax_  = uAbsMax_  > uAbs ? uAbsMax_ : uAbs;
	    vAbsMax_  = vAbsMax_  > vAbs ? vAbsMax_ : vAbs;
	    
	    freq.uvrMax_   = freq.uvrMax_   > data.r_ ? freq.uvrMax_   : data.r_;
	    stokes.uvrMax_ = stokes.uvrMax_ > data.r_ ? stokes.uvrMax_ : data.r_;
	    group.uvrMax_  = group.uvrMax_  > data.r_ ? group.uvrMax_  : data.r_;

	    //------------------------------------------------------------
	    // And increment the number of visibilities in this
	    // VisFreqData container
	    //------------------------------------------------------------

	    freq.nVis_++;
	    nVis++;

	    //------------------------------------------------------------
	    // Accumulate the means and weights for this VisFreqData
	    // object.  These will be used to check any user-specified
	    // weight scale against the weights implied by the data
	    //------------------------------------------------------------

	    freq.reMean_ += (data.re_ - freq.reMean_) / (freq.nVis_);
	    freq.imMean_ += (data.im_ - freq.imMean_) / (freq.nVis_);

	    if(ngoodvis == 0) {
	      wtMin_ = data.wt_;
	      wtMax_ = data.wt_;
	    } else {
	      wtMin_ = data.wt_ < wtMin_ ? data.wt_ : wtMin_;
	      wtMax_ = data.wt_ > wtMax_ ? data.wt_ : wtMax_;
	    }

	    ngoodvis++;
	    wtSumTotal_ += data.wt_;

	  } else {
	    ++badVis;
	  } // End if(goodVis)

	} // End loop over Stokes parameters
      } // End loop over IFs

      //------------------------------------------------------------
      // And increment the count of this baseline group, but only if
      // some unflagged data were found for this group
      //------------------------------------------------------------

      if(nVis != 0)
	group.nBaseline_++;
      else
	++nBadGroup;

      //------------------------------------------------------------
      // If this baseline was bas (or excluded), add it to the count
      // of unusable visibilities
      //------------------------------------------------------------

    } else {
      badVis += (obs_.nStokes_ * obs_.nFreq_);
    }

    if(iGroup % 1000 == 0) {
      std::cout << "\rReading..." << (100*iGroup)/obs_.nGroup_ << "%";
      fflush(stdout);
    }

  } // End loop over groups in file

  //------------------------------------------------------------
  // Finally, resize the internal bitmask to accommodate all axes of
  // this data set
  //------------------------------------------------------------

  synchronizer_.resize(obs_.getNumberOfFrequencies() * obs_.getNumberOfStokesParameters() * baselineGroups_.size());

  //------------------------------------------------------------
  // Close the file and report statistics
  //------------------------------------------------------------

  closeFileReader();

  std::cout << "\rReading...100%\r                   ";

  COUTCOLOR(std::endl << "UT ranged from " << Date::mjdToCal(obs_.mjdMin_) << " to " << Date::mjdToCal(obs_.mjdMax_)
	    << " (MJD " << obs_.mjdMin_ << " to " << obs_.mjdMax_ << ")", "cyan");

  COUTCOLOR(std::endl << "Estimated integration time per visibility: " << avIntTime << "s", "cyan");
  COUTCOLOR(             "Estimated total integration time:          " << setprecision(1) << std::fixed << ((obs_.nGroup_/obs_.nBaseline_) * avIntTime)/3600 << "h", "cyan");

  COUTCOLOR(std::endl << "Weights ranged from " << wtMin_ << " to " << wtMax_ << " (rms = " << sqrt(1.0/wtMax_) << " to " << sqrt(1.0/wtMin_) << " Jy)", "cyan");

  COUTCOLOR(String::formatHumanReadableInteger(badVis) << " visibilities were flagged as unusable ("
	    << std::setprecision(4) << (double)(badVis)/(obs_.nGroup_ * obs_.nStokes_ * obs_.nFreq_) * 100 <<
	    "% of the total)", "cyan");

  COUTCOLOR(String::formatHumanReadableInteger(nBadGroup) << " groups were flagged as unusable", "cyan");
}

/**.......................................................................
 * Load data from multiple datasets into this one
 */
void VisDataSet::loadDataMultiple()
{
  //------------------------------------------------------------
  // Now initialize gridders and images to match the requested percent
  // correlation or the requested image scale
  //------------------------------------------------------------

  if(usePerc_) {
    initializeVisibilityArrays(percentCorrelation_);
  } else {
    initializeVisibilityArrays(image_);
  }

  //------------------------------------------------------------
  // Now read the data from each dataset in turn, dispatching the
  // resulting data into this object's arrays
  //------------------------------------------------------------

  initializeForMomentAccumulation(true);

  for(unsigned i=0; i < datasets_.size(); i++) {
    VisDataSet* dataset = datasets_[i];

    //------------------------------------------------------------
    // Accumulate data for this dataset into our gridders
    //------------------------------------------------------------

    dataset->accumulateMoments(fileList_[i], true, true, *this, VisDataSet::dispatchExternal);

    //------------------------------------------------------------
    // And store the weight sums associated with this dataset -- these
    // will be used later to construct an appropriate primary beam
    //------------------------------------------------------------
    
    storeWtSums(dataset->xShift_, dataset->yShift_);
  }
  
  if(storeDataInternally_) {
    obs_.visibilities_.resize(lastGroup_+1);
    obs_.setNumberOfTimestamps(lastGroup_+1);
  }

  calculateErrorInMean();
}

void VisDataSet::addDataSet(std::string file)
{
  ThrowError("Inheritor has not defined this method");
}

/**.......................................................................
 * Parse any parameters needed to read the data in
 */
void VisDataSet::parseParameters()
{
  //------------------------------------------------------------
  // Check whether we are gridding to match a target correlation
  // percentage, or an image scale.  If neither is specified, default
  // to correlation percentage
  //------------------------------------------------------------

  //------------------------------------------------------------
  // Sanity checks for gridding.  First check if a percentage was
  // specified.
  //------------------------------------------------------------

  if(getParameter("perc", false)->data_.hasValue()) {

    if(getParameter("npix", false)->data_.hasValue() || getParameter("size", false)->data_.hasValue())
      ThrowColorError("You can specify an image size (" << name_ << ".size) and resolution ('" << name_ << ".npix) "
		      "or a correlation percentage (" << name_ << ".perc), but not both", "red");

    percentCorrelation_ = getDoubleVal("perc");
    usePerc_ = true;
    
    //------------------------------------------------------------
    // If not, if either a size or resolution was specified, make sure
    // both were
    //------------------------------------------------------------

  } else if(getParameter("npix", false)->data_.hasValue() || getParameter("size", false)->data_.hasValue()) {

    if(!getParameter("npix", false)->data_.hasValue() || !getParameter("size", false)->data_.hasValue())
      ThrowSimpleColorError("You must specify both the image size (" << name_ << ".size) and resolution (" << name_ <<".npix)", "red");

    unsigned npix;
    npix = getUintVal("npix");

    Angle size;
    size.setName(name_, "size");
    size.setVal(getDoubleVal("size"), getParameter("size", true)->units_);

    image_.xAxis().setNpix(npix);
    image_.xAxis().setAngularSize(size);
    
    image_.yAxis().setNpix(npix);
    image_.yAxis().setAngularSize(size);
    usePerc_ = false;
  }

  //------------------------------------------------------------
  // Sanity check for plotting
  //------------------------------------------------------------

  if(getParameter("clean", false)->data_.hasValue() && getBoolVal("clean")) {
    if(!getParameter("npix", false)->data_.hasValue() || !getParameter("size", false)->data_.hasValue())
      ThrowSimpleColorError("If displaying CLEAN images (" << name_ << ".clean = true), you must specify an image size (" 
			    << name_ << ".size) and resolution (" << name_ <<".npix) to force the data to a fixed resolution", "red");
  }

  //------------------------------------------------------------
  // Check for other externally-specified variables
  //------------------------------------------------------------

  if(getParameter("uvmin", false)->data_.hasValue())
    uvMin_ = getDoubleVal("uvmin");

  if(getParameter("uvmax", false)->data_.hasValue())
    uvMax_ = getDoubleVal("uvmax");

  if(getParameter("intscale", false)->data_.hasValue())
    intScale_ = getDoubleVal("intscale");

  if(getParameter("wtscale", false)->data_.hasValue()) {
    String wtScaleStr(getStringVal("wtscale"));

    if(wtScaleStr.contains("auto")) {
      wtScale_   = 1.0;
      autoScale_ = true;
    } else {
      wtScale_ = wtScaleStr.toDouble();
      autoScale_ = false;
    }
  }

  if(getParameter("forcewt", false)->data_.hasValue()) {
    forceWt_ = getBoolVal("forcewt");
  }

  if(getParameter("store", false)->data_.hasValue()) {
    bool store = getBoolVal("store");
    storeDataInternallyOnReadin(store);
    releaseDataAfterReadin(!store);
  }

  // Override the user's choice if we were requested to output data

  if(getParameter("datauvf", false)->data_.hasValue() ||
     getParameter("modeluvf", false)->data_.hasValue() ||
     getParameter("resuvf", false)->data_.hasValue()) 
  {
    storeDataInternallyOnReadin(true);
    releaseDataAfterReadin(true);
  }
     
  if(getParameter("useanttypes", false)->data_.hasValue()) {
    initializeIncludedAntennaTypes(getStringVal("useanttypes"));
  }

  if(getParameter("excludeants", false)->data_.hasValue()) {
    initializeExcludedAntennaNumbers(getStringVal("excludeants"));
  }

  if(getParameter("excludeifs", false)->data_.hasValue()) {
    initializeExcludedIfNumbers(getStringVal("excludeifs"));
  }

  if(getParameter("includeifs", false)->data_.hasValue()) {
    initializeIncludedIfNumbers(getStringVal("includeifs"));
  }

  if(getParameter("uvtaper", false)->data_.hasValue()) {
    initializeTaper(getStringVal("uvtaper"));
  }
}

/**.......................................................................
 * Load data, checking the weights for sanity.  If weights don't match
 * the variance in the data, this function may re-load the data with
 * weight scaling.
 */
void VisDataSet::loadDataWithChecks()
{
  try {
    loadData();
  } catch(...) {
    return;
  }

  if(reload())
    loadData();
}

/**.......................................................................
 * Load data from one or more data sets
 */
void VisDataSet::loadData()
{
  loadDataMultiple();
}

/**.......................................................................
 * Return true if data should be reloaded
 */
bool VisDataSet::reload()
{
  return reloadMultiple();
}

/**.......................................................................
 * Print warnings and return true if the data should be reloaded with
 * different scaling
 */
bool VisDataSet::reloadSingle(std::string namePrefix, std::string file)
{
  std::string name;
  if(namePrefix == "")
    name = name_;
  else {
    std::ostringstream nameStr;
    nameStr << namePrefix << "." << name_;
    name = nameStr.str();
  }

  //------------------------------------------------------------
  // Check the weights, but only if no taper was requested --
  // otherwise we can't reliably estimate the variance from the
  // data and estimated or calculated weight scalings won't make
  // any sense
  //------------------------------------------------------------

  if(!taper_) {

    double wtScale = estimateWtScale();

    if(autoScale_) {

      wtScale_ = wtScale;

      // Zero estimated chisq 

      storeEstChisqAsWtScale();
      initializeEstChisq();

      if(abs(1.0 - wtScale_) > 0.05 || forceWt_) {

	if(forceWt_) {
	  COUTCOLOR(std::endl << "Warning: forcing the weights to match the rms estimated from the data." << std::endl
		    << "You can use '" << name << ".forcewt = false' to disable this behavior. " << std::endl
		    << "But your confidence intervals won't make sense if the reported errors don't match the actual data rms " 
		    << std::endl << " (reduced chi-squared will not be close to 1)"
		    , "yellow");
	} else {
	  COUTCOLOR(std::endl << "Warning: scaling the weights from reported values by " << std::setprecision(2) << wtScale_ 
		    << ", as estimated from the data." << std::endl
		    << "You can use '" << name << ".wtscale = 1' to disable this behavior. " << std::endl 
		    << "But your confidence intervals won't make sense if the reported errors don't match the actual data rms " 
		    << std::endl << " (reduced chi-squared will not be close to 1)"
		    , "yellow");
	}

	if(namePrefix != "")
	  COUTCOLOR("(While loading data from file: " << file << ")", "yellow");

	return true;

      }

    } else {

      if(abs(1.0 - wtScale) > 0.05) {
	COUTCOLOR(std::endl << "Warning: you are scaling weights by " << std::setprecision(2) << wtScale_ 
		  << " but I estimate from the data that the weights should be scaled by " 
		  << std::setprecision(2) << (wtScale_ * wtScale) << "." << std::endl
		  << "You can use '" << name << ".wtscale = " << std::setprecision(2) << (wtScale_ * wtScale) << "' or '" 
		  << name << ".wtscale = auto' (or " << name << ".wtscale = auto, " << name << ".forcewt = true) to fix this." << std::endl
		  << "Your confidence intervals won't make sense until you do " << std::endl 
		  << "(reduced chi-squared will not be close to 1)", "red");
	if(namePrefix != "")
	  COUTCOLOR("(While loading data from file: " << file << ")", "red");
      }
    }
  }

  return false;
}

/**.......................................................................
 * Print warnings and return true if any data set should be reloaded with
 * different scaling
 */
bool VisDataSet::reloadMultiple()
{
  bool doReload = false;
  for(unsigned iDataSet=0; iDataSet < datasets_.size(); iDataSet++) {
    VisDataSet* vds = (VisDataSet*)datasets_[iDataSet];

    bool shift = false;
    Angle xShift, yShift;
    std::string file = parseFileName(fileList_[iDataSet], shift, xShift, yShift);

    if(vds->reloadSingle(name_, file))
      doReload = true;
  }

  return doReload;
}

/**.......................................................................
 * Method for counting data and initializing arrays
 */
void VisDataSet::initializeAndCountData(bool simulate)
{
  if(simulate)
    initializeAndCountDataSingle("");
  else
    initializeAndCountDataMultiple();
}

/**.......................................................................
 * Method for counting data and initializing arrays from a single file
 */
void VisDataSet::initializeAndCountDataSingle(std::string fileName)
{
  //------------------------------------------------------------
  // Parse file parameters
  //------------------------------------------------------------

  std::string file;

  if(!obsWasSet_)
    file = parseFileName(fileName, shiftRequested_, xShift_, yShift_);

  //------------------------------------------------------------
  // Parse parameters needed for reading in data
  //------------------------------------------------------------

  parseParameters();

  //------------------------------------------------------------
  // Now load the data from a file.  
  //------------------------------------------------------------

  if(!obsWasSet_) {
    initializeFromFile(file);
    countData(file);
    
    //------------------------------------------------------------
    // Else if obs was set, then we are
    // simulating, and no file has been specified.  In this case, just
    // initialize internal arrays to the size specified, etc.
    //------------------------------------------------------------
    
  } else {

    determineUniqueBaselineGroupings(obs_);
    
    if(usePerc_) {
      determineUvMax(obs_);
      initializeVisibilityArrays(percentCorrelation_);
    } else {
      initializeVisibilityArrays(image_);
    }

  }
}

/**.......................................................................
 * Method for counting data and initializing arrays from a single file
 */
void VisDataSet::initializeAndCountDataMultiple()
{
  //------------------------------------------------------------
  // Parse our own parameters
  //------------------------------------------------------------

  parseParameters();

  //------------------------------------------------------------
  // Parse each file and count data
  //------------------------------------------------------------

  std::map<std::string, std::string> exc;

  std::string file;
  bool shiftRequested;
  Angle xShift, yShift;

  // We exclude the store parameter from the instantiated datasets --
  // data read from multiple datasets will be stored in this object,
  // not the component datasets

  exc["store"] = "store";

  //------------------------------------------------------------
  // Iterate over each dataset, initializing it only to the point
  // where we can determine relevant information about the contents
  // (and not allocate internal memory needed to actually read the
  // data)
  //------------------------------------------------------------

  for(unsigned i=0; i < fileList_.size(); i++) {

    VisDataSet* dataset = datasets_[i];

    //------------------------------------------------------------
    // Copy any parameters that weren't otherwise specified (false
    // argument)
    //------------------------------------------------------------

    dataset->copyParameters(this, exc, false);
    dataset->initializeCommonParameters();
    dataset->initializeAndCountDataSingle(fileList_[i]);

    //------------------------------------------------------------
    // Now initialize internal arrays describing unique baseline
    // groupings in this file
    //------------------------------------------------------------

    dataset->determineUniqueBaselineGroupings(dataset->obs_);

    //------------------------------------------------------------
    // Update information about the source from the first file
    // encountered
    //------------------------------------------------------------

    if(i==0)
      updateObservationInformation(dataset);
  }

  //------------------------------------------------------------
  // Determine internal arrays from all baseline groups and
  // frequencies that were encountered
  //------------------------------------------------------------

  buildInternalMap(datasets_);
}

/**.......................................................................
 * Main method for loading data into this VisDataSet
 */
void VisDataSet::loadData(bool simulate)
{
  //------------------------------------------------------------
  // Count the data and initialize internal arrays
  //------------------------------------------------------------

  initializeAndCountData(simulate);

  //------------------------------------------------------------
  // Load data, checking for chisq that matches the weights
  //------------------------------------------------------------

  loadDataWithChecks();

  //------------------------------------------------------------
  // If debuggin was requested, print stats now
  //------------------------------------------------------------

  if(debug_) {
    printOccupiedIndices();
    printReImRms();
  }

  //------------------------------------------------------------
  // Compute primary beams
  //------------------------------------------------------------

  computePrimaryBeams();

  //------------------------------------------------------------
  // Compute the combined synthesized beam solid angle
  //------------------------------------------------------------

  computeGlobalSynthesizedBeam();

  //------------------------------------------------------------
  // And get the estimated synthesized beam parameters for all
  // VisFreqData subsets
  //------------------------------------------------------------

  estimateSynthesizedBeams();
}

/**.......................................................................
 * Load data from the specified file
 */
void VisDataSet::loadDataSingle(std::string fileName) 
{
  //------------------------------------------------------------
  // Initialize gridding arrays as specified
  //------------------------------------------------------------

  if(usePerc_) {
    initializeVisibilityArrays(percentCorrelation_);
  } else {
    initializeVisibilityArrays(image_);
  }

  //------------------------------------------------------------
  // Now, accumulate first moments
  //------------------------------------------------------------

  accumulateMoments(fileName, true, true, *this, VisDataSet::dispatchInternal);

  //------------------------------------------------------------
  // If we are estimating the error in the mean from the data
  // themselves, accumulate second moments now
  //------------------------------------------------------------

  if(estimateErrInMeanFromData_)
    accumulateMoments(fileName, false);

  //------------------------------------------------------------
  // Regardless of the method of estimating errors, convert from
  // variance to error in mean
  //------------------------------------------------------------

  calculateErrorInMean();

  //------------------------------------------------------------
  // Store weight sums and shifts
  //------------------------------------------------------------

  storeWtSums(xShift_, yShift_);

  //------------------------------------------------------------
  // If we were storing the data internally on read-in, release the
  // memory now (but only if we were told to; if simulation features
  // are being used, these may be wanted later).
  //------------------------------------------------------------

  if(releaseDataAfterReadin_)
    obs_.visibilities_.resize(0);
}

/**.......................................................................
 * Accumulate first or second moments, either by reading directly from
 * the file (slower), or by iterating through the internal copy of the
 * data that has been previously stored (faster, but can use a large
 * memory footprint).
 *
 * Note that the ability to do the latter requires a call to
 * storeDataInternallyOnReadin() prior to countData().
 */
void VisDataSet::accumulateMoments(std::string fileName, bool first)
{
  //------------------------------------------------------------
  // First, initialize all internal arrays needed for moment
  // calculation. 
  //------------------------------------------------------------

  initializeForMomentAccumulation(first);

  //------------------------------------------------------------
  // Open the file reader if we have not stored the data internally
  //------------------------------------------------------------

  if(!storeDataInternally_) {
    std::string file = parseFileName(fileName, shiftRequested_, xShift_, yShift_);
    openFileReader(file);
  }

  ObsInfo::Vis visFromFile;
  ObsInfo::Vis& vis = visFromFile;

  //------------------------------------------------------------
  // Iterate over all groups in the UVF file
  //------------------------------------------------------------

  Wavelength wavelength;

  for(unsigned iGroup=0; iGroup < obs_.nGroup_; iGroup++) {

    if(storeDataInternally_) {
      vis = obs_.visibilities_[iGroup];
    } else {
      getGroup(iGroup, vis);
    }

    //------------------------------------------------------------
    // Find the baseline grouping to which this group belongs
    //------------------------------------------------------------

    int baselineGroupIndex = aipsBaselineIndexToGroupIndexMap_[vis.baseline_];

    // Determine whether this visibility is valid

    bool goodBaseline = false;
    if(baselineGroupIndex >= 0) {
      goodBaseline = isfinite(vis.u_) && isfinite(vis.v_) && isfinite(vis.w_);
    } else {
      goodBaseline = false;
    }

    if(goodBaseline) {

      VisBaselineGroup& group = baselineGroups_[baselineGroupIndex];

      //------------------------------------------------------------
      // Now iterate over Stokes parameters and IFs to put the data in
      // the right place
      //------------------------------------------------------------

      Delay delay;
      VisData data;
      unsigned iVis=0;

      for(unsigned iStokes=0; iStokes < obs_.nStokes_; iStokes++) {
	VisStokesData& stokes = group.stokesData_[iStokes];

	for(unsigned iIf=0; iIf < obs_.nFreq_; iIf++) {

	  wavelength.setFrequency(obs_.frequencies_[iIf]);

	  VisFreqData& freq     = stokes.freqData_[iIf];

	  //------------------------------------------------------------
	  // Don't proceed if this IF number is excluded
	  //------------------------------------------------------------

	  delay.setDelayInSeconds(vis.u_);
	  data.u_  = (reverseDelays_ ? 1 : -1) * delay.cm() / wavelength.cm();

	  delay.setDelayInSeconds(vis.v_);
	  data.v_  = (reverseDelays_ ? 1 : -1) * delay.cm() / wavelength.cm();

	  delay.setDelayInSeconds(vis.w_);
	  data.w_  = (reverseDelays_ ? 1 : -1) * delay.cm() / wavelength.cm();
	  
	  data.r_ = sqrt(data.u_ * data.u_ + data.v_ * data.v_);
	  
	  data.re_ = vis.re_[iVis];
	  data.im_ = vis.im_[iVis];

	  //------------------------------------------------------------
	  // By convention, we assume that the weights represent the
	  // Re/Im variance, or half the variance on the total
	  // intensity
	  //
	  // We also apply any user-supplied scaling and taper at this
	  // point.
	  //
	  // If forceWt_ = true, we use the estimated weight scaling for
	  // each frequency band.  Else, we use the user-supplied or 
	  // estimated weight scaling for the whole data set
	  //------------------------------------------------------------

	  if(forceWt_)
	    data.wt_ = vis.wt_[iVis] * freq.wtScale_ * taper(data.u_, data.v_);
	  else
	    data.wt_ = vis.wt_[iVis] *      wtScale_ * taper(data.u_, data.v_);
	  
	  iVis++;

	  bool goodVis = isfinite(data.re_) && isfinite(data.im_) && isfinite(data.wt_) && data.wt_ > 0.0;

	  //------------------------------------------------------------
	  // Only grid data that fall within the specified UV min/max 
	  //------------------------------------------------------------

	  goodVis &= (uvMin_ < 0.0 || data.r_ > uvMin_);
	  goodVis &= (uvMax_ < 0.0 || data.r_ < uvMax_);

	  //------------------------------------------------------------
	  // If this data point should be gridded, grid it now
	  //------------------------------------------------------------

	  if(goodVis) {
	    
	    //------------------------------------------------------------
	    // Load this data point into the data gridder for the
	    // current VisFreqData object
	    //------------------------------------------------------------
	    
	    if(first) {
	      freq.griddedData_.accumulateFirstMoments( data.u_, data.v_, data.re_, data.im_, data.wt_);
	    } else {
	      freq.griddedData_.accumulateSecondMoments(data.u_, data.v_, data.re_, data.im_, data.wt_);
	    }

	    freq.iVis_++;

	    if(useIfNumber(freq.ifNo_))
	      freq.nVisUsed_++;

	    //------------------------------------------------------------
	    // Accumulate the variance for this VisFreqData object
	    //------------------------------------------------------------

	    double reCurr = data.re_ - freq.reMean_;
	    double imCurr = data.im_ - freq.imMean_;

	    reCurr *= reCurr;
	    imCurr *= imCurr;

	    double chisq;

	    chisq = reCurr * data.wt_;
	    freq.estChisq_ += (chisq - freq.estChisq_) / freq.iVis_;

	    chisq = imCurr * data.wt_;
	    freq.estChisq_ += (chisq - freq.estChisq_) / freq.iVis_;

	    freq.wtSumTotal_ += data.wt_;
	  }

	}
      }
    }

    if(iGroup % 1000 == 0) {
      std::cout << "\rReading..." << (100*iGroup)/obs_.nGroup_ << "%";
      fflush(stdout);
    }
  }

  std::cout << "\rReading...100%\r                   ";

  //------------------------------------------------------------
  // If data were read in from the file, close the file now
  //------------------------------------------------------------

  if(!storeDataInternally_)
    closeFileReader();
}

/**.......................................................................
 * Method to accumulate visibility moments, dispatching the data via
 * an externally supplied function
 */
void VisDataSet::accumulateMoments(std::string fileName, bool first, bool init, VisDataSet& vds, ACC_DISPATCH_FN(*dispatchFn))
{
  //------------------------------------------------------------
  // First, initialize all internal arrays needed for moment
  // calculation. 
  //------------------------------------------------------------

  initializeForMomentAccumulation(first);

  //------------------------------------------------------------
  // Open the file reader if we have not stored the data internally
  //------------------------------------------------------------

  if(!storeDataInternally_) {
    std::string file = parseFileName(fileName, shiftRequested_, xShift_, yShift_);
    openFileReader(file);
  }

  ObsInfo::Vis visFromFile;
  ObsInfo::Vis& vis = visFromFile;

  //------------------------------------------------------------
  // Iterate over all groups in the UVF file
  //------------------------------------------------------------

  Wavelength wavelength;

  for(unsigned iGroup=0; iGroup < obs_.nGroup_; iGroup++) {

    if(storeDataInternally_) {
      vis = obs_.visibilities_[iGroup];
    } else {
      getGroup(iGroup, vis);
    }

    //------------------------------------------------------------
    // Find the baseline grouping to which this group belongs
    //------------------------------------------------------------

    int baselineGroupIndex = aipsBaselineIndexToGroupIndexMap_[vis.baseline_];

    //------------------------------------------------------------
    // Determine whether this visibility is valid
    //------------------------------------------------------------

    bool goodBaseline = false;
    if(baselineGroupIndex >= 0) {
      goodBaseline = isfinite(vis.u_) && isfinite(vis.v_) && isfinite(vis.w_);
    } else {
      goodBaseline = false;
    }

    if(goodBaseline) {

      VisBaselineGroup& group = baselineGroups_[baselineGroupIndex];

      //------------------------------------------------------------
      // Now iterate over Stokes parameters and IFs to put the data in
      // the right place
      //------------------------------------------------------------

      Delay delay;
      VisData data;
      unsigned iVis=0;

      for(unsigned iStokes=0; iStokes < obs_.nStokes_; iStokes++) {
	VisStokesData& stokes = group.stokesData_[iStokes];

	for(unsigned iIf=0; iIf < obs_.nFreq_; iIf++) {

	  wavelength.setFrequency(obs_.frequencies_[iIf]);

	  VisFreqData& freq     = stokes.freqData_[iIf];

	  //------------------------------------------------------------
	  // Don't proceed if this IF number is excluded
	  //------------------------------------------------------------

	  delay.setDelayInSeconds(vis.u_);
	  data.u_  = (reverseDelays_ ? 1 : -1) * delay.cm() / wavelength.cm();

	  delay.setDelayInSeconds(vis.v_);
	  data.v_  = (reverseDelays_ ? 1 : -1) * delay.cm() / wavelength.cm();

	  delay.setDelayInSeconds(vis.w_);
	  data.w_  = (reverseDelays_ ? 1 : -1) * delay.cm() / wavelength.cm();
	  
	  data.r_ = sqrt(data.u_ * data.u_ + data.v_ * data.v_);
	  
	  data.re_ = vis.re_[iVis];
	  data.im_ = vis.im_[iVis];

	  //------------------------------------------------------------
	  // By convention, we assume that the weights represent the
	  // Re/Im variance, or half the variance on the total
	  // intensity
	  //
	  // We also apply any user-supplied scaling and taper at this
	  // point.
	  //
	  // If forceWt_ = true, we use the estimated weight scaling for
	  // each frequency band.  Else, we use the user-supplied or 
	  // estimated weight scaling for the whole data set
	  //------------------------------------------------------------

	  if(forceWt_)
	    data.wt_ = vis.wt_[iVis] * freq.wtScale_ * taper(data.u_, data.v_);
	  else
	    data.wt_ = vis.wt_[iVis] *      wtScale_ * taper(data.u_, data.v_);

	  //	  COUT("Reading re = " << data.re_ << " im = " << data.im_ << " wt = " << data.wt_ << " baseline = " << vis.baseline_ << " iGroup = " << iGroup << " jd = " << data.jd_);

	  data.jd_       = vis.jd_;
	  data.baseline_ = vis.baseline_;

	  iVis++;

	  bool goodVis = isfinite(data.re_) && isfinite(data.im_) && isfinite(data.wt_) && data.wt_ > 0.0;

	  //------------------------------------------------------------
	  // Only grid data that fall within the specified UV min/max 
	  //------------------------------------------------------------

	  goodVis &= (uvMin_ < 0.0 || data.r_ > uvMin_);
	  goodVis &= (uvMax_ < 0.0 || data.r_ < uvMax_);

	  //------------------------------------------------------------
	  // If this data point should be gridded, do it now
	  //------------------------------------------------------------

	  if(goodVis) {
	    dispatchFn(data, first, *this, vds, group, stokes, freq, vis);
	  }
	}
      }
    }

    if(iGroup % 1000 == 0) {
      std::cout << "\rReading..." << (100*iGroup)/obs_.nGroup_ << "%";
      fflush(stdout);
    }

  }

  std::cout << "\rReading...100%\r                   ";

  //------------------------------------------------------------
  // If data were read in from the file, close the file now
  //------------------------------------------------------------

  if(!storeDataInternally_)
    closeFileReader();
}

/**.......................................................................
 * Convert from first and second moments to mean and error in mean
 */
void VisDataSet::initializeForMomentAccumulation(bool first)
{
  //------------------------------------------------------------
  // Before reading in data, zero all relevant data and weight arrays
  //------------------------------------------------------------

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	// Tell the UV gridders whether or not the error in the mean
	// is being estimataed from the data, or by using the data
	// weights

	freqData.griddedData_.estimateErrorInMeanFromData(estimateErrInMeanFromData_);

	if(first)
	  freqData.griddedData_.initializeForFirstMoments();
	else
	  freqData.griddedData_.initializeForSecondMoments();

      }
    }
  }
}

/**.......................................................................
 * Convert from first and second moments to mean and error in mean
 */
void VisDataSet::calculateErrorInMean()
{
  //------------------------------------------------------------
  // Now that the data have been read in, convert second moments to
  // error in mean
  //------------------------------------------------------------

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	//------------------------------------------------------------
	// Convert internal weight sums to error in mean and mark
	// populated indices
	//------------------------------------------------------------

	freqData.griddedData_.calculateErrorInMean();

	//------------------------------------------------------------
	// Now that the gridded data index array has been populated,
	// copy it into the model component objects too
	//------------------------------------------------------------

	freqData.fourierModelComponent_.assignPopulatedIndicesFrom(freqData.griddedData_);
	freqData.compositeFourierModelDft_.assignPopulatedIndicesFrom(freqData.griddedData_);
      }
    }
  }
}

/**.......................................................................
 * Examine what we know so far about the antennas to determine how
 * many unique baseline groupings there are.
 * 
 * A baseline grouping is defined as all baselines of the same pair
 * of antenna types; i.e., for an array with 2 different antenna
 * types, call them type 1 and type 2, there would be 3 distinct
 * baseline groupings: 1-1, 2-2, and 1-2.
 */
void VisDataSet::determineUniqueBaselineGroupings(ObsInfo& obs)
{
  //------------------------------------------------------------
  // First check that antenna information has been set up before loading
  // visibility data.  This will affect how we group data on read-in
  //------------------------------------------------------------

  if(obs.antennas_.size() == 0) {
    ThrowError("You must specify antenna types before reading the visibility data, using the "
	       << "setAntennaType() or setAntennaDiameter() commands. Otherwise, I don't know "
	       << "how to group visibilities by antenna diameters");
  }

  //------------------------------------------------------------
  // Now iterate through antennas, checking for unique types
  //------------------------------------------------------------

  std::vector<Antenna> uniqueAnts;
  for(unsigned i=0; i < obs.antennas_.size(); i++) {
    checkAntenna(obs.antennas_[i], uniqueAnts);
  }
  
  //------------------------------------------------------------
  // Now we have unique tag IDs.  Maximum number of unique baseline
  // groupings will be (number of unique antenna types choose 2) +
  // number of ant types
  //------------------------------------------------------------

  unsigned nUAnt = uniqueAnts.size();

  std::ostringstream osAnts;
  getAntOs(osAnts, uniqueAnts);
  COUTCOLOR(osAnts.str(), "cyan");

  unsigned nIncludedAnt = (includedAntTypes_.size() > 0 ? includedAntTypes_.size() : nUAnt);

  //------------------------------------------------------------
  // Trap the user specifying more antennas than there actually are
  //------------------------------------------------------------

  if(nIncludedAnt > nUAnt)
    nIncludedAnt = nUAnt;

  unsigned nUBase = (nIncludedAnt * (nIncludedAnt-1))/2 + nIncludedAnt;

  baselineGroups_.resize(nUBase);

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& group = baselineGroups_[iGroup];
    group.initialize(obs.nStokes_, obs.frequencies_, obs.bandwidths_, debug_);
    group.dataset_ = this;
  }

  //------------------------------------------------------------
  // Iterate through unique antenna types, using the antenna idTags to
  // create a unique baseline grouping tag.  We then associate each
  // baseline tag with a baseline group.
  //------------------------------------------------------------

  for(unsigned iAnt1=0, iBase=0; iAnt1 < uniqueAnts.size(); iAnt1++) {
    Antenna& ant1 = uniqueAnts[iAnt1];

    for(unsigned iAnt2=iAnt1; iAnt2 < uniqueAnts.size(); iAnt2++) {
      Antenna& ant2 = uniqueAnts[iAnt2];
      unsigned baselineTag = ant1.idTag_ | ant2.idTag_;

      if(antennaTypeIsIncluded(ant1) && antennaTypeIsIncluded(ant2)) {
	baselineTagToGroupIndexMap_[baselineTag]  = iBase;

	//------------------------------------------------------------
	// And initialize the antennas for this baseline group
	//------------------------------------------------------------

	baselineGroups_[iBase].antennaPair_.first  = ant1;
	baselineGroups_[iBase].antennaPair_.second = ant2;

	iBase++;

      } else {
	baselineTagToGroupIndexMap_[baselineTag]  = -1;
      }
    }
  }

  std::ostringstream osBase;
  getBaseOs(osBase);
  COUTCOLOR(osBase.str(), "cyan");

  //------------------------------------------------------------
  // Now that we have the association between baseline tag and
  // baseline grouping index, iterate through all possible baselines,
  // constructing the map of AIPS baseline index to baseline group
  //------------------------------------------------------------

  for(unsigned iAnt1=0, iBase=0; iAnt1 < obs.antennas_.size()-1; iAnt1++) {
    Antenna& ant1 = obs.antennas_[iAnt1];

    for(unsigned iAnt2=iAnt1+1; iAnt2 < obs.antennas_.size(); iAnt2++) {
      Antenna& ant2 = obs.antennas_[iAnt2];

      //------------------------------------------------------------
      // Calculate the AIPS-standard baseline index
      //------------------------------------------------------------

      unsigned aipsBaselineIndex  = (iAnt1+1) * 256 + (iAnt2+1);

      // Construct a bitmask of the two antenna types

      unsigned baselineTag        = ant1.idTag_ | ant2.idTag_;

      //------------------------------------------------------------
      // Only proceed if this antenna type and antenna number is
      // included by the current user selection
      //------------------------------------------------------------

      if(antennaIsIncluded(ant1) && antennaIsIncluded(ant2)) {

	// Get the baseline group index that corresponds to these
	// two antenna types
	
	if(baselineTagToGroupIndexMap_[baselineTag] < 0) {
	  aipsBaselineIndexToGroupIndexMap_[aipsBaselineIndex] = -1;
	  aipsBaselineIndexCount_[aipsBaselineIndex] = 0;
	  continue;
	}

	unsigned baselineGroupIndex = baselineTagToGroupIndexMap_[baselineTag];
	
	aipsBaselineIndexToGroupIndexMap_[aipsBaselineIndex] = baselineGroupIndex;
	aipsBaselineIndexCount_[aipsBaselineIndex] = 0;
	
	// For convenience, add this baseline to the list of baselines
	// associated with this group.  These will be used for simulation
	// purposes only
	
	VisBaselineGroup& group = baselineGroups_[baselineGroupIndex];
	group.addBaseline(aipsBaselineIndex, ant1, ant2);

	++iBase;

	//------------------------------------------------------------
	// If either antenna is excluded, set this baseline to point
	// to no valid baselinegroup
	//------------------------------------------------------------

      } else {
	aipsBaselineIndexToGroupIndexMap_[aipsBaselineIndex] = -1;
	aipsBaselineIndexCount_[aipsBaselineIndex] = 0;
      }
    }
  }

  //------------------------------------------------------------
  // Finally, initialize multithreading data in the VisFreqData
  // objects, regardless of whether or not it is used
  //------------------------------------------------------------

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& group = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < group.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = group.stokesData_[iStokes];

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	if(freqData.execData_) {
	  delete freqData.execData_;
	  freqData.execData_ = 0;
	}

	freqData.execData_ = new VisExecData(this, &freqData, iGroup, iStokes, iFreq);

	//------------------------------------------------------------
	// Also calculate the estimated primary beam for this group
	// and frequency, and store the max over all baseline groups
	// and frequencies
	//------------------------------------------------------------

	Angle pbHw = group.estimatePrimaryBeamHalfWidth(freqData.frequency_);

	if(pbHw > maxPrimaryBeamHalfwidth_)
	  maxPrimaryBeamHalfwidth_ = pbHw;
      }
    }
  }
}

/**.......................................................................
 * Check if the current antenna type has been encountered before, or
 * if this is a new antenna type.
 */
void VisDataSet::checkAntenna(Antenna& ant, std::vector<Antenna>& uniqueAnts)
{
  if(uniqueAnts.size() == 0) {
    ant.idTag_ = 1;
    uniqueAnts.push_back(ant);
  } else {
    
    unsigned i=0;
    for(i=0; i < uniqueAnts.size(); i++) {
      Antenna& currAnt = uniqueAnts[i];
      
      //------------------------------------------------------------
      // If this antenna was not specified by diameter, and it
      // matches the type of the current antenna, it is the same
      // type of antenna
      //------------------------------------------------------------

      if(ant.type_ != Antenna::ANT_DIAM && ant.type_ == currAnt.type_) {
	ant.idTag_ = currAnt.idTag_;
	break;
      }
      
      //------------------------------------------------------------
      // If this antenna was specified by diameter, and it has the
      // same diameter as the current antenna, treat it as the
      // same type of antenna
      //------------------------------------------------------------
      
      if(ant.type_ == Antenna::ANT_DIAM && currAnt.type_ && ant.diameter_ == currAnt.diameter_) {
	ant.idTag_ = currAnt.idTag_;
	break;
      }
    }
    
    //------------------------------------------------------------
    // If no match, this is a new unique antenna type
    //------------------------------------------------------------
    
    if(i == uniqueAnts.size()) {
      ant.idTag_ = (uniqueAnts[uniqueAnts.size()-1].idTag_ << 1);
      uniqueAnts.push_back(ant);
    }
  }
}

/**.......................................................................
 * Guess at the type of antennas, based on the name of the telescope
 * from the UVF header.
 */
void VisDataSet::guessAtAntennaType()
{
  String telStr(obs_.getTelescopeName());

  if(obs_.antennas_.size() == 0) {
    COUT("The number and type of antennas is unknown (no AIPS AN table is present in the UVF file)");
    return;
  }

  bool unknownAnts = false;
  for(unsigned i=0; i < obs_.antennas_.size(); i++) {
    if(obs_.antennas_[i].type_ == Antenna::ANT_UNKNOWN) {
      unknownAnts = true;
      break;
    }
  }

  if(unknownAnts) {

    //------------------------------------------------------------
    // Check if this is the SZA
    //------------------------------------------------------------

    if((telStr.contains("CARMA") || telStr.contains("SZA")) && obs_.antennas_.size() == 8) {
      obs_.setAntennaTypeIfUnknown(Antenna::ANT_SZA);
      COUTCOLOR("There are " << obs_.antennas_.size() << " antennas, some of which have not been specified.  " << std::endl << "Judging from the telescope name, "
		"I'm going to assume they are SZA (3.5m) antennas", "cyan");

      //------------------------------------------------------------
      // Check if OVRO
      //------------------------------------------------------------

    } else if((telStr.contains("CARMA") || telStr.contains("OVRO")) && obs_.antennas_.size() == 6) {
      obs_.setAntennaTypeIfUnknown(Antenna::ANT_OVRO);
      COUTCOLOR("There are " << obs_.antennas_.size() << " antennas, some of which have not been specified.  " << std::endl << "Judging from the telescope name, "
		"I'm going to assume they are OVRO (10.4m) antennas", "cyan");

      //------------------------------------------------------------
      // Check if BIMA
      //------------------------------------------------------------

    } else if((telStr.contains("CARMA") || telStr.contains("BIMA")) && obs_.antennas_.size() == 9) {
      obs_.setAntennaTypeIfUnknown(Antenna::ANT_BIMA);
      COUTCOLOR("There are " << obs_.antennas_.size() << " antennas, some of which have not been specified.  " << std::endl << "Judging from the telescope name, "
		"I'm going to assume they are BIMA (6.1m) antennas", "cyan");

      //------------------------------------------------------------
      // Check if CARMA15 (also called OVRO)
      //------------------------------------------------------------

    } else if((telStr.contains("CARMA") || telStr.contains("OVRO") || telStr.contains("BIMA") || telStr.contains("SZA")) && obs_.antennas_.size() == 15) {

      COUTCOLOR("There are " << obs_.antennas_.size() << " antennas, some of which have not been specified.  " << std::endl << "Judging from the telescope name, "
		"I'm going to assume this is CARMA data." << std::endl, "cyan");

      for(unsigned iAnt=0; iAnt < 6; iAnt++)
	obs_.setAntennaTypeIfUnknown(Antenna::ANT_OVRO, iAnt);

      for(unsigned iAnt=6; iAnt < 15; iAnt++)
	obs_.setAntennaTypeIfUnknown(Antenna::ANT_BIMA, iAnt);

      //------------------------------------------------------------
      // Check if CARMA23 (also called OVRO)
      //------------------------------------------------------------

    } else if((telStr.contains("CARMA") || telStr.contains("OVRO") || telStr.contains("BIMA") 
	       || telStr.contains("SZA")) && obs_.antennas_.size() == 23) {

      COUTCOLOR("There are " << obs_.antennas_.size() << " antennas, some of which have not been specified.  " << std::endl << "Judging from the telescope name, "
		"I'm going to assume this is CARMA data." << std::endl, "cyan");

      std::ostringstream os;
      for(unsigned iAnt=0; iAnt < 6; iAnt++) {
	os.str("");
	os << "ant" << iAnt+1 << ".type";

	if(!obs_.getParameter(os.str(), false)->data_.hasValue()) {
	  obs_.setAntennaTypeIfUnknown(Antenna::ANT_OVRO, iAnt);
	}
      }

      for(unsigned iAnt=6; iAnt < 15; iAnt++) {
	os.str("");
	os << "ant" << iAnt+1 << ".type";
      
	if(!obs_.getParameter(os.str(), false)->data_.hasValue()) {
	  obs_.setAntennaTypeIfUnknown(Antenna::ANT_BIMA, iAnt);
	}
      }

      for(unsigned iAnt=15; iAnt < 23; iAnt++) {
	os.str("");
	os << "ant" << iAnt+1 << ".type";
      
	if(!obs_.getParameter(os.str(), false)->data_.hasValue()) {
	  obs_.setAntennaTypeIfUnknown(Antenna::ANT_SZA, iAnt);
	}
      }

      //------------------------------------------------------------
      // Check if VLA
      //------------------------------------------------------------

    } else if(telStr.contains("VLA")) {
      obs_.setAntennaTypeIfUnknown(Antenna::ANT_VLA);
      COUTCOLOR("There are " << obs_.antennas_.size() << " antennas, some of which have not been specified.  " << std::endl << "Judging from the telescope name, "
		"I'm going to assume they are VLA (27m) antennas", "cyan");
    }
  }

  //-------------------------------------------------------------
  // Plot antennas, if requested
  //-------------------------------------------------------------

  bool antplot = false;
  
  if(getParameter("antplot", false)->data_.hasValue())
    antplot = getBoolVal("antplot");
    
  if(antplot)
    obs_.plotAntennas();
  
  if(getParameter("reversedelays", false)->data_.hasValue())
    reverseDelays_ = getBoolVal("reversedelays");

  if(getParameter("reversedisplay", false)->data_.hasValue()) {
    if(getBoolVal("reversedisplay")) {
      reverseDisplay_ = !reverseDisplay_;
      COUT("Reversing display: " << reverseDisplay_);
    }
  }

  return;
}

/**.......................................................................
 * Print information about the file just read in
 */
void VisDataSet::printFileStats(std::string fileName)
{
  std::ostringstream os;
  os << std::endl << "UVF File " << fileName << " contains observations of: " << std::endl << std::endl
    
     << "  Object     = " << obs_.getSourceName() << std::endl
     << "  RA         = " << " " << obs_.getObsRa()      << std::endl
     << "  DEC        = " << obs_.getObsDec()     << std::endl 
     << "  Equinox    = " << obs_.getObsEquinox() << std::endl << std::endl

     << "with: " << std::endl << std::endl

     << "  Telescope  = " << obs_.getTelescopeName()  << std::endl
     << "  Instrument = " << obs_.getInstrumentName() << std::endl << std::endl

     << "There are "   << String::formatHumanReadableInteger(obs_.nGroup_ * obs_.nStokes_ * obs_.nFreq_) << " total visibilities, " << "comprising: " << std::endl << std::endl;

  if(obs_.antennas_.size() > 0) {
    os << "  " << std::setw(12) << std::right << String::formatHumanReadableInteger(obs_.nGroup_/obs_.nBaseline_) << " timestamps" << std::endl
       << " x" << std::setw(12) << std::right << obs_.nBaseline_         << " baselines"  << std::endl;
  } else {
    os << "  " << std::setw(12) << std::right << obs_.nGroup_ << " timestamps x baselines" << std::endl;
  }

  os << " x" << std::setw(12) << std::right << obs_.nStokes_            
     << (obs_.nStokes_ > 1 ? " Stokes parameters" : " Stokes parameter") << std::endl
     << " x" << std::setw(12) << std::right << obs_.nFreq_ 
     << (obs_.nFreq_ > 1 ? " IFs: " : " IF: ") << std::endl;

  for(unsigned i=0; i < obs_.nFreq_; i++) {
    os << std::setw(25) << std::right << std::setprecision(3) << (i+1) 
       << ") " << std::setprecision(1) << std::fixed << std::setfill(' ') << obs_.frequencies_[i].GHz() << " GHz" << std::setw(5) << " " 
       << std::setw(6) << std::setprecision(1) << std::fixed << std::right << std::setfill(' ') << obs_.bandwidths_[i].MHz() << " MHz";
    if(i < obs_.nFreq_ - 1) {
      os << std::endl;
    }
  }

  std::ostringstream osIf;
  getIfOs(osIf);
  os << osIf.str();
  
  os << std::endl;

  COUTCOLOR(os.str(), "cyan");

#if 0
  for(unsigned i=0; i < obs_.antennas_.size(); i++) {
    COUT("Ant " << i << " = " << obs_.antennas_[i].antNo_);
  }
#endif
}

/**.......................................................................
 * Remove the composite model from the data
 */
void VisDataSet::remModel()
{
  transformModels();

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];
    
    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];
	freqData.remModel();
      }
    }
  }

  clearModel();
}

/**.......................................................................
 * Add a model component to the composite image-plane model
 */
void VisDataSet::addModel(Model& model)
{
  //------------------------------------------------------------
  // Only proceed if this model applies to this data set
  //------------------------------------------------------------

  if(applies(model)) {

    //------------------------------------------------------------
    // Dynamic cast will throw if this model cannot be cast to the
    // right type, but the call to applies() should protect against
    // this
    //------------------------------------------------------------

    Generic2DAngularModel& model2D = dynamic_cast<Generic2DAngularModel&>(model);

    initWait();

    for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
      VisBaselineGroup& groupData = baselineGroups_[iGroup];
      
      for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
	VisStokesData& stokesData = groupData.stokesData_[iStokes];
	
	for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	  VisFreqData& freqData = stokesData.freqData_[iFreq];
	  addModelMultiThread(freqData, model2D, iGroup, iStokes, iFreq);
	}
      }
    }

    waitUntilDone();
  }
}

/**.......................................................................
 * Multi-thread-aware version of addModel
 */
void VisDataSet::addModelMultiThread(VisFreqData& vfd, Generic2DAngularModel& model, 
				     unsigned iGroup, unsigned iStokes, unsigned iFreq)
{
  if(!pool_) {
    vfd.addModel(model);
  } else {
    VisExecData* ved = vfd.execData_;
    ved->initialize(&model);
    registerPending(iGroup, iStokes, iFreq);
    pool_->execute(&execAddModel, ved);
  }
}

/**.......................................................................
 * Static method which can be passed to a thread pool, to add a model to a 
 * single VisFreqData
 */
EXECUTE_FN(VisDataSet::execAddModel)
{
  VisExecData*           ved   = (VisExecData*)args;
  VisDataSet*            vds   = ved->vds_;
  VisFreqData*           vfd   = ved->vfd_;
  Generic2DAngularModel* model = ved->model_;

  vfd->addModel(*model);
  vds->registerDone(ved->iGroup_, ved->iStokes_, ved->iFreq_);
}

/**.......................................................................
 * Clear all model components
 */
void VisDataSet::clearModel()
{
  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	freqData.clearModel();
      }
    }
  }
}

/**.......................................................................
 * Compute the chi-squared of this dataset with the given model
 */
ChisqVariate VisDataSet::computeChisq()
{
  //------------------------------------------------------------
  // Transform models
  //------------------------------------------------------------

  transformModels();

  //------------------------------------------------------------
  // Accumulate chi-squared over all frequencies
  //------------------------------------------------------------

  ChisqVariate chisq = accumulateChisq();

  //------------------------------------------------------------
  // Now clear models so that the next call to addModel()
  // re-initializes the model (rather than adds to it).
  //------------------------------------------------------------

  clearModel();

  //------------------------------------------------------------
  // Return the chi-squared
  //------------------------------------------------------------

  return chisq;
}

/**.......................................................................
 * Transform models
 */
void VisDataSet::transformModels()
{
  initWait();

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	transformModelMultiThread(freqData, iGroup, iStokes, iFreq);
      }
    }
  }

  waitUntilDone();
}

/**.......................................................................
 * Inverse transform data
 */
void VisDataSet::inverseTransformData()
{
  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];
	freqData.inverseTransformData();
      }
    }
  }
}

/**.......................................................................
 * Static method which can be passed to a thread pool, to transform a
 * single VisFreqData image
 */
EXECUTE_FN(VisDataSet::execTransformModel)
{
  VisExecData* ved = (VisExecData*)args;
  VisDataSet*  vds = ved->vds_;
  VisFreqData* vfd = ved->vfd_;

  vfd->transformModel();
  vds->registerDone(ved->iGroup_, ved->iStokes_, ved->iFreq_);
}

/**.......................................................................
 * Multi-thread-aware version of transformModel
 */
void VisDataSet::transformModelMultiThread(VisFreqData& vfd, unsigned iGroup, unsigned iStokes, unsigned iFreq)
{
  if(!pool_) {
    vfd.transformModel();
  } else {
    VisExecData* ved = vfd.execData_;
    registerPending(iGroup, iStokes, iFreq);
    pool_->execute(&execTransformModel, ved);
  }
}

/**.......................................................................
 * Compute the chi-squared of this dataset with the given model
 */
ChisqVariate VisDataSet::accumulateChisq()
{
  ChisqVariate chisq;

  initWait();

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	if(freqData.hasData())
	  computeChisqMultiThread(freqData, chisq, iGroup, iStokes, iFreq);
      }
    }
  }
  
  waitUntilDone();

  return chisq;
}

/**.......................................................................
 * Accumulate data over all baselines, and frequencies, into internal
 * gridders.
 */
void VisDataSet::accumulate(VisDataSet::AccumulatorType type)
{
  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	if(freqData.hasData()) {
	  freqData.utilityGridder_.initializeForFirstMoments();
	  freqData.accumulate(freqData.utilityGridder_, type);
	  freqData.utilityGridder_.shift();
	  freqData.utilityGridder_.computeInverseTransform();
	}
      }
    }
  }
}

/**.......................................................................
 * Accumulate data over all baselines, and frequencies, into the
 * passed data gridder.
 */
void VisDataSet::accumulate(VisDataSet::AccumulatorType type, UvDataGridder& gridder)
{
  gridder.initializeForFirstMoments();

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	if(freqData.hasData()) {
	  freqData.accumulate(gridder, type);
	}
      }
    }
  }

  gridder.shift();
  gridder.computeInverseTransform();
}

/**.......................................................................
 * Accumulate a dirty map, corrected for the primary beam
 */
void VisDataSet::accumulateBeamCorrectedImage(Image& image, Image& wtimage)
{
  double power = 0.1;

  if(getParameter("power", false)->data_.hasValue())
    power = getDoubleVal("power");

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    Length diam1 = groupData.antennaPair_.first.getDiameter();
    Length diam2 = groupData.antennaPair_.second.getDiameter();

    // For this group, store the smaller of the diameters

    Length diameter = diam1 < diam2 ? diam1 : diam2;

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	if(freqData.hasData()) {

	  Wavelength wave(freqData.frequency_);

	  // Estimate the FWHM of this primary beam

	  Angle sigma(Angle::Radians(), (1.22 * wave.meters() / diameter.meters()) / sqrt(8*log(2.0)));

	  Image& pb = freqData.primaryBeam_;
	  double wt = freqData.griddedData_.wtSumTotal_;
	  Image currImage = freqData.utilityGridder_.getImage(false);

	  // Zero this image beyond the specified radius

	  Image incrNum = (pb * currImage) * wt;
	  Image incrDen = (pb * pb) * wt;

	  incrNum.setRaDecFft(ra_, dec_);
	  incrDen.setRaDecFft(ra_, dec_);

	  incrNum.invalidateBeyondRadius(sigma * sqrt(-2*log(power)));;
	  incrDen.invalidateBeyondRadius(sigma * sqrt(-2*log(power)));;

	  image   += incrNum;
	  wtimage += incrDen;
	}
      }
    }
  }
}

/**.......................................................................
 * Static method which can be passed to a thread pool, to accumulate
 * chisq for a single VisFreqData image
 */
EXECUTE_FN(VisDataSet::execComputeChisq)
{
  VisExecData*  ved   = (VisExecData*)args;
  VisDataSet*   vds   = ved->vds_;
  VisFreqData*  vfd   = ved->vfd_;
  ChisqVariate* chisq = ved->chisq_;
  ChisqVariate  chisqTmp;

  chisqTmp = vfd->computeChisq();

  ved->dataAccessGuard_.lock();
  (*chisq) += chisqTmp;
  ved->dataAccessGuard_.unlock();

  vds->registerDone(ved->iGroup_, ved->iStokes_, ved->iFreq_);
}

/**.......................................................................
 * Multi-thread-aware version of computeChisq
 */
void VisDataSet::computeChisqMultiThread(VisFreqData& vfd, ChisqVariate& chisq, unsigned iGroup, unsigned iStokes, unsigned iFreq)
{
  if(!pool_) {
    chisq += vfd.computeChisq();
  } else {
    VisExecData* ved = vfd.execData_;
    ved->initialize(&chisq);
    registerPending(iGroup, iStokes, iFreq);
    pool_->execute(&execComputeChisq, ved);
  }
}

/**.......................................................................
 * Check the validity of a baseline
 */
void VisDataSet::checkBaselineGroupIndex(unsigned iGroup)
{
  if(iGroup > baselineGroups_.size()-1) {
    ThrowError("Invalid group index: " << iGroup << ".  Should be < " << baselineGroups_.size());
  }
}

void VisDataSet::plotUv(int iGroup, int iStokes, int iFreq)
{
  std::vector<float> x;
  std::vector<float> y;

  getVis(Dft2d::DATA_UV, x, y, iGroup, iStokes, iFreq);

  PgUtil::linePlot(x.size(), &x[0], &y[0], 0, "U", "V", "UV plot", false, false);
}

void VisDataSet::plotReal(int iGroup, int iStokes, int iFreq)
{
  std::vector<float> x;
  std::vector<float> y;

  getVis(Dft2d::DATA_REAL, x, y, iGroup, iStokes, iFreq);

  PgUtil::linePlot(x.size(), &x[0], &y[0], 0, "UV Radius", "Re", "Real plot", false, false);
}

void VisDataSet::plotImag(int iGroup, int iStokes, int iFreq)
{
  std::vector<float> x;
  std::vector<float> y;

  getVis(Dft2d::DATA_IMAG, x, y, iGroup, iStokes, iFreq);

  PgUtil::linePlot(x.size(), &x[0], &y[0], 0, "UV Radius", "Im", "Imag plot", false, false);
}

void VisDataSet::plotAbs(int iGroup, int iStokes, int iFreq)
{
  std::vector<float> x;
  std::vector<float> y;

  getVis(Dft2d::DATA_ABS, x, y, iGroup, iStokes, iFreq);

  PgUtil::linePlot(x.size(), &x[0], &y[0], 0, "UV Radius", "Abs", "Abs plot", false, false);
}

/**.......................................................................
 * Extract the gridded visibility data from this object
 */
void VisDataSet::getVis(Dft2d::DataType type,
			std::vector<float>& x, std::vector<float>& y, 
			int iGroup, int iStokes, int iFreq)
{
  unsigned groupStart,  groupStop;
  unsigned stokesStart, stokesStop;
  unsigned freqStart,   freqStop;
  unsigned nVis = 0;

  if(iGroup < 0) {
    groupStart = 0;
    groupStop  = baselineGroups_.size()-1;
  } else {
    checkBaselineGroupIndex(iGroup);
    groupStart = groupStop = iGroup;
  }

  for(unsigned iGroup=groupStart; iGroup <= groupStop; iGroup++) {
    VisBaselineGroup& group = baselineGroups_[iGroup];
    
    if(iStokes < 0) {
      stokesStart = 0;
      stokesStop  = group.stokesData_.size()-1;
    } else {
      group.checkStokesIndex(iStokes);
      stokesStart = stokesStop = iStokes;
    }

    for(unsigned iStokes=stokesStart; iStokes <= stokesStop; iStokes++) {
      VisStokesData& stokesData = group.stokesData_[iStokes];
      
      if(iFreq < 0) {
	freqStart = 0;
	freqStop  = stokesData.freqData_.size()-1;
      } else {
	stokesData.checkFrequencyIndex(iFreq);
	freqStart = freqStop = iFreq;
      }

      for(unsigned iFreq=freqStart; iFreq <= freqStop; iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];
	nVis += freqData.griddedData_.populatedIndices_.size();
      }
    }
  }

  x.resize(nVis);
  y.resize(nVis);

  unsigned iDat = 0;
  for(unsigned iGroup=groupStart; iGroup <= groupStop; iGroup++) {
    VisBaselineGroup& group = baselineGroups_[iGroup];

    if(iStokes < 0) {
      stokesStart = 0;
      stokesStop  = group.stokesData_.size()-1;
    } else {
      stokesStart = stokesStop = iStokes;
    }

    for(unsigned iStokes=stokesStart; iStokes <= stokesStop; iStokes++) {
      VisStokesData& stokesData = group.stokesData_[iStokes];
    
      if(iFreq < 0) {
	freqStart = 0;
	freqStop  = stokesData.freqData_.size()-1;
      } else {
	freqStart = freqStop = iFreq;
      }

      for(unsigned iFreq=freqStart; iFreq <= freqStop; iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	for(unsigned iVis=0; iVis < freqData.griddedData_.populatedIndices_.size(); iVis++) {

	  // Extract the (u, v) coordinate, and the requested value
	  // from the data grid container

	  double u, v, val, r;
	  freqData.griddedData_.getUVData(freqData.griddedData_.populatedIndices_[iVis], type, u, v, val);

	  // For convenience, precompute the uv radius

	  r = sqrt(u*u + v*v);

	  switch (type) {
	  case Dft2d::DATA_UV:
	    x[iDat] = u;
	    y[iDat] = v;
	    break;
	  case Dft2d::DATA_REAL:
	  case Dft2d::DATA_IMAG:
	  case Dft2d::DATA_ABS:
	    x[iDat] = r;
	    y[iDat] = val;
	    break;
	  default:
	    break;
	  }

	  ++iDat;
	}
      }
    }
  }
}

/**.......................................................................
 * Once data have been sorted into unique baseline pairings (groups),
 * purge any groups that don't contain any data.
 */
void VisDataSet::purgeZeroSizedGroups()
{
  bool found=false,anyFound=false;
  unsigned eraseIndex;

  do {
    found = false;

    for(unsigned i=0; i < baselineGroups_.size(); i++) {

      if(baselineGroups_[i].nBaseline_ == 0) {
	found = true;
	anyFound = true;
	baselineGroups_.erase(baselineGroups_.begin() + i);
	break;
      }

    }

  } while(found);

  if(anyFound) {
    COUTCOLOR("There are now only " << baselineGroups_.size() 
	      << " distinct baseline groupings with valid data and/or that you have requested", "cyan");
  }

}

/**.......................................................................
 * Initialize estimated chisq 
 */
void VisDataSet::initializeEstChisq()
{
  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& group = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < group.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = group.stokesData_[iStokes];

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];
	freqData.iVis_     = 0;
	freqData.estChisq_ = 0.0;
      }
    }
  }
}

/**.......................................................................
 * Initialize estimated chisq 
 */
void VisDataSet::storeEstChisqAsWtScale()
{
  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& group = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < group.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = group.stokesData_[iStokes];

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];
	freqData.wtScale_ = 1.0/freqData.estChisq_;
      }
    }
  }
}

/**.......................................................................
 * Resize internal Dfts to the size needed to grid the data to the
 * requested percent correlation
 */
void VisDataSet::initializeVisibilityArrays(double percentCorrelation)
{
  UvDataGridder* max = 0;
  
  Length diameter1, diameter2;
  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& group = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < group.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = group.stokesData_[iStokes];

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	//------------------------------------------------------------
	// Now check what percent correlation we want for combining data
	//------------------------------------------------------------

	if(percentCorrelation < 1.0) {

	  //------------------------------------------------------------
	  // Get the correlation length corresponding to this baseline
	  // group, this frequency
	  //------------------------------------------------------------
	  
	  diameter1 = group.antennaPair_.first.getDiameter();
	  diameter2 = group.antennaPair_.second.getDiameter();
	  double correlationLength = 
	    Dft2d::correlationLength(diameter1, diameter2,
				     freqData.frequency_, percentCorrelation);
	  
	  //------------------------------------------------------------
	  // Resize images and dfts to match the correlation length
	  //------------------------------------------------------------

	  freqData.resize(percentCorrelation, correlationLength, obsWasSet_);

	  UvDataGridder* curr = &freqData.griddedData_;

	  if(max == 0 
	     || curr->xAxis().getNpix() > max->xAxis().getNpix()
	     || curr->yAxis().getNpix() > max->yAxis().getNpix()) {
	    max = curr;
	  }

	} else {
	  ThrowColorError(std::endl << "I won't let you specify a correlation percentage of 1.  " << std::endl
			  << "This corresponds to a correlation length of zero, which would require an infinitely large image. " << std::endl
			  << "Use " << name_ << ".size and " << name_ << ".npix if you want to manually control the gridding", "red");
	}

	freqData.iVis_     = 0;
	freqData.estChisq_ = 0.0;
      }
    }
  }

  initializeGlobalGridders(max);
}

/**.......................................................................
 * Resize internal Dfts to the size needed to grid the data to the
 * requested percent correlation
 */
void VisDataSet::initializeVisibilityArrays2(double percentCorrelation)
{
  bool first = true;
  Angle maxRes;
  Angle maxSize;
  unsigned maxNpix;

  Length diameter1, diameter2;
  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& group = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < group.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = group.stokesData_[iStokes];

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	//------------------------------------------------------------
	// Now check what percent correlation we want for combining data
	//------------------------------------------------------------

	if(percentCorrelation < 1.0) {

	  // Get the correlation length corresponding to this baseline
	  // group, this frequency
	  
	  diameter1 = group.antennaPair_.first.getDiameter();
	  diameter2 = group.antennaPair_.second.getDiameter();
	  double correlationLength = 
	    Dft2d::correlationLength(diameter1, diameter2,
				     obs_.frequencies_[iFreq], percentCorrelation);
	  
	  // Resize images and dfts to match the correlation length

	  freqData.resize(percentCorrelation, correlationLength, obsWasSet_);

	  UvDataGridder* curr = &freqData.griddedData_;
	  Angle xRes  = curr->xAxis().getAngularResolution();
	  Angle yRes  = curr->yAxis().getAngularResolution();
	  Angle xSize = curr->xAxis().getAngularSize();
	  Angle ySize = curr->yAxis().getAngularSize();
	  unsigned xNpix = curr->xAxis().getNpix();
	  unsigned yNpix = curr->yAxis().getNpix();

	  if(first) {
	    maxRes = xRes;
	    maxSize = xSize;
	    maxNpix = xNpix;
	    first = false;
	  }

	  maxRes  = (maxRes < xRes) ? maxRes : xRes;
	  maxRes  = (maxRes < yRes) ? maxRes : yRes;
	  maxSize = (maxSize > xSize) ? maxSize : xSize;
	  maxSize = (maxSize > ySize) ? maxSize : ySize;
	  maxNpix = (maxNpix > xNpix) ? maxNpix : xNpix;
	  maxNpix = (maxNpix > yNpix) ? maxNpix : yNpix;
	  
	} else {
	  ThrowColorError(std::endl << "I won't let you specify a correlation percentage of 1.  " << std::endl
			  << "This corresponds to a correlation length of zero, which would require an infinitely large image. " << std::endl
			  << "Use " << name_ << ".size and " << name_ << ".npix if you want to manually control the gridding", "red");
	}
      }
    }
  }

  UvDataGridder max;
  unsigned nPix = Dft2d::nearestPowerOf2NotLessThan(maxSize / maxRes);

  maxSize.setDegrees(0.2);
  maxNpix = 128;
  max.initializeForVis(maxSize, maxSize, maxNpix, maxNpix);

  initializeGlobalGridders(&max);
}

/**.......................................................................
 * Return the image required to grid the data to the requested percent
 * correlation
 */
Image VisDataSet::getRequiredImage(double percentCorrelation)
{
  UvDataGridder* max = 0;

  unsigned nPixMax      = 0;
  double   maxSizeInRad = 0.0;
  bool first = true;

  Length diameter1, diameter2;
  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& group = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < group.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = group.stokesData_[iStokes];

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	//------------------------------------------------------------
	// Now check what percent correlation we want for combining data
	//------------------------------------------------------------

	if(percentCorrelation < 1.0) {

	  diameter1 = group.antennaPair_.first.getDiameter();
	  diameter2 = group.antennaPair_.second.getDiameter();

	  double correlationLength = 
	    Dft2d::correlationLength(diameter1, diameter2,
				     obs_.frequencies_[iFreq], percentCorrelation);

	  // Find the power-of-2 size of the array that will grid the
	  // data _at least_ this finely
	  
	  double absMax = freqData.uAbsMax_ > freqData.vAbsMax_ ? freqData.uAbsMax_ : freqData.vAbsMax_;
	  
	  unsigned nPix = 
	    Dft2d::nearestPowerOf2NotLessThan(sqrt(2.0) * absMax / correlationLength);
	  
	  double spatialFrequencyResolution = (absMax / (nPix/2));
	  double sizeInRad = 1.0/spatialFrequencyResolution;

	  nPixMax      = nPixMax > nPix ? nPixMax : nPix;

	  if(first) {
	    maxSizeInRad = sizeInRad;
	    first = false;
	  } else {
	    maxSizeInRad = maxSizeInRad < sizeInRad ? maxSizeInRad : sizeInRad;
	  }

	} else {
	  ThrowError("I won't let you specify a correlation percentage of 1.  "
		     << " This corresponds to a correlation length of zero, "
		     << "which would require an infinitely large image.");
	}
      }
    }
  }

  Angle size(Angle::Radians(), maxSizeInRad);
  Image retVal(nPixMax, size);

  return retVal;
}

/**.......................................................................
 * Initialize data gridders use to store data for all frequencies (for
 * displaying image-plane maps and residuals)
 */
void VisDataSet::initializeGlobalGridders(UvDataGridder* max)
{
  Image image = max->getImage(false);
  
  utilityGridder_.initializeForVis(image);

  // If simulating, set populated indices to all

  if(obsWasSet_)
    utilityGridder_.initializePopulatedIndicesToAll();

  // Piggy-back initializing storage parameters onto this call

  initializeStoreParameters();
}

void VisDataSet::initializeStoreParameters()
{
  TimeVal tVal;
  tVal.setToCurrentTime();

  lastGroup_    = -1;
  lastDate_     = -1;
  lastDataSet_  =  0;
  lastJd_       =  0.0;
  lastBaseline_ =  0;
  writeJd_      =  tVal.getMjd() + 2400000.5;
}

/**.......................................................................
 * Resize internal Dfts to the size needed to grid the data to match the 
 * passed image
 */
void VisDataSet::initializeVisibilityArrays(Image& image)
{
  UvDataGridder* max = 0;
  UvDataGridder* planGridder = 0;

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& group = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < group.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = group.stokesData_[iStokes];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	//------------------------------------------------------------
	// Resize images and dfts to match the passed image
	//------------------------------------------------------------

	freqData.resize(image, obsWasSet_, &planGridder);

	UvDataGridder* curr = &freqData.griddedData_;

	if(max == 0 
	   || curr->xAxis().getNpix() > max->xAxis().getNpix()
	   || curr->yAxis().getNpix() > max->yAxis().getNpix()) {
	  max = curr;
	}

	freqData.iVis_     = 0;
	freqData.estChisq_ = 0.0;
      }
    }
  }

  initializeGlobalGridders(max);
}

/**.......................................................................
 * Compute the primary beam (power pattern) for this baseline
 */
void VisDataSet::computePrimaryBeams()
{
  std::ostringstream os;
  PgUtil::setInteractive(false);

  initWait();

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& group = baselineGroups_[iGroup];

    Antenna& ant1 = group.antennaPair_.first;
    Antenna& ant2 = group.antennaPair_.second;
    
    for(unsigned iStokes=0; iStokes < group.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = group.stokesData_[iStokes];

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	// Set this first, so that freqData.hasData() will return the
	// right thing if we are simulating

	freqData.generatingFakeData_ = obsWasSet_;

	if(freqData.hasData()) {
	  computePrimaryBeamMultiThread(freqData, ant1, ant2, iGroup, iStokes, iFreq);

	  std::cout << "\rComputing beams for group " << group << " frequency = " 
		    << setw(5) << std::setprecision(2) << std::fixed << std::left << std::setfill('0') 
		    << freqData.frequency_.GHz() << " GHz";

	  fflush(stdout);
	}
      }
    }
  }

  COUT("\r                                                                                    ");

  waitUntilDone();
}

/**.......................................................................
 * Compute the synthesized beam for this dataset
 */
void VisDataSet::computeGlobalSynthesizedBeam()
{
  Image image = getImage(ACC_BEAM, IMG_DIRTY);

  Angle fwhmMin, fwhmMaj;
  double fwhmToSigma = 1.0/(2*sqrt(2*log(2)));

  utilityGridder_.getEstimatedSynthesizedBeam(fwhmMaj, fwhmMin, synthBeamRotAngle_);

  COUTCOLOR("Estimated synthesized beam: angle =  " << std::setprecision(2) << std::fixed << synthBeamRotAngle_.degrees()  << " deg,"
	    << " Maj (fwhm) = " << fwhmMaj.arcsec() << "\"," 
	    << " Min (fwhm) = " << fwhmMin.arcsec() << "\"" << std::endl, "cyan");

  //------------------------------------------------------------
  // Store the parameters as sigmas
  //------------------------------------------------------------

  synthBeamMajSig_.setRadians(fwhmMin.radians() * fwhmToSigma);
  synthBeamMinSig_.setRadians(fwhmMaj.radians() * fwhmToSigma);

  SolidAngle sa(fwhmMin, fwhmMaj);

  //------------------------------------------------------------
  // Now for convenience, set this beam as the global synthesized beam
  // for all baseline subsets
  //------------------------------------------------------------

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];
	freqData.estimatedGlobalSynthesizedBeam_ = sa;
      }
    }
  }

  storeSynthesizedBeamModelForPlots(pgManager_);
}

/**.......................................................................
 * Shift the gridded data if requested to do so
 */
void VisDataSet::shiftIfRequested()
{
  if(!shiftRequested_)
    return;

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	freqData.griddedData_.shiftBy(xShift_, yShift_);
      }
    }
  }
}

/**.......................................................................
 * Estimate synthesized beams for all VisFreqData subsets
 */
void VisDataSet::estimateSynthesizedBeams()
{
  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	if(freqData.hasData()) {
	  Angle fwhmMin, fwhmMaj;
	  
	  freqData.griddedData_.getEstimatedSynthesizedBeam(fwhmMaj, fwhmMin, freqData.synthBeamRotAngle_);
	  
	  // Store the parameters as sigmas
	  
	  freqData.synthBeamMajSig_.setRadians(fwhmMin.radians() / 2.35);
	  freqData.synthBeamMinSig_.setRadians(fwhmMaj.radians() / 2.35);
	}

      }
    }
  }
}

/**.......................................................................
 * Multi-thread-aware version of computePrimaryBeam
 */
void VisDataSet::computePrimaryBeamMultiThread(VisFreqData& vfd, Antenna& ant1, Antenna& ant2, 
					       unsigned iGroup, unsigned iStokes, unsigned iFreq)
{
  if(!pool_) {

    //------------------------------------------------------------
    // Compute the weighted sum of the primary beam over all shifted
    // datasets
    //------------------------------------------------------------

    double wtSumTotal = 0;

    for(unsigned i=0; i < vfd.wtSums_.size(); i++) {

      //      COUT(" size = " << vfd.wtSums_.size() << " iGroup = " << iGroup << " freq = " << vfd.frequency_.GHz() << " wt = " << vfd.wtSums_[i] << " pb has data " << vfd.primaryBeam_.hasData());

      if(!vfd.primaryBeam_.hasData()) {

	if(vfd.wtSums_[i] > 0.0) {
	  vfd.primaryBeam_  = ant1.getRealisticApertureField(vfd.primaryBeam_, vfd.frequency_, vfd.xShifts_[i], vfd.yShifts_[i]);
	  vfd.primaryBeam_ *= ant2.getRealisticApertureField(vfd.primaryBeam_, vfd.frequency_, vfd.xShifts_[i], vfd.yShifts_[i]);
	  wtSumTotal = vfd.wtSums_[i];
	}

      } else {

	Image pbPrev = vfd.primaryBeam_;
	Image pbCurr;

	double wtCurr = vfd.wtSums_[i];
	double wtPrev = wtSumTotal;

	pbCurr  = ant1.getRealisticApertureField(vfd.primaryBeam_, vfd.frequency_, vfd.xShifts_[i], vfd.yShifts_[i]);
	pbCurr *= ant1.getRealisticApertureField(vfd.primaryBeam_, vfd.frequency_, vfd.xShifts_[i], vfd.yShifts_[i]);

	pbPrev *= wtPrev;
	pbCurr *= wtCurr;

	vfd.primaryBeam_ = (pbCurr + pbPrev) / (wtCurr + wtPrev);

	wtSumTotal += wtCurr;
      }
    }

    if(debug_) {
      PgUtil::setInteractive(true);
      vfd.primaryBeam_.display();
      COUT("wtSumtotal =  " << wtSumTotal << " cf " << vfd.wtSumTotal_ << std::endl);
    }

  } else {
    VisExecData* ved = vfd.execData_;
    ved->initialize(&ant1, &ant2);
    registerPending(iGroup, iStokes, iFreq);
    pool_->execute(&execComputePrimaryBeam, ved);
  }
}

/**.......................................................................
 * Static method which can be passed to a thread pool, to compute the
 * primary beam for a single VisFreqData object
 */
EXECUTE_FN(VisDataSet::execComputePrimaryBeam)
{
  VisExecData* ved  = (VisExecData*)args;
  VisDataSet*  vds  = ved->vds_;
  VisFreqData* vfd  = ved->vfd_;
  Antenna*     ant1 = ved->ant1_;
  Antenna*     ant2 = ved->ant2_;

  //------------------------------------------------------------
  // Compute the weighted sum of the primary beam over all shifted
  // datasets
  //------------------------------------------------------------
  
  double wtSumTotal = 0;

  for(unsigned i=0; i < vfd->wtSums_.size(); i++) {
    if(i==0) {
      vfd->primaryBeam_  = ant1->getRealisticApertureField(vfd->primaryBeam_, vfd->frequency_, vfd->xShifts_[i], vfd->yShifts_[i]);
      vfd->primaryBeam_ *= ant2->getRealisticApertureField(vfd->primaryBeam_, vfd->frequency_, vfd->xShifts_[i], vfd->yShifts_[i]);
      wtSumTotal = vfd->wtSums_[i];
    } else {
      Image pbPrev = vfd->primaryBeam_;
      Image pbCurr;
      
      double wtCurr = vfd->wtSums_[i];
      double wtPrev = wtSumTotal;
      
      pbCurr  = ant1->getRealisticApertureField(vfd->primaryBeam_, vfd->frequency_, vfd->xShifts_[i], vfd->yShifts_[i]);
      pbCurr *= ant1->getRealisticApertureField(vfd->primaryBeam_, vfd->frequency_, vfd->xShifts_[i], vfd->yShifts_[i]);
      
      pbPrev *= wtPrev;
      pbCurr *= wtCurr;
      
      vfd->primaryBeam_ = (pbCurr + pbPrev) / (wtCurr + wtPrev);
      
      wtSumTotal += wtCurr;
    }
  }
  
  vds->registerDone(ved->iGroup_, ved->iStokes_, ved->iFreq_);
}

/**.......................................................................
 * Perform simulated observations matching the current array
 * configuration and visibilities.  
 *
 * (We decouple transforming the images from computing the
 * visibilities, so that these can be separately multi-threaded)
 */
void VisDataSet::observe()
{
  //------------------------------------------------------------
  // First transform images
  //------------------------------------------------------------

  transformImages();

  //------------------------------------------------------------
  // Then replace existing visibilities with simulated ones
  //------------------------------------------------------------

  replaceVisibilities();
}

/**.......................................................................
 * Perform simulated observations from an ObsInfo object specifying
 * the array parameters.
 */
void VisDataSet::observe(ObsInfo& obs)
{
  try {

    //------------------------------------------------------------
    // First transform images
    //------------------------------------------------------------

    transformImages();

    //------------------------------------------------------------
    // Then simulate visibilities matching the obs specification
    //------------------------------------------------------------

    calculateVisibilities(obs);

  } catch(Exception& err) {
    COUT(err.what());
    throw err;
  }
}

/**.......................................................................
 * Transform all simulated images prior to calculating simulated
 * visibilities
 */
void VisDataSet::transformImages()
{
  //------------------------------------------------------------
  // Iterate over all unique baseline groups, Stokes parameters and
  // frequencies, observing the source
  //------------------------------------------------------------

  initWait();

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];
    
    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];
	transformImageMultiThread(freqData, iGroup, iStokes, iFreq);
      }
    }
  }

  waitUntilDone();
}

/**.......................................................................
 * Initialize waiting for a multi-threaded process to complete.  If not
 * running in a multi-threaded context, this is a no-op
 */
void VisDataSet::initWait()
{
  if(pool_) {
    synchronizer_.reset();
    synchronizer_.initWait();
  }
}

/**.......................................................................
 * Wait until we are signalled that a transaction has completed.  If
 * not running in multi-threaded context, this is a no-op
 */
void VisDataSet::waitUntilDone()
{
  if(pool_) {
    synchronizer_.wait();
  }
}

/**.......................................................................
 * Static method which can be passed to a thread pool, to transform a
 * single VisFreqData image
 */
EXECUTE_FN(VisDataSet::execTransformImage)
{
  VisExecData* ved = (VisExecData*)args;
  VisDataSet*  vds = ved->vds_;
  VisFreqData* vfd = ved->vfd_;

  vfd->transformImage();
  vds->registerDone(ved->iGroup_, ved->iStokes_, ved->iFreq_);
}

/**.......................................................................
 * Multi-thread-aware version of transformImage
 */
void VisDataSet::transformImageMultiThread(VisFreqData& vfd, unsigned iGroup, unsigned iStokes, unsigned iFreq)
{
  if(!pool_) {
    vfd.transformImage();
  } else {
    VisExecData* ved = vfd.execData_;
    registerPending(iGroup, iStokes, iFreq);
    pool_->execute(&execTransformImage, ved);
  }
}

/**.......................................................................
 * Register a single iteration of a transaction as done.  If this iteration
 * completes the entire transaction, then any waiting threads will be signalled
 */
void VisDataSet::registerDone(unsigned iGroup, unsigned iStokes, unsigned iFreq)
{
  unsigned nStokes = obs_.getNumberOfStokesParameters();
  unsigned nFreq   = obs_.getNumberOfFrequencies();
  unsigned nGroup  = baselineGroups_.size();

  unsigned iBit = (iGroup * nStokes + iStokes) * nFreq + iFreq;

  synchronizer_.registerDone(iBit, nGroup * nStokes * nFreq); 
}

/**.......................................................................
 * Register a single iteration of a transaction as pending.
 */
void VisDataSet::registerPending(unsigned iGroup, unsigned iStokes, unsigned iFreq)
{
  unsigned nStokes = obs_.getNumberOfStokesParameters();
  unsigned nFreq   = obs_.getNumberOfFrequencies();
  unsigned nGroup  = baselineGroups_.size();

  unsigned iBit = (iGroup * nStokes + iStokes) * nFreq + iFreq;

  synchronizer_.registerPending(iBit);
}

/**.......................................................................
 * Calculate simulated visibilities from previously-transformed images
 */
void VisDataSet::calculateVisibilities(ObsInfo& obs)
{
  //------------------------------------------------------------
  // Initialize the visibility array
  //------------------------------------------------------------
  
  obs.initializeSimulationVisibilityArray();

  //------------------------------------------------------------
  // Now fill the visibility array
  //------------------------------------------------------------

  fillSimulationVisibilityArray(obs);
}

/**.......................................................................
 * Fill the internal visibility array with simulated values
 */
void VisDataSet::
fillSimulationVisibilityArray(ObsInfo& obs)
{
#ifdef SIM_TIMER_TEST
  Timer t1, t2, t3, t4, t5, t6;
  double t1time = 0.0;
  double t2time = 0.0;
  double t3time = 0.0;
  double t4time = 0.0;
  double t5time = 0.0;
  double t6time = 0.0;
  t1.start();
#endif

  std::ofstream fout;
  if(debug_) {
    fout.open("visout.txt", ios::out);
  }

  //------------------------------------------------------------
  // First, initialize all internal arrays needed for moment
  // calculation.  This is to mimic what is done when we are reading
  // in data, so that the data/model and residuals can be displayed
  // even when simulating data, if desired
  //------------------------------------------------------------

  initializeForMomentAccumulation(true);
  utilityGridder_.initializeForFirstMoments();

  //------------------------------------------------------------
  // Now start simulating data
  //------------------------------------------------------------

  Geoid geoid;
  unsigned iVisGroup = 0;
  Flux noiseRms0, noiseRms, reNoise, imNoise;
  Frequency bw0;

  obs.calculateStartJd();

  //------------------------------------------------------------  
  // Iterate over timestamps
  //------------------------------------------------------------

  HourAngle startHa = obs.getStartHa();
  HourAngle stopHa  = obs.getStopHa();
  HourAngle deltaHa = obs.getDeltaHa();
  HourAngle ha;
  Time dt;

  gcp::util::Lla lla = obs.getArrayLocation();
  Declination dec = obs.obsDec_;
  gcp::util::PolarLengthVector azel;

  dt.setSeconds(deltaHa.seconds() * Constants::utSecPerSiderealSec_);

  unsigned nHa = (unsigned)((stopHa - startHa) / deltaHa);

  for(unsigned iHa=0; iHa < nHa; iHa++) {

    std::cout << "\rComputing iHa = " << iHa << " of " << nHa;
    fflush(stdout);

    //------------------------------------------------------------
    // Compute the current HA
    //------------------------------------------------------------

    ha = (startHa) + (deltaHa * iHa) + (deltaHa/2);

    //------------------------------------------------------------
    // Compute the current Az/El.  Although this is technically
    // different for different antennas, to speed up the calculation
    // we just compute it for the array center and assume all antennas
    // have the same Az/El
    //------------------------------------------------------------

    azel = geoid.geodeticLlaAndHaDecToAzEl(lla, ha, dec);
    
    //------------------------------------------------------------
    // Iterate over all baseline groups for this timestamp
    //------------------------------------------------------------

    for(unsigned iBaseGroup=0; iBaseGroup < baselineGroups_.size(); iBaseGroup++) {

      VisBaselineGroup& group = baselineGroups_[iBaseGroup];
      unsigned nBaseline = group.baselines_.size();

      //------------------------------------------------------------
      // Iterate through all baselines of this group
      //------------------------------------------------------------

      for(unsigned iBase=0; iBase < nBaseline; iBase++) {

	VisBaseline& baseline = group.baselines_[iBase];

	Antenna* ant1 = baseline.ant1_;
	Antenna* ant2 = baseline.ant2_;
	
	//------------------------------------------------------------
	// Store the XYZ coordinates of this baseline
	//------------------------------------------------------------
	
#ifdef SIM_TIMER_TEST
	t4.start();
#endif

	LengthTriplet xyz1 = ant1->getXyz();
	LengthTriplet xyz2 = ant2->getXyz();
	LengthTriplet dxyz = xyz2 - xyz1;

#ifdef SIM_TIMER_TEST
	t4.stop();
	t4time += t4.deltaInSeconds();
#endif
	//------------------------------------------------------------
	// Compute the UVW coordinate that corresponds to the
	// current HA
	//------------------------------------------------------------
	    
	LengthTriplet uvw = geoid.haDecAndXyzToUvw(ha, obs.obsDec_, dxyz);

	//------------------------------------------------------------
	// Iterate over all Stokes parameters for this baseline
	//
	// Store the first bandwidth we encounter
	//------------------------------------------------------------

	bw0 = group.stokesData_[0].freqData_[0].bandwidth_;

	//------------------------------------------------------------
	// And generate the noise rms for this bandwidth -- we can
	// simply scale from this for other bandwidths
	//------------------------------------------------------------

#ifdef SIM_TIMER_TEST
	t2.start();
#endif

	obs.getNoiseRms(noiseRms0, &ha, &obs.obsDec_, &bw0, &dt, ant1, ant2, &azel);

	unsigned nStokes = group.stokesData_.size();
	
#ifdef SIM_TIMER_TEST
	t2.stop();
	t2time += t2.deltaInSeconds();
#endif

	for(unsigned iStokes=0; iStokes < nStokes; iStokes++) {

	  VisStokesData& stokesData = group.stokesData_[iStokes];
	  unsigned nFreq = stokesData.freqData_.size();
      
	  //------------------------------------------------------------
	  // Iterate over all frequencies of this Stokes parameter
	  //------------------------------------------------------------
	  
	  for(unsigned iFreq=0; iFreq < nFreq; iFreq++) {
	    
#ifdef SIM_TIMER_TEST
	    t5.start();
#endif

	    VisFreqData& freqData = stokesData.freqData_[iFreq];
	    Frequency& freq = freqData.frequency_;
	    
	    // Convert to u and v in units of wavelength

	    double u = uvw.u_.meters() / freq.meters();
	    double v = uvw.v_.meters() / freq.meters();
	    
	    // Finally, interpolate the visibility
	    
	    double re=0.0, im=0.0;
	    double reIp=0.0, imIp=0.0;
	    double reFp=0.0, imFp=0.0;
	    bool validIp=true, validFp=true, valid;
	    
	    //------------------------------------------------------------
	    // Interpolate data from the composite Image and
	    // Fourier-plane models
	    //------------------------------------------------------------

	    if(freqData.compositeImageModelDft_.hasData_) 
	      freqData.compositeImageModelDft_.interpolateReImData(u, v, reIp, imIp, validIp);
	    
	    if(freqData.compositeFourierModelDft_.hasData_) 
	      freqData.compositeFourierModelDft_.interpolateReImData(u, v, reFp, imFp, validFp);

#ifdef SIM_TIMER_TEST
	    t5.stop();
	    t5time += t5.deltaInSeconds();

	    t6.start();
#endif
	    re = reIp + reFp;
	    im = imIp + imFp;
	    valid = validIp && validFp;

	    // Now do something with it!
	    
	    ObsInfo::Vis& vis = obs.visibilities_[iVisGroup];

	    // Store the uv as light travel time, in seconds.  We take
	    // the negative, so that when we reverse the sign on
	    // read-in, we recover the correct UV

	    vis.u_ = -uvw.u_.meters() / Constants::lightSpeed_.metersPerSec();
	    vis.v_ = -uvw.v_.meters() / Constants::lightSpeed_.metersPerSec();
	    vis.w_ = -uvw.w_.meters() / Constants::lightSpeed_.metersPerSec();

	    // Install the AIPS baseline code

	    vis.baseline_ = baseline.aipsBaselineIndex_;

	    // And come up with a JD date

	    vis.jd_ = obs.startJd_ + iHa * obs.deltaHa_.hours()/24;
	    
	    // Note that when images were installed, they were already
	    // converted to Jy, so no need to convert here
	    
	    unsigned visInd = iStokes * nFreq + iFreq;
	    
	    double wt = 1.0;

	    //------------------------------------------------------------
	    // Adjust the noise rms to match the current bandwidth.
	    // 
	    // We also allow adjusting by the intScale_ so that we can
	    // effectively simulate multiple observations over the
	    // same HA range, if intScale_ > 1.0
	    //------------------------------------------------------------

	    noiseRms.setJy((noiseRms0.Jy() * sqrt(bw0.MHz() / freqData.bandwidth_.MHz())) / sqrt(intScale_));

	    //------------------------------------------------------------
	    // And generate noise samples from it.  
	    //------------------------------------------------------------

#ifdef SIM_TIMER_TEST
	    t3.start();
#endif

	    obs.generateNoise(noiseRms, reNoise, imNoise, wt);

#ifdef SIM_TIMER_TEST
	    t3.stop();
	    t3time += t3.deltaInSeconds();
#endif

	    vis.re_[visInd] = re + reNoise.Jy();
	    vis.im_[visInd] = im + imNoise.Jy();
	    vis.wt_[visInd] = valid ? wt : 0.0;

	    if(debug_) {
	      fout << u << " " << v << " " << vis.re_[visInd] << " " << vis.im_[visInd] << " " << vis.wt_[visInd] << std::endl;
	    }

	    //------------------------------------------------------------
	    // Accumulate the results into our containers too, so we
	    // can display it if we want to.
	    //------------------------------------------------------------

	    {
	      double accu  = uvw.u_.meters() / freqData.frequency_.meters();
	      double accv  = uvw.v_.meters() / freqData.frequency_.meters();
	      double accwt = vis.wt_[visInd] * wtScale_ * taper(accu, accv);
	      double accre = vis.re_[visInd];
	      double accim = vis.im_[visInd];

	      double r = sqrt(accu * accu + accv * accv);

	      bool goodVis = (uvMin_ < 0.0 || r > uvMin_);
	      goodVis &= (uvMax_ < 0.0 || r < uvMax_);

	      if(goodVis) {
		freqData.griddedData_.accumulateFirstMoments(accu, accv, accre, accim, accwt);
	      }
	    }

#ifdef SIM_TIMER_TEST
	    t6.stop();
	    t6time += t6.deltaInSeconds();
#endif
	  }
	}

	iVisGroup++;
      }
    }
  }

  //------------------------------------------------------------
  // For display of simulated data, call calculateErrorInMean, as we
  // do when reading in data, to correctly initialize populated
  // indices in the gridders, weight sums, etc.
  //------------------------------------------------------------

  calculateErrorInMean();

  COUT("\r                                                                ");
  fflush(stdout);

#ifdef SIM_TIMER_TEST
  t1.stop();
  t1time += t1.deltaInSeconds();

  COUT("T1 = " << t1time);
  COUT("T2 = " << t2time);
  COUT("T3 = " << t3time);
  COUT("T4 = " << t4time);
  COUT("T5 = " << t5time);
  COUT("T6 = " << t6time);

  obs.printTime();
#endif

  if(debug_) {
    fout.close();
  }
}

/**.......................................................................
 * Fill the internal visibility array with simulated values
 */
void VisDataSet::
determineUvMax(ObsInfo& obs)
{
  //------------------------------------------------------------
  // Now start simulating data
  //------------------------------------------------------------

  Geoid geoid;

  obs.calculateStartJd();

  //------------------------------------------------------------  
  // Iterate over timestamps
  //------------------------------------------------------------

  HourAngle startHa = obs.getStartHa();
  HourAngle stopHa  = obs.getStopHa();
  HourAngle deltaHa = obs.getDeltaHa();
  HourAngle ha;
  Time dt;

  dt.setSeconds(deltaHa.seconds() * Constants::utSecPerSiderealSec_);

  unsigned nHa = (unsigned)((stopHa - startHa) / deltaHa);

  for(unsigned iHa=0; iHa < nHa; iHa++) {

    std::cout << "\rComputing iHa = " << iHa << " of " << nHa;
    fflush(stdout);

    //------------------------------------------------------------
    // Compute the current HA
    //------------------------------------------------------------

    ha = (startHa) + (deltaHa * iHa) + (deltaHa/2);

    //------------------------------------------------------------
    // Iterate over all baseline groups for this timestamp
    //------------------------------------------------------------

    for(unsigned iBaseGroup=0; iBaseGroup < baselineGroups_.size(); iBaseGroup++) {

      VisBaselineGroup& group = baselineGroups_[iBaseGroup];
      unsigned nBaseline = group.baselines_.size();

      //------------------------------------------------------------
      // Iterate through all baselines of this group
      //------------------------------------------------------------

      for(unsigned iBase=0; iBase < nBaseline; iBase++) {

	VisBaseline& baseline = group.baselines_[iBase];

	Antenna* ant1 = baseline.ant1_;
	Antenna* ant2 = baseline.ant2_;
	
	//------------------------------------------------------------
	// Store the XYZ coordinates of this baseline
	//------------------------------------------------------------
	
	LengthTriplet xyz1 = ant1->getXyz();
	LengthTriplet xyz2 = ant2->getXyz();
	LengthTriplet dxyz = xyz2 - xyz1;

	//------------------------------------------------------------
	// Compute the UVW coordinate that corresponds to the
	// current HA
	//------------------------------------------------------------
	    
	LengthTriplet uvw = geoid.haDecAndXyzToUvw(ha, obs.obsDec_, dxyz);

	unsigned nStokes = group.stokesData_.size();
	
	for(unsigned iStokes=0; iStokes < nStokes; iStokes++) {

	  VisStokesData& stokesData = group.stokesData_[iStokes];
	  unsigned nFreq = stokesData.freqData_.size();
      
	  //------------------------------------------------------------
	  // Iterate over all frequencies of this Stokes parameter
	  //------------------------------------------------------------
	  
	  for(unsigned iFreq=0; iFreq < nFreq; iFreq++) {
	    
	    VisFreqData& freqData = stokesData.freqData_[iFreq];
	    Frequency& freq = freqData.frequency_;
	    
	    // Convert to u and v in units of wavelength

	    double u = uvw.u_.meters() / freq.meters();
	    double v = uvw.v_.meters() / freq.meters();

	    double uAbs = fabs(u);
	    double vAbs = fabs(v);

	    freqData.uAbsMax_  = freqData.uAbsMax_  > uAbs ? freqData.uAbsMax_     : uAbs;
	    freqData.vAbsMax_  = freqData.vAbsMax_  > vAbs ? freqData.vAbsMax_     : vAbs;
	  }
	}
      }
    }
  }
}

/**.......................................................................
 * Replace internal visibilities with simulated visibilities
 */
void VisDataSet::
replaceVisibilities()
{
  unsigned nGroupsBetweenUpdates = obs_.visibilities_.size() / 10;
  bool first = true;

  Flux reNoise, imNoise;

  //------------------------------------------------------------  
  // Iterate over all groups in the file
  //------------------------------------------------------------

  for(unsigned iGroup=0; iGroup < obs_.visibilities_.size(); iGroup++) {

    // Get the associated visibility

    ObsInfo::Vis& vis = obs_.visibilities_[iGroup];

    // Convert from light-travel time in seconds to meters

    double uMeters = vis.u_ * Constants::lightSpeed_.metersPerSec();
    double vMeters = vis.v_ * Constants::lightSpeed_.metersPerSec();

    // Find the baseline grouping to which this baseline belongs

    unsigned baselineGroupIndex = aipsBaselineIndexToGroupIndexMap_[vis.baseline_];
    VisBaselineGroup& group = baselineGroups_[baselineGroupIndex];

    // Iterate over Stokes and frequency for this visibility

    unsigned nStokes = group.stokesData_.size();
    for(unsigned iStokes=0; iStokes < nStokes; iStokes++) {

      VisStokesData& stokesData = group.stokesData_[iStokes];

      unsigned nFreq = stokesData.freqData_.size();
      for(unsigned iFreq=0; iFreq < nFreq; iFreq++) {
	
	VisFreqData& freqData = stokesData.freqData_[iFreq];
	Frequency& freq = freqData.frequency_;

	// Convert from uv in meters to units of wavelength

	double u = uMeters / freq.meters();
	double v = vMeters / freq.meters();

	// Finally, interpolate the visibility

	double re, im;
	bool valid;
	    
	freqData.compositeImageModelDft_.interpolateReImData(u, v, re, im, valid);

	// Replace the visibility with the interpolated values

	unsigned visInd = iStokes * nFreq + iFreq;
	    
	double wt = 1.0;
	obs_.generateNoise(reNoise, imNoise, wt);

	vis.re_[visInd] = re + reNoise.Jy();
	vis.im_[visInd] = im + imNoise.Jy();
	
	// Only replace visibilities that weren't already flagged

	if(goodWt(vis.wt_[visInd])) {
	  vis.wt_[visInd] = valid ? wt : 0.0;
	}

      }
    }

    if(iGroup % nGroupsBetweenUpdates == 0) {
      std::cout << "\rReplacing visibilities..." << (100*iGroup)/obs_.visibilities_.size() << "%";
      fflush(stdout);
    }

  }

  COUT("\rReplacing visibilities...100%");
}

/**.......................................................................
 * Replace internal visibilities with simulated visibilities
 */
void VisDataSet::
replaceVisibilities(VisDataSet::AccumulatorType type)
{
  if(!storeDataInternally_)
    ThrowSimpleColorError("Visibility data must be stored on read-in to write out data.  Use " << name_ << ".store = true", "red");

  unsigned nGroupsBetweenUpdates = obs_.visibilities_.size() / 10;
  bool first = true;

  Flux reNoise, imNoise;
  
  std::string visType;

  switch (type) {
  case VisDataSet::ACC_MODEL:
    visType = "model";
    break;
  case VisDataSet::ACC_RES:
    visType = "residual";
    break;
  default:
    visType = "data";
    break;
  }

  //------------------------------------------------------------  
  // Iterate over all groups in the file
  //------------------------------------------------------------

  for(unsigned iGroup=0; iGroup < obs_.visibilities_.size(); iGroup++) {

    //------------------------------------------------------------
    // Get the associated visibility
    //------------------------------------------------------------

    ObsInfo::Vis& vis = obs_.visibilities_[iGroup];

    //------------------------------------------------------------
    // Convert from light-travel time in seconds to meters
    //------------------------------------------------------------

    double uMeters = vis.u_ * Constants::lightSpeed_.metersPerSec();
    double vMeters = vis.v_ * Constants::lightSpeed_.metersPerSec();

    //------------------------------------------------------------
    // Find the baseline grouping to which this baseline belongs
    //------------------------------------------------------------

    unsigned baselineGroupIndex = aipsBaselineIndexToGroupIndexMap_[vis.baseline_];
    VisBaselineGroup& group = baselineGroups_[baselineGroupIndex];

    //------------------------------------------------------------
    // Iterate over Stokes and frequency for this visibility
    //------------------------------------------------------------

    unsigned nStokes = group.stokesData_.size();
    for(unsigned iStokes=0; iStokes < nStokes; iStokes++) {

      VisStokesData& stokesData = group.stokesData_[iStokes];

      unsigned nFreq = stokesData.freqData_.size();
      for(unsigned iFreq=0; iFreq < nFreq; iFreq++) {
	
	VisFreqData& freqData = stokesData.freqData_[iFreq];
	Frequency& freq = freqData.frequency_;

	//------------------------------------------------------------
	// Convert from uv in meters to units of wavelength
	//------------------------------------------------------------

	double u = uMeters / freq.meters();
	double v = vMeters / freq.meters();

	//------------------------------------------------------------
	// Finally, interpolate the visibility, replacing the
	// visibility with the interpolated values
	//------------------------------------------------------------

	double re, im, wt;
	bool valid = true;

	unsigned visInd = iStokes * nFreq + iFreq;

	re = vis.re_[visInd];
	im = vis.im_[visInd];
	wt = vis.wt_[visInd];

	switch (type) {
	case VisDataSet::ACC_MODEL:
	  try {
	    freqData.compositeImageModelDft_.interpolateReImData(u, v, re, im, valid);
	  } catch(...) {
	    valid = false;
	  }
	  break;
	case VisDataSet::ACC_RES:
	  try {
	    freqData.compositeImageModelDft_.interpolateReImData(u, v, re, im, valid);
	    re = vis.re_[visInd] - re;
	    im = vis.im_[visInd] - im;
	  } catch(...) {
	    valid = false;
	  }
	  break;
	default:
	  break;
	}

	//------------------------------------------------------------
	// Replace the visibility with the requested values
	//------------------------------------------------------------

	vis.re_[visInd] = re;
	vis.im_[visInd] = im;
	vis.wt_[visInd] = valid ? wt : 0.0;
      }
    }

    if(iGroup % nGroupsBetweenUpdates == 0) {
      std::cout << "\rWriting " << visType << " visibilities..." << (100*iGroup)/obs_.visibilities_.size() << "%";
      fflush(stdout);
    }

  }

  std::cout << "\rWriting " << visType << " visibilities...100%\r                                       \r";
}

/**.......................................................................
 * Calculate the UVW coordinates of the observation
 */
void VisDataSet::calculateSimulatedUvw()
{
  Geoid geoid;

  unsigned nHa   = (unsigned)((obs_.stopHa_ - obs_.startHa_) / obs_.deltaHa_);
  unsigned nAnt  = obs_.antennas_.size();
  unsigned nBase = (nAnt * (nAnt-1))/2;

  std::vector<float> u(nHa * nBase);
  std::vector<float> v(nHa * nBase);

  std::vector<float> up(nHa * nBase);
  std::vector<float> vp(nHa * nBase);
  std::vector<float> mup(nHa * nBase);
  std::vector<float> mvp(nHa * nBase);

  //------------------------------------------------------------
  // First calculate XYZ for all antennas
  //------------------------------------------------------------

  std::vector<LengthTriplet> xyz(obs_.antennas_.size());

  Lla lla = obs_.getArrayLocation();
  for(unsigned iAnt=0; iAnt < obs_.antennas_.size(); iAnt++) {
    LengthTriplet enu = obs_.antennas_[iAnt].getEnu();
    xyz[iAnt] = geoid.geodeticLlaAndEnuToXyz(lla, enu);
  }

  PgUtil::open("2/xs");
  PgUtil::setOverplot(false);
  PgUtil::setXmin(-2e4);
  PgUtil::setXmax( 2e4);
  PgUtil::setYmin(-2e4);
  PgUtil::setYmax( 2e4);
  PgUtil::setUsedefs(true);

  //------------------------------------------------------------
  // Now iterate over timestamps
  //------------------------------------------------------------

  HourAngle ha;
  HourAngle startHa = obs_.getStartHa();
  HourAngle deltaHa = obs_.getDeltaHa();

  for(unsigned iHa=0; iHa < nHa; iHa++) {

    ha = startHa + (deltaHa * iHa);

    for(unsigned iAnt1=0, iBase=0; iAnt1 < obs_.antennas_.size()-1; iAnt1++) {
      LengthTriplet& xyz1 = xyz[iAnt1];
      for(unsigned iAnt2=iAnt1+1; iAnt2 < obs_.antennas_.size(); iAnt2++, iBase++) {
	LengthTriplet& xyz2 = xyz[iAnt2];
	LengthTriplet xyz = xyz2 - xyz1;
	LengthTriplet uvw = geoid.haDecAndXyzToUvw(ha, obs_.getObsDec(), xyz);
	  
	u[iHa * nBase + iBase]  =  uvw.u_.meters();
	v[iHa * nBase + iBase]  =  uvw.v_.meters();
      }
    }
  }

  //------------------------------------------------------------
  // Now display, iterating over frequency
  //------------------------------------------------------------
  
  for(unsigned iFreq=0; iFreq < obs_.nFreq_; iFreq++) {
    Frequency& freq = obs_.frequencies_[iFreq];
    
    for(unsigned iUv=0; iUv < u.size(); iUv++) {

      up[iUv]  =   u[iUv] / freq.meters();
      vp[iUv]  =   v[iUv] / freq.meters();
      mup[iUv] = -up[iUv];
      mvp[iUv] = -vp[iUv];

    }

    PgUtil::linePlot(up.size(), &up[0], &vp[0], 0, "U", "V", "", false, false);
    PgUtil::setOverplot(true);
    PgUtil::linePlot(mup.size(), &mup[0], &mvp[0], 0, "U", "V", "", false, false);
  }
}

/**.......................................................................
 * Simulation only: Add an image to the image to be observed by the
 * array
 */
void VisDataSet::addImage(Image& image, int iFreq, int iStokes)
{
  UvDataGridder* planGridder = 0;

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];
    
    unsigned iStokesStart = (iStokes < 0) ? 0 : (unsigned) iStokes;
    unsigned iStokesStop  = (iStokes < 0) ? groupData.stokesData_.size() : (unsigned) iStokes;

    for(unsigned iStokes=iStokesStart; iStokes < iStokesStop; iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      unsigned iFreqStart = (iFreq < 0) ? 0 : (unsigned) iFreq;
      unsigned iFreqStop  = (iFreq < 0) ? stokesData.freqData_.size() : (unsigned) iFreq;

      for(unsigned iFreq=iFreqStart; iFreq < iFreqStop; iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];
	freqData.addImage(image, &planGridder);
      }
    }
  }
}

/**.......................................................................
 * Return true if we have enough information to simulate an array, and
 * we currently have anything to observe
 */
bool VisDataSet::canSimulate()
{
  return obs_.canSimulate() && haveImages();
}

bool VisDataSet::haveImages()
{
  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];
    
    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      if(!stokesData.hasImage())
	return false;
    }
  }

  return true;
}

//=======================================================================
// Methods of VisDataSet::VisBaselineGroup
//=======================================================================

/**.......................................................................
 * Initialize the internal arrays of this object
 */
void VisDataSet::
VisBaselineGroup::initialize(unsigned nStokes, std::vector<Frequency>& frequencies, std::vector<Frequency>& bandwidths, bool debug)
{
  stokesData_.resize(nStokes);

  for(unsigned iStokes=0; iStokes < stokesData_.size(); iStokes++) {
    VisStokesData& stokesData = stokesData_[iStokes];
    stokesData.stokes_ = Stokes::STOKES_NONE;

    stokesData.freqData_.resize(frequencies.size());

    for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
      VisFreqData& freqData = stokesData.freqData_[iFreq];
      freqData.frequency_ = frequencies[iFreq];
      freqData.ifNo_      = iFreq+1;
      freqData.bandwidth_ = bandwidths[iFreq];
      freqData.execData_  = 0;
      freqData.debug_     = debug;
      freqData.group_     = this;
      freqData.stokes_    = &stokesData;
    }
  }
}

/**.......................................................................
 * Initialize the internal arrays of this object from a per-Stokes map
 */
void VisDataSet::
VisBaselineGroup::initialize(std::map<gcp::util::Stokes::Param, std::map<double, VisDataSet::VisFreqData*> >& stokesMap)
{
  stokesData_.resize(stokesMap.size());

  unsigned iStokes=0;
  for(std::map<gcp::util::Stokes::Param, std::map<double, VisFreqData*> >::iterator stokesIter=stokesMap.begin(); 
      stokesIter != stokesMap.end(); ++stokesIter, ++iStokes) {

    VisStokesData& stokesData = stokesData_[iStokes];
    stokesData.stokes_ = stokesIter->first;

    std::map<double, VisFreqData*>& freqMap = stokesIter->second;

    stokesData.freqData_.resize(freqMap.size());

    unsigned iFreq=0;
    for(std::map<double, VisFreqData*>::iterator freqIter=freqMap.begin(); freqIter != freqMap.end(); ++freqIter, ++iFreq) {

      VisFreqData& freqDataDest = stokesData.freqData_[iFreq];
      VisFreqData* freqDataSrc  = freqIter->second;

      freqDataDest.frequency_ = freqDataSrc->frequency_;
      freqDataDest.ifNo_      = iFreq+1;
      freqDataDest.bandwidth_ = freqDataSrc->bandwidth_;
      freqDataDest.execData_  = 0;
      freqDataDest.debug_     = freqDataSrc->debug_;
      freqDataDest.group_     = this;
      freqDataDest.stokes_    = &stokesData;
      freqDataDest.uAbsMax_   = freqDataSrc->uAbsMax_;
      freqDataDest.vAbsMax_   = freqDataSrc->vAbsMax_;
    }
  }
}

void VisDataSet::
VisBaselineGroup::checkStokesIndex(unsigned iStokes)
{
  if(iStokes > stokesData_.size()-1) {
    ThrowError("Invalid stokes index: " << iStokes << ".  Should be < " << stokesData_.size());
  }
}

void VisDataSet::
VisBaselineGroup::addBaseline(unsigned aipsBaselineIndex, Antenna& ant1, Antenna& ant2)
{
  baselines_.push_back(VisBaseline(aipsBaselineIndex, ant1, ant2));
}

/**.......................................................................
 * Return the approximate PB halfwidth for this antenna pair at the
 * requested frequency
 */
Angle VisDataSet::
VisBaselineGroup::estimatePrimaryBeamHalfWidth(Frequency& freq)
{
  Angle resHw;
  Angle pbSigma1 = Image::gaussianSigma(antennaPair_.first.diameter_, freq);
  Angle pbSigma2 = Image::gaussianSigma(antennaPair_.second.diameter_, freq);

  resHw.setRadians(sqrt(pbSigma1.radians() * pbSigma2.radians()) * 2.35 / 2.0);

  return resHw;
}

/**.......................................................................
 * Return true if this objects represents the same group
 */
bool VisDataSet::VisBaselineGroup::operator==(VisDataSet::VisBaselineGroup& group)
{
  return (group.antennaPair_.first == antennaPair_.first && group.antennaPair_.second == antennaPair_.second) || 
    (group.antennaPair_.first == antennaPair_.second && group.antennaPair_.second == antennaPair_.first);
}

//=======================================================================
// Methods of VisDataSet::VisStokesData
//=======================================================================

void VisDataSet::VisStokesData::checkFrequencyIndex(unsigned iFreq)
{
  if(iFreq > freqData_.size()-1) {
    ThrowError("Invalid frequency index: " << iFreq << ".  Should be < " << freqData_.size());
  }
}

bool VisDataSet::VisStokesData::hasImage()
{
  bool freqsHaveImage = true;
  for(unsigned iFreq=0; iFreq < freqData_.size(); iFreq++) {
    if(!freqData_[iFreq].hasImage()) {
      freqsHaveImage = false;
      break;
    }
  }

  return freqsHaveImage;
}

/**.......................................................................
 * Return true if this objects represents the same Stokes parameter
 */
bool VisDataSet::VisStokesData::operator==(VisDataSet::VisStokesData& stokes)
{
  return stokes_ == stokes.stokes_;
}

std::ostream& gcp::datasets::operator<<(std::ostream& os, const VisDataSet::VisStokesData& stokes)
{
  os << stokes.stokes_;
  return os;
}

//=======================================================================
// Methods of VisDataSet::VisFreqData
//=======================================================================

/**.......................................................................
 * Clear the composite model(s)
 */
void VisDataSet::VisFreqData::clearModel()
{
  if(compositeImageModel_.hasData_) {
    compositeImageModel_.zero();
    compositeImageModel_.hasData_            = false;

    compositeImageModelDft_.zero();
    compositeImageModelDft_.hasData_         = false;
    compositeImageModelDft_.isTransformed_   = false;
  }

  if(compositeFourierModelDft_.hasData_) {
    compositeFourierModelDft_.zero();
    compositeFourierModelDft_.hasData_       = false;
    compositeFourierModelDft_.isTransformed_ = false;
  }
}

/**.......................................................................
 * Add a model component to the composite model
 */
void VisDataSet::VisFreqData::addModel(Generic2DAngularModel& model)
{
  if(hasData()) {
    if(isImagePlaneModel(model)) {
      addImagePlaneModel(model);
    } else {
      addFourierPlaneModel(model);
    }
  }
}

/**.......................................................................
 * Remove the current composite model from the data
 */
void VisDataSet::VisFreqData::remModel()
{
  unsigned dftInd;

  if(hasData()) {
    for(unsigned i=0; i < griddedData_.populatedIndices_.size(); i++) {

      dftInd  = griddedData_.populatedIndices_[i];

      //------------------------------------------------------------
      // The model we compare to is the sum of the composite
      // Image-plane and Fourier-plane models
      //------------------------------------------------------------
      
      double reModel = compositeImageModelDft_.out_[dftInd][0] + compositeFourierModelDft_.out_[dftInd][0];
      double imModel = compositeImageModelDft_.out_[dftInd][1] + compositeFourierModelDft_.out_[dftInd][1];

      griddedData_.out_[dftInd][0] -= reModel;
      griddedData_.out_[dftInd][1] -= imModel;
    }
  }
}

/**.......................................................................
 * Add a Fourier-plane model to this dataset
 */
void VisDataSet::VisFreqData::addFourierPlaneModel(Generic2DAngularModel& model)
{
  //------------------------------------------------------------
  // Load the model component into our temporary array
  //------------------------------------------------------------
  
  gcp::models::PtSrcModel::UvParams params;
  params.beam_ = &primaryBeam_;
  params.freq_ = &frequency_;
  
  model.fillUvData(DataSetType::DATASET_RADIO, fourierModelComponent_, &params);

  //------------------------------------------------------------
  // Now convert to Jy.  If we are fitting components in Jy/bm, we
  // must use a single global beam width, else we would be converting
  // to different intensity units for each VisFreqData set
  //------------------------------------------------------------
  
  fourierModelComponent_.convertToJy(frequency_, estimatedGlobalSynthesizedBeam_);

  //------------------------------------------------------------
  // Finally, add it to the composite model
  //------------------------------------------------------------
  
  if(!compositeFourierModelDft_.hasData()) {
    compositeFourierModelDft_.assignDataFrom(fourierModelComponent_);
  } else {
    compositeFourierModelDft_ += fourierModelComponent_;
  }
}

/**.......................................................................
 * Add an image-plane model to this data set
 */
void VisDataSet::VisFreqData::addImagePlaneModel(Generic2DAngularModel& model)
{
  //------------------------------------------------------------
  // Load the model component into our temporary array
  //------------------------------------------------------------
  
#if 0
  addmodeltimer1.start();
#endif

  model.fillImage(DataSetType::DATASET_RADIO, imageModelComponent_, &frequency_);

#if 0
  addmodeltimer1.stop();
  amt1 += addmodeltimer1.deltaInSeconds();
  addmodeltimer2.start();
#endif

  //------------------------------------------------------------
  // Now convert to Jy.  If we are fitting components in Jy/bm, we
  // must use a single global beam width, else we would be converting
  // to different intensity units for each VisFreqData set
  //------------------------------------------------------------

  imageModelComponent_.convertToJy(frequency_, estimatedGlobalSynthesizedBeam_);

  //------------------------------------------------------------
  // Finally, add it to the composite model
  //------------------------------------------------------------
  
  if(!compositeImageModel_.hasData()) {
    compositeImageModel_.assignDataFrom(imageModelComponent_);
  } else {
    compositeImageModel_ += imageModelComponent_;
  }
}

/**.......................................................................
 * Simulation only: add an image to the image to be observed
 */
void VisDataSet::VisFreqData::addImage(Image& image, UvDataGridder** planGridder)
{
  // If no image has been installed, initialize containers to match
  // the new image size

  if(!hasImage_) {
    resize(image, true, planGridder);
  } else if(!primaryBeam_.axesAreEquivalent(image)) {
    ThrowError("Attempt to add an Image that is a different size than the previously added image");
  }

  // Initialize the temporary model component from the image

  imageModelComponent_ = image;

  //------------------------------------------------------------
  // Now convert to Jy.  If we are fitting components in Jy/bm, we
  // must use a single global beam width, else we would be converting
  // to different intensity units for each VisFreqData set
  //------------------------------------------------------------
    
  imageModelComponent_.convertToJy(frequency_, estimatedGlobalSynthesizedBeam_);

  // And add it to the composite model.  If no image has been added,
  // set the composite model equal to this component.  If an image has
  // already been added, add the image to what's already there.
    
  if(hasImage_)
    compositeImageModel_ += imageModelComponent_;
  else
    compositeImageModel_.assignDataFrom(imageModelComponent_);

  hasImage_ = true;
}

bool VisDataSet::VisFreqData::hasImage()
{
  return hasImage_;
}

/**.......................................................................
 * Transform this VisFreqData's composite image-plane model (if it
 * needs to be transformed)
 */
void VisDataSet::VisFreqData::transformModel()
{
  if(hasData() && compositeImageModel_.hasData() && !compositeImageModelDft_.isTransformed_) {

    //------------------------------------------------------------
    // The composite model should already be in units of Jy.  All we
    // have to do is apply the primary beam here.
    //------------------------------------------------------------
    
    compositeImageModel_ *= primaryBeam_;
    
#if DO_INTERP
    //------------------------------------------------------------
    // We multiply the image by a function that mimicks the data
    // convolution in the Fourier plane
    //------------------------------------------------------------
    
    unsigned nx = compositeImageModel_.xAxis().getNpix();
    unsigned ny = compositeImageModel_.yAxis().getNpix();
    
    double sigx = (double)(nx)/(2*M_PI*Dft2d::convSigInPixels_);
    double sigy = (double)(ny)/(2*M_PI*Dft2d::convSigInPixels_);
    
    Image convCorrection;
    convCorrection.createGaussianImage(nx, ny, sigx, sigy);
    compositeImageModel_ *= convCorrection;
#endif

    //------------------------------------------------------------
    // And transform the image
    //------------------------------------------------------------
    
    compositeImageModelDft_.setInput(compositeImageModel_);
    compositeImageModelDft_.computeForwardTransform();
    compositeImageModelDft_.shift();

    compositeImageModelDft_.isTransformed_ = true;
  }
}

/**.......................................................................
 * Compute chi-square for a single frequency
 */
ChisqVariate VisDataSet::VisFreqData::computeChisq()
{
  bool first = true;

  //------------------------------------------------------------
  // Iterate over just the Fourier components that were populated with
  // data, computing chi-squared from the contributions of both real
  // and imaginary components
  //------------------------------------------------------------
    
  ChisqVariate chisq;
    
  unsigned dftInd;
  double reData, reErr, imData, imErr, reModel, imModel;
  double reCont, imCont;

  for(unsigned i=0; i < griddedData_.populatedIndices_.size(); i++) {
      
#ifdef TIMER_TEST
    double chisqtest=0.0;
    cc1.start();
#endif

    dftInd  = griddedData_.populatedIndices_[i];
 
    reData  = griddedData_.out_[dftInd][0];
    imData  = griddedData_.out_[dftInd][1];
 
    reErr   = griddedData_.errorInMean_[dftInd][0];
    imErr   = griddedData_.errorInMean_[dftInd][1];
      
    //------------------------------------------------------------
    // The model we compare to is the sum of the composite Image-plane
    // and Fourier-plane models
    //------------------------------------------------------------
      
    reModel = compositeImageModelDft_.out_[dftInd][0] + compositeFourierModelDft_.out_[dftInd][0];
    imModel = compositeImageModelDft_.out_[dftInd][1] + compositeFourierModelDft_.out_[dftInd][1];
 
    //------------------------------------------------------------
    // Get the contribution to chisq of the real data
    //------------------------------------------------------------
      
    reCont = (reData - reModel) / reErr;
    reCont *= reCont;
      
    //------------------------------------------------------------
    // And the contribution to chisq of the imaginary data
    //------------------------------------------------------------
      
    imCont = (imData - imModel) / imErr;
    imCont *= imCont;

#ifdef TIMER_TEST
    cc1.stop();
    cct1 += cc1.deltaInSeconds();
    cc2.start();
#endif

    //------------------------------------------------------------
    // And co-add both contributions to chi-square (but do this
    // separately for each component so that chisq registers the
    // correct number of degrees of freedom)
    //------------------------------------------------------------
      
    chisq.directAdd(reCont+imCont, 2);

#ifdef TIMER_TEST
    cc2.stop();
    cct2 += cc2.deltaInSeconds();
    cc3.start();
#endif

#ifdef TIMER_TEST
    chisqtest += reCont;
    chisqtest += imCont;
#endif

#ifdef TIMER_TEST
    cc3.stop();
    cct3 += cc3.deltaInSeconds();
#endif
  }

  return chisq;
}

/**.......................................................................
 * Store the weight sum and shifts of the last dataset that was added
 * to this object.  These will be used when creating primary beams for
 * this object.
 */
void VisDataSet::VisFreqData::storeWtSum(gcp::util::Angle& xShift, gcp::util::Angle& yShift)
{
  if(wtSums_.size() == 0) {
    wtSums_.push_back(wtSumTotal_);
  } else {
    double wtSum = 0;
    for(unsigned i=0; i < wtSums_.size(); i++)
      wtSum += wtSums_[i];
    wtSums_.push_back(wtSumTotal_ - wtSum);
  }

  xShifts_.push_back(xShift);
  yShifts_.push_back(yShift);
}

/**.......................................................................
 * Accumulate visibility data into this object
 */
void VisDataSet::VisFreqData::accumulateMoments(bool first, VisDataSet::VisData& data, gcp::util::Angle& xShift, gcp::util::Angle& yShift)
{
  double xRad = xShift.radians();
  double yRad = yShift.radians();

  double cosFac = cos(-2*M_PI*(data.u_ * xRad + data.v_ * yRad));
  double sinFac = sin(-2*M_PI*(data.u_ * xRad + data.v_ * yRad));

  double re = data.re_ * cosFac - data.im_ * sinFac;
  double im = data.re_ * sinFac + data.im_ * cosFac;

  //------------------------------------------------------------
  // Load this data point into the data gridder for the
  // destination VisFreqData object
  //------------------------------------------------------------
  
  if(first) {
    griddedData_.accumulateFirstMoments( data.u_, data.v_, re, im, data.wt_);
  } else {
    griddedData_.accumulateSecondMoments(data.u_, data.v_, re, im, data.wt_);
  }
  
  accumulateVarianceStats(data);
}

/**.......................................................................
 * Accumulate the estimated chisq and wtsum for this VisFreqData
 * object
 */
void VisDataSet::VisFreqData::accumulateVarianceStats(VisDataSet::VisData& data)
{
  iVis_++;
  nVisUsed_++;
  
  double reCurr = data.re_ - reMean_;
  double imCurr = data.im_ - imMean_;
  
  reCurr *= reCurr;
  imCurr *= imCurr;
  
  double chisq;
  
  chisq = reCurr * data.wt_;
  estChisq_ += (chisq - estChisq_) / iVis_;
  
  chisq = imCurr * data.wt_;
  estChisq_ += (chisq - estChisq_) / iVis_;
  
  wtSumTotal_ += data.wt_;
}

/**.......................................................................
 * Return true if this object represents the same frequency
 */
bool VisDataSet::VisFreqData::operator==(VisDataSet::VisFreqData& freq)
{
  double dFreq = fabs(frequency_.GHz() - freq.frequency_.GHz())/frequency_.GHz();
  return dFreq < 0.01;
}

std::ostream& gcp::datasets::operator<<(std::ostream& os, const VisDataSet::VisFreqData& freq)
{
  return operator<<(os, (VisDataSet::VisFreqData&) freq);
}

std::ostream& gcp::datasets::operator<<(std::ostream& os, VisDataSet::VisFreqData& freq)
{
  os << freq.frequency_.GHz() << " GHz";
  return os;
}

/**.......................................................................
 * Merge two VisFreqData objects
 */
void VisDataSet::VisFreqData::mergeData(VisDataSet::VisFreqData& freq, 
					gcp::util::Angle& xShift, gcp::util::Angle& yShift)
{
  //------------------------------------------------------------
  // Co-add the visibility data
  //------------------------------------------------------------

  griddedData_.mergeData(freq.griddedData_);

  //------------------------------------------------------------
  // Form a weighted mean of the primary beams
  //------------------------------------------------------------

  double wt1 =      wtSumTotal_;
  double wt2 = freq.wtSumTotal_;

  Image pb1  =      primaryBeam_*wt1;
  Image pb2  = freq.primaryBeam_*wt2;

  //------------------------------------------------------------
  // If the existing images don't match, we have to recalculate the
  // (possibly shifted) beam onto our grid
  //------------------------------------------------------------

  if(!pb1.axesAreEquivalent(pb2)) {

    ReportSimpleColorError("Images are on different grids... recomputing for this grid", "yellow");

    pb2   = group_->antennaPair_.first. getRealisticApertureField(pb1, frequency_, xShift, yShift);
    pb2  *= group_->antennaPair_.second.getRealisticApertureField(pb1, frequency_, xShift, yShift);

    pb2 *= wt2;
  }

  primaryBeam_ = (pb1 + pb2) / (wt1 + wt2);

  //------------------------------------------------------------
  // Update internal parameters of this object
  //------------------------------------------------------------

  wtSumTotal_ = wt1 + wt2;

  nVis_       =     nVis_ + freq.nVis_;
  nVisUsed_   = nVisUsed_ + freq.nVisUsed_;

  uvrMax_     =  uvrMax_ > freq.uvrMax_  ?  uvrMax_ : freq.uvrMax_;
  uAbsMax_    = uAbsMax_ > freq.uAbsMax_ ? uAbsMax_ : freq.uAbsMax_;
  vAbsMax_    = vAbsMax_ > freq.vAbsMax_ ? vAbsMax_ : freq.vAbsMax_;

}

void VisDataSet::VisFreqData::operator=(VisDataSet::VisFreqData& data) 
{
  compositeImageModel_      = data.compositeImageModel_;

  imageModelComponent_      = data.imageModelComponent_;
  compositeImageModelDft_   = data.compositeImageModelDft_;

  fourierModelComponent_    = data.fourierModelComponent_;
  compositeFourierModelDft_ = data.compositeFourierModelDft_;

  griddedData_              = data.griddedData_;
  utilityGridder_           = data.utilityGridder_;

  primaryBeam_              = data.primaryBeam_;
  frequency_                = data.frequency_;
  ifNo_                     = data.ifNo_;
  bandwidth_                = data.bandwidth_;
  primaryBeamHalfWidth_     = data.primaryBeamHalfWidth_;

  uvrMax_                   = data.uvrMax_;

  uvrMax_                   = data.uvrMax_;
  uAbsMax_                  = data.uAbsMax_;
  vAbsMax_                  = data.vAbsMax_;
  reMean_                   = data.reMean_;
  imMean_                   = data.imMean_;
  estChisq_                 = data.estChisq_;
  wtScale_                  = data.wtScale_;
  wtSumTotal_               = data.wtSumTotal_;

  nVis_                     = data.nVis_;
  iVis_                     = data.iVis_;
  nVisUsed_                 = data.nVisUsed_;

  group_                    = 0;

  hasImage_                 = data.hasImage_;
  generatingFakeData_       = data.generatingFakeData_;

  execData_                 = 0;

  estimatedGlobalSynthesizedBeam_ = data.estimatedGlobalSynthesizedBeam_;

  synthBeamMajSig_          = data.synthBeamMajSig_;
  synthBeamMinSig_          = data.synthBeamMinSig_;
  synthBeamRotAngle_        = data.synthBeamRotAngle_;

  debug_                    = data.debug_;
}

void VisDataSet::VisFreqData::duplicate(VisDataSet::VisFreqData& data) 
{
  operator=(data);

  fourierModelComponent_.sizeToMatch(data.fourierModelComponent_);
  compositeFourierModelDft_.sizeToMatch(data.compositeFourierModelDft_);

  griddedData_.duplicate(data.griddedData_);
  utilityGridder_.sizeToMatch(data.utilityGridder_);
}

/**.......................................................................
 * Debugging use only: display primary beams
 */
void VisDataSet::displayPrimaryBeams()
{
  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	freqData.primaryBeam_.display();
      }
    }
  }
}

/**.......................................................................
 * Return an image of the model for this dataset, convolved with a
 * Gaussian approximation of the synthesized beam.
 */
Image VisDataSet::getCleanImage(VisDataSet::AccumulatorType type)
{
  Image ret;
  bool first=true;
  double wt, wtSum = 0.0;

  transformModels();

  //------------------------------------------------------------
  // Iterate over VisFreqData sets, forming the weighted sum of the
  // requested image
  //------------------------------------------------------------

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	if(freqData.hasData()) {

	  wt = freqData.griddedData_.wtSumTotal_;

	  Image curr = freqData.getCleanImage(type);

	  if(first) {
	    ret = curr;
	    ret.zero();
	    first = false;
	  } else {

	    if(!(ret.xAxis() == curr.xAxis() && ret.yAxis() == curr.yAxis()))
	      ThrowSimpleColorError("Image dimensions don't match", "red");
	  }

	  for(unsigned i=0; i < ret.data_.size(); i++)
	    ret.data_[i] += (curr.data_[i] - ret.data_[i]) * wt / (wtSum + wt);

	  wtSum += wt;
	}
      }
    }
  }

  return ret;
}

/**.......................................................................
 * Return an image of the model for this dataset, convolved with a
 * gaussian approximation of the synthesized beam.
 */
Image VisDataSet::VisFreqData::getCleanImage(VisDataSet::AccumulatorType type)
{
  utilityGridder_.initializeForFirstMoments();
  accumulate(utilityGridder_, type);
  utilityGridder_.shift();
  utilityGridder_.computeInverseTransform();

  return utilityGridder_.getImage();
}

/**.......................................................................
 * Accumulate data for a single frequency
 */
void VisDataSet::VisFreqData::accumulate(UvDataGridder& gridder, VisDataSet::AccumulatorType type)
{
  if(type == VisDataSet::ACC_CLEAN || type == VisDataSet::ACC_CLEANBEAM || type == VisDataSet::ACC_CLEANMODEL) {
    return accumulateClean(gridder, type);
  } else
    return accumulateDirty(gridder, type);
}

/**.......................................................................
 * Accumulate data for a single frequency
 */
void VisDataSet::VisFreqData::accumulateDirty(UvDataGridder& gridder, VisDataSet::AccumulatorType type)
{
  //------------------------------------------------------------
  // Now iterate over just the Fourier components that were populated
  // with data, accumulating whatever was requested.
  //------------------------------------------------------------

  unsigned dftInd;
  double reData, imData, err, reModel, imModel, wt;

  for(unsigned i=0; i < griddedData_.populatedIndices_.size(); i++) {

    dftInd  = griddedData_.populatedIndices_[i];

    //------------------------------------------------------------
    // Get the data
    //------------------------------------------------------------

    reData  = griddedData_.out_[dftInd][0];
    imData  = griddedData_.out_[dftInd][1];

    //------------------------------------------------------------
    // Use real weights from the data
    //------------------------------------------------------------

    err = griddedData_.errorInMean_[dftInd][0];
    wt = 1.0/(err*err);

    //------------------------------------------------------------
    // Get the model component for this index
    //------------------------------------------------------------

    reModel = compositeImageModelDft_.out_[dftInd][0] + compositeFourierModelDft_.out_[dftInd][0];
    imModel = compositeImageModelDft_.out_[dftInd][1] + compositeFourierModelDft_.out_[dftInd][1];

    //------------------------------------------------------------
    // Construct a running mean of the model by co-adding this data to
    // any data that already exists for this uv cell
    //------------------------------------------------------------

    switch (type) {
    case VisDataSet::ACC_DATA:
      gridder.accumulateFirstMoments(dftInd, reData, imData, wt);
      break;
    case VisDataSet::ACC_RES:
      gridder.accumulateFirstMoments(dftInd, reData - reModel, imData - imModel, wt);
      break;
    case VisDataSet::ACC_MODEL:
      gridder.accumulateFirstMoments(dftInd, reModel, imModel, wt);
      break;
    case VisDataSet::ACC_BEAM:
      gridder.accumulateFirstMoments(dftInd, 1.0, 0.0, wt);
      break;
    default:
      break;
    }
  }
}

/**.......................................................................
 * Accumulate the transform of the clean model for a single frequency
 */
void VisDataSet::VisFreqData::accumulateClean(gcp::util::UvDataGridder& gridder, AccumulatorType type)
{
  //------------------------------------------------------------
  // For the clean model, we are approximating a filled aperture, so
  // we iterate over all indices, tapering by an approximation of the
  // synthesized beam for this VisFreqData set
  //------------------------------------------------------------

  unsigned dftInd;
  double reData, imData, err, reModel, imModel;

  //------------------------------------------------------------
  // Precompute parameters of the synthesized beam
  //------------------------------------------------------------

  double cRotAng = cos(synthBeamRotAngle_.radians());
  double sRotAng = sin(synthBeamRotAngle_.radians());
  double invMajSigma = 1.0/(2*M_PI*synthBeamMajSig_.radians());
  double invMinSigma = 1.0/(2*M_PI*synthBeamMinSig_.radians());

  double u,v,ur,vr,val;
 
  for(unsigned dftInd=0; dftInd < griddedData_.nOutZeroPad_; dftInd++) {

    //------------------------------------------------------------
    // Get the model component for this index
    //------------------------------------------------------------

    reModel = compositeImageModelDft_.out_[dftInd][0] + compositeFourierModelDft_.out_[dftInd][0];
    imModel = compositeImageModelDft_.out_[dftInd][1] + compositeFourierModelDft_.out_[dftInd][1];

    //------------------------------------------------------------
    // Get the UV coordinate of this point
    //------------------------------------------------------------

    griddedData_.getUVData(dftInd, Dft2d::DATA_UV, u, v, val);

    //------------------------------------------------------------
    // Construct the rotated coordinate
    //------------------------------------------------------------

    ur =   u * cRotAng + v * sRotAng;
    vr = - u * sRotAng + v * cRotAng;

    double fac = ur*ur/(2*invMajSigma*invMajSigma) + vr*vr/(2*invMinSigma*invMinSigma);
    fac = exp(-fac);

    //------------------------------------------------------------
    // Construct a running mean of the model by co-adding this data to
    // any data that already exists for this uv cell
    //------------------------------------------------------------

    switch(type) {
    case ACC_CLEANBEAM:
      gridder.accumulateFirstMoments(dftInd, 1.0, 0.0, fac);
      break;
    case ACC_CLEANMODEL:
      gridder.accumulateFirstMoments(dftInd, reModel, imModel, 1.0);
      break;
    default:
      gridder.accumulateFirstMoments(dftInd, reModel, imModel, fac);
      break;
    }
  }

  //------------------------------------------------------------
  // Check if any weights are zero for a dft pixel -- the values will
  // be nans and should be set to zero
  //------------------------------------------------------------

  for(unsigned dftInd=0; dftInd < griddedData_.nOutZeroPad_; dftInd++) {

    if(!(gridder.wtSum_[dftInd] > 0.0)) {
      gridder.out_[dftInd][0] = 0.0;
      gridder.out_[dftInd][1] = 0.0;
    }

  }
}

/**.......................................................................
 * Take the composite image, multiply by the primary beam, and transform
 * it
 */
void VisDataSet::VisFreqData::inverseTransformData()
{
  compositeImageModelDft_.computeInverseTransform();
}

/**.......................................................................
 * Take the composite image, multiply by the primary beam, and transform
 * it
 */
void VisDataSet::VisFreqData::transformImage()
{
  if(compositeImageModel_.hasData()) {

    //------------------------------------------------------------
    // The composite model should already be in units of Jy.  All we
    // have to do is apply the primary beam...
    //------------------------------------------------------------
    
    compositeImageModel_ *= primaryBeam_;
    
    //------------------------------------------------------------
    // We divide the image by a function that corrects for the
    // convolution in the Fourier plane
    //------------------------------------------------------------
    
    unsigned nx = compositeImageModel_.xAxis().getNpix();
    unsigned ny = compositeImageModel_.yAxis().getNpix();
    
    double convSigInPixels = Dft2d::convSigInPixels_;

    double sigx = (double)(nx)/(2*M_PI*convSigInPixels);
    double sigy = (double)(ny)/(2*M_PI*convSigInPixels);
    
    Image convCorrection;
    convCorrection.createGaussianImage(nx, ny, sigx, sigy);
    
    compositeImageModel_ /= convCorrection;
    
    //------------------------------------------------------------
    // Finally, transform the image
    //------------------------------------------------------------
    
    compositeImageModelDft_.zeropad(true);
    compositeImageModelDft_.initialize(compositeImageModel_);
    compositeImageModelDft_.computeForwardTransform();
    compositeImageModelDft_.shift();

    compositeImageModelDft_.isTransformed_ = true;

  }
}

/**.......................................................................
 * Resize images and dfts to the requested size
 */
void VisDataSet::VisFreqData::resize(double percentCorrelation, double correlationLength, bool isSim)
{
  // Find the power-of-2 size of the array that will grid the data _at
  // least_ this finely.  We divide by correlationLength/sqrt(2) since
  // the diagonal is the widest separation in UV we will tolerate.
  
  unsigned nx = 
    Dft2d::nearestPowerOf2NotLessThan(sqrt(2.0) * uAbsMax_ / correlationLength);
  
  unsigned ny = 
    Dft2d::nearestPowerOf2NotLessThan(sqrt(2.0) * vAbsMax_ / correlationLength);

  // First resize the data grid and model dft to match.  n/4 instead
  // of n/2 to account for zeropadding

  Angle xSize(Angle::Radians(), (double)(nx/4) / uAbsMax_);
  Angle ySize(Angle::Radians(), (double)(ny/4) / vAbsMax_);

  griddedData_.initializeForVis(xSize, ySize, nx, ny);
  utilityGridder_.initializeForVis(xSize, ySize, nx, ny);

  if(debug_) {
    COUT("Correlation percentage of: " << percentCorrelation * 100 << "% for group " << *group_
	 << " corresponds to correlation length = " << setw(8) << setprecision(2) 
	 << std::fixed << std::right << correlationLength << " at freq " << frequency_.GHz() << " GHz. " 
	 << "Data will be gridded to: " << setw(4) << std::right << nx << " x " << setw(4) << std::right << ny
	 << " du = " << (uAbsMax_ / nx) << " dv = " << (vAbsMax_ / ny));
  }
  
  //------------------------------------------------------------
  // Resize components needed for manipulating image models
  //------------------------------------------------------------
  
  // Resize the primary beam to match

  primaryBeam_.initialize(xSize, ySize, nx, ny);
  imageModelComponent_.initialize(xSize, ySize, nx, ny);
  compositeImageModel_.initialize(xSize, ySize, nx, ny);

  compositeImageModelDft_.initializeForVis(xSize, ySize, nx, ny);

  //------------------------------------------------------------
  // Resize components needed for manipulating fourier-plane models
  //------------------------------------------------------------

  fourierModelComponent_.initializeForVis(xSize, ySize, nx, ny);
  compositeFourierModelDft_.initializeForVis(xSize, ySize, nx, ny);

  //------------------------------------------------------------
  // If simulating, set populated indices to all
  //------------------------------------------------------------

  if(isSim) {
    fourierModelComponent_.initializePopulatedIndicesToAll();
    compositeFourierModelDft_.initializePopulatedIndicesToAll();
  }
}

/**.......................................................................
 * Resize images and dfts to the requested size
 */
void VisDataSet::VisFreqData::resize(Image& image, bool isSim, UvDataGridder** planGridder)
{
  //------------------------------------------------------------
  // First resize the data grid to match
  //------------------------------------------------------------
  
  if(*planGridder == 0) {
    griddedData_.initializeForVis(image);
    *planGridder = &griddedData_;
  } else {
    griddedData_.setPlan((*planGridder)->forwardPlan_, (*planGridder)->inversePlan_);
    griddedData_.initializeForVis(image);
  }

  utilityGridder_.setPlan((*planGridder)->forwardPlan_, (*planGridder)->inversePlan_);
  utilityGridder_.initializeForVis(image);

  //------------------------------------------------------------
  // Resize components needed for manipulating image models
  //------------------------------------------------------------
  
  // Resize the primary beam to match

  primaryBeam_.initialize(image);
  imageModelComponent_.initialize(image);
  compositeImageModel_.initialize(image);

  compositeImageModelDft_.setPlan((*planGridder)->forwardPlan_, (*planGridder)->inversePlan_);
  compositeImageModelDft_.initializeForVis(image);

  //------------------------------------------------------------
  // Resize components needed for manipulating fourier-plane models
  //------------------------------------------------------------

  fourierModelComponent_.setPlan((*planGridder)->forwardPlan_, (*planGridder)->inversePlan_);
  fourierModelComponent_.initializeForVis(image);

  compositeFourierModelDft_.setPlan((*planGridder)->forwardPlan_, (*planGridder)->inversePlan_);
  compositeFourierModelDft_.initializeForVis(image);

  //------------------------------------------------------------
  // If simulating, set populated indices to all
  //------------------------------------------------------------

  if(isSim) {
    fourierModelComponent_.initializePopulatedIndicesToAll();
    compositeFourierModelDft_.initializePopulatedIndicesToAll();
  }
}

//=======================================================================
// Methods of VisDataSet
//=======================================================================

void VisDataSet::estimateErrorInMeanFromData(bool estimate)
{
  estimateErrInMeanFromData_ = estimate;
}

//-----------------------------------------------------------------------
// Methods for writing FITS files
//-----------------------------------------------------------------------

/**.......................................................................
 * Write data in the visibilities_ array to a FITS file
 */
void VisDataSet::writeUvfFile(std::string fileName)
{
  FitsIoHandler fitsio;
  fitsio.setFirstTelescopeNum(1);
  fitsio.writeUvfFile(fileName, obs_);
}

bool VisDataSet::goodWt(double wt)
{
  return isfinite(wt) && wt > 0.0;
}

//-----------------------------------------------------------------------
// Methods for accumulating and displaying images from this data set
//-----------------------------------------------------------------------

/**.......................................................................
 * Display the data.  For visibility data, we display the dirty map of
 * the composite data (all frequencies, all baselines)
 */
void VisDataSet::display()
{
  PgUtil::setInteractive(interactive_);
  display(ACC_DATA);
}

/**.......................................................................
 * Display the synthesized beam.  For visibility data, we display the
 * dirty map of the composite data (all frequencies, all baselines)
 */
void VisDataSet::displayBeam()
{
  PgUtil::setInteractive(interactive_);
  display(ACC_BEAM);
}

/**.......................................................................
 * Display the composite model.  For visibility data we will transform
 * to the Fourier plane
 */
void VisDataSet::displayCompositeModel()
{
  PgUtil::setInteractive(interactive_);

  bool doClean = false;

  if(getParameter("clean", false)->data_.hasValue())
    doClean = getBoolVal("clean");

  if(doClean) {

    if(getStringVal("cleantype") == "model") {
      display(VisDataSet::ACC_CLEAN);
    } else if(getStringVal("cleantype") == "delta") {
      clean();
    } else {
      ThrowSimpleColorError("Unrecognized cleantype: " << getStringVal("cleantype"), "red");
    }

  } else {
    display(ACC_MODEL);
  }
}

/**.......................................................................
 * Display residuals.  We display the dirty map of the composite data
 * (all frequencies, all baselines).
 */
void VisDataSet::displayResiduals()
{
  PgUtil::setInteractive(interactive_);
  display(ACC_RES);
}

/**.......................................................................
 * Define what it means to display data for a visibility data set
 */
void VisDataSet::display(VisDataSet::AccumulatorType type)
{
  //------------------------------------------------------------
  // Get the requested image
  //------------------------------------------------------------

  Image image;

  if(type == ACC_CLEAN) {
    image     = getCleanImage(ACC_CLEAN);
    Image res = getImage(ACC_RES);
    image    += res;
  } else if(type == ACC_CLEANMODEL) {
    image = getCleanImage(ACC_CLEANMODEL);
  } else if(type == ACC_CLEANBEAM) {
    image = getCleanImage(ACC_CLEANBEAM);
  } else {
    image = getImage(type, IMG_DIRTY);
  }

  //------------------------------------------------------------
  // Adjust the window to display the real data aspect ratio
  //------------------------------------------------------------

  PgUtil::setInteractive(interactive_);
  PgUtil::setWnad(true);
  PgUtil::insertPgManager(pgManager_);

  //------------------------------------------------------------
  // If a colormap was specified, set it now
  //------------------------------------------------------------

  if(getParameter("cmap", false)->data_.hasValue()) {
    PgUtil::setColormap(getStringVal("cmap"));
  } else {
    PgUtil::setColormap("grey");
  }

  //------------------------------------------------------------
  // Default to displaying the full data range, unless a range was
  // specified.  If displaying the beam, ignore user-specified range
  //------------------------------------------------------------

  if(type == ACC_DATA) {
    zmin_ = image.min();
    zmax_ = image.max();
  }

  if(type != ACC_BEAM) {
    if(getParameter("zmin", false)->data_.hasValue() && getParameter("zmax", false)->data_.hasValue()) {
      zmin_ = getDoubleVal("zmin");
      zmax_ = getDoubleVal("zmax");
    }
  } else {
    zmin_ = 0.0;
    zmax_ = 0.0;
    image.units_    = Unit::UNITS_NONE;
    image.hasUnits_ = false;
  }

  PgUtil::setZmin(zmin_);
  PgUtil::setZmax(zmax_);

  VisDataSet::setWedge(type, zmin_, zmax_, getParameter("zmin", false)->data_.hasValue() && getParameter("zmax", false)->data_.hasValue());

  //------------------------------------------------------------
  // Finally, display it
  //------------------------------------------------------------

  PgUtil::useHeader(true);
  PgUtil::setHeader(displayHeaderString(type), PgUtil::JUST_LEFT);

  if(reverseDisplay_) {
    image.display();
  } else {
    image.difmapDisplay();
  }

  PgUtil::useHeader(false);
  PgUtil::clearPgManager();

  //------------------------------------------------------------
  // If requested, write out the image too
  //------------------------------------------------------------

  switch(type) {
  case VisDataSet::ACC_DATA:
    if(getParameter("dataimage", false)->data_.hasValue()) {
      image.writeToFitsFile(getStringVal("dataimage"));
    }
    break;
  case VisDataSet::ACC_MODEL:
  case VisDataSet::ACC_CLEAN:
    if(getParameter("modelimage", false)->data_.hasValue()) {
      image.writeToFitsFile(getStringVal("modelimage"));
    }
    break;
  case VisDataSet::ACC_RES:
    if(getParameter("resimage", false)->data_.hasValue()) {
      image.writeToFitsFile(getStringVal("resimage"));
    }
    break;
  default:
    break;
  }
}

/**.......................................................................
 * Write out data as requsted
 */
void VisDataSet::writeData()
{
  transformModels();

  if(getParameter("datauvf", false)->data_.hasValue())
    writeDataToFile(getStringVal("datauvf"), VisDataSet::ACC_DATA);

  if(getParameter("modeluvf", false)->data_.hasValue())
    writeDataToFile(getStringVal("modeluvf"), VisDataSet::ACC_MODEL);
  
  if(getParameter("resuvf", false)->data_.hasValue())
    writeDataToFile(getStringVal("resuvf"), VisDataSet::ACC_RES);
}

/**.......................................................................
 * Return the requested image (data, model or residuals, displayed as
 * beam-corrected dirty image, or snr, depending on imgType)
 */
Image VisDataSet::getImage(AccumulatorType accType, ImageType imgType) 
{
  //------------------------------------------------------------
  // Dirty images (uncorrected for primary beam) don't require special
  // treatment
  //------------------------------------------------------------

  if(imgType == IMG_DIRTY) {

    //------------------------------------------------------------
    // For 'clean' images, we convolve the model with an approximation
    // to the synthesized beam, and add the residuals
    //------------------------------------------------------------

    if(accType == ACC_CLEAN) {

      Image modelImage = getImage(ACC_CLEAN);
      Image resImage   = getImage(ACC_RES);

      return modelImage + resImage;

      //------------------------------------------------------------
      // Otherwise, we are just calculating dirty images
      //------------------------------------------------------------

    } else {
      return getImage(accType);
    }
  }

  //------------------------------------------------------------
  // Otherwise we are requesting beam-corrected images, which must be
  // handled differently
  //------------------------------------------------------------

  Image image, wtimage;
  getImage(image, wtimage, accType);

  switch (imgType) {
  case IMG_SKY:
    {
      image /= wtimage;
      image.setInvalidPixelsTo(0.0);
      image.setUnits(Unit::UNITS_JYBEAM);
    }
    break;
  case IMG_SNR:
    {
      Image wtsqrt = wtimage.getSqrt();
      if(getParameter("wtmin", false)->data_.hasValue()) {
	wtsqrt.invalidateDataLessThan(getDoubleVal("wtmin"));
	wtsqrt.setInvalidPixelsTo(0.0);
      }
      
      image /= wtsqrt;
      image.setInvalidPixelsTo(0.0);
      image.setUnits(Unit::UNITS_SNR);
    }
    break;
  default:
    ThrowError("Unrecognized image type: " << imgType);
    break;
  }

  return image;
}

/**.......................................................................
 * Return the requested dirty mosaicked image
 */
void VisDataSet::getImage(Image& image, Image& wtimage, AccumulatorType type) 
{
  //------------------------------------------------------------
  // If accumulating model or residuals, we need to transform image
  // models first
  //------------------------------------------------------------

  if(type != ACC_DATA && type != ACC_BEAM)
    transformModels();

  //------------------------------------------------------------
  // Now we need to accumulate the requested data
  //------------------------------------------------------------
  
  accumulate(type);
  
  //------------------------------------------------------------
  // Now correct for the beams
  //------------------------------------------------------------
  
  accumulateBeamCorrectedImage(image, wtimage);
}

/**.......................................................................
 * Return the requested average image
 */
Image VisDataSet::getImage(AccumulatorType type) 
{
  //------------------------------------------------------------
  // If accumulating model or residuals, we need to transform image
  // models first
  //------------------------------------------------------------

  if(type != ACC_DATA && type != ACC_BEAM) {
    transformModels();
  }

  //------------------------------------------------------------
  // Now we need to accumulate the requested data
  //------------------------------------------------------------
  
  accumulate(type, utilityGridder_);

  //------------------------------------------------------------
  // Now return the image
  //------------------------------------------------------------

  return utilityGridder_.getImage(false);
}

/**.......................................................................
 * For now, we treat only point sources in the Fourier plane
 */
bool VisDataSet::VisFreqData::isImagePlaneModel(gcp::util::Generic2DAngularModel& model)
{
  return !(model.dataSetType_ & DataSetType::DATASET_PTSRC);
}

/**.......................................................................
 * Inherited method called when checkPosition() is invoked.
 * checkPosition() checks if any position was specified via an
 * external parameter first, so if we get here and
 * hasAbsolutePosition_ is true, then it means a parameter has been
 * set via an initialization script.  In that case, we ignore any
 * position that was specified in the data file.
 */
void VisDataSet::initializePositionDependentData()
{
  //------------------------------------------------------------
  // Set the position
  //------------------------------------------------------------
  
  if(!ra_.hasValue_ && obs_.obsRa_.hasValue_)
    setRa(obs_.obsRa_);

  if(!dec_.hasValue_ && obs_.obsDec_.hasValue_)
    setDec(obs_.obsDec_);

  if(hasAbsolutePosition_) {

    //------------------------------------------------------------
    // And now that the position has been set, set it in all of our
    // model components as well.
    //------------------------------------------------------------
    
    for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
      VisBaselineGroup& groupData = baselineGroups_[iGroup];
      
      for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
	VisStokesData& stokesData = groupData.stokesData_[iStokes];
	
	for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	  VisFreqData& freqData = stokesData.freqData_[iFreq];

	  freqData.imageModelComponent_.setRaDecFft(ra_, dec_);
	  freqData.fourierModelComponent_.setRaDec(ra_, dec_);
	  freqData.griddedData_.setRaDec(ra_, dec_);
	}
      }
    }

    //------------------------------------------------------------
    // Also in the global gridder
    //------------------------------------------------------------

    utilityGridder_.setRaDec(ra_, dec_);
  }
}

/**.......................................................................
 * Return an estimate of the largest primary beam size in this dataset
 */
Angle VisDataSet::estimateLargestPrimaryBeamFwhm()
{
  bool first=true;
  Wavelength wave;
  Frequency freq;
  Length diameter;
  Angle fwhm;

  //------------------------------------------------------------
  // Iterate over distinct baseline groups
  //------------------------------------------------------------

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    Length diam1 = groupData.antennaPair_.first.getDiameter();
    Length diam2 = groupData.antennaPair_.second.getDiameter();

    //------------------------------------------------------------
    // For this group, store the smaller of the diameters
    //------------------------------------------------------------

    diameter = diam1 < diam2 ? diam1 : diam2;

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      freq = stokesData.freqData_[0].frequency_;
      wave.setFrequency(freq);

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];
	
	//------------------------------------------------------------
	// And store the lowest frequency for this group
	//------------------------------------------------------------

	if(freqData.frequency_ < freq)
	  wave.setFrequency(freqData.frequency_);
      }
    }

    //------------------------------------------------------------
    // Now convert to radians, and store the largest size over all
    // groups
    //------------------------------------------------------------

    double rad = 1.22 * wave.meters() / diameter.meters();

    if(first) {
      fwhm.setRadians(rad);
      first = false;
    } else {
      if(rad > fwhm.radians()) {
	fwhm.setRadians(rad);
      }
    }
    
  }

  return fwhm;
}

void VisDataSet::loadDataFromObs()
{
  unsigned npix;
  npix = getUintVal("npix");
  
  Angle size;
  size.setName(name_, "size");
  size.setVal(getDoubleVal("size"), getParameter("size", true)->units_);
  
  image_.xAxis().setNpix(npix);
  image_.xAxis().setAngularSize(size);
  
  image_.yAxis().setNpix(npix);
  image_.yAxis().setAngularSize(size);

  storeDataInternallyOnReadin(true);
  usePerc_ = false;

  // Now load the data as if reading in

  loadDataSingle("");
}

/**.......................................................................
 * Overloaded interface from DataSet.  We ignore the passed noise
 * parameter, because the type and magnitude of the noise to generate
 * is specified in the ObsInfo object that has previously been loaded.
 *
 * Here we simply call observe() to generate simulated visibilities
 * from the installed models
 */
void VisDataSet::simulateData(double sigma)
{
  observe(getObs());
}

void VisDataSet::writeCompositeModelToFile(std::string fileName, double sigma)
{
  ThrowError("No writeCompositeModelToFile() method defined by this inheritor");
}

void VisDataSet::initializeIncludedAntennaTypes(std::string includedTypes)
{
  String antStr(includedTypes);

  bool cont=true;
  while(cont) {
    String typeStr;
    if(antStr.remainder().contains(",")) {
      typeStr = antStr.findNextInstanceOf(",", false, ",", true, true);
    } else {
      typeStr = antStr.remainder();
      cont = false;
    }

    Antenna::AntennaType type = Antenna::getType(typeStr.str());
    includedAntTypes_[type] = type;
  }
}

void VisDataSet::initializeExcludedAntennaNumbers(std::string excludedNos)
{
  String antStr(excludedNos);

  std::vector<double> vals = antStr.parseRange();

  for(unsigned iVal=0; iVal < vals.size(); iVal++) {
    unsigned antNo = (unsigned)(vals[iVal]);
    excludedAntNos_[antNo] = antNo;
  }
}

void VisDataSet::initializeExcludedIfNumbers(std::string excludedIfs)
{
  RangeParser parser;

  bool cont=true;

  String excStr(excludedIfs);

  while(cont) {

    String str;

    if(excStr.remainder().contains(",")) {
      str = excStr.findNextInstanceOf(",", false, ",", true, true);
    } else {
      str = excStr.remainder();
      cont = false;
    }

    std::ostringstream os;
    String ifStr(str);

    if(!ifStr.contains("*")) {
      os << "[" << ifStr.str() << "]";
      ifStr = os.str();
    }

    std::vector<unsigned> inds = parser.extractIndexRange(ifStr);
  
    for(unsigned i=0; i < inds.size(); i++) {
      excludedIfNos_[inds[i]] = inds[i];
    }
  }
}

void VisDataSet::initializeIncludedIfNumbers(std::string includedIfs)
{
  RangeParser parser;

  bool cont=true;

  String incStr(includedIfs);

  while(cont) {

    String str;

    if(incStr.remainder().contains(",")) {
      str = incStr.findNextInstanceOf(",", false, ",", true, true);
    } else {
      str = incStr.remainder();
      cont = false;
    }

    std::ostringstream os;
    String ifStr(str);

    if(!ifStr.contains("*")) {
      os << "[" << ifStr.str() << "]";
      ifStr = os.str();
    }

    std::vector<unsigned> inds = parser.extractIndexRange(ifStr);
  
    for(unsigned i=0; i < inds.size(); i++) {
      includedIfNos_[inds[i]] = inds[i];
    }
  }
}

void VisDataSet::initializeTaper(std::string taperSpec)
{
  String taperStr(taperSpec);
  String valStr, valRadStr;
  double val, rad;
  
  valStr    = taperStr.findNextInstanceOf(",", false, ",", true, true);
  valRadStr = taperStr.remainder();

  if(valStr.isEmpty() || valRadStr.isEmpty()) {
    ThrowError("Invalid uvtaper string: '" << taperSpec << "'.  Should be 'val, uvrad'");
  }

  taper_ = true;

  val = valStr.toFloat();
  rad = valRadStr.toFloat();

  if(val < 0.0) {
    taperInvert_ = true;
    val = -val;
  }

  taperSigma_ = sqrt((rad*rad) / (-2.0 * log(val)));
}

/**.......................................................................
 * See if this antenna is included
 */
bool VisDataSet::antennaIsIncluded(Antenna& ant)
{
  return antennaTypeIsIncluded(ant) && antennaNumberIsIncluded(ant);
}

/**.......................................................................
 * See if this antenna is included
 */
bool VisDataSet::antennaTypeIsIncluded(Antenna& ant)
{
  bool inc = true;
  bool exc = false;

  if(includedAntTypes_.size() == 0)
    inc = true;
  else {
    inc = includedAntTypes_.find(ant.type_) != includedAntTypes_.end();
  }

  return inc;
}

bool VisDataSet::antennaNumberIsIncluded(Antenna& ant)
{
  bool exc = false;

  if(excludedAntNos_.size() == 0)
    exc = false;
  else {
    exc = excludedAntNos_.find(ant.antNo_) != excludedAntNos_.end();
  }

  return !exc;
}

bool VisDataSet::useIfNumber(unsigned ifNo)
{
  return ifNumberIsIncluded(ifNo) && !ifNumberIsExcluded(ifNo);
}

bool VisDataSet::ifNumberIsExcluded(unsigned ifNo)
{
  bool exc = false;

  if(excludedIfNos_.size() == 0)
    exc = false;
  else {
    exc = excludedIfNos_.find(ifNo) != excludedIfNos_.end();
  }

  return exc;
}

bool VisDataSet::ifNumberIsIncluded(unsigned ifNo)
{
  bool inc = false;

  if(includedIfNos_.size() == 0)
    inc = true;
  else {
    inc = includedIfNos_.find(ifNo) != includedIfNos_.end();
  }

  return inc;
}

/**.......................................................................
 * See if this baseline is included
 */
bool VisDataSet::isIncluded(VisBaselineGroup& group)
{
  Antenna& ant1 = group.antennaPair_.first;
  Antenna& ant2 = group.antennaPair_.second;

  return antennaIsIncluded(ant1) && antennaIsIncluded(ant2);
}

/**.......................................................................
 * Create an output listing of the estimated rms for each frequency of
 * each baseline grouping
 */
void VisDataSet::printReImRms()
{
  unsigned nTotal = 0;

  double estimatedWtScale = 0;
  unsigned nEstimate= 0;

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	if(freqData.hasData()) {
	  COUT("Chisq for group " << groupData << " freq = " << freqData.frequency_ 
	       << " (ifno = " << std::setw(2) << std::setfill(' ') << std::right << freqData.ifNo_ << ") = " 
	       << std::fixed << std::setprecision(3) << freqData.estChisq_ 
	       << " (wtscale = " << 1.0/freqData.estChisq_ << ")"
	       << "   Estimated map rms = " 
	       << 1e3*sqrt(1.0/freqData.wtSumTotal_) << " mJy/bm from " << std::setw(8) << std::right << freqData.nVis_ << " visibilities" << " # populated = " << freqData.griddedData_.populatedIndices_.size());
	}
      }
    }
  }
}

/**.......................................................................
 * Compute the estimated weight scaling for this data set
 */
double VisDataSet::estimateWtScale()
{
  double estimatedWtScale = 0;
  unsigned nEstimate= 0;

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	if(freqData.hasData())
	  estimatedWtScale += (freqData.estChisq_ - estimatedWtScale) / (++nEstimate);
      }
    }
  }

  //------------------------------------------------------------
  // Estimate the scale factor by which the noise is wrong.  We are
  // not trying to compensate for slightly inaccurate noise estimates
  // here, but for gross errors in the weight scaling (because, for
  // example, we disagree on the convention for UVF weights -- re/im
  // weights, vs. weights on the total power).  So assume that the
  // scaling can be off by whole numbers (in either direction)
  // ------------------------------------------------------------

  if(debug_) {
    COUT("Estimate wt scale to be (1): " << estimatedWtScale);
  }

  if(estimatedWtScale > 1.0)
    estimatedWtScale = rint(estimatedWtScale);
  else {
    estimatedWtScale = 1.0/rint(1.0/estimatedWtScale);
  }

  if(debug_) {
    COUT("Estimate wt scale to be (2): " << estimatedWtScale);
  }

  return 1.0/estimatedWtScale;
}

/**.......................................................................
 * Print a listing of occupied indices
 */
void VisDataSet::printOccupiedIndices()
{
  unsigned nTotal = 0;
  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];
	nTotal += freqData.griddedData_.populatedIndices_.size();

	COUT("Gridding for group " << groupData << " freq = " << freqData.frequency_ 
	     << " = " << freqData.griddedData_.xAxis().getNpix()
	     << " x " << freqData.griddedData_.yAxis().getNpix());
      }
    }
  }

  COUT("Total occupied size = " << nTotal);
}

/**.......................................................................
 * Calculate the taper to apply to this weight
 */
double VisDataSet::taper(double u, double v)
{
  if(taper_) {
    double taperVal = exp(- (u*u + v*v) / (2*taperSigma_*taperSigma_));
    return taperInvert_ ? 1.0 - taperVal : taperVal;
  } else {
    return 1.0;
  }
}

std::ostream& gcp::datasets::operator<<(std::ostream& os, const VisDataSet::VisBaselineGroup& group)
{
  Antenna& ant1 = (Antenna&)group.antennaPair_.first;
  Antenna& ant2 = (Antenna&)group.antennaPair_.second;

  os << std::setw(4) << std::right << std::setfill(' ') << ant1.typeStr() 
     << " - "
     << std::setw(4) << std::right << std::setfill(' ') << ant2.typeStr();

  return os;
}

void VisDataSet::debugPrint()
{
#if 0
  COUTCOLOR("Amt1 = " << amt1               << "s", "yellow");
  COUTCOLOR("Amt2 = " << amt2               << "s", "yellow");
  COUTCOLOR("Amt3 = " << amt3               << "s", "yellow");
#endif
#ifdef TIMER_TEST
  COUTCOLOR("cct1 = " << cct1               << "s", "yellow");
  COUTCOLOR("cct2 = " << cct2               << "s", "yellow");
  COUTCOLOR("cct3 = " << cct3               << "s", "yellow");
#endif
}

/**.......................................................................
 * Return stats about which IFs are currently included by user
 * selection
 */
void VisDataSet::getIfOs(std::ostream& os)
{
  os << std::endl << std::endl << "of which you have included " 
     << (includedIfNos_.size() > 0 ? includedIfNos_.size() : obs_.frequencies_.size())
     << " (";

  if(includedIfNos_.size() > 0) {
    unsigned iIf=0;
    for(std::map<unsigned, unsigned >::iterator iter = includedIfNos_.begin(); 
	iter != includedIfNos_.end(); iter++) {
      os << iter->second;
      if(iIf < includedIfNos_.size()-1)
	os << ", ";
      ++iIf;
    }
  } else {
    for(unsigned iIf=0; iIf < obs_.frequencies_.size(); iIf++) {
      os << (iIf+1);
      if(iIf < obs_.frequencies_.size()-1)
	os << ", ";
    }
  }

  os << ")";
}

/**.......................................................................
 * Return stats about which antennas are currently included by user
 * selection
 */
void VisDataSet::getAntOs(std::ostream& os, std::vector<Antenna>& uniqueAnts)
{
  if(uniqueAnts.size() > 1)
    os << "There are " << uniqueAnts.size() << " distinct antenna types (";
  else
    os << "There is " << uniqueAnts.size() << " distinct antenna type (";

  for(unsigned iAnt=0; iAnt < uniqueAnts.size(); iAnt++) {
    os << uniqueAnts[iAnt].type_;
    if(iAnt < uniqueAnts.size()-1)
      os << ", ";
  }
  os << "), of which you have included " 
     << (includedAntTypes_.size() > 0 ? includedAntTypes_.size() : uniqueAnts.size())
     << " (";

  if(includedAntTypes_.size() > 0) {
    unsigned iAnt=0;
    for(std::map<gcp::util::Antenna::AntennaType, gcp::util::Antenna::AntennaType>::iterator iter = includedAntTypes_.begin(); 
	iter != includedAntTypes_.end(); iter++) {
      os << iter->second;
      if(iAnt < includedAntTypes_.size()-1)
	os << ", ";
      ++iAnt;
    }
  } else {
    for(unsigned iAnt=0; iAnt < uniqueAnts.size(); iAnt++) {
      os << uniqueAnts[iAnt].type_;
      if(iAnt < uniqueAnts.size()-1)
	os << ", ";
    }
  }

  os << ")" << std::endl;

  os << "There are " << obs_.antennas_.size() << " antennas";
  os << ", of which you have excluded " << excludedAntNos_.size();

  if(excludedAntNos_.size() > 0) {
    os << " (";
    unsigned iAnt=0;
    for(std::map<unsigned, unsigned>::iterator iter = excludedAntNos_.begin(); 
	iter != excludedAntNos_.end(); iter++) {
      os << iter->second;
      if(iAnt < excludedAntNos_.size()-1)
	os << ", ";
      ++iAnt;
    }
    os << ")";
  }
}

/**.......................................................................
 * Return a string describing this baseline group
 */
void VisDataSet::getBaseOs(std::ostream& os)
{
  if(baselineGroups_.size() > 1)
    os << "There are now " << baselineGroups_.size() << " distinct baseline groupings (";
  else
    os << "There is now 1 distinct baseline grouping (";

  for(unsigned iBase=0; iBase < baselineGroups_.size(); iBase++) {
    VisBaselineGroup& group = baselineGroups_[iBase];
    Antenna& ant1 = group.antennaPair_.first;
    Antenna& ant2 = group.antennaPair_.second;

    os << ant1.type_ << "-" << ant2.type_;
    if(iBase < baselineGroups_.size()-1) {
      os << ", ";
    }
  }

  os << ")";
}

/**.......................................................................
 * Method by which this dataset is displayed
 */
void VisDataSet::displayIfRequested()
{
  double vpmin = 0.1;
  double vpmax = 0.9;
  double vsep  = 0.15;
  double dvp   = (vpmax-vpmin-vsep)/2;
  double vp11  = vpmin;
  double vp12  = vpmin+dvp;
  double vp21  = vpmin+dvp+vsep;
  double vp22  = vpmax;

  PgUtil::initialize();

  if(getParameter("display", false)->data_.hasValue()) {
    if(getBoolVal("display")) {
      
      if(getParameter("dev", false)->data_.hasValue())
	PgUtil::open(getStringVal("dev"));
      else
	PgUtil::open("/xs");

      PgUtil::setVp(false);
      PgUtil::setOverplot(true);

      cpgsvp(vp11, vp12, vp21, vp22);
      display();

      cpgsvp(vp21, vp22, vp21, vp22);
      radPlot();

      cpgsvp(vp11, vp12, vp11, vp12);
      uvPlot();

      cpgsvp(vp21, vp22, vp11, vp12);
      displayBeam();

      PgUtil::close();

      // Reset any legacy display settings

      PgUtil::setZmin(0.0);
      PgUtil::setZmax(0.0);
    }
  }
}

/**.......................................................................
 * Plot the uv sampling for this dataset
 */
void VisDataSet::uvPlot()
{
  std::vector<double> u;
  std::vector<double> v;

  unsigned nVis=0;

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];
	if(freqData.hasData()) {
	  for(unsigned i=0; i < freqData.griddedData_.populatedIndices_.size(); i++) {
	    unsigned ind = freqData.griddedData_.populatedIndices_[i];
	    double uVal, vVal, val;
	    freqData.griddedData_.getUVData(ind, Dft2d::DATA_UV, uVal, vVal, val);
	    u.push_back( uVal/1000);
	    v.push_back( vVal/1000);
	    u.push_back(-uVal/1000);
	    v.push_back(-vVal/1000);
	  }
	}
      }
    }
  }

  PgUtil::setWnad(true);
  PgUtil::linePlot(u, v, "U (k\\gl)", "V (k\\gl)", "", false);
}

/**.......................................................................
 * Compute a radial plot of amplitudes for this dataset
 */
void VisDataSet::radPlot()
{
  std::vector<double> r;
  std::vector<double> amp;
  std::vector<double> phase;

  unsigned nVis=0;

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];
	if(freqData.hasData()) {
	  for(unsigned i=0; i < freqData.griddedData_.populatedIndices_.size(); i++) {
	    unsigned ind = freqData.griddedData_.populatedIndices_[i];
	    double uVal, vVal, reVal, imVal;
	    freqData.griddedData_.getUVData(ind, Dft2d::DATA_REAL, uVal, vVal, reVal);
	    freqData.griddedData_.getUVData(ind, Dft2d::DATA_IMAG, uVal, vVal, imVal);

	    uVal /= 1000;
	    vVal /= 1000;

	    r.push_back(sqrt((uVal*uVal) + (vVal*vVal)));
	    amp.push_back(sqrt((reVal*reVal + imVal*imVal)));
	  }
	}
      }
    }
  }

  PgUtil::setWnad(false);
  PgUtil::setReverseX(false);
  PgUtil::linePlot(r, amp, "UV Radius (k\\gl)", "Amp (Jy)", "", false);
}

/**.......................................................................
 * Construct a header string for displays of this dataset
 */
std::string VisDataSet::displayHeaderString(VisDataSet::AccumulatorType type)
{
  std::ostringstream os;
  String source = obs_.getSourceName();
  source.strip(' ');
  source.strip('\r');
  source.strip('\n');
  source.strip('"');
  source.strip('\0');

  os << typeString(type) << " (" << source << ")";

  return os.str();
}

std::string VisDataSet::typeString(VisDataSet::AccumulatorType type)
{
  switch (type) {
  case ACC_BEAM:
    return "Synthesized Beam";
    break;
  case ACC_RES:
    return "Residuals";
    break;
  case ACC_MODEL:
    return "Model";
    break;
  case ACC_CLEAN:
    return "Clean image";
    break;
  case ACC_CLEANMODEL:
    return "Clean model";
    break;
  case ACC_CLEANBEAM:
    return "Clean beam";
    break;
  default:
    return "Data";
    break;
  }
}

/**.......................................................................
 * Inherited method to setup for markov display.  This will be called
 * before models are asserted, so that internals can be set up
 * accordingly.
 *
 * For this dataset, currently the only thing that has to happen is
 * that if clean images were requested, we need to initialize our
 * model component DFTs' populated indices to all, since we will
 * iterate over all Fourier components, not just those sampled by the
 * data
 */
void VisDataSet::initializeForDisplay()
{
  bool clean = false;

  if(getParameter("clean", false)->data_.hasValue())
    clean = getBoolVal("clean");

  if(clean) {

    //------------------------------------------------------------
    // Now that the data have been read in, convert second moments to
    // error in mean
    //------------------------------------------------------------
    
    for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
      VisBaselineGroup& groupData = baselineGroups_[iGroup];
      
      for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
	VisStokesData& stokesData = groupData.stokesData_[iStokes];
	
	for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	  VisFreqData& freqData = stokesData.freqData_[iFreq];
	  
	  freqData.fourierModelComponent_.initializePopulatedIndicesToAll();
	  freqData.compositeFourierModelDft_.initializePopulatedIndicesToAll();
	}
      }
    }
  }
}

/**.......................................................................
 * Perform an iterative image-plane CLEAN of the image
 */
void VisDataSet::clean()
{
  //------------------------------------------------------------
  // First, clear the model to get rid of any legacy models that might
  // have been installed before us
  //------------------------------------------------------------

  clearModel();

  //------------------------------------------------------------
  // Iterate, removing delta functions
  //------------------------------------------------------------

  double val, flux, fluxTotal = 0.0;
  Angle xOff, yOff;
  unsigned iMax;
  Image image;
  std::vector<Model*> deltaComponents;

  image = getImage(ACC_RES);

  //------------------------------------------------------------
  // Parse window specifications now.  Pass in the image so that if
  // the user specified 'peak' or '-peak', we can replace with the
  // offsets corresponding to those
  //------------------------------------------------------------

  std::vector<gcp::util::Image::Window> windows = VisDataSet::parseWindows(getStringVal("cleanwindow"), image);
  
  double cleanGain   = getDoubleVal("cleangain");
  double cleanCutoff = getDoubleVal("cleancutoff");
  unsigned cleanIter = getUintVal("cleaniter");

  if(cleanGain > 1.0)
    ThrowSimpleColorError("CLEAN gain should be <= 1.0", "red");

#if 1
  COUT("Cleaning with niter = " << cleanIter << " gain = " << cleanGain << " and windows: ");
  for(unsigned iWin=0; iWin < windows.size(); iWin++) 
    COUT(windows[iWin] << std::endl);
#endif

  //------------------------------------------------------------
  // Iterate, removing delta functions
  //------------------------------------------------------------

  unsigned ndigit = ceil(log10((double)cleanIter));

  for(unsigned i=0; i < cleanIter; i++) {

    COUTCOLORNNL("\rClean iteration: " << std::setw(ndigit) << std::right << (i+1) << "/" 
		 << cleanIter << " (" << std::right << setprecision(0) << std::fixed << (100*(double)(i+1)/cleanIter) << "%)", "green");

    //------------------------------------------------------------
    // Get the image and display prior to removal of the delta component
    //------------------------------------------------------------

    image = getImage(ACC_RES);
    image.getMax(val, xOff, yOff, iMax, windows, Image::TYPE_ABS);

    if(fabs(image.data_[iMax]) < cleanCutoff)
      break;

    flux = cleanGain*image.data_[iMax];
    Model* model = addDeltaFunctionModel(deltaComponents, flux, xOff, yOff);
    fluxTotal += flux;

    //------------------------------------------------------------
    // Now remove the model from the dataset
    //------------------------------------------------------------
    
    addModel(*model);
  }

  //------------------------------------------------------------
  // Done with clean -- displaying restored image
  //------------------------------------------------------------

  COUTCOLOR("\rAdded models with combined flux of: " << fluxTotal, "green");
  display(VisDataSet::ACC_CLEAN);

  for(unsigned iModel=0; iModel < deltaComponents.size(); iModel++) {
    delete deltaComponents[iModel];
    deltaComponents[iModel] = 0;
  }
}

/**.......................................................................
 * Install a PgModel to display the synthesized beam
 */
void VisDataSet::insertSynthesizedBeamModelForPlots(PgModelManager& pgManager)
{
  Image image = utilityGridder_.getImage();

  double xImWidth = image.xAxis().getAngularSize().degrees();
  double yImWidth = image.yAxis().getAngularSize().degrees();

  gcp::models::Generic2DGaussian model;

  model.getVar("xoff")->setVal(xImWidth/2 - (xImWidth)*0.1, "deg");
  model.getVar("xoff")->setUnits("deg");
  model.getVar("xoff")->wasSpecified_ = true;
  
  model.getVar("yoff")->setVal(-yImWidth/2 + (yImWidth)*0.1, "deg");
  model.getVar("yoff")->setUnits("deg");
  model.getVar("yoff")->wasSpecified_ = true;
  
  model.getVar("rotang")->setVal(synthBeamRotAngle_.degrees(), "deg");
  model.getVar("rotang")->setUnits("deg");
  model.getVar("rotang")->wasSpecified_ = true;
  
  model.getVar("majSigma")->setVal(synthBeamMajSig_.degrees(), "deg");
  model.getVar("majSigma")->setUnits("deg");
  model.getVar("majSigma")->wasSpecified_ = true;

  model.getVar("axialRatio")->setVal(synthBeamMinSig_.degrees()/synthBeamMajSig_.degrees(), "");
  model.getVar("axialRatio")->setUnits("");
  model.getVar("axialRatio")->wasSpecified_ = true;

  PgModel pgModel = model.pgModel();
  pgModel.drawCenter_ = false;
  pgModel.fill_       = true;
  pgModel.type_       = PgModel::TYPE_BEAM;

  pgManager.models_.push_back(pgModel);
  pgManager.display_ = true;
  PgUtil::insertPgManager(pgManager);
}

/**.......................................................................
 * Install a PgModel to display the synthesized beam
 */
void VisDataSet::storeSynthesizedBeamModelForPlots(PgModelManager& pgManager)
{
  Image image = utilityGridder_.getImage();

  double xImWidth = image.xAxis().getAngularSize().degrees();
  double yImWidth = image.yAxis().getAngularSize().degrees();

  gcp::models::Generic2DGaussian model;

  model.getVar("xoff")->setVal(xImWidth/2 - (xImWidth)*0.1, "deg");
  model.getVar("xoff")->setUnits("deg");
  model.getVar("xoff")->wasSpecified_ = true;
  
  model.getVar("yoff")->setVal(-yImWidth/2 + (yImWidth)*0.1, "deg");
  model.getVar("yoff")->setUnits("deg");
  model.getVar("yoff")->wasSpecified_ = true;
  
  model.getVar("rotang")->setVal(synthBeamRotAngle_.degrees(), "deg");
  model.getVar("rotang")->setUnits("deg");
  model.getVar("rotang")->wasSpecified_ = true;
  
  model.getVar("majSigma")->setVal(synthBeamMajSig_.degrees(), "deg");
  model.getVar("majSigma")->setUnits("deg");
  model.getVar("majSigma")->wasSpecified_ = true;

  model.getVar("axialRatio")->setVal(synthBeamMinSig_.degrees()/synthBeamMajSig_.degrees(), "");
  model.getVar("axialRatio")->setUnits("");
  model.getVar("axialRatio")->wasSpecified_ = true;

  PgModel pgModel = model.pgModel();
  pgModel.drawCenter_ = false;
  pgModel.fill_       = true;
  pgModel.type_       = PgModel::TYPE_BEAM;

  pgManager.models_.push_back(pgModel);
  pgManager.display_ = true;
}

Model* VisDataSet::addDeltaFunctionModel(std::vector<Model*>& models, double fluxJy, Angle& xOff, Angle& yOff)
{
  Model* model = new gcp::models::PtSrcModel();

  model->getVar("Sradio")->setVal(fluxJy, "Jy");
  model->getVar("Sradio")->wasSpecified_ = true;
  model->getVar("Sradio")->setUnits("Jy");

  model->getVar("xoff")->setVal(xOff.degrees(), "deg");
  model->getVar("xoff")->setUnits("deg");
  model->getVar("xoff")->wasSpecified_ = true;

  model->getVar("yoff")->setVal(yOff.degrees(), "deg");
  model->getVar("yoff")->setUnits("deg");
  model->getVar("yoff")->wasSpecified_ = true;

  model->getVar("spectralType")->setVal("alpha");
  model->getVar("spectralType")->wasSpecified_ = true;

  model->getVar("spectralIndex")->setVal(0.0, "");
  model->getVar("spectralIndex")->wasSpecified_ = true;

  models.push_back(model);

  return model;
}

/**.......................................................................
 * Parse window specifications of the format [xmin:xmax, ymin:ymax, units]
 */
std::vector<gcp::util::Image::Window> 
VisDataSet::parseWindows(std::string windowSpec, Image& image)
{
  String nextWindow, xRange, yRange, units;
  std::vector<gcp::util::Image::Window> windows;

  String windowSpecStr(windowSpec);

  do {
    nextWindow = windowSpecStr.findNextInstanceOf("[", true, "]", true, false);
    nextWindow.strip(' ');

    if(!nextWindow.isEmpty()) {

      xRange = nextWindow.findNextInstanceOf(" ", false, ",", true,  false);
      yRange = nextWindow.findNextInstanceOf(",", true,  ",", true,  false);
      units  = nextWindow.findNextInstanceOf(",", true,  " ", false, false);

      Image::Window window = image.getWindow(xRange, yRange, units);
      windows.push_back(window);

      addDisplayWindow(window);
    }

  } while(!nextWindow.isEmpty());

  return windows;
}

void VisDataSet::addDisplayWindow(Image::Window& win)
{
  PgModel pgModel;
  pgModel.drawCenter_ = false;
  pgModel.fill_       = false;
  pgModel.type_       = PgModel::TYPE_BOX;

  pgModel.xMin_        = win.xMin_.degrees();
  pgModel.xMax_        = win.xMax_.degrees();
  pgModel.yMin_        = win.yMin_.degrees();
  pgModel.yMax_        = win.yMax_.degrees();

  PgUtil::addModel(pgModel);
}

void VisDataSet::setWedge(VisDataSet::AccumulatorType type, double zmin, double zmax, bool specified)
{
  //------------------------------------------------------------
  // Always display a wedge if displaying the data or the beam
  //------------------------------------------------------------

  switch (type) {
  case VisDataSet::ACC_BEAM:
  case VisDataSet::ACC_DATA:
    PgUtil::setWedge(true);
    break;

    // Else display a wedge only if zmin/zmax were explicitly specified and equal to each other

  default:
    PgUtil::setWedge(specified && (zmin == zmax));
    break;
  }
}

/**.......................................................................
 * Parse a filename specification, which we allow to be of the form
 *
 *   'file = name'
 *
 * or
 * 
 *   'file = name, shift = dx, dy units'
 */
std::string VisDataSet::parseFileName(std::string fileName, bool& shiftRequested, Angle& xoff, Angle& yoff)
{
  String fileStr(fileName);
  fileStr.strip(" ");

  std::string file = fileStr.findNextInstanceOf(" ", false, ",", false, true).str();
  String rem = fileStr.remainder();
  
  if(rem.isEmpty()) {
    shiftRequested = false;
    return file;
  }
    
  if(!rem.contains("shift="))
    ThrowSimpleColorError("Unable to parse file string: '" << fileStr << "'", "red");
  
  String coord = rem.findNextInstanceOf("shift=", true, " ", false, true);
  
  double xVal = coord.findNextInstanceOf(" ", false, ",", true, true).toDouble();
  double yVal = coord.findNextNumericString().toDouble();
  
  String units = coord.remainder();
  
  xoff.setVal(xVal, units.str());
  yoff.setVal(yVal, units.str());

  shiftRequested = true;

  return file;
}

/**.......................................................................
 * Co-add two VisDataSet objects
 */
void VisDataSet::operator+=(VisDataSet& vds)
{
  //------------------------------------------------------------
  // Iterate over all VisFreqData objects in the passed datasets,
  // adding them to ours.  We don't assume that these objects
  // necessarily match; if a match is found, we coadd, if no match is
  // found, we will add new objects as needed
  //------------------------------------------------------------

  for(unsigned iGroup=0; iGroup < vds.baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& group = vds.baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < group.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = group.stokesData_[iStokes];

      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	if(freqData.hasData())
	  mergeData(vds, group, stokesData, freqData);
      }
    }
  }

  //------------------------------------------------------------
  // Re-compute the combined synthesized beam solid angle
  //------------------------------------------------------------

  computeGlobalSynthesizedBeam();

  //------------------------------------------------------------
  // And get the estimated synthesized beam parameters for all
  // VisFreqData subsets
  //------------------------------------------------------------

  estimateSynthesizedBeams();
}

/**.......................................................................
 * Merge a frequency group into this object
 */
void VisDataSet::mergeData(VisDataSet& vds, VisBaselineGroup& findGroup, VisStokesData& findStokes, VisFreqData& findFreq)
{
  bool groupMatch  = false;
  bool stokesMatch = false;
  bool freqMatch   = false;

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& group = baselineGroups_[iGroup];

    if(group == findGroup) {
      groupMatch = true;

      for(unsigned iStokes=0; iStokes < group.stokesData_.size(); iStokes++) {
	VisStokesData& stokesData = group.stokesData_[iStokes];
	
	if(stokesData == findStokes) {
	  stokesMatch = true;

	  for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	    VisFreqData& freqData = stokesData.freqData_[iFreq];

	    if(freqData == findFreq) {
	      freqMatch = true;

	      //------------------------------------------------------------
	      // Call VisFreqData increment operator
	      //------------------------------------------------------------
	      
	      freqData.mergeData(findFreq, vds.xShift_, vds.yShift_);
	    }
	  }
	}
      }
    }
  }

  if(!(groupMatch && stokesMatch && freqMatch)) {
    std::ostringstream os;
    os << "Dumping data for: " << findGroup << " " << findStokes << " " << findFreq << "(";
    if(!groupMatch) {
      os << "Found no matching group)";
    } else {
      os << "Found matching group";
      if(!stokesMatch) {
	os << ", but no matching Stokes parameter)";
      } else {
	os << ", and matching Stokes parameter, but no matching frequency)";
      }
    }

    COUT(os.str());
  }
}

/**.......................................................................
 * Set a parameter
 */
void VisDataSet::setParameter(std::string name, std::string val, std::string units)
{
  String nameStr(name);

  //------------------------------------------------------------
  // If this is the file keyword, add a dataset corresponding to it
  //------------------------------------------------------------

  if(name == "file") {
    fileList_.push_back(val);
    addDataSet(val);
  }

  //------------------------------------------------------------
  // If this is a parameter for a dataset, set its parameter instead
  //------------------------------------------------------------

  if(nameStr.contains("dataset")) {
    String dataSetName = nameStr.findNextInstanceOf("", false, ".", true, true);
    String parName     = nameStr.remainder();

    DataSet* dataSet = getDataSet(dataSetName.str());
    dataSet->setParameter(parName.str(), val, units);

    //------------------------------------------------------------
    // Else set our parameter
    //------------------------------------------------------------

  } else {
    gcp::util::ParameterManager::setParameter(name, val, units);
  }
}

/**.......................................................................
 * Increment a parameter
 */
void VisDataSet::incrementParameter(std::string name, std::string val, std::string units)
{
  String nameStr(name);

  //------------------------------------------------------------
  // If incrementing the file list, add a new dataset
  //------------------------------------------------------------

  if(name == "file") {
    fileList_.push_back(val);
    addDataSet(val);
  }

  //------------------------------------------------------------
  // If this is a parameter for a dataset, increment its parameter instead
  //------------------------------------------------------------

  if(nameStr.contains("dataset")) {
    String dataSetName = nameStr.findNextInstanceOf("", false, ".", true, true);
    String parName     = nameStr.remainder();

    DataSet* dataSet = getDataSet(dataSetName.str());
    dataSet->incrementParameter(parName.str(), val, units);

    //------------------------------------------------------------
    // Else increment our parameter
    //------------------------------------------------------------

  } else {
    gcp::util::ParameterManager::incrementParameter(name, val, units);
  }
}

/**.......................................................................
 * Method used for loading data into this object internally
 */
ACC_DISPATCH_FN(VisDataSet::dispatchInternal)
{
  if(src.useIfNumber(freq.ifNo_))
    freq.accumulateMoments(first, data, src.xShift_, src.yShift_);
}

//=======================================================================
// Methods used for loading data into this object from an external
// VisDataSet
//=======================================================================

ACC_DISPATCH_FN(VisDataSet::dispatchExternal)
{
  if(src.useIfNumber(freq.ifNo_)) {

    //------------------------------------------------------------
    // Find the VisFreqData object that this freq object maps to, and
    // accumulate the data into that object instead
    //------------------------------------------------------------

    VisFreqData* destFreq = dest.freqMap_[&freq];
    destFreq->accumulateMoments(first, data, src.xShift_, src.yShift_);

    //------------------------------------------------------------
    // But accumulate variance stats into the parent object, in case
    // we want to check weight scalings on data load
    //------------------------------------------------------------

    freq.accumulateVarianceStats(data);

    //------------------------------------------------------------
    // If requested to store the data, put it in the right place
    //------------------------------------------------------------

    destFreq->group_->dataset_->storeDataIfRequested(data, destFreq, src, vis);
  }
}

/**.......................................................................
 * Store data if requested to store on read-in
 */
void VisDataSet::storeDataIfRequested(VisData& data, VisFreqData* destFreq, VisDataSet& src, ObsInfo::Vis& vis)
{
  if(!storeDataInternally_)
    return;

  unsigned iGroup, iVis, iDate;
  getStoreIndices(data, destFreq, src, iGroup, iVis, iDate);

  //  COUT("iGroup = " << iGroup);
  ObsInfo::Vis& visOut = obs_.visibilities_[iGroup];

  visOut.re_[iVis] = data.re_;
  visOut.im_[iVis] = data.im_;
  visOut.wt_[iVis] = data.wt_;

  visOut.u_  = vis.u_;
  visOut.v_  = vis.v_;
  visOut.w_  = vis.w_;

  //------------------------------------------------------------
  // If we are concatentatin multiple datasets, fake the timestamps by
  // adding iDate sec to the start JD.  (This is because programs such
  // as difmap will choke if dates are out of order)
  //------------------------------------------------------------

  if(fileList_.size() > 1)
    visOut.jd_ = writeJd_ + iDate * 1.0/(24*3600);
  else
    visOut.jd_ = vis.jd_;

  visOut.baseline_ = getAipsBaselineIndex(src, data.baseline_);

  //  COUT("SToring data: iGroup = " << iGroup << " iVis = " << iVis << " iDate = " << iDate << " baseline = " << data.baseline_ << " u = " << data.u_ << " v = " << data.v_ << " re = " << data.re_ << " im = " << data.im_ << " wt = " << data.wt_ << " my base = " << visOut.baseline_ << " jd = " << setprecision(12) << visOut.jd_ << " orig jd = " << data.jd_);
}

unsigned VisDataSet::getAipsBaselineIndex(VisDataSet& dataset, unsigned aipsIndex)
{
  // Get the native antenna numbering of this baseline

  unsigned antNo1 = aipsIndex / 256;
  unsigned antNo2 = aipsIndex % 256;

  // Now find the antenna numbers that correspond to these in the
  // current antenna array

  std::map<unsigned, unsigned>& indMap = datasetAntMap_[&dataset];

  //  COUT("Before antNo1 = " << antNo1 << " antNo2 = " << antNo2);

  antNo1 = indMap[antNo1];
  antNo2 = indMap[antNo2];

  //  COUT("after antNo1 = " << antNo1 << " antNo2 = " << antNo2);

  return antNo1 * 256 + antNo2;
}

/**.......................................................................
 * Return the index into our stored array of this group & vis
 */
void VisDataSet::getStoreIndices(VisData& data, VisFreqData* destFreq, VisDataSet& src, unsigned& iGroup, unsigned& iVis, unsigned& iDate)
{
  //------------------------------------------------------------
  // We increment the date index if the date has changed
  //------------------------------------------------------------

  if(data.jd_ != lastJd_)
    ++lastDate_;

  //------------------------------------------------------------
  // We increment the group index if either the date or the baseline
  // has changed
  //------------------------------------------------------------

  if(data.jd_ != lastJd_ || data.baseline_ != lastBaseline_)
    ++lastGroup_;

  lastJd_       = data.jd_;
  lastBaseline_ = data.baseline_;

  //------------------------------------------------------------
  // Get the global stokes and frequency index that this object
  // belongs to
  //------------------------------------------------------------

  unsigned iStokes = destFreq->stokes_->globalStokesIndex_;
  unsigned iFreq   = destFreq->globalIfIndex_;

  unsigned nStokes = allStokes_.size();
  unsigned nFreq   = allFreqs_.size();

  iGroup = lastGroup_;
  iVis   = iStokes * nStokes + iFreq;
  iDate  = lastDate_;
}

/**.......................................................................
 * Find the VisFreqData object that matches the passed group, stokes
 * and freq
 */
VisDataSet::VisFreqData& VisDataSet::findMatch(VisBaselineGroup& groupMatch, VisStokesData& stokesMatch, VisFreqData& freqMatch)
{
  for(unsigned iBase=0; iBase < baselineGroups_.size(); iBase++) {
    VisBaselineGroup& groupData = baselineGroups_[iBase];

    if(groupData == groupMatch) {
      for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
	VisStokesData& stokesData = groupData.stokesData_[iStokes];
	
	if(stokesData == stokesMatch) {
	  for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	    VisFreqData& freqData = stokesData.freqData_[iFreq];

	    if(freqData == freqMatch) {
	      return freqData;
	    }
	  }
	} // End stokes match

      }
    } // End group match

  } // End iteration over groups

  ThrowError("No match found for group " << groupMatch << " stokes " << stokesMatch << " freq " << freqMatch);
}

/**.......................................................................
 * Build a single unique antenna map from all datasets
 *
 * At the end we want a single array that encompasses the maximum
 * number of antennas of a given type over all datasets
 *
 * We additionally want to store a map of antenna number in each
 * dataset to the corresponding antenna in our ordering scheme
 */
void VisDataSet::buildAntennaMap(std::vector<VisDataSet*>& datasets)
{
  for(unsigned i=0; i < datasets.size(); i++) {
    VisDataSet* dataset = datasets[i];
    for(unsigned iAnt=0; iAnt < dataset->obs_.antennas_.size(); iAnt++) {
      //      COUTCOLOR("ant " << iAnt << " = " << dataset->obs_.antennas_[iAnt].antNo_, "green");
      obs_.antennas_ = obs_.mergeAnts(dataset->obs_);
      updateAntennaMap(dataset);
    }
  }

#if 0
  // Print the map

  for(std::map<VisDataSet*, std::map<unsigned, unsigned> >::iterator iter=datasetAntMap_.begin(); iter != datasetAntMap_.end(); iter++) {
    COUTCOLOR("Dataset = " << iter->first << ": ", "green");
    std::map<unsigned,unsigned>& indMap = iter->second;
    
    for(std::map<unsigned, unsigned>::iterator indIter=indMap.begin(); indIter != indMap.end(); indIter++) {
      COUTCOLOR("antNo " << indIter->first << " maps to " << indIter->second, "green");
    } 
  }
#endif

}

/**.......................................................................
 * Update our antenna map from this VisDataSet
 */
void VisDataSet::updateAntennaMap(VisDataSet* dataset)
{
  //------------------------------------------------------------
  // Iterate through each antenna in the dataset, finding the
  // corresponding entry in our map
  //------------------------------------------------------------

  for(unsigned iAnt=0; iAnt < obs_.antennas_.size(); iAnt++)
    obs_.antennas_[iAnt].antNo_ = iAnt + 1;

  std::map<Antenna*, std::vector<Antenna*> > thisVecMap = obs_.getAntVecMap();
  std::map<Antenna*, std::vector<Antenna*> > thatVecMap = dataset->obs_.getAntVecMap();

  std::map<unsigned, unsigned>& indMap = datasetAntMap_[dataset];

  for(std::map<Antenna*, std::vector<Antenna*> >::iterator iter=thatVecMap.begin(); iter != thatVecMap.end(); iter++) {

    Antenna* key = obs_.getAntMapKey(iter->first, thisVecMap);

    std::vector<Antenna*>& thisVec = thisVecMap[key];
    std::vector<Antenna*>& thatVec = thatVecMap[iter->first];

    for(unsigned iAnt=0; iAnt < thatVec.size(); iAnt++) {
      unsigned thisAntNo = thisVec[iAnt]->antNo_;
      unsigned thatAntNo = thatVec[iAnt]->antNo_;
      indMap[thatAntNo] = thisAntNo;
    }
  }
}

/**.......................................................................
 * Iterate through all datasets, constructing a map of unique baseline
 * groups, stokes parameters and frequencies
 */
void VisDataSet::buildInternalMap(std::vector<VisDataSet*>& datasets)
{
  //------------------------------------------------------------
  // First iterate over all datasets, and construct a map of all
  // antennas
  //------------------------------------------------------------
  
  buildAntennaMap(datasets);
  
  //------------------------------------------------------------
  // Now use that antenna map to determine unique baseline groupings
  //------------------------------------------------------------
  
  determineUniqueBaselineGroupings(obs_);
  
  //------------------------------------------------------------
  // Now determine which stokes parameters and frequencies are present
  // for each baseline group
  //------------------------------------------------------------
  
  initializeBaselineGroups(datasets);

  //------------------------------------------------------------
  // Construct a map from VisFreqData objects in the dataset
  // vector to VisFreqData objects in this dataset
  //------------------------------------------------------------

  mapFrequencies(datasets);

  //------------------------------------------------------------
  // Remap stokes and frequencies across all baseline groups
  // to reflect the total number of stokes and freqs encountered
  //------------------------------------------------------------

  remapStokesAndFreqs();

  //------------------------------------------------------------
  // Install store arrays if requested to store the data on read-in
  //------------------------------------------------------------

  initializeStoreArrays(datasets);
}

/**.......................................................................
 * Iterate through all datasets, constructing a map of unique Stokes
 * parameters and frequencies for each group, then initialize the
 * group to reflect that structure
 */
void VisDataSet::mapFrequencies(std::vector<VisDataSet*>& datasets)
{
  //------------------------------------------------------------
  // Iterate over all VisFreqData objects in all datasets,
  // constructing a map of matching VisFreqData objects in this object
  //------------------------------------------------------------

  for(unsigned i=0; i < datasets.size(); i++) {
    VisDataSet* dataset = datasets[i];
    
    for(unsigned iGroup=0; iGroup < dataset->baselineGroups_.size(); iGroup++) {
      VisBaselineGroup& groupData = dataset->baselineGroups_[iGroup];
      
      for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
	VisStokesData& stokesData = groupData.stokesData_[iStokes];
	
	for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	  VisFreqData& freqData = stokesData.freqData_[iFreq];
	  
	  if(dataset->useIfNumber(freqData.ifNo_))
	    freqMap_[&freqData] = &findMatch(groupData, stokesData, freqData);
	}
      }
    }
  }
}

/**.......................................................................
 * Iterate through all datasets, constructing a map of unique Stokes
 * parameters and frequencies for each group, then initialize the
 * group to reflect that structure
 */
void VisDataSet::initializeStoreArrays(std::vector<VisDataSet*>& datasets)
{
  if(!storeDataInternally_)
    return;

  unsigned nGroup  = 0;
  unsigned nStokes = allStokes_.size();
  unsigned nFreq   = allFreqs_.size();

  //------------------------------------------------------------
  // Iterate over all datasets, storing the maximum number of groups
  // we might have to store
  //------------------------------------------------------------

  for(unsigned i=0; i < datasets.size(); i++) {
    VisDataSet* dataset = datasets[i];
    nGroup += dataset->obs_.nGroup_;
  }

  //------------------------------------------------------------
  // And resize our internal arrays to match
  //------------------------------------------------------------

  obs_.visibilities_.resize(nGroup);
  obs_.setNumberOfAntennas(obs_.antennas_.size());
  obs_.setNumberOfBaselines(1);
  obs_.setNumberOfTimestamps(nGroup);
  obs_.markAntennaLocationsAsReceived();
  obs_.setNumberOfStokesParameters(allStokes_.size());

  std::vector<Frequency> freqs(allFreqs_.size());
  std::vector<Frequency>   bws(allFreqs_.size());

  unsigned iFreq=0;
  for(std::map<double, unsigned>::iterator iter=allFreqs_.begin(); iter != allFreqs_.end(); iter++, iFreq++) {
    freqs[iFreq].setGHz(iter->first);
    bws[iFreq].setGHz(allBws_[iter->first]);
  }

  obs_.setFrequencyInformation(freqs, bws);

  for(unsigned iGroup=0; iGroup < obs_.visibilities_.size(); iGroup++) {

    // Each group will be the maximum size, even if no data exists for
    // it

    ObsInfo::Vis& vis = obs_.visibilities_[iGroup];
    vis.initialize(nStokes * nFreq);
  }
}

/**.......................................................................
 * Iterate through all datasets, constructing a map of unique Stokes
 * parameters and frequencies for each group, then initialize the
 * group to reflect that structure
 */
void VisDataSet::remapStokesAndFreqs()
{
  //------------------------------------------------------------
  // Assign sequential indices in the maps of all Stokes and
  // Frequencies
  //------------------------------------------------------------
  
  unsigned iStokes=0;
  for(std::map<gcp::util::Stokes::Param, unsigned>::iterator iter=allStokes_.begin(); iter != allStokes_.end(); iter++, iStokes++)
    iter->second = iStokes;
  
  unsigned iFreq=0;
  for(std::map<double, unsigned>::iterator iter=allFreqs_.begin(); iter != allFreqs_.end(); iter++, iFreq++)
    iter->second = iFreq;

  //------------------------------------------------------------
  // Now iterate over all VisFreqData objects in this dataset,
  // assigning indices to them from the maps
  //------------------------------------------------------------

  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];

    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      stokesData.globalStokesIndex_ = allStokes_[stokesData.stokes_];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];
	freqData.globalIfIndex_ = allFreqs_[freqData.frequency_.GHz()];
      }
    }
  }
}

/**.......................................................................
 * Iterate through all datasets, constructing a map of unique Stokes
 * parameters and frequencies for each group, then initialize the
 * group to reflect that structure
 */
void VisDataSet::initializeBaselineGroups(std::vector<VisDataSet*>& datasets)
{
  //------------------------------------------------------------
  // Iterate over each unique baseline group in this dataset.  For
  // each one, we search for any baseline group that matches in the
  // list of datasets, and find all Stokes parameters and frequencies
  // that occur for this group
  //------------------------------------------------------------

  for(unsigned iBase=0; iBase < baselineGroups_.size(); iBase++) {

    VisBaselineGroup& group = baselineGroups_[iBase];
    std::map<gcp::util::Stokes::Param, std::map<double, VisFreqData*> > uniqueStokes;

    for(unsigned i=0; i < datasets.size(); i++) {
      VisDataSet* dataset = datasets[i];
      
      for(unsigned iGroup=0; iGroup < dataset->baselineGroups_.size(); iGroup++) {
	VisBaselineGroup& groupData = dataset->baselineGroups_[iGroup];
	
	//------------------------------------------------------------
	// If this matches the current group, proceed
	//------------------------------------------------------------
	
	if(groupData == group) {
	  
	  for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
	    VisStokesData& stokesData = groupData.stokesData_[iStokes];
	    
	    for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	      VisFreqData& freqData = stokesData.freqData_[iFreq];
	      
	      //------------------------------------------------------------
	      // Don't add this frequency in if we were told to exclude
	      // it
	      //------------------------------------------------------------
	      
	      if(!dataset->useIfNumber(freqData.ifNo_))
		continue;
	      
	      std::map<double, VisFreqData*>& freqMap = uniqueStokes[stokesData.stokes_];
	      
	      VisFreqData* freqPtr = freqMap[freqData.frequency_.GHz()];
	      
	      // Store any stokes and frequency we encounter

	      allStokes_[stokesData.stokes_]       = 0;
	      allFreqs_[freqData.frequency_.GHz()] = 0;
	      allBws_[freqData.frequency_.GHz()]   = freqData.bandwidth_.GHz();

	      //------------------------------------------------------------
	      // Store the maximum UV over all frequency objects
	      // matching this frequency
	      //------------------------------------------------------------
	      
	      if(freqPtr != 0) {
		freqData.uAbsMax_ = (freqData.uAbsMax_ > freqPtr->uAbsMax_ ? freqData.uAbsMax_ : freqPtr->uAbsMax_);
		freqData.vAbsMax_ = (freqData.vAbsMax_ > freqPtr->vAbsMax_ ? freqData.vAbsMax_ : freqPtr->vAbsMax_);
	      }
	      
	      freqMap[freqData.frequency_.GHz()] = &freqData;
	    }
	  }
	}
      }
    }

    //------------------------------------------------------------
    // And initialize this group from the stokes map
    //------------------------------------------------------------
    
    group.initialize(uniqueStokes);
    group.dataset_ = this;
  }
}

/**.......................................................................
 * Iterate over all VisFreqData objects in this dataset, storing the
 * contribution to the weight sum from the last dataset added
 */
void VisDataSet::storeWtSums(gcp::util::Angle& xShift, gcp::util::Angle& yShift)
{
  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];
    
    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];
	freqData.storeWtSum(xShift, yShift);
      }
    }
  }
}

void VisDataSet::printReIm()
{
  for(unsigned iGroup=0; iGroup < baselineGroups_.size(); iGroup++) {
    VisBaselineGroup& groupData = baselineGroups_[iGroup];
    
    for(unsigned iStokes=0; iStokes < groupData.stokesData_.size(); iStokes++) {
      VisStokesData& stokesData = groupData.stokesData_[iStokes];
      
      for(unsigned iFreq=0; iFreq < stokesData.freqData_.size(); iFreq++) {
	VisFreqData& freqData = stokesData.freqData_[iFreq];

	for(unsigned iPop=0; iPop < freqData.griddedData_.populatedIndices_.size(); iPop++) {
	  unsigned ind = freqData.griddedData_.populatedIndices_[iPop];
	  COUT("re = " << freqData.griddedData_.out_[ind][0] << " im = " << freqData.griddedData_.out_[ind][1] << " reerr = " 
	       << freqData.griddedData_.errorInMean_[ind][0] << " imerr = " << freqData.griddedData_.errorInMean_[ind][1]);
	}
      }
    }
  }
}

void VisDataSet::writeDataToFile(std::string fileName, VisDataSet::AccumulatorType type)
{
  replaceVisibilities(type);
  writeUvfFile(fileName);
}

DataSet* VisDataSet::getDataSet(std::string name)
{
  for(unsigned i=0; i < datasets_.size(); i++) {
    if(datasets_[i]->name_ == name)
      return (DataSet*)datasets_[i];
  }

  ThrowSimpleColorError("No dataset '" << name << "' has been defined", "red");
  return 0;
}
