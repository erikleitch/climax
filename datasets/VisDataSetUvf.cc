#include "gcp/datasets/VisDataSetUvf.h"

#include "gcp/util/FitsBinTableReader.h"

#include <vector>

using namespace gcp::datasets;
using namespace gcp::util;
using namespace std;

/**.......................................................................
 * Constructors.
 */
VisDataSetUvf::VisDataSetUvf(ThreadPool* pool) :
  VisDataSet(pool) {}

/**.......................................................................
 * Destructor.
 */
VisDataSetUvf::~VisDataSetUvf() {}

/**.......................................................................
 * Initialize the file reader
 */
void VisDataSetUvf::openFileReader(std::string fileName)
{
  reader_.openFile(fileName);
}

/**.......................................................................
 * Close the file readear
 */
void VisDataSetUvf::closeFileReader()
{
  reader_.close();
}

/**.......................................................................
 * Method to get the next group of visibilities from a file
 */
void VisDataSetUvf::getGroup(unsigned iGroup, ObsInfo::Vis& vis)
{
  FitsUvfReader::Vis uvfVis;
  reader_.readData(iGroup, uvfVis);
  vis = uvfVis;
}

/**.......................................................................
 * Initialize antenna/baseline information from the antenna table (if any) 
 */
void VisDataSetUvf::initializeAntennaInformation(std::string fileName)
{
  //------------------------------------------------------------
  // Initialize antenna information from the FITS binary table
  //------------------------------------------------------------

  FitsBinTableReader tableReader(fileName);

  //------------------------------------------------------------
  // Iterate over all tables in the file until we find an antenna
  // table, or we reach the end
  //------------------------------------------------------------

  bool stop=false;
  while(!stop) {

    try {

      tableReader.getNextTableInfo();

      //------------------------------------------------------------
      // If this is an AIPS AN table, the number of antennas should be
      // the number of rows in the table
      //------------------------------------------------------------

      if(tableReader.tableName() == "AIPS AN") {

	//------------------------------------------------------------
	// Only resize the antenna array if it hasn't already been set
	// up manually
	//------------------------------------------------------------

        if(obs_.antennas_.size() == 0) {

	  std::ostringstream os;
	  os << tableReader.nRow();
	  obs_.setParameter("nant", os.str());

	  //------------------------------------------------------------
	  // Initialize the antenna numbers
	  //------------------------------------------------------------

	  for(unsigned iAnt=0; iAnt < obs_.antennas_.size(); iAnt++) {
	    obs_.antennas_[iAnt].antNo_ = iAnt+1;
	  }

	  //------------------------------------------------------------
	  // If the antennas have already been specified, check that
	  // that specification agrees with the AIPS AN table
	  //------------------------------------------------------------

	} else if(obs_.antennas_.size() != tableReader.nRow()) {
	  ReportError("Warning: The previously specified number of antennas (" << obs_.antennas_.size() << ")"
		      << " doesn't match what's in the AIPS AN table (" << tableReader.nRow() << ")");
	}

	//------------------------------------------------------------
	// See if the antenna locations were specified
	//------------------------------------------------------------

	unsigned nCol = tableReader.nCol();

	for(unsigned iCol = 0; iCol < nCol; iCol++) {

	  if(tableReader.colName(iCol) == "LX" || tableReader.colName(iCol) == "LY" 
	     || tableReader.colName(iCol) == "LZ") {

	    std::vector<double> vals = tableReader.getDoubleData(iCol);

	    if(vals.size() != obs_.antennas_.size()) {
	      ThrowError("Column size is > number of antennas");
	    }

	    for(unsigned iAnt=0; iAnt < obs_.antennas_.size(); iAnt++) {

	      Length length(Length::Meters(), vals[iAnt]);

	      if(tableReader.colName(iCol)        == "LX") {
		obs_.setAntennaX(length, iAnt);
	      } else if(tableReader.colName(iCol) == "LY") {
		obs_.setAntennaY(length, iAnt);
	      } else if(tableReader.colName(iCol) == "LZ") {
		obs_.setAntennaZ(length, iAnt);
	      }

	    }

	  }


	  if(tableReader.colName(iCol) == "ANT NO.") {

	    std::vector<short> vals = tableReader.getShortData(iCol);

	    if(vals.size() != obs_.antennas_.size()) {
	      ThrowError("Column size is > number of antennas");
	    }

	    for(unsigned iAnt=0; iAnt < obs_.antennas_.size(); iAnt++) {
	      obs_.antennas_[iAnt].antNo_ = vals[iAnt];
	    }

	  }

	}

	stop = true;
      }
    } catch(Exception& err) {
      stop = true;
    } catch(...) {
      stop = true;
    }
  }

  //------------------------------------------------------------
  // Now that antennas have been initialized, setup an array for
  // mapping from consecutive groups to aips baseline indices, in case
  // these are not specified with the groups
  //------------------------------------------------------------

  obs_.initializeBaselineIndices();
}

/**.......................................................................
 * Initialize frequency information from the FQ table (if any) 
 */
void VisDataSetUvf::initializeFrequencyInformation(std::string fileName)
{
  //------------------------------------------------------------
  // Initialize frequency information from the FITS binary table
  //------------------------------------------------------------

  FitsBinTableReader tableReader(fileName);

  //------------------------------------------------------------
  // Iterate over all tables in the file until we find a frequency
  // table, or we reach the end.
  //------------------------------------------------------------

  bool stop=false;
  while(!stop) {
    try {
      tableReader.getNextTableInfo();

      //------------------------------------------------------------
      // If this is an AIPS FQ table, the number of antennas should be
      // given by the NO_IF keyword
      //------------------------------------------------------------

      if(tableReader.tableName() == "AIPS FQ") {

	unsigned nIf = tableReader.getLongKey("NO_IF");

	std::vector<Frequency> freqs(nIf);
	std::vector<Frequency> bws(nIf);

	//------------------------------------------------------------
	// Read the frequency info
	//------------------------------------------------------------

	unsigned nCol = tableReader.nCol();

	ifOffsetHz_ = tableReader.getDoubleData("IF FREQ");

	bwHz_               = tableReader.getDoubleData("TOTAL BANDWIDTH");
	std::string bwUnit  = tableReader.getUnit("TOTAL BANDWIDTH"); 

	//------------------------------------------------------------
	// Make sure these are actually in Hz
	//------------------------------------------------------------

	Frequency bw;
	for(unsigned i=0; i < bwHz_.size(); i++) {
	  bw.setVal(bwHz_[i], bwUnit);
	  bwHz_[i] = bw.Hz();
	}

	stop = true;
      }
    } catch(...) {
      stop = true;
    }
  }
}

/**.......................................................................
 * Update frequency information from the UVF file
 */
void VisDataSetUvf::updateFrequencyInformation()
{
  double freqsHz0 = reader_.freqs()[0];
  std::vector<double> freqsHz = reader_.ifs();
  double bwHz = reader_.ifDelta();

  std::vector<Frequency> freqs;
  std::vector<Frequency> bws;

  freqs.resize(freqsHz.size());
  bws.resize(freqsHz.size());

  //------------------------------------------------------------
  // If an FQ table was encountered, use the frequency offsets
  // relative to the first IF frequency in the header
  //------------------------------------------------------------

  if(ifOffsetHz_.size() == freqsHz.size()) {
    for(unsigned i=0; i < freqs.size(); i++) {
      freqs[i].setHz(freqsHz0 + ifOffsetHz_[i]);
      bws[i].setHz(bwHz_[i]);
    }

    //------------------------------------------------------------
    // Else use the frequency information from the axis headers
    //------------------------------------------------------------

  } else {
    for(unsigned i=0; i < freqs.size(); i++) {
      freqs[i].setHz(freqsHz[i]);
      bws[i].setHz(bwHz);
    }
  }

  obs_.setFrequencyInformation(freqs, bws);
}

/**.......................................................................
 * Update information about the observation from a dataset
 */
void VisDataSetUvf::updateObservationInformation(gcp::datasets::VisDataSet* dataset)
{
  VisDataSetUvf* uvf = (VisDataSetUvf*)dataset;
  updateObservationInformation(uvf->reader_);
}

/**.......................................................................
 * Update information about the observation from a FITS reader
 */
void VisDataSetUvf::updateObservationInformation(gcp::util::FitsUvfReader& reader)
{
  obs_.setSourceName(    reader.object());
  obs_.setObsRa(         reader.obsra());
  obs_.setObsDec(        reader.obsdec());
  obs_.setObsEquinox(    reader.equinox());
  obs_.setTelescopeName( reader.telescope());
  obs_.setInstrumentName(reader.instrument());
}

/**.......................................................................
 * Update information about the observation
 */
void VisDataSetUvf::updateObservationInformation()
{
  updateObservationInformation(reader_);
}

/**.......................................................................
 * Update visibility information
 */
void VisDataSetUvf::updateVisibilityInformation()
{
  obs_.setNumberOfGroups(reader_.nGroup());
  obs_.setNumberOfStokesParameters(reader_.nStokes());
}

/**.......................................................................
 * For UVF datasets, we just call the underlying writeUvfFile() method
 * to write our composite model to a file
 */
void VisDataSetUvf::writeCompositeModelToFile(std::string fileName, double sigma)
{
  writeUvfFile(fileName);
}

/**.......................................................................
 * Add a dataset to the vector of datasets we will stack into this one
 */
void VisDataSetUvf::addDataSet(std::string file)
{
  VisDataSet* dataset = new VisDataSetUvf();
  datasets_.push_back(dataset);

  std::ostringstream name;
  unsigned iDataSet = datasets_.size();
  name << "dataset" << iDataSet;

  dataset->setName(name.str());

  dataSetMap_[name.str()] = dataset;
  COUT(this << " Calling addDataSet with dataset = " << dataset << " map size = " << dataSetMap_.size() << " file = " << file);

  //------------------------------------------------------------
  // Add a 'parameter' corresponding to this dataset, so it appears
  // in the parameter listing
  //------------------------------------------------------------
  
  std::ostringstream os;
  os << "The dataset corresponding to file " << file;
  addParameter(name.str(),           DataType::STRING, os.str());
}

