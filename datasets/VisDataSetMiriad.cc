#include "gcp/datasets/VisDataSetMiriad.h"

#include "gcp/util/FitsBinTableReader.h"
#include "gcp/util/MiriadIo.h"

using namespace gcp::datasets;
using namespace gcp::util;
using namespace std;

/**.......................................................................
 * Constructor.
 */
VisDataSetMiriad::VisDataSetMiriad(ThreadPool* pool) :
  VisDataSet(pool) {}

/**.......................................................................
 * Destructor.
 */
VisDataSetMiriad::~VisDataSetMiriad() {}

/**.......................................................................
 * Initialize the file reader
 */
void VisDataSetMiriad::openFileReader(std::string fileName)
{
  reader_.openFileForRead(fileName);
}

/**.......................................................................
 * Close the file readear
 */
void VisDataSetMiriad::closeFileReader()
{
  reader_.closeFile();
}

/**.......................................................................
 * Method to get the next group of visibilties from a file
 */
void VisDataSetMiriad::getGroup(unsigned iGroup, ObsInfo::Vis& vis)
{
}

/**.......................................................................
 * Initialize antenna/baseline information from the antenna table (if any) 
 */
void VisDataSetMiriad::initializeAntennaInformation(std::string fileName)
{
  // Initialize antenna information from the FITS binary table

  gcp::util::MiriadIo reader;
  reader.openFileForRead(fileName);

  // Only resize the antenna array if it hasn't already been set
  // up manually
  
  if(obs_.antennas_.size() == 0) {
    obs_.setNumberOfAntennas(reader.getNumberOfTelescopes());

    // If the antennas have already been specified, check that
    // that specification agrees with the file

  } else if(obs_.antennas_.size() != reader.getNumberOfTelescopes()) {
    ReportError("Warning: The previously specified number of antennas (" << obs_.antennas_.size() << ")"
		<< " doesn't match what's in the Miriad file (" << reader.getNumberOfTelescopes() << ")");
  }

  //------------------------------------------------------------
  // See if the antenna locations were specified
  //------------------------------------------------------------
  
  std::vector<LengthTriplet> xyz = reader.getTelescopeLocations();
  for(unsigned iAnt=0; iAnt < obs_.antennas_.size(); iAnt++) {
    obs_.setAntennaLocation(xyz[iAnt], iAnt);
  }
}

/**.......................................................................
 * Update frequency information from the UVF file
 */
void VisDataSetMiriad::updateFrequencyInformation()
{
}

/**.......................................................................
 * Update RA/DEC information
 */
void VisDataSetMiriad::updateObservationInformation()
{
  obs_.setSourceName(reader_.getSourceName());
  obs_.setObsRa(reader_.getRaApp());
  obs_.setObsDec(reader_.getDecApp());
  obs_.setObsEquinox(2000);
  obs_.setTelescopeName(reader_.getTelescopeName());
  obs_.setInstrumentName("Unknown");
}

/**.......................................................................
 * Update visibility information
 */
void VisDataSetMiriad::updateVisibilityInformation()
{
}

void VisDataSetMiriad::initializeFileReader(std::string fileName)
{
}
