
#include "gcp/matlab/UvfWriter.h"

#include "gcp/fftutil/FitsIoHandler.h"

#include "gcp/util/Astrometry.h"
#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/MiriadIo.h"

using namespace std;
using namespace gcp::matlab;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
UvfWriter::UvfWriter(unsigned nTelescope) :
  nTelescope_(nTelescope) {
  //  dfreq_.setGHz(0.5);
  fixedBaselines_=false;
  coordsAreJ2000_=false;
  firstTelescopeNum_=0;
}

/**.......................................................................
 * Destructor.
 */
UvfWriter::~UvfWriter() {}

/**.......................................................................
*  Set parameters appropriate for writing out Miriad data to FITS
*/
void UvfWriter::setCoordsToJ2000(const mxArray* array) {
  mp_.setTo(array);
  
  if(mp_.getNumberOfDimensions() != 2 || 
     mp_.getDimension(0) != 1 || mp_.getDimension(1) != 1) {
    ThrowError("CoordsAreJ2000 must be of size 1 x 1");
  }
  
  bool* coordsAreJ2000 = mp_.getLogicalData();

  coordsAreJ2000_=coordsAreJ2000[0];
}

void UvfWriter::setNumberOfTelescopes(const mxArray* array) {
  mp_.setTo(array);
  
  if(mp_.getNumberOfDimensions() != 2 || 
     mp_.getDimension(0) != 1 || mp_.getDimension(1) != 1) {
    ThrowError("Number of telescopes must be of size 1 x 1");
  }
  
  unsigned* ntel = mp_.getUintData();

  nTelescope_=ntel[0];
}

void UvfWriter::setFirstTelescopeNum(const mxArray* array) {
  mp_.setTo(array);
  
  if(mp_.getNumberOfDimensions() != 2 || 
     mp_.getDimension(0) != 1 || mp_.getDimension(1) != 1) {
    ThrowError("Number of first telescope must be of size 1 x 1");
  }
  
  unsigned* itel = mp_.getUintData();

  firstTelescopeNum_=itel[0];
}

void UvfWriter::setDeltaIfFrequencies(const mxArray* array)
{
  mp_.setTo(array);
  
  if(mp_.getNumberOfDimensions() != 2 || 
     mp_.getDimension(0) != 1 || 
     mp_.getDimension(1) != nFrequency_)
    ThrowError("Delta frequency must be of size 1 x x Nfrequency (=" << nFrequency_ << ")");
  
  double* dfreq = mp_.getDoubleData();
  
  dfreq_.resize(nFrequency_);
  
  for(unsigned iFreq=0; iFreq < nFrequency_; iFreq++)
    dfreq_[iFreq].setHz(dfreq[iFreq]);
}

void UvfWriter::setBaselines(const mxArray* array)
{
  mp_.setTo(array);
  
  if(mp_.getNumberOfDimensions() != 2 || 
     mp_.getDimension(0) != 1)
    ThrowError("Baselines parameter must be of size 1 x nBaseline)");
  
  unsigned* baselines = mp_.getUintData();
  
  baselines_=baselines;
  fixedBaselines_=true;
}

/**.......................................................................
 * Set a file name
 */
void UvfWriter::setFileName(const mxArray* fileName, const mxArray* openMode)
{
  fileName_ =  mxArrayToString(fileName);

  if(openMode != 0)
    openMode_ =  mxArrayToString(openMode);
  else
    openMode_ = "new";
}

/**.......................................................................
 * Set a src name
 */
void UvfWriter::setSourceName(const mxArray* array)
{
  if(!MexParser::isString(array)) {
    ThrowError("Source name array does not represent a string");
  }
  srcName_ =  mxArrayToString(array);
}

/**.......................................................................
 * Set a pointer to the data array
 */
void UvfWriter::setMiriadData(const mxArray* array)
{
  mp_.setTo(array);

  if(mp_.getNumberOfDimensions() != 4 || 
     mp_.getDimension(3) != 2) {
    ThrowError("Visibility data must be of size Nframe x Nbaseline x Nfrequency x 2");
  }

  nFrame_     = mp_.getDimension(0);
  nBaseline_  = mp_.getDimension(1);
  nFrequency_ = mp_.getDimension(2);

  data_ = mp_.getDoubleData();
}

/**.......................................................................
 * Set a pointer to the data array
 */
void UvfWriter::setUvfData(const mxArray* array)
{
  mp_.setTo(array);

  if(mp_.getNumberOfDimensions() != 4 || 
     mp_.getDimension(3) != 3) {
    ThrowError("Visibility data must be of size Nframe x Nbaseline x Nfrequency x 3");
  }

  nFrame_     = mp_.getDimension(0);
  nBaseline_  = mp_.getDimension(1);
  nFrequency_ = mp_.getDimension(2);

  data_ = mp_.getDoubleData();
}

/**.......................................................................
 * Set a pointer to the rms array (not currently used -- rms is the
 * 3rd dimension of the uvf data array)
 */
void UvfWriter::setRms(const mxArray* array)
{
  mp_.setTo(array);
  
  if(mp_.getNumberOfDimensions() != 3 || 
     mp_.getDimension(1) != 8 ||
     mp_.getDimension(2) != 16)
    ThrowError("Rms array must be of size N x 8 x 16");
  
  rms_ = mp_.getDoubleData();
}

/**.......................................................................
 * Set a pointer to the dates
 */
void UvfWriter::setDate(const mxArray* array)
{
  mp_.setTo(array);
  
  if(mp_.getNumberOfDimensions() != 2 || 
     mp_.getDimension(0) != nFrame_ ||
     mp_.getDimension(1) != 2)
    ThrowError("Date array must be of size Nframe (=" << nFrame_ << ") x 2");
  
  date_ = mp_.getDoubleData();
}

/**.......................................................................
 * Set the observing RA/DEC
 */
void UvfWriter::setRefCoord(const mxArray* array)
{
  mp_.setTo(array);
  
  if(mp_.getNumberOfDimensions() != 2 || 
     mp_.getDimension(0) != 1 || 
     mp_.getDimension(1) != 2) {
    ThrowError("Reference coordinate array must be of size 1 x 2");
  }
  
  double* coord = mp_.getDoubleData();

  refRa_.setHours(coord[0]);
  refDec_.setDegrees(coord[1]);
}

/**.......................................................................
 * Set the reference position RA/DEC
 */
void UvfWriter::setCoord(const mxArray* array)
{
  mp_.setTo(array);
  
  if(mp_.getNumberOfDimensions() != 2 || 
     mp_.getDimension(0) != 1 || 
     mp_.getDimension(1) != 3) {
    ThrowError("Reference coordinate array must be of size 1 x 3");
  }
  
  double* coord = mp_.getDoubleData();

  obsRa_.setHours(coord[0]);
  obsDec_.setDegrees(coord[1]);
  obsMjd_.setMjd(coord[2]);
}

/**.......................................................................
 * Set the frequencies
 */
void UvfWriter::setFreq(const mxArray* array)
{
  mp_.setTo(array);
  
  if(mp_.getNumberOfDimensions() != 2 || 
     mp_.getDimension(0) != 1 || 
     mp_.getDimension(1) != nFrequency_)
    ThrowError("Frequency array must be of size 1 x Nfrequency (=" << nFrequency_ << ")");
  
  double* freq = mp_.getDoubleData();
  
  frequencies_.resize(nFrequency_);
  
  for(unsigned iFreq=0; iFreq < nFrequency_; iFreq++)
    frequencies_[iFreq].setHz(freq[iFreq]);
}

/**.......................................................................
 * Set the XYZ coordinates of each antenna
 */
void UvfWriter::setXyz(const mxArray* array)
{
  mp_.setTo(array);
  
  if(mp_.getNumberOfDimensions() != 2 || 
     mp_.getDimension(0) != 3 || 
     mp_.getDimension(1) != nTelescope_)
    ThrowError("Xyz array must be of size 3 x nTelescope");
  
  double* xyz = mp_.getDoubleData();

  xyz_.resize(nTelescope_);
  
  unsigned index=0;
  for(unsigned iTel=0; iTel < nTelescope_; iTel++) {
    xyz_[iTel].X_.setMeters(xyz[0]);
    xyz_[iTel].Y_.setMeters(xyz[1]);
    xyz_[iTel].Z_.setMeters(xyz[2]);
  }
}

/**.......................................................................
 * Set the UVW coordinates of each point
 */
void UvfWriter::setUvw(const mxArray* array)
{
  mp_.setTo(array);
  
  if(mp_.getNumberOfDimensions() != 3 || 
     mp_.getDimension(0) != nFrame_ || 
     mp_.getDimension(1) != nBaseline_ ||
     mp_.getDimension(2) != 3) {
    ThrowError("Uvw array must be of size Nframe (=" << nFrame_ << ") x Nbaseline (=" << nBaseline_ << ") x 3");
  }
  
  uvw_ = mp_.getDoubleData();
}

void UvfWriter::writeUvfFile()
{
  FitsIoHandler fitsio;
  HourAngle meanRa;
  Declination meanDec;

  // Need to convert from epoch of observation to standard epoch
  
  if (!coordsAreJ2000_) {
    gcp::util::Astrometry::apparentToJ2000Place(obsRa_, obsDec_, obsMjd_, 
		  			        meanRa, meanDec);
  } else {
    meanRa=obsRa_;
    meanDec=obsDec_;
  }

  fitsio.setTelescopeName("CARMA");
  fitsio.setInstrumentName("CARMA");

  fitsio.setRa(meanRa);
  fitsio.setDec(meanDec);

  fitsio.setNumberOfChannelsPerIf(1);
  fitsio.setNumberOfStokesParameters(1);

  fitsio.setSourceName(srcName_);
  fitsio.setTelescopeLocations(xyz_);

  fitsio.setIfFrequencies(frequencies_);
  fitsio.setDeltaIfFrequencies(dfreq_);

  fitsio.setDeltaIfFrequency(dfreq_[0]);
  fitsio.setDeltaChannelFrequency(dfreq_[0]);

  fitsio.setNumberOfTelescopes(nTelescope_);
  fitsio.setNumberOfBaselines(nBaseline_);

  if(fixedBaselines_) 
    fitsio.setBaselines(baselines_);

  fitsio.setFirstTelescopeNum(firstTelescopeNum_);
  fitsio.setNumberOfFrames(nFrame_);
  fitsio.setDate();

  fitsio.openFile(fileName_);

  // Write the fits file

  fitsio.installVisibilityData(data_, date_, uvw_);
  fitsio.writeUvfFile();
}

void UvfWriter::writeFakeUvfFile()
{
  FitsIoHandler fitsio;
  HourAngle meanRa;
  Declination meanDec;

  // Need to convert from epoch of observation to standard epoch
  
  if (!coordsAreJ2000_) {
    gcp::util::Astrometry::apparentToJ2000Place(obsRa_, obsDec_, obsMjd_, 
		  			        meanRa, meanDec);
  } else {
    meanRa=obsRa_;
    meanDec=obsDec_;
  }

  fitsio.setRa(meanRa);
  fitsio.setDec(meanDec);

  fitsio.setSourceName(srcName_);
  fitsio.setTelescopeLocations(xyz_);
  fitsio.setIfFrequencies(frequencies_);
  fitsio.setNumberOfBaselines(nBaseline_);
  fitsio.setNumberOfFrames(nFrame_);
  fitsio.setDate();

  fitsio.openFile(fileName_);

  // Write the fits file

  fitsio.writeFakeUvfFile(data_, date_, uvw_);
}

void UvfWriter::writeMiriadFile()
{
  MiriadIo miriadio;
  HourAngle meanRa;
  Declination meanDec;

  // Need to convert from epoch of observation to standard epoch
  
  if (!coordsAreJ2000_) {
    gcp::util::Astrometry::apparentToJ2000Place(obsRa_, obsDec_, obsMjd_, 
		  			        meanRa, meanDec);
  } else {
    meanRa=obsRa_;
    meanDec=obsDec_;
  }

  miriadio.setRa(meanRa);
  miriadio.setDec(meanDec);

  miriadio.setRaApp(obsRa_);
  miriadio.setDecApp(obsDec_);

  miriadio.setRaRef(refRa_);
  miriadio.setDecRef(refDec_);

  miriadio.setSourceName(srcName_);
  miriadio.setTelescopeLocations(xyz_);
  miriadio.setIfFrequencies(frequencies_);
  miriadio.setNumberOfBaselines(nBaseline_); // Was 28
  miriadio.setFirstTelescopeNum(firstTelescopeNum_);
  if (fixedBaselines_) miriadio.setBaselines(baselines_);
  miriadio.setNumberOfFrames(nFrame_);

  // Write the miriad file

  miriadio.openFile(fileName_, openMode_);
  miriadio.writeFile(data_, date_, uvw_, rms_);
  miriadio.closeFile();
}

void UvfWriter::writeFakeMiriadFile()
{
  MiriadIo miriadio;
  HourAngle meanRa;
  Declination meanDec;

  // Need to convert from epoch of observation to standard epoch
  
  if (!coordsAreJ2000_) {
    gcp::util::Astrometry::apparentToJ2000Place(obsRa_, obsDec_, obsMjd_, 
		  			        meanRa, meanDec);
  } else {
    meanRa=obsRa_;
    meanDec=obsDec_;
  }

  miriadio.setRa(meanRa);
  miriadio.setDec(meanDec);

  miriadio.setRa(obsRa_);
  miriadio.setDec(obsDec_);

  miriadio.setSourceName(srcName_);
  miriadio.setTelescopeLocations(xyz_);
  miriadio.setIfFrequencies(frequencies_);
  miriadio.setNumberOfBaselines(28);
  miriadio.setNumberOfFrames(nFrame_);

  // Write the fits file

  miriadio.openFile(fileName_, openMode_);
  miriadio.writeFakeFile(data_, date_, uvw_, rms_);
  miriadio.closeFile();
}

void UvfWriter::writeMiriadMosaic()
{
  MiriadIo miriadio;
  HourAngle meanRa;
  Declination meanDec;

  HourAngle refRa;
  Declination refDec;

  refRa.setHours("0:17:59.999");
  refDec.setDegrees("16:22:29.992");

  miriadio.setRaRef(refRa);
  miriadio.setDecRef(refDec);


  // Need to convert from epoch of observation to standard epoch
  
  if (!coordsAreJ2000_) {
    gcp::util::Astrometry::apparentToJ2000Place(obsRa_, obsDec_, obsMjd_, 
		  			        meanRa, meanDec);
  } else {
    meanRa=obsRa_;
    meanDec=obsDec_;
  }

  miriadio.setTelescopeLocations(xyz_);
  miriadio.setIfFrequencies(frequencies_);
  miriadio.setNumberOfBaselines(28);
  miriadio.setNumberOfFrames(nFrame_);

  // Write the fits file

  gcp::util::Declination ddec;
  gcp::util::HourAngle dra;

  // Center pointing

  miriadio.setRa(meanRa);
  miriadio.setDec(meanDec);

  miriadio.setRaApp(meanRa);
  miriadio.setDecApp(meanDec);

  ddec.setRadians(0.0);
  dra.setRadians(0.0);
  miriadio.setDRaApp(dra);
  miriadio.setDDecApp(ddec);

  miriadio.setSourceName("mos_0");

  miriadio.openFile("mos0", openMode_);
  miriadio.writeFakeFile(data_, date_, uvw_, rms_);
  miriadio.closeFile();

  // Offset right 

  ddec.setArcMinutes(6.0);
  miriadio.setDDecApp(ddec);

  miriadio.openFile("mos1", openMode_);
  miriadio.writeFakeFile(data_, date_, uvw_, rms_);
  miriadio.closeFile();

  // Offset left

  ddec.setArcMinutes(-6.0);
  miriadio.setDDecApp(ddec);

  miriadio.openFile("mos2", openMode_);
  miriadio.writeFakeFile(data_, date_, uvw_, rms_);
  miriadio.closeFile();
}

void UvfWriter::writeMiriadMosaicTest()
{
  MiriadIo miriadio;
  HourAngle meanRa;
  Declination meanDec;

  HourAngle refRa;
  Declination refDec;

  refRa.setHours("0:17:59.999");
  refDec.setDegrees("16:22:29.992");

  miriadio.setRaRef(refRa);
  miriadio.setDecRef(refDec);

  // Need to convert from epoch of observation to standard epoch
  
  if (!coordsAreJ2000_) {
    gcp::util::Astrometry::apparentToJ2000Place(obsRa_, obsDec_, obsMjd_, 
		  			        meanRa, meanDec);
  } else {
    meanRa=obsRa_;
    meanDec=obsDec_;
  }

  miriadio.setTelescopeLocations(xyz_);
  miriadio.setIfFrequencies(frequencies_);
  miriadio.setNumberOfBaselines(28);
  miriadio.setNumberOfFrames(nFrame_);

  // Write the fits file

  gcp::util::Declination ddec;
  gcp::util::HourAngle dra;

  // Center pointing

  meanRa.setHours("0:17:59.999");
  meanDec.setDegrees("16:22:29.992");

  miriadio.setRa(meanRa);
  miriadio.setDec(meanDec);
  miriadio.setSourceName("mos_0");
  miriadio.openFile("mos0", openMode_);
  miriadio.writeFakeFile(data_, date_, uvw_, rms_);
  miriadio.closeFile();

  // Offset right 

  meanRa.setHours("0:17:59.999");
  meanDec.setDegrees("16:28:29.992");

  miriadio.setRa(meanRa);
  miriadio.setDec(meanDec);
  miriadio.openFile("mos1", openMode_);
  miriadio.writeFakeFile(data_, date_, uvw_, rms_);
  miriadio.closeFile();

  // Offset left

  meanRa.setHours("0:17:59.999");
  meanDec.setDegrees("16:14:29.992");

  miriadio.setRa(meanRa);
  miriadio.setDec(meanDec);
  miriadio.openFile("mos2", openMode_);
  miriadio.writeFakeFile(data_, date_, uvw_, rms_);
  miriadio.closeFile();
}
