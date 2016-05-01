/**.......................................................................
 * MATLAB Mex file for reading from UVF files
 *
 * Use like:
 *
 * d=gcpMatReadUvf({'array.frame.record','corr.band0.usb[0][0]',
 *                     'antenna*.tracker.actual double', 'antenna*.tracker.source string'},
 *                     '06-jan-2005:15','06-jan-2005:16',
 *                     '/data/gcpdaq/arc','/home/gcpdaq/carma_unstable/gcp/array/conf/cal');
 *
 */
#include "gcp/util/Geoid.h"
#include "gcp/util/MiriadIo.h"
#include "gcp/util/FitsUvfReader.h"
#include "gcp/util/FitsBinTableReader.h"

#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/MexParser.h"

#include "mex.h"
#include "matrix.h"

#include <iostream>
#include <math.h>

using namespace std;
using namespace gcp::util;
using namespace gcp::matlab;

#define ASSIGN_VALS(cast)				\
  {							\
    cast* mptr = (cast*)data;				\
    unsigned iM, iC;					\
    unsigned nEl0 = reader.nRow();			\
    unsigned nEl1 = reader.colRepeat(iCol);		\
    for(unsigned iEl1=0; iEl1 < nEl1; iEl1++)		\
      for(unsigned iEl0=0; iEl0 < nEl0; iEl0++) {	\
	iC = iEl0*nEl1 + iEl1;				\
	iM = iEl1*nEl0 + iEl0;				\
	*(mptr + iM) = vals[iC];			\
      }							\
  }

class AntennaInfo {
public:
  std::vector<Length> east_;
  std::vector<Length> up_;
  std::vector<Length> north_;

  std::vector<Length> elAxisHeight_;
  std::vector<Length> diameter_;
  std::vector<Length> sweptVolumeDiameter_;

  AntennaInfo();
  virtual ~AntennaInfo();

  bool isSzaAntShadowed(unsigned iSzaAnt,
			Angle& az, Angle& el, bool doPrint=false);

  bool isAntShadowed(Angle& az, Angle& el,
		     unsigned iAntTarget, unsigned iAntTest, bool doPrint=false);

  static const unsigned nAnt_;
  static const unsigned nOvro_;
  static const unsigned nBima_;
  static const unsigned nSza_;
};


mxArray* readMiriadFile(const mxArray* prhs[]);
void addAxes(mxArray* hdr, FitsUvfReader& reader);
void addTables(mxArray* hdr, std::string fileName);
void addTable(mxArray* tables, FitsBinTableReader& tableReader, std::string extname, unsigned index);
void assignTableData(FitsBinTableReader& reader, unsigned iCol, void* data);
void addData(mxArray* data, std::string fName, unsigned nGroup, unsigned maxVisPerRecord, mxArray* hdr);
void addHeaderData(MiriadIo& io, mxArray* hdr, unsigned nAnt, AntennaInfo& antennaInfo);


/**.......................................................................
 * Entry point from the matlab environment
 */
void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  gcp::util::Logger::installStdoutPrintFn(MexHandler::stdoutPrintFn);
  gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);
  
  gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
  gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
  gcp::util::ErrHandler::installLogFn(MexHandler::logFn);
  
  plhs[0] = readMiriadFile(prhs);
}
   
/**.......................................................................
 * Read the Miriad file
 */
mxArray* readMiriadFile(const mxArray* prhs[])
{
  mxArray* ptr=0;

  // Open the file

  char* cptr = mxArrayToString(prhs[0]);
  std::string fName(cptr);

  gcp::util::MiriadIo io;

  // Get stats we need to read the file

  unsigned nGroup, minVisPerRecord, maxVisPerRecord;
  io.getStats(fName, nGroup, minVisPerRecord, maxVisPerRecord);

  COUT("Found nGroup = " << nGroup);

  // Create a top-level struct
  
  int dims[2]={1,1};
  ptr = mxCreateStructArray(2, dims, 0, NULL);

  // Create two substructures

  mxArray* hdr  = MexHandler::addNamedStructField(ptr, "header");

  //------------------------------------------------------------
  // Add header fields
  //------------------------------------------------------------

  MexHandler::addNamedStringStructField(hdr,       "file", fName);

  unsigned int* ngroup = MexHandler::addNamedUintStructField(hdr,       "ngroup", 1);
  *ngroup = nGroup;

  // Finally, read the data!

  addData(ptr, fName, nGroup, maxVisPerRecord, hdr);

  return ptr;
}

void addData(mxArray* data, std::string fName, unsigned nGroup, unsigned maxVisPerRecord, mxArray* hdr)
{
  FitsUvfReader::Vis vis;
  std::vector<unsigned> dims;
  bool first=true;

  dims.resize(1);

  // Arrays for managing information about sources                                                                                                                                                             
  std::map<std::string, std::vector<unsigned>* > sourceIndMap;
  std::map<std::string, std::string>             sourcePurposeMap;
  std::map<std::string, double>                  sourceObsRaMap;
  std::map<std::string, double>                  sourceObsDecMap;
  std::map<std::string, double>                  sourceRaMap;
  std::map<std::string, double>                  sourceDecMap;
  std::map<std::string, double>                  sourceEpochMap;

  MiriadIo io;
  io.openFile(fName, "old");

  AntennaInfo antennaInfo;

  // Add the groups data field

  mxArray*  groups = MexHandler::addNamedStructField(data, "groups", nGroup);

  // Now iterate over all groups in the file, reading each one into
  // the corresponding element of the matlab array

  bool haveBfMask  = true;
  bool haveAxisOff = true;

  for(unsigned iGroup=0; iGroup < nGroup; iGroup++) {

    // Add singleton data fields for each of the following

    float*       u  = MexHandler::addNamedFloatStructField(groups,      "u", 1, iGroup);
    float*       v  = MexHandler::addNamedFloatStructField(groups,      "v", 1, iGroup);
    unsigned* ant1  = MexHandler::addNamedUintStructField(groups,    "ant1", 1, iGroup);
    unsigned* ant2  = MexHandler::addNamedUintStructField(groups,    "ant2", 1, iGroup);
    double*    mjd  = MexHandler::addNamedDoubleStructField(groups,   "mjd", 1, iGroup);
    double*    lst  = MexHandler::addNamedDoubleStructField(groups,   "lst", 1, iGroup);
    float*  inttime = MexHandler::addNamedFloatStructField(groups,"inttime", 1, iGroup);

    unsigned* ant1Shadowed = MexHandler::addNamedUintStructField(groups,    "ant1Shadowed", 1, iGroup);
    unsigned* ant2Shadowed = MexHandler::addNamedUintStructField(groups,    "ant2Shadowed", 1, iGroup);

    // Now read the next group of visibility data

    io.readNextRecord();

    //------------------------------------------------------------
    // Read the number of antennas
    //------------------------------------------------------------

    unsigned nAnt   = io.readIntVar("nants");

    // If this is the first group, write header fields

    if(iGroup == 0) {
      addHeaderData(io, hdr, nAnt, antennaInfo);
    }

    *inttime = io.readFloatVar("inttime");

    //------------------------------------------------------------
    // Get the source
    //------------------------------------------------------------

    std::string source  = io.readStringVar("source");
    std::string purpose = io.readStringVar("purpose");
    double obsra        = io.readDoubleVar("obsra");
    double obsdec       = io.readDoubleVar("obsdec");
    double ra           = io.readDoubleVar("ra");
    double dec          = io.readDoubleVar("dec");
    float epoch         = io.readFloatVar("epoch");

    // See if we have already encountered this source, or it is a new
    // source

    if(sourceIndMap.find(source) == sourceIndMap.end()) {
      sourceIndMap[source] = new std::vector<unsigned>;
    }

    sourceIndMap[source]->push_back(iGroup+1);
    sourcePurposeMap[source] = purpose;
    sourceObsRaMap[source]   = obsra;
    sourceObsDecMap[source]  = obsdec;
    sourceRaMap[source]      = ra;
    sourceDecMap[source]     = dec;
    sourceEpochMap[source]   = epoch;

    //------------------------------------------------------------
    // Get the baseline index & decompose into antennas
    //------------------------------------------------------------

    unsigned base = (unsigned)io.preamble_[3];

    // base = 256*(ant1) + (ant2)

    *ant1 = base / 256;
    *ant2 = base % 256;

    unsigned iAnt1 = *ant1-1;
    unsigned iAnt2 = *ant2-1;

    //------------------------------------------------------------
    // Add an array of the antenna az & els
    //------------------------------------------------------------

    dims.resize(1);
    dims[0] = nAnt;

    mxArray*   azArr = MexHandler::addNamedStructField(groups,   "az", dims.size(), (const int*)&dims[0], 
						       DataType::DOUBLE, iGroup);

    mxArray*   elArr = MexHandler::addNamedStructField(groups,   "el", dims.size(), (const int*)&dims[0], 
						       DataType::DOUBLE, iGroup);

    mxArray*   axisDelayArr = MexHandler::addNamedStructField(groups,   "axisDelay", dims.size(), (const int*)&dims[0], 
						       DataType::FLOAT, iGroup);

    double* azPtr     = (double*)mxGetData(azArr);
    double* elPtr     = (double*)mxGetData(elArr);
    float*  axisDelayPtr = (float*)mxGetData(axisDelayArr);

    std::vector<double> az = io.readDoubleVar("antaz", nAnt);
    std::vector<double> el = io.readDoubleVar("antel", nAnt);

    for(unsigned i=0; i < nAnt; i++) {
      azPtr[i]        = az[i];
      elPtr[i]        = el[i];
    }

    std::vector<float>  axisDelay;

    if(haveAxisOff) {
      try {
	axisDelay = io.readFloatVar("axisoff", nAnt);
      } catch(...) {
	haveAxisOff = false;
      }
    }

    for(unsigned i=0; i < nAnt; i++) {
      if(axisDelay.size() > 0) {
	axisDelayPtr[i] = axisDelay[i];
      } else {
	axisDelayPtr[i] = 0.0;
      }
    }

    // Calculate if either of the antennas of this baseline were
    // shadowed.

    Angle azAnt, elAnt;
    azAnt.setDegrees(azPtr[iAnt1]);
    elAnt.setDegrees(elPtr[iAnt1]);

    *ant1Shadowed = (unsigned)antennaInfo.isSzaAntShadowed(iAnt1, azAnt, elAnt);

    azAnt.setDegrees(azPtr[iAnt2]);
    elAnt.setDegrees(elPtr[iAnt2]);

    *ant2Shadowed = (unsigned)antennaInfo.isSzaAntShadowed(iAnt2, azAnt, elAnt);

    //------------------------------------------------------------
    // Read other pertinent variables
    //------------------------------------------------------------

    unsigned nSpect = io.readIntVar("nspect");

    //------------------------------------------------------------
    // Add the bfmask array -- if it is present.  We set 01 Mar 2012
    // as the first date to check for this, since it was not really
    // present before
    //------------------------------------------------------------

    double mjdTest = io.preamble_[2] - 2400000.5;
    
    if(haveBfMask) {
      if(mjdTest > 55987) {
	
	try {
	  std::vector<int> bfmaskArr(nSpect);
	  bfmaskArr = io.readIntVar("bfmask", nSpect);
	  
	  unsigned* bfmask = MexHandler::addNamedUintStructField(groups,  "bfmask", nSpect, iGroup);
	  
	  for(unsigned i=0; i < nSpect; i++) 
	    *(bfmask + i) = bfmaskArr[i];
	} catch(...) {
	  haveBfMask = false;
	}
      }
    }

    //------------------------------------------------------------
    // Add arrays that record sfreq, ischan, sdf
    //------------------------------------------------------------

    dims.resize(1);
    dims[0] = nSpect;

    mxArray*   sfreqArr  = MexHandler::addNamedStructField(groups, "sfreq",  dims.size(), (const int*)&dims[0], 
							   DataType::DOUBLE, iGroup);
    
    mxArray*   sdfArr    = MexHandler::addNamedStructField(groups, "sdf",    dims.size(), (const int*)&dims[0], 
							   DataType::DOUBLE, iGroup);
    
    mxArray*   ischanArr = MexHandler::addNamedStructField(groups, "ischan", dims.size(), (const int*)&dims[0], 
							   DataType::INT, iGroup);

    double* sfreqPtr  = (double*) mxGetData(sfreqArr);
    double* sdfPtr    = (double*) mxGetData(sdfArr);
    int*    ischanPtr = (int*)    mxGetData(ischanArr);

    std::vector<double> sfreq = io.readDoubleVar("sfreq",  nSpect);
    std::vector<double>   sdf = io.readDoubleVar("sdf",    nSpect);
    std::vector<int>   ischan = io.readIntVar(   "ischan", nSpect);
    std::vector<int>   nschan = io.readIntVar(   "nschan", nSpect);

    unsigned nSpectValid = 0;
    for(unsigned iSpect=0; iSpect < nSpect; iSpect++) {
      sfreqPtr[iSpect]  = sfreq[iSpect];
      sdfPtr[iSpect]    = sdf[iSpect];
      ischanPtr[iSpect] = ischan[iSpect];

      if(nschan[iSpect] > 0) {
	++nSpectValid;
      }

    }

    //------------------------------------------------------------
    // Now THIS is stupid.  If spectral windows are missing (let's say
    // 2), miriad still writes the same number of spectral windows
    // (ie, 32 instead of 30), but there will be zeros in the sfreq
    // array for these windows.  However, the number of visibilities
    // returned will only be 30 windows long!  This means that we
    // cannot rely on the nspect read from the miriad file to
    // determine the number of channels
    //------------------------------------------------------------

    // Number of channels can change (if pipeline is configured to
    // return all 17 channels, for example)

    unsigned nVis   = io.nVisLastRead_;
    unsigned nChan  = nVis / nSpectValid;

    //------------------------------------------------------------
    // Add an array that is 2 antennas x nspec long, to hold the sys
    // temps for the ants comprising each baseline
    //------------------------------------------------------------

    dims.resize(2);
    dims[0] = 2;
    dims[1] = nSpect;

    mxArray* systempArr  = 
      MexHandler::addNamedStructField(groups, "systemp",  dims.size(), (const int*)&dims[0], 
				      DataType::FLOAT, iGroup);

    float* systempPtr = (float*)mxGetData(systempArr);

    // Get the system temps

    std::vector<float> systemps = io.readFloatVar("systemp", nSpect * nAnt);

    // Iterate over the return array, to extract the temps for our two
    // antennas

    for(unsigned iSpect=0; iSpect < nSpect; iSpect++) {
      
      // Index in fortran order for miriad
	
      unsigned mirInd1 = iSpect * nAnt + (*ant1-1);
      unsigned mirInd2 = iSpect * nAnt + (*ant2-1);
      unsigned matInd1 = iSpect * 2 + 0;
      unsigned matInd2 = iSpect * 2 + 1;
      
      systempPtr[matInd1] = systemps[mirInd1];
      systempPtr[matInd2] = systemps[mirInd2];
    }

    // Add an array with the same dimensionality as the group array in
    // the Miriad file

    dims.resize(2);
    dims[0] = nSpect;
    dims[1] = nChan;

    mxArray* visArr  = MexHandler::addNamedStructField(groups, "vis",  dims.size(), 
						       (const int*)&dims[0], DataType::COMPLEX_FLOAT, iGroup);

    mxArray* flagArr = MexHandler::addNamedStructField(groups, "flag", dims.size(), 
						       (const int*)&dims[0], DataType::UINT,          iGroup);

#if 0
    dims.resize(1);
    dims[0] = nVis;
    mxArray* visTestArr  = MexHandler::addNamedStructField(groups, "visTest",  dims.size(), 
							   (const int*)&dims[0], DataType::COMPLEX_FLOAT, iGroup);
#endif

    MexHandler::addNamedStringStructField(groups, "source",  source, iGroup);
    double*  raPtr  = MexHandler::addNamedDoubleStructField(groups,"ra", 1, iGroup);
    double*  decPtr = MexHandler::addNamedDoubleStructField(groups,"dec", 1, iGroup);
    double*  raApp  = MexHandler::addNamedDoubleStructField(groups,"raApp", 1, iGroup);
    double*  decApp = MexHandler::addNamedDoubleStructField(groups,"decApp", 1, iGroup);

    *raPtr  = io.readDoubleVar("ra");
    *decPtr = io.readDoubleVar("dec");
    *raApp  = io.readDoubleVar("obsra");
    *decApp = io.readDoubleVar("obsdec");

    float* rePtr = (float*)mxGetData(visArr);
    float* imPtr = (float*)mxGetImagData(visArr);
    unsigned int* flagPtr = (unsigned int*)mxGetData(flagArr);
    
#if 0
    float* reTestPtr = (float*)mxGetData(visTestArr);
    float* imTestPtr = (float*)mxGetImagData(visTestArr);
#endif

    // Set the group data

    *u   = io.preamble_[0];
    *v   = io.preamble_[1];
    *mjd = io.preamble_[2] - 2400000.5;

    *lst = io.readDoubleVar("lst");

    // Now iterate over the vis array, 

    unsigned iVis=0;
    unsigned nSpectSkip=0;

    // Iterate over ALL reported spectral windows, as reported in the
    // miriad file, even if there is no data for them

    for(unsigned iSpect=0; iSpect < nSpect; iSpect++) {

      if(!(sfreq[iSpect] > 0.0) || nschan[iSpect] == 0) {
	++nSpectSkip;
      }

      for(unsigned iChan=0; iChan < nChan; iChan++, iVis++) {

	// Modify the miriad indices by nSpsctSkip to reflect the fact
	// that miriad skips over visibility data for 'missing'
	// spectral windows

	unsigned reInd   = 2*((iSpect-nSpectSkip) * nChan + iChan);
	unsigned imInd   = 2*((iSpect-nSpectSkip) * nChan + iChan) + 1;
	unsigned flagInd =    (iSpect-nSpectSkip) * nChan + iChan;

	// But don't modify the matlab index -- we want to correctly
	// reflect in the matlab data when data for a spectral window
	// are gone

	unsigned matInd  =    iChan * nSpect + (iSpect);

	// If this was a valid spectral window, for which there is
	// actually data in the miriad file, then fill it from the
	// read arrays

	if(sfreq[iSpect] > 0.0 && nschan[iSpect] > 0) {
	  rePtr[matInd]   = io.readVisData_[reInd];
	  imPtr[matInd]   = io.readVisData_[imInd];
	  flagPtr[matInd] = io.readFlags_[flagInd];
	} else {
	  rePtr[matInd]   = 0.0;
	  imPtr[matInd]   = 0.0;
	  flagPtr[matInd] = 0; // Apparently the miriad convention is that
			       // bad data has flag=0, good data have
			       // flag=1
	}

#if 0
	reTestPtr[iVis] = io.readVisData_[2*iVis];
	imTestPtr[iVis] = io.readVisData_[2*iVis+1];
#endif
      }
    }

  }

  //------------------------------------------------------------ 
  // Close the file now that we're done reading
  //------------------------------------------------------------

  io.closeFile();

  //------------------------------------------------------------ 
  // And create structures to reflect the sources encountered in the
  // data
  //------------------------------------------------------------

  mxArray* srcs = MexHandler::addNamedStructField(hdr, "sources", sourceIndMap.size());

  unsigned iSrc=0;

  std::map<std::string, std::vector<unsigned>* >::iterator iter;
  for(iter=sourceIndMap.begin(); iter != sourceIndMap.end(); iter++, iSrc++) {

    std::string src         = iter->first;
    std::vector<unsigned>* indArr = iter->second;

    COUT("Found source: " << src << " with inds = " << indArr->size());

    MexHandler::addNamedStringStructField(srcs, "name",   src, iSrc);

    double* raptr     = MexHandler::addNamedDoubleStructField(srcs, "ra",     1, iSrc);
    double* decptr    = MexHandler::addNamedDoubleStructField(srcs, "dec",    1, iSrc);
    double* obsraptr  = MexHandler::addNamedDoubleStructField(srcs, "obsra",  1, iSrc);
    double* obsdecptr = MexHandler::addNamedDoubleStructField(srcs, "obsdec", 1, iSrc);
    double* epochptr  = MexHandler::addNamedDoubleStructField(srcs, "epoch",  1, iSrc);

    *raptr     = sourceRaMap[src];
    *decptr    = sourceDecMap[src];
    *obsraptr  = sourceObsRaMap[src];
    *obsdecptr = sourceObsDecMap[src];
    *epochptr  = sourceEpochMap[src];

    std::string intent = sourcePurposeMap[src];
    MexHandler::addNamedStringStructField(srcs, "intent", intent, iSrc);

    unsigned int* inds = MexHandler::addNamedUintStructField(srcs, "inds", indArr->size(), iSrc);

    for(unsigned iInd=0; iInd < indArr->size(); iInd++) {
      inds[iInd] = (*indArr)[iInd];
    }
  }

  // Now release any memory that was allocated                                                                                                                                                                 

  for(iter=sourceIndMap.begin(); iter != sourceIndMap.end(); iter++, iSrc++) {
    delete iter->second;
  }
}

void addHeaderData(MiriadIo& io, mxArray* hdr, unsigned nAnt, AntennaInfo& antennaInfo)
{
  std::vector<unsigned> dims;

  double* lat = MexHandler::addNamedDoubleStructField(hdr, "latitude",  1);
  double* lng = MexHandler::addNamedDoubleStructField(hdr, "longitude", 1);
  double* alt = MexHandler::addNamedDoubleStructField(hdr, "altitude",  1);

  *lat = io.readDoubleVar("latitud");
  *lng = io.readDoubleVar("longitu");
  *alt = 2196.265;

  Lla lla;
  lla.latitude_.setRadians(*lat);
  lla.longitude_.setRadians(*lng);
  lla.altitude_.setMeters(2196.265);

  //------------------------------------------------------------
  // Add Jy Per K
  //------------------------------------------------------------

    dims.resize(1);
    dims[0] = nAnt;

    mxArray* jyperkaArr = MexHandler::addNamedStructField(hdr,   "jyperka", dims.size(), (const int*)&dims[0], 
							  DataType::FLOAT);

    float* jyperkaPtr = (float*)mxGetData(jyperkaArr);

    std::vector<float>  jyperka = io.readFloatVar("jyperka", nAnt);

    for(unsigned i=0; i < nAnt; i++) {
      jyperkaPtr[i] = jyperka[i];
    }

  //------------------------------------------------------------
  // Add XYZ
  //------------------------------------------------------------

  dims.resize(1);
  dims[0] = nAnt;
    
  mxArray*   XArr  = MexHandler::addNamedStructField(hdr,   "X", dims.size(), 
						     (const int*)&dims[0], 
						     DataType::DOUBLE);
    
  mxArray*   YArr  = MexHandler::addNamedStructField(hdr,   "Y", dims.size(), 
						     (const int*)&dims[0], 
						     DataType::DOUBLE);
    
  mxArray*   ZArr  = MexHandler::addNamedStructField(hdr,   "Z", dims.size(), 
						     (const int*)&dims[0], 
						     DataType::DOUBLE);

  mxArray*   eastArr  = MexHandler::addNamedStructField(hdr,   "east", dims.size(), 
						     (const int*)&dims[0], 
						     DataType::DOUBLE);
    
  mxArray*   northArr  = MexHandler::addNamedStructField(hdr,   "north", dims.size(), 
						     (const int*)&dims[0], 
						     DataType::DOUBLE);
    
  mxArray*   upArr  = MexHandler::addNamedStructField(hdr,   "up", dims.size(), 
						     (const int*)&dims[0], 
						     DataType::DOUBLE);

  double* XPtr     = (double*)mxGetData(XArr);
  double* YPtr     = (double*)mxGetData(YArr);
  double* ZPtr     = (double*)mxGetData(ZArr);

  double* eastPtr  = (double*)mxGetData(eastArr);
  double* northPtr = (double*)mxGetData(northArr);
  double* upPtr    = (double*)mxGetData(upArr);

  // Get the antenna positions

  std::vector<double> antpos = io.readDoubleVar("antpos", 3 * nAnt);

  unsigned ind;
  Time XTime, YTime, ZTime;
  Length X, Y, Z;
  LengthTriplet xyz, enu;
  Geoid geoid;

  for(unsigned iAnt=0; iAnt < nAnt; iAnt++) {

    ind = 0 * nAnt + iAnt;

    // Convert from light travel time in nanoseconds to length

    XTime.setNanoSeconds(antpos[ind]);
    X.setLightTravelTime(XTime);
    xyz.X_ = X;
    XPtr[iAnt] = X.meters();

    ind = 1 * nAnt + iAnt;
    YTime.setNanoSeconds(antpos[ind]);
    Y.setLightTravelTime(YTime);
    xyz.Y_ = Y;    
    YPtr[iAnt] = Y.meters();

    ind = 2 * nAnt + iAnt;
    ZTime.setNanoSeconds(antpos[ind]);
    Z.setLightTravelTime(ZTime);
    xyz.Z_ = Z;
    ZPtr[iAnt] = Z.meters();

    // Finally, convert to ENU coordinates

    enu = geoid.geodeticLlaAndXyzToEnu(lla, xyz);

    // And store the result

    eastPtr[iAnt]  = enu.east_.meters();
    northPtr[iAnt] = enu.north_.meters();
    upPtr[iAnt]    = enu.up_.meters();
    
    // Store the location information in the AntennaInfo object as
    // well.  This will be used later for shadowing calculation

    antennaInfo.east_[iAnt]  = enu.east_;
    antennaInfo.north_[iAnt] = enu.north_;
    antennaInfo.up_[iAnt]    = enu.up_;
  }
}

//=======================================================================
// AntennaInfo class
//=======================================================================

const unsigned AntennaInfo::nAnt_  = 23;
const unsigned AntennaInfo::nOvro_ = 6;
const unsigned AntennaInfo::nBima_ = 9;
const unsigned AntennaInfo::nSza_  = 8;

AntennaInfo::AntennaInfo()
{
  east_.resize(nAnt_);
  north_.resize(nAnt_);
  up_.resize(nAnt_);

  diameter_.resize(nAnt_);
  sweptVolumeDiameter_.resize(nAnt_);
  elAxisHeight_.resize(nAnt_);

  for(unsigned iAnt=0; iAnt < nAnt_; iAnt++) {

    // OVRO

    if(iAnt < nOvro_) {

      diameter_[iAnt].setMeters(10.4);
      sweptVolumeDiameter_[iAnt].setMeters(2*7.043);
      elAxisHeight_[iAnt].setMeters(5.435);

      // BIMA

    } else if(iAnt < nOvro_ + nBima_) {

      diameter_[iAnt].setMeters(6.1);
      sweptVolumeDiameter_[iAnt].setMeters(2*5.6388);
      elAxisHeight_[iAnt].setMeters(5.198);

      // SZA
    } else {

      diameter_[iAnt].setMeters(3.5);
      sweptVolumeDiameter_[iAnt].setMeters(4.41442);
      elAxisHeight_[iAnt].setMeters(2.7526);
    }
  }
}

AntennaInfo::~AntennaInfo()
{
}

bool AntennaInfo::isSzaAntShadowed(unsigned iSzaAnt, Angle& az, Angle& el, bool doPrint)
 {
   if(iSzaAnt < nOvro_ + nBima_) {
     return false;
   }

   // Only iterate over non-SZA antennas

   bool isShadowed = false;
   for(unsigned iAnt=0; iAnt < nOvro_ + nBima_; iAnt++) {
     isShadowed = isShadowed | isAntShadowed(az, el, iSzaAnt, iAnt, doPrint);
   }

   return isShadowed;
 }

bool AntennaInfo::isAntShadowed(Angle& az0, Angle& el0,
				unsigned iAntTarget, unsigned iAntTest, bool doPrint)
{
  bool shadowed = false;
  double percShadowed = 0.0;

  double dU1  = 
    (up_[iAntTest].meters()   + elAxisHeight_[iAntTest].meters()) - 
    (up_[iAntTarget].meters() + elAxisHeight_[iAntTarget].meters());

  double dN1  = north_[iAntTest].meters() - north_[iAntTarget].meters();
  double dE1  = east_[iAntTest].meters()  - east_[iAntTarget].meters();

  double mag  = sqrt(dE1*dE1 + dN1*dN1 + dU1*dU1);

  // Precompute cos/sin terms

  double saz = sin(az0.radians());
  double caz = cos(az0.radians());
  double sel = sin(el0.radians());
  double cel = cos(el0.radians());

  // Now get the coordinates (direction cosines) of the pointing vector

  double dU0 = mag*sel;
  double dE0 = mag*cel*saz;
  double dN0 = mag*cel*caz;

  double cdang = (dU0*dU1 + dE0*dE1 + dN0*dN1)/(mag*mag);
  double dang  = fabs(acos(cdang));

  double anglim = fabs(asin(((diameter_[iAntTarget].meters() + 
			      sweptVolumeDiameter_[iAntTest].meters())/2)/mag));

  // We check that the projected separation of the dish center from
  // the center of the swept-volume sphere, orthogonal to the pointing
  // direction, is greater than the sum of the radii of the dish and
  // swept volume sphere

  if(dang < anglim) {
    shadowed = true;
  }

  return shadowed;
}
