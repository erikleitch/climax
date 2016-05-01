#include "gcp/util/Constants.h"
#include "gcp/util/Exception.h"
#include "gcp/util/FileHandler.h"
#include "gcp/util/OsInfo.h"
#include "gcp/util/String.h"

#include "gcp/pgutil/UsgsReader.h"

#include <fstream>
#include <sstream>

#include <fcntl.h>

using namespace std;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
UsgsReader::UsgsReader() 
{
  isInitialized_ = false;
  byteOrder_ = LSB_FIRST;
  type_ = TYPE_UNKNOWN;
}

/**.......................................................................
 * Destructor.
 */
UsgsReader::~UsgsReader() {}

bool UsgsReader::initialized()
{
  return isInitialized_;
}

/**.......................................................................
 * Constructor.
 */
void UsgsReader::setTo(std::string directory, std::string prefix) 
{
  std::ostringstream os;

  os.str("");
  os << directory << "/" << prefix << "/" << prefix << ".hdr";
  hdrFileName_ = os.str();

  imageFileName_ = getDataFile(directory, prefix);

  COUT("Done wqith setTo typ = " << type_);
  isInitialized_ = true;
}

std::string UsgsReader::getDataFile(std::string directory, std::string prefix)
{
  std::ostringstream os, osPref;

  osPref.str("");
  osPref << directory << "/" << prefix << "/" << prefix;

  os.str("");
  os << osPref.str() << ".flt";
  if(FileHandler::fileExists(os.str())) {
    COUT("float file exists");
    type_ = TYPE_FLT;
    return os.str();
  }

  os.str("");
  os << osPref.str() << ".dem";
  if(FileHandler::fileExists(os.str())) {
    COUT("dem file exists");
    type_ = TYPE_DEM;
    return os.str();
  }

  ThrowSimpleColorError("Unable to determine the type of dataset", "red");

  return os.str();
}

void UsgsReader::parseHeaderInfo()
{
  if(!initialized()) {
    ThrowError("No image data have been specified: use UsgsReader::setTo()");
  }

  switch (type_) {
  case TYPE_FLT:
    COUT("Parsing float header");
    parseFltHeaderInfo();
    break;
  case TYPE_DEM:
    COUT("Parsing dem header");
    parseDemHeaderInfo();
    break;
  default:
    ThrowError("Unrecognized type");
    break;
  }
}

void UsgsReader::parseFltHeaderInfo()
{
  ifstream inFile(hdrFileName_.c_str());

  if (!inFile) {
    ThrowError("Can't open input file " << hdrFileName_);
  }

  std::ostringstream os;

  std::string s;
  while(getline(inFile, s)) {
    os << s << std::endl;
  }

  inFile.close();

  String str(os.str());

  unsigned nCols, nRows;
  Angle cell;
  Angle lat, lon;

  str.findNextStringSeparatedByChars(" \n");
  nCols = str.findNextStringSeparatedByChars(" \n").toInt();

  str.findNextStringSeparatedByChars(" \n");
  nRows = str.findNextStringSeparatedByChars(" \n").toInt();

  str.findNextStringSeparatedByChars(" \n");
  lon.setDegrees(str.findNextStringSeparatedByChars(" \n").toFloat());

  str.findNextStringSeparatedByChars(" \n");
  lat.setDegrees(str.findNextStringSeparatedByChars(" \n").toFloat());

  str.findNextStringSeparatedByChars(" \n");
  cell.setDegrees(str.findNextStringSeparatedByChars(" \n").toFloat());

  str.findNextStringSeparatedByChars(" \n");
  noDataVal_ = str.findNextStringSeparatedByChars(" \n").toFloat();

  str.findNextStringSeparatedByChars(" \n");
  byteOrder_ = (str.findNextStringSeparatedByChars(" \n").contains("LSBFIRST")) ? LSB_FIRST : USB_FIRST; 

  Angle xSize, ySize;
  xSize.setDegrees(nCols * cell.degrees());
  ySize.setDegrees(nRows * cell.degrees());

  image_.initialize(xSize, ySize, nCols, nRows);
  COUT("Setting lat =" << lat << " lon = " << lon);
  image_.setLatLng(lat, lon);
}

void UsgsReader::parseDemHeaderInfo()
{
  ifstream inFile(hdrFileName_.c_str());

  if(!inFile) {
    ThrowError("Can't open input file " << hdrFileName_);
  }

  std::ostringstream os;

  std::string s;
  while(getline(inFile, s)) {
    os << s << std::endl;
  }

  inFile.close();

  String str(os.str());

  unsigned nCols, nRows;
  Angle xCell, yCell;
  Angle lat, lng;

  String tok, val;

  do {
    getNextKeywordPair(str, tok, val);

    if(tok.contains("NROWS")) {
      nRows = val.toInt();
    } else if(tok.contains("NCOLS")) {
      nCols = val.toInt();
    } else if(tok.contains("NBITS")) {
      nBitPerPixel_ = val.toInt();
      COUT("Reading nBitPerPipxel = " << nBitPerPixel_);

    } else if(tok.contains("NODATA")) {
      noDataVal_ = val.toInt();
    } else if(tok.contains("ULXMAP")) {
      lng.setDegrees(val.toDouble());
    } else if(tok.contains("ULYMAP")) {
      lat.setDegrees(val.toDouble());
    } else if(tok.contains("XDIM")) {
      xCell.setDegrees(val.toDouble());
    } else if(tok.contains("YDIM")) {
      yCell.setDegrees(val.toDouble());
    } else if(tok.contains("BYTEORDER")) {
      byteOrder_ = val.contains("M") ? USB_FIRST : LSB_FIRST;
    }

  } while(!tok.isEmpty());

  COUT("Image has size: " << nCols << " x " << nRows);
  Angle xSize, ySize;
  xSize.setDegrees(nCols * xCell.degrees());
  ySize.setDegrees(nRows * yCell.degrees());

  image_.initialize(xSize, ySize, nCols, nRows);
  image_.setLatLng(lat, lng);
}

void UsgsReader::readImageData()
{
  COUT("Inside rid with type = " << type_);
  if(!initialized()) {
    ThrowError("No image data have been specified: use UsgsReader::setTo()");
  }

  int fd = -1;

  fd = ::open(imageFileName_.c_str(), O_RDONLY);

  if(fd < 0) {
    ThrowSysError("in open(): ");
  }

  std::vector<float> floatImage;
  std::vector<short> shortImage;
  std::vector<int> intImage;

  floatImage.resize(image_.data_.size());

  if(type_ == TYPE_FLT) {
    ::read(fd, (unsigned char*)&floatImage[0], 4 * image_.data_.size());
  } else {

    bool swap = needsSwap();

    if(nBitPerPixel_ == 16) {

      shortImage.resize(image_.data_.size());
      ::read(fd, (unsigned char*)&shortImage[0], 2 * image_.data_.size());
      
      for(unsigned i=0; i < image_.data_.size(); i++) {

	if(swap)
	  reverseBytes(shortImage[i]);

	floatImage[i] = (float) shortImage[i];
      }
    } else     if(nBitPerPixel_ == 32) {
      intImage.resize(image_.data_.size());
      ::read(fd, (unsigned char*)&intImage[0], 4 * image_.data_.size());
      for(unsigned i=0; i < image_.data_.size(); i++)
	floatImage[i] = (float) intImage[i];
    } else {
      ThrowError("Unsupported nBitPerPixel_");
    }
  }

  if(fd > 0)
    ::close(fd);

  std::valarray<float>& imageData = image_.data_;
  unsigned nx = image_.xAxis().getNpix();
  unsigned ny = image_.yAxis().getNpix();

  double dVal;

  unsigned indexFrom, indexTo;
  for(unsigned ix=0; ix < nx; ix++) {
    for(unsigned iy=0; iy < ny; iy++) {
      
      indexFrom = iy * nx + ix;
      indexTo   = (ny - 1 - iy) * nx + ix;

      imageData[indexTo] = floatImage[indexFrom] * (40.0/12)/1000;

      dVal = floatImage[indexFrom];

      if(dVal == noDataVal_)
	image_.valid_[indexTo] = false;
      else
	image_.valid_[indexTo] = true;

    }
  }
  
  imageReceived_ = true;
  image_.hasData_ = true;
  image_.setUnits(Unit::UNITS_KFEET);

  image_.setInvalidPixelsTo(0.0);

  return;
}

void UsgsReader::display()
{
  PgUtil::setWnad(true);
  image_.xAxis().setScale(Constants::milesPerDegree_, "Miles");
  image_.yAxis().setScale(Constants::milesPerDegree_, "Miles");
  image_.display();
}

void UsgsReader::getNextKeywordPair(String& str, String& tok, String& val)
{
  tok = str.findNextStringSeparatedByChars(" \n");
  val = str.findNextStringSeparatedByChars(" \n");
}

void UsgsReader::reverseBytes(short& s)
{
  short sc = s;
  char* sp = (char*)&s;
  char* scp = (char*)&sc;

  sp[0] = scp[1];
  sp[1] = scp[0];
}

bool UsgsReader::needsSwap()
{
  return (byteOrder_ == USB_FIRST && OsInfo::isLittleEndian()) ||
    (byteOrder_ == LSB_FIRST && OsInfo::isBigEndian());
}
