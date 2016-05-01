#include "gcp/util/Exception.h"
#include "gcp/util/FitsImageReader.h"
#include "gcp/util/FitsImageWriter.h"
#include "gcp/util/FileHandler.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
FitsImageWriter::FitsImageWriter() 
{
  initialize();
}

FitsImageWriter::FitsImageWriter(std::string fileName) 
{
  initialize();
  open(fileName);
}

void FitsImageWriter::initialize()
{
  status_ = 0;
  fptr_   = 0;
  wasOpened_ = false;
  hasBitPix_ = false;
  hasAxes_ = false;
}

void FitsImageWriter::open(std::string fileName)
{
  fileName_ = fileName;

  if(FileHandler::fileExists(fileName)) {
    ThrowSimpleColorError(std::endl << "Unable to open file " << fileName << ":" 
			  << std::endl << std::endl << "file already exists", "red");
  }

  if(ffinit(&fptr_, fileName.c_str(), &status_) > 0) {
    COUTCOLOR("Unable to open file: " << fileName, "red");
    ThrowSimpleFitsErrorColor("Unable to open file: " << fileName, "red");
  }

  wasOpened_ = true;
}

void FitsImageWriter::close()
{
  // Only close the file if it was opened by this object
  
  if(wasOpened_) {
    if(fptr_) {
      if(ffclos(fptr_, &status_) > 0)
	ReportFitsError("Error closing file");
      fptr_ = 0;
    }
  }
}

/**.......................................................................
 * Destructor.
 */
FitsImageWriter::~FitsImageWriter() 
{
  FitsImageWriter::close();
}

void FitsImageWriter::writeStandardKeys()
{
  if(!(hasBitPix_)) {
    ThrowError("Bits per pixel has not been specified: use setBitsPerPixel()");
  }

  if(!(hasAxes_)) {
    ThrowError("Axes have not been specified: use setAxes()");
  }

  writeStandardKeys(bitPix_, axisLengths_);
}

void FitsImageWriter::writeStandardKeys(int bitsPerPixel, std::vector<long>& axisLengths)
{
  if(ffphps(fptr_, bitsPerPixel, axisLengths.size(), &axisLengths[0], &status_)) {
    ThrowFitsError("Error writing standard keys");
  }
}

void FitsImageWriter::writeImageKeys(int bitsPerPixel, std::vector<long>& axisLengths)
{
  if(fits_create_img(fptr_, bitsPerPixel, axisLengths.size(), &axisLengths[0], &status_)) {
    ThrowFitsError("Error writing standard keys");
  }
}

void FitsImageWriter::putKey(std::string name, std::string val, std::string comment)
{
  if(ffpky(fptr_, TSTRING, (char*)name.c_str(), (void*)val.c_str(), (char*)comment.c_str(), &status_)) {
    ThrowFitsError("Error writing key: " << name);
  }
}

void FitsImageWriter::putKey(std::string name, bool val, std::string comment)
{
  if(ffpky(fptr_, TLOGICAL, (char*)name.c_str(), (void*)&val, (char*)comment.c_str(), &status_)) {
    ThrowFitsError("Error writing key: " << name);
  }
}

void FitsImageWriter::putKey(std::string name, char val, std::string comment)
{
  if(ffpky(fptr_, TBYTE, (char*)name.c_str(), (void*)&val, (char*)comment.c_str(), &status_)) {
    ThrowFitsError("Error writing key: " << name);
  }
}

void FitsImageWriter::putKey(std::string name, short val, std::string comment)
{
  if(ffpky(fptr_, TSHORT, (char*)name.c_str(), (void*)&val, (char*)comment.c_str(), &status_)) {
    ThrowFitsError("Error writing key: " << name);
  }
}

void FitsImageWriter::putKey(std::string name, int val, std::string comment)
{
  if(ffpky(fptr_, TINT, (char*)name.c_str(), (void*)&val, (char*)comment.c_str(), &status_)) {
    ThrowFitsError("Error writing key: " << name);
  }
}

void FitsImageWriter::putKey(std::string name, long val, std::string comment)
{
  if(ffpky(fptr_, TLONG, (char*)name.c_str(), (void*)&val, (char*)comment.c_str(), &status_)) {
    ThrowFitsError("Error writing key: " << name);
  }
}

void FitsImageWriter::putKey(std::string name, float val, std::string comment)
{
  if(ffpky(fptr_, TFLOAT, (char*)name.c_str(), (void*)&val, (char*)comment.c_str(), &status_)) {
    ThrowFitsError("Error writing key: " << name);
  }
}

void FitsImageWriter::putKey(std::string name, double val, std::string comment)
{
  if(ffpky(fptr_, TDOUBLE, (char*)name.c_str(), (void*)&val, (char*)comment.c_str(), &status_)) {
    ThrowFitsError("Error writing key: " << name);
  }
}


void FitsImageWriter::setBitsPerPixel(unsigned bitPix)
{
  bitPix_ = bitPix;
  hasBitPix_ = true;
}

void  FitsImageWriter::setAxes(std::vector<long>& axisLengths)
{
  axisLengths_ = axisLengths;
  hasAxes_ = true;
}

void FitsImageWriter::writeData(std::vector<float>& data)
{
  if(fits_write_img(fptr_, TFLOAT, 1, data.size(), &data[0], &status_))
    ThrowFitsError("Error writing image");
}

void FitsImageWriter::writeData(std::valarray<float>& data)
{
  if(fits_write_img(fptr_, TFLOAT, 1, data.size(), &data[0], &status_))
    ThrowFitsError("Error writing image");
}
