#include "gcp/util/Exception.h"
#include "gcp/util/FitsImageReader.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/String.h"

#include <sstream>
#include <iomanip>
#include <string.h>

using namespace std;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
FitsImageReader::FitsImageReader()
{
  initialize();
}

FitsImageReader::FitsImageReader(std::string file, unsigned iHdu) 
{
  setTo(file, iHdu) ;
}

void FitsImageReader::setTo(std::string file, unsigned iHdu) 
{
  initialize();

  open(file);

  //------------------------------------------------------------
  // Get the number of keywords
  //------------------------------------------------------------

  int morekeys;
  if(ffghsp(fptr_, &nKeyword_, &morekeys, &status_) > 0)
    ThrowFitsError("Unable to get header information");

  getAxes(iHdu);

  //------------------------------------------------------------
  // Now read all relevant header information
  //------------------------------------------------------------

  getHeaderInfo();
}

void FitsImageReader::initialize()
{
  units_    = "Unknown";
  nKeyword_ = 0;
}

/**.......................................................................
 * Destructor.
 */
FitsImageReader::~FitsImageReader() {}

/**.......................................................................
 * Initialize axis pointers to NULL
 */
void FitsImageReader::initAxes()
{}

/**.......................................................................
 * Read generic header information
 */
void FitsImageReader::getHeaderInfo()
{
  std::string name, val, comment;
  
  name.resize(80);
  val.resize(80);
  comment.resize(80);

  try {
    FitsReader::getHeaderInfo();
  } catch(...) {
    status_ = 0;
  }

  for(unsigned iKey=1; iKey <= nKeyword_; iKey++) {


    for(unsigned i=0; i < 80; i++) {
      name[i]    = '\0';
      val[i]     = '\0';
      comment[i] = '\0';
    }

    if(ffgkyn(fptr_, iKey, &name[0], &val[0], &comment[0], &status_))
      ThrowFitsError("Unable to get header information");
    
    HeaderCard card(name, val, comment);

    if(name.compare(0, strlen("COMMENT"), "COMMENT")==0) {
      commentCards_.push_back(card);
    } else if(name.compare(0, strlen("HISTORY"),"HISTORY")==0) {
      historyCards_.push_back(card);
    } else {
      headerCards_.push_back(card);

      // Store header values of interest to us

      if(name.compare(0, strlen("BUNIT"),"BUNIT")==0) {
	String valStr(val);
	valStr.strip("'");
	valStr.strip('\0');
	units_ = valStr.str();
      }
    }
  }

  double ra;
  if(ffgkyd(fptr_, "OBSRA", &ra, &comment[0], &status_)) {
    ra = ra_->refVal_;
    status_ = 0;
    obsRa_.setDegrees(ra);
    infoMask_ |= ObsParameter::INFO_OBS_RA;
  }

  double dec;
  if(ffgkyd(fptr_, "OBSDEC", &dec, &comment[0], &status_)) {
    dec = dec_->refVal_;
    status_ = 0;
    obsDec_.setDegrees(dec);
    infoMask_ |= ObsParameter::INFO_OBS_DEC;
  }
}

/**.......................................................................
 * Read axis descriptions from the main HDU
 */
void FitsImageReader::getAxes(unsigned iHdu)
{
  moveToHduOffsetBy(iHdu);

  initAxes();

  std::string comment;
  comment.resize(80);

  //------------------------------------------------------------
  // Read the NAXIS keyword
  //------------------------------------------------------------

  int nfound;

  if(ffgkyj(fptr_, "NAXIS", &nAxis_, &comment[0], &status_) )
    ThrowFitsError("Unable to get number of axes");

  axes_.resize(nAxis_);

  ostringstream os;

  //------------------------------------------------------------
  // For each axis in the header, read its size and type
  //------------------------------------------------------------

  for(int iAxis=0; iAxis < nAxis_; iAxis++) {

    os.str("");

    //------------------------------------------------------------
    // Get the size
    //------------------------------------------------------------

    os << "NAXIS" << (iAxis+1);

    if(ffgkyj(fptr_, (char*)os.str().c_str(), &axes_[iAxis].n_, 
	      &comment[0], &status_) )
      ThrowFitsError("Unable to get keyword '" << os.str() << "'"); 

    //------------------------------------------------------------
    // Get the type
    //------------------------------------------------------------

    os.str("");
    os << "CTYPE" << (iAxis+1);

    axes_[iAxis].type_.resize(70);
    axes_[iAxis].typeComment_.resize(80);

    if(ffgkys(fptr_, (char*)os.str().c_str(), &axes_[iAxis].type_[0], 
	      &axes_[iAxis].typeComment_[0], &status_)) {
      axes_[iAxis].type_ = "Unknown";
      ReportSimpleFitsError("Unable to get column type for column:     " << os.str());
      status_ = 0;
    }

    //------------------------------------------------------------
    // Get the reference value
    //------------------------------------------------------------

    os.str("");
    os << "CRVAL" << (iAxis+1);

    if(ffgkyd(fptr_, (char*)os.str().c_str(), &axes_[iAxis].refVal_, 
	      &comment[0], &status_)) {
      axes_[iAxis].refVal_ = 1;
      ReportSimpleFitsError("Unable to get reference value for column: " << os.str());
    }

    //------------------------------------------------------------
    // Get the reference "pixel"
    //------------------------------------------------------------

    os.str("");
    os << "CRPIX" << (iAxis+1);

    if(ffgkyd(fptr_, (char*)os.str().c_str(), &axes_[iAxis].refPix_, 
	      &comment[0], &status_)) {
      axes_[iAxis].refPix_ = 1;
      ReportSimpleFitsError("Unable to get reference pixel for column: " << os.str());
      status_ = 0;
    }

    //------------------------------------------------------------
    // Get the pixel delta
    //------------------------------------------------------------

    os.str("");
    os << "CDELT" << (iAxis+1);

    if(ffgkyd(fptr_, (char*)os.str().c_str(), &axes_[iAxis].delta_, 
	      &comment[0], &status_)) {

      status_ = 0;
      os.str("");

      os << "CD" << iAxis+1 << "_" << iAxis+1;

      if(ffgkyd(fptr_, (char*)os.str().c_str(), &axes_[iAxis].delta_, 
		&comment[0], &status_)) {
	axes_[iAxis].delta_ = 1;
	ReportSimpleFitsError("Unable to get delta for column:           " << os.str());
	status_ = 0;
      } else {
	ReportFitsError("Found a rotation matrix -- interpreting as delta for column: " << iAxis << " but I'm not really set up to handle this!");
      }
    }

    //------------------------------------------------------------
    // Store the column number associated with this column
    //------------------------------------------------------------

    axes_[iAxis].iCol_ = iAxis+1;

    if(strstr(axes_[iAxis].type_.c_str(),"RA") != 0)
      ra_ = &axes_[iAxis];
    else if(strstr(axes_[iAxis].type_.c_str(),"DEC") != 0)
      dec_ = &axes_[iAxis];
  }

  //------------------------------------------------------------
  // Calculate the total number of pixels
  //------------------------------------------------------------

  nPixel_ = axes_[0].n_;

  for(unsigned iAxis=1; iAxis < nAxis_; iAxis++)
    nPixel_ *= axes_[iAxis].n_;
}

/**.......................................................................
 * Read all data for the current HDU
 */
void FitsImageReader::readData(Image& image)
{
  // And resize the output image

  image.data_.resize(nPixel_);

  // Now read the image data

  float nullval = 0.0;
  int anynull=0;

  if(fits_read_img(fptr_, TFLOAT, 1, nPixel_, &nullval,
		   &image.data_[0], &anynull, &status_))
    ThrowFitsError("Error reading image");
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, FitsImageReader::Image& image)
{
  return os;
}

void FitsImageReader::getAxisData(double* vals, unsigned iAxis)
{
  Axis* axis = &axes_[iAxis];

  for(unsigned i=0; i < axis->n_; i++) 
    vals[i] = axis->refVal_ + (i+1 - axis->refPix_) * axis->delta_;
}

void FitsImageReader::printAxes()
{
  COUT("Axes size = " << axes_.size());
  for(unsigned iAxis=0; iAxis < axes_.size(); iAxis++) {
    COUT(axes_[iAxis]);
  }
}

void FitsImageReader::printHeaderCards()
{
  COUT("Headers size = " << headerCards_.size());
  for(unsigned iHdr=0; iHdr < headerCards_.size(); iHdr++) {
    COUT(headerCards_[iHdr]);
  }
}
