#include "gcp/util/Exception.h"
#include "gcp/util/FirstReader.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/String.h"

#include <iomanip>

using namespace std;
using namespace gcp::util;

Angle FirstReader::eps_ = 
Angle(Angle::ArcSec(), 0.1);

// The restoring beam for FIRST is 6.4x5.4 arcseconds (in the north),
// according to the FIRST web pages

Angle FirstReader::northResMaj_ = 
Angle(Angle::ArcSec(), 5.4);

Angle FirstReader::northResMin_ = 
Angle(Angle::ArcSec(), 5.4);

// Below +4:33:21, the beam is elliptical, 6.4x5.4

Angle FirstReader::medianDecLimit_ = 
Angle("4:33:21");

Angle FirstReader::medianResMaj_ = 
Angle(Angle::ArcSec(), 6.4);

Angle FirstReader::medianResMin_ = 
Angle(Angle::ArcSec(), 5.4);

// Below -2:30:25, the beam is elliptical, 6.4x5.4

Angle FirstReader::southDecLimit_ = 
Angle("-2:30:25");

Angle FirstReader::southResMaj_ = 
Angle(Angle::ArcSec(), 6.8);

Angle FirstReader::southResMin_ = 
Angle(Angle::ArcSec(), 5.4);

/**.......................................................................
 * Constructor.
 */
FirstReader::FirstReader() {}

/**.......................................................................
 * Destructor.
 */
FirstReader::~FirstReader() {}

/**.......................................................................
 * Apply FIRST corrections
 */
void FirstReader::applyCorrections(PtSrcReader::Source& src)
{
  // Set synthesized beam parameters

  setRestoringBeam(src);

  // Now set up errors

  // Empirical expression for the SNR (from FIRST catalog web pages)

  double snr    = (src.peak_.mJy() - 0.25) / src.rms_.mJy();

  // Construct the average beam for uncertainty calculations

  double avBeam = sqrt(src.resMin_.arcsec() * src.resMaj_.arcsec());

  // FIRST gives an empirical expression for the positional error
  // along the axes of the fitted ellipse.  Since these will all be
  // tiny anyway, I'll just use a single positional error estimated
  // from the mean of the two axes.

  setPositionErrors(src);

  // Error in the sizes is given by the following empirical formula

  if(src.decMin_ < eps_) {
    src.decMin_ = src.fitMin_;
    src.decMinErr_.setArcSec(-1);
  } else {
    setSizeError(src, src.decMin_, src.decMinErr_);
  }

  if(src.decMaj_ < eps_) {
    src.decMaj_ = src.fitMaj_;
    src.decMajErr_.setArcSec(-1);
  } else {
    setSizeError(src, src.decMaj_, src.decMajErr_);
  }

  // No sensible errors are given for either the peak value or the
  // integrated flux.  However, the FIRST catalog description states
  // that "For bright sources the accuracies of Fpeak and Fint are
  // limited to about 5% by systematic effects".  I will take this
  // latter as the error on both, as a conservative estimate of the
  // error on the flux values.

  src.peakErr_.setJy(sqrt((0.05 * 0.05) * src.peak_.Jy() * src.peak_.Jy() +
			 src.rms_.Jy() * src.rms_.Jy()));

  src.intErr_.setJy(sqrt((0.05 * 0.05) * src.int_.Jy() * src.int_.Jy() +
			 src.rms_.Jy() * src.rms_.Jy()));

  // Set all other errors to zero - we don't know what they are:

  src.fitMinErr_.setArcSec(0.0);
  src.fitMajErr_.setArcSec(0.0);
  src.fitPaErr_.setDegrees(0.0);
}

/**.......................................................................
 * Set the size error for an axis
 */
void FirstReader::
setSizeError(PtSrcReader::Source& src, Angle& axis, Angle& error)
{
  double snr    = getSnr(src);
  Angle avBeam  = getAvBeam(src);
  double psiF   = sqrt(avBeam.arcsec() * avBeam.arcsec() + axis.arcsec() * axis.arcsec());
  double sigmaF = psiF * (1.0/snr + 1.0/75);
  double rat    = axis.arcsec()/avBeam.arcsec();
  double rat2   = rat * rat;
  double fac    = avBeam.arcsec()/(2 * sigmaF);
  double prefac = sqrt(avBeam.arcsec() * sigmaF/2);

  error.setArcSec(prefac * (1.0 + rat2 + sqrt(fac) * rat2*rat2)/(1.0 + fac * rat2 * rat2));
}

/**.......................................................................
 * Set the appropriate restoring beam
 */
void FirstReader::setRestoringBeam(PtSrcReader::Source& src)
{
  if(src.dec_ < southDecLimit_) {
    src.resMaj_ = southResMaj_;
    src.resMin_ = southResMin_;
  } else if(src.dec_ < medianDecLimit_) {
    src.resMaj_ = medianResMaj_;
    src.resMin_ = medianResMin_;
  } else {
    src.resMaj_ = northResMaj_;
    src.resMin_ = northResMin_;
  }
}

/**.......................................................................
 * Get the SNR for this source
 */
double FirstReader::getSnr(PtSrcReader::Source& src)
{
  return (src.peak_.mJy() - 0.25) / src.rms_.mJy();
}

/**.......................................................................
 * Construct the average beam for uncertainty calculations
 */
Angle FirstReader::getAvBeam(PtSrcReader::Source& src)
{
  Angle avBeam;
  avBeam.setArcSec(sqrt(src.resMin_.arcsec() * src.resMaj_.arcsec()));
  return avBeam;
}

/**.......................................................................
 * Set position errors
 */
void FirstReader::setPositionErrors(PtSrcReader::Source& src)
{
  Angle avBeam = getAvBeam(src);

  double snr = getSnr(src);
  
  double minPosErr = sqrt(avBeam.arcsec() * avBeam.arcsec() + 
			  src.decMin_.arcsec() * src.decMin_.arcsec()) * 
    (1.0/snr + 1.0/20);

  double majPosErr = sqrt(avBeam.arcsec() * avBeam.arcsec() + 
			  src.decMaj_.arcsec() * src.decMaj_.arcsec()) * 
    (1.0/snr + 1.0/20);
}

