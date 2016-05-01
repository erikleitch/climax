#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/NvssReader.h"

#include <iomanip>

using namespace std;
using namespace gcp::util;

// The latitude of the VLA

Angle NvssReader::vlaLat_ = 
Angle(Angle::Radians(), 0.5948);

// NVSS position biases

Angle NvssReader::raBias_ = 
Angle(Angle::ArcSec(), -0.025);

Angle NvssReader::decBias_ = 
Angle(Angle::ArcSec(),  0.113);

// NVSS flux biases
//
// DnC bias values

Flux NvssReader::dncBiasAv_ = 
Flux(Flux::Jansky(), 0.00030);

Flux NvssReader::dncBiasErr_ = 
Flux(Flux::Jansky(), 0.00030);

// D bias values

Flux NvssReader::dBiasAv_ = 
Flux(Flux::Jansky(), 0.00020);

Flux NvssReader::dBiasErr_ = 
Flux(Flux::Jansky(), 0.00020);

// RA calibration error (arcseconds)

Angle NvssReader::calRaErr_ = 
Angle(Angle::ArcSec(), 0.45);

// DEC calibration error (arcseconds)

Angle NvssReader::calDecErr_ = 
Angle(Angle::ArcSec(), 0.56);

// Amplitude calibration error (%)

double NvssReader::calAmpErr_ = 0.03;

/**.......................................................................
 * Constructor.
 */
NvssReader::NvssReader() : PtSrcFitsReader() {}

/**.......................................................................
 * Constructor.
 */
NvssReader::NvssReader(std::string catalogFile) :
  PtSrcFitsReader(catalogFile) {}

/**.......................................................................
 * Destructor.
 */
NvssReader::~NvssReader() {}

/**.......................................................................
 * Read the next entry from the catalog file
 */
PtSrcReader::Source NvssReader::parseData()
{
  // Calculate the index into the current chunk containing this source

  unsigned iSrc = iRow_ % chunkSize_;
  
  // Stick the relevant data into a source structure and return it.
  
  // Store the flux, rms and axes - we will need these for
  // error calculations below

  PtSrcReader::Source src;

  src.ra_.setDegrees(ras_[iSrc]);
  src.dec_.setDegrees(decs_[iSrc]);

  src.rawPeak_.setJy(peakFluxes_[iSrc]);
  src.peak_.setJy(peakFluxes_[iSrc]);
  src.fitMaj_.setDegrees(majorAxes_[iSrc]);
  src.fitMin_.setDegrees(minorAxes_[iSrc]);

  src.fitPa_.setDegrees(positionAngles_[iSrc]);
  src.rms_.setJy(rmsFluxes_[iSrc]);

  // Set the source name

  src.name_ = sourceNames_[iSrc];

  // NVSS catalog assumes every source is real

  src.warn_ = false;

  // Return the source

  return src;
}

/**.......................................................................
 * Read the next chunk of data from the FITS file
 *
 * Data structure from the NVSS catalog.ps file:   
 *
 *   1 RA(2000)     J2000 Right Ascension (degrees)
 *   2 DEC(2000)    J2000 declination (degrees)
 *   3 PEAK INT     Peak Stokes I (Jy/beam)
 *   4 MAJOR AX     Fitted major axis (degrees)
 *   5 MINOR AX     Fitted minior axis (degrees)
 *   6 POSANGLE     Fitted position angle of the major axis (degrees)
 *   7 Q CENTER     Interpolated Q value at position of I peak (Jy/beam)
 *   8 U CENTER     Interpolated U value at position of I peak (Jy/beam)
 *   9 P FLUX       Integrated polarized (linear) flux (Jy)
 *  10 I RMS        RMS noise in Stokes I image (Jy/beam)
 *  11 POL RMS      RMS noise in Stokes Q,U image (Jy/beam)
 *  12 RES RMS      RMS Stokes I residual after fit
 *  13 RES PEAK     Peak Stokes I residual (Jy/beam)
 *  14 RES FLUX     Integrated Stokes I residual (Jy)
 *  15 CENTER X     X pixel in FIELD of center (pixels)
 *  16 CENTER Y     Y pixel in FIELD of center (pixels)
 *  17 FIELD        Name of the 4x4 degree field (string)
 *  18 JD PROCESSED Julian date on which this entry was derived from an image (day)

 */
void NvssReader::readFitsData(long startRow, long nElement)
{
  float floatnull = 0;
  static char strnull[10] = {" "};
  int anynull;

  // Now read the data

  fits_read_col(fitsFile_, TDOUBLE,  1, startRow, 1, nElement, &floatnull,  ras_,
		&anynull, &status_);
	
  fits_read_col(fitsFile_, TDOUBLE,  2, startRow, 1, nElement, &floatnull,  decs_,
		&anynull, &status_);

  fits_read_col(fitsFile_, TFLOAT,   3, startRow, 1, nElement, &floatnull,  peakFluxes_,
		&anynull, &status_);

  fits_read_col(fitsFile_, TFLOAT,   4, startRow, 1, nElement, &floatnull,  majorAxes_,
		&anynull, &status_);

  fits_read_col(fitsFile_, TFLOAT,   5, startRow, 1, nElement, &floatnull,  minorAxes_,
		&anynull, &status_);

  fits_read_col(fitsFile_, TFLOAT,   6, startRow, 1, nElement, &floatnull,  positionAngles_,
		&anynull, &status_);

  fits_read_col(fitsFile_, TFLOAT,  10, startRow, 1, nElement, &floatnull,  rmsFluxes_,
		&anynull, &status_);

  fits_read_col(fitsFile_, TSTRING, 17, startRow, 1, nElement, &strnull,    sourceNames_,
		&anynull, &status_);
}

/**.......................................................................
 * NVSS notes say that to correct for confusion bias, the following
 * correction should be applied:
 */
void NvssReader::correctForConfusionBias(PtSrcReader::Source& src)
{
  // Separate CLEAN bias for D and DnC

  double decDeg = src.dec_.degrees();
  double biasMilliJy;
  if ((decDeg < 77.994) && (decDeg > -10.125)) {
    biasMilliJy = dBiasAv_.mJy();
  } else {
    biasMilliJy = dncBiasAv_.mJy();
  }

  double peakMilliJy = src.peak_.mJy();
  double rmsMilliJy  = src.rms_.mJy();

  peakMilliJy += (biasMilliJy - (rmsMilliJy * rmsMilliJy) / peakMilliJy);

  // And correct the peak

  src.peak_.setMilliJy(peakMilliJy);
}

/**........................................................................
 * Equation (1) from catalog.ps, for the effective SNR in a given
 * parameter
 */
double NvssReader::effSnr2(PtSrcReader::Source& src,
			   double alphaMajor, double alphaMinor)
{
  // Get the FWHMs in the same units

  double fitMinFwhm = src.fitMin_.arcsec();
  double fitMajFwhm = src.fitMaj_.arcsec();
  double resMajFwhm = src.resMaj_.arcsec();
  double resMinFwhm = src.resMin_.arcsec();

  // Now take ratios that will be reused below

  double majRat  = fitMajFwhm/resMajFwhm;
  double minRat  = fitMinFwhm/resMinFwhm;
  double fluxRat = src.rawPeak_.mJy()/src.rms_.mJy();

  double preFac  = (fitMinFwhm/resMajFwhm) * (fitMajFwhm/resMajFwhm)/4;

  return preFac * 
    pow(1.0 + 1.0/(majRat * majRat), alphaMajor) * 
    pow(1.0 + 1.0/(minRat * minRat), alphaMinor) * 
    (fluxRat * fluxRat);
}

/**.......................................................................
 * Set the error in the fitted peak flux value (see catalog.ps,
 * section 6.1)
 */
void NvssReader::setPeakError(PtSrcReader::Source& src)
{
  // Separate CLEAN bias for D and DnC

  double biasErrJy = getBiasErrInJy(src);
  double snrAmp2   = effSnr2(src, 3.0/2, 3.0/2);
  double fluxJy    = src.peak_.Jy();

  double err2      = (fluxJy * fluxJy) * (2.0/snrAmp2 + (calAmpErr_ * calAmpErr_)) + (biasErrJy * biasErrJy);

  // And set the (sqrt) value in the passed Flux object

  src.peakErr_.setJy(sqrt(err2));
}

/**.......................................................................
 * Set the error in the fitted major axis size (see catalog.ps,
 * section 6.2)
 */
void NvssReader::setMajAxisError(PtSrcReader::Source& src)
{
  double snrMaj2    = effSnr2(src, 5.0/2, 1.0/2);
  double majFwhm    = src.fitMaj_.arcsec();
  double resMajFwhm = src.resMaj_.arcsec();
  double err2       = (2.0 * majFwhm * majFwhm / snrMaj2) + (0.02 * 0.02 * resMajFwhm * resMajFwhm);

  // And set the (sqrt) value in the passed Angle object

  src.fitMajErr_.setArcSec(sqrt(err2));
}

/**.......................................................................
 * Set the error in the fitted minor axis size (see catalog.ps,
 * section 6.2)
 */
void NvssReader::setMinAxisError(PtSrcReader::Source& src)
{
  double snrMin2    = effSnr2(src, 1.0/2, 5.0/2);
  double minFwhm    = src.fitMin_.arcsec();
  double resMinFwhm = src.resMin_.arcsec();
  double err2       = (2.0 * minFwhm * minFwhm / snrMin2) + (0.02 * 0.02 * resMinFwhm * resMinFwhm);

  // And set the (sqrt) value in the passed Angle object

  src.fitMinErr_.setArcSec(sqrt(err2));
}

/**.......................................................................
 * Set the error in the fitted minor axis size (see catalog.ps,
 * section 6.4)
 */
void NvssReader::setPositionAngleError(PtSrcReader::Source& src)
{
  double snrMin2 = effSnr2(src, 1.0/2, 5.0/2);
  double majFwhm = src.fitMaj_.degrees();
  double minFwhm = src.fitMin_.degrees();
  double preFac  = (majFwhm * minFwhm);
  double sqFac   = 1.0/(majFwhm * majFwhm - minFwhm * minFwhm + 1e-20);

  double err2    = preFac * (sqFac * sqFac) * (4.0/snrMin2);

  // And set the (sqrt) value in the passed Angle object

  double deg = sqrt(err2);

  if(deg > 90)
    deg = 90;

  src.fitPaErr_.setDegrees(deg);

  // Also put the position angle into the range +- 90 degrees

  deg = src.fitPa_.degrees();

  if(deg < -90)
    deg += 180;
  if(deg > 90)
    deg -= 180;

  src.fitPa_.setDegrees(deg);
}

/**.......................................................................
 * Set the error in the fitted RA and DEC (see catalog.ps, section 6.8
 * & 6.9)
 */
void NvssReader::setPositionErrors(PtSrcReader::Source& src)
{
  double snrMaj2  = effSnr2(src, 5.0/2, 1.0/2);
  double snrMin2  = effSnr2(src, 1.0/2, 5.0/2);

  double majFwhm = src.fitMaj_.degrees();
  double minFwhm = src.fitMin_.degrees();
  double raErr   = calRaErr_.degrees();
  double decErr  = calDecErr_.degrees();

  double sinPos  = sin(src.fitPa_.radians());
  double cosPos  = cos(src.fitPa_.radians());

  double ln2     = log(2.0);

  double errX2   = (2.0 * 2.0 * majFwhm * majFwhm) / (8 * ln2 * snrMaj2);
  double errY2   = (2.0 * 2.0 * minFwhm * minFwhm) / (8 * ln2 * snrMin2);

  double raErr2    = (sinPos * sinPos * errX2 + cosPos * cosPos * errY2) + raErr * raErr;

  // And set the (sqrt) value in the passed HourAngle object

  src.raErr_.setDegrees(sqrt(raErr2)/cos(src.dec_.radians()));
  
  double decErr2   = (cosPos * cosPos * errX2 + sinPos * sinPos * errY2) + decErr * decErr;

  // And set the (sqrt) value in the passed HourAngle object

  src.decErr_.setDegrees(sqrt(decErr2));
}

/**.......................................................................
 * Set the error in the fitted RA and DEC (see catalog.ps, section 6.8
 * & 6.9)
 */
void NvssReader::
setGenericPositionErrors(PtSrcReader::Source& src, Angle& calRaErr, Angle& calDecErr)
{
  double snrMaj2  = effSnr2(src, 5.0/2, 1.0/2);
  double snrMin2  = effSnr2(src, 1.0/2, 5.0/2);

  double majFwhm = src.fitMaj_.degrees();
  double minFwhm = src.fitMin_.degrees();
  double raErr   = calRaErr.degrees();
  double decErr  = calDecErr.degrees();

  double sinPos  = sin(src.fitPa_.radians());
  double cosPos  = cos(src.fitPa_.radians());

  double ln2     = log(2.0);

  // Note that this differs from the NVSS case because they seem to
  // have inserted an ad hoc extra factor of 2 on both errX2 and errY2
  // to "bring the NVSS positions into agreement with the more
  // accurate FIRST positions".  For the general case here, I'm going
  // to stick with Condon's original derivation in PASP 109, 166.

  double errX2   = 2 * (majFwhm * majFwhm) / (8 * ln2 * snrMaj2);
  double errY2   = 2 * (minFwhm * minFwhm) / (8 * ln2 * snrMin2);

  double raErr2  = (sinPos * sinPos * errX2 + cosPos * cosPos * errY2) + raErr * raErr;

  // And set the (sqrt) value in the passed HourAngle object

  src.raErr_.setDegrees(sqrt(raErr2)/cos(src.dec_.radians()));

  double decErr2   = (cosPos * cosPos * errX2 + sinPos * sinPos * errY2) + decErr * decErr;

  // And set the (sqrt) value in the passed HourAngle object

  src.decErr_.setDegrees(sqrt(decErr2));

}

/**.......................................................................
 * NVSS main correction functions (translated from NVSSlist.f)
 */ 
void NvssReader::applyCorrections(PtSrcReader::Source& src)
{
  // Calculate the analytic expression for the point-spread function
  // (synthesized beam)

  setNvssPsFn(src);

  // Correct for positional biases

  correctForPositionBias(src);

  // Correct for confusion bias

  correctForConfusionBias(src);

  // Set the error in the peak value
  
  setPeakError(src);

  // Set the error in the major axis value

  setMajAxisError(src);

  // Set the error in the minor axis value

  setMinAxisError(src);

  // Set the error in the position angle

  setPositionAngleError(src);

  // Set the error in the position

  setPositionErrors(src);

  // Convert from fitted position to deconvolved position

  deconvolve(src);

  // And calculate the integrated flux and errors

  setIntegratedFlux(src);
}

/**.......................................................................
 * NVSS point-spread function calculation (translated from NVSSlist.f)
 */ 
void NvssReader::setNvssPsFn(PtSrcReader::Source& src)
{
  double psMin, psMaj, psPa;
  double zaRad, omega, r;
  double decRad = src.dec_.radians();

  // First calculate dirty-beam parameters for the 1 mJy/beam
  // uncleaned "pedestal".  DnC configuration:

  if(decRad <= -0.176719 || decRad > 1.361255) {

    // minor axis = 51 arcsec, in units of 45 arcsec restoring beam

    psMin = 1.133333;

    // zenith angle in radians

    zaRad = decRad - vlaLat_.radians();

    // dirty beam solid angle from NVSS raw fits
    // (/aten/ATEN_1/isot.dir/beamsize.sm, BEAMSIZE.OUT)

    omega = 1.405 + 0.065 / cos(zaRad);
    psMaj = omega / psMin;
    psPa  = 0.0;

  } else {

    // D configuration: minor axis = 54 arcsec, in units of 45 arcsec
    // restoring beam

    psMin = 1.2000;
    zaRad = decRad - vlaLat_.radians();

    // beam is elongated north-south near transit

    psPa = 0.0;

    // "zenith" band 6, observed off transit:

    if(decRad > 0.437220 && decRad <= 0.768818) {

      // average zenith angle is 27 deg

      zaRad = 0.471;

      // beam is elongated east-west instead of north-south

      psPa = 90.0;

    }
    omega = 0.91 + 0.61 / cos(zaRad);
    psMaj = omega / psMin;
  }

  // next calculate point-source response for whole source, cleaned
  // down to 1 mJy/beam equation derived in notes 010201

  r = src.peak_.mJy() - 1.0;

  psMin = (r + psMin*psMin*psMin) / (r + psMin);
  psMin = sqrt(psMin) * 45.0;

  psMaj = (r + psMaj*psMaj*psMaj) / (r + psMaj);
  psMaj = sqrt(psMaj) * 45.0;

  // Set the values in the Source object

  src.resPa_.setDegrees(psPa);
  src.resMaj_.setArcSec(psMaj);
  src.resMin_.setArcSec(psMin);

   // Force fitted value to be at least as large as the PSF

  if(src.resMaj_ > src.fitMaj_)
    src.fitMaj_ = src.resMaj_;

  if(src.resMin_ > src.fitMin_)
    src.fitMin_ = src.resMin_;
}

/**.......................................................................
 * Correct for positional biases
 */
void NvssReader::correctForPositionBias(PtSrcReader::Source& src)
{
  double decRad = src.dec_.radians();

  src.ra_.setRadians(src.ra_.radians() + raBias_.radians()/cos(decRad));
  src.dec_.setRadians(decRad + decBias_.radians());
}

/**-----------------------------------------------------------------------
 * This subroutine deconvolves a gaussian point-source response from a
 * fitted elliptical gaussian to yield the gaussian source parameters.
 *
 * This calculation is from AJ, 109, 2318 and Condon's notes of
 * 010209.  deconvolve() also determines whether each source axis is
 * significantly (98% confidence = 2.33 sigma) resolved.  
 *
 * If so, the rms error in that axis is calculated.  If not, the 98%
 * confidence upper limit to that axis is calculated and its rms error
 * is set to -1.  These error calculations are based on the old DECONV
 * subroutine and the NVSS paper (AJ, 115, 1693).
 *
 * Input arguments:
 *
 *   ptsmaj  = major-axis of pt src response (arcsec)
 *   ptsmin  = minor-axis of pt src response (arcsec)
 *   ptspa   = position angle of pt src response (DEG)
 *             (MUST BE IN RANGE -180 DEG TO +180 DEG)
 *   fitmaj  = major-axis of raw fit (arcsec)
 *   ufitmaj = rms error in fitted major axis (arcsec)
 *   fitmin  = minor-axis of raw fit (arcsec)
 *   ufitmin = rms error in fitted minor axis (arcsec)
 *   fitpa   = position angle of raw fit (DEG)
 *             (MUST BE IN RANGE -180 DEG TO +180 DEG)
 * 
 * Output arguments:
 *
 *   srcmaj  = major-axis of deconv source (arcsec)
 *             or 98% confidence upper limit
 *   usrcmaj = rms error in srcmaj if resolved,
 *             or -1 if srcmaj is an upper limit
 *   srcmin  = minor-axis of deconv source (arcsec)
 *             or 98% confidence upper limit
 *   usrcmin = rms error in srcmin if resolved,
 *             or -1 if srcmaj in an upper limit
 *   srcpa   = position angle of deconv source (deg)
 *             (RETURNED IN RANGE -90 TO +90 deg)
 */
void NvssReader::deconvolve(PtSrcReader::Source& src)
{
  double ptsmaj = src.resMaj_.arcsec();
  double ptsmin = src.resMin_.arcsec();
  double ptspa  = src.resPa_.degrees();
  double fitmaj = src.fitMaj_.arcsec();
  double ufitmj = src.fitMajErr_.arcsec();
  double fitmin = src.fitMin_.arcsec();
  double ufitmn = src.fitMinErr_.arcsec();
  double fitpa  = src.fitPa_.degrees();
  double ufitpa = src.fitPaErr_.degrees();
   
  double srcmaj, usrcmj, srcmin, usrcmn, srcpa;

  double temp,   ptprad, fitprd, phirad;
  double ptmjsq, ptmnsq, ftmjsq, ftmnsq, deltar, scmjsq;
  double scmnsq, phird2, cs2phi, sn2phi, denom,  xnum;
  double cosphi, sinphi, srcpar, ptspar, radarg, scmjmx;
  double umsmin, scmnmx, test,   umsmaj, upsmaj, upsmin;
  
  double degrad = 180.0/M_PI;
  double sigmax = 2.33;

  // fix any reversed input major and minor axe;
    
  if(ptsmaj < ptsmin)  {
    temp = ptsmaj;
    ptsmaj = ptsmin;
    ptsmin = temp;
    ptspa = ptspa + 90;
    
    if(ptspa > 180) 
      ptspa = ptspa - 360;
  }

  if (fitmaj < fitmin) {
    temp = fitmaj;
    fitmaj = fitmin;
    fitmin = temp;
    fitpa = fitpa + 90;
      
    if (fitpa > 180)
      fitpa = fitpa - 360;
  }
      
  // convert pa's to radians, in range 0 to pi

  ptprad = ptspa / degrad;
  
  if(ptprad < 0) 
    ptprad = ptprad + M_PI;
  
  fitprd = fitpa / degrad;
  
  if (fitprd < 0) 
    fitprd = fitprd + M_PI;
  
  // calculate pa difference phirad (radians) between raw fit and
  // point-source response
  
  phirad = fitprd - ptprad;
  
  // make range of phirad from -pi/2 to +pi/2
  
  if(phirad > M_PI/2)
    phirad = phirad - M_PI;
  
  if(phirad < -M_PI/2) 
    phirad = phirad + M_PI;
  
  cosphi = cos (phirad);
  sinphi = sin (phirad);
  
  // calculate squares of fwhm sizes
  
  ptmjsq = ptsmaj * ptsmaj;
  ptmnsq = ptsmin * ptsmin;
  ftmjsq = fitmaj * fitmaj;
  ftmnsq = fitmin * fitmin;
  
  //  do various special cases for which the general formula is
  //  not needed or is not well defined (division by zero)
  //
  // ptsmajsq = ptsminsq? (circular pt-src response)
  
  if(abs(ptmjsq - ptmnsq) < 1.e-3) {
    
    //  delta is the position angle difference
    //  between the deconvolved src major axis and
    //  the pt src response major axis
    
    deltar = phirad;
    scmjsq = ftmjsq - ptmjsq;
    scmnsq = ftmnsq - ptmnsq;
    
    //  fitmajsq = fitminsq? (circular raw fit) 
    
  } else if (abs (ftmjsq - ftmnsq) < 1.e-3) {
    
    deltar = M_PI / 2;
    scmjsq = ftmjsq - ptmnsq;
    scmnsq = ftmnsq - ptmjsq;
    
    // phirad = 0? (fit parallel to pt-src response)
    
  } else if(abs (phirad) < 1.e-5) {
    
    scmjsq = ftmjsq - ptmjsq;
    scmnsq = ftmnsq - ptmnsq;
    
    if(scmjsq >= scmnsq) {
      deltar = 0;
    } else {
      temp = scmjsq;
      scmjsq = scmnsq;
      scmnsq = temp;
      deltar = M_PI / 2;
    }
    
    // phirad = +- pi/2? (fit perpendicular to beam)
    
  } else if(abs (phirad - M_PI/2) < 1.e-5 ||
	    abs (phirad + M_PI/2) < 1.e-5) {
    
    scmjsq = ftmjsq - ptmnsq;
    scmnsq = ftmnsq - ptmjsq;
    deltar = phirad;
    
  } else {
    
    // end of special cases calculate deltarad
    
    phird2 = 2 * phirad;
    cs2phi = cos (phird2);
    sn2phi = sin (phird2);
    denom = (ftmjsq - ftmnsq) * cs2phi - (ptmjsq - ptmnsq);
    xnum = (ftmjsq - ftmnsq) * sn2phi;
    
    // calculate deltarad
    
    if (abs(denom) >= 1.e-5) {
      deltar = 0.5 * atan (xnum / denom);
    } else {
      if (xnum > 0) 
	deltar = M_PI / 4;
      if (xnum <= 0) 
	deltar = -M_PI / 4;
    }
    
    // range is now +- pi/4.  resolve ambiguities
    // to make range +- pi/2.
    
    if (denom < 0) {
      if (xnum > 0) 
	deltar = deltar + M_PI / 2;
      if (xnum <= 0) 
	deltar = deltar - M_PI / 2;
    }
    
    // calculate srcmajsq
    
    scmjsq = (ftmjsq - ftmnsq) * cosphi * sinphi;
    scmjsq = scmjsq / (cos (deltar) * sin (deltar));
    scmjsq = 0.5 * (scmjsq + ftmjsq + ftmnsq - (ptmjsq + ptmnsq));
    
    // calculate srcminsq
    
    scmnsq = ftmjsq + ftmnsq - (ptmjsq + ptmnsq) - scmjsq;
  }
  
  // srcmajsq < srcminsq?
  
  if (scmjsq < scmnsq) {
    temp = scmjsq;
    scmjsq = scmnsq;
    scmnsq = temp;
    deltar = deltar + M_PI / 2;
  }
  
  // if deconvolution fails (that is, the
  // square of the source size is negative,
  // set source size negative
  
  if (scmjsq > 0) {
    srcmaj = sqrt (scmjsq);
  } else {
    srcmaj = -sqrt(-scmjsq);
  }
  
  if(scmnsq > 0) {
    srcmin = sqrt (scmnsq);
  } else {
    srcmin = -sqrt(-scmnsq);
  }
  
  // calculate source position angle
  
  srcpar = deltar + ptprad;
  srcpa = srcpar * degrad;
  
  // make sure srcpa in range -90 to +90 deg
  
  if (srcpa < -90)
    srcpa = srcpa + 180;
  if (srcpa >= 90)
    srcpa = srcpa - 180;
  
  // end of fit calculation next calculate fit uncertainties
  //
  // test for significant major-axis resolution
  // (see section 5.2.4 in aj, 115, 1693)
  
  ptspar = sqrt (ptsmaj*ptsmaj*cosphi*cosphi + ptsmin*ptsmin*sinphi*sinphi);
  
  if (fitmaj < ptspar) {
    temp = ptspar;
  } else {
    temp = fitmaj;
  } 
  
  // is fit major axis > 2.33*sigma + beam
  // projected onto fit major axis?
  
  test = temp - ufitmj * sigmax - ptspar;
  
  if (test > 0) {
    
    // src major axis is significantly
    // resolved calculate rms error in src
    // major axis nvss paper, eq. 32a
    
    radarg = (temp - ufitmj) * (temp - ufitmj) - ptspar * ptspar;
    umsmaj = srcmaj - sqrt(radarg);
    radarg = (temp + ufitmj) * (temp + ufitmj) - ptspar * ptspar;
    upsmaj = sqrt(radarg) - srcmaj;
    usrcmj = 0.5 * (umsmaj + upsmaj);
    
  } else {
    
    // src major axis is not significantly resolved
    // calculate source major axis
    // 98% confidence upper limit
    
    scmjmx = temp + sigmax * ufitmj;
    radarg = scmjmx * scmjmx - ptspar * ptspar;
    scmjmx = sqrt (radarg);
    
    // set source size limit to srcmajmax
    
    srcmaj = scmjmx;
    usrcmj = -1;
  }
  
  // test for significant minor-axis resolution
  
  ptspar = sqrt (ptsmaj*ptsmaj*sinphi*sinphi + ptsmin*ptsmin*cosphi*cosphi);
  
  if (fitmin < ptspar) {
    temp = ptspar;
  } else {
    temp = fitmin;
  }
  
  test = temp - ufitmn * sigmax - ptspar;
  
  if (test > 0) {
    
    // src minor axis is significantly resolved
    // calculate rms error in src minor axis
    
    radarg = (temp - ufitmn) * (temp - ufitmn) - ptspar * ptspar;
    umsmin = srcmin - sqrt(radarg);
    radarg = (temp + ufitmn) * (temp + ufitmn) - ptspar * ptspar;
    upsmin = sqrt(radarg) - srcmin;
    usrcmn = 0.5 *(upsmin + umsmin);
  } else {
  
    // src minor axis is not significantly resolved
    // calculate source minor axis upper limt
    
    scmnmx = temp + sigmax * ufitmn;
    radarg = scmnmx * scmnmx - ptspar * ptspar;
    scmnmx = sqrt (radarg);
    
    // set source size limit to srcminmax
    
    srcmin = scmnmx;
    usrcmn = -1;
  }

  // make sure larger upper limit is called srcmaj

  if (usrcmn == -1. && usrcmj == -1.) {
    if (srcmin > srcmaj) {
      temp = srcmin;
      srcmin = srcmaj;
      srcmaj = temp;
    }
  }

  // set output values

  src.decMaj_.setArcSec(srcmaj);
  src.decMajErr_.setArcSec(usrcmj);
  src.decMin_.setArcSec(srcmin);
  src.decMinErr_.setArcSec(usrcmn);
  src.decPa_.setDegrees(srcpa);

  // If both axes were deconvolved:

  if(usrcmj > 0.0 && usrcmn > 0.0) {

    // Set the error on the decovolved pa

    ufitpa = (ufitpa > 90 ? 90 : ufitpa);
    src.decPaErr_.setDegrees(ufitpa);

    // Else can't estimate an error

  } else {
    src.decPaErr_.setDegrees(-1);
  }
}

/**.......................................................................
 * Calculate the integrated flux
 */    
void NvssReader::setIntegratedFlux(PtSrcReader::Source& src)
{
  // Store the deconvolved major and minor axis "errors" -- we will use
  // these to determine if the source is resolved

  double majerr = src.decMajErr_.arcsec();
  double minerr = src.decMinErr_.arcsec();

  // We have three cases: 
  //
  // 1) Source was resolved on one axis only
  // 2) Source was r esolved on bothe axes
  // 3) Source was unresolved

  // Source was resolved on one or more axes

  if(majerr > 0.0 || minerr > 0.0) {

    // Source was resolved on both axes

    if(majerr > 0.0 && minerr > 0.0) {
      setResolvedIntegratedFlux(src);

      // Source was resolved on one axis only

    } else {
      setPartiallyResolvedIntegratedFlux(src);
    }

    // Source was unresolved on both axes

  } else {
    setUnresolvedIntegratedFlux(src); 
  }
}

/**.......................................................................
 * Set the integrated flux and error for an unresolved source
 */
void NvssReader::setUnresolvedIntegratedFlux(PtSrcReader::Source& src)
{
  // Calculate the ratio of the fitted beams to the calculated ps
  // response
  
  double beamRat = (src.fitMaj_.arcsec()/src.resMaj_.arcsec()) * (src.fitMin_.arcsec()/src.resMin_.arcsec());
    
  // Trap under size beams

  if(beamRat < 1.0)
    beamRat = 1.0;
  
  src.int_.setJy(src.peak_.Jy() * sqrt(beamRat));

  // And set the error

  setUnresolvedIntegratedFluxError(src);
}

/**.......................................................................
 * Set the error in the integrated flux for an unresolved source (see
 * catalog.ps, section 6.6)
 */
void NvssReader::
setUnresolvedIntegratedFluxError(PtSrcReader::Source& src)
{
  // Separate CLEAN bias for D and DnC

  double biasErrJy = getBiasErrInJy(src);
  double snrAmp2   = effSnr2(src, 3.0/2, 3.0/2);

  // Use the "integrated" flux here 

  double fluxJy    = src.int_.Jy();

  // Note that this differs from the peak error by a factor of 2.0 on the snrAmp term

  double err2      = (fluxJy * fluxJy) * (
					  (1.0/snrAmp2) + 
					  (calAmpErr_ * calAmpErr_)
					  ) + (biasErrJy * biasErrJy);

  // And set the (sqrt) value in the passed Flux object

  src.intErr_.setJy(sqrt(err2));
}

/**.......................................................................
 * Set the integrated flux and error for an unresolved source
 */
void NvssReader::setResolvedIntegratedFlux(PtSrcReader::Source& src)
{
  // Calculate the ratio of the fitted beams to the calculated ps
  // response
  
  double beamRat = (src.fitMaj_.arcsec()/src.resMaj_.arcsec()) * (src.fitMin_.arcsec()/src.resMin_.arcsec());
    
  // Trap under size beams

  if(beamRat < 1.0)
    beamRat = 1.0;

  // Note that this differs from the unresolved case by a factor of
  // sqrt(beamRat)
 
  src.int_.setJy(src.peak_.Jy() * beamRat);

  // And set the error

  setResolvedIntegratedFluxError(src);
}

/**.......................................................................
 * Set the error in the integrated flux for an unresolved source (see
 * catalog.ps, section 6.5)
 */
void NvssReader::
setResolvedIntegratedFluxError(PtSrcReader::Source& src)
{
  double fluxJy  = src.int_.Jy();

  double peak    = src.peak_.Jy();
  double peakErr = src.peakErr_.Jy();

  double maj     = src.fitMaj_.arcsec();
  double majErr  = src.fitMajErr_.arcsec();

  double min     = src.fitMin_.arcsec();
  double minErr  = src.fitMinErr_.arcsec();

  // Calculate the ratio of the fitted beams to the calculated ps
  // response
  
  double beamRat = (maj/src.resMaj_.arcsec()) * (min/src.resMin_.arcsec());

  // Trap under size beams

  if(beamRat < 1.0)
    beamRat = 1.0;

  double err2    = (fluxJy * fluxJy) * (
					(peakErr/peak)*(peakErr/peak) + 
	
					(1.0/beamRat) * (
							 (majErr/maj)*(majErr/maj) + 
							 (minErr/min)*(minErr/min)
							 )
					);

  // And set the (sqrt) value in the passed Flux object

  src.intErr_.setJy(sqrt(err2));
}

/**.......................................................................
 * Set the integrated flux and error for a partially resolved source
 */
void NvssReader::
setPartiallyResolvedIntegratedFlux(PtSrcReader::Source& src)
{
  // Which axis was resolved?

  double pseudoPeakJy;
  double intFluxJy;

  // If resolved on the major axis

  if(src.decMajErr_.arcsec() > 0.0) {
    pseudoPeakJy = src.peak_.Jy() * sqrt(src.fitMin_.arcsec()/src.resMin_.arcsec());
    intFluxJy    = pseudoPeakJy * (src.fitMaj_.arcsec()/src.resMaj_.arcsec()); 
  } else {
    pseudoPeakJy = src.peak_.Jy() * sqrt(src.fitMaj_.arcsec()/src.resMaj_.arcsec());
    intFluxJy    = pseudoPeakJy * (src.fitMin_.arcsec()/src.resMin_.arcsec()); 
  }
  
  // Set the flux

  src.int_.setJy(intFluxJy);

  // And set the error

  setPartiallyResolvedIntegratedFluxError(src, pseudoPeakJy);
}

/**.......................................................................
 * Set the error in the integrated flux for an unresolved source (see
 * catalog.ps, section 6.5)
 */
void NvssReader::
setPartiallyResolvedIntegratedFluxError(PtSrcReader::Source& src, double pseudoPeakJy)
{
  // Store the integrated flux

  double intFluxJy  = src.int_.Jy();
  double biasErrJy  = getBiasErrInJy(src);

  // Error in the pseudo-peak intensity

  double snrAmp2    = effSnr2(src, 3.0/2, 3.0/2);
  double pseudoErr2 = (pseudoPeakJy * pseudoPeakJy) * (
						       (3.0/(2 * snrAmp2)) + 
						       (calAmpErr_ * calAmpErr_)
						       ) + (biasErrJy * biasErrJy);
  
  // Propagated error in the integrated flux

  double maj    = src.fitMaj_.arcsec();
  double resMaj = src.resMaj_.arcsec();
  double majErr = src.fitMajErr_.arcsec();

  double err2   = (intFluxJy * intFluxJy) * (
					     (pseudoErr2/(pseudoPeakJy * pseudoPeakJy)) + 
					     ((resMaj/maj) * ((majErr/maj) * (majErr/maj)))
					     );

  // And set the (sqrt) value in the passed Flux object

  src.intErr_.setJy(sqrt(err2));
}

/**.......................................................................
 * Return the appropriate bias error for D or DnC configuration
 */
double NvssReader::getBiasErrInJy(PtSrcReader::Source& src)
{
  // Separate CLEAN bias for D and DnC
  
  double decDeg = src.dec_.degrees();

  if ((decDeg < 77.994) && (decDeg > -10.125)) {
    return dBiasErr_.Jy();
  } else {
    return dncBiasErr_.Jy();
  }
}
