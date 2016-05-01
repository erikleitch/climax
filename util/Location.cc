#include "gcp/util/Location.h"
#include "gcp/util/Coordinates.h"
#include "gcp/util/Debug.h"

using namespace std;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Location::Location() 
{
  lacking_ = ALL;
  changed_ = true;
  
  fiducialAltitude_ = 0.0;
  
  uen_.resize(3);
  
  uen_[0] = uen_[1] = uen_[2] = 0.0;
  
  actualAltitude_ = 0.0;
  
  geocentricXyz_.resize(3);
  
  geocentricXyz_[0] = 0.0;
  geocentricXyz_[1] = 0.0;
  geocentricXyz_[2] = 0.0;
  
  topocentricXyz_.resize(3);
  
  topocentricXyz_[0] = 0.0;
  topocentricXyz_[1] = 0.0;
  topocentricXyz_[2] = 0.0;
}

Location::Location(Location& location) 
{
}

Location::Location(const Location& location) 
{
}

/**.......................................................................
 * Destructor.
 */
Location::~Location() 
{
}

/**.......................................................................
 * A method to set a fiducial LLA point
 */
void Location::setFiducialSite(Angle longitude, Angle latitude, double altitude)
{
  fiducialLongitude_ = longitude;
  fiducialLatitude_  = latitude;
  fiducialAltitude_  = altitude;
  
  lacking_ &= ~SITE;
  changed_  = true;
  
  // If the site has changed, we need to update coordinates
  
  updateCoordinates();
}

void Location::setUp(Length up)
{
  uen_[0] = up.meters();

  lacking_ &= ~UP;
  changed_  = true;
  
  // If the offset has changed, we need to update coordinates
  
  updateCoordinates();
}

void Location::setEast(Length east)
{
  uen_[1] = east.meters();

  lacking_ &= ~EAST;
  changed_  = true;
  
  // If the offset has changed, we need to update coordinates
  
  updateCoordinates();
}

void Location::setNorth(Length north)
{
  uen_[1] = north.meters();

  lacking_ &= ~NORTH;
  changed_  = true;
  
  // If the offset has changed, we need to update coordinates
  
  updateCoordinates();
}

/**.......................................................................
 * A method to set an UEN offset relative to the fiducial
 */
void Location::setOffset(double up, double east, double north)
{
  uen_[0] = up;
  uen_[1] = east;
  uen_[2] = north;
  
  lacking_ &= ~LOCATION;
  changed_  = true;
  
  // If the offset has changed, we need to update coordinates
  
  updateCoordinates();
}

/**.......................................................................
 * Return true if a site and offset have ben set for this
 * location
 */
bool Location::canLocate()
{
  return lacking_ == NONE;
}

/**.......................................................................
 * Update coordinate representations of this location
 */
void Location::updateCoordinates()
{
  // Update the (L, L, A) representation of this location
  
  Vector<double> lla = Coordinates::
    llaAndUenToLla(longitude(true), latitude(true), altitude(true), 
		   up(), east(), north());

  actualLongitude_.setRadians(lla[0]);
  actualLatitude_.setRadians(lla[1]);
  actualAltitude_ = lla[2];
  
  // Needed for tropospheric calculation

  Length alt;
  alt.setMeters(lla[2]);

  // Needed for refraction correction for tropospheric calculation

  Angle lat;
  lat.setRadians(lla[1]);

  // Update the geocentric (X, Y, Z) representation of this location
  
  geocentricXyz_  = Coordinates::laAndUenToXyz(latitude(true), altitude(true), 
					       up(), east(), north(), true);
  
  // Update the topocentric (X, Y, Z) representation of this location
  
  topocentricXyz_ = Coordinates::laAndUenToXyz(latitude(true), altitude(true), 
					       up(), east(), north(), false);
}

/**.......................................................................
 * Get geometric delay for an Ha Dec source position, in
 * nanoseconds.
 */
Delay Location::geometricDelay(Location* refDLoc,
			       bool doMotionCorrection)
{
  // Pass earth-centered coordinates, e.g., X(true), in case we are
  // doing the motion correction
  
  return Coordinates::
    getGeometricDelay(ha_, dec_,
		      X(true), Y(true), Z(true),
		      refDLoc->X(true), refDLoc->Y(true), refDLoc->Z(true), 
		      doMotionCorrection);
}

/**.......................................................................
 * Get geometric delay for an Az El source position, in
 * nanoseconds
 */
Delay Location::geometricDelay(Angle az, Angle el,
			       Location* refDLoc)
{
  // If we have a valid reference location, it doesn't matter which
  // coordinates we pass, earth-centered or surface-relative, since
  // all positions are relative.  Note however that if the reference
  // location hasn't been specified, it will default to (0,0,0), in
  // which case we want to make sure that the antenna coordinates are
  // surface-relative.
  //
  // Make sure we use the actual latitude (latitude(false)) of this
  // location, and not the fiducial latitude.
  
  return Coordinates::
    getGeometricDelay(latitude(false), az, el, 
		      X(false)-refDLoc->X(false), 
		      Y(false)-refDLoc->Y(false), 
		      Z(false)-refDLoc->Z(false));
}

/**.......................................................................
 * Convert mjd to lst for the actual location of this antenna
 */
HourAngle Location::getLst(double mjdUtc)
{
  // Use the actual position of this antenna

  return astrom_.mjdUtcToLst(mjdUtc, longitude(false));
}

/**.......................................................................
 * Convert to Ha for the actual location of this antenna
 */
HourAngle Location::getHa(double mjdUtc, HourAngle ra)
{
  return getLst(mjdUtc) - ra;
}


