#include "gcp/util/Exception.h"
#include "gcp/util/String.h"

#include "gcp/fftutil/Antenna.h"
#include "gcp/fftutil/CarmaConfig.h"
#include "gcp/fftutil/Dft2d.h"

using namespace std;

using namespace gcp::util;

unsigned  Antenna::nxCached_ = 0;
unsigned  Antenna::nyCached_ = 0;
fftw_plan Antenna::forwardPlan_;
fftw_plan Antenna::inversePlan_;
bool      Antenna::planCached_ = false;

/**.......................................................................
 * Constructor.
 */
Antenna::Antenna() 
{
  type_  = ANT_UNKNOWN;
  idTag_ = 0;
  antNo_ = 0;

  parameterMask_    = PARAM_NONE;

  enu_.coordSystem_ = COORD_ENU;
  xyz_.coordSystem_ = COORD_XYZ;

  initializeAntennaDiameterMap();
}

/**.......................................................................
 * Destructor.
 */
Antenna::~Antenna() {}

Antenna::Antenna(const Antenna& ant)
{
  *this = (Antenna&)ant;
  initializeAntennaDiameterMap();
}

Antenna::Antenna(Antenna& ant)
{
  *this = ant;
}

void Antenna::operator=(const Antenna& ant)
{
  *this = (Antenna&)ant;
}

void Antenna::operator=(Antenna& ant)
{
  diameter_            = ant.diameter_;
  apertureEfficiency_  = ant.apertureEfficiency_;
  receiverTemperature_ = ant.receiverTemperature_;
  groundSpillover_     = ant.groundSpillover_;

  type_  = ant.type_;
  idTag_ = ant.idTag_;
  antNo_ = ant.antNo_;
  enu_   = ant.enu_;
  xyz_   = ant.xyz_;

  parameterMask_  = ant.parameterMask_;
  refGeodeticLla_ = ant.refGeodeticLla_;
  antGeodeticLla_ = ant.antGeodeticLla_;
}

/**.......................................................................
 * Return true if these antennas are the same type
 */
bool Antenna::operator==(Antenna& ant)
{
  if(type_ == ANT_UNKNOWN || ant.type_ == ANT_UNKNOWN)
    return diameter_ == ant.diameter_;
  else
    return type_ == ant.type_;
}

Antenna::AntennaType Antenna::getType(std::string type)
{
  String typeStr(type);
  typeStr = typeStr.toLower();

  if(typeStr == "ami") {
    return ANT_AMI;
  } else if(typeStr == "bima") {
    return ANT_BIMA;
  } else if(typeStr == "cbi") {
    return ANT_CBI;
  } else if(typeStr == "dasi") {
    return ANT_DASI;
  } else if(typeStr == "sza") {
    return ANT_SZA;
  } else if(typeStr == "ovro") {
    return ANT_OVRO;
  } else if(typeStr == "vla") {
    return ANT_VLA;
  } else {
    ThrowSimpleError("Unrecognized antenna type: " << type);
  }

  return ANT_UNKNOWN;
}

std::string Antenna::typeStr()
{
  return getType(type_);
}

std::string Antenna::getType(Antenna::AntennaType type)
{
  switch (type) {
  case ANT_AMI:
    return "AMI";
    break;
  case ANT_BIMA:
    return "BIMA";
    break;
  case ANT_CBI:
    return "CBI";
    break;
  case ANT_DASI:
    return "DASI";
    break;
  case ANT_SZA:
    return "SZA";
    break;
  case ANT_OVRO:
    return "OVRO";
    break;
  case ANT_VLA:
    return "VLA";
    break;
  default:
    ThrowSimpleError("Unrecognized antenna type: " << type);
    break;
  }

  return "";
}

void Antenna::setType(std::string type)
{
  setType(getType(type));
}

void Antenna::setType(AntennaType type)
{
  type_     = type; // Antenna was specified by type

  try {
    diameter_ = diameterMap_[type];
  } catch(Exception& err) {
    COUTCOLOR(err.what(), "red");
    throw err;
  }

  parameterMask_ |= PHYS_DIAM;
}

void Antenna::setDiameter(Length& diameter)
{
  type_ = ANT_DIAM; // Antenna was specified by diameter

  try {
    diameter_ = diameter;
  } catch(Exception& err) {
    COUTCOLOR(err.what(), "red");
    throw err;
  }

  parameterMask_ |= PHYS_DIAM;
}

Length Antenna::getRadius()
{
  if(!(parameterMask_ & PHYS_DIAM)) {
    ThrowError("No diameter has been specified for this antenna");
  }

  Length radius;
  radius = diameter_ / 2.0;

  return radius;
}

Length Antenna::getDiameter()
{
  if(!(parameterMask_ & PHYS_DIAM)) {
    ThrowError("No diameter has been specified for this antenna");
  }

  return diameter_;
}

void Antenna::setApertureEfficiency(Percent& apEff)
{
  apertureEfficiency_ = apEff;
  parameterMask_ |= PHYS_APEFF;
}

Percent Antenna::getApertureEfficiency()
{
  if(!(parameterMask_ & PHYS_APEFF)) {
    ThrowError("No aperture efficiency has been specified for this antenna");
  }

  return apertureEfficiency_;
}

/**.......................................................................
 * Initialize the map of known antenna diameters
 */
void Antenna::initializeAntennaDiameterMap()
{
  diameterMap_[ANT_AMI]           = Length(Length::Meters(),       4.00);
  diameterMap_[ANT_BIMA]          = Length(Length::Feet(),        20.00);
  diameterMap_[ANT_CBI]           = Length(Length::Meters(),       1.00);
  diameterMap_[ANT_DASI]          = Length(Length::Centimeters(), 25.00);
  diameterMap_[ANT_SZA]           = Length(Length::Meters(),       3.50);
  diameterMap_[ANT_OVRO]          = Length(Length::Meters(),      10.40);
  diameterMap_[ANT_VLA]           = Length(Length::Meters(),      25.00);

  secondaryDiameterMap_[ANT_BIMA] = Length(Length::Feet(),         2.00);
  secondaryDiameterMap_[ANT_OVRO] = Length(Length::Feet(),         2.00);
  secondaryDiameterMap_[ANT_SZA]  = Length(Length::Meters(),       0.35);
}

std::ostream& gcp::util::operator<<(std::ostream& os, const Antenna& antenna)
{
  return operator<<(os, (Antenna&)antenna);
}

std::ostream& gcp::util::operator<<(std::ostream& os, Antenna& antenna)
{
  os << "Type     = "   << antenna.type_     << std::endl
     << "Diameter = "   << antenna.diameter_ << std::endl
     << "Diameter was specified = " << antenna.parameterIsSet(Antenna::PHYS_DIAM) << std::endl
     << "Id Tag   = "   << antenna.idTag_    << std::endl
     << "Ant #    = "   << antenna.antNo_    << std::endl
     << "ENU Location is: " << std::endl << antenna.enu_ << std::endl
     << "XYZ Location is: " << std::endl << antenna.xyz_ << std::endl
     << "Apeff is: " << std::endl << antenna.apertureEfficiency_.percentMax1() << std::endl
     << "Trx is: " << std::endl << antenna.receiverTemperature_ << std::endl
     << "LLA is: " << std::endl << antenna.refGeodeticLla_ << std::endl;

  return os;
}

std::ostream& gcp::util::operator<<(std::ostream& os, Antenna::AntennaType& type)
{
  switch (type) {
  case Antenna::ANT_UNKNOWN:
    os << "Unknown";
    break;
  case Antenna::ANT_AMI:
    os << "AMI";
    break;
  case Antenna::ANT_BIMA:
    os << "BIMA";
    break;
  case Antenna::ANT_CBI:
    os << "CBI";
    break;
  case Antenna::ANT_DASI:
    os << "DASI";
    break;
  case Antenna::ANT_SZA:
    os << "SZA";
    break;
  case Antenna::ANT_OVRO:
    os << "OVRO";
    break;
  case Antenna::ANT_VLA:
    os << "VLA";
    break;
  default:
    os << "Unknown";
    break;
  }

  return os;
}

/**.......................................................................
 * Return the Fourier-space correlation length for this antenna
 */
double Antenna::getFourierCorrelationLength(Frequency& freq, double correlationCoeff)
{
  return Dft2d::correlationLength(diameter_, freq, correlationCoeff);
}

/**.......................................................................
 * Return the gaussian approximation to the primary beam of this antenna
 */
Image Antenna::getGaussianPrimaryBeam(unsigned nPix, Angle& size, Frequency& freq)
{
  Image beam;

  beam.setNpix(nPix);
  beam.setAngularSize(size);

  beam.createGaussianPrimaryBeam(diameter_, freq);

  return beam;
}

/**.......................................................................
 * Return a realistic aperture field for this antenna (power pattern =
 * primary beam = the square of the aperture field)
 */
Image Antenna::getRealisticApertureField(Image& image, Frequency& freq, Angle xoff, Angle yoff)
{
  if(type_ == ANT_DIAM) {
    return getGenericRealisticApertureField(image, freq, xoff, yoff);
  } else {
    return getSpecificRealisticApertureField(image, freq, xoff, yoff);
  }
}

/**.......................................................................
 * Return a generic aperture field for an antenna of this diameter,
 * assuming a generic current grading.
 */
Image Antenna::getGenericRealisticApertureField(Image& image, Frequency& freq, Angle xoff, Angle yoff)
{
  Dft2d dft;
  bool cached = planIsCached(image);

  //------------------------------------------------------------
  // If we already have a plan cached for a transform of this size,
  // don't recompute it, but set it externally instead
  //------------------------------------------------------------

  if(cached)
    dft.setPlan(forwardPlan_, inversePlan_);

  dft.zeropad(true, 4);

  dft.xAxis().setNpix(image.xAxis().getNpix());
  dft.yAxis().setNpix(image.yAxis().getNpix());

  dft.xAxis().setAngularSize(image.xAxis().getAngularSize());
  dft.yAxis().setAngularSize(image.yAxis().getAngularSize());

  //------------------------------------------------------------
  // And if no plan was cached, cache it now
  //------------------------------------------------------------

  if(!cached)
    cachePlan(image, dft.forwardPlan_, dft.inversePlan_);

  wave_.setFrequency(freq);

  // At this point, I'm considering a J0 Bessel current grading across
  // the primary, with the 20 dB point of the bessel main lobe at the
  // edge of the dish, to be 'generic'.

  dft.createJ0BesselFnDft(wave_, diameter_, 0.01);
  dft.shiftBy(xoff, yoff);

  // Compute the inverse transform to get the far-field voltage
  // pattern of this dish

  dft.computeInverseTransform();
  Image apfield = dft.getImage();

  // And return the normalized pattern

  return apfield / apfield.max();
}

/**.......................................................................
 * Return the aperture field for an antenna of this diameter, assuming
 * a current grading specific to this antenna
 */
Image Antenna::getSpecificRealisticApertureField(Image& image, Frequency& freq, Angle xoff, Angle yoff)
{
  Dft2d dft;
  bool cached = planIsCached(image);

  //------------------------------------------------------------
  // If we already have a plan cached for a transform of this size,
  // don't recompute it, but set it externally instead
  //------------------------------------------------------------

  if(cached)
    dft.setPlan(forwardPlan_, inversePlan_);

  dft.zeropad(true, 4);

  dft.xAxis().setNpix(image.xAxis().getNpix());
  dft.yAxis().setNpix(image.yAxis().getNpix());

  dft.xAxis().setAngularSize(image.xAxis().getAngularSize());
  dft.yAxis().setAngularSize(image.yAxis().getAngularSize());

  //------------------------------------------------------------
  // And if no plan was cached, cache it now
  //------------------------------------------------------------

  if(!cached)
    cachePlan(image, dft.forwardPlan_, dft.inversePlan_);

  wave_.setFrequency(freq);

  switch(type_) {

    // SZA, OVRO and BIMA antennas are all pretty much the same.  J0
    // Bessel fn current grading, with the 20 dB point at the edge of
    // the dish

  case ANT_SZA:
  case ANT_OVRO:
  case ANT_BIMA:
    innerDiameter_ = secondaryDiameterMap_[type_];
    dft.createBlockedApertureJ0BesselFnDft(wave_, innerDiameter_, diameter_, 0.01);
    break;

    // DASI was an unblocked feed horn

  case ANT_DASI:
    dft.createUniformDiskDft((diameter_ / wave_)/2);
    break;

    // Others we just assume J0 Bessel current grading with unblocked
    // apertures, since we don't in this case know the size of the
    // secondary

  default:
    dft.createJ0BesselFnDft(wave_, diameter_, 0.01);
    break;
  }

  dft.shiftBy(xoff, yoff);

  // Now get the inverse transform to get the far-field voltage
  // pattern of this antenna

  dft.computeInverseTransform();
  Image apfield = dft.getImage();

#if 0
  PgUtil::setInteractive(true);
  apfield.display();
#endif

  // And return the normalized pattern

  return apfield / apfield.max();
}

/**.......................................................................
 * Return a gaussian primary beam for this antenna
 */
Image Antenna::getGaussianPrimaryBeam(Image& image, Frequency& freq)
{
  Image beam;

  beam.xAxis().setNpix(image.xAxis().getNpix());
  beam.xAxis().setAngularSize(image.xAxis().getAngularSize());

  beam.yAxis().setNpix(image.yAxis().getNpix());
  beam.yAxis().setAngularSize(image.yAxis().getAngularSize());
  
  beam.createGaussianPrimaryBeam(diameter_, freq);

  return beam;
}

/**.......................................................................
 * Return a more realistic primary beam for this antenna
 */
Image Antenna::getRealisticPrimaryBeam(Image& image, Frequency& freq, Angle xoff, Angle yoff)
{
  if(type_ == ANT_DIAM) {
    return getGenericRealisticPrimaryBeam(image, freq, xoff, yoff);
  } else {
    return getSpecificRealisticPrimaryBeam(image, freq, xoff, yoff);
  }
}

/**.......................................................................
 * Return a generic primary beam for an antenna of this diameter,
 * assuming a generic uniform aperture field.
 */
Image Antenna::getGenericRealisticPrimaryBeam(Image& image, Frequency& freq, Angle xoff, Angle yoff)
{
  Image beam = getGenericRealisticApertureField(image, freq, xoff, yoff);

  // Square it to get the power pattern

  beam *= beam;

  // And return the normalized beam

  return beam / beam.max();
}

/**.......................................................................
 * Return a primary beam for an antenna of this diameter, assuming an
 * aperture field specific to this antenna
 */
Image Antenna::getSpecificRealisticPrimaryBeam(Image& image, Frequency& freq, Angle xoff, Angle yoff)
{
  Image beam = getSpecificRealisticApertureField(image, freq, xoff, yoff);

  // Square it to get the power pattern

  beam *= beam;

  // And return the normalized beam

  return beam / beam.max();
}

//-----------------------------------------------------------------------
// Location set methods
//-----------------------------------------------------------------------

void Antenna::setLocation(LengthTriplet& location)
{
  if(location.coordSystem_ == COORD_ENU) {

    setEast(location.east_);
    setNorth(location.north_);
    setUp(location.up_);

  } else if(location.coordSystem_ == COORD_XYZ) {

    setX(location.X_);
    setY(location.Y_);
    setZ(location.Z_);

  } else {
    ThrowError("Unrecognized coordinate system");
  }
}

void Antenna::setEast(Length& east)
{
  enu_.east_ = east;
  markParameterAsSet(LOC_EAST);
}

void Antenna::setNorth(Length& north)
{
  enu_.north_ = north;
  markParameterAsSet(LOC_NORTH);
}

void Antenna::setUp(Length& up)
{
  enu_.up_ = up;
  markParameterAsSet(LOC_UP);
}

void Antenna::setX(Length& X)
{
  xyz_.X_ = X;
  markParameterAsSet(LOC_X);
}

void Antenna::setY(Length& Y)
{
  xyz_.Y_ = Y;
  markParameterAsSet(LOC_Y);
}

void Antenna::setZ(Length& Z)
{
  xyz_.Z_ = Z;
  markParameterAsSet(LOC_Z);
}

void Antenna::setReferenceLla(Lla& refGeodeticLla)
{
  refGeodeticLla_ = refGeodeticLla;
  markParameterAsSet(LOC_REF_LLA);
}

Lla Antenna::getReferenceLla()
{
  if(!parameterIsSet(LOC_REF_LLA)) {
    ThrowError("Reference LLA has not been specified");
  }

  return refGeodeticLla_;
}

Lla Antenna::getAntennaLla()
{
  if(parameterIsSet(LOC_ANT_LLA)) {
  } else {
    Lla refLla = getReferenceLla();
    LengthTriplet enu = getEnu();
    antGeodeticLla_ = geoid_.geodeticLlaAndEnuToGeodeticLla(refLla, enu);
    markParameterAsSet(LOC_ANT_LLA);
  }

  return antGeodeticLla_;
}

LengthTriplet Antenna::getEnu()
{
  if(parameterIsSet(LOC_ENU_ALL)) {
  } else if(parameterIsSet(LOC_XYZ_ALL) && parameterIsSet(LOC_REF_LLA)) {
    enu_ = geoid_.geodeticLlaAndXyzToEnu(refGeodeticLla_, xyz_);
    markParameterAsSet(LOC_ENU_ALL);
  } else {
    ThrowError("Not enough information to compute ENU coordinates");
  }

  return enu_;
}

LengthTriplet Antenna::getXyz()
{
  if(parameterIsSet(LOC_XYZ_ALL)) {
  } else if(parameterIsSet(LOC_ENU_ALL) && parameterIsSet(LOC_REF_LLA)) {
    xyz_ = geoid_.geodeticLlaAndEnuToXyz(refGeodeticLla_, enu_);
    markParameterAsSet(LOC_XYZ_ALL);
  } else {
    ThrowError("Not enough information to compute XYZ coordinates");
  }

  return xyz_;
}

bool Antenna::hasLocation()
{
  return (parameterMask_ & LOC_ENU_ALL)==LOC_ENU_ALL || (parameterMask_ & LOC_XYZ_ALL)==LOC_XYZ_ALL;
}

void Antenna::setReceiverTemperature(Temperature& trx)
{
  receiverTemperature_ = trx;
  markParameterAsSet(TEMP_RX);
}

Temperature Antenna::getReceiverTemperature()
{
  if(!parameterIsSet(TEMP_RX))
    ThrowError("No receiver temperature has been specified for this antenna");
    
  return receiverTemperature_;
}

void Antenna::setGroundSpillover(Temperature& tgnd)
{
  groundSpillover_ = tgnd;
  markParameterAsSet(TEMP_FIXED);
}

Temperature Antenna::getGroundSpillover()
{
  if(!parameterIsSet(TEMP_FIXED)) 
    ThrowError("No ground spillover temperature has been specified for this antenna");

  return groundSpillover_;
}

bool Antenna::parameterIsSet(unsigned param)
{
  return (parameterMask_ & param)==param;
}

void Antenna::markParameterAsSet(unsigned param)
{
  parameterMask_ |= param;
}

//-----------------------------------------------------------------------
// External calling interface
//-----------------------------------------------------------------------

void Antenna::addParameters()
{
  addParameter("diameter", DataType::DOUBLE,  "The diameter of this antenna, if not determinable from the type");
  addParameter("type",     DataType::STRING,  "The type of this antenna");
  addParameter("apeff",    DataType::DOUBLE,  "The aperture efficiency");

  addParameter("trx",      DataType::DOUBLE,  "The receiver temperature to use when simulating realistic noise");
  addParameter("tfixed",   DataType::DOUBLE,  "A fixed addition to Tsys when simulating realistic noise");

  addParameter("east",     DataType::DOUBLE,  "The east coordinate of this antenna");
  addParameter("north",    DataType::DOUBLE,  "The north coordinate of this antenna");
  addParameter("up",       DataType::DOUBLE,  "The up coordinate of this antenna");
  addParameter("pad",      DataType::STRING,  "A pad id, for recognized arrays (carma32, for example)");

  addParameter("X",        DataType::DOUBLE,  "The X coordinate of this antenna");
  addParameter("Y",        DataType::DOUBLE,  "The Y coordinate of this antenna");
  addParameter("Z",        DataType::DOUBLE,  "The Z coordinate of this antenna");
}

/**.......................................................................
 * Add a parameter to our map.  If resizable=true, we create a new
 * vector of the specified type.
 */
void Antenna::setParameter(std::string name, std::string val, std::string units)
{
  ParameterManager::setParameter(name, val, units);

  if(name == "diameter") {
    Length diam;
    diam.setVal(getDoubleVal("diameter"), getParameter("diameter", true)->units_);
    setDiameter(diam);
  } else if(name == "type") {
    setType(getStringVal("type"));
  } else if(name == "apeff") {
    Percent perc;
    perc.setPercentMax1(getDoubleVal("apeff"));
    setApertureEfficiency(perc);
  } else if(name == "trx") {
    Temperature temp;
    temp.setVal(getDoubleVal("trx"), getParameter("trx", true)->units_);
    setReceiverTemperature(temp);
  } else if(name == "tfixed") {
    Temperature temp;
    temp.setVal(getDoubleVal("tfixed"), getParameter("tfixed", true)->units_);
    setGroundSpillover(temp);
  } else if(name == "east") {
    Length east;
    east.setVal(getDoubleVal("east"), getParameter("east", true)->units_);
    setEast(east);
  } else if(name == "north") {
    Length north;
    north.setVal(getDoubleVal("north"), getParameter("north", true)->units_);
    setNorth(north);
  } else if(name == "up") {
    Length up;
    up.setVal(getDoubleVal("up"), getParameter("up", true)->units_);
    setUp(up);
  } else if(name == "X") {
    Length X;
    X.setVal(getDoubleVal("X"), getParameter("X", true)->units_);
    setX(X);
  } else if(name == "Y") {
    Length Y;
    Y.setVal(getDoubleVal("Y"), getParameter("Y", true)->units_);
    setY(Y);
  } else if(name == "Z") {
    Length Z;
    Z.setVal(getDoubleVal("Z"), getParameter("Z", true)->units_);
    setZ(Z);

    // Pad is a special one -- it refers to pre-determined pad
    // locations.  Check for known arrays (currently only CARMA)

  } else if(name == "pad") {
    String padStr = getStringVal("pad");
    padStr = padStr.toLower();

    if(padStr.contains("carma")) {
      CarmaConfig cc;
      unsigned padNo = padStr.findNextInstanceOf("carma", true, " ", false, true).toInt();
      CarmaConfig::PadLocation pad = cc.getPadByNumber(padNo);
      setNorth(pad.north_);
      setEast(pad.east_);
      setUp(pad.up_);
    }

  }
}

bool Antenna::canComputeNoise()
{
  return parameterIsSet(NOISE_ALL);
}

Flux Antenna::jyPerK()
{
  if(!parameterIsSet(PHYS_DIAM | PHYS_APEFF)) {
    ThrowError("Not enough information has been specified to calculate jyPerK for this antenna");
  }

  double radius = diameter_.centimeters() / 2;
  double effectiveArea = apertureEfficiency_.percentMax1() * (M_PI * radius * radius);

  Flux gain;
  gain.setJy((2 * Constants::kBoltzCgs_ * Constants::JyPerCgs_) / effectiveArea);

  return gain;
}

PolarLengthVector Antenna::getAzEl(HourAngle* ha, Declination* dec)
{
  lla_ = getAntennaLla();
  return geoid_.geodeticLlaAndHaDecToAzElTest(lla_, *ha, *dec);
}

void Antenna::cachePlan(Image& image, fftw_plan forwardPlan, fftw_plan inversePlan)
{
  nxCached_    = image.xAxis().getNpix();
  nyCached_    = image.yAxis().getNpix();
  forwardPlan_ = forwardPlan;
  inversePlan_ = inversePlan;
  planCached_  = true;
}

bool Antenna::planIsCached(Image& image)
{
  if(planCached_ && nxCached_ == image.xAxis().getNpix() && nyCached_ == image.yAxis().getNpix()) 
    return true;
  else {
    planCached_ = false;
    return false;
  }
}
