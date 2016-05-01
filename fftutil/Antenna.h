// $Id: Antenna.h,v 1.5 2012/05/09 23:52:39 eml Exp $

#ifndef GCP_UTIL_ANTENNA_H
#define GCP_UTIL_ANTENNA_H

/**
 * @file Antenna.h
 * 
 * Tagged: Thu Jun 17 16:37:20 PDT 2010
 * 
 * @version: $Revision: 1.5 $, $Date: 2012/05/09 23:52:39 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Frequency.h"
#include "gcp/util/Geoid.h"
#include "gcp/util/Length.h"
#include "gcp/util/ParameterManager.h"
#include "gcp/util/Percent.h"
#include "gcp/util/Wavelength.h"

#include "gcp/fftutil/Image.h"

#include <map>

#include <fftw3.h>

namespace gcp {

  namespace datasets {
    class VisDataSet;
  }

  namespace util {

    class Antenna : public ParameterManager {
    public:

      //------------------------------------------------------------
      // Enumerate antenna types we know about
      //------------------------------------------------------------

      enum AntennaType {

	// These are used internally as defaults

	ANT_UNKNOWN,
	ANT_DIAM,
	ANT_SAME,

	// The following are recognized antenna types

	ANT_AMI,
	ANT_BIMA,
	ANT_CBI,
	ANT_DASI,
	ANT_SZA,
	ANT_OVRO,
	ANT_VLA
      };
     
      //------------------------------------------------------------
      // Enumerate parameters that can be specified for this antenna
      //------------------------------------------------------------

      enum {
	PARAM_NONE  = 0x0,

	LOC_EAST    = 0x1,
	LOC_NORTH   = 0x2,
	LOC_UP      = 0x4,

	LOC_X       = 0x8,
	LOC_Y       = 0x10,
	LOC_Z       = 0x20,

	LOC_REF_LLA = 0x40,
	LOC_ANT_LLA = 0x80,

	LOC_ENU_ALL = LOC_EAST | LOC_NORTH | LOC_UP,
	LOC_XYZ_ALL = LOC_X | LOC_Y | LOC_Z,

	TEMP_RX     = 0x100,
	TEMP_FIXED  = 0x200,
	PHYS_DIAM   = 0x400,
	PHYS_APEFF  = 0x800,

	// All parameters needed to calculate noise in Jy for this
	// antenna

	NOISE_ALL  = TEMP_RX | TEMP_FIXED | PHYS_DIAM | PHYS_APEFF,
      };

      /**
       * Constructor.
       */
      Antenna();

      /**
       * Destructor.
       */
      virtual ~Antenna();

      Antenna(const Antenna& ant);
      Antenna(Antenna& ant);
      void operator=(const Antenna& ant);
      void operator=(Antenna& ant);
      bool operator==(Antenna& ant);

      static AntennaType getType(std::string type);
      static std::string getType(Antenna::AntennaType type);

      std::string typeStr();
      void setType(AntennaType);
      void setType(std::string typeString);

      void setDiameter(Length& diameter);
      Length getDiameter();
      Length getRadius();

      void setApertureEfficiency(Percent& apeff);
      Percent getApertureEfficiency();

      void setReceiverTemperature(Temperature& trx);
      Temperature getReceiverTemperature();

      void setGroundSpillover(Temperature& trx);
      Temperature getGroundSpillover();

      bool canComputeNoise();

      Flux jyPerK();

      friend std::ostream& operator<<(std::ostream& os, const Antenna& antenna);
      friend std::ostream& operator<<(std::ostream& os, Antenna& antenna);
      friend std::ostream& operator<<(std::ostream& os, Antenna::AntennaType& type);

      double getFourierCorrelationLength(Frequency& freq, double correlationCoeff);

      Image getRealisticApertureField(Image& image, Frequency& freq,
				      Angle xoff = Angle(Angle::Degrees(), 0.0), Angle yoff = Angle(Angle::Degrees(), 0.0));

      Image getGenericRealisticApertureField(Image& image, Frequency& freq,
				    Angle xoff = Angle(Angle::Degrees(), 0.0), Angle yoff = Angle(Angle::Degrees(), 0.0));

      Image getSpecificRealisticApertureField(Image& image, Frequency& freq,
				    Angle xoff = Angle(Angle::Degrees(), 0.0), Angle yoff = Angle(Angle::Degrees(), 0.0));

      Image getGaussianPrimaryBeam(unsigned nPix, Angle& size, Frequency& freq);

      Image getGaussianPrimaryBeam(Image& image, Frequency& freq); 

      Image getRealisticPrimaryBeam(Image& image, Frequency& freq,
				    Angle xoff = Angle(Angle::Degrees(), 0.0), Angle yoff = Angle(Angle::Degrees(), 0.0));

      Image getGenericRealisticPrimaryBeam(Image& image, Frequency& freq,
				    Angle xoff = Angle(Angle::Degrees(), 0.0), Angle yoff = Angle(Angle::Degrees(), 0.0));

      Image getSpecificRealisticPrimaryBeam(Image& image, Frequency& freq,
				    Angle xoff = Angle(Angle::Degrees(), 0.0), Angle yoff = Angle(Angle::Degrees(), 0.0));
      
      // Methods to specify an antenna location

      void setLocation(LengthTriplet& enu);
      void setReferenceLla(Lla& geodeticLla);

      Lla getReferenceLla();
      Lla getAntennaLla();

      // Method to separately specify ENU coordinates

      void setEast(Length& east);
      void setNorth(Length& north);
      void setUp(Length& up);

      // Method to separately specify XYZ coordinates

      void setX(Length& X);
      void setY(Length& Y);
      void setZ(Length& Z);

      bool hasLocation();

      LengthTriplet getEnu();
      LengthTriplet getXyz();

      //------------------------------------------------------------
      // Parsing interface
      //------------------------------------------------------------

      void addParameters();
      void setParameter(std::string name, std::string val, std::string units="");

      //------------------------------------------------------------
      // Internal parameter mask handling
      //------------------------------------------------------------

      bool parameterIsSet(unsigned param);
      void markParameterAsSet(unsigned param);

      PolarLengthVector getAzEl(HourAngle* ra, Declination* dec);

      //------------------------------------------------------------
      // Plan caching
      //------------------------------------------------------------

      bool planIsCached(Image& image);
      void cachePlan(Image& image, fftw_plan forwardPlan, fftw_plan inversePlan);

    public:

      static unsigned nxCached_;
      static unsigned nyCached_;
      static bool planCached_;
      static fftw_plan forwardPlan_;
      static fftw_plan inversePlan_;

      friend class gcp::datasets::VisDataSet;

      Length diameter_;
      Percent apertureEfficiency_;
      Temperature receiverTemperature_;
      Temperature groundSpillover_;

      AntennaType type_;
      std::map<AntennaType, Length> diameterMap_;
      std::map<AntennaType, Length> secondaryDiameterMap_;

      unsigned idTag_;
      unsigned antNo_;

      LengthTriplet enu_;
      LengthTriplet xyz_;

      unsigned parameterMask_;

      Lla refGeodeticLla_;
      Lla antGeodeticLla_;

      void initializeAntennaDiameterMap();

      //------------------------------------------------------------
      // Used for efficient calculation
      //------------------------------------------------------------
      
      Wavelength wave_;
      Length innerDiameter_;
      Geoid geoid_;
      Lla lla_;

    }; // End class Antenna

    std::ostream& operator<<(std::ostream& os, const Antenna& antenna);
    std::ostream& operator<<(std::ostream& os, Antenna& antenna);
    std::ostream& operator<<(std::ostream& os, Antenna::AntennaType& type);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_ANTENNA_H
