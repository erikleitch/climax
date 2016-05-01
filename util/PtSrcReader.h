// $Id: PtSrcReader.h,v 1.12 2011/04/20 22:19:21 eml Exp $

#ifndef GCP_UTIL_PTSRCREADER_H
#define GCP_UTIL_PTSRCREADER_H

/**
 * @file PtSrcReader.h
 * 
 * Tagged: Tue Aug 15 17:53:08 PDT 2006
 * 
 * @version: $Revision: 1.12 $, $Date: 2011/04/20 22:19:21 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Angle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Exception.h"
#include "gcp/util/Flux.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/String.h"

#include "fitsio.h"

#include <vector>

namespace gcp {
  namespace util {

    class PtSrcReader {
    public:

      /**
       * Structure to encapsulate information about a source
       */
      struct Source {

	enum Type {
	  TYPE_J2000,
	  TYPE_EPHEM
	};

	// The type of this source

	Type type_;

	// The name of the source

	String name_;

	// RA of the source, and estimated error

	HourAngle ra_;
	HourAngle raErr_;

	// DEC of the source, and estimated error

	Declination dec_;
	Declination decErr_;

	// The distance of this source from a field center

	Angle distance_;

	// Specific to FIRST -- will be set to 1 (true) if the source
	// may be a sidelobe of a bright nearby source

	bool warn_;

	// Peak flux and error in the flux (technically this is an
	// intensity, in flux/beam, and not a flux)

	Flux rawPeak_;

	Flux peak_;
	Flux peakErr_;

	// Rms on the peak flux

	Flux rms_;

	// Integrated flux and error

	Flux int_;
	Flux intErr_;

	//------------------------------------------------------------
	// Fitted parameters

	// Fitted major axis (before deconvolution with restoring
	// beam)

	Angle fitMaj_;
	Angle fitMajErr_;

	// Fitted minor axis (before deconvolution with restoring
	// beam)

	Angle fitMin_;
	Angle fitMinErr_;

	// Fitted position angle (before deconvolution with restoring
	// beam)

	Angle fitPa_;
	Angle fitPaErr_;

	//------------------------------------------------------------
	// Deconvolved parameters

	// Deconvolved major axis

	Angle decMaj_;
	Angle decMajErr_;

	// Deconvolved minor axis

	Angle decMin_;
	Angle decMinErr_;

	// Deconvolved position angle

	Angle decPa_;
	Angle decPaErr_;

	//------------------------------------------------------------
	// Restoring beam parameters for this source

	Angle resMaj_;
	Angle resMin_;
	Angle resPa_;

	//------------------------------------------------------------
	// For some sources, this will record a spectral index

	double specInd_;

	//------------------------------------------------------------
	// For some sources, we will record a survey name

	std::string survey_;

	Source() {};

	Source(const Source& src) {
	  *this = (Source&)src;
	};
	
	Source(Source& src) {
	  *this = src;
	};
	
	void operator=(const Source& src) {
	  *this = (Source&)src;
	};
	  
	void operator=(Source& src) {
	  type_ = src.type_;
	  name_ = src.name_;

	  ra_ = src.ra_;
	  raErr_ = src.raErr_;

	  dec_ = src.dec_;
	  decErr_ = src.decErr_;
	  
	  distance_ = src.distance_;

	  warn_ = src.warn_;

	  rawPeak_ = src.rawPeak_;

	  peak_ = src.peak_;
	  peakErr_ = src.peakErr_;

	  rms_ = src.rms_;

	  int_ = src.int_;
	  intErr_ = src.intErr_;

	  fitMaj_ = src.fitMaj_;
	  fitMajErr_ = src.fitMajErr_;

	  fitMin_ = src.fitMin_;
	  fitMinErr_ = src.fitMinErr_;

	  fitPa_ = src.fitPa_;
	  fitPaErr_ = src.fitPaErr_;

	  resMaj_ = src.resMaj_;
	  resMin_ = src.resMin_;
	  resPa_ = src.resPa_;

	  specInd_ = src.specInd_;
	  survey_ = src.survey_;
	};

	static bool isLessThan(PtSrcReader::Source& src1, PtSrcReader::Source& src2) 
	{
	  return (src1.survey_ < src2.survey_) || (src1.survey_ == src2.survey_ && src1.name_.str() <= src2.name_.str());
	}

	static bool isEqualTo(PtSrcReader::Source& src1, PtSrcReader::Source& src2) 
	{
	  return (src1.survey_ == src2.survey_ && src1.name_.str() == src2.name_.str());
	}

	static bool isCloserThan(const PtSrcReader::Source& src1, const PtSrcReader::Source& src2) 
	{
	  return (Angle&)src1.distance_ < (Angle&)src2.distance_;
	}

	Angle distance(PtSrcReader::Source& src) {
	  Angle distance;

	  double sd1 = sin(dec_.radians());
	  double cd1 = cos(dec_.radians());

	  double sd2 = sin(src.dec_.radians());
	  double cd2 = cos(src.dec_.radians());
	    
	  double sr1 = sin(ra_.radians());
	  double cr1 = cos(ra_.radians());
	    
	  double sr2 = sin(src.ra_.radians());
	  double cr2 = cos(src.ra_.radians());
	    
	  double arg = cr1*cd1*cr2*cd2 + sr1*cd1*sr2*cd2 + sd1*sd2;

	  if(arg > 1.0) {
	    arg = 1.0;
	  }

	  double rad = acos(arg);
	    
	  distance.setRadians(rad);

	  return distance;
	};

      };

      // Constructor.

      PtSrcReader(std::string catalogFile);
      PtSrcReader();
      void initialize();

      // Destructor.

      virtual ~PtSrcReader();

      void setCatalogFile(std::string catalogFile);

      // Find sources within radius of the requested position

      std::vector<PtSrcReader::Source> findSources(HourAngle ra, Declination dec, Angle radius, 
						   Flux fMin=minFlux_, Flux fMax=maxFlux_, bool doPrint=true);

      // Return a list of sources that match the passed regexp string
      
      std::vector<PtSrcReader::Source> findSources(std::string regExpStr, bool exact=false, bool caseSensitive=true);

      // Return the number of sources

      unsigned countSources(HourAngle ra, Declination dec, Angle radius, 
			    Flux fMin=minFlux_, Flux fMax=maxFlux_);

      // Index the list of sources
      
      void indexSources();

      // Check if a position if within a given radius of the passed ra and dec

      bool checkAngle(PtSrcReader::Source& src, HourAngle& ra, Declination& dec, Angle& radius);

      // Friend funtion to print out source information

      friend std::ostream& operator<<(std::ostream& os, PtSrcReader::Source& src);

      void printHeader(std::ostream& os);

      // Catalog-specific functions to open/cose the catalog file

      virtual void openCatalogFile()  = 0;
      virtual void closeCatalogFile() = 0;

      // Catalog-specific function to read the next entry from a file

      virtual Source readNextEntry() = 0;

      // Return true if we are at the end of the catalog file

      virtual bool eof() = 0;

      // Apply any corrections to convert the catalog values

      virtual void applyCorrections(Source& src);

      // Calculate the minimal RA range we need to search

      virtual void setRaRange(HourAngle& ra, Declination& dec, Angle& radius);

      // Defaults for flux searching when none are specified

      static Flux minFlux_;
      static Flux maxFlux_;

      HourAngle raMin_;
      HourAngle raMax_;

    protected: 

      // The name of the catalog File

      std::string catalogFile_;

      //
      long sourceIndices_[25];

      // The number of sources in a catalog

      unsigned nSrc_;

      // Report an error generated by the cfitsio library
      
      void throwCfitsioError(int status);

    }; // End class PtSrcReader

    // A predicate for testing if a src equals another

    class Src_eq : public std::unary_function<PtSrcReader::Source, bool> {
      PtSrcReader::Source* src_;
    public:
      explicit Src_eq(PtSrcReader::Source* src) : src_(src) {}
	bool operator() (const PtSrcReader::Source* src) const {return (src_ != src && src_->survey_ == src->survey_ && src_->name_ == src->name_);}
    };

    // A predicate for testing if a src is lexically less than another

    struct Src_lt : public std::binary_function<PtSrcReader::Source, PtSrcReader::Source, bool> {
      bool operator() (PtSrcReader::Source& src1, PtSrcReader::Source& src2) const 
      {
	return (src1.survey_ < src2.survey_ && src1.name_.str() < src2.name_.str());
      }
    };


    std::ostream& operator<<(std::ostream& os, PtSrcReader::Source& src);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PTSRCREADER_H
