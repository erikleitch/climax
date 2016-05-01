// $Id: FitsImageReader.h,v 1.3 2012/05/09 21:17:48 eml Exp $

#ifndef GCP_UTIL_FITSIMAGEREADER_H
#define GCP_UTIL_FITSIMAGEREADER_H

/**
 * @file FitsImageReader.h
 * 
 * Tagged: Tue May 29 14:24:38 PDT 2007
 * 
 * @version: $Revision: 1.3 $, $Date: 2012/05/09 21:17:48 $
 * 
 * @author username: Erik M. Leitch
 */
#include <iostream>
#include <string>
#include <vector>

#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/FitsReader.h"

namespace gcp {
  namespace util {

    class FitsImageReader : public FitsReader {
    public:

      struct Image {

	std::vector<unsigned> dims_;

	double jd_;
	HourAngle ra_;
	Declination dec_;
	Angle dx_;
	Angle dy_;

	std::vector<float> data_;

	friend std::ostream& operator<<(std::ostream& os, Image& image);
      };

      /**
       * Constructor.
       */
      FitsImageReader();
      FitsImageReader(std::string file, unsigned iHdu=0);

      /**
       * Destructor.
       */
      virtual ~FitsImageReader();

      void setTo(std::string file, unsigned iHdu=0);
      void initialize();

      std::vector<double> getData(unsigned groupPar);

      std::string units() {
	return units_;
      }

      void getHeaderInfo();
      void listHeader();
      void initAxes();
      void getAxes(unsigned iHdu);
      void readData(Image& image);

      void printAxes();
      void printHeaderCards();

      //------------------------------------------------------------
      // Image-specific
      //------------------------------------------------------------

      unsigned nPixel() {
	return nPixel_;
      }


      //------------------------------------------------------------
      // Axis methods
      //------------------------------------------------------------

      unsigned nAxis() {
	return nAxis_;
      }

      std::string axisName(unsigned iAxis) {
	return axes_[iAxis].type_;
      }

      std::string axisComment(unsigned iAxis) {
	return axes_[iAxis].typeComment_;
      }

      unsigned axisSize(unsigned iAxis) {
	return axes_[iAxis].n_;
      }

      double axisRefVal(unsigned iAxis) {
	return axes_[iAxis].refVal_;
      }

      double axisRefPix(unsigned iAxis) {
	return axes_[iAxis].refPix_;
      }

      double axisDelta(unsigned iAxis) {
	return axes_[iAxis].delta_;
      }

      std::vector<int> axisDims() {
	std::vector<int> dims;
	dims.resize(nAxis_);

	for(unsigned iAxis=0; iAxis < nAxis_; iAxis++)
	  dims[iAxis] = axes_[iAxis].n_;

	return dims;
      }

      // Return the data for axis iAxis

      void getAxisData(double* vals, unsigned iAxis);

      enum UVWUnit {
	SECONDS,
	CM
      };

    public:

      struct HeaderCard {
	std::string name_;
	std::string val_;
	std::string comment_;

	HeaderCard(std::string& name, std::string& val, std::string& comment) {
	  name_    = name;
	  val_     = val;
	  comment_ = comment;
	}

	HeaderCard(HeaderCard& card) {
	  *this = card;
	}

	HeaderCard(const HeaderCard& card) {
	  *this = card;
	}

	void operator=(const HeaderCard& card) {
	  *this = (HeaderCard&) card;
	}

	void operator=(HeaderCard& card) {
	  name_    = card.name_;
	  val_     = card.val_;
	  comment_ = card.comment_;
	}


	friend std::ostream& operator<<(std::ostream& os, HeaderCard& card) {
	  os << "Name    = '" << card.name_    << "'" << std::endl
	     << "Value   = '" << card.val_     << "'" << std::endl
	     << "Comment = '" << card.comment_ << "'" << std::endl;
	  return os;
	}

      };

      std::vector<HeaderCard> headerCards_;
      std::vector<HeaderCard> commentCards_;
      std::vector<HeaderCard> historyCards_;

      struct Par {
	std::string type_;
	unsigned iPar_;
	double scale_;
	double zero_;
      };

      long nAxis_;
      unsigned nPixel_;

      std::vector<Axis> axes_;

      Axis* complex_;
      Axis* stokes_;
      Axis* freq_;
      Axis* if_;
      Axis* ra_;
      Axis* dec_;

      std::vector<Par> pars_;
      
      Par* u_;
      Par* v_;
      Par* w_;
      Par* baseline_;
      Par* date_;

      UVWUnit uUnit_;
      UVWUnit vUnit_;
      UVWUnit wUnit_;

      std::string units_;

      int nKeyword_;
      
    }; // End class FitsImageReader

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_FITSIMAGEREADER_H
