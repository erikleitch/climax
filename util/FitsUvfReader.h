// $Id: FitsUvfReader.h,v 1.4 2012/05/09 21:17:48 eml Exp $

#ifndef GCP_UTIL_FITSUVFREADER_H
#define GCP_UTIL_FITSUVFREADER_H

/**
 * @file FitsUvfReader.h
 * 
 * Tagged: Tue May 29 14:24:38 PDT 2007
 * 
 * @version: $Revision: 1.4 $, $Date: 2012/05/09 21:17:48 $
 * 
 * @author Erik Leitch.
 */
#include <iostream>
#include <string>
#include <vector>

#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/FitsReader.h"

namespace gcp {
  namespace util {

    class FitsUvfReader : public FitsReader {
    public:

      struct Vis {

	std::vector<unsigned> dims_;

	double u_;
	double v_;
	double w_;
	unsigned int baseline_;
	double jd_;
	unsigned int source_;

	std::vector<unsigned> index_;
	std::vector<double> re_;
	std::vector<double> im_;
	std::vector<double> wt_;
	std::vector<double> stokes_;
	std::vector<double> freq_;
	std::vector<double> if_;
	std::vector<double> ra_;
	std::vector<double> dec_;

	friend std::ostream& operator<<(std::ostream& os, Vis& vis);
      };

      /**
       * Constructor.
       */
      FitsUvfReader();
      FitsUvfReader(std::string file);

      /**
       * Destructor.
       */
      virtual ~FitsUvfReader();

      void openFile(std::string file);

      std::vector<double> getData(unsigned groupPar);

      unsigned nGroup() {
	return (unsigned)nGroup_;
      }

      unsigned nPar() {
	return (unsigned)nPar_;
      }

      void initAxes();
      void getAxes();
      void initPars();
      void getPars();

      void readData(long group, Vis& vis);
      void readGroup(long group, Vis& vis);
      void readPrimaryData(long group, Vis& vis);

      virtual void getHeaderInfo();

      std::vector<double> freqs();
      void setFreqs(double* vals);
      double freqDelta();

      std::vector<double> ifs();
      void setIfs(double* vals);
      double ifDelta();

      void setStokes(double* vals);

      std::string axisName(unsigned iAxis) {
	return axes_[iAxis].type_;
      }

      std::string axisComment(unsigned iAxis) {
	return axes_[iAxis].typeComment_;
      }

      unsigned axisSize(unsigned iAxis) {
	return axes_[iAxis].n_;
      }

      // Return the group axis dimensions

      std::vector<unsigned> groupAxisDimsMatlabOrder();

      void getAxisData(double* vals, unsigned iAxis);

      unsigned nRa() {
	return ra_ ? ra_->n_ : 0;
      }

      unsigned nDec() {
	return dec_ ? dec_->n_ : 0;
      }

      unsigned nFreq() {
	return freq_ ? freq_->n_ : 0;
      }

      unsigned nIf() {
	return if_ ? if_->n_ : 0;
      }

      unsigned nStokes() {
	return stokes_ ? stokes_->n_ : 0;
      }

      unsigned nAxis() {
	return nAxis_;
      }

      std::string raComment() {
	return ra_ ? ra_->typeComment_ : "unknown";
      }

      std::string decComment() {
	return dec_ ? dec_->typeComment_ : "unknown";
      }

      std::string freqComment() {
	return freq_ ? freq_->typeComment_ : "unknown";
      }

      std::string stokesComment() {
	return stokes_ ? stokes_->typeComment_ : "unknown";
      }

    private:

      enum UVWUnit {
	SECONDS,
	CM
      };

      struct Axis {
	long n_;
	std::string type_;
	std::string typeComment_;
	double refVal_;
	double refPix_;
	double delta_;
	unsigned iCol_;
	bool isPresent_;
      };

      struct Par {
	std::string type_;
	unsigned iPar_;
	double scale_;
	double zero_;
      };

      long nAxis_;

      std::vector<Axis> axes_;

      Axis ifDefault_;

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
      Par* source_;

      UVWUnit uUnit_;
      UVWUnit vUnit_;
      UVWUnit wUnit_;

      int iu_;
      int iv_;
      int iw_;
      int ibaseline_;
      int idate_;
      int isource_;

      int nKeyword_;
      long int nGroup_;
      long int nPar_;

    }; // End class FitsUvfReader

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_FITSUVFREADER_H
