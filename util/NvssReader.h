// $Id: NvssReader.h,v 1.7 2006/09/11 21:20:58 gcpdaq Exp $

#ifndef GCP_UTIL_NVSSREADER_H
#define GCP_UTIL_NVSSREADER_H

/**
 * @file NvssReader.h
 * 
 * Tagged: Tue Aug 15 18:00:57 PDT 2006
 * 
 * @version: $Revision: 1.7 $, $Date: 2006/09/11 21:20:58 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/PtSrcFitsReader.h"

#include "fitsio.h"

namespace gcp {
  namespace util {

    class NvssReader : public PtSrcFitsReader {
    public:

      /**
       * Constructor.
       */
      NvssReader(std::string catalogFile);
      NvssReader();

      /**
       * Destructor.
       */
      virtual ~NvssReader();

    private:

      static Angle vlaLat_;
      static Angle raBias_;
      static Angle decBias_;
      static Flux dBiasAv_;
      static Flux dBiasErr_;
      static Flux dncBiasAv_;
      static Flux dncBiasErr_;
      static Angle calRaErr_;
      static Angle calDecErr_;
      static double calAmpErr_;

    public:

      //------------------------------------------------------------
      // Beam functions
      //------------------------------------------------------------

      // NVSS point-spread function calculation

      void setNvssPsFn(PtSrcReader::Source& src);

      // Set the error in the fitted major axis size
      
      void setMajAxisError(PtSrcReader::Source& src);

      // Set the error in the fitted minor axis size

      void setMinAxisError(PtSrcReader::Source& src);
      
      // Set the error in the fitted minor axis size

      void setPositionAngleError(PtSrcReader::Source& src);

      //------------------------------------------------------------
      // Positional routines
      //------------------------------------------------------------

      // Correct fitted positions for NVSS position bias

      void correctForPositionBias(PtSrcReader::Source& src);
      
      // Set the error in the fitted RA

      void setPositionErrors(PtSrcReader::Source& src);

      // Generic function, not just specific to NVSS

      static void setGenericPositionErrors(PtSrcReader::Source& src, Angle& calRaErr, Angle& calDecErr);

      //------------------------------------------------------------
      // Peak flux routines
      //------------------------------------------------------------

      // Correct a fitted flux for confusion bias

      void correctForConfusionBias(PtSrcReader::Source& src);

      // Set the error in the fitted peak flux value
      
      void setPeakError(PtSrcReader::Source& src);


      //------------------------------------------------------------
      // Integrated flux routines
      //------------------------------------------------------------

      // Calculate the integrated flux and error

      void setIntegratedFlux(PtSrcReader::Source& src);
      
      // Set the integrated flux and error for an unresolved source
      
      void setUnresolvedIntegratedFlux(PtSrcReader::Source& src);

      // Set the error in the integrated flux for an unresolved source
      // (see catalog.ps, section 6.6)

      void setUnresolvedIntegratedFluxError(PtSrcReader::Source& src);

      
      // Set the integrated flux and error for an unresolved source

      void setResolvedIntegratedFlux(PtSrcReader::Source& src);

      // Set the error in the integrated flux for an unresolved source

      void setResolvedIntegratedFluxError(PtSrcReader::Source& src);

      // Set the integrated flux and error for a partially resolved
      // source

      void setPartiallyResolvedIntegratedFlux(PtSrcReader::Source& src);

      // Set the error in the integrated flux for a partially-resolved source
      
      void setPartiallyResolvedIntegratedFluxError(PtSrcReader::Source& src, double pseudoPeakJy);
      
      //------------------------------------------------------------
      // Utility functions
      //------------------------------------------------------------

      // Calculate the effective SNR in a given parameter

      static double effSnr2(PtSrcReader::Source& src,
			    double alphaMajor, double alphaMinor);

      // Return the appropriate bias error for D or DnC configuration

      double getBiasErrInJy(PtSrcReader::Source& src);

      // NVSS main correction functions (translated from NVSSlist.f)

      void applyCorrections(PtSrcReader::Source& src);

      // Deconvolve the calculated point-source response from the
      // fitted parameters

      void deconvolve(PtSrcReader::Source& src);

      //------------------------------------------------------------
      // Overloaded base-class methods

      PtSrcReader::Source parseData();

      void readFitsData(long startRow, long nElement);

    }; // End class NvssReader

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_NVSSREADER_H
