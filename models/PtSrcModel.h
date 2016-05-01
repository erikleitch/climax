// $Id: PtSrcModel.h,v 1.5 2012/05/31 16:22:56 eml Exp $

#ifndef GCP_MODELS_PTSRCMODEL_H
#define GCP_MODELS_PTSRCMODEL_H

/**
 * @file PtSrcModel.h
 * 
 * Tagged: Wed Aug 18 16:57:59 PDT 2010
 * 
 * @version: $Revision: 1.5 $, $Date: 2012/05/31 16:22:56 $
 * 
 * @author tcsh: Erik Leitch
 */
#include "gcp/fftutil/Generic2DAngularModel.h"

#include "gcp/util/Variate.h"

namespace gcp {

  namespace util {
    class UvDataGridder;
  }

  namespace models {

    class PtSrcModel : public gcp::util::Generic2DAngularModel {
    public:

      struct UvParams {
	gcp::util::Image*     beam_;
	gcp::util::Frequency* freq_;
      };

      /**
       * Constructor.
       */
      PtSrcModel();

      /**
       * Destructor.
       */
      virtual ~PtSrcModel();

      void fillSzImage(gcp::util::Image& image, gcp::util::Frequency& frequency);
      void fillXrayImage(gcp::util::Image& image, gcp::util::Frequency& frequency);
      void fillImage(unsigned type, gcp::util::Image& image, void* params);

      void fillUvData(unsigned type, gcp::util::UvDataGridder& gridder, void* params=0);

      // The flux of the source

      void setFlux(gcp::util::Flux& flux);
      void setJy(double fluxJy);

      // The frequency of the source

      void setFrequency(gcp::util::Frequency& freq);
      void setGHz(double freqGHz);

      // The spectral index of this source

      void setSpectralIndex(double spectralIndex);

      // An offset for this source

      void setOffset(gcp::util::Angle& xOff, gcp::util::Angle& yOff);

      void addComponents();

      void debugPrint();

      static double cos(double argRadians);
      static double sin(double argRadians);
      static void sincos(double argRad, double& sVal, double& cVal);

    private:

      void initializeLookupTable();
      void initializeLookupTable(double precision);

      // The flux of this source

      gcp::util::Flux flux_;

      // The frequency at which the source flux was specified

      gcp::util::Frequency freq_;

      // The spectral index of this source

      gcp::util::Variate spectralIndex_;

      // A lookup table for sin/cos

      static std::valarray<double> sinLookup_;
      static double lookupResRadians_;
      static unsigned lookupNpt_;

    }; // End class PtSrcModel

  } // End namespace models
} // End namespace gcp



#endif // End #ifndef GCP_MODELS_PTSRCMODEL_H
