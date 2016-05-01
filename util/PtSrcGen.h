// $Id: PtSrcGen.h,v 1.1.1.1 2010/07/13 17:56:53 eml Exp $

#ifndef GCP_UTIL_PTSRCGEN_H
#define GCP_UTIL_PTSRCGEN_H

/**
 * @file PtSrcGen.h
 * 
 * Tagged: Fri Mar 27 13:26:18 PDT 2009
 * 
 * @version: $Revision: 1.1.1.1 $, $Date: 2010/07/13 17:56:53 $
 * 
 * @author Erik Leitch.
 */
#include "gcp/util/Angle.h"
#include "gcp/util/Flux.h"
#include "gcp/util/SolidAngle.h"

#include <vector>

namespace gcp {
  namespace util {

    class PtSrcGen {
    public:

      /**
       * Constructor.
       */
      PtSrcGen();

      /**
       * Destructor.
       */
      virtual ~PtSrcGen();

      //------------------------------------------------------------
      // Set the dN/dS to use, as a power law.
      //
      // Note I use the convention:
      //
      //      dN          -gamma 
      //   --------- = k S 
      //   dS dOmega

      void setDnDs(double k, double gamma, Flux fu, SolidAngle au);

      // Set the dN/dS to use, as a user-specified function x = flux, y = num

      void setDnDs(std::vector<double> flux, std::vector<double> num, 
		   const Flux::Jansky& fluxUnit, 
		   const SolidAngle::Steradians& angleUnit);

      unsigned getNSrc(Flux& fluxMin, SolidAngle& area, bool doRand=true);
      unsigned getNSrc(Flux& fluxMin, Flux& fluxMax, SolidAngle& area, bool doRand=true);

      // Generate a list of sources drawn from the specified
      // distribution, within the specified solid angle

      std::vector<Flux> generateSources(Flux fluxMin, SolidAngle sr, bool doRand=true);
      std::vector<Flux> generateSources(Flux fluxMin, Flux fluxMax, SolidAngle sr, bool doRand=true);

      // Generate a list of source fluxes and positions within the
      // specified x/y box
 
      void generateSources(Flux fluxMin, Angle x, Angle y, 
			   std::vector<Flux>&  srcFlux, 
			   std::vector<Angle>& srcX, 
			   std::vector<Angle>& srcY,
			   bool doRand=true);

      void generateSources(Flux fluxMin, Flux fluxMax, Angle x, Angle y, 
			   std::vector<Flux>&  srcFlux, 
			   std::vector<Angle>& srcX, 
			   std::vector<Angle>& srcY,
			   bool doRand=true);


      // Return the confusion limit for the requested solid angle

      Flux confusionLimit(SolidAngle omega);

    private:

      bool dNdSIsSet_;
      double k_;
      double gamma_;
      
      SolidAngle au_;
      Flux       fu_;

    }; // End class PtSrcGen

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PTSRCGEN_H
