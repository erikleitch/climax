// $Id: $

#ifndef GCP_UTIL_PTSRCFINDER_H
#define GCP_UTIL_PTSRCFINDER_H

/**
 * @file PtSrcFinder.h
 * 
 * Tagged: Fri May  2 13:22:08 PDT 2014
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include "gcp/util/Angle.h"
#include "gcp/util/Constants.h"
#include "gcp/util/Cosmology.h"
#include "gcp/util/Declination.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/PtSrcReader.h"
#include "gcp/util/FirstFitsReader.h"
#include "gcp/util/NvssReader.h"
#include "gcp/util/Flux.h"
#include "gcp/util/Mass.h"

namespace gcp {
  namespace util {

    class PtSrcFinder {
    public:

      /**
       * Constructor.
       */
      PtSrcFinder();

      /**
       * Destructor.
       */
      virtual ~PtSrcFinder();

      void findSources(PtSrcReader* reader, 
		       HourAngle& ra, Declination& dec, Angle& radius, 
		       Flux& fMin, Flux& fMax,
		       std::vector<HourAngle>& ras, std::vector<Declination>& decs);

      void findSources(std::string dir, std::string cat, 
		       std::string raStr, std::string decStr, std::string radStr, 
		       std::vector<HourAngle>& ras, std::vector<Declination>& decs);
      

      void findSources(std::string dir, std::string cat, 
		       HourAngle ra, Declination dec, Angle radius,
		       std::vector<HourAngle>& ras, std::vector<Declination>& decs);
      
    }; // End class PtSrcFinder

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PTSRCFINDER_H
