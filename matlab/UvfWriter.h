#ifndef GCP_MATLAB_UVFWRITER_H
#define GCP_MATLAB_UVFWRITER_H

/**
 * @file UvfWriter.h
 * 
 * Tagged: Thu May  5 08:15:06 PDT 2005
 * 
 * @author Erik Leitch
 */
#include <vector>

#include "mex.h"
#include "matrix.h"

#include "gcp/matlab/MexParser.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Geoid.h"

namespace gcp {
  namespace matlab {
    
    class UvfWriter {
    public:
      
      /**
       * Constructor.
       */
      UvfWriter(unsigned nTelescope=8);
      
      /**
       * Destructor.
       */
      virtual ~UvfWriter();
      
      void setUvfData(const mxArray* array);
      void setMiriadData(const mxArray* array);
      void setDate(const mxArray* array);
      void setCoord(const mxArray* array);
      void setRefCoord(const mxArray* array);
      void setFreq(const mxArray* array);
      void setDeltaFreq(const mxArray* array);
      void setXyz(const mxArray* array);
      void setUvw(const mxArray* array);
      void setRms(const mxArray* array);
      void setSourceName(const mxArray* array);
      void setFileName(const mxArray* array, const mxArray* openMode=0);
      
      void setNumberOfTelescopes(const mxArray* array);
      void setDeltaIfFrequencies(const mxArray* array);
      void setBaselines(const mxArray* array);
      void setCoordsToJ2000(const mxArray* array);
      void setFirstTelescopeNum(const mxArray* array);
      
      void writeUvfFile();
      void writeFakeUvfFile();
      void writeMiriadFile();
      void writeFakeMiriadFile();
      void writeMiriadMosaic();
      void writeMiriadMosaicTest();

    private:
      
      MexParser mp_;

      unsigned nFrame_;
      unsigned nBaseline_;
      unsigned nFrequency_;
      unsigned nTelescope_;
      unsigned firstTelescopeNum_;

      gcp::util::HourAngle obsRa_;
      gcp::util::Declination  obsDec_;
      gcp::util::HourAngle refRa_;
      gcp::util::Declination  refDec_;
      gcp::util::TimeVal   obsMjd_;
      std::string fileName_;
      std::string openMode_;
      std::string srcName_;
      
      double* data_;
      double* date_;
      std::vector<gcp::util::Frequency> frequencies_;
      std::vector<gcp::util::Frequency> dfreq_;
      std::vector<gcp::util::LengthTriplet> xyz_;
      double* uvw_;
      double* rms_;
      unsigned* baselines_;
      bool fixedBaselines_;
      bool coordsAreJ2000_;

    }; // End class UvfWriter
    
  } // End namespace matlab
} // End namespace gcp



#endif // End #ifndef GCP_MATLAB_UVFWRITER_H
