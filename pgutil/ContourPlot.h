// $Id: $

#ifndef GCP_UTIL_CONTOURPLOT_H
#define GCP_UTIL_CONTOURPLOT_H

/**
 * @file ContourPlot.h
 * 
 * Tagged: Wed Apr  9 13:38:47 PDT 2014
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include "gcp/pgutil/PgUtil.h"
#include "gcp/pgutil/Trans.h"
#include "gcp/pgutil/PlotBound.h"
#include "gcp/pgutil/Plot2D.h"
#include "gcp/pgutil/PlotAxisRange.h"

namespace gcp {
  namespace util {

    class ContourPlot : public Plot2D {
    public:

      /**
       * Constructor.
       */
      ContourPlot();

      ContourPlot(PgUtil* parent, 
		  int ndata, float *zdata, int nx,int ny, 
		  float xmin, float xmax, float ymin, float ymax, 
		  float *flag,float z1, float z2, int ncontour, bool wasSigned,
		  std::string xlab, std::string ylab, std::string title, std::string unit, bool expand);

      ContourPlot(Plot2D& plot2d, int ncontour, bool wasSigned);
      
      // Destructor.

      virtual ~ContourPlot();

      void display(bool erase=false);
      void setupContours();

    private:

      std::vector<float> contours_;

      int nContour_;
      bool wasSigned_;

    }; // End class ContourPlot

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_CONTOURPLOT_H
