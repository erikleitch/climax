// $Id: $

#ifndef GCP_UTIL_LINEPLOT_H
#define GCP_UTIL_LINEPLOT_H

/**
 * @file LinePlot.h
 * 
 * Tagged: Thu Apr 10 09:46:48 PDT 2014
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include "gcp/pgutil/Plot.h"

#include <string>
#include <vector>

namespace gcp {
  namespace util {

    class PgUtil;

    class LinePlot : public Plot {
    public:

      // Constructor.
      
      LinePlot(PgUtil* parent,
	       int narr, float* xarr, float* yarr, float* earr,
	       std::string xlab, std::string ylab, std::string title, bool doLine, bool doErr);

      // Destructor.

      virtual ~LinePlot();

      void initialize();
      void plot();
      void display(bool erase=false);
      void printValNearestToPoint(float x, float y);
      void cacheHighestPoint();
      void cacheNHighestPeaks();
      void cachePoint(float x1, float x2, float y1, float y2, bool doLine);
      void labelPoints(bool erase);
      void clearLabeledPoints();

    private:
      
      std::vector<float> xs_;
      std::vector<float> ys_;

      std::vector<float> yDataErrHi_;
      std::vector<float> yDataErrLo_;

      float* xdata_;
      float* edata_;

      bool doLine_;
      bool doErr_;
      
      std::vector<Plot::LabeledPoint> labeledPoints_;
     
    }; // End class LinePlot

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_LINEPLOT_H
