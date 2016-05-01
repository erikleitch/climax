// $Id: $

#ifndef GCP_UTIL_PLOT2D_H
#define GCP_UTIL_PLOT2D_H

/**
 * @file Plot2D.h
 * 
 * Tagged: Wed Apr  9 13:38:47 PDT 2014
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include "gcp/pgutil/PgUtil.h"
#include "gcp/pgutil/Trans.h"
#include "gcp/pgutil/Plot.h"
#include "gcp/pgutil/PlotBound.h"
#include "gcp/pgutil/PlotAxisRange.h"

namespace gcp {
  namespace util {

    class Plot2D : public Plot {
    public:

      // Constructor.

      Plot2D();

      Plot2D(PgUtil* parent, 
	     int ndata, float *zdata, int nx,int ny, 
	     float xmin, float xmax, float ymin, float ymax, 
	     float *flag,float z1, float z2, 
	     std::string xlab, std::string ylab, std::string title, std::string unit, bool expand=true);

      Plot2D(Plot2D& plot2d);
      
      // Destructor.

      virtual ~Plot2D();

      void plot();
      virtual void display(bool erase=false);
      void printStats();
      void getInfo();
      void setupColormap();
      void resetContrast();
      void linePlot();
      void multiLinePlot();

      void getLine(std::vector<float>& xdata, std::vector<float>& ydata, Point& pt1, Point& pt2, bool& first);

    protected:

      virtual void initialize(bool expand);

      Cmap* grey_;
      Cmap* rain_;
      Cmap* heat_;
      Cmap* cmap_;

      unsigned cursor_;

      std::vector<float> zeros_;

      bool read_;

      Trans trans_;

      std::string unit_;

      float bright_;
      float contrast_;

      float tr_[6];
      int i1_;
      int i2_;
      int j1_;
      int j2_;

    }; // End class Plot2D

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PLOT2D_H
