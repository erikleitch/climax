// $Id: $

#ifndef GCP_UTIL_GREYSCALEPLOT_H
#define GCP_UTIL_GREYSCALEPLOT_H

/**
 * @file GreyscalePlot.h
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

    class GreyscalePlot : public Plot2D {
    public:

      /**
       * Constructor.
       */
      GreyscalePlot(PgUtil* parent, 
		    int ndata, float *zdata, int nx,int ny, 
		    float xmin, float xmax, float ymin, float ymax, 
		    float *flag,float z1, float z2, 
		    std::string xlab, std::string ylab, std::string title, std::string unit);
      /**
       * Destructor.
       */
      virtual ~GreyscalePlot();

      void display(bool erase=false);

    }; // End class GreyscalePlot

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_GREYSCALEPLOT_H
