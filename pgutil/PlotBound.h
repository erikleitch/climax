// $Id: $

#ifndef GCP_UTIL_PLOTBOUND_H
#define GCP_UTIL_PLOTBOUND_H

/**
 * @file PlotBound.h
 * 
 * Tagged: Wed Apr  9 13:45:54 PDT 2014
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include "gcp/pgutil/PlotAxisRange.h"

#include <iostream>

namespace gcp {
  namespace util {

    class PlotBound {
    public:

      /**
       * Constructor.
       */
      PlotBound();

      /**
       * Destructor.
       */
      virtual ~PlotBound();

      PlotAxisRange xrng_;
      PlotAxisRange yrng_;

      void setTo(double xmin, double xmax, double ymin, double ymax, 
		 bool reverseX=false, bool reverseY=false);
      void expandPerc(double perc);
      void expandAbs(double delta);

      // True if this bound contains the requested point

      bool contains(double x, double y);
      bool containsX(double x);
      bool containsY(double y);
      double maxRadialExtent();
 
      double absXmin();
      double absXmax();
      double absYmin();
      double absYmax();

      double plotXmin();
      double plotXmax();
      double plotYmin();
      double plotYmax();

      double dataXmin();
      double dataXmax();
      double dataYmin();
      double dataYmax();

      double xRange();
      double yRange();

      double range();

      friend std::ostream& operator<<(std::ostream& os, PlotBound& bound);
      
    }; // End class PlotBound

    std::ostream& operator<<(std::ostream& os, PlotBound& bound);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PLOTBOUND_H
