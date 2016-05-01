// $Id: $

#ifndef GCP_UTIL_PLOTAXISRANGE_H
#define GCP_UTIL_PLOTAXISRANGE_H

/**
 * @file PlotAxisRange.h
 * 
 * Tagged: Wed Apr  9 13:50:08 PDT 2014
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include <iostream>

namespace gcp {
  namespace util {

    class PlotAxisRange {
    public:

      /**
       * Constructor.
       */
      PlotAxisRange();

      /**
       * Destructor.
       */
      virtual ~PlotAxisRange();

      // Assignment specifying the sense of the axes

      void setTo(double min, double max, bool reverse);
      void setTo(double min, double max, unsigned ndata, float* data, bool reverse);

      // Assignment preserving the current sense of this axis

      void setTo(double min, double max);
      void setTo(double min, double max, unsigned ndata, float* data);

      void expandPerc(double perc);
      void expandAbs(double val);

      double midPoint();
      double range();

      double absMin();
      double absMax();

      double plotMin();
      double plotMax();

      double dataMin();
      double dataMax();

      bool contains(double val);

      friend std::ostream& operator<<(std::ostream& os, PlotAxisRange& par);

    public:

      void setMinMax(double min, double max);

      double min_;
      double max_;
      double dataMin_;
      double dataMax_;
      bool reverse_;
      bool isSpecified_;

    }; // End class PlotAxisRange

    std::ostream& operator<<(std::ostream& os, PlotAxisRange& par);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PLOTAXISRANGE_H
