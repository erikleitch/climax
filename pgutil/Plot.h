// $Id: $

#ifndef GCP_UTIL_PLOT_H
#define GCP_UTIL_PLOT_H

/**
 * @file Plot.h
 * 
 * Tagged: Thu Apr 10 09:50:30 PDT 2014
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
#include "gcp/pgutil/PlotAxisRange.h"
#include "gcp/pgutil/PlotBound.h"
#include "cpgplot.h"

#include <string>
#include <vector>

namespace gcp {
  namespace util {

    class PgUtil;

    class Plot {
    public:

      //------------------------------------------------------------
      // A point
      //------------------------------------------------------------

      struct Point {
	double x_;
	double y_;
      };

      struct LabeledPoint {
	PlotBound bound_;
	std::string label_;
	bool doLine_;
      };

      //------------------------------------------------------------
      // Type of plot
      //------------------------------------------------------------

      enum PlotType {
	PLOT_PLOT1D,
	PLOT_PLOT2D
      };

      //------------------------------------------------------------
      // Function keys in the graphical interface
      //------------------------------------------------------------

      enum {                     
	G_CUR   = 'A',
	G_FLG   = 'B',           
	G_CUT   = 'C',
	G_CONT  = 'C',
	G_CAN   = 'D',
	G_DEF   = 'D',

	G_FID   = 'F',
	G_GREY  = 'G',
	G_HELP  = 'H',
	G_INS   = 'I',
	G_FIT   = 'J',
	G_RAD   = 'K',
	G_DIS   = 'L',

	G_OVRPLT= 'O',
	G_MODEL = 'M',
	G_NXT   = 'N',
	G_HEAT  = 'Q',
	G_RAIN  = 'R',
	G_STAT  = 'S',

	G_BEAM  = 'T',   

	G_HORI  = 'U',
	G_VERT  = 'V',

	G_QUIT  = 'X',
	G_YSC   = 'Y',
	G_ZOOM  = 'Z',
	G_CROS  = '+',
      };	    

      //------------------------------------------------------------
      // Types of zoom
      //------------------------------------------------------------

      enum ZoomType {
	ZOOM_UNKNOWN = 0x0,
	ZOOM_X = 0x1,
	ZOOM_Y = 0x2,
	ZOOM_BOTH = ZOOM_X | ZOOM_Y,
	ZOOM_LINE = 0x4
      };

      /**
       * Constructor.
       */
      Plot();
      Plot(PgUtil* parent, unsigned ndata, std::string xlab, std::string ylab, std::string title);

      /**
       * Destructor.
       */
      virtual ~Plot();

      void initialize();
      void zoom(char keyFull, ZoomType type);
      void getRange(char keyFull, ZoomType type, PlotBound& bound, Point& pt1, Point& pt2);

      void getPoint(ZoomType type, Point& pt1);
      void getLine(ZoomType type, Point& pt1, Point& pt2);

      void getEndOfRange(char keyFull, ZoomType type, PlotBound& bound, Point& pt1, Point& pt2);
      void setupPlotBoundaries();

      virtual void display(bool erase=false) = 0;
      void getNumber(int& num, bool& wasSigned);
      void computeGreyscaleRange();
      std::vector<float> getCurrentlyDisplayedData(unsigned& nxSub, unsigned& nySub);

    protected:

      float xpos_[2],ypos_[2];

      // Store the last position of the cursor when the user entered a key

      Point anchor_;

      PlotType plotType_;
      PgUtil* parent_;
      std::string xlab_;
      std::string ylab_;
      std::string title_;

      float xvp2_;
      float xvp1_;
      float yvp2_;
      float yvp1_;
      float xwin2_;
      float xwin1_;
      float ywin2_;
      float ywin1_;

      unsigned ndata_;
      float* data_;

      PlotBound bound_;
      PlotBound boundCurr_;

      bool docurs_;
      bool wasopen_;

      bool dofull_;
      bool cancel_;
      char key_;
      bool autoscale_;

      double dx_;
      double dy_;

      unsigned nx_;
      unsigned ny_;

      PlotAxisRange zrng_;
      PlotAxisRange zrngCurr_;

    }; // End class Plot

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PLOT_H
