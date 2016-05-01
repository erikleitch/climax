// $Id: PgUtil.h,v 1.6 2012/05/16 18:00:50 eml Exp $

#ifndef GCP_UTIL_PGUTIL_H
#define GCP_UTIL_PGUTIL_H

/**
 * @file PgUtil.h
 * 
 * Tagged: Fri Aug 22 11:03:53 PDT 2008
 * 
 * @version: $Revision: 1.6 $, $Date: 2012/05/16 18:00:50 $
 * 
 * @author Erik Leitch.
 */
#include "gcp/pgutil/color_tab.h"
#include "gcp/pgutil/PgModelManager.h"

#include <string>
#include <vector>
#include <valarray>

#define PGUTIL_UNIT_CALLBACK(fn)  void (fn)(double x, double y, std::string& xstr, std::string& ystr, void* args)
#define PGUTIL_COORD_CALLBACK(fn) void (fn)(double x, double y, std::string& xstr, std::string& ystr, void* args)
#define PGUTIL_PRIOR_CALLBACK(fn) std::string (fn)(double x, double y, void* args)

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

namespace gcp {
  namespace util {

    enum CursorMode {
      B_NORM=0,
      B_LINE=1,
      B_RECT=2,
      B_YRNG=3,
      B_XRNG=4,
      B_YVAL=5,
      B_XVAL=6,
      B_CROSS=7
    };

    class PgUtil {
    public:

      enum Stat {
	STAT_NONE,
	STAT_MAXL,
	STAT_UPPER,
	STAT_LOWER,
      };

      enum Just {
	JUST_LEFT,
	JUST_RIGHT,
	JUST_CENTER,
      };

      struct Header {
	std::string text_;
	Just just_;
      };

      /**
       * Constructor.
       */
      PgUtil();

      /**
       * Destructor.
       */
      virtual ~PgUtil();

      enum {
	BLACK        = 0,
	WHITE        = 1,
	RED          = 2,
	GREEN        = 3,
	DEEP_BLUE    = 4,
	PALE_BLUE    = 5,
	MAGENTA      = 6,
	YELLOW       = 7,
	ORANGE       = 8,
	YELLOW_GREEN = 9,
	SAGE_GREEN   = 10,
	SLATE_BLUE   = 11,
	PURPLE       = 12,
	PINK         = 13,
	DARK_GRAY    = 14,
	DARK_GREY    = 14,
	LIGHT_GRAY   = 15,
	LIGHT_GREY   = 15,
      };

      static int open(std::string device);
      static void close();
      static void subplot(int nx, int ny);
      static void advance();

      static void setInteractive(bool inter);

      static void contour(std::vector<double>& zdata, int nx, int ny, 
			  double xmina=0, double xmaxa=1, double ymina=0, double ymaxa=1, 
			  double *flag=0,double z1=0, double z2=0, unsigned ncontour=10,
			  std::string xlab=std::string(""), std::string ylab=std::string(""), std::string title=std::string(""), std::string unit=std::string(""));
      
      static void contour(int ndata, double *zdata, int nx, int ny, 
			  double xmina=0, double xmaxa=1, double ymina=0, double ymaxa=1, 
			  double *flag=0,double z1=0, double z2=0, unsigned ncontour=10,
			  std::string xlab=std::string(""), std::string ylab=std::string(""), std::string title=std::string(""), std::string unit=std::string(""));
      
      static void contour(int ndata, float *zdata, int nx, int ny, 
			  float xmina=0, float xmaxa=1, float ymina=0, float ymaxa=1, 
			  float *flag=0,float z1=0, float z2=0, unsigned ncontour=10, 
			  std::string xlab=std::string(""), std::string ylab=std::string(""), std::string title=std::string(""), std::string unit=std::string(""));
      
      static void greyScale(std::vector<double>& zdata, int nx, int ny, 
			    double xmina=0, double xmaxa=1, double ymina=0, double ymaxa=1, 
			    double *flag=0,double z1=0, double z2=0, 
			    std::string xlab=std::string(""), std::string ylab=std::string(""), std::string title=std::string(""), std::string unit=std::string(""));
      
      static void greyScale(int ndata, double *zdata, int nx, int ny, 
			    double xmina=0, double xmaxa=1, double ymina=0, double ymaxa=1, 
			    double *flag=0,double z1=0, double z2=0, 
			    std::string xlab=std::string(""), std::string ylab=std::string(""), std::string title=std::string(""), std::string unit=std::string(""));
      
      static void greyScale(int ndata, float *zdata, int nx, int ny, 
			    float xmina=0, float xmaxa=1, float ymina=0, float ymaxa=1, 
			    float *flag=0,float z1=0, float z2=0, 
			    std::string xlab=std::string(""), std::string ylab=std::string(""), std::string title=std::string(""), std::string unit=std::string(""));
      
      static void grayScale(int ndata, float *zdata, int nx,int ny, 
			    float xmina=0, float xmaxa=1, float ymina=0, float ymaxa=1, 
			    float *flag=0,float z1=0, float z2=0, 
			    std::string xlab=std::string(""), std::string ylab=std::string(""), std::string title=std::string(""), std::string unit=std::string(""));
      
      static void errPlot(std::vector<double>& xarr, std::vector<double>& yarr, std::vector<double>& earr,
			  std::string xlab=std::string(""), std::string ylab=std::string(""), std::string title=std::string(""));

      static void errPlot(std::valarray<double>& xarr, std::valarray<double>& yarr, std::valarray<double>& earr,
			  std::string xlab=std::string(""), std::string ylab=std::string(""), std::string title=std::string(""));

      // Plot a simple line plot.
      
      static void linePlot(std::vector<double>& xarr, std::vector<double>& yarr, 
			   std::string xlab=std::string(""), std::string ylab=std::string(""), std::string title=std::string(""),
			   bool doLine=true);

      static void linePlot(std::vector<float>& xarr, std::vector<float>& yarr, 
			   std::string xlab=std::string(""), std::string ylab=std::string(""), std::string title=std::string(""),
			   bool doLine=true);

      static void linePlot(int narr, double* xarr, double* yarr, double* earr,
			   std::string xlab=std::string(""), std::string ylab=std::string(""), std::string title=std::string(""),
			   bool doLine=true, bool doErr=false);

      static void linePlot(int narr, float* xarr, float* yarr, float* earr, 
			   std::string xlab=std::string(""), std::string ylab=std::string(""), std::string title=std::string(""),
			   bool doLine=true, bool doErr=false);

      static void linePlot(std::vector<double>& yarr, bool doLine=true);
 
      static void linePlot(std::vector<float>& yarr, bool doLine=true);
      static void linePlot(std::vector<unsigned>& yarr, bool doLine=true);

      static void plotPowerSpectrum(std::vector<float>& yarr, bool doLine=true);
      static void plotPowerSpectrum(std::vector<double>& yarr, bool doLine=true);

      static void setPrompt(bool prompt);

      static void useHeader(bool use) {
	useHeader_ = use;
      }

      static void setHeader(std::string text, Just just) {
	header_.text_ = text;
	header_.just_ = just;
      }

      static void setNcontour(unsigned nCont) {
	nContour_ = nCont;
      }

      static void setReverseX(bool reverseX) {
	reverseX_ = reverseX;
      }

      static void setReverseY(bool reverseY) {
	reverseY_ = reverseY;
      }

      static void setOverplot(bool overplot) {
	overplot_ = overplot;
      }

      static void setLogPlot(bool logPlot) {
	logPlot_ = logPlot;
      }

      static void setUnitCallback(PGUTIL_UNIT_CALLBACK(*fn), void* args)
      {
	unitCallback_ = fn;
	unitCallbackArgs_ = args;
      }

      static void setCoordCallback(PGUTIL_COORD_CALLBACK(*fn), void* args)
      {
	coordCallback_ = fn;
	coordCallbackArgs_ = args;
      }

      static void setPriorCallback(PGUTIL_PRIOR_CALLBACK(*fn), void* args)
      {
	priorCallback_ = fn;
	priorCallbackArgs_ = args;
      }

      static void setWnad(bool wnad) {
	wnad_ = wnad;
      }

      static void setWin(bool win) {
	win_ = win;
      }

      static void setBox(bool box) {
	box_ = box;
      }

      static void setVp(bool vp) {
	vp_ = vp;
      }

      static void setTick(bool tick) {
	tick_ = tick;
      }

      static void setXTick(bool xTick) {
	xTick_ = xTick;
      }

      static void setYTick(bool yTick) {
	yTick_ = yTick;
      }

      static void setXTickLabelAtBottom(bool bottom) {
	xTickLabelAtBottom_ = bottom;
      }

      static void setYTickLabelAtLeft(bool left) {
	yTickLabelAtLeft_ = left;
      }

      static void setXTickLabeling(bool xTickLabeling) {
	xTickLabeling_ = xTickLabeling;
      }

      static void setYTickLabeling(bool yTickLabeling) {
	yTickLabeling_ = yTickLabeling;
      }

      static void setLabel(bool label) {
	label_ = label;
      }

      static void setTitle(bool title) {
	title_ = title;
      }

      static void setXLabel(bool xLabel) {
	xLabel_ = xLabel;
      }

      static void setYLabel(bool yLabel) {
	yLabel_ = yLabel;
      }

      static void setXLabelString(std::string xLabel) {
	xLabelString_ = xLabel;
      }

      static void setYLabelString(std::string yLabel) {
	yLabelString_ = yLabel;
      }

      static void setWedge(bool wedge) {
	wedge_ = wedge;
      }

      static void setUsedefs(bool usedefs) {
	usedefs_ = usedefs;
      }

      static void setXmin(float xmin) {
	xmin_ = xmin;
      }

      static void setXmax(float xmax) {
	xmax_ = xmax;
      }

      static void setYmin(float ymin) {
	ymin_ = ymin;
      }

      static void setYmax(float ymax) {
	ymax_ = ymax;
      }

      static void setColormap(std::string cmap);

      static void setPostscriptFontName(std::string fontName);

      static void setCharacterHeight(double ch) {
	ch_ = ch;
      }

      static void queryDevice(bool& wasopen);
      static bool haveCursor();
      static bool haveOpenDevice();

      static void wnad(float xmin, float xmax, float ymin, float ymax);

      static void useXNlab(bool xnlab) {
	xnlab_ = xnlab;
      }

      static void drawMean(bool mean) {
	drawMean_ = mean;
      }

      static void draw1SigmaConfidenceInterval(bool draw) {
	draw1SigmaConfidenceInterval_ = draw;
      }

      static void setStat(PgUtil::Stat stat) {
	stat_ = stat;
      }

      static void setNsigma(double nSigma) {
	nSigma_ = nSigma;
      }

      static void setBoxColor(int ci) {
	boxCi_ = ci;
	useBoxCi_ = true;
      }

      static void clearBoxColor() {
	useBoxCi_ = false;
      }

      static void setTraceColor(int ci) {
	traceCi_ = ci;
	useTraceCi_ = true;
      }

      static void clearTraceColor() {
	useTraceCi_ = false;
      }

      static void setZmin(float zmin) {
	zmin_ = zmin;
      }

      static void setZmax(float zmax) {
	zmax_ = zmax;
      }

      static PGUTIL_UNIT_CALLBACK(*unitCallback_);
      static void* unitCallbackArgs_;
      static PGUTIL_COORD_CALLBACK(*coordCallback_);
      static void* coordCallbackArgs_;
      static PGUTIL_PRIOR_CALLBACK(*priorCallback_);
      static void* priorCallbackArgs_;
      static PgModelManager pgManager_;

      static void clearPgManager();
      static void insertPgManager(PgModelManager& pm);

      static void addModel(PgModel& model);

    public:

      static Stat stat_;
      static double nSigma_;
      static bool useBoxCi_;
      static unsigned boxCi_;
      static bool useTraceCi_;
      static unsigned traceCi_;
      static bool draw1SigmaConfidenceInterval_;
      static bool drawMean_;
      static bool xnlab_;
      static bool logPlot_;
      static double ch_;
      static bool interactive_;
      static Cmap* cmap_;
      static bool overplot_;
      static bool reverseX_;
      static bool reverseY_;
      static unsigned nContour_;
      static bool vp_;
      static bool win_;
      static bool box_;
      static bool wnad_;
      static bool tick_;
      static bool xTick_;
      static bool yTick_;
      static bool xTickLabeling_;
      static bool yTickLabeling_;
      static bool xTickLabelAtBottom_;
      static bool yTickLabelAtLeft_;
      static bool label_;
      static bool xLabel_;
      static bool yLabel_;
      static bool title_;
      static bool wedge_;
      static bool usedefs_;
      static float xmin_, xmax_;
      static float ymin_, ymax_;
      static float zmin_, zmax_;
      static std::string xLabelString_;
      static std::string yLabelString_;
      static Header header_;
      static bool useHeader_;

      // Make a contour plot

      static int v_contour(int ndata, float *zdata, int nx,int ny, 
			   float xmina, float xmaxa, float ymina, float ymaxa, 
			   float *flag, float z1, float z2, unsigned ncontour,
			   char *xlab, char *ylab, char *title, char *unit);

      // Make a gray scale plot

      int v_grey2(int ndata, float *zdata, int nx,int ny, 
		  float xmina=0, float xmaxa=1, float ymina=0, float ymaxa=1, 
		  float* flag=0,float z1=0, float z2=0, 
		  std::string xlab=std::string(""),  std::string ylab=std::string(""),
                  std::string title=std::string(""), std::string unit=std::string(""));

      // Plot a simple line plot.
      
      static int v_lplot(int narr, float* xarr, float* yarr, float* earr,
			 std::string xlab=std::string(""), std::string ylab=std::string(""), std::string title=std::string(""),
			 bool doLine=true, bool doErr=false);

      static int v_radplot(float data[],int nbin, float rmin, float rmax, 
			   float xmin, float xmax, float ymin, float ymax, int nx, int ny);


      // Module for v_lplot.

      static int v_lnewx(float xmins, float xmaxs, float ymins, float ymaxs, 
			 float *xmin, float *xmax, float *ymin, float *ymax, 
			 int narr,float xarr[], float yarr[]);
	
      // Module for v_lplot.

      static int v_lwrite(int narr, float xarr[], float yarr[], float xmin,
			  float xmax, float ymin, float ymax);

      // Module for v_lplot.

      static int v_lnum(int narr, float xarr[], float yarr[], float xmin,
			float xmax, float ymin, float ymax);

      // Module for v_lplot.
      
      static int v_lten(int narr, float xarr[], float yarr[], float xmin,
			float xmax, float ymin, float ymax);

      // Module for v_hist.

      static int v_lzoom(float xmins, float xmaxs, float ymins, float ymaxs, 
			 float *xmin, float *xmax, float *ymin, float *ymax);

      static void draw1DHistogramMean(float* hist, unsigned nbin, float min, float max);
      static void draw1DHistogramConfidenceIntervalAboutMean(float* hist, unsigned nbin, float min, float max);
      static void draw1DHistogramMax(float* hist, unsigned nbin, float min, float max);
      static void drawRefMarker(float* hist, unsigned nbin, float min, float max, double refXVal);
      static void draw1DHistogramConfidenceIntervalAboutMax(float* hist, unsigned nbin, float min, float max);

      static double drawConfidenceInterval(float* hist, unsigned nbin, float min, float max);

      // Module for v_lplot -- Draw the line plot.
      
      static int v_ldraw(int narr, float xarr[], float yarr[], float earr[], 
			 std::string xlab, std::string ylab, std::string title,
			 bool doLine, bool doErr,
			 double xmin, double xmax, double ymin, double ymax);
      
      static void drawHeader();

      void drawBox(unsigned ci=1, bool override=false);
      void setupPlotBoundaries(double xmin, double xmax, double ymin, double ymax);
      void drawWedge(double zmin, double zmax, std::string unit);
      void drawLabels(std::string xlab, std::string ylab, std::string title);

    public:

      static void histogram(std::vector<double>* dptr,   unsigned nbin, std::vector<unsigned>* multiplicity=0);
      static void histogram(unsigned ndata, float* data, unsigned nbin, std::vector<unsigned>* multiplicity=0);

      static void binPlot(std::vector<double>& x, std::vector<double>& y, unsigned nbin);

      static void binData(std::vector<double>& x, std::vector<double>& y, unsigned nbin,
			  std::vector<double>& xret, std::vector<double>& yret);

      static void histogram2D(unsigned ndata, float* data1, float* data2, unsigned nbin1, unsigned nbin2, std::vector<unsigned>* multiplicity=0);
      static void histogram2D(std::vector<double>* dptr1, std::vector<double>* dptr2, unsigned nbin, std::vector<unsigned>* multiplicity=0);

      static int v_hist(float* hist, float xmins, float xmaxs, int nbin, std::string title);

      // Draw the histogram box and envelope

      static int v_hdraw(float* hist, int nbin, float xmins, float xmaxs, 
			 float ymins, float ymaxs, float xmin, float xmax, 
			 float ymin,float dx, std::string ylab, std::string title);

      static int v_hzoom(float xmins, float xmaxs, float ymins, float ymaxs, 
			 float dx);

      static int v_hnewx(float xmins, float xmaxs, float ymins, float ymaxs, 
			 float dx);

      static void substituteForPgplot(std::string& str);

      static void label(std::string xlabel, std::string ylabel, std::string title);

      static std::vector<std::string> fontNames_;

      static std::vector<std::string> getFontNames();

      void displayColors();

      static void plotPoint(double x, double y, int marker);

      static int  indexx(int npts, float arrin[], int **indx);

      static void text(std::string text, double x, double y, Angle angle=Angle(Angle::Degrees(), 0.0), Just just=JUST_LEFT);

      static void initialize();

    private:

      static void insert_in_heap(float arrin[], int indx[], int node, int num_node, int new_el);
      static float pgplotJust(Just just);

    }; // End class PgUtil

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PGUTIL_H
