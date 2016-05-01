#include "gcp/util/Angle.h"
#include "gcp/util/Exception.h"

#include "gcp/pgutil/ContourPlot.h"
#include "gcp/pgutil/GreyscalePlot.h"
#include "gcp/pgutil/LinePlot.h"
#include "gcp/pgutil/PgModelManager.h"
#include "gcp/pgutil/PgModel.h"
#include "gcp/pgutil/PgUtil.h"
#include "gcp/pgutil/Trans.h"

#include "gcp/util/Stats.h"
#include "gcp/util/String.h"

#include "gcp/fftutil/Dft1d.h"

#include <cmath>
#include <cstring>
#include <cstdio>
#include <cctype>

#include <stdlib.h>

#include "cpgplot.h"

using namespace std;

using namespace gcp::util;

#define FNINT(f) floor((f)+0.5f)
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

/*
 * Define the selection keys.
 */
enum {
  KEY_NONE='\0',  /* Null key press */
  KEY_CUR ='A',   /* Key for cursor position input */
  KEY_BIN ='B',   /* Key to introduce binning */
  KEY_CUT ='C',   /* Key to introduce cutting */
  KEY_CAN ='D',   /* Key to cancel incomplete select range */
  KEY_ERR ='E',   /* Toggle display of error bars */
  KEY_FLG ='F',   /* Toggle display of flagged data */
  KEY_LINE='G',   /* Toggle connect pts for binned plots */
  KEY_HELP='H',   /* Key to list usage information */
  KEY_INS ='I',   /* Inspect the value of the point nearest the cursor */
  KEY_FIT ='J',   /* Fit a gaussian to the selected data */
  KEY_KEEP='K',   /* Key to keep default flag displays */
  KEY_DIS ='L',   /* Key to redisplay the current plot */
  KEY_MENU='M',   /* Key to toggle menu display */
  KEY_NXT ='N',   /* Key to display the next sub-plot(s) */
  KEY_OUTL='O',   /* Perform outlier editing from the interface */
  KEY_PREV='P',   /* Key to display the previous sub-plot(s) */
  KEY_TEST='Q',   /* Test function key */
  KEY_REST='R',   /* Key to introduce boxed-restore */
  KEY_STAT='S',   /* Compute simple statistics for the selected points */
  KEY_BEAM='T',   /* Draw the beam */
  KEY_UT  ='U',   /* Key to select UT input mode */
  KEY_Y   ='V',   /* Key to select Y input mode */
  KEY_WRAP='W',   /* Key to wrap/unwrap a phase point */
  KEY_QUIT='X',   /* Key to quit from this function */
  KEY_YSC ='Y',   /* Toggle autoscale of y-axis */
  KEY_ZOOM='Z',   /* Select a zoom box */
  KEY_CROS='+',   /* Cross hair cursor */
  KEY_EDIT=' '    /* Toggle editing type */
};

//------------------------------------------------------------
// Initialize static variables
//------------------------------------------------------------

bool PgUtil::useBoxCi_ = false;
unsigned PgUtil::boxCi_ = 10;

bool PgUtil::useTraceCi_ = false;
unsigned PgUtil::traceCi_ = 10;

Cmap* PgUtil::cmap_ = &std_cmaps[4];
bool PgUtil::overplot_ = false;
bool PgUtil::reverseX_ = false;
bool PgUtil::reverseY_ = false;
unsigned PgUtil::nContour_ = 0;
bool PgUtil::vp_ = true;
bool PgUtil::win_ = true;
bool PgUtil::box_ = true;
bool PgUtil::wnad_ = false;
bool PgUtil::tick_ = true;
bool PgUtil::xTick_ = true;
bool PgUtil::yTick_ = true;
bool PgUtil::xTickLabeling_ = true;
bool PgUtil::yTickLabeling_ = true;
bool PgUtil::xTickLabelAtBottom_ = true;
bool PgUtil::yTickLabelAtLeft_ = true;
bool PgUtil::label_ = true;
bool PgUtil::xLabel_ = true;
bool PgUtil::yLabel_ = true;
bool PgUtil::title_ = true;
bool PgUtil::wedge_ = true;
bool PgUtil::usedefs_ = false;
float PgUtil::xmin_ = 0.0;
float PgUtil::xmax_ = 0.0;
float PgUtil::ymin_ = 0.0;
float PgUtil::ymax_ = 0.0;
double PgUtil::ch_ = 1.0;
bool PgUtil::interactive_ = true;
bool PgUtil::logPlot_     = false;
bool PgUtil::xnlab_       = true;
bool PgUtil::drawMean_    = true;
bool PgUtil::draw1SigmaConfidenceInterval_ = false;
float PgUtil::zmin_ = 0.0;
float PgUtil::zmax_ = 0.0;

PGUTIL_UNIT_CALLBACK(*PgUtil::unitCallback_)   = 0;
void* PgUtil::unitCallbackArgs_  = 0;
PGUTIL_COORD_CALLBACK(*PgUtil::coordCallback_) = 0;
void* PgUtil::coordCallbackArgs_ = 0;
PGUTIL_PRIOR_CALLBACK(*PgUtil::priorCallback_) = 0;
void* PgUtil::priorCallbackArgs_ = 0;

std::string PgUtil::xLabelString_ = "";
std::string PgUtil::yLabelString_ = "";

std::vector<std::string> PgUtil::fontNames_ = PgUtil::getFontNames();

PgUtil::Stat  PgUtil::stat_ = PgUtil::STAT_NONE;
double PgUtil::nSigma_ = 1.0;

PgUtil::Header  PgUtil::header_;
bool  PgUtil::useHeader_ = false;

PgModelManager PgUtil::pgManager_;

/**.......................................................................
 * Constructor.
 */
PgUtil::PgUtil() {}

/**.......................................................................
 * Destructor.
 */
PgUtil::~PgUtil() {}

void PgUtil::setInteractive(bool inter)
{
  interactive_ = inter;
}

void PgUtil::contour(std::vector<double>& zdata, int nx,int ny, 
		     double xmina, double xmaxa, double ymina, double ymaxa, 
		     double *flag, double z1, double z2, unsigned ncontour,
		     std::string xlab, std::string ylab, std::string title, std::string unit)
{
  return contour(zdata.size(), &zdata[0], nx, ny,
		 xmina, xmaxa, ymina, ymaxa, 
		 flag, z1, z2, ncontour,
		 xlab, ylab, title, unit);
}

void PgUtil::contour(int ndata, double *zdata, int nx,int ny, 
		     double xmina, double xmaxa, double ymina, double ymaxa, 
		     double *flag, double z1, double z2, unsigned ncontour,
		     std::string xlab, std::string ylab, std::string title, std::string unit)
{
  std::vector<float> fdata(ndata);
  float fflag = (flag ? *flag : 0);
  
  for(unsigned i=0; i < ndata; i++) {
    fdata[i] = zdata[i];
  }

  return contour(ndata, &fdata[0], nx, ny,
		 xmina, xmaxa, ymina, ymaxa, 
		 &fflag, z1, z2, ncontour,
		 xlab, ylab, title, unit);
}

void PgUtil::contour(int ndata, float *zdata, int nx,int ny, 
		     float xmina, float xmaxa, float ymina, float ymaxa, 
		     float *flag, float z1, float z2, unsigned ncontour,
		     std::string xlab, std::string ylab, std::string title, std::string unit)
{
#if 0
  int status = v_contour(ndata, zdata, nx, ny,
			 xmina, xmaxa, ymina, ymaxa,
			 flag, z1, z2, ncontour,
			 (char*)xlab.c_str(), (char*)ylab.c_str(), (char*)title.c_str(), (char*)unit.c_str());

  if(status) {
    ThrowError("Error occurred in contour");
  }
#else
  PgUtil pgUtil;
  ContourPlot gp(&pgUtil,
		 ndata, zdata, nx, ny,
		 xmina, xmaxa, ymina, ymaxa,
		 flag, z1, z2, ncontour, false,
		 (char*)xlab.c_str(), (char*)ylab.c_str(), (char*)title.c_str(), (char*)unit.c_str(), true);
  
  gp.plot();
#endif
}

/**.......................................................................
 * Grid data & make a grayscale map 
 *
 * Input: 
 *  ferret    Ferret *  The ferret to be plotted.
 *  xmem        Dmem *  The xmember.
 *  ymem        Dmem *  The ymember.
 *  zmem        Dmem *  The zmember.
 *  nx           int    The number of points to use in x.
 *  ny           int    The number of points to use in y.
 *  z1         float    The foreground greyscale.
 *  z2         float    The background greyscale.
 * Output:
 *  return    int       0 - OK.
 */
void PgUtil::greyScale(std::vector<double>& zdata, int nx,int ny, 
		       double xmina, double xmaxa, double ymina, double ymaxa, 
		       double *flag, double z1, double z2, 
		       std::string xlab, std::string ylab, std::string title, std::string unit)
{
  return greyScale(zdata.size(), &zdata[0], nx, ny,
		   xmina, xmaxa, ymina, ymaxa, 
		   flag, z1, z2, 
		   xlab, ylab, title, unit);
}

void PgUtil::greyScale(int ndata, double *zdata, int nx,int ny, 
		       double xmina, double xmaxa, double ymina, double ymaxa, 
		       double *flag, double z1, double z2, 
		       std::string xlab, std::string ylab, std::string title, std::string unit)
{
  std::vector<float> fdata(ndata);
  float fflag = (flag ? *flag : 0);

  for(unsigned i=0; i < ndata; i++) {
    fdata[i] = zdata[i];
  }

  return greyScale(ndata, &fdata[0], nx, ny,
		   xmina, xmaxa, ymina, ymaxa, 
		   &fflag, z1, z2, 
		   xlab, ylab, title, unit);
}

void PgUtil::grayScale(int ndata, float *zdata, int nx,int ny, 
		       float xmina, float xmaxa, float ymina, float ymaxa, 
		       float *flag, float z1, float z2, 
		       std::string xlab, std::string ylab, std::string title, std::string unit)
{
  greyScale(ndata, zdata, nx, ny,
	    xmina, xmaxa, ymina, ymaxa,
	    flag, z1, z2,
	    xlab, ylab, title, unit);
}

void PgUtil::greyScale(int ndata, float *zdata, int nx,int ny, 
		       float xmina, float xmaxa, float ymina, float ymaxa, 
		       float *flag, float z1, float z2, 
		       std::string xlab, std::string ylab, std::string title, std::string unit)
{
  PgUtil pgUtil;

  GreyscalePlot gp(&pgUtil,
		   ndata, zdata, nx, ny,
		   xmina, xmaxa, ymina, ymaxa,
		   flag, z1, z2,
		   (char*)xlab.c_str(), (char*)ylab.c_str(), (char*)title.c_str(), (char*)unit.c_str());

  gp.plot();
}

int PgUtil::v_grey2(int ndata, float *zdata, int nx,int ny, 
		    float xmina, float xmaxa, float ymina, float ymaxa, 
		    float *flag,float z1, float z2, 
		    char *xlab, char *ylab, char *title, char *unit)
{
  PgModelManager modelManager = pgManager_;

  float xmins,xmaxs,xmin,xmax,dx;
  float ymins, ymaxs,ymin,ymax,dy;
  float zmin, zmax;
  int i,j;
  char answer[100];
  bool docurs=0,wasopen=0;
  int first = 1;
  float tr[6];
  int slen;
  int i1,i2,j1,j2,nbin;
  float xtemp,ytemp,xpos[2],ypos[2];
  static string mess1("\n %c - Select start and end vertices using this key\n");
  static string mess2(" %c - Select the whole plot\n");
  static string mess3(" %c - Abort selection\n");
  int cancel,dofull,accepted;
  int iter;
  int oldcol;

  float bright   =  0.5;
  float contrast = -1.0;

  bool read = true;
  Cmap *grey=0, *rain=0, *cmap=0, *heat=0;
  int n,xind,yind,exc;
  float mean,sd,min,max;
  float x,y,rad,rad0,r1,r2,r,h;
  float xmid,ymid,xmid0,ymid0;
  float sig;
  int autoscale=(z1==z2);

  enum {                      /* Function keys in the gaphical interface */
    G_CUR   = 'A',
    G_FLG   = 'B',            /* Toggle display flagged data. */
    G_CUT   = 'C',
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
    G_HEAT  = 'Q',
    G_RAIN  = 'R',
    G_STAT  = 'S',

    G_HORI  = 'U',
    G_VERT  = 'V',

    G_QUIT  = 'X',
    G_YSC   = 'Y',
    G_ZOOM  = 'Z',
    G_CROS  = '+',
  };	    

  char key = G_DIS;

  for(i=0;i < n_std_cmap;i++) {
    if(strcmp(std_cmaps[i].name.c_str(),"grey")==0) 
      grey = &std_cmaps[i];
    if(strcmp(std_cmaps[i].name.c_str(),"rainbow")==0) 
      rain = &std_cmaps[i];
    if(strcmp(std_cmaps[i].name.c_str(),"heat")==0) 
      heat = &std_cmaps[i];
  };

  if(autoscale) {
    zmin = zmax = zdata[0];
    for(i=0;i < ndata;i++) {
      zmin = MIN(zmin, zdata[i]);
      zmax = MAX(zmax, zdata[i]);
    }
  
    if(zmin==zmax) {
      zmin -= 0.1*zmin;
      zmax += 0.1*zmax;
    }
  } else {
    zmin=z1;
    zmax=z2;
  }
  
  xmins = xmina;
  xmaxs = xmaxa;
  ymins = ymina;
  ymaxs = ymaxa;

  xmid = (xmaxs+xmins)/2;
  ymid = (ymaxs+ymins)/2;

  // Store these for later use.

  xmid0 = xmid;
  ymid0 = ymid;
  rad0 = sqrt((xmaxs-xmid)*(xmaxs-xmid) + (ymaxs-ymid)*(ymaxs-ymid));

  dx = (xmaxs-xmins)/(nx-1);
  dy = (ymaxs-ymins)/(ny-1);

  // Set the transformation matrix for the data array.
  
  if(reverseX_) {
    tr[0]=xmaxs+dx;
    tr[1]=-dx;
  } else {
    tr[0]=xmins-dx;
    tr[1]=dx;
  }

  tr[2]=0.0;
  tr[3]=ymins-dy;
  tr[4]=0.0;
  tr[5]=dy;

  i1 = j1 = 1;
  i2 = nx;
  j2 = ny;

  queryDevice(wasopen);
  docurs = haveCursor();

  cpgsch(ch_);

  // Now that the transformation matrix has been set up, reset the display limits if requested

  if(usedefs_) {
    xmins = (xmin_ == xmax_) ? xmins : xmin_;
    xmaxs = (xmin_ == xmax_) ? xmaxs : xmax_;
    ymins = (ymin_ == ymax_) ? ymins : ymin_;
    ymaxs = (ymin_ == ymax_) ? ymaxs : ymax_;
  }

  if(win_) {
    cpgswin(xmins,xmaxs,ymins,ymaxs);
  }

  cmap = cmap_;

  // Expand the plot limits

  xmins -= dx/2;
  xmaxs += dx/2;
  ymins -= dy/2;
  ymaxs += dy/2;

  xmin = xmins;
  xmax = xmaxs;
  ymin = ymins;
  ymax = ymaxs;

  unsigned cursor = B_NORM;
  Trans trans;

  if(docurs) {

    fprintf(stdout,"For HELP, hit the \'%c\' key on your keyboard\n", G_HELP);

    do {
      cancel = 0;

      trans.dx_    = fabs(dx);
      trans.dy_    = fabs(dy);
      trans.nx_    = nx;
      trans.ny_    = ny;
      trans.ndata_ = ndata;
      trans.xmins_ = xmins;
      trans.ymins_ = ymins;
      trans.xmaxs_ = xmaxs;
      trans.ymaxs_ = ymaxs;
      trans.zdata_ = zdata;

      switch(key) {

	//------------------------------------------------------------
	// Make a radial plot of the image.
	//------------------------------------------------------------

      case G_RAD:
	cpgqci(&oldcol);
	cpgsci(1);
	dofull = 0;
	cancel = 0;
	for(iter = 0;iter<2 && !dofull && !cancel;iter++) {
	  do {
	    accepted = 0;
	    cpgband(B_LINE, 0, xmid, ymid, &xtemp,&ytemp, &key);
	    if(islower((int) key))
	      key = (char) toupper((int) key);
	    xpos[iter] = xtemp;
	    ypos[iter] = ytemp;
	    switch(key) {
	    case G_CAN:      /* Abort box selection */
	      accepted = 1;
	      exc = 0;
	      break;
	    case G_QUIT:     /* Quit now */
	      accepted = cancel = 1;
	      break;
	    case G_CUR:             /* Accept the selected start vertex */
	      accepted=1;
	      exc = 1;
	      break;
	    default:            /* Unexpected cursor input key - show usage */
	      fprintf(stdout,mess1.c_str(), G_CUR);
	      fprintf(stdout,mess2.c_str(), G_ZOOM);
	      fprintf(stdout,mess3.c_str(), G_CAN);
	      break;
	    };
	  } while(!accepted);
	};
	if(xpos[0] > xpos[1]) {
	  xtemp = xpos[1];
	  xpos[1] = xpos[0];
	  xpos[0] = xtemp;
	}
	if(ypos[0] > ypos[1]) {
	  ytemp = ypos[1];
	  ypos[1] = ypos[0];
	  ypos[0] = ytemp;
	}
	
	// Get the length of this radial segment.

	rad = sqrt((xpos[1]-xpos[0])*(xpos[1]-xpos[0])+(ypos[1]-ypos[0])*(ypos[1]-ypos[0]));
	
	// Use a number of radial bins proportional to the distance
	// (max will be the native image resolution, ie, ngrid/2

	nbin = (int)(nx/2*rad/rad0);

	r1 = sqrt((xpos[0]-xmid0)*(xpos[0]-xmid0)+(ypos[0]-ymid0)*(ypos[0]-ymid0));
	r2 = sqrt((xpos[1]-xmid0)*(xpos[1]-xmid0)+(ypos[1]-ymid0)*(ypos[1]-ymid0));

  	v_radplot(zdata,nbin,r1<r2?r1:r2,r1>r2?r1:r2,xmins,xmaxs,ymins,ymaxs,nx,ny);

	break;

	//------------------------------------------------------------
	// Zoom the plot
	//------------------------------------------------------------

      case G_ZOOM:
	cpgqci(&oldcol);
	cpgsci(5);
	dofull = 0;
	cancel = 0;
	for(iter = 0;iter<2 && !dofull && !cancel;iter++) {
	  do {
	    accepted = 0;
	    cpgband((iter==0) ? B_NORM : B_RECT, 0, xtemp, ytemp, &xtemp,
		    &ytemp, &key);
	    if(islower((int) key))
	      key = (char) toupper((int) key);
	    xpos[iter] = xtemp;
	    ypos[iter] = ytemp;
	    switch(key) {
	    case G_ZOOM:
	      accepted = dofull = 1;
	      break;
	    case G_CAN:      /* Abort box selection */
	      accepted = cancel = 1;
	      break;
	    case G_QUIT:     /* Quit now */
	      if(!wasopen)
		cpgend();
	      return 0;
	      break;
	    case G_CUR:             /* Accept the selected start vertex */
	      accepted=1;
	      break;
	    default:            /* Unexpected cursor input key - show usage */
	      fprintf(stdout,mess1.c_str(), G_CUR);
	      fprintf(stdout,mess2.c_str(), G_ZOOM);
	      fprintf(stdout,mess3.c_str(), G_CAN);
	      break;
	    };
	  } while(!accepted);
	};

	if(dofull) {
	  xmin = xmins;
	  ymin = ymins;
	  xmax = xmaxs;
	  ymax = ymaxs;
	} else {
	  
	  // Only reverse the boundaries if the "min" and "max" are
	  // contrary to the sense of dx and dy.

	  if((xpos[0] < xpos[1] && dx > 0.0) || 
	     (xpos[0] > xpos[1] && dx < 0.0)) {
	    xmin = xpos[0];
	    xmax = xpos[1];
	  }
	  else {
	    xmin = xpos[1];
	    xmax = xpos[0];
	  }
	  if((ypos[0] < ypos[1] && dy > 0.0) || 
	     (ypos[0] > ypos[1] && dy < 0.0)){
	    ymin = ypos[0];
	    ymax = ypos[1];
	  }
	  else {
	    ymin = ypos[1];
	    ymax = ypos[0];
	  }
	}
	
	// Recompute the midpoint of the displayed image.

	xmid = (xmin+xmax)/2;
	ymid = (ymin+ymax)/2;
	
	// Compute the new greyscale boundaries

	if(autoscale || !dofull){
	  float x,y;
	  int first=1,ind;

	  for(i=0;i < nx;i++)
	    for(j=0;j < ny;j++) {
	      x = xmins+dx*i;
	      y = ymins+dy*j;
	      ind = i+j*nx;

	      if(x >= xmin && x <= xmax && y >= ymin && y <= ymax) {
		if(first) {
		  zmin = zmax = zdata[ind];
		  first = 0;
		}
		zmin = MIN(zmin, zdata[ind]);
		zmax = MAX(zmax, zdata[ind]);
	      }
	    }
	} else {
	  zmin = z1;
	  zmax = z2;
	}

	cpgswin(xmin+dx/2,xmax+dx/2,ymin+dy/2,ymax+dy/2);
	cpgsci(oldcol);

      case G_DIS:
	if(!cancel) {

	  setupPlotBoundaries(xmin, xmax, ymin, ymax);

	  cpgctab(cmap->l,cmap->r,cmap->g,cmap->b,cmap->n,contrast,bright);

	  // Override min/max greyscale with user-specified values
	  
	  if(zmin_ != zmax_) {
	    zmin = zmin_;
	    zmax = zmax_;
	  }

	  cpgimag(zdata,nx,ny,i1,i2,j1,j2,zmax,zmin,tr);
	  cpgsci(1);

	  //------------------------------------------------------------
	  // Draw models on top, if displaying models
	  //------------------------------------------------------------
	  
	  if(modelManager.display_) {
	    modelManager.display();
	  }

	  drawBox(1);
	  drawLabels(xlab, ylab, title);
	  drawWedge(zmin, zmax, unit);

	  // Draw a header if using a header

	  drawHeader();
	};
	break;
	
	// Compute statistics on a selected region of the plot.

      case G_STAT:
	cpgqci(&oldcol);
	cpgsci(5);
	dofull = 0;
	cancel = 0;
	for(iter = 0;iter<2 && !dofull && !cancel;iter++) {
	  do {
	    accepted = 0;
	    cpgband((iter==0) ? B_NORM : B_RECT, 0, xtemp, ytemp, &xtemp,
		    &ytemp, &key);
	    if(islower((int) key))
	      key = (char) toupper((int) key);
	    xpos[iter] = xtemp;
	    ypos[iter] = ytemp;
	    switch(key) {
	    case G_STAT:
	      accepted = dofull = 1;
	      break;
	    case G_CAN:      /* Abort box selection */
	      accepted = cancel = 1;
	      break;
	    case G_QUIT:     /* Quit now */
	      /* 
	       * Close PGPLOT only if it was opened in this function 
	       */
	      if(!wasopen)
		cpgend(); 
	      return 0;
	      break;
	    case G_CUR:             /* Accept the selected start vertex */
	      accepted=1;
	      break;
	    default:            /* Unexpected cursor input key - show usage */
	      fprintf(stdout,mess1.c_str(), G_CUR);
	      fprintf(stdout,mess2.c_str(), G_ZOOM);
	      fprintf(stdout,mess3.c_str(), G_CAN);
	      break;
	    };
	  } while(!accepted);
	};

	if(dofull) {
	  xpos[0] = xmins;
	  ypos[0] = ymins;
	  xpos[1] = xmaxs;
	  ypos[1] = ymaxs;
	}



#if 0
	/*
	 * Here we want xpos[0] and xpos[1], etc. to be the absolute 
	 * minimum and maximum, since we test if a data point falls between 
	 * these values.
	 */
	if(xpos[0] > xpos[1]) {
	  xtemp = xpos[1];
	  xpos[1] = xpos[0];
	  xpos[0] = xtemp;
	}
	if(ypos[0] > ypos[1]) {
	  ytemp = ypos[1];
	  ypos[1] = ypos[0];
	  ypos[0] = ytemp;
	}
	mean = 0.0;
	n = 0;
	for(i=0;i < ndata;i++) {
	  yind = i/nx;
	  xind = i - yind*nx;
	  xtemp = xmins + dx/2 + xind*dx;
	  ytemp = ymins + dy/2 + yind*dy;
	  if(xtemp >= xpos[0] && xtemp <= xpos[1] && ytemp >= ypos[0] && 
	     ytemp <= ypos[1]) {
	    if(first) {
	      min = max = zdata[i];
	      first = 0;
	    }
	    min = MIN(zdata[i],min);
	    max = MAX(zdata[i],max);
	    mean += (zdata[i] - mean)/(n+1);
	    ++n;
	  }
	}
	sd = 0.0;
	n = 0;
	for(i=0;i < ndata;i++) {
	  yind = i/nx;
	  xind = i - yind*nx;
	  xtemp = xmins + dx/2 + xind*dx;
	  ytemp = ymins + dy/2 + yind*dy;
	  if(xtemp >= xpos[0] && xtemp <= xpos[1] && ytemp >= ypos[0] && 
	     ytemp <= ypos[1]) {
	    sd += ((zdata[i] - mean)*(zdata[i]-mean) - sd)/(n+1);
	    ++n;
	  }
	}
	if(n > 1)
	  sd = sqrt(sd*n/(n-1));
	else sd = 0.0f;
	fprintf(stdout, "\n\n\t\tmean\t=\t%g\n\t\tsd\t=\t%g\n\t\tmin\t=\t%g\n\t\tmax\t=\t%g\n\t\tnpts\t=\t%d\n", mean, sd, min, max, n);
	first = 1;
#else
	trans.printStats(xpos[0], xpos[1], ypos[0], ypos[1]);
#endif
	cpgsci(oldcol);
	break;

	//------------------------------------------------------------
	// Define a new model component
	//------------------------------------------------------------

      case G_MODEL:
	modelManager.getModel(xpos[0], ypos[0], read, key, unit, trans);
	break;
	
	//------------------------------------------------------------
	// Print information about the nearest point
	//------------------------------------------------------------

      case G_INS:
	/*
	 * Find the nearest grid point
	 */
	{
#if 0
	  bool first=true;
	  float dist,distmin,xtmp,ytmp,val;
	  int xmin,ymin;
	  for(i=0;i < ndata;i++) {
	    yind = i/nx;
	    xind = i - yind*nx;
	    xtmp = xmins + dx/2 + xind*dx;
	    ytmp = ymins + dy/2 + yind*dy;
	    dist = sqrt((xtmp-xpos[0])*(xtmp-xpos[0]) + 
			(ytmp-ypos[0])*(ytmp-ypos[0]));
	    if(dist < distmin || first) {
	      xmin = xind;
	      ymin = yind;
	      distmin = dist;
	      first = false;
	      val = zdata[i];
	    }
	  }
#else
	  float val = trans.valNearestToPoint(xpos[0], ypos[0]);
#endif

	  if(coordCallback_) {
	    std::string xstr, ystr;
	    (*coordCallback_)(xpos[0], ypos[0], xstr, ystr, coordCallbackArgs_);
	    fprintf(stdout,"Pixel value at (%s, %s) is %f %s\n", xstr.c_str(), ystr.c_str(), val, unit);
	  } else {
	    fprintf(stdout,"Pixel value at (%lf, %lf) is %f %s\n", xmin, ymin, val, unit);
	  }

	}
	break;
      case G_HELP:     /* Print usage info */
	fprintf(stdout,"\nYou requested help by pressing \'%c\'.\n", G_HELP);
	fprintf(stdout,"All cursor positions are entered with \'%c\' key (Left mouse button)\n", G_CUR);
	fprintf(stdout,"\n %c - Select a sub-image to be displayed\n", G_ZOOM);
	fprintf(stdout," %c - Calculate statistics on a subset of data\n", G_STAT);
	fprintf(stdout," %c - Redisplay current plot\n", G_DIS);
	fprintf(stdout," %c - Add/remove/display models\n", G_MODEL);
	fprintf(stdout," %c - Fiddle contrast & brightness\n", G_FID);
	fprintf(stdout," %c - Use greyscale\n", G_GREY);
	fprintf(stdout," %c - Use rainbow colormap\n", G_RAIN);
	fprintf(stdout," %c - Use heat colormap\n", G_HEAT);
	fprintf(stdout," %c - Toggle crosshair cursor\n", G_CROS);
	fprintf(stdout,"\nTo end this session hit the \'%c\' key (Right mouse button)\n", G_QUIT);
	fprintf(stdout,"\n");
	break;
      default :
	break;
      case G_GREY:
	cmap = grey;
	read = false;
	key = G_DIS;
	break;
      case G_RAIN:
	cmap = rain;
	read = false;
	key = G_DIS;
	break;
      case G_HEAT:
	cmap = heat; 
	read = false;
	key = G_DIS;
	break;
      case G_CROS:
	if(cursor == B_CROSS)
	  cursor = B_NORM;
	else
	  cursor = B_CROSS;
	break;
      case G_FID:
	{
	  contrast = 5.0 * (ypos[0]-ymid)/(ypos[0] < ymid ? (ymin-ymid) : -(ymax-ymid));
	  bright = 0.5 + 1.0 * (fabs(contrast)+1.0)*((xpos[0] - xmax)/(xmin - xmax) - 0.5);
	  cpgctab(cmap->l,cmap->r,cmap->g,cmap->b,cmap->n,contrast,bright);
	  
	  cpgband(cursor, 0, xpos[0], ypos[0], &xpos[0], &ypos[0], &key);
	  key = G_DIS;
	  read = false;
#if 0
	  do {
	    contrast = 5.0 * (ypos[0]-ymid)/(ypos[0] < ymid ? (ymin-ymid) : -(ymax-ymid));
	    bright = 0.5 + 1.0 * (fabs(contrast)+1.0)*((xpos[0] - xmax)/(xmin - xmax) - 0.5);
	    cpgctab(cmap->l,cmap->r,cmap->g,cmap->b,cmap->n,contrast,bright);
	    
	    cpgband(cursor, 0, xpos[0], ypos[0], &xpos[0], &ypos[0], &key);
	    if(islower((int) key))
	      key = (char) toupper((int) key);
	  } while(key == G_FID);
	  read = false;
#endif
	}


      } // End of switch

      if(read) 
	cpgband(cursor, 0, xpos[0], ypos[0], &xpos[0], &ypos[0], &key);

      read = true;

      if(islower((int) key))
	key = (char) toupper((int) key);

    } while(key != G_QUIT);
  }
  /*
   * Else if no cursor, just plot and exit.
   */
  else {
    
    if(!overplot_)
      cpgpage();

    if(vp_) {
      cpgvstd();
    }

    
    if(wnad_) {
      if(win_)
	cpgswin(0,1,0,1);
      cpgwnad(0,1,0,1); 
    }

    if(win_) {
      cpgswin(xmins,xmaxs,ymins,ymaxs);
    }

    if(wnad_) {
      cpgwnad(xmins,xmaxs,ymins,ymaxs);
    }

    if(!autoscale){
      zmin=z1;
      zmax=z2;
    }

    cpgctab(cmap->l,cmap->r,cmap->g,cmap->b,cmap->n,contrast,bright);

    // Override min/max greyscale with user-specified values

    if(zmin_ != zmax_) {
      zmin = zmin_;
      zmax = zmax_;
    }

    cpgimag(zdata,nx,ny,i1,i2,j1,j2,zmax,zmin,tr);

    if(modelManager.display_) {
      modelManager.display();
    }

    cpgsci(1);

    std::string xOpts, yOpts;

    xOpts = (tick_ && xTick_) ? (xTickLabeling_ ? "BCNST" : "BCST") : "BC";
    yOpts = (tick_ && yTick_) ? (yTickLabeling_ ? "BCNST" : "BCST") : "BC";

    if(box_) {

      if(useBoxCi_)
	cpgsci(boxCi_);
      else
	cpgsci(1);

      cpgbox(xOpts.c_str(), 0.0,0, yOpts.c_str(),0.0,0);
    }
    
    std::string xLab, yLab, tLab;

    xLab = (label_ && xLabel_) ? xlab : "";
    yLab = (label_ && yLabel_) ? ylab : "";
    tLab = title_ ? title : "";
    
    cpglab(xLab.c_str(), yLab.c_str(), tLab.c_str());
	  
    /*
     * Draw a Wedge on the side.
     */
    if(wedge_)
      cpgwedg("RI",0,4,zmax,zmin,unit);

    // Draw a header if using a header

    drawHeader();
  };
  /* 
   * Close PGPLOT only if it was opened in this function 
   */
  if(!wasopen) {
    cpgend(); 
  }

  return 0;
}

void PgUtil::drawHeader()
{
  if(!useHeader_)
    return;

  float x1,x2,y1,y2;
  cpgqwin(&x1,&x2,&y1,&y2);

  cpgptxt(x1, y2 + 0.05*(y2-y1), 0.0, 0.0, header_.text_.c_str());
}

void PgUtil::text(std::string text, double x, double y, Angle angle, Just just)
{
  cpgptxt(x, y, angle.degrees(), pgplotJust(just), text.c_str());
}

float PgUtil::pgplotJust(Just just)
{
  switch(just)  {
  case JUST_LEFT:
    return 0.0;
    break;
  case JUST_CENTER:
    return 0.5;
    break;
  default:
    return 1.0;
    break;
  }
}

int PgUtil::v_contour(int ndata, float *zdata, int nx,int ny, 
		      float xmina, float xmaxa, float ymina, float ymaxa, 
		      float *flag, float z1, float z2, unsigned ncontour,
		      char *xlab, char *ylab, char *title, char *unit)
{
  float xmins,xmaxs,xmin,xmax,dx;
  float ymins, ymaxs,ymin,ymax,dy;
  float zmin, zmax;
  int i,j;
  char answer[100];
  bool docurs=0,wasopen=0;
  int first = 1;
  float tr[6];
  int slen;
  int i1,i2,j1,j2,nbin;
  float xtemp,ytemp,xpos[2],ypos[2];
  static string mess1("\n %c - Select start and end vertices using this key\n");
  static string mess2(" %c - Select the whole plot\n");
  static string mess3(" %c - Abort selection\n");
  int cancel,dofull,accepted;
  int iter;
  char key = 'L';
  int oldcol;

  float bright   =  0.5;
  float contrast = -1.0;

  int read = 1;
  Cmap *grey=0, *rain=0, *cmap=0, *heat=0;
  int n,xind,yind,exc;
  float mean,sd,min,max;
  float x,y,rad,rad0,r1,r2,r,h;
  float xmid,ymid,xmid0,ymid0;
  float sig;
  int autoscale=(z1==z2);

  enum {                      /* Function keys in the gaphical interface */
    G_CUR   = 'A',
    G_FLG   = 'B',            /* Toggle display flagged data. */
    G_CUT   = 'C',
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
    G_COPY  = 'P',
    G_HEAT  = 'Q',
    G_RAIN  = 'R',
    G_STAT  = 'S',

    G_HORI  = 'U',
    G_VERT  = 'V',

    G_QUIT  = 'X',
    G_YSC   = 'Y',
    G_ZOOM  = 'Z',
  };	    

  for(i=0;i < n_std_cmap;i++) {
    if(strcmp(std_cmaps[i].name.c_str(),"grey")==0) 
      grey = &std_cmaps[i];
    if(strcmp(std_cmaps[i].name.c_str(),"rainbow")==0) 
      rain = &std_cmaps[i];
    if(strcmp(std_cmaps[i].name.c_str(),"heat")==0) 
      heat = &std_cmaps[i];
  };

  if(autoscale) {
    zmin = zmax = zdata[0];
    for(i=0;i < ndata;i++) {
      zmin = MIN(zmin, zdata[i]);
      zmax = MAX(zmax, zdata[i]);
    }
  
    if(zmin==zmax) {
      zmin -= 0.1*zmin;
      zmax += 0.1*zmax;
    }
  } else {
    zmin=z1;
    zmax=z2;
  }
  
  if(nContour_ > 0) {
    ncontour = nContour_;
  }

  std::vector<float> contours(ncontour);
  float dz = (zmax - zmin)/ncontour;

  for(unsigned icont=0; icont < ncontour; icont++) {
    contours[icont] = zmin + icont * dz;
  }
  
  xmins = xmina;
  xmaxs = xmaxa;
  ymins = ymina;
  ymaxs = ymaxa;

  xmid = (xmaxs+xmins)/2;
  ymid = (ymaxs+ymins)/2;

  /**
   * Store these for later use.
   */
  xmid0 = xmid;
  ymid0 = ymid;
  rad0 = sqrt((xmaxs-xmid)*(xmaxs-xmid) + (ymaxs-ymid)*(ymaxs-ymid));

  dx = (xmaxs-xmins)/(nx-1);
  dy = (ymaxs-ymins)/(ny-1);

  /**
   * Set the transformation matrix for the data array.
   */
  tr[0]=xmins-dx;
  tr[1]=dx;
  tr[2]=0.0;
  tr[3]=ymins-dy;
  tr[4]=0.0;
  tr[5]=dy;
  

  i1 = j1 = 1;
  i2 = nx;
  j2 = ny;

  queryDevice(wasopen);
  docurs = haveCursor();

  cpgswin(xmins,xmaxs,ymins,ymaxs);

  cmap = cmap_;
  
  // Expand the plot limits

  xmins -= dx/2;
  xmaxs += dx/2;
  ymins -= dy/2;
  ymaxs += dy/2;

  xmin = xmins;
  xmax = xmaxs;
  ymin = ymins;
  ymax = ymaxs;

  if(docurs) {
    fprintf(stdout,"For HELP, hit the \'%c\' key on your keyboard\n", G_HELP);
    do {
      cancel = 0;
      switch(key) {
      case G_CUR:
	cpgqci(&oldcol);
	cpgsci(1);
	dofull = 0;
	cancel = 0;
	for(iter = 0;iter<2 && !dofull && !cancel;iter++) {
	  do {
	    accepted = 0;
	    cpgband(B_LINE, 0, xmid, ymid, &xtemp,&ytemp, &key);
	    if(islower((int) key))
	      key = (char) toupper((int) key);
	    xpos[iter] = xtemp;
	    ypos[iter] = ytemp;
	    switch(key) {
	    case G_CAN:      /* Abort box selection */
	      accepted = 1;
	      exc = 0;
	      break;
	    case G_QUIT:     /* Quit now */
	      accepted = cancel = 1;
	      break;
	    case G_CUR:             /* Accept the selected start vertex */
	      accepted=1;
	      exc = 1;
	      break;
	    default:            /* Unexpected cursor input key - show usage */
	      fprintf(stdout,mess1.c_str(), G_CUR);
	      fprintf(stdout,mess2.c_str(), G_ZOOM);
	      fprintf(stdout,mess3.c_str(), G_CAN);
	      break;
	    };
	  } while(!accepted);
	};
	/*
	 * Distance from the origin.
	 */
	rad = sqrt((xmid-xpos[0])*(xmid-xpos[0])+(ymid-ypos[0])*(ymid-ypos[0]));
	/*
	 * Width of the filter.
	 */
	sig = 0.5*sqrt((xpos[0]-xpos[1])*(xpos[0]-xpos[1])+(ypos[0]-ypos[1])*(ypos[0]-ypos[1]));
	/*
	 * Zap the requested points.
	 */
	if(flag!=NULL) {
	  first = 1;
	  for(j=0;j < ny;j++) 
	    for(i=0;i < nx;i++) {
	      x = xmins+dx*i;
	      y = ymins+dy*j;
	      r = sqrt((xmid-x)*(xmid-x)+(ymid-y)*(ymid-y));
	      h = 1.0/(1+(r*sig/(r*r-rad*rad))*(r*sig/(r*r-rad*rad)));
	      flag[i+nx*j] = exc ? h : -(h-1);
	      zdata[i+nx*j] *= flag[i+nx*j];
	      if(autoscale) {
		if(first) {
		  zmin = zmax = zdata[i+nx*j];
		  first = 0;
		}
		zmin = MIN(zdata[i+nx*j],zmin);
		zmax = MAX(zdata[i+nx*j],zmax);
	      }
	    }  
	}
	break;
      case G_CUT:
	cpgqci(&oldcol);
	cpgsci(5);
	dofull = 0;
	cancel = 0;
	for(iter = 0;iter<2 && !dofull && !cancel;iter++) {
	  do {
	    accepted = 0;
	    cpgband((iter==0 ? B_NORM : B_LINE), 0, xtemp, ytemp, &xtemp, &ytemp, &key);
	    if(islower((int) key))
	      key = (char) toupper((int) key);
	    xpos[iter] = xtemp;
	    ypos[iter] = ytemp;
	    switch(key) {
	    case G_CAN:      /* Abort box selection */
	      accepted = 1;
	      exc = 0;
	      break;
	    case G_QUIT:     /* Quit now */
	      accepted = cancel = 1;
	      break;
	    case G_CUR:             /* Accept the selected start vertex */
	      accepted=1;
	      exc = 1;
	      break;
	    default:            /* Unexpected cursor input key - show usage */
	      fprintf(stdout,mess1.c_str(), G_CUR);
	      fprintf(stdout,mess2.c_str(), G_ZOOM);
	      fprintf(stdout,mess3.c_str(), G_CAN);
	      break;
	    };
	  } while(!accepted);
	};
	sig = 2*sqrt((xpos[0]-xpos[1])*(xpos[0]-xpos[1])+(ypos[0]-ypos[1])*(ypos[0]-ypos[1]));
	/*
	 * Zap the requested points.
	 */
	if(flag!=NULL) {
	  first = 1;
	  for(j=0;j < ny;j++) 
	    for(i=0;i < nx;i++) {
	      x = xmins+dx*i;
	      y = ymins+dy*j;
	      r1 = sqrt((x-xpos[0])*(x-xpos[0])+(y-ypos[0])*(y-ypos[0]));
	      r2 = sqrt((2*xmid-x-xpos[0])*(2*xmid-x-xpos[0])+(2*ymid-y-ypos[0])*(2*ymid-y-ypos[0]));
	      h = 1.0/(1+(sig/r1))*1.0/(1+(sig/r2));
	      flag[i+nx*j] = h;
	      zdata[i+nx*j] *= flag[i+nx*j];
	      if(autoscale) {
		if(first) {
		  zmin = zmax = zdata[i+nx*j];
		  first = 0;
		}
		zmin = MIN(zdata[i+nx*j],zmin);
		zmax = MAX(zdata[i+nx*j],zmax);
	      }
	    } 
	}
	break;
	/*
	 * Make a radial plot of the image.
	 */
      case G_RAD:
	cpgqci(&oldcol);
	cpgsci(1);
	dofull = 0;
	cancel = 0;
	for(iter = 0;iter<2 && !dofull && !cancel;iter++) {
	  do {
	    accepted = 0;
	    cpgband(B_LINE, 0, xmid, ymid, &xtemp,&ytemp, &key);
	    if(islower((int) key))
	      key = (char) toupper((int) key);
	    xpos[iter] = xtemp;
	    ypos[iter] = ytemp;
	    switch(key) {
	    case G_CAN:      /* Abort box selection */
	      accepted = 1;
	      exc = 0;
	      break;
	    case G_QUIT:     /* Quit now */
	      accepted = cancel = 1;
	      break;
	    case G_CUR:             /* Accept the selected start vertex */
	      accepted=1;
	      exc = 1;
	      break;
	    default:            /* Unexpected cursor input key - show usage */
	      fprintf(stdout,mess1.c_str(), G_CUR);
	      fprintf(stdout,mess2.c_str(), G_ZOOM);
	      fprintf(stdout,mess3.c_str(), G_CAN);
	      break;
	    };
	  } while(!accepted);
	};
	if(xpos[0] > xpos[1]) {
	  xtemp = xpos[1];
	  xpos[1] = xpos[0];
	  xpos[0] = xtemp;
	}
	if(ypos[0] > ypos[1]) {
	  ytemp = ypos[1];
	  ypos[1] = ypos[0];
	  ypos[0] = ytemp;
	}
	/*
	 * Get the length of this radial segment.
	 */
	rad = sqrt((xpos[1]-xpos[0])*(xpos[1]-xpos[0])+(ypos[1]-ypos[0])*(ypos[1]-ypos[0]));
	/*
	 * Use a number of radial bins proportional to the distance (max will
	 * be the native image resolution, ie, ngrid/2
	 */
	nbin = (int)(nx/2*rad/rad0);

	r1 = sqrt((xpos[0]-xmid0)*(xpos[0]-xmid0)+(ypos[0]-ymid0)*(ypos[0]-ymid0));
	r2 = sqrt((xpos[1]-xmid0)*(xpos[1]-xmid0)+(ypos[1]-ymid0)*(ypos[1]-ymid0));

  	v_radplot(zdata,nbin,r1<r2?r1:r2,r1>r2?r1:r2,xmins,xmaxs,ymins,ymaxs,nx,ny);

	break;
      case G_YSC:
	cpgqci(&oldcol);
	cpgsci(5);
	dofull = 0;
	cancel = 0;
	for(iter = 0;iter<1 && !dofull && !cancel;iter++) {
	  do {
	    accepted = 0;
	    cpgband(B_YRNG, 0, xtemp, ytemp, &xtemp,
		    &ytemp, &key);
	    if(islower((int) key))
	      key = (char) toupper((int) key);
	    xpos[iter] = xtemp;
	    ypos[iter] = ytemp;
	    switch(key) {
	    case G_CAN:      /* Abort box selection */
	      accepted = 1;
	      exc = 0;
	      break;
	    case G_QUIT:     /* Quit now */
	      accepted = cancel = 1;
	      break;
	    case G_CUR:             /* Accept the selected start vertex */
	      accepted=1;
	      break;
	    default:            /* Unexpected cursor input key - show usage */
	      fprintf(stdout,mess1.c_str(), G_CUR);
	      fprintf(stdout,mess2.c_str(), G_ZOOM);
	      fprintf(stdout,mess3.c_str(), G_CAN);
	      break;
	    };
	  } while(!accepted);
	};
	sig = fabs(ypos[0]-ymid);
	/*
	 * Zap the requested points.
	 */
	if(flag!=NULL) {
	  first = 1;
	  for(j=0;j < ny;j++) 
	    for(i=0;i < nx;i++) {
	      y = ymins+dy*j;
	      if(fabs(y-ymid) <= sig)
		flag[i+nx*j] = 0;
	      zdata[i+nx*j] *= flag[i+nx*j];
	      if(first) {
		zmin = zmax = zdata[i+nx*j];
		first = 0;
	      }
	      zmin = MIN(zdata[i+nx*j],zmin);
	      zmax = MAX(zdata[i+nx*j],zmax);
	    }  
	}
	break;
      case G_HORI:
	cpgqci(&oldcol);
	cpgsci(5);
	dofull = 0;
	cancel = 0;
	for(iter = 0;iter<1 && !dofull && !cancel;iter++) {
	  do {
	    accepted = 0;
	    cpgband(B_XRNG, 0, xtemp, ytemp, &xtemp,
		    &ytemp, &key);
	    if(islower((int) key))
	      key = (char) toupper((int) key);
	    xpos[iter] = xtemp;
	    ypos[iter] = ytemp;
	    switch(key) {
	    case G_CAN:      /* Abort box selection */
	      accepted = 1;
	      exc = 0;
	      break;
	    case G_QUIT:     /* Quit now */
	      accepted = cancel = 1;
	      break;
	    case G_CUR:             /* Accept the selected start vertex */
	      accepted=1;
	      break;
	    default:            /* Unexpected cursor input key - show usage */
	      fprintf(stdout,mess1.c_str(), G_CUR);
	      fprintf(stdout,mess2.c_str(), G_ZOOM);
	      fprintf(stdout,mess3.c_str(), G_CAN);
	      break;
	    };
	  } while(!accepted);
	};
	sig = fabs(xmid-xpos[0]);
	/*
	 * Zap the requested points.
	 */
	if(flag!=NULL) {
	  first = 1;
	  for(j=0;j < ny;j++) 
	    for(i=0;i < nx;i++) {
	      x = xmins+dx*i;
	      if(fabs(x-xmid) <= sig)
		flag[i+nx*j] = 0;
	      zdata[i+nx*j] *= flag[i+nx*j];
	      if(first) {
		zmin = zmax = zdata[i+nx*j];
		first = 0;
	      }
	      zmin = MIN(zdata[i+nx*j],zmin);
	      zmax = MAX(zdata[i+nx*j],zmax);
	    }  
	}
	break;
      case G_FIT:
	cpgqci(&oldcol);
	cpgsci(5);
	dofull = 0;
	cancel = 0;
	for(iter = 0;iter<1 && !dofull && !cancel;iter++) {
	  do {
	    accepted = 0;
	    cpgband(B_LINE, 0, xmid, ymid, &xtemp, &ytemp, &key);
	    if(islower((int) key))
	      key = (char) toupper((int) key);
	    xpos[iter] = xtemp;
	    ypos[iter] = ytemp;
	    switch(key) {
	    case G_CAN:      /* Abort box selection */
	      accepted = 1;
	      exc = 0;
	      break;
	    case G_QUIT:     /* Quit now */
	      accepted = cancel = 1;
	      break;
	    case G_CUR:             /* Accept the selected start vertex */
	      accepted=1;
	      exc = 1;
	      break;
	    default:            /* Unexpected cursor input key - show usage */
	      fprintf(stdout,mess1.c_str(), G_CUR);
	      fprintf(stdout,mess2.c_str(), G_ZOOM);
	      fprintf(stdout,mess3.c_str(), G_CAN);
	      break;
	    };
	  } while(!accepted);
	};
	rad = sqrt((xpos[0]-xmid)*(xpos[0]-xmid)+(ypos[0]-ymid)*(ypos[0]-ymid));
	/*
	 * Zap the requested points.
	 */
	if(flag!=NULL) {
	  first = 1;
	  for(j=0;j < ny;j++) 
	    for(i=0;i < nx;i++) {
	      x = xmins+dx*i;
	      y = ymins+dy*j;
	      r = sqrt((x-xmid)*(x-xmid)+(y-ymid)*(y-ymid));
	      h = 1.0/(1+(rad/r)*(rad/r));
	      flag[i+nx*j] = exc ? h : -(h-1);
	      zdata[i+nx*j] *= flag[i+nx*j];
	      if(first) {
		zmin = zmax = zdata[i+nx*j];
		first = 0;
	      }
	      zmin = MIN(zdata[i+nx*j],zmin);
	      zmax = MAX(zdata[i+nx*j],zmax);
	    }  
	}
	break;
	/*
	 * Zoom the plot
	 */
      case G_ZOOM:
	cpgqci(&oldcol);
	cpgsci(5);
	dofull = 0;
	cancel = 0;
	for(iter = 0;iter<2 && !dofull && !cancel;iter++) {
	  do {
	    accepted = 0;
	    cpgband((iter==0) ? B_NORM : B_RECT, 0, xtemp, ytemp, &xtemp,
		    &ytemp, &key);
	    if(islower((int) key))
	      key = (char) toupper((int) key);
	    xpos[iter] = xtemp;
	    ypos[iter] = ytemp;
	    switch(key) {
	    case G_ZOOM:
	      accepted = dofull = 1;
	      break;
	    case G_CAN:      /* Abort box selection */
	      accepted = cancel = 1;
	      break;
	    case G_QUIT:     /* Quit now */
	      if(!wasopen)
		cpgend();
	      return 0;
	      break;
	    case G_CUR:             /* Accept the selected start vertex */
	      accepted=1;
	      break;
	    default:            /* Unexpected cursor input key - show usage */
	      fprintf(stdout,mess1.c_str(), G_CUR);
	      fprintf(stdout,mess2.c_str(), G_ZOOM);
	      fprintf(stdout,mess3.c_str(), G_CAN);
	      break;
	    };
	  } while(!accepted);
	};
	if(dofull) {
	  xmin = xmins;
	  ymin = ymins;
	  xmax = xmaxs;
	  ymax = ymaxs;
	}
	else {
	  /*
	   * Only reverse the boundaries if the "min" and "max" are contrary
	   * to the sense of dx and dy.
	   */
	  if((xpos[0] < xpos[1] && dx > 0.0) || 
	     (xpos[0] > xpos[1] && dx < 0.0)) {
	    xmin = xpos[0];
	    xmax = xpos[1];
	  }
	  else {
	    xmin = xpos[1];
	    xmax = xpos[0];
	  }
	  if((ypos[0] < ypos[1] && dy > 0.0) || 
	     (ypos[0] > ypos[1] && dy < 0.0)){
	    ymin = ypos[0];
	    ymax = ypos[1];
	  }
	  else {
	    ymin = ypos[1];
	    ymax = ypos[0];
	  }
	}
	/*
	 * Recompute the midpoint of the displayed image.
	 */
	xmid = (xmin+xmax)/2;
	ymid = (ymin+ymax)/2;
	/*
	 * Compute the new greyscale boundaries
	 */
	if(autoscale || !dofull){
	  float x,y;
	  int first=1,ind;

	  for(i=0;i < nx;i++)
	    for(j=0;j < ny;j++) {
	      x = xmins+dx*i;
	      y = ymins+dy*j;
	      ind = i+j*nx;

	      if(x >= xmin && x <= xmax && y >= ymin && y <= ymax) {
		if(first) {
		  zmin = zmax = zdata[ind];
		  first = 0;
		}
		zmin = MIN(zmin, zdata[ind]);
		zmax = MAX(zmax, zdata[ind]);
	      }
	    }
	}
	else {
	  zmin = z1;
	  zmax = z2;
	}

	cpgswin(xmin+dx/2,xmax+dx/2,ymin+dy/2,ymax+dy/2);
	cpgsci(oldcol);
      case G_DIS:
	if(!cancel) {

	  if(!overplot_) {
	    cpgpage();
	  }

	  if(vp_) {
	    cpgvstd();
	  }

	  if(wnad_) {
	    if(win_)
	      cpgswin(0,1,0,1);
	    cpgwnad(0,1,0,1); 
	  }

	  if(win_)
	    cpgswin(xmin-dx/2,xmax+dx/2,ymin-dy/2,ymax+dy/2);
	  
	  if(wnad_)
	    cpgwnad(xmin-dx/2,xmax+dx/2,ymin-dy/2,ymax+dy/2);

	  //	  cpgbbuf();

	  cpgctab(cmap->l,cmap->r,cmap->g,cmap->b,cmap->n,contrast,bright);

	  // Override min/max greyscale with user-specified values
	  
	  if(zmin_ != zmax_) {
	    zmin = zmin_;
	    zmax = zmax_;

	    float dz = (zmax - zmin)/ncontour;

	    for(unsigned icont=0; icont < ncontour; icont++) {
	      contours[icont] = zmin + icont * dz;
	    }
	  }

	  cpgcont(zdata,nx,ny,i1,i2,j1,j2, &contours[0], ncontour,tr);
	  cpgsci(1);

	  std::string xOpts, yOpts;

	  if(tick_) { 
	    xOpts = xTick_ ? (xTickLabeling_ ? "BCNST" : "BCST") : "BC";
	    yOpts = yTick_ ? (yTickLabeling_ ? "BCNST" : "BCST") : "BC";
	  }

	  if(box_)
	    cpgbox(xOpts.c_str(), 0.0,0, yOpts.c_str(),0.0,0);

	  std::string xLab, yLab, tLab;

	  xLab = (label_ && xLabel_) ? xlab : "";
	  yLab = (label_ && yLabel_) ? ylab : "";
	  tLab = title_ ? title : "";

	  cpglab(xLab.c_str(), yLab.c_str(), tLab.c_str());

	  // Draw a ramp on the side.

	  if(wedge_)
	    cpgwedg("RI",0,4,zmax,zmin,unit); 
	  //	  cpgebuf();
	};
	break;
	/*
	 * Compute statistics on a selected region of the plot.
	 */
      case G_STAT:
	cpgqci(&oldcol);
	cpgsci(5);
	dofull = 0;
	cancel = 0;
	for(iter = 0;iter<2 && !dofull && !cancel;iter++) {
	  do {
	    accepted = 0;
	    cpgband((iter==0) ? B_NORM : B_RECT, 0, xtemp, ytemp, &xtemp,
		    &ytemp, &key);
	    if(islower((int) key))
	      key = (char) toupper((int) key);
	    xpos[iter] = xtemp;
	    ypos[iter] = ytemp;
	    switch(key) {
	    case G_STAT:
	      accepted = dofull = 1;
	      break;
	    case G_CAN:      /* Abort box selection */
	      accepted = cancel = 1;
	      break;
	    case G_QUIT:     /* Quit now */
	      /* 
	       * Close PGPLOT only if it was opened in this function 
	       */
	      if(!wasopen)
		cpgend(); 
	      return 0;
	      break;
	    case G_CUR:             /* Accept the selected start vertex */
	      accepted=1;
	      break;
	    default:            /* Unexpected cursor input key - show usage */
	      fprintf(stdout,mess1.c_str(), G_CUR);
	      fprintf(stdout,mess2.c_str(), G_ZOOM);
	      fprintf(stdout,mess3.c_str(), G_CAN);
	      break;
	    };
	  } while(!accepted);
	};
	if(dofull) {
	  xpos[0] = xmins;
	  ypos[0] = ymins;
	  xpos[1] = xmaxs;
	  ypos[1] = ymaxs;
	}
	/*
	 * Here we want xpos[0] and xpos[1], etc. to be the absolute 
	 * minimum and maximum, since we test if a data point falls between 
	 * these values.
	 */
	if(xpos[0] > xpos[1]) {
	  xtemp = xpos[1];
	  xpos[1] = xpos[0];
	  xpos[0] = xtemp;
	}
	if(ypos[0] > ypos[1]) {
	  ytemp = ypos[1];
	  ypos[1] = ypos[0];
	  ypos[0] = ytemp;
	}
	mean = 0.0;
	n = 0;
	for(i=0;i < ndata;i++) {
	  yind = i/nx;
	  xind = i - yind*nx;
	  xtemp = xmins + dx/2 + xind*dx;
	  ytemp = ymins + dy/2 + yind*dy;
	  if(xtemp >= xpos[0] && xtemp <= xpos[1] && ytemp >= ypos[0] && 
	     ytemp <= ypos[1]) {
	    if(first) {
	      min = max = zdata[i];
	      first = 0;
	    }
	    min = MIN(zdata[i],min);
	    max = MAX(zdata[i],max);
	    mean += (zdata[i] - mean)/(n+1);
	    ++n;
	  }
	}
	sd = 0.0;
	n = 0;
	for(i=0;i < ndata;i++) {
	  yind = i/nx;
	  xind = i - yind*nx;
	  xtemp = xmins + dx/2 + xind*dx;
	  ytemp = ymins + dy/2 + yind*dy;
	  if(xtemp >= xpos[0] && xtemp <= xpos[1] && ytemp >= ypos[0] && 
	     ytemp <= ypos[1]) {
	    sd += ((zdata[i] - mean)*(zdata[i]-mean) - sd)/(n+1);
	    ++n;
	  }
	}
	if(n > 1)
	  sd = sqrt(sd*n/(n-1));
	else sd = 0.0f;
	fprintf(stdout, "\n\n\t\tmean\t=\t%g\n\t\tsd\t=\t%g\n\t\tmin\t=\t%g\n\t\tmax\t=\t%g\n\t\tnpts\t=\t%d\n", mean, sd, min, max, n);
	first = 1;
	cpgsci(oldcol);
	break;
	/*
	 * Print information about the nearest point
	 */
      case G_INS:
	/*
	 * Find the nearest grid point
	 */
	{
	  bool first=true;
	  float dist,distmin,xtmp,ytmp,val;
	  int xmin,ymin;
	  for(i=0;i < ndata;i++) {
	    yind = i/nx;
	    xind = i - yind*nx;
	    xtmp = xmins + dx/2 + xind*dx;
	    ytmp = ymins + dy/2 + yind*dy;
	    dist = sqrt((xtmp-xpos[0])*(xtmp-xpos[0]) + 
			(ytmp-ypos[0])*(ytmp-ypos[0]));
	    if(dist < distmin || first) {
	      xmin = xind;
	      ymin = yind;
	      distmin = dist;
	      first = false;
	      val = zdata[i];
	    }
	  }

	  if(coordCallback_) {
	    std::string xstr, ystr;
	    (*coordCallback_)(xpos[0], ypos[0], xstr, ystr, coordCallbackArgs_);
	    fprintf(stdout,"Pixel value at (%s, %s) is %f %s\n", xstr.c_str(), ystr.c_str(), val, unit);
	  } else {
	    fprintf(stdout,"Pixel value at (%f, %f) is %f %s\n", xpos[0], ypos[0], val, unit);
	  }

	}
	break;
      case G_HELP:     /* Print usage info */
	fprintf(stdout,"\nYou requested help by pressing \'%c\'.\n", G_HELP);
	fprintf(stdout,"All cursor positions are entered with \'%c\' key (Left mouse button)\n", G_CUR);
	fprintf(stdout,"\n %c - Select a sub-image to be displayed.\n", G_ZOOM);
	fprintf(stdout," %c - Calculate statistics on a subset of data.\n", G_STAT);
	fprintf(stdout," %c - Redisplay current plot.\n", G_DIS);
	fprintf(stdout," %c - Fiddle contrast & brightness.\n", G_FID);
	fprintf(stdout," %c - Use greyscale.\n", G_GREY);
	fprintf(stdout," %c - Use rainbow colormap.\n", G_RAIN);
	fprintf(stdout," %c - Use heat colormap.\n", G_HEAT);
	fprintf(stdout,"\nTo end this session hit the \'%c\' key (Right mouse button)\n", G_QUIT);
	fprintf(stdout,"\n");
	break;
      default :
	break;
      case G_GREY:
	cmap = grey;
	break;
      case G_RAIN:
	cmap = rain;
	break;
      case G_HEAT:
	cmap = heat; 
	break;
      case G_FID:
	do {
	  contrast = 5.0 * (ypos[0]-ymid)/(ypos[0] < ymid ? (ymin-ymid) : -(ymax-ymid));
	  bright = 0.5 + 1.0 * (fabs(contrast)+1.0)*((xpos[0] - xmax)/(xmin - xmax) - 0.5);
	  cpgctab(cmap->l,cmap->r,cmap->g,cmap->b,cmap->n,contrast,bright);
	  
	  cpgband(B_NORM, 0, xpos[0], ypos[0], &xpos[0], &ypos[0], &key);
	  if(islower((int) key))
	    key = (char) toupper((int) key);
	} while(key == G_FID);
	read = 0;
      }
      if(read) 
	cpgband(B_NORM, 0, xpos[0], ypos[0], &xpos[0], &ypos[0], &key);
      read = 1;
      if(islower((int) key))
	key = (char) toupper((int) key);
    } while(key != G_QUIT);
  }
  /*
   * Else if no cursor, just plot and exit.
   */
  else {
    
    if(!overplot_)
      cpgpage();

    if(vp_) {
      cpgvstd();
    }

    if(wnad_) {
      if(win_)
	cpgswin(0,1,0,1);
      cpgwnad(0,1,0,1); 
    }

    if(win_)
      cpgswin(xmin-dx/2,xmax+dx/2,ymin-dy/2,ymax+dy/2);

    if(wnad_)
      cpgwnad(xmin-dx/2,xmax+dx/2,ymin-dy/2,ymax+dy/2);

    //    cpgbbuf();
    if(!autoscale){
      zmin=z1;
      zmax=z2;
    }

    //    fprintf(stdout,"zmin = %g\tzmax = %g\n",zmin,zmax);
    cpgctab(cmap->l,cmap->r,cmap->g,cmap->b,cmap->n,contrast,bright);

    // Override min/max greyscale with user-specified values

    if(zmin_ != zmax_) {
      zmin = zmin_;
      zmax = zmax_;

      float dz = (zmax - zmin)/ncontour;

      for(unsigned icont=0; icont < ncontour; icont++) {
	contours[icont] = zmin + icont * dz;
      }
    }

    cpgcont(zdata,nx,ny,i1,i2,j1,j2,&contours[0], ncontour,tr);
    cpgsci(1);

    std::string xOpts, yOpts;

    xOpts = (tick_ && xTick_) ? (xTickLabeling_ ? "BCNST" : "BCST") : "BC";
    yOpts = (tick_ && yTick_) ? (yTickLabeling_ ? "BCNST" : "BCST") : "BC";

    cpgbox(xOpts.c_str(), 0.0,0, yOpts.c_str(),0.0,0);
    
    std::string xLab, yLab, tLab;

    xLab = (label_ && xLabel_) ? xlab : "";
    yLab = (label_ && yLabel_) ? ylab : "";
    tLab = title_ ? title : "";
    
    cpglab(xLab.c_str(), yLab.c_str(), tLab.c_str());
	  
    /*
     * Draw a Wedge on the side.
     */
    if(wedge_)
      cpgwedg("RI",0,4,zmax,zmin,unit);
    //    cpgebuf();
  };
  /* 
   * Close PGPLOT only if it was opened in this function 
   */
  if(!wasopen)
    cpgend(); 

  return 0;
}

/*.......................................................................
 * Make a radial plot of an array
 */
int PgUtil::v_radplot(float data[],int nbin, float rmin, float rmax, 
		      float xmin, float xmax, float ymin, float ymax, int nx, int ny)
{
  int i,j,ind,waserr=0;
  float dr,r,min,max;
  float *rxs=NULL,*rys=NULL;
  int *rns=NULL;
  float x,dx = (xmax-xmin)/(nx-1);
  float y,dy = (ymax-ymin)/(ny-1);
  float xmid = (xmax+xmin)/2;
  float ymid = (ymax+ymin)/2;

  waserr = (rns = (int *)malloc(nbin*sizeof(int)))==NULL;
  waserr |= (rxs = (float *)malloc(nbin*sizeof(float)))==NULL;
  waserr |= (rys = (float *)malloc(nbin*sizeof(float)))==NULL;

  if(!waserr) {
    /*
     * Use nbin bins between rmin and rmax
     */
    dr = (rmax-rmin)/(nbin-1);

    for(i=0;i < nbin;i++) {
      rns[i] = 0;
      rxs[i] = rmin + dr*i;
      rys[i] = 0.0;
    }

    for(i=0;i < nx;i++) {
      for(j=0;j < ny;j++) {
	x = xmin + i*dx;
	y = ymin + j*dy;
	r = sqrt((x-xmid)*(x-xmid)+(y-ymid)*(y-ymid));
	if(r >= rmin && r <= rmax) {
	  ind = (int)floor((r-rmin)/dr);
	  /*
	   * Keep a running mean for each bin.
	   */
	  rys[ind] += (data[i+j*nx] - rys[ind])/(rns[ind]+1);
	  ++rns[ind];
	}
      }
    }
    min = max = rys[0];
    for(i=0;i < nbin;i++) {
      max = MAX(max,rys[i]);
      min = MIN(min,rys[i]);
    }
    cpgask(0);
    cpgpage();
    cpgvstd();
    //    cpgbbuf();
    cpgsci(1);
    {
      float range=rmax-rmin;
      cpgswin(rmin-0.1*range,rmax+0.1*range,min - (max-min)*0.1,max+(max-min)*0.1);
    }
    cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
    cpgsci(10);
    cpgsls(1);
    cpgline(nbin, rxs, rys);
    //    cpgebuf(); 

  }
  /*
   * Free any allocated memory,
   */
  if(rns)
    free(rns);
  if(rxs)
    free(rxs);
  if(rys)
    free(rys);
  return waserr;
}

/*.......................................................................
 * Plot a simple line plot.
 *
 * Input:
 *  ndata   int     The number of data points.
 *  data  float  *  The array of values to plot.
 *  xmin  float     The minimum xvalu.e
 *  xmax  float     The maximum xvalue.
 *  nbin    int     The number of bins to use.
 */
void PgUtil::linePlot(std::vector<float>& yfarr, bool doLine)
{
  std::vector<double> ydarr;
  ydarr.resize(yfarr.size());

  for(unsigned i=0; i < yfarr.size(); i++) {
    ydarr[i] = yfarr[i];
  }
  
  return linePlot(ydarr, doLine);
}

void PgUtil::linePlot(std::vector<unsigned>& yfarr, bool doLine)
{
  std::vector<double> ydarr;
  ydarr.resize(yfarr.size());

  for(unsigned i=0; i < yfarr.size(); i++) {
    ydarr[i] = yfarr[i];
  }
  
  return linePlot(ydarr, doLine);
}

void PgUtil::plotPowerSpectrum(std::vector<float>& yarr, bool doLine)
{
  std::vector<double> dyarr(yarr.size());

  for(unsigned i=0; i < yarr.size(); i++)
    dyarr[i] = yarr[i];

  return plotPowerSpectrum(dyarr, doLine);
}

void PgUtil::plotPowerSpectrum(std::vector<double>& yarr, bool doLine)
{
  Dft1d dft(yarr.size(), true);

  for(unsigned i=0; i < yarr.size(); i++)
    dft.pushSample(yarr[i]);

  std::vector<double> y = dft.abs();
  std::vector<double> xp(y.size()-1);
  std::vector<double> yp(y.size()-1);

  double dx = 1.0;
  unsigned n = xp.size();

  for(unsigned i=0; i < n; i++) {
    xp[i] = (i+1)*2*M_PI/n;
    yp[i] = y[i+1];
  }

  setLogPlot(true);
  linePlot(xp, yp, "", "", "", doLine);
}

void PgUtil::linePlot(std::vector<double>& yarr, bool doLine)
{
  std::vector<double> xarr;
  xarr.resize(yarr.size());

  for(unsigned i=0; i < yarr.size(); i++) {
    xarr[i] = i;
  }

  return linePlot(xarr, yarr, "", "", "", doLine);
}

void PgUtil::errPlot(std::vector<double>& xarr, std::vector<double>& yarr, std::vector<double>& earr,
		     std::string xlab, std::string ylab, std::string title)
{
  return linePlot(xarr.size(), &xarr[0], &yarr[0], &earr[0],
		  xlab, ylab, title, false, true);
}

void PgUtil::errPlot(std::valarray<double>& xarr, std::valarray<double>& yarr, std::valarray<double>& earr,
		     std::string xlab, std::string ylab, std::string title)
{
  return linePlot(xarr.size(), &xarr[0], &yarr[0], &earr[0],
		  xlab, ylab, title, false, true);
}

void PgUtil::linePlot(std::vector<float>& xarr, std::vector<float>& yarr,
		      std::string xlab, std::string ylab, std::string title, bool doLine) 
{
  std::vector<double> xd(xarr.size());
  std::vector<double> yd(yarr.size());

  for(unsigned i=0; i < xd.size(); i++) {
    xd[i] = xarr[i];
    yd[i] = yarr[i];
  }

  return linePlot(xd, yd, xlab, ylab, title, doLine);
}

void PgUtil::linePlot(std::vector<double>& xarr, std::vector<double>& yarr,
		      std::string xlab, std::string ylab, std::string title, bool doLine) 
{
  return linePlot(xarr.size(), &xarr[0], &yarr[0], 0,
		  xlab, ylab, title, doLine, false);
}

void PgUtil::linePlot(int narr, double* xarr, double* yarr, double* earr,
		      std::string xlab, std::string ylab, std::string title, bool doLine, bool doErr) 
{
  std::vector<float> fxarr(narr);
  std::vector<float> fyarr(narr);
  std::vector<float> fearr(narr);

  for(unsigned i=0; i < narr; i++) {
    fxarr[i] = xarr[i];
    fyarr[i] = yarr[i];
    fearr[i] = earr ? earr[i] : 0.0;
  }

  return linePlot(narr, &fxarr[0], &fyarr[0], &fearr[0],
		  xlab, ylab, title, doLine, (earr ? doErr : false));
}

void PgUtil::linePlot(int narr, float* xarr, float* yarr, float* earr,
		      std::string xlab, std::string ylab, std::string title, bool doLine, bool doErr) 
{
#if 1
  PgUtil pgUtil;
  LinePlot lp(&pgUtil, narr, xarr, yarr, earr, xlab, ylab, title, doLine, doErr);
  lp.plot();
#else
  v_lplot(narr, xarr, yarr, earr, xlab, ylab, title, doLine, doErr);
#endif
}

int PgUtil::v_lplot(int narr, float* xarr, float* yarr, float* earr,
		    std::string xlab, std::string ylab, std::string title, bool doLine, bool doErr) 
{
  int i;
  char value[20];
  int length;
  float xmins,xmaxs,ymins,ymaxs; /* Limits for the whole plot */
  float xmin,xmax,ymin,ymax;     /* temporary limits (zoomed, etc.) */
  char answer[10];
  int slen,ix,iy,indx,indy,index;
  int i1,i2,j1,j2;
  float xtemp,ytemp,xpos[2],ypos[2];
  static string mess1("\n %c - Select start and end vertices using this key\n");
  static string mess2(" %c - Select the whole plot\n");
  static string mess3(" %c - Abort selection\n");
  int cancel,dofull,accepted;
  int iter;
  char key = 'L';
  bool docurs = false;;
  bool wasopen = false;
  int read=1;
  int oldcol;
  int redisp;
  double ylo, yhi;

  std::vector<float> xs(narr);
  std::vector<float> ys(narr);

  xmins = xmaxs = xarr[0];
  ymins = (earr ? yarr[0] - earr[0] : yarr[i]);
  ymaxs = (earr ? yarr[0] - earr[0] : yarr[i]);

  for(i=0;i < narr;i++) {
    xmins = (xmins < xarr[i]) ? xmins : xarr[i];
    xmaxs = (xmaxs > xarr[i]) ? xmaxs : xarr[i];
    ylo = earr ? yarr[i]-earr[i] : yarr[i];
    yhi = earr ? yarr[i]+earr[i] : yarr[i];
    ymins = (ymins < ylo) ? ymins : ylo;
    ymaxs = (ymaxs > yhi) ? ymaxs : yhi;
  }

#if 1
  if(logPlot_) {
    xmins = log10(xmins);
    xmaxs = log10(xmaxs);
    ymins = log10(ymins);
    ymaxs = log10(ymaxs);
  }
#endif

  xmins -= 0.1*(xmaxs-xmins);
  xmaxs += 0.1*(xmaxs-xmins);

  ymins -= 0.1*(ymaxs-ymins);
  ymaxs += 0.1*(ymaxs-ymins);

  if(ymins == ymaxs) {
    ymins -= 0.1*ymins;
    ymaxs += 0.1*ymaxs;
  }

  if(usedefs_) {

    if(xmin_ != xmax_) {
      xmins = logPlot_ ? log10(xmin_) : xmin_;
      xmaxs = logPlot_ ? log10(xmax_) : xmax_;
    }

    if(ymin_ != ymax_) {
      ymins = logPlot_ ? log10(ymin_) : ymin_;
      ymaxs = logPlot_ ? log10(ymax_) : ymax_;
    }
  }

  for(unsigned i=0; i < narr; i++) {
    if(logPlot_) {
      xs[i] = log10(xarr[i]);
      ys[i] = log10(yarr[i]);
    } else {
      xs[i] = xarr[i];
      ys[i] = yarr[i];
    }
  }

  // If no PGPLOT device has been opened, open one.

  queryDevice(wasopen);

  // Do we have a cursor?

  docurs = haveCursor();

  xmin = xmins;
  xmax = xmaxs;
  ymin = ymins;
  ymax = ymaxs;

  cpgsch(ch_);

  unsigned cursor;
  if(docurs) {
    cursor = B_NORM;
    printf("For HELP, hit the \'%c\' key on your keyboard\n", KEY_HELP);
    do {
      redisp = 0;
      cancel = 0;

      switch(key) {
      case KEY_ZOOM:
	v_lzoom(xmins, xmaxs, ymins, ymaxs, &xmin, &xmax, &ymin, &ymax);
	redisp = 1;
	cursor = B_NORM;
	break;
      case KEY_UT:
	redisp=v_lnewx(xmins,xmaxs,ymins,ymaxs,&xmin,&xmax,&ymin,&ymax,narr,&xs[0],&ys[0]);
	break;
      case KEY_FLG:
	v_lwrite(narr, xarr, yarr, xmin, xmax, ymin,ymax);
	redisp = 0;
	break;
      case KEY_NXT:
	v_lnum(narr, xarr, yarr, xmin, xmax, ymin,ymax);
	redisp = 0;
	break;
      case KEY_BEAM: /* Label the ten highest peaks */
	v_lten(narr, &xs[0], &ys[0], xmin, xmax, ymin,ymax);
	redisp = 0;
	break;
      case KEY_DIS:
	redisp = 1;
	break;
      case KEY_CROS:
	if(cursor == B_CROSS)
	  cursor = B_NORM;
	else
	  cursor = B_CROSS;
	break;
      case KEY_HELP:     /* Print usage info */
	printf("\nYou requested help by pressing \'%c\'.\n", KEY_HELP);
	printf("All cursor positions are entered with \'%c\' key (Left mouse button)\n", KEY_CUR);
	printf(" %c - Fit a gaussian/polynomial to selected data.\n", KEY_FIT);
	printf(" %c - Label highest peak\n", KEY_FLG);
	printf(" %c - Label highest peak with just a number\n", KEY_NXT);
	printf(" %c - Label highest peaks over a range\n", KEY_BEAM);
	printf(" %c - Select X range to be displayed (hit %c twice for full range)\n", KEY_UT, KEY_UT);
	printf(" %c - Select a sub-image to be displayed.\n", KEY_ZOOM);
	printf(" %c - Redisplay current plot.\n", KEY_DIS);
	printf("\nTo end this session hit the \'%c\' key (Right mouse button)\n", KEY_QUIT);
	printf("\n");
	break;
      default :
	break;
      }

      if(redisp)
      	if(!cancel) 
	  v_ldraw(narr,&xs[0],&ys[0],earr,xlab,ylab,title,doLine,doErr,xmin, xmax, ymin, ymax);  

      if(read) {
	cpgband(cursor, 0, xpos[0], ypos[0], &xpos[0], &ypos[0], &key);
      }

      read = 1;

      if(islower((int) key))
	key = (char) toupper((int) key);

    } while(key != KEY_QUIT);
  }
  else 
    v_ldraw(narr,&xs[0],&ys[0],earr,xlab,ylab,title,doLine,doErr, xmin, xmax, ymin, ymax);  
  return 0;
}

/*.......................................................................
 * Module for v_lplot.
 */
int PgUtil::v_lnewx(float xmins, float xmaxs, float ymins, float ymaxs, 
		    float *xmin, float *xmax, float *ymin, float *ymax, 
		    int narr,float xarr[], float yarr[])
{
  int iter;
  int oldcol;
  int dofull;
  int cancel;
  int accepted;
  float xpos[2],ypos[2];
  float xtemp,ytemp;
  char key;
  static string mess1("\n %c - Select start and end vertices using this key\n");
  static string mess2(" %c - Select the whole plot\n");
  static string mess3(" %c - Abort selection\n");

  cpgqci(&oldcol);
  cpgsci(5);

  dofull = 0;
  cancel = 0;
  for(iter = 0;iter<2 && !dofull && !cancel;iter++) {
    do {
      accepted = 0;
      cpgband((iter==0) ? B_XVAL : B_XRNG, 0, xtemp, ytemp, &xtemp,
	      &ytemp, &key);
      if(islower((int) key))
	key = (char) toupper((int) key);
      xpos[iter] = xtemp;
      ypos[iter] = ytemp;
      switch(key) {
      case KEY_UT:
	accepted = dofull = 1;
	break;
      case KEY_QUIT: case KEY_CAN:      /* Abort box selection */
	return 0;
	break;
      case KEY_CUR:             /* Accept the selected start vertex */
	accepted=1;
	break;
      default:            /* Unexpected cursor input key - show usage */
	printf(mess1.c_str(), KEY_CUR);
	printf(mess2.c_str(), KEY_UT);
	printf(mess3.c_str(), KEY_CAN);
	break;
      };
    } while(!accepted);
  };
  if(dofull) {
    *xmin = xmins;
    *ymin = ymins;
    *xmax = xmaxs;
    *ymax = ymaxs;
  }
  else {
    *xmin = (xpos[0] < xpos[1]) ? xpos[0] : xpos[1];
    *xmax = (xpos[0] > xpos[1]) ? xpos[0] : xpos[1];
  }
  {
    int first,i;

    for(i=0;i < narr;i++) 
      if(xarr[i] < (*xmax) && xarr[i] > (*xmin)) {
	if(first) {
	  *ymin=*ymax=yarr[i];
	  first=0;
	}
	*ymin = MIN(*ymin, yarr[i]);
	*ymax = MAX(*ymax, yarr[i]);
      }
  }
  *ymax += 0.1*(*ymax-*ymin);
  *ymin -= 0.1*(*ymax-*ymin);
  cpgswin(*xmin,*xmax,*ymin,*ymax);
  cpgsci(oldcol);

  return 1;
}
/*.......................................................................
 * Module for v_lplot.
 */
int PgUtil::v_lwrite(int narr, float xarr[], float yarr[], float xmin,
		     float xmax, float ymin, float ymax)
{
  int i,imax,iter;
  int oldcol;
  int dofull;
  int cancel;
  int accepted;
  float xpos[2],ypos[2];
  float x1,x2,y1,y2;
  float xtemp,ytemp;
  float max;
  char num[10];
  char key;
  int first=1;
  static string mess1("\n %c - Select start and end vertices using this key\n");
  static string mess2(" %c - Select the whole plot\n");
  static string mess3(" %c - Abort selection\n");

  cpgqci(&oldcol);
  cpgsci(5);

  dofull = 0;
  cancel = 0;
  for(iter = 0;iter<2 && !dofull && !cancel;iter++) {
    do {
      accepted = 0;
      cpgband((iter==0) ? B_XVAL : B_XRNG, 0, xtemp, ytemp, &xtemp,
	      &ytemp, &key);
      if(islower((int) key))
	key = (char) toupper((int) key);
      xpos[iter] = xtemp;
      ypos[iter] = ytemp;
      switch(key) {
      case KEY_FLG:
	accepted = dofull = 1;
	break;
      case KEY_CAN:      /* Abort box selection */
	accepted = cancel = 1;
	break;
      case KEY_QUIT:     /* Quit now */
	cpgend();
	return 0;
	break;
      case KEY_CUR:             /* Accept the selected start vertex */
	accepted=1;
	break;
      default:            /* Unexpected cursor input key - show usage */
	printf(mess1.c_str(), KEY_CUR);
	printf(mess2.c_str(), KEY_FLG);
	printf(mess3.c_str(), KEY_CAN);
	break;
      };
    } while(!accepted);
  };
  if(dofull) {
    x1 = xmin;
    x2 = xmax;
  }
  else {
    x1 = (xpos[0] < xpos[1]) ? xpos[0] : xpos[1];
    x2 = (xpos[0] > xpos[1]) ? xpos[0] : xpos[1];
  }
  for(i=0;i < narr;i++) {
    if(xarr[i] > x1 && xarr[i] < x2) {
      if(first) {
	imax = i;
	max = yarr[i];
	first = 0;
      }
      if(yarr[i] > max) {
	imax = i;
	max = yarr[i];
      }
    }
  }
  /*
   * Draw a line on the plot.
   */
  y1 = yarr[imax]+(ymax-ymin)/32;
  y2 = yarr[imax]+5*(ymax-ymin)/32;
  cpgmove(xarr[imax], y1);
  cpgdraw(xarr[imax], y2);
  sprintf(num,"%.3g",xarr[imax]);
  cpgptxt(xarr[imax],y2+(ymax-ymin)/64,0,0.5,num);

  cpgsci(oldcol);
  return 1;
}
/*.......................................................................
 * Module for v_lplot.
 */
int PgUtil::v_lnum(int narr, float xarr[], float yarr[], float xmin,
		   float xmax, float ymin, float ymax)
{
  int i,imax,iter;
  int oldcol;
  int dofull;
  int cancel;
  int accepted;
  float xpos[2],ypos[2];
  float x1,x2,y1,y2;
  float xtemp,ytemp;
  float max;
  char num[10];
  char key;
  int first=1;
  static string mess1("\n %c - Select start and end vertices using this key\n");
  static string mess2(" %c - Select the whole plot\n");
  static string mess3(" %c - Abort selection\n");

  cpgqci(&oldcol);
  cpgsci(5);

  dofull = 0;
  cancel = 0;
  for(iter = 0;iter<2 && !dofull && !cancel;iter++) {
    do {
      accepted = 0;
      cpgband((iter==0) ? B_XVAL : B_XRNG, 0, xtemp, ytemp, &xtemp,
	      &ytemp, &key);
      if(islower((int) key))
	key = (char) toupper((int) key);
      xpos[iter] = xtemp;
      ypos[iter] = ytemp;
      switch(key) {
      case KEY_NXT:
	accepted = dofull = 1;
	break;
      case KEY_CAN:      /* Abort box selection */
	accepted = cancel = 1;
	break;
      case KEY_QUIT:     /* Quit now */
	cpgend();
	return 0;
	break;
      case KEY_CUR:             /* Accept the selected start vertex */
	accepted=1;
	break;
      default:            /* Unexpected cursor input key - show usage */
	printf(mess1.c_str(), KEY_CUR);
	printf(mess2.c_str(), KEY_FLG);
	printf(mess3.c_str(), KEY_CAN);
	break;
      };
    } while(!accepted);
  };
  if(dofull) {
    x1 = xmin;
    x2 = xmax;
  }
  else {
    x1 = (xpos[0] < xpos[1]) ? xpos[0] : xpos[1];
    x2 = (xpos[0] > xpos[1]) ? xpos[0] : xpos[1];
  }
  for(i=0;i < narr;i++) {
    if(xarr[i] > x1 && xarr[i] < x2) {
      if(first) {
	imax = i;
	max = yarr[i];
	first = 0;
      }
      if(yarr[i] > max) {
	imax = i;
	max = yarr[i];
      }
    }
  }
  /*
   * Draw the number on the plot.
   */
  y1 = yarr[imax]+(ymax-ymin)/32;
  y2 = yarr[imax]+5*(ymax-ymin)/32;
  sprintf(num,"%.3g",xarr[imax]);
  cpgptxt(xarr[imax],y1,0,0.5,num);

  cpgsci(oldcol);
  return 1;
}
/*.......................................................................
 * Module for v_lplot.
 */
int PgUtil::v_lten(int narr, float xarr[], float yarr[], float xmin,
		   float xmax, float ymin, float ymax)
{
  int imax,iter;
  int oldcol;
  int dofull;
  int cancel;
  int accepted;
  float xpos[2],ypos[2];
  float x1,x2,y1,y2;
  float xtemp,ytemp;
  float max;
  char num[10];
  char key;
  int first=1;
  static string mess1("\n %c - Select start and end vertices using this key\n");
  static string mess2(" %c - Select the whole plot\n");
  static string mess3(" %c - Abort selection\n");
  int *index=NULL;
  float xval,yval;
  int n=10;
  char nstring[5];
  int keymax=5;
  char *ptr=NULL;
  int i=0,j=0;

  cpgqci(&oldcol);
  cpgsci(5);
  /*
   * Get a numerical string.
   */
  do {
    accepted = 0;
    cpgband(B_NORM, 0, xpos[0], ypos[0], &xpos[0], &ypos[0], &key);
    switch(key) {
    case KEY_QUIT: case KEY_CAN:     /* Abort box selection */
      return 0;
      break;
    case KEY_HELP:
      fprintf(stdout,"\nEnter a number, followed by a carriage return.\n");
      break;
    case '\n': case '\r':
      nstring[i] = '\0';
      break;
    default:
      nstring[i] = key;
      ++i;
      break;
    } 
  } while(!(key == '\n' || key == '\r') && i < keymax);;
  /* 
   * If a number
   */
  n = (int)strtod(nstring,&ptr);
  if(*ptr != '\0') {
    fprintf(stderr,"Not a number.\n");
    return 0;
  }
  /*
   * Now select the range over which we want to flag peaks
   */
  dofull = 0;
  cancel = 0;
  for(iter = 0;iter<2 && !dofull && !cancel;iter++) {
    do {
      accepted = 0;
      cpgband((iter==0) ? B_XVAL : B_XRNG, 0, xtemp, ytemp, &xtemp,
	      &ytemp, &key);
      if(islower((int) key))
	key = (char) toupper((int) key);
      xpos[iter] = xtemp;
      ypos[iter] = ytemp;
      switch(key) {
      case KEY_NXT:
	accepted = dofull = 1;
	break;
      case KEY_CAN:      /* Abort box selection */
	accepted = cancel = 1;
	break;
      case KEY_QUIT:     /* Quit now */
	cpgend();
	return 0;
	break;
      case KEY_CUR:             /* Accept the selected start vertex */
	accepted=1;
	break;
      default:            /* Unexpected cursor input key - show usage */
	printf(mess1.c_str(), KEY_CUR);
	printf(mess2.c_str(), KEY_FLG);
	printf(mess3.c_str(), KEY_CAN);
	break;
      };
    } while(!accepted);
  };
  if(dofull) {
    x1 = xmin;
    x2 = xmax;
  }
  else {
    x1 = (xpos[0] < xpos[1]) ? xpos[0] : xpos[1];
    x2 = (xpos[0] > xpos[1]) ? xpos[0] : xpos[1];
  }
  /*
   * Sort the arrays
   */
  indexx(narr, yarr, &index);
  for(j=0,i=narr-1;i > 0;i--) {
    xval = xarr[index[i]];
    yval = yarr[index[i]];
    if(j < n && (xval >= x1 && xval <= x2)) {
      /*
       * Draw the number on the plot.
       */
      y1 = yval +(ymax-ymin)/32;
      sprintf(num,"%.3g",xval);
      cpgptxt(xval,y1,0,0.5,num);
      j++;
    }
  }
  cpgsci(oldcol);
  if(index)
    free(index);
  return 1;
}

/**.......................................................................
 * Module for v_hist.
 */
int PgUtil::v_lzoom(float xmins, float xmaxs, float ymins, float ymaxs, 
		    float *xmin, float *xmax, float *ymin, float *ymax)
{
  int iter;
  int oldcol;
  int dofull;
  int cancel;
  int accepted;
  float xpos[2],ypos[2];
  float xtemp,ytemp;
  char key;
  static string mess1("\n %c - Select start and end vertices using this key\n");
  static string mess2(" %c - Select the whole plot\n");
  static string mess3(" %c - Abort selection\n");

  cpgqci(&oldcol);
  cpgsci(5);

  dofull = 0;
  cancel = 0;
  for(iter = 0;iter<2 && !dofull && !cancel;iter++) {
    do {
      accepted = 0;
      cpgband((iter==0) ? B_NORM : B_RECT, 0, xtemp, ytemp, &xtemp,
	      &ytemp, &key);
      if(islower((int) key))
	key = (char) toupper((int) key);
      xpos[iter] = xtemp;
      ypos[iter] = ytemp;
      switch(key) {
      case KEY_ZOOM:
	accepted = dofull = 1;
	break;
      case KEY_CAN:      /* Abort box selection */
	accepted = cancel = 1;
	break;
      case KEY_QUIT:     /* Quit now */
	cpgend();
	return 0;
	break;
      case KEY_CUR:             /* Accept the selected start vertex */
	accepted=1;
	break;
      default:            /* Unexpected cursor input key - show usage */
	printf(mess1.c_str(), KEY_CUR);
	printf(mess2.c_str(), KEY_ZOOM);
	printf(mess3.c_str(), KEY_CAN);
	break;
      };
    } while(!accepted);
  };
  if(dofull) {
    *xmin = xmins;
    *ymin = ymins;
    *xmax = xmaxs;
    *ymax = ymaxs;
  }
  else {
    *xmin = (xpos[0] < xpos[1]) ? xpos[0] : xpos[1];
    *ymin = (ypos[0] < ypos[1]) ? ypos[0] : ypos[1];
    *xmax = (xpos[0] > xpos[1]) ? xpos[0] : xpos[1];
    *ymax = (ypos[0] > ypos[1]) ? ypos[0] : ypos[1];
  }
  cpgswin(*xmin,*xmax,*ymin,*ymax);
  cpgsci(oldcol);

  return 1;
}
/*.......................................................................
 * Module for v_lplot -- Draw the line plot.
 *
 */
int PgUtil::v_ldraw(int narr, float xarr[], float yarr[], float earr[],
		    std::string xlab, std::string ylab, std::string title,
		    bool doLine, bool doErr,
		    double xmin, double xmax, double ymin, double ymax)
{
  int i;
  float x1,x2;
  float ypts[5],xpts[5];

  if(!overplot_)
    cpgpage();

  if(vp_) {
    cpgvstd();
  }

  if(win_)
    cpgswin(xmin,xmax,ymin,ymax);
  
  if(wnad_) {
    cpgwnad(xmin,xmax,ymin,ymax);
  }


  //  cpgbbuf();

  cpgsci(1);

  //------------------------------------------------------------
  // Label the plot
  //------------------------------------------------------------

  std::string xOpts, yOpts;

  if(logPlot_) {
    xOpts = (tick_ && xTick_) ? (xTickLabeling_ ? (xTickLabelAtBottom_ ? "BCNLST" : "BCMLST") : "BCLST")  : "BLC";
  } else {
    xOpts = (tick_ && xTick_) ? (xTickLabeling_ ? (xTickLabelAtBottom_ ? "BCNST" : "BCMST") : "BCST")   : "BC";
  }

  if(logPlot_) {
    yOpts = (tick_ && yTick_) ? (yTickLabeling_ ? (yTickLabelAtLeft_ ? "BCNLST" : "BCMLST") : "BCLST") : "BLC";
  } else {
    yOpts = (tick_ && yTick_) ? (yTickLabeling_ ? (yTickLabelAtLeft_ ? "BCNST" : "BCMST") : "BCST")  : "BC";
  }

  // Draw the box only after we're done, so the histogram envelope
  // lines don't overlap it

  if(useBoxCi_)
    cpgsci(boxCi_);
  else
    cpgsci(1);

  cpgbox(xOpts.c_str(), 0.0,0, yOpts.c_str(), 0.0,0);
  
  std::string xLab, yLab, tLab;
  
  xLab = (label_ && xLabel_) ? xlab : "";
  yLab = (label_ && yLabel_) ? ylab : "";
  tLab = title_ ? title : "";
  
  cpglab(xLab.c_str(), yLab.c_str(), tLab.c_str());

  if(useTraceCi_)
    cpgsci(traceCi_);
  else
    cpgsci(10);

  if(doLine) {
    cpgline(narr, xarr, yarr);
  } else {

    if(doErr) {
      std::vector<float> y1, y2;
      y1.resize(narr);
      y2.resize(narr);
      for(unsigned i=0; i < narr;i++) {
	y1[i] = yarr[i] - earr[i];
	y2[i] = yarr[i] + earr[i];
      }

      cpgerry(narr, xarr, &y1[0], &y2[0], 1.0);

    } else {
      cpgpt(narr, xarr, yarr, -1);
    }
  }

  // And reset the color to the box color

  if(useBoxCi_)
    cpgsci(boxCi_);
  else
    cpgsci(1);
  
  return 1;
}

void PgUtil::setColormap(std::string cmap)
{
  for(unsigned i=0; i < n_std_cmap; i++) {
    if(strcmp(cmap.c_str(), std_cmaps[i].name.c_str())==0) {
      cmap_ = &std_cmaps[i];
      return;
    }
  }

  ThrowSimpleColorError("Unrecognized colormap name: '" << cmap << "'.  Should be one of 'grey', 'heat', or 'rainbow'", "red");
}

/*.......................................................................
  This function takes a data array, arrin[] with npts elements and
  returns, via the argument list, an index array *indx. The sort method
  is the heap-sort, an NlogN process.  IMPORTANT: The indx arrays is
  allocated in this function, so the calling routine MUST deallocate it
  after use. On error in this routine, -1 is returned and the indx
  array is automatically zapped.
*/
int PgUtil::indexx(int npts, float arrin[], int **indx)
{
  int indxt,i;
  /*
    Allocate memory for the index array.
  */
  if( (*indx = (int *) calloc(npts+1, sizeof(int))) == NULL) {
    fprintf(stderr, "sort: Memory allocation of index array failed.\n");
    return -1;
  };
  /*
    Initialize the index array elements with their element numbers.
    The data array will be indexed through this array for the
    comparisons between data array elements.
  */

  for (i=0; i<npts; i++) (*indx)[i]=i;
  /*
    The algorithm fails for npts=1. However no sorting is required in this case
    anyway so return with no error.
  */
  if(npts==1)
    return 0;
  /*
    Re-arrange the index array to index the data array as a heap.
    Start at the bottom node of the heap tree, i=npts/2-1 and
    work upwards, rearanging the values in the tree until all
    branches have lower values than there nodes.
  */
  for (i = npts/2 - 1 ; i >= 0; i--)
    insert_in_heap(arrin, *indx, i, npts, (*indx)[i]);
  /*
    Now parse the tree, from the top node downwards, removing the
    value of the top node, and recursively promoting values in sub-nodes to
    fill the gap. The removed value is placed at the top of the index array,
    in the element vacated by the last promotion. 'i' keeps a record
    of the number of nodes still to be removed.
  */
  for(i = npts-1; i > 0; i--) {
    /*
      Remove the root value and copy to its resting place at the end
      of the un-treated portion of the array. This overwrites the last
      element, so keep a temporary record of it in indxt.
    */
    indxt = (*indx)[i];
    (*indx)[i] = (*indx)[0];
    /*
      Now follow down the branches, promoting the highest valued branches.
    */
    insert_in_heap(arrin, *indx, 0, i, indxt);
  };
  return 0;
}

/*.......................................................................
  This function is called by indexx. Given the element number, new_el of
  the data array, arrin[], search arrin[] as a heap, from the node at
  element number, node, to find the correct position for the new
*/
void PgUtil::insert_in_heap(float arrin[], int indx[], int node, int num_node, int new_el)
{
  int branch, right;
  float temp_value;
  /*
    node holds a record of the current node being treated. At the start
    of each iteration, branch points to the leftmost branch of that node.
  */
  branch = node + node + 1;
  /*
    Keep a record of the value pointed to by new_el in the data array.
  */
  temp_value = arrin[new_el];
  /*
    Follow the tree branches, promoting branch values where necessary,
    until the appropriate node is found for the index to the new value
    recorded in new_value.
  */
  while (branch < num_node) {
    /*
      Make 'branch' point at the branch with the highest value in it.
      Initially it points at the left branch.
    */
    right = branch+1;
    if (right < num_node && arrin[indx[branch]] < arrin[indx[right]]) 
      branch = right;
    /*
      Stop looking when the root value exceeds the values of both branches.
    */
    if(temp_value >= arrin[indx[branch]])
      break;
    /*
      Otherwise promote the highest branch value and continue the search
      down that branch.
    */
    indx[node] = indx[branch];
    node = branch;
    branch += branch+1;
  };
  /*
    Install the index of the temp_value at the newly found location in
    the index array.
  */
  indx[node]=new_el;
  return;
}


int PgUtil::open(std::string device)
{
  int pgid = cpgopen(device.c_str());
  cpgslct(pgid);
  return pgid;
}

void PgUtil::close()
{
  cpgclos();
}

void PgUtil::subplot(int nx, int ny)
{
  cpgsubp(nx, ny);
}

void PgUtil::advance()
{
  cpgpage();
}

bool PgUtil::haveOpenDevice()
{
  char answer[10];
  int slen = sizeof(answer)-1;
  cpgqinf("STATE", answer, &slen);
  bool wasopen = strncmp(answer,"OPEN",4)==0;

  return wasopen;
}

void PgUtil::queryDevice(bool& wasopen)
{
  wasopen = haveOpenDevice();
  if(!wasopen) {
    if(cpgbeg(0,"?",1,1)!=1) {
      ThrowError("Error in cpgbeg()");
    }
    cpgask(0);
  }
}

bool PgUtil::haveCursor()
{
  if(!interactive_)
    return false;

  int slen;
  char answer[10];
  
  // Do we have a cursor?

  slen = sizeof(answer)-1;
  cpgqinf("CURSOR", answer, &slen);
  return (strncmp(answer,"YES",3) == 0);
}

void PgUtil::wnad(float xmin, float xmax, float ymin, float ymax)
{
  cpgwnad(xmin, xmax, ymin, ymax);
}

void PgUtil::binPlot(std::vector<double>& x, std::vector<double>& y, unsigned nbin)
{
  std::vector<float> hist(nbin);
  double min, max, mid;

  for(unsigned i=0; i < nbin; i++)
    hist[i] = 0.0;
  
  // Get the min and max of the data.

  for(unsigned i=0; i < x.size(); i++) {

    if(i==0) 
      min = max = x[i];

    min = (x[i] < min ? x[i] : min);
    max = (x[i] > max ? x[i] : max);
  }

  // And construct the histogram.

  double dx = (max-min)/(nbin-1);
  mid = (max + mid)/2;

  if(fabs(dx/mid) < 1e-12) {
    min = 0.9*min;
    max = 1.1*max;
    dx = (max-min)/(nbin-1);
  }
    
  for(unsigned i=0; i < x.size(); i++) {
    int ind = (int)floor((x[i]-min)/dx + 0.5);
    hist[ind] += y[i];
  }
  
  v_hist(&hist[0], min, max, nbin, "");
}

void PgUtil::binData(std::vector<double>& x, std::vector<double>& y, unsigned nbin,
		     std::vector<double>& xret, std::vector<double>& yret)
{
  xret.resize(nbin);
  yret.resize(nbin);
  double min, max;

  // Get the min and max of the data.

  for(unsigned i=0; i < x.size(); i++) {

    if(i==0) 
      min = max = x[i];

    min = (x[i] < min ? x[i] : min);
    max = (x[i] > max ? x[i] : max);
  }

  // And construct the histogram.

  double dx = (max-min)/(nbin-1);

  for(unsigned i=0; i < nbin; i++) {
    xret[i] = min + dx*i;
    yret[i] = 0.0;
  }

  for(unsigned i=0; i < x.size(); i++) {
    int ind = (int)floor((x[i]-min)/dx + 0.5);
    yret[ind] += y[i];
  }
}

void PgUtil::histogram(std::vector<double>* dptr, unsigned nbin, std::vector<unsigned>* multiplicity)
{
  std::vector<float> fvec(dptr->size());

  for(unsigned i=0; i < fvec.size(); i++) {
    fvec[i] = dptr->at(i);
  }

  return histogram(fvec.size(), &fvec[0], nbin, multiplicity);
}

void PgUtil::histogram2D(std::vector<double>* dptr1, std::vector<double>* dptr2, unsigned nbin, std::vector<unsigned>* multiplicity)
{
  std::vector<float> fvec1(dptr1->size());
  std::vector<float> fvec2(dptr2->size());

  for(unsigned i=0; i < fvec1.size(); i++) {
    fvec1[i] = dptr1->at(i);
    fvec2[i] = dptr2->at(i);
  }

  return histogram2D(fvec1.size(), &fvec1[0], &fvec2[0], nbin, nbin, multiplicity);
}

void PgUtil::histogram(unsigned ndata, float* data, unsigned nbin, std::vector<unsigned>* multiplicity)
{
  int waserr = 0;
  int first = 1;
  int i,ind;
  double min,max,mid,dx;

  try {

    std::vector<float> hist(nbin);

    for(i=0;i < nbin;i++)
      hist[i] = 0.0;
  
    // Get the min and max of the data.

    if(!waserr) {
      for(i=0;i < ndata;i++) {

	if(!finite(data[i])) {
	  ThrowColorError("\nData not finite", "red");
	}

	if(first) {
	  min = max = data[i];
	  first = 0;
	}
	min = (data[i] < min ? data[i] : min);
	max = (data[i] > max ? data[i] : max);
      }
    }

    // And construct the histogram.

    dx = (max-min)/(nbin-1);
    mid = (max + min)/2;

    if(fabs(dx/mid) < 1e-12) {
      min = 0.9*min;
      max = 1.1*max;
      dx = (max-min)/(nbin-1);
    }

    for(i=0;i < ndata;i++) {
      ind = (int)floor((data[i]-min)/dx + 0.5);
      if(ind < 0 || ind > nbin-1) {
	hist[0]   += 1.0 * (multiplicity ? multiplicity->at(i) : 1);
      } else {
	hist[ind] += 1.0 * (multiplicity ? multiplicity->at(i) : 1);
      }
    }

    v_hist(&hist[0], min, max, nbin, "");

  } catch(Exception& err) {
    ReportError(err.what());
  }
}

void PgUtil::histogram2D(unsigned ndata, float* data1, float* data2, unsigned nbin1, unsigned nbin2, std::vector<unsigned>* multiplicity)
{
  int waserr = 0;
  int first = 1;
  int i,ind1, ind2;
  double min1,max1,dx1,mid1;
  double min2,max2,dx2,mid2;
  
  std::vector<float> hist2D(nbin1 * nbin2);

  try {
    for(i=0;i < nbin1*nbin2; i++)
      hist2D[i] = 0.0;
  
    // Get the min and max of the data.

    if(!waserr) {

      for(i=0; i < ndata; i++) {

	if(!finite(data1[i]) || !finite(data2[i])) {
	  ThrowColorError("\nData not finite", "red");
	}

	if(first) {
	  min1 = max1 = data1[i];
	  min2 = max2 = data2[i];
	  first = 0;
	}

	min1 = (data1[i] < min1 ? data1[i] : min1);
	max1 = (data1[i] > max1 ? data1[i] : max1);

	min2 = (data2[i] < min2 ? data2[i] : min2);
	max2 = (data2[i] > max2 ? data2[i] : max2);
      }

    }
  
    // And construct the histogram.

    dx1 = (max1-min1)/(nbin1-1);
    dx2 = (max2-min2)/(nbin2-1);

    mid1 = (max1 + min1)/2;
    mid2 = (max2 + min2)/2;

    if(fabs(dx1/mid1) < 1e-12) {
      min1 = 0.9*min1;
      max1 = 1.1*max1;
      dx1 = (max1-min1)/(nbin1-1);
    }

    if(fabs(dx2/mid2) < 1e-12) {
      min2 = 0.9*min2;
      max2 = 1.1*max2;
      dx2 = (max2-min2)/(nbin2-1);
    }

    for(i=0; i < ndata; i++) {

      ind1 = (int)floor((data1[i]-min1)/dx1 + 0.5);
      ind2 = (int)floor((data2[i]-min2)/dx2 + 0.5);

      if(ind1 < 0 || ind1 > nbin1-1 || ind2 < 0 || ind2 > nbin2-1) {
	hist2D[0] += 1.0 * (multiplicity ? multiplicity->at(i) : 1);
      } else {
	hist2D[nbin1 * ind2 + ind1] += 1.0 * (multiplicity ? multiplicity->at(i) : 1);
      }
    }

    double xrange = max1-min1;
    double yrange = max2-min2;

    greyScale(nbin1 * nbin2, &hist2D[0], nbin1, nbin2, 
	      min1, max1, 
	      min2, max2);

  } catch(Exception& err) {
    ReportError(err.what());
  }
}

int PgUtil::v_hist(float* hist, float xmins, float xmaxs, int nbin, std::string title)
{
  float dx,ymaxs,ymins=0.0;
  int i;
  float xmin,xmax,ymin,ymax;
  char answer[10];
  int slen;
  float xpos[2],ypos[2];
  int cancel;
  char key = 'L';
  int docurs;
  bool wasopen;
  int read=1;
  int redisp;

  for(ymaxs=0.0,i=0;i < nbin;i++) 
    ymaxs = (ymaxs > hist[i]) ? ymaxs : hist[i];   

  double xrange = xmaxs - xmins;

  dx = (xmaxs-xmins)/(nbin-1);

  ymaxs *= 1.1;
  
  // If no PGPLOT device has been opened, open one.

  cpgsch(ch_);

  queryDevice(wasopen);
  docurs = haveCursor();

  //------------------------------------------------------------
  // See if plot limits were specified by the user
  //------------------------------------------------------------

  if(usedefs_) {

    if(xmin_ == xmax_) {
      xmin = logPlot_ ? (xmins <= 0.0 ? 0.0 : log10(xmins)) : xmins;
      xmax = logPlot_ ? (xmaxs <= 0.0 ? 0.0 : log10(xmaxs)) : xmaxs;
    } else {

      if(logPlot_ && xmin_ < 0.0)
	xmin = 0.0;
      else
	xmin = (logPlot_ ? log10(xmin_) : xmin_);

      if(logPlot_ && xmax_ < 0.0)
	xmax = 0.0;
      else
	xmax = (logPlot_ ? log10(xmax_) : xmax_);
    };

    if(ymin_ == ymax_) {
      ymin = logPlot_ ? (ymins <= 0.0 ? 0.0 : log10(ymins)) : ymins;
      ymax = logPlot_ ? (ymaxs <= 0.0 ? 0.0 : log10(ymaxs)) : ymaxs;
    } else {

      if(logPlot_ && ymin_ < 0.0)
	ymin = 0.0;
      else
	ymin = (logPlot_ ? log10(ymin_) : ymin_);

      if(logPlot_ && ymax_ < 0.0)
	ymax = 0.0;
      else
	ymax = (logPlot_ ? log10(ymax_) : ymax_);
    };

    //------------------------------------------------------------
    // Plot limits weren't specified by the user
    //------------------------------------------------------------

  } else {
    xmin = logPlot_ ? (xmins <= 0.0 ? 0.0 : log10(xmins)) : xmins;
    xmax = logPlot_ ? (xmaxs <= 0.0 ? 0.0 : log10(xmaxs)) : xmaxs;

    ymin = logPlot_ ? (ymins <= 0.0 ? 0.0 : log10(ymins)) : ymins;
    ymax = logPlot_ ? (ymaxs <= 0.0 ? 0.0 : log10(ymaxs)) : ymaxs;
  }

  if(!usedefs_) {
    PlotAxisRange xrng;
    xrng.setTo(xmins, xmaxs);
    xrng.expandAbs(1.5*dx);
    double xminTmp = xrng.plotMin();
    double xmaxTmp = xrng.plotMax();
    double plxmin = logPlot_ ? ((xminTmp) <= 0.0 ? 0.0 : log10(xminTmp)) : (xminTmp);
    double plxmax = logPlot_ ? ((xmaxTmp) <= 0.0 ? 0.0 : log10(xmaxTmp)) : (xmaxTmp);

    if(win_)
      cpgswin(plxmin, plxmax, ymin, ymax);

  } else {
    if(win_)
      cpgswin(xmin, xmax, ymin, ymax);
  }

  if(docurs) {
    printf("For HELP, hit the \'%c\' key on your keyboard\n", KEY_HELP);
    do {
      redisp = 0;
      cancel = 0;
      switch(key) {
      case KEY_ZOOM:
	if(v_hzoom(xmins,xmaxs,ymins,ymaxs,dx))
	  redisp = 1;
	break;
      case KEY_UT:
	if(v_hnewx(xmins,xmaxs,ymin,ymax,dx))
	  redisp = 1;
	break;
      case KEY_DIS:
	redisp = 1;
	break;
      case KEY_HELP:     /* Print usage info */
	printf("\nYou requested help by pressing \'%c\'.\n", KEY_HELP);
	printf("All cursor positions are entered with the \'%c\' key (Left mouse button)\n", KEY_CUR);
	printf("\n %c - Redisplay current plot.\n", KEY_DIS);
	printf(" %c - Fit a gaussian/polynomial to selected data.\n", KEY_FIT);
	printf(" %c - Select X range to be displayed (hit %c twice for full range)\n", KEY_UT, KEY_UT);
	printf(" %c - Select a sub-image to be displayed.\n", KEY_ZOOM);
	printf("\nTo end this session hit the \'%c\' key (Right mouse button)\n", KEY_QUIT);
	printf("\n");
	break;
      default :
	break;
	
      }
      if(redisp)
      	if(!cancel) 
	  v_hdraw(&hist[0],nbin,xmins,xmaxs,ymins,ymaxs,xmin,xmax,ymin,dx,"N", title);
      if(read) 
	cpgband(B_NORM, 0, xpos[0], ypos[0], &xpos[0], &ypos[0], &key);
      read = 1;
      if(islower((int) key))
	key = (char) toupper((int) key);
    } while(key != KEY_QUIT);
  } else {
    v_hdraw(&hist[0],nbin,xmins,xmaxs,ymins,ymaxs,xmin,xmax,ymin,dx,"N",title.c_str());
  }
  
  return 0;
}

/**.......................................................................
 * Draw the histogram box and envelope
 */
int PgUtil::v_hdraw(float* hist, int nbin, float xmins, float xmaxs, 
		    float ymins, float ymaxs, float xmin, float xmax, 
		    float ymin,float dx, std::string ylab, std::string title)
{
  int i;
  float x1,x2;
  float ypts[5],xpts[5];

  if(!overplot_) 
    cpgpage();

  if(vp_)
    cpgvstd();

  cpgbbuf(); 

  if(useTraceCi_)
    cpgsci(traceCi_);
  else
    cpgsci(10);

  //------------------------------------------------------------
  // If drawing the confidence interval, do it first
  //------------------------------------------------------------

  double refXVal;
  if(draw1SigmaConfidenceInterval_) {
    //    draw1DHistogramConfidenceIntervalAboutMax(hist, nbin, xmins, xmaxs);
    refXVal = drawConfidenceInterval(hist, nbin, xmins, xmaxs);
  }

  // Draw the histogram envelope.

  for(i=0;i < nbin;i++) {
    
    x1 = xmins+(i-0.5)*dx;
    x2 = xmins+(i+0.5)*dx;
    
    xpts[0] = logPlot_ ? (x1 <= 0.0 ? 0.0 : log10(x1)) : x1;
    xpts[1] = logPlot_ ? (x1 <= 0.0 ? 0.0 : log10(x1)) : x1;
    xpts[2] = logPlot_ ? (x2 <= 0.0 ? 0.0 : log10(x2)) : x2;
    
    if(i==0) {
      ypts[0] = logPlot_ ? (ymins <= 0.0 ? 0.0 : log10(ymins)) : ymins;
      ypts[1] = logPlot_ ? (hist[i] <= 0.0 ? 0.0 : log10(hist[i])) : hist[i];
      ypts[2] = logPlot_ ? (hist[i] <= 0.0 ? 0.0 : log10(hist[i])) : hist[i];
      
      cpgline(3,xpts,ypts);

    } else if(i==nbin-1) {

      xpts[3] = logPlot_ ? (x2 <= 0.0 ? 0.0 : log10(x2)) : x2;
      
      ypts[0] = logPlot_ ? (hist[i-1] <= 0.0 ? 0.0 : log10(hist[i-1])) : hist[i-1];
      ypts[1] = logPlot_ ? (hist[i] <= 0.0 ? 0.0 : log10(hist[i])): hist[i];
      ypts[2] = logPlot_ ? (hist[i] <= 0.0 ? 0.0 : log10(hist[i])) : hist[i];
      ypts[3] = logPlot_ ? (ymins <= 0.0 ? 0.0 : log10(ymins)) : ymins;
      
      cpgline(4,xpts,ypts);
    } else {
      ypts[0] = logPlot_ ? (hist[i-1] <= 0.0 ? 0.0 : log10(hist[i-1])) : hist[i-1];
      ypts[1] = logPlot_ ? (hist[i] <= 0.0 ? 0.0 : log10(hist[i])) : hist[i];
      ypts[2] = logPlot_ ? (hist[i] <= 0.0 ? 0.0 : log10(hist[i])) : hist[i];
      
      cpgline(3,xpts,ypts);    
    }
  }

  // If requested to draw the mean, compute and draw it now

  if(drawMean_) {
    //    draw1DHistogramMax(hist, nbin, xmins, xmaxs);
    drawRefMarker(hist, nbin, xmins, xmaxs, refXVal);
  }

  cpgsci(1);
  
  // Redraw the bottom axis.

  xpts[0] = logPlot_ ? log10(xmin) : xmin;
  xpts[1] = logPlot_ ? log10(xmax) : xmax;
  ypts[0] = ypts[1] = logPlot_ ? log10(ymin) : ymin;

  cpgsls(1);
  cpgline(2, xpts, ypts);
  
  // Label the plot

  std::string xOpts, yOpts;

  if(logPlot_) {
    xOpts = (tick_ && xTick_) ? (xTickLabeling_ ? (xTickLabelAtBottom_ ? "BCNLST" : "BCMLST") : "BCLST")  : "BLC";
  } else {
    xOpts = (tick_ && xTick_) ? (xTickLabeling_ ? (xTickLabelAtBottom_ ? "BCNST" : "BCMST") : "BCST")   : "BC";
  }

  if(logPlot_) {
    yOpts = (tick_ && yTick_) ? (yTickLabeling_ ? (yTickLabelAtLeft_ ? "BCNLST" : "BCMLST") : "BCLST") : "BLC";
  } else {
    yOpts = (tick_ && yTick_) ? (yTickLabeling_ ? (yTickLabelAtLeft_ ? "BCNST" : "BCMST") : "BCST")  : "BC";
  }

  // Draw the box only after we're done, so the histogram envelope
  // lines don't overlap it

  if(useBoxCi_)
    cpgsci(boxCi_);
  else
    cpgsci(1);

  cpgbox(xOpts.c_str(), 0.0,0, yOpts.c_str(), 0.0,0);

  std::string xLab, yLab, tLab;

  xLab = (label_ && xLabel_) ? xLabelString_ : "";
  yLab = (label_ && yLabel_) ? yLabelString_ : "";
  tLab = title_ ? title : "";

  cpglab(xLab.c_str(), yLab.c_str(), tLab.c_str());

  cpgebuf(); 

  return 1;
}

void PgUtil::draw1DHistogramConfidenceIntervalAboutMean(float* hist, unsigned nbin, float min, float max)
{
  int ci;
  cpgqci(&ci);

  cpgscr(16, 0.35, 0.35, 0.35);
  cpgsci(16);

  // Calculate the mean of this histogram

  double msum = 0.0, psum = 0.0;
  double dx = (max-min)/(nbin-1);

  for(unsigned i=0; i < nbin; i++) {
    double x = min + (i + 0.5) * dx;

    msum +=  x * hist[i];
    psum +=      hist[i];
  }

  double mean = msum / psum;

  int indMean = (int)floor((mean-min)/dx + 0.5);
  double frac = hist[indMean];

  float xpts[4], ypts[4];
  bool stop=false;
  bool highStop=false;
  bool lowStop=false;

  int iLow=indMean, iHigh=indMean;

  for(int i=1; frac < 0.68 * psum && !stop; i++) {

    iHigh = indMean + i;

    if(iHigh < nbin-1) {

      frac += hist[iHigh];

      //------------------------------------------------------------
      // Drawed filled rectangles within the 1-sigma confidence interval
      //------------------------------------------------------------

      xpts[0] = min + (iHigh - 0.5) * dx;
      xpts[1] = xpts[0];
      xpts[2] = min + (iHigh + 0.5) * dx;
      xpts[3] = xpts[2];
      
      ypts[0] = 0.0;
      ypts[1] = hist[iHigh];
      ypts[2] = hist[iHigh];
      ypts[3] = 0.0;
      
      cpgpoly(4, xpts, ypts);
    } else {
      highStop = true;
    }

    iLow = indMean - i;

    if(iLow >= 0) {

      frac += hist[iLow];

      //------------------------------------------------------------
      // Drawed filled rectangles within the 1-sigma confidence interval
      //------------------------------------------------------------
      
      xpts[0] = min + (iLow - 0.5) * dx;
      xpts[1] = xpts[0];
      xpts[2] = min + (iLow + 0.5) * dx;
      xpts[3] = xpts[2];
      
      ypts[0] = 0.0;
      ypts[1] = hist[iLow];
      ypts[2] = hist[iLow];
      ypts[3] = 0.0;
      
      cpgpoly(4, xpts, ypts);
    } else {
      lowStop = true;
    }

    stop = lowStop && highStop;
  }

  cpgsci(ci);
}

void PgUtil::draw1DHistogramMean(float* hist, unsigned nbin, float min, float max)
{
  // Calculate the mean of this histogram

  double msum = 0.0, psum = 0.0;
  double dx = (max-min)/(nbin-1);

  for(unsigned i=0; i < nbin; i++) {
    double x = min + (i + 0.5) * dx;
    msum += x*hist[i];
    psum += hist[i];
  }

  double mean = msum / psum;

  float ymin, ymax, xmin, xmax;
  cpgqwin(&xmin, &xmax, &ymin, &ymax);
  int ci;
  cpgqci(&ci);

  cpgsci(2);

  cpgmove(mean, ymin);
  cpgdraw(mean, ymax);

  cpgsci(ci);
}

void PgUtil::draw1DHistogramConfidenceIntervalAboutMax(float* hist, unsigned nbin, float min, float max)
{
  int ci;
  cpgqci(&ci);

  cpgscr(16, 0.35, 0.35, 0.35);
  cpgsci(16);

  // Calculate the mean of this histogram

  double dx = (max-min)/(nbin-1);
  double hmax = hist[0];
  double xmax = min + 0.5 * dx;
  int indMax = 0;

  double msum = 0.0, psum = 0.0;

  for(unsigned i=0; i < nbin; i++) {
    double x = min + (i + 0.5) * dx;

    msum += x* hist[i];
    psum +=    hist[i];

    if(hist[i] > hmax) {
      xmax = x;
      hmax = hist[i];
      indMax = i;
    }
  }

  double frac = hist[indMax];

  float xpts[4], ypts[4];
  bool stop=false;
  bool highStop=false;
  bool lowStop=false;

  int iLow=indMax, iHigh=indMax;

  //------------------------------------------------------------
  // Draw the first bin (max value)
  //------------------------------------------------------------

  {
    xpts[0] = min + (indMax - 0.5) * dx;
    xpts[1] = xpts[0];
    xpts[2] = min + (indMax + 0.5) * dx;
    xpts[3] = xpts[2];
    
    ypts[0] = 0.0;
    ypts[1] = hist[indMax];
    ypts[2] = hist[indMax];
    ypts[3] = 0.0;
    
    cpgpoly(4, xpts, ypts);
  }

  //------------------------------------------------------------
  // Now iterate drawing bins within the 1sigma confidence interval    
  //------------------------------------------------------------

  for(int i=1; frac < 0.68 * psum && !stop; i++) {

    iHigh = indMax + i;

    if(iHigh < nbin-1) {

      frac += hist[iHigh];

      //------------------------------------------------------------
      // Drawed filled rectangles within the 1-sigma confidence interval
      //------------------------------------------------------------

      xpts[0] = min + (iHigh - 0.5) * dx;
      xpts[1] = xpts[0];
      xpts[2] = min + (iHigh + 0.5) * dx;
      xpts[3] = xpts[2];
      
      ypts[0] = 0.0;
      ypts[1] = hist[iHigh];
      ypts[2] = hist[iHigh];
      ypts[3] = 0.0;
      
      cpgpoly(4, xpts, ypts);
    } else {
      highStop = true;
    }

    iLow = indMax - i;

    if(iLow >= 0) {

      frac += hist[iLow];

      //------------------------------------------------------------
      // Drawed filled rectangles within the 1-sigma confidence interval
      //------------------------------------------------------------

      xpts[0] = min + (iLow - 0.5) * dx;
      xpts[1] = xpts[0];
      xpts[2] = min + (iLow + 0.5) * dx;
      xpts[3] = xpts[2];
      
      ypts[0] = 0.0;
      ypts[1] = hist[iLow];
      ypts[2] = hist[iLow];
      ypts[3] = 0.0;
      
      cpgpoly(4, xpts, ypts);
    } else {
      lowStop = true;
    }

    stop = lowStop && highStop;
  }

  cpgsci(ci);
}

double PgUtil::drawConfidenceInterval(float* hist, unsigned nbin, float min, float max)
{
  int ci;
  cpgqci(&ci);

  cpgscr(16, 0.35, 0.35, 0.35);
  cpgsci(16);

  double perc = 1.0 - 2*Sampler::gaussCdf(-nSigma_, 0.0, 1.0);

  //------------------------------------------------------------
  // Calculate the mean of this histogram
  //------------------------------------------------------------

  double dx = (max-min)/(nbin-1);
  double hmax = hist[0];
  double xmax = min + 0.5 * dx;

  int indMode = 0;
  int indMin  = 0;
  int indMax  = nbin-1;

  int indRef = 0;

  double msum = 0.0, psum = 0.0;

  for(unsigned i=0; i < nbin; i++) {
    double x = min + (i + 0.5) * dx;

    msum += x* hist[i];
    psum +=    hist[i];

    if(hist[i] > hmax) {
      xmax = x;
      hmax = hist[i];
      indMode = i;
    }
  }

  switch(stat_) {
  case PgUtil::STAT_MAXL:
    indRef = indMode;
    break;
  case PgUtil::STAT_UPPER:
    indRef = indMin;
    break;
  default:
    indRef = indMax;
    break;
  }

  double frac = hist[indRef];

  float xpts[4], ypts[4];
  bool stop=false;
  bool highStop=false;
  bool lowStop=false;
  double refXVal = min + indRef * dx;
  double lastLowXVal  = min + indRef * dx;
  double lastHighXVal = min + indRef * dx;
  int iLow=indRef, iHigh=indRef;

  //------------------------------------------------------------
  // Draw the first bin (max value)
  //------------------------------------------------------------

  {
    xpts[0] = min + (indRef - 0.5) * dx;
    xpts[1] = xpts[0];
    xpts[2] = min + (indRef + 0.5) * dx;
    xpts[3] = xpts[2];
    
    ypts[0] = 0.0;
    ypts[1] = hist[indRef];
    ypts[2] = hist[indRef];
    ypts[3] = 0.0;
    
    cpgpoly(4, xpts, ypts);
  }

  //------------------------------------------------------------
  // Now iterate drawing bins within the 1sigma confidence interval    
  //------------------------------------------------------------

  for(int i=1; frac < perc * psum && !stop; i++) {

    iHigh = indRef + i;

    if(iHigh < nbin-1) {

      frac += hist[iHigh];

      //------------------------------------------------------------
      // Drawed filled rectangles within the 1-sigma confidence interval
      //------------------------------------------------------------

      lastHighXVal = min + iHigh * dx;

      xpts[0] = min + (iHigh - 0.5) * dx;
      xpts[1] = xpts[0];
      xpts[2] = min + (iHigh + 0.5) * dx;
      xpts[3] = xpts[2];
      
      ypts[0] = 0.0;
      ypts[1] = hist[iHigh];
      ypts[2] = hist[iHigh];
      ypts[3] = 0.0;
      
      cpgpoly(4, xpts, ypts);
    } else {
      highStop = true;
    }

    iLow = indRef - i;

    if(iLow >= 0) {

      frac += hist[iLow];

      //------------------------------------------------------------
      // Drawed filled rectangles within the 1-sigma confidence interval
      //------------------------------------------------------------

      lastLowXVal = min + iLow * dx;

      xpts[0] = min + (iLow - 0.5) * dx;
      xpts[1] = xpts[0];
      xpts[2] = min + (iLow + 0.5) * dx;
      xpts[3] = xpts[2];
      
      ypts[0] = 0.0;
      ypts[1] = hist[iLow];
      ypts[2] = hist[iLow];
      ypts[3] = 0.0;
      
      cpgpoly(4, xpts, ypts);
    } else {
      lowStop = true;
    }

    stop = lowStop && highStop;
  }

  cpgsci(ci);

  switch (stat_) {
  case PgUtil::STAT_MAXL:
    return refXVal;
    break;
  case PgUtil::STAT_UPPER:
    return lastHighXVal;
    break;
  default:
    return lastLowXVal;
    break;
  }
}

void PgUtil::draw1DHistogramMax(float* hist, unsigned nbin, float min, float max)
{
  // Calculate the mean of this histogram

  double dx = (max-min)/(nbin-1);
  double hymax = hist[0];
  double hxmax = min + 0.5 * dx;

  for(unsigned i=0; i < nbin; i++) {
    double x = min + (i + 0.5) * dx;
    if(hist[i] > hymax) {
      hxmax = x;
      hymax = hist[i];
    }
  }

  int indMax = (int)floor((max-min)/dx + 0.5);
  double frac = hist[indMax];

  float ymin, ymax, xmin, xmax;
  cpgqwin(&xmin, &xmax, &ymin, &ymax);
  int ci;
  cpgqci(&ci);

  cpgsci(2);

  cpgmove(hxmax, ymin);
  cpgdraw(hxmax, ymax);

  cpgsci(ci);
}

void PgUtil::drawRefMarker(float* hist, unsigned nbin, float min, float max, double refXVal)
{
  // Calculate the mean of this histogram

  double dx = (max-min)/(nbin-1);
  double hymax = hist[0];
  double hxmax = min + 0.5 * dx;
  double hy = 0.0;

  for(unsigned i=0; i < nbin; i++) {
    double x = min + (i + 0.5) * dx;
    if(hist[i] > hymax) {
      hxmax = x;
      hymax = hist[i];
    }
  }

  float ymin, ymax, xmin, xmax;
  cpgqwin(&xmin, &xmax, &ymin, &ymax);
  int ci;
  cpgqci(&ci);

  cpgsci(6);

  cpgmove(refXVal, ymin);
  cpgdraw(refXVal, ymax);

  cpgsci(ci);
}

/*.......................................................................
 * Module for v_hist.
 */
int PgUtil::v_hzoom(float xmins, float xmaxs, float ymins, float ymaxs, float dx)
{
  int iter;
  int oldcol;
  int dofull;
  int cancel;
  int accepted;
  float xpos[2],ypos[2];
  float xtemp,ytemp;
  char key;
  float xmin,xmax,ymin,ymax;
  static string mess1("\n %c - Select start and end vertices using this key\n");
  static string mess2(" %c - Select the whole plot\n");
  static string mess3(" %c - Abort selection\n");

  cpgqci(&oldcol);
  cpgsci(5);

  dofull = 0;
  cancel = 0;
  for(iter = 0;iter<2 && !dofull && !cancel;iter++) {
    do {
      accepted = 0;
      cpgband((iter==0) ? B_NORM : B_RECT, 0, xtemp, ytemp, &xtemp,
	      &ytemp, &key);
      if(islower((int) key))
	key = (char) toupper((int) key);
      xpos[iter] = xtemp;
      ypos[iter] = ytemp;
      switch(key) {
      case KEY_ZOOM:
	accepted = dofull = 1;
	break;
      case KEY_CAN:      /* Abort box selection */
	accepted = cancel = 1;
	return 0;
	break;
      case KEY_QUIT:     /* Quit now */
	cpgend();
	return 0;
	break;
      case KEY_CUR:             /* Accept the selected start vertex */
	accepted=1;
	break;
      default:            /* Unexpected cursor input key - show usage */
	printf(mess1.c_str(), KEY_CUR);
	printf(mess2.c_str(), KEY_ZOOM);
	printf(mess3.c_str(), KEY_CAN);
	break;
      };
    } while(!accepted);
  };
  if(dofull) {
    xmin = xmins;
    ymin = ymins;
    xmax = xmaxs;
    ymax = ymaxs;
  }
  else {
    xmin = (xpos[0] < xpos[1]) ? xpos[0] : xpos[1];
    ymin = (ypos[0] < ypos[1]) ? ypos[0] : ypos[1];
    xmax = (xpos[0] > xpos[1]) ? xpos[0] : xpos[1];
    ymax = (ypos[0] > ypos[1]) ? ypos[0] : ypos[1];
    if(ymin < 0.0)
      ymin = 0.0;
  }
  if(!usedefs_)
    cpgswin(xmin-dx,xmax+dx,ymin,ymax);
  cpgsci(oldcol);

  return 1;
}
/*.......................................................................
 * Module for v_hist.
 */
int PgUtil::v_hnewx(float xmins, float xmaxs, float ymins, float ymaxs, float dx)
{
  int iter;
  int oldcol;
  int dofull;
  int cancel;
  int accepted;
  float xpos[2],ypos[2];
  float xtemp,ytemp;
  char key;
  float xmin,xmax,ymin,ymax;
  static string mess1("\n %c - Select start and end vertices using this key\n");
  static string mess2(" %c - Select the whole plot\n");
  static string mess3(" %c - Abort selection\n");

  /*
   * Print help messages
   */
  printf(mess1.c_str(), KEY_CUR);
  printf(mess2.c_str(), KEY_UT);
  printf(mess3.c_str(), KEY_CAN);

  cpgqci(&oldcol);
  cpgsci(5);

  dofull = 0;
  cancel = 0;
  for(iter = 0;iter<2 && !dofull && !cancel;iter++) {
    do {
      accepted = 0;
      cpgband((iter==0) ? B_XVAL : B_XRNG, 0, xtemp, ytemp, &xtemp,
	      &ytemp, &key);
      if(islower((int) key))
	key = (char) toupper((int) key);
      xpos[iter] = xtemp;
      ypos[iter] = ytemp;
      switch(key) {
      case KEY_UT:
	accepted = dofull = 1;
	break;
      case KEY_CAN:      /* Abort box selection */
	accepted = cancel = 1;
	return 0;
	break;
      case KEY_QUIT:     /* Quit now */
	cpgend();
	return 0;
	break;
      case KEY_CUR:             /* Accept the selected start vertex */
	accepted=1;
	break;
      default:            /* Unexpected cursor input key - show usage */
	printf(mess1.c_str(), KEY_CUR);
	printf(mess2.c_str(), KEY_UT);
	printf(mess3.c_str(), KEY_CAN);
	break;
      };
    } while(!accepted);
  };
  if(!cancel) {
    if(dofull) {
      xmin = xmins;
      ymin = ymins;
      xmax = xmaxs;
      ymax = ymaxs;
    }
    else {
      xmin = (xpos[0] < xpos[1]) ? xpos[0] : xpos[1];
      ymin = ymins;
      xmax = (xpos[0] > xpos[1]) ? xpos[0] : xpos[1];
      ymax = ymaxs;
    }
    if(!usedefs_)
      cpgswin(xmin-dx,xmax+dx, logPlot_ ? log10(ymin):ymin, logPlot_ ? log(ymax):ymax);
  }
  cpgsci(oldcol);

  return 1;
}

void PgUtil::setPrompt(bool prompt)
{
  if(haveOpenDevice()) {
    cpgask(prompt);
  }
}
     
void PgUtil::substituteForPgplot(std::string& str)
{
  String pgStr(str);

  //------------------------------------------------------------
  // Replace special characters that won't render correctly with their
  // Hershey symbol equivalents
  //------------------------------------------------------------

  pgStr.superscriptNumbersForPgplot();
  pgStr.replace("^", "\\(0756)");
  pgStr.replace("solar", "\\d\\(2281)\\u");
  pgStr.replace("muK", "\\(2138)K");

  str = pgStr.str();
}

void PgUtil::label(std::string xlabel, std::string ylabel, std::string title)
{
  substituteForPgplot(xlabel);
  substituteForPgplot(ylabel);
  substituteForPgplot(title);

  cpglab(xlabel.c_str(), ylabel.c_str(), title.c_str());
}

void PgUtil::setPostscriptFontName(std::string fontName)
{
  String str(fontName);
  String strLow = str.toLower();

  bool match = false;
  for(unsigned i=0; i < fontNames_.size(); i++) {
    String str1(fontNames_[i]);
    String str1Low = str1.toLower();
    if(str1Low.contains(strLow.str())) {
      match = true;
      break;
    }
  }

  if(match)
    setenv("PGPLOT_PS_FONT", fontName.c_str(), 1);
  else {
    XtermManip xtm;
    std::ostringstream os;
    for(unsigned i=0; i < fontNames_.size(); i++) {
      os << fontNames_[i] << std::endl;
    }

    ThrowError(COLORIZE(xtm, "red", std::endl << "Font " << fontName << " is not supported" << std::endl << std::endl) << 
	       COLORIZE(xtm, "green", "Supported fonts are: " << std::endl << std::endl << os.str() << std::endl));
  }
}

std::vector<std::string> PgUtil::getFontNames()
{
  std::vector<std::string> fn;

  fn.push_back("Symbol");
  fn.push_back("Times-Roman");
  fn.push_back("Times-Italic");
  fn.push_back("Times-Bold");
  fn.push_back("Times-BoldItalic");
  fn.push_back("Helvetica");
  fn.push_back("Helvetica-Oblique");
  fn.push_back("Helvetica-Bold");
  fn.push_back("Helvetica-BoldOblique");
  fn.push_back("Courier");
  fn.push_back("Courier-Oblique");
  fn.push_back("Courier-Bold");
  fn.push_back("Courier-BoldOblique");
  fn.push_back("NewCenturySchlbk-Roman");
  fn.push_back("NewCenturySchlbk-Italic");
  fn.push_back("NewCenturySchlbk-Bold");
  fn.push_back("NewCenturySchlbk-BoldItalic");
  fn.push_back("ZapfChancery");
  fn.push_back("ZapfChancery-Oblique");
  fn.push_back("ZapfChancery-Bold");
  fn.push_back("URWGroteskT-Bold");
  fn.push_back("URWAntiquaT-RegularCondensed");
  fn.push_back("Cyrillic");
  fn.push_back("ZapfDingbats");

  return fn;
}

void PgUtil::displayColors()
{
  open("/xs");
  cpgvstd();
  cpgbox("BCNST", 0, 0, "BCNST", 0, 0);
  cpgswin(-1, 16, -1, 1);

  for(unsigned i=0; i < 15; i++) {
    std::ostringstream os;
    os << i;
    cpgsci(i);
    cpgtext(i, 0, os.str().c_str());
  }

  close();
}

void PgUtil::clearPgManager()
{
  pgManager_.models_.resize(0);
}

void PgUtil::insertPgManager(PgModelManager& pm)
{
  pgManager_ = pm;
}

void PgUtil::drawBox(unsigned ci, bool override)
{
  std::string xOpts, yOpts;

  if(logPlot_) {
    xOpts = (tick_ && xTick_) ? (xTickLabeling_ ? (xTickLabelAtBottom_ ? "BCNLST" : "BCMLST") : "BCLST")  : "BLC";
  } else {
    xOpts = (tick_ && xTick_) ? (xTickLabeling_ ? (xTickLabelAtBottom_ ? "BCNST" : "BCMST") : "BCST")   : "BC";
  }

  if(logPlot_) {
    yOpts = (tick_ && yTick_) ? (yTickLabeling_ ? (yTickLabelAtLeft_ ? "BCNLST" : "BCMLST") : "BCLST") : "BLC";
  } else {
    yOpts = (tick_ && yTick_) ? (yTickLabeling_ ? (yTickLabelAtLeft_ ? "BCNST" : "BCMST") : "BCST")  : "BC";
  }

  if(box_) {
    if(useBoxCi_ && !override)
      cpgsci(boxCi_);
    else
      cpgsci(ci);
    
    cpgbox(xOpts.c_str(), 0.0,0, yOpts.c_str(),0.0,0);
  }
}

void PgUtil::drawLabels(std::string xlab, std::string ylab, std::string title)
{
  std::string xLab, yLab, tLab;
  
  xLab = (label_ && xLabel_) ? xlab : "";
  yLab = (label_ && yLabel_) ? ylab : "";
  tLab = title_ ? title : "";
  
  cpglab(xLab.c_str(), yLab.c_str(), tLab.c_str());
}

void PgUtil::drawWedge(double zmin, double zmax, std::string unit)
{
  // Draw a ramp on the side.
  
  if(wedge_)
    cpgwedg("RI",0,4,zmax,zmin,unit.c_str()); 
}

void PgUtil::setupPlotBoundaries(double xmin, double xmax, double ymin, double ymax)
{
  float x1, x2, y1, y2;

  cpgqvp(0, &x1, &x2, &y1, &y2);

  if(!overplot_)
    cpgpage();
  
  if(vp_)
    cpgvstd();
  
  if(win_)
    cpgswin(xmin,xmax,ymin,ymax);
  
  if(wnad_)
    cpgwnad(xmin,xmax,ymin,ymax);
}


void PgUtil::plotPoint(double x, double y, int marker)
{
  int ci;
  cpgqci(&ci);
  cpgsci(5);
  cpgpt1(x, y, marker);
  cpgsci(ci);
}

void PgUtil::addModel(PgModel& model)
{
  pgManager_.addModel(model);
}

void PgUtil::initialize()
{
  useBoxCi_ = false;
  boxCi_ = 10;
  useTraceCi_ = false;
  traceCi_ = 10;
  cmap_ = &std_cmaps[4];
  overplot_ = false;
  reverseX_ = false;
  reverseY_ = false;
  nContour_ = 0;
  vp_ = true;
  win_ = true;
  box_ = true;
  wnad_ = false;
  tick_ = true;
  xTick_ = true;
  yTick_ = true;
  xTickLabeling_ = true;
  yTickLabeling_ = true;
  xTickLabelAtBottom_ = true;
  yTickLabelAtLeft_ = true;
  label_ = true;
  xLabel_ = true;
  yLabel_ = true;
  title_ = true;
  wedge_ = true;
  usedefs_ = false;
  xmin_ = 0.0;
  xmax_ = 0.0;
  ymin_ = 0.0;
  ymax_ = 0.0;
  ch_ = 1.0;
  interactive_ = true;
  logPlot_     = false;
  xnlab_       = true;
  drawMean_    = true;
  draw1SigmaConfidenceInterval_ = false;
  zmin_ = 0.0;
  zmax_ = 0.0;
  unitCallback_   = 0;
  unitCallbackArgs_  = 0;
  coordCallback_ = 0;
  coordCallbackArgs_ = 0;
  priorCallback_ = 0;
  priorCallbackArgs_ = 0;

  xLabelString_ = "";
  yLabelString_ = "";

  fontNames_ = PgUtil::getFontNames();

  stat_ = PgUtil::STAT_NONE;
  nSigma_ = 1.0;

  useHeader_ = false;
}
