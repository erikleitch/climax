#include "gcp/util/Exception.h"
#include "gcp/util/Fitter.h"
#include "gcp/util/Stats.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/fftutil/Dft1d.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

using namespace gcp::util;

#define EPS 1e-10
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

/**.......................................................................
 * Constructor.
 */
Fitter::Fitter() {}

/**.......................................................................
 * Destructor.
 */
Fitter::~Fitter() {}

/**.......................................................................
 * Take a chain of accepted samples, compute the power spectrum of the
 * chain and fit a model to the power spectrum
 */
void Fitter::fitPowerSpectrum(std::vector<double>& yarr, 
                              std::vector<unsigned>* multiplicity, 
                              unsigned nFit,
                              double& p0, double& ks, double& alpha,
                              unsigned& nSamp)
{
  //------------------------------------------------------------
  // Get the power spectrum and inital guesses for fit parameters
  //------------------------------------------------------------

  Fitter::Fit_scan scan = 
    getPowerSpectrumFitDataNew(yarr, multiplicity, nFit, p0, ks, alpha, nSamp);

  //------------------------------------------------------------
  // Now fit
  //------------------------------------------------------------

  std::vector<float> pars(3);
  float chisq;

  pars[0] = p0;
  pars[2] = alpha;
  pars[1] = ks;

  //------------------------------------------------------------
  // We fit twice, since this is more robust to bad initial guesses
  // for the parameters
  //------------------------------------------------------------

  fit(&scan, &pars[0], 3, Fitter::F_POWSPEC, &chisq, 0);
  fit(&scan, &pars[0], 3, Fitter::F_POWSPEC, &chisq, 0);

  p0    = pars[0];
  ks    = pars[1];
  alpha = pars[2];

  double r = p0/nSamp;
  double js = ks/(2*M_PI)*nSamp;
}

/**.......................................................................
 * Take a chain of accepted samples, compute the power spectrum of the
 * chain, fit a model to the power spectrum and plot it
 */
void Fitter::plotPowerSpectrum(std::vector<double>& yarr, std::vector<unsigned>* multiplicity, unsigned nFit,
                               double& p0, double& ks, double& alpha)
{
  try {
    //------------------------------------------------------------
    // Get the power spectrum and inital guesses for fit parameters
    //------------------------------------------------------------
    
    unsigned nSamp;
    
    Fitter::Fit_scan scan = getPowerSpectrumFitData(yarr, multiplicity, p0, ks, alpha, nSamp);
    
    //------------------------------------------------------------
    // Now fit and plot
    //------------------------------------------------------------
    
    std::vector<float> pars(3);
    float chisq;
    
    pars[0] = p0;
    pars[2] = alpha;
    pars[1] = ks;
    
    PgUtil::setLogPlot(true);
    PgUtil::linePlot(scan.x, scan.y, "", "", "", true);
    
    //------------------------------------------------------------
    // We fit twice, since this is more robust to bad initial guesses
    // for the parameters
    //------------------------------------------------------------
    
    fit(&scan, &pars[0], 3, Fitter::F_POWSPEC, &chisq, 0);
    fit(&scan, &pars[0], 3, Fitter::F_POWSPEC, &chisq, 0);
    
    p0    = pars[0];
    ks    = pars[1];
    alpha = pars[2];
    
    double r = p0/nSamp;
    double js = ks/(2*M_PI)*nSamp;
    
    //  COUT("nSamp = " << nSamp << " p0 = " << p0 << " ks = " << ks << " alpha = " << alpha << " r = " << p0/nSamp << " j* = " << ks/(2*M_PI)*nSamp);
    //  COUT("Converged = " << ((r < 0.01 && js > 20) ? "true" : "false"));
    
    for(unsigned i=0; i < scan.x.size(); i++) {
      double k = scan.x[i];
      double r = ks/k;
      double f = powf(r, alpha);
      scan.y[i] = p0 * f / (f+1);
    }
    
    PgUtil::setOverplot(true);
    PgUtil::setWin(false);
    PgUtil::setBox(false);
    PgUtil::setTraceColor(6);
    
    PgUtil::linePlot(scan.x, scan.y, "", "", "", true);
    
    PgUtil::setWin(true);
    PgUtil::setBox(true);
    PgUtil::clearTraceColor();
  } catch(...) {
  }
}

/**.......................................................................
 * Fit some data
 */
int Fitter::fit(Fit_scan *fitscan, float apar[], int ma, Fitter::Func_type type, 
		float *chisq_fin, int printchisq)
{
  float *a=0,**covar=0,**alpha=0,chisq,alamda,chisqb;
  int i,j,waserr=0,mfit,*lista=0;
  void (*funcs)(float x, float a[], float *y, float dyda[], int na);
  void (*funcs2d)(float x1, float x2, float a[], float *y, float dyda[], int na);

  bool is2d = false;

  if(fitscan==0)
    return 1;
  
  //------------------------------------------------------------
  // First read through the fitscan to check if the data look ok
  //------------------------------------------------------------

  if(fitscan->x.size() > 0) {
    for(i=0;i < fitscan->npt;i++) 
      if(!isfinite(fitscan->x[i])) {
	waserr |= 1;
	break;
      }
  } else
    waserr = 1;

  if(fitscan->y.size() > 0) {
    for(i=0;i < fitscan->npt;i++) 
      if(!isfinite(fitscan->y[i])) {
	waserr |= 1;
	break;
      }
  } else 
    waserr = 1;

  if(fitscan->sd.size() > 0)
    for(i=0;i < fitscan->npt;i++) 
      if(!isfinite(fitscan->sd[i])) {
	waserr |= 1;
	break;
      }
  
  //------------------------------------------------------------
  // Allocate the matrices we need to do the fitting
  //------------------------------------------------------------

  if(!waserr)
    waserr |= (covar=matrix(ma,ma))==0;
  if(!waserr)
    waserr |= (alpha=matrix(ma,ma))==0;
  if(!waserr)
    waserr |= (a = (float *)malloc(ma*sizeof(float)))==0;
  if(!waserr)
    waserr |= (lista = (int *)malloc(ma*sizeof(int)))==0;

  if(!waserr) {
    for(i=0;i < ma;i++) {
      a[i] = apar[i];
      lista[i] = 1;
    }
    
    mfit = ma;

    switch(type) {
    case(F_GAUSS):
      funcs = gauss_nobase;
      break;
    case(F_POLY):
      funcs = poly;
      break;
    case(F_POWER):
      funcs = power;
      break;
    case(F_POWSPEC):
      funcs = powspec;
      break;
    case(F_GAUSS2D):
      funcs2d = gauss2d;
      is2d = true;
      break;
    case(F_GAUSSTEST):
      funcs = gausstest;
      break;
    default:
      funcs = gauss_nobase;
      break;
    };
    
    //------------------------------------------------------------
    // Now we've initialized the parameters and are ready to fit.
    //------------------------------------------------------------

    alamda = -1.0;
    chisqb =  1.0;
    chisq  =  0.0;
    
    //------------------------------------------------------------
    // Main loop -- try to reduce chisq until the new value differs
    // from the old by < EPS.
    //------------------------------------------------------------

    //    COUT("chisq = " << chisq << " chisqb = " << chisqb << " EPS = " << EPS << " waserr = " << waserr);

    while (fabs(chisqb-chisq) > EPS && !waserr) {

      chisqb = chisq;

      if(is2d) {
	waserr |= mrqmin2d(fitscan->x,fitscan->x2,fitscan->y,fitscan->sd,a,lista,
			   ma, covar,alpha,&chisq,funcs2d,&alamda);
      } else {
	waserr |= mrqmin(fitscan->x,fitscan->y,fitscan->sd,a,lista,
			 ma, covar,alpha,&chisq,funcs,&alamda);
      }

      //      COUT("chisq = " << chisq << " chisqb = " << chisqb << " waserr = " << waserr);
    }

    if(!waserr) {
      alamda = 0.0;

      if(is2d) {
	waserr |= mrqmin2d(fitscan->x,fitscan->x2,fitscan->y,fitscan->sd,a,lista,
			   ma,covar,alpha,&chisq,funcs2d,&alamda);
      } else {
	waserr |= mrqmin(fitscan->x,fitscan->y,fitscan->sd,a,lista,
			 ma,covar,alpha,&chisq,funcs,&alamda);
      }
      
      for(i=0;i < ma;i++)
	apar[i] = a[i];

      //------------------------------------------------------------
      // Convert back to gaussian amplitude
      //------------------------------------------------------------

      if(type == F_GAUSS) {

	//	apar[1] = a[1]/sqrt(2*M_PI*a[0]*a[0]);

	//------------------------------------------------------------
	// Return the positive solution for sigma only.
	//------------------------------------------------------------

	apar[1] = (apar[1] < 0) ? -apar[1] : apar[1];
      }

      //------------------------------------------------------------
      // Convert back to gaussian amplitude
      //------------------------------------------------------------

      if(type == F_GAUSS2D) {

	//------------------------------------------------------------
	// Return the positive solution for sigma only.
	//------------------------------------------------------------

	apar[2] = (apar[2] < 0) ? -apar[2] : apar[2];
	apar[4] = (apar[4] < 0) ? -apar[4] : apar[4];
      }

      *chisq_fin = chisq;
      
      //------------------------------------------------------------
      // Print the covariance matrix before exiting
      //------------------------------------------------------------

      if(printchisq) {

	fprintf(stdout,"\n");  
	for(i=0;i < ma;i++) {
	  fprintf(stdout,"%8.2g ",apar[i]);
	}

	fprintf(stdout,"\n");  
	for(i=0;i < ma;i++) {
	  fprintf(stdout,"\n\t");
	  for(j=0;j < ma;j++)
	    fprintf(stdout,"%8.2g ",covar[i][j]);
	}

	fprintf(stdout,"\n");
      }
    }
  }
  
  //------------------------------------------------------------
  // Free any allocated memory
  //------------------------------------------------------------

  if(a)
    free(a);
  if(lista)
    free(lista);
  if(covar)
    free_matrix(covar);
  if(alpha)
    free_matrix(alpha);

  return waserr;
}

int Fitter::mrqcof(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, float a[], 
		   int ia[], int ma, float **alpha, float beta[], float *chisq,
		   FITTER_EVAL_FN(*funcs))
{
  int i,j,k,l,m,mfit=0,waserr=0;
  float ymod,wt,sig2i,dy,*dyda=0;
  unsigned ndata = x.size();

  waserr = (dyda=(float *)malloc(ma*sizeof(float)))==0;

  if(!waserr) {
    for (j=0;j<ma;j++)
      if (ia[j]) mfit++;
    for (j=0;j<mfit;j++) {
      for (k=0;k<=j;k++) alpha[j][k]=0.0;
      beta[j]=0.0;
    }
    *chisq=0.0;

    for (i=0;i<ndata;i++) {
      (*funcs)(x[i],a,&ymod,dyda,ma);
      sig2i=1.0/(sig[i]*sig[i]);
      dy=y[i]-ymod;
      for (j=-1,l=0;l<ma;l++) {
	if (ia[l]) {
	  wt=dyda[l]*sig2i;
	  for (j++,k=-1,m=0;m<=l;m++)
	    if (ia[m]) alpha[j][++k] += wt*dyda[m];
	  beta[j] += dy*wt;
	}
      }

      *chisq += dy*dy*sig2i;
    }

    for (j=1;j<mfit;j++)
      for (k=0;k<j;k++) alpha[k][j]=alpha[j][k];
  }
  
  // Free allocated memory.

  if(dyda)
    free(dyda);

  return waserr;
}

int Fitter::mrqcof2d(std::vector<double>& x1, std::vector<double>& x2, std::vector<double>& y, 
		     std::vector<double>& sig, float a[], 
		     int ia[], int ma, float **alpha, float beta[], float *chisq,
		     FITTER_EVAL_FN2D(*funcs))
{
  int i,j,k,l,m,mfit=0,waserr=0;
  float ymod,wt,sig2i,dy,*dyda=0;
  unsigned ndata = x1.size();

  waserr = (dyda=(float *)malloc(ma*sizeof(float)))==0;

  if(!waserr) {
    for (j=0;j<ma;j++)
      if (ia[j]) mfit++;
    for (j=0;j<mfit;j++) {
      for (k=0;k<=j;k++) alpha[j][k]=0.0;
      beta[j]=0.0;
    }
    *chisq=0.0;
    for (i=0;i<ndata;i++) {
      (*funcs)(x1[i],x2[i],a,&ymod,dyda,ma);
      sig2i=1.0/(sig[i]*sig[i]);
      dy=y[i]-ymod;
      for (j=-1,l=0;l<ma;l++) {
	if (ia[l]) {
	  wt=dyda[l]*sig2i;
	  for (j++,k=-1,m=0;m<=l;m++)
	    if (ia[m]) alpha[j][++k] += wt*dyda[m];
	  beta[j] += dy*wt;
	}
      }
      *chisq += dy*dy*sig2i;
    }
    for (j=1;j<mfit;j++)
      for (k=0;k<j;k++) alpha[k][j]=alpha[j][k];
  }
  
  // Free allocated memory.

  if(dyda)
    free(dyda);

  return waserr;
}

/*.......................................................................
 * Levenberg-Marquardt non-linear least squares.
 */
int Fitter::mrqmin2d(std::vector<double>& x1, std::vector<double>& x2, std::vector<double>& y, 
		     std::vector<double>& sig, float a[], 
		     int ia[], int ma, float **covar, float **alpha, float *chisq,
		     FITTER_EVAL_FN2D(*funcs),
		     float *alamda)
{
  int j,k,l,waserr=0;
  static int mfit;
  static float ochisq,*atry=0,*beta=0,*da=0,**oneda=0;
  unsigned ndata = x1.size();

  if (*alamda < 0.0) {

    waserr  = (atry=(float *)malloc(ma*sizeof(float)))==0;
    waserr |= (beta=(float *)malloc(ma*sizeof(float)))==0;
    waserr |= (da=(float *)malloc(ma*sizeof(float)))==0;

    if(!waserr) {
      for (mfit=0,j=0;j < ma;j++)
	if (ia[j]) mfit++;
      oneda=matrix(mfit,1);
      *alamda=0.001;
      waserr |= mrqcof2d(x1,x2,y,sig,a,ia,ma,alpha,beta,chisq,funcs);
      ochisq=(*chisq);
      for (j=0;j<ma;j++) atry[j]=a[j];
    }
  }
  if(!waserr) {
    for (j=0;j<mfit;j++) {
      for (k=0;k<mfit;k++) 
	covar[j][k]=alpha[j][k];
      covar[j][j]=alpha[j][j]*(1.0+(*alamda));
      oneda[j][0]=beta[j];
    }

    waserr |= gaussj(covar,mfit,oneda,1);

    if(!waserr) {
      for (j=0;j<mfit;j++) da[j]=oneda[j][0];
      
      if (*alamda == 0.0) {
	covsrt(covar,ma,ia,mfit);
	
	// If we are done, free all allocated memory.

	free_matrix(oneda);

	if(da)
	  free(da);
	if(beta)
	  free(beta);
	if(atry)
	  free(atry);

	return waserr;
      }
      for (j=0,l=0;l<ma;l++)
	if (ia[l]) 
	  atry[l]=a[l]+da[j++];
      
      waserr |= mrqcof2d(x1,x2,y,sig,atry,ia,ma,covar,da,chisq,funcs);
     
      if(!waserr) {
	if (*chisq < ochisq) {
	  *alamda *= 0.1;
	  ochisq=(*chisq);
	  for (j=0;j<mfit;j++) {
	    for (k=0;k<mfit;k++) 
	      alpha[j][k]=covar[j][k];
	    beta[j]=da[j];
	  }
	  for (l=0;l<ma;l++) a[l]=atry[l];
	} 
	else {
	  *alamda *= 10.0;
	  *chisq=ochisq;
	}
      }
    }
  }
  
  // If we encountered an error, make sure any allocated memory is
  // freed before we exit.

  if(waserr) {
    free_matrix(oneda);
    if(da)
      free(da);
    if(beta)
      free(beta);
    if(atry)
      free(atry);
  }
  return waserr;
}

/*.......................................................................
 * Levenberg-Marquardt non-linear least squares.
 */
int Fitter::mrqmin(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, float a[], 
		   int ia[], int ma, float **covar, float **alpha, float *chisq,
		   FITTER_EVAL_FN(*funcs),
		   float *alamda)
{
  int j,k,l,waserr=0;
  static int mfit;
  static float ochisq,*atry=0,*beta=0,*da=0,**oneda=0;
  unsigned ndata = x.size();

  if (*alamda < 0.0) {

    waserr  = (atry=(float *)malloc(ma*sizeof(float)))==0;
    waserr |= (beta=(float *)malloc(ma*sizeof(float)))==0;
    waserr |= (da=(float *)malloc(ma*sizeof(float)))==0;

    if(!waserr) {
      for (mfit=0,j=0;j < ma;j++)
	if (ia[j]) mfit++;
      oneda=matrix(mfit,1);
      *alamda=0.001;
      waserr |= mrqcof(x,y,sig,a,ia,ma,alpha,beta,chisq,funcs);
      ochisq=(*chisq);
      for (j=0;j<ma;j++) atry[j]=a[j];
    }
  }
  if(!waserr) {
    for (j=0;j<mfit;j++) {
      for (k=0;k<mfit;k++) 
	covar[j][k]=alpha[j][k];
      covar[j][j]=alpha[j][j]*(1.0+(*alamda));
      oneda[j][0]=beta[j];
    }

    waserr |= gaussj(covar,mfit,oneda,1);

    if(!waserr) {
      for (j=0;j<mfit;j++) da[j]=oneda[j][0];
      
      if (*alamda == 0.0) {
	covsrt(covar,ma,ia,mfit);
	
	// If we are done, free all allocated memory.

	free_matrix(oneda);

	if(da)
	  free(da);
	if(beta)
	  free(beta);
	if(atry)
	  free(atry);

	return waserr;
      }
      for (j=0,l=0;l<ma;l++)
	if (ia[l]) 
	  atry[l]=a[l]+da[j++];
      
      waserr |= mrqcof(x,y,sig,atry,ia,ma,covar,da,chisq,funcs);
     
      if(!waserr) {
	if (*chisq < ochisq) {
	  *alamda *= 0.1;
	  ochisq=(*chisq);
	  for (j=0;j<mfit;j++) {
	    for (k=0;k<mfit;k++) 
	      alpha[j][k]=covar[j][k];
	    beta[j]=da[j];
	  }
	  for (l=0;l<ma;l++) a[l]=atry[l];
	} 
	else {
	  *alamda *= 10.0;
	  *chisq=ochisq;
	}
      }
    }
  }
  
  // If we encountered an error, make sure any allocated memory is
  // freed before we exit.

  if(waserr) {
    free_matrix(oneda);
    if(da)
      free(da);
    if(beta)
      free(beta);
    if(atry)
      free(atry);
  }
  return waserr;
}

/*.......................................................................
 * Allocate a two-dimensional float array.  The entire matrix can be 
 * deleted by deleting the first row.
 */
float** Fitter::matrix(long nrow, long ncol)
{
  int i,waserr=0;
  float **m=0,*temp=0;
  
  // First allocate the array.

  if((m=(float **) malloc((size_t)((nrow)*sizeof(float*))))==0) {
    fprintf(stderr,"Unable to allocate matrix.\n");
    waserr = 1;
  }
  
  // Next allocate the columns as one block.

  if(!waserr) {
    if((temp=(float *)malloc(sizeof(float)*nrow*ncol))==0) {
      fprintf(stderr,"Unable to allocate columns.\n");
      return 0;
    }
    if(!waserr)
      for(i=0;i < nrow;i++) 
	m[i] = &temp[i*ncol];
  }
  
  // Delete the (possibly) allocated memory if there was an error.

  if(waserr) {
    if(temp)
      free(temp);
    if(m)
      free(m);
  }
  return waserr ? 0 : m;
}
/**.......................................................................
 * Free a two-dimensional float array allocated by matrix().  
 * The matrix should be deleted by deleting the first row only.
 */
float** Fitter::free_matrix(float **m)
{
  // Check that the array hasn't already been deleted.

  if(m==0)
    return m;

  if(m[0])
    free(m[0]);
  free(m);

  return (float **) 0;
}
/*.......................................................................
 * Reorder the parameter array.
 */
void Fitter::covsrt(float **covar, int ma, int ia[], int mfit)
{
  int i,j,k;
  float temp;
  
  for (i=mfit;i<ma;i++)
    for (j=0;j<i;j++) covar[i][j]=covar[j][i]=0.0;
  k=mfit-1;
  for (j=ma-1;j>=0;j--) {
    if (ia[j]) {
      for (i=0;i<ma;i++) 
	SWAP(covar[i][k],covar[i][j]);
      for (i=0;i<ma;i++) 
	SWAP(covar[k][i],covar[j][i]);
      k--;
    }
  }
  return;
}

/*.......................................................................
 * Allocate a new Fit_scan descriptor
 */
Fitter::Fit_scan::Fit_scan(int n)
{
  x.resize(n);
  x2.resize(n);
  y.resize(n);
  sd.resize(n);
  npt = n;
}

/*.......................................................................
 * Delete a Fit_scan descriptor.
 */
Fitter::Fit_scan::~Fit_scan() {}

/**.......................................................................
 * This function assumes: a[0] = x0
 *                        a[1] = sig
 *                        a[2] = amp
 */
FITTER_EVAL_FN(Fitter::gauss_nobase)
{
  float arg,fac,ex;
  
  //  COUT("Inside gauss nobase with a[0] = " << a[0] << " a[1] = " << a[1] << " a[2] = " << a[2]);

  float amp = a[2];
  float sig = a[1];
  float x0  = a[0];

  float dx = x - x0;

  ex  = amp * exp(-(dx*dx)/(2*sig*sig));
  *y = ex;
  
  // df/dx0

  dyda[0] = ex * dx / (sig*sig);

  // df/dsig

  dyda[1] = ex * dx*dx / (sig*sig*sig);

  // df/damp

  dyda[2] = ex / amp;
  
  return;
}

/**.......................................................................
 * This function assumes: a[0] = amp
 *                        a[1] = x0
 *                        a[2] = sigx
 *                        a[3] = y0
 *                        a[4] = sigy
 */
FITTER_EVAL_FN2D(Fitter::gauss2d)
{
  //  COUT("Inside gauss nobase with a[0] = " << a[0] << " a[1] = " << a[1] << " a[2] = " << a[2]);

  float amp  = a[0];
  float x0   = a[1];
  float sigx = a[2];
  float y0   = a[3];
  float sigy = a[4];

  float dx = x1-x0;
  float dy = x2-y0;

  float xfac = exp(-(dx*dx)/(2*sigx*sigx));
  float yfac = exp(-(dy*dy)/(2*sigy*sigy));

  float ex = amp * xfac * yfac;
  *y = ex;
  
  // df/damp

  dyda[0] = ex / amp;

  // df/dx0

  dyda[1] = ex * dx/(sigx*sigx);

  // df/dsigx

  dyda[2] = ex * dx*dx/(sigx*sigx*sigx);

  // df/dy0

  dyda[3] = ex * dy/(sigy*sigy);

  // df/dsigy

  dyda[4] = ex * dy*dy/(sigy*sigy*sigy);
  
  return;
}

/**.......................................................................
 * This function assumes: a[0] = x0
 *                        a[1] = sig
 *                        a[2] = amp
 */
FITTER_EVAL_FN(Fitter::gausstest)
{
#if 0
  float arg,fac,ex;
  
  //  COUT("Inside gauss nobase with a[0] = " << a[0] << " a[1] = " << a[1] << " a[2] = " << a[2]);

  float dx = x-a[0];
  float arg = -(dx*dx)/(2*a[1]*a[1]);
  *y = a[2] * exp(arg);

  dyda[1] = ex*(x-a[0])/(a[1]*a[1]);
  dyda[2] = ex*(powf((x-a[0]),2.0)/powf(a[1],3.0)-1/a[1]);
  dyda[3] = ex/a[2];
#endif
  return;
}

FITTER_EVAL_FN(Fitter::poly)
{
  int i;
  float sum = 0.0;
  for(i=0;i < npar;i++) { 
    sum += a[i]*((i == 0) ? 1.0 : powf(x,(float)(i)));
    dyda[i] = ((i == 0) ? 1.0 : powf(x,(float)(i)));
  }
  *y = sum;

  return;
}

FITTER_EVAL_FN(Fitter::powspec)
{
  // This defines the order of the parameters

  double P0    = a[0]; 
  double ks    = a[1]; 
  double alpha = a[2];

  double r    = ks/x;
  double f    = powf(r, alpha);
  double fac  = 1.0/(1 + f);
  double dpdf = P0 * fac * (1.0 - f*fac);

  double dfdks    = f * (alpha/ks);
  double dfdalpha = f * log(r);
  
  // Evaluate the function

  *y = P0 * f * fac;

  // Now evaluate the derivatives

  dyda[0] = f * fac;
  dyda[1] = dpdf * dfdks; 
  dyda[2] = dpdf * dfdalpha; 

  return;
}

FITTER_EVAL_FN(Fitter::power)
{
  // This defines the order of the parameters

  double norm  = a[0]; 
  double alpha = a[1]; 

  // Evaluate the function

  *y = norm * powf(x, alpha);

  // Now evaluate the derivatives

  dyda[0] = powf(x, alpha);
  dyda[1] = (*y) * log(x);

  //  COUT("Inside power with norm = " << std::setprecision(12) << norm << " alpha = " << std::setprecision(12) << alpha);
  return;
}

/*.......................................................................
 * Based on the NR routine for gauss-jordan elimination.
 */
int Fitter::gaussj(float **a, int n, float **b, int m)
{
  int *indxc=NULL,*indxr=NULL,*ipiv=NULL;
  int i,icol,irow,j,k,l,ll;
  float big,dum,pivinv,temp;
  int waserr=0;
  
  waserr = (indxc=(int *)malloc(sizeof(int)*n))==NULL;
  waserr |= (indxr=(int *)malloc(sizeof(int)*n))==NULL;
  waserr |= (ipiv=(int *)malloc(sizeof(int)*n))==NULL;
  
  if(!waserr) {
    for (j=0;j<n;j++) ipiv[j]=0;
    for (i=0;!waserr && i < n;i++) {
      big=0.0;
      for (j=0;!waserr && j < n;j++)
	if (ipiv[j] != 1)
	  for (k=0; !waserr && k < n;k++) {
	    if (ipiv[k] == 0) {
	      if (fabs(a[j][k]) >= big) {
		big=fabs(a[j][k]);
		irow=j;
		icol=k;
	      }
	    } 
	    else if (ipiv[k] > 1) {
	      //fprintf(stderr,"gaussj: Singular Matrix-1.\n");
	      waserr = 1;
	    }
	  }
      if(!waserr) {
	++(ipiv[icol]);
	if (irow != icol) {
	  for(l=0;l < n;l++)  SWAP(a[irow][l],a[icol][l]);
	  for(l=0;l < m;l++) SWAP(b[irow][l],b[icol][l]);
	}
	indxr[i]=irow;
	indxc[i]=icol;
	if (a[icol][icol] == 0.0) {
	  //	  fprintf(stderr,"gaussj: Singular Matrix-2.\n");
	  waserr = 1;
	}
	
	if(!waserr) {
	  pivinv=1.0/a[icol][icol];
	  a[icol][icol]=1.0;
	  for (l=0;l<n;l++) a[icol][l] *= pivinv;
	  for (l=0;l<m;l++) b[icol][l] *= pivinv;
	  for (ll=0;ll<n;ll++)
	    if (ll != icol) {
	      dum=a[ll][icol];
	      a[ll][icol]=0.0;
	      for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
	      for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
	    }
	}
      }
    }
    if(!waserr) {
      for (l=n-1;l>=0;l--) {
	if (indxr[l] != indxc[l])
	  for (k=0;k<n;k++)
	    SWAP(a[k][indxr[l]],a[k][indxc[l]]);
      }
    }
  }
  /*
   * Free any allocated memory
   */  
  if(ipiv)
    free(ipiv);
  if(indxr)
    free(indxr);
  if(indxc)
    free(indxc);
  return waserr;
}

/**.......................................................................
 * Return the indices into the power spectrum that correspond to the
 * iInterval th log interval of fractional width frac
 */
void Fitter::getIndices(unsigned n, unsigned iInterval, double frac, int& iStart, int& iStop)
{
  //------------------------------------------------------------
  // Get the total log interval spanned by this array
  //------------------------------------------------------------

  double lnMin = log(2*M_PI/n);
  double lnMax = log((n/2+1)*2*M_PI/n);

  double dLnTotal = lnMax - lnMin;

  unsigned nInterval = (unsigned)floor(1.0/frac);

  if(iInterval > nInterval-1)
    ThrowError("You have requested the " << iInterval << "th interval of fractional width " << frac << " but there can be only " << nInterval);

  //------------------------------------------------------------
  // Now get the midpoint of the iInterval th fractional interval
  //------------------------------------------------------------

  double intLnMin = lnMin + iInterval * dLnTotal*frac;
  double intLnMax = intLnMin + dLnTotal*frac;
  double xMin = exp(intLnMin);
  double xMax = exp(intLnMax);

  //  COUT("xMin = " << xMin << " xMax = " << xMax);

  iStart = xMin/(2*M_PI) * n;
  iStop  = xMax/(2*M_PI) * n;

  if(iStart < 0)
    iStart = 0;

  if(iStop > n/2-1)
    iStop = n/2-1;
}


Fitter::Fit_scan Fitter::getPowerSpectrumFitDataNew(std::vector<double>& yarr, std::vector<unsigned>* multiplicity, unsigned nMax,
						    double& p0, double& ks, double& alpha, unsigned& nSamp)
{
  int subChainLen = 2000;

  unsigned nData = (nMax < yarr.size() ? nMax : yarr.size());

  //------------------------------------------------------------
  // First get the power spectrum of the chain, normalized to unit
  // variance (hence the divide by the rms).  We also divide by
  // sqrt(N) so that the power spectrum returns the variance, and not
  // N*Var
  //------------------------------------------------------------

  std::vector<double>norm = yarr;
  double rms = Stats::rms(yarr, multiplicity, 0, nData-1);

  unsigned nSubChain;

  if(nData > subChainLen) {
    nSamp = subChainLen;
    nSubChain = (unsigned) floor(nData / subChainLen);
  } else {
    nSamp = nData;
    nSubChain = 1;
  }

  double nsqrt = sqrt((double)nSamp);

  Dft1d dft(nSamp, true);
  dft.setAverage(true);
  dft.setVectorAverage(false);

  //------------------------------------------------------------
  // Accumulate the power spectrum (possibly over multiple subchains)
  //------------------------------------------------------------

  for(unsigned iSubChain=0; iSubChain < nSubChain; iSubChain++)
    for(unsigned iSamp=0; iSamp < nSamp; iSamp++) 
      dft.pushSample(yarr[iSubChain * nSamp + iSamp]/(nsqrt*rms));

  //------------------------------------------------------------
  // Now get the (possibly averaged) power spectrum, and create a
  // Fit_scan object from it
  //------------------------------------------------------------

  std::vector<double> y = dft.abs2();
  unsigned n = y.size();
  unsigned nFit = n-1;

  Fitter::Fit_scan scan(nFit);

  try {
    //------------------------------------------------------------
    // Estimate sigma from the first third (logarithmic interval) of the
    // samples
    //------------------------------------------------------------

    int iPowStart, iPowStop;
    getIndices(nSamp, 0, 1.0/3, iPowStart, iPowStop);

    double sig = Stats::rms(y, 0, iPowStart, iPowStop);

    for(unsigned i=0; i < nFit; i++) {
      scan.x[i]  = (i+1)*2*M_PI/nSamp;
      scan.y[i]  = y[i+1];
      scan.sd[i] = sig;
    }

    //------------------------------------------------------------
    // Estimate the white-noise power from this section of the samples
    // too
    //------------------------------------------------------------

    p0 = Stats::mean(y, 0, iPowStart, iPowStop);

    //------------------------------------------------------------
    // Estimate the power-law drop off from the last two thirds of the
    // samples
    //------------------------------------------------------------

    int iPowStart1, iPowStop1;
    int iPowStart2, iPowStop2;

    getIndices(nSamp, 1, 1.0/3, iPowStart1, iPowStop1);
    getIndices(nSamp, 2, 1.0/3, iPowStart2, iPowStop2);

    double p2 = Stats::mean(scan.y, 0, iPowStart2, iPowStop2);
    double p1 = Stats::mean(scan.y, 0, iPowStart1, iPowStop1);

    double k2 = Stats::mean(scan.x, 0, iPowStart2, iPowStop2);
    double k1 = Stats::mean(scan.x, 0, iPowStart1, iPowStop1);

    alpha = log(p1/p2) / log(k2/k1);

    //------------------------------------------------------------
    // Finally, estimate k* by walking in intervals along the power
    // spectrum until the power differs by more than sigma/3
    //------------------------------------------------------------

    double nInt = (nSubChain > 1 ? 6 : 5);
    double p;

    for(unsigned i=0; i < nInt; i++) {
      double frac = 1.0/nInt;
      getIndices(nSamp, i, frac, iPowStart, iPowStop);

      p = Stats::mean(scan.y, 0, iPowStart, iPowStop);

      if(p0/p > 1.5)
	break;
    }

    ks = iPowStart * (2*M_PI)/nSamp;

  } catch(...) {
  }

  return scan;
  
}

Fitter::Fit_scan Fitter::getPowerSpectrumFitData(std::vector<double>& yarr, std::vector<unsigned>* multiplicity, 
						 double& p0, double& ks, double& alpha, unsigned& nSamp)
{
  int subChainLen = 2000;

  //------------------------------------------------------------
  // First get the power spectrum of the chain, normalized to unit
  // variance (hence the divide by the rms).  We also divide by
  // sqrt(N) so that the power spectrum returns the variance, and not
  // N*Var
  //------------------------------------------------------------

  std::vector<double>norm = yarr;
  double rms = Stats::rms(yarr, multiplicity, 0, yarr.size()-1);

  unsigned nSubChain;

  if(yarr.size() > subChainLen) {
    nSamp = subChainLen;
    nSubChain = (unsigned) floor(yarr.size() / subChainLen);
  } else {
    nSamp = yarr.size();
    nSubChain = 1;
  }

  double nsqrt = sqrt((double)nSamp);

  Dft1d dft(nSamp, true);
  dft.setAverage(true);
  dft.setVectorAverage(false);

  //------------------------------------------------------------
  // Accumulate the power spectrum (possibly over multiple subchains)
  //------------------------------------------------------------

  for(unsigned iSubChain=0; iSubChain < nSubChain; iSubChain++)
    for(unsigned iSamp=0; iSamp < nSamp; iSamp++) 
      dft.pushSample(yarr[iSubChain * nSamp + iSamp]/(nsqrt*rms));

  //------------------------------------------------------------
  // Now get the (possibly averaged) power spectrum, and create a
  // Fit_scan object from it
  //------------------------------------------------------------

  std::vector<double> y = dft.abs2();
  unsigned n = y.size();
  unsigned nFit = n-1;

  Fitter::Fit_scan scan(nFit);

  try {
    //------------------------------------------------------------
    // Estimate sigma from the first third (logarithmic interval) of the
    // samples
    //------------------------------------------------------------

    int iPowStart, iPowStop;
    getIndices(nSamp, 0, 1.0/3, iPowStart, iPowStop);

    double sig = Stats::rms(y, 0, iPowStart, iPowStop);

    for(unsigned i=0; i < nFit; i++) {
      scan.x[i]  = (i+1)*2*M_PI/nSamp;
      scan.y[i]  = y[i+1];
      scan.sd[i] = sig;
    }

    //------------------------------------------------------------
    // Estimate the white-noise power from this section of the samples
    // too
    //------------------------------------------------------------

    p0 = Stats::mean(y, 0, iPowStart, iPowStop);

    //------------------------------------------------------------
    // Estimate the power-law drop off from the last two thirds of the
    // samples
    //------------------------------------------------------------

    int iPowStart1, iPowStop1;
    int iPowStart2, iPowStop2;

    getIndices(nSamp, 1, 1.0/3, iPowStart1, iPowStop1);
    getIndices(nSamp, 2, 1.0/3, iPowStart2, iPowStop2);

    double p2 = Stats::mean(scan.y, 0, iPowStart2, iPowStop2);
    double p1 = Stats::mean(scan.y, 0, iPowStart1, iPowStop1);

    double k2 = Stats::mean(scan.x, 0, iPowStart2, iPowStop2);
    double k1 = Stats::mean(scan.x, 0, iPowStart1, iPowStop1);

    alpha = log(p1/p2) / log(k2/k1);

    //------------------------------------------------------------
    // Finally, estimate k* by walking in intervals along the power
    // spectrum until the power differs by more than sigma/3
    //------------------------------------------------------------

    double nInt = (nSubChain > 1 ? 6 : 5);
    double p;

    for(unsigned i=0; i < nInt; i++) {
      double frac = 1.0/nInt;
      getIndices(nSamp, i, frac, iPowStart, iPowStop);

      p = Stats::mean(scan.y, 0, iPowStart, iPowStop);

      if(p0/p > 1.5)
	break;
    }

    ks = iPowStart * (2*M_PI)/nSamp;

  } catch(...) {
  }

  return scan;
  
}
