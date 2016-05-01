#include <iostream>
#include <fstream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/models/GaussianClusterModel.h"

#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/IoLock.h"
#include "gcp/util/GaussianVariate.h"
#include "gcp/util/Timer.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/fftutil/DataSet.h"
#include "gcp/fftutil/Dft2d.h"
#include "gcp/fftutil/Image.h"
#include "gcp/fftutil/Model.h"

#include "gcp/util/String.h"

#include "gcp/models/ModelManager.h"
#include "gcp/datasets/DataSetManager.h"

using namespace std;
using namespace gcp::datasets;
using namespace gcp::util;
using namespace gcp::program;
using namespace gcp::models;

KeyTabEntry Program::keywords[] = {
  { "file",     "",               "s", "FITS file to read in"},
  { "nburn",    "1000",           "i", "Number of burn-in samples"},
  { "ntry",     "10000",          "i", "Number of trials"},

  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

void parseFile(std::string modelFile);
void getDataSetNameAndType(String& line, std::string& name, std::string& type);
void getModelNameAndType(String& line, std::string& name, std::string& type);
void getTokVal(String& line, String& tok, String& val);
void parseCorrelationVariables(String& token, std::string& varName1, std::string& varName2);
void parseVarname(String& mVarName, String& modelName, String& varName);
void parseUniformPrior(String& val, String& min, String& max, String& units);
void parseGaussianPrior(String& val, String& mean, String& sigma, String& units);
void parseVal(String& val, String& value, String& units);
void parseDataSetVariableAssignment(DataSet* dataSet, String& varName, String& valStr);
void parseModelVariableAssignment(String& tok, String& valStr);
void parseVariableAssignment(String& tok, String& val);
void displayModel();
void checkIfNameAlreadyExists(std::string name);
void runMarkov(unsigned nBurn, unsigned nTry);
void updateHessian(ChisqVariate& chisq);

ModelManager   mm_;
DataSetManager dm_;

void generate1DSample();
void generate2DSample();
void multivariateSample();
void multivariateSample2();

int Program::main()
{
  multivariateSample2();
  return 0;
}
		 
void generate1DSample()
{
  double logpFirst,  logpNext,  logpCurr;
  double dlogpFirst, dlogpNext, dlogpCurr, dlogpLast;
  double xFirst, xNext, xCurr;
  double d2;
  double x1, x2;

  for(unsigned i=0; i < 10; i++) {

    xFirst   = xNext;
    xNext    = xCurr;
    xCurr    = Sampler::generateGaussianSample(0.2);

    logpFirst = logpNext;
    logpNext  = logpCurr;
    logpCurr  = Sampler::lnGaussPdf(xCurr, 0.0, 0.2);

    
    x2 = (xNext + xFirst)/2;
    x1 = (xCurr + xNext)/2;

    dlogpFirst = (logpNext - logpFirst) / (xNext - xFirst);
    dlogpLast  = (logpCurr - logpNext)  / (xCurr - xNext);

    d2 = (dlogpLast - dlogpFirst) / (x1 - x2);

    COUT("xCurr = " << xCurr << " logpCurr = " << logpCurr  << " d2 = " << d2);
  }

}  

void generate2DSample()
{
  unsigned nVar = 2;

  Vector<double> xCurr(nVar);
  xCurr[0] = 0.1; xCurr[1] = 0.1;

  Vector<double> mean(nVar);
  mean[0] = 0.0; mean[1] = 0.0;

  Matrix<double> sigma = Matrix<double>::identity(nVar);
  sigma[0][0] = 1.0;
  sigma[1][1] = 1.0;

  Matrix<double> correl = Matrix<double>::identity(nVar);

  Matrix<double> cov = sigma * correl * sigma;
  Matrix<double> invCov = cov.inverse();
  double detC = cov.determinant();
  
  double logpFirst,  logpNext,  logpCurr;

  Vector<double> dlogpFirst(nVar), dlogpNext(nVar), dlogpCurr(nVar), dlogpLast(nVar);
  Vector<double> xFirst(nVar), xNext(nVar);
  Vector<double> d2(nVar);
  Vector<double> x1(nVar), x2(nVar);

  COUT("cov = " << std::endl << cov);

  for(unsigned i=0; i < 10; i++) {

    xFirst   = xNext;
    xNext    = xCurr;
    //    xCurr    = Sampler::generateMultiVariateGaussianSample(xCurr, mean, invCov);
    xCurr[0] = Sampler::generateGaussianSample(1.0);
    xCurr[1] = Sampler::generateGaussianSample(1.0);

    logpFirst = logpNext;
    logpNext  = logpCurr;
    //    logpCurr  = Sampler::lnMultiVariateGaussPdf(xCurr, mean, invCov, detC);
    logpCurr  = Sampler::lnGaussPdf(xCurr[0], 0.0, 1.0) + Sampler::lnGaussPdf(xCurr[1], 0.0, 1.0);

    for(unsigned iVar=0; iVar < nVar; iVar++) {

      x2[iVar] = (xNext[iVar] + xFirst[iVar])/2;
      x1[iVar] = (xCurr[iVar] + xNext[iVar])/2;

      dlogpFirst[iVar] = (logpNext - logpFirst) / (xNext[iVar] - xFirst[iVar]);
      dlogpLast[iVar]  = (logpCurr - logpNext)  / (xCurr[iVar] - xNext[iVar]);

      d2[iVar] = (dlogpLast[iVar] - dlogpFirst[iVar]) / (x1[iVar] - x2[iVar]);
    }

    COUT("xCurr = " << xCurr << " logpCurr = " << logpCurr  << " d2 = " << d2);
  }

}  

void multivariateSample()
{
  unsigned nVar = 3;
  Vector<double> mean(nVar);
  Vector<double> val(nVar);
  Matrix<double> sigma = Matrix<double>::identity(nVar);

  for(unsigned iVar=0; iVar < nVar; iVar++) {
    val[iVar]  = 0.0;
    mean[iVar] = 0.0;
    if(iVar==0)
      sigma[iVar][iVar] = 1.0;
    if(iVar==1)
      sigma[iVar][iVar] = 0.2;
    if(iVar==2)
      sigma[iVar][iVar] = 0.5;
  }

  Matrix<double> correl = Matrix<double>::identity(nVar);
  Matrix<double> cov = sigma * correl * sigma;
  Matrix<double> invCov = cov.inverse();
  double detC = cov.determinant();

  // Ok -- now we're set up to generate multivariate samples. Now
  // generate samples in turn from each variate's conditional
  // distribution and use them to calculate the Hessian matrix

  val = Sampler::generateMultiVariateGaussianSample(val, mean, invCov);

  // And store it as the starting point for each variate 

  Vector<double> storedVal = val;

  Vector<double> x0(nVar), x1(nVar), x2(nVar), tmp(nVar);
  double lp0, lp1, lp2, dlp1, dlp2, d2;
  double val1, val2;

  for(unsigned iVar=0; iVar < nVar; iVar++) {

    // Generate a new random point from this variate's conditional
    // distribution

    tmp = storedVal;

    x0 = tmp;
    Sampler::multiVariateSampleIterator(iVar, tmp, mean, invCov);
    x1 = tmp;
    Sampler::multiVariateSampleIterator(iVar, tmp, mean, invCov);
    x2 = tmp;

    lp0 = Sampler::lnMultiVariateGaussPdf(x0, mean, invCov, detC);
    lp1 = Sampler::lnMultiVariateGaussPdf(x1, mean, invCov, detC);
    lp2 = Sampler::lnMultiVariateGaussPdf(x2, mean, invCov, detC);
    
    val1 = (x1[iVar] + x0[iVar])/2;
    val2 = (x2[iVar] + x1[iVar])/2;

    // Calculate dlogp/dx at each point

    dlp1 = (lp1 - lp0) / (x1[iVar] - x0[iVar]);
    dlp2 = (lp2 - lp1) / (x2[iVar] - x1[iVar]);

    // Finally, calculate d2logp/dx2

    d2 = (dlp2 - dlp1) / (val2 - val1);

    COUT("ivar = " << iVar << " d2 = " << d2);
  }
}



void multivariateSample2()
{
  unsigned nVar = 3;
  Vector<double> mean(nVar);
  Vector<double> val(nVar);
  Matrix<double> sigma = Matrix<double>::identity(nVar);

  for(unsigned iVar=0; iVar < nVar; iVar++) {
    val[iVar]  = 0.0;

    if(iVar==0) {
      sigma[iVar][iVar] = 1.0;
      mean[iVar] = 0.0;
    }

    if(iVar==1) {
      sigma[iVar][iVar] = 0.2;
      mean[iVar] = 1.0;
    }

    if(iVar==2) {
      sigma[iVar][iVar] = 0.5;
      mean[iVar] = 2.0;
    }
  }

  Matrix<double> correl = Matrix<double>::identity(nVar);
  Matrix<double> cov = sigma * correl * sigma;
  Matrix<double> invCov = cov.inverse();
  double detC = cov.determinant();

  // Ok -- now we're set up to generate multivariate samples. Now
  // generate samples in turn from each variate's conditional
  // distribution and use them to calculate the Hessian matrix

  //  val = Sampler::generateMultiVariateGaussianSample(val, mean, invCov);
  val[0] = 80;
  val[1] = 90000;
  val[2] = 0.9;

  // And store it as the starting point for each variate 

  Vector<double> storedVal = val;

  Vector<double> xlo(nVar), xmid(nVar), xhi(nVar), tmp(nVar);
  double lplo, lpmid, lphi, dlp0, dlp1, dlp2, d2;
  double val1, val2;

  unsigned iIter=0, nIter=10;
  bool stop = false;

  do {

    ++iIter;

    if(iIter > nIter)
      stop = true;

    for(unsigned iVar=0; iVar < nVar; iVar++) {
      
      // Generate a new random point from this variate's conditional
      // distribution
      
      tmp = val;
      
      xlo  = tmp;
      xmid = tmp;
      xhi  = tmp;
      
      xlo[iVar] = xmid[iVar] - 0.1*sigma[iVar][iVar]; 
      xhi[iVar] = xmid[iVar] + 0.1*sigma[iVar][iVar]; 
      
      lplo  = Sampler::lnMultiVariateGaussPdf(xlo,  mean, invCov, detC);
      lpmid = Sampler::lnMultiVariateGaussPdf(xmid, mean, invCov, detC);
      lphi  = Sampler::lnMultiVariateGaussPdf(xhi,  mean, invCov, detC);
      
      val1 = (xmid[iVar] +  xlo[iVar])/2;
      val2 = (xhi[iVar]  + xmid[iVar])/2;
      
      // Calculate dlogp/dx at each point
      
      dlp1 = (lpmid - lplo) / (xmid[iVar] -  xlo[iVar]);
      dlp2 = (lphi - lpmid) / (xhi[iVar]  - xmid[iVar]);
      
      // Estimate dlogp/dx at the midpoint by interpolating between
      // these two
      
      dlp0 = dlp1 + (dlp2 - dlp1) / (val2 - val1) * (xmid[iVar] - val1);
      
      // Finally, calculate d2logp/dx2
      
      d2 = (dlp2 - dlp1) / (val2 - val1);
      
      double s2  = -1.0/d2;
      
      // And the estimated mean of the distribution is offset from the
      // current sample position by the first derivative
      
      // Iterate to the estimated mean

      val[iVar] = xmid[iVar] + dlp0 * s2;
      
      COUT("iVar = " << iVar << " est sigma = " << sqrt(s2) << " est mean = " << val[iVar]);
    }
  } while(!stop);
}


