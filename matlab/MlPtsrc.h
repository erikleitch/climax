// $Id: MlPtsrc.h,v 1.1.1.1 2010/07/13 17:56:46 eml Exp $

#ifndef GCP_MATLAB_MLPTSRC_H
#define GCP_MATLAB_MLPTSRC_H

/**
 * @file MlPtsrc.h
 * 
 * Tagged: Fri Apr 17 15:40:33 PDT 2009
 * 
 * @version: $Revision: 1.1.1.1 $, $Date: 2010/07/13 17:56:46 $
 * 
 * @author Erik Leitch.
 */
namespace gcp {
  namespace matlab {

    class MlPtsrc {
    public:

      static void kernelint(unsigned nsrc, double* sobs, double* retVals, double gamma, double sigma, double nsig, unsigned nbin);
      static void dkernelint(unsigned nsrc, double* sobs, double* retVals, double gamma, double sigma, double nsig, unsigned nbin);
      static void d2kernelint(unsigned nsrc, double* sobs, double* retVals, double gamma, double sigma, double nsig, unsigned nbin);
      static void kernelparams(double& peak, double& width, double sobs, double gamma, double sigma);

      /**
       * Constructor.
       */
      MlPtsrc();

      /**
       * Destructor.
       */
      virtual ~MlPtsrc();

    private:
    }; // End class MlPtsrc

  } // End namespace matlab
} // End namespace gcp



#endif // End #ifndef GCP_MATLAB_MLPTSRC_H
