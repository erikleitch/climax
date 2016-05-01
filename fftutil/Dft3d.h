// $Id: Dft3d.h,v 1.1.1.1 2010/07/13 17:56:32 eml Exp $

#ifndef GCP_UTIL_DFT3D_H
#define GCP_UTIL_DFT3D_H

/**
 * @file Dft3d.h
 * 
 * Tagged: Thu Jun  3 11:41:14 PDT 2010
 * 
 * @version: $Revision: 1.1.1.1 $, $Date: 2010/07/13 17:56:32 $
 * 
 * @author Erik Leitch
 */
#include <fftw3.h>

#include "gcp/fftutil/Image.h"

namespace gcp {
  namespace util {

    class Dft3d {
    public:

      // Constructor.

      Dft3d(bool optimize=true);
      Dft3d(int nx, int ny, int nz, bool optimize=true);
      Dft3d(Image& image, bool optimize=true);

      // Destructor.

      virtual ~Dft3d();

      void resize(unsigned int nx, unsigned int ny, unsigned nz);
      void initialize(Image& image);

      // Compute the transforms

      void computeForwardTransform();
      void computeInverseTransform();

    public:

      double* in_;            // The input array to be transformed
      fftw_complex* out_;     // The output of the transformed

      fftw_plan forwardPlan_; // Instructions to fftw for the best forward method
      fftw_plan inversePlan_; // Instructions to fftw for the best inverse method
			      // to FFT
      bool optimize_;         // True if we want fftw to perform expensive
			      // tests to determine the optimal FFT
			      // method.  Note that these tests will only
			      // be done on calls to resize()

      bool normalize_;
      bool zeropad_;

      unsigned nIn_;
      unsigned nOut_;

      unsigned nx_;
      unsigned ny_;
      unsigned nz_;

      // Compute a plan for this fft

      void computePlan(unsigned flag); 
      void initialize();

    }; // End class Dft3d

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_DFT3D_H
