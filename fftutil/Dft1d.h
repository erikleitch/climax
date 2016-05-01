// $Id: emacs_macros,v 1.1.1.1.4.1 2006/01/25 08:01:42 bicep Exp $

#ifndef GCP_UTIL_DFT1D_H
#define GCP_UTIL_DFT1D_H

/**
 * @file Dft1d.h
 * 
 * Tagged: Tue Jan  9 03:18:37 CST 2007
 * 
 * @version: $Revision: 1.1.1.1.4.1 $, $Date: 2006/01/25 08:01:42 $
 * 
 * @author Erik Leitch
 */

#define COMPLEX_AV

#include <iostream>
#include <vector>

#include "gcp/util/CircBuf.h"
#include "gcp/util/Frequency.h"
#include "gcp/util/TimeVal.h"

#include <fftw3.h>

#define APOD_FN(fn) double (fn)(unsigned iSamp, unsigned n)

namespace gcp {
  namespace util {

    class Dft1d {
    public:

      enum Apodization {
	APOD_RECTANGLE = 0, // Rectangular (Dirichlet) window
	APOD_TRIANGLE  = 1, // Triangle apodization (goes to zero at edges)
	APOD_HAMMING   = 2, // Hamming window
	APOD_HANN      = 3, // Hann window
	APOD_COS       = 4, // Cosine window
	APOD_SINC      = 5, // Sinc window
      };

      /**
       * Constructor for a real transform n points long
       */
      Dft1d(int n, bool optimize=true, Apodization apod=APOD_RECTANGLE);

      /**
       * Copy Constructor.
       */
      Dft1d(const Dft1d& objToBeCopied);

      /**
       * Copy Constructor.
       */
      Dft1d(Dft1d& objToBeCopied);

      /**
       * Const Assignment Operator.
       */
      void operator=(const Dft1d& objToBeAssigned);

      /**
       * Assignment Operator.
       */
      void operator=(Dft1d& objToBeAssigned);

      /**
       * Output Operator.
       */
      friend std::ostream& operator<<(std::ostream& os, Dft1d& obj);

      /**
       * Destructor.
       */
      virtual ~Dft1d();

      //-----------------------------------------------------------------------
      // Set methods
      //-----------------------------------------------------------------------

      /**
       * Resize this object to accomodate transforms of length n.  If
       * the pushSample() method is called on this object, it will
       * automatically transform the data when n samples have been
       * pushed into the buffer.
       */
      void resize(unsigned int);

      /**
       * Set the time resolution of the x-axis.  This will be used to
       * calculate frequency span and resolution of the transform.
       */
      void setTimeRes(TimeVal tVal);

      /**
       * If true, we will store on-the-fly averages of the transforms
       */
      void setAverage(bool doAv);

      /**
       * If integration is turned on, set the type of integration to do
       */
      void setVectorAverage(bool vecAv);

      /**
       * Set the apodization type (stages the change until the start of the next transform)
       */
      void setApodizationType(Apodization apod);

      //-----------------------------------------------------------------------
      // Query methods
      //-----------------------------------------------------------------------

      /**
       * Return the size of the input array
       */
      unsigned inputSize();

      /**
       * Return the size of the transformed array
       */
      unsigned transformSize();

      /**
       * Return the size of the transformed array
       */
      static unsigned transformSize(unsigned n);

      /**
       * Get the frequency resolution of the transform
       */
      Frequency getFrequencyResolution();

      // Static method to return the same

      static Frequency getFrequencyResolution(unsigned npt, TimeVal timeRes);

      /**
       * Return the power spectrum
       */
      std::vector<double> powerSpectrum(double* data);

      /**
       * Return the absolute value of the output array
       */
      std::vector<double> abs();
      std::vector<double> abs2();

      /**
       * Return the absolute value of the power spectrum, filling an
       * external array
       */
      void abs(float* outputArray, 
	       bool& first, double& min, double& max, bool doLog=false);

      /**
       * Update the min max
       */
      void updateMinMax(unsigned i, double& val, bool doLog, 
			bool& first, double& min, double& max);

      //-----------------------------------------------------------------------
      // Manual methods for using this class
      //-----------------------------------------------------------------------

      /**
       * Fill the input data array
       */
      void fillInputArray(double* data);

      /**
       * Compute the transform
       */
      void computeTransform();

      /**
       * Return a pointer to the input data
       */
      double* getInputData();

      /**
       * Return a pointer to the transformed data
       */
      fftw_complex* getTransform();

      //-----------------------------------------------------------------------
      // Automatic methods for computing transforms
      //-----------------------------------------------------------------------

      /**
       * Push a sample onto the input array.  If this push fills the
       * input buffer, the transform of the input buffer will be
       * computed.
       */
      void pushSample(double sample); 

      /**
       * Will return true after a new transform has been computed, and
       * before the next call to pushSample()
       */
      bool transformIsReady();

      //-----------------------------------------------------------------------
      // Apodization functions
      //-----------------------------------------------------------------------

      static APOD_FN(apodRectangle);
      static APOD_FN(apodTriangle);
      static APOD_FN(apodHamming);
      static APOD_FN(apodHann);
      static APOD_FN(apodCos);
      static APOD_FN(apodSinc);

    private:

      /**
       * Compute a plan for this fft
       */
      void computePlan(unsigned flag); 

      /**
       * Add a transform to the averaege
       */
      void addToAverage();

      // Remove the mean from the input array

      void removeMean();

      // Assert a pending apodization type

      void assertApodizationType(Apodization apod);

      // Return the apodization coefficient for the current sample

      double apodizationCoefficient(unsigned iSamp);

      void copySamples();

    private:

      // Fundamental variables

      int n_;                       // The size transform we will compute
      TimeVal tRes_;                // The time resolution of the samples
      bool haveRes_;                // True once the resolution is set
      Apodization currentApodType_; // The current apodization type
      Apodization pendingApodType_; // The pending apodization type
      bool apodTypePending_;        // An apodization command was received,
			            // that should be put into effect the next
			            // time the sample counter rolls over
      APOD_FN(*apodFn_);            // The current apodization function

      double mean_;
      bool subtractMean_;

      // FFTW specific variables

      double* in_;          // The input array to be transformed
      fftw_complex* out_;   // The output of the transformed
      fftw_plan plan_;      // Instructions to fftw for the best method
			    // to FFT
      bool optimize_;       // True if we want fftw to perform expensive
			    // tests to determine the optimal FFT
			    // method.  Note that these tests will only
			    // be done on calls to resize()

      // Averaging variables

      unsigned nAv_;        // The number of transforms currently in
			    // the running average
      bool doAv_;           // If true, average the transforms

      bool vecAv_;          // If true, vector-average the transforms.
			    // If false, rms-average the power spectra
      
      bool linear_;         // If true, output data in linear (rather
			    // than log) units

      fftw_complex* outAv_; // The averaged output of the transform
      std::valarray<double> absBuf_; // The average output of the
				     // power spectra

      // Auto-transform variables

      std::valarray<double> sampBuf_; // Buffer into which we store
				      // samples
      std::valarray<double> apodBuf_; // Buffer into which we store
				      // apodization coefficients
#if 0
      CircBuf<double> buf_;     // A circular buffer for storing
			        // on-the-fly samples
      CircBuf<double> apodBuf_; // A circular buffer for storing
			        // on-the-fly samples
#endif
      unsigned nSinceLastTransform_;
      bool transformIsReady_;

      
    }; // End class Dft1d

  } // End namespace util
} // End namespace gcp


#endif // End #ifndef GCP_UTIL_DFT_H
