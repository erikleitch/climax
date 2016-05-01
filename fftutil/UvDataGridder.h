// $Id: UvDataGridder.h,v 1.4 2012/05/31 16:22:56 eml Exp $

#ifndef GCP_UTIL_UVDATAGRIDDER_H
#define GCP_UTIL_UVDATAGRIDDER_H

/**
 * @file UvDataGridder.h
 * 
 * Tagged: Thu Jun  3 11:41:14 PDT 2010
 * 
 * @version: $Revision: 1.4 $, $Date: 2012/05/31 16:22:56 $
 * 
 * @author Erik Leitch
 */
#include "gcp/fftutil/Dft2d.h"

#define DO_INTERP 0

namespace gcp {
  namespace util {

    class UvDataGridder : public Dft2d {
    public:

      UvDataGridder(bool optimize=true);
      virtual ~UvDataGridder();

      UvDataGridder(const UvDataGridder& gridder);
      UvDataGridder(UvDataGridder& gridder);

      void operator=(const UvDataGridder& gridder);
      void operator=(UvDataGridder& gridder);

      void assignDataFrom(const UvDataGridder& gridder);
      void assignDataFrom(UvDataGridder& gridder);

      void sizeToMatch(UvDataGridder& gridder);
      void duplicate(UvDataGridder& gridder);

      void operator+=(const UvDataGridder& gridder);
      void operator+=(UvDataGridder& gridder);

      void mergeData(UvDataGridder& gridder);

      virtual void operator*=(double mult);

      void estimateErrorInMeanFromData(bool estimate);

      // Reinitialize this object for first-moment accumulation

      void initializeForFirstMoments();

      // Accumulate the weighted mean of visibilities

      void accumulateFirstMoments(double u, double v, double re, double im, double wt);
      void accumulateFirstMomentsTest(double u, double v, double re, double im, double wt);
      void accumulateFirstMomentsWithInterpolation(double u, double v, double re, double im, double wt);
      void accumulateFirstMoments(unsigned dftInd, double re, double im, double wt);

      // Reinitialize this object for second-moment accumulation

      void initializeForSecondMoments();

      // Accumulate moments needed for calculation of the error in the
      // mean

      void accumulateSecondMoments(double u, double v, double re, double im, double wt);
      void accumulateSecondMomentsWithInterpolation(double u, double v, double re, double im, double wt);

      // After visibilities have been added using addVis(), we will
      // have stored weighted first and second moments.  Call
      // calculateErrorInMean() to convert those moments into an error
      // in the mean.

      void calculateErrorInMean();

      std::vector<double> getPopulatedIndices();

      // Explicitly zero the data array

      void zero();

      void assignPopulatedIndicesFrom(UvDataGridder& gridder);

      // Return an estimate of the map pixel rms

      double estimateMapPixelRms();

      // Initialize populated indices to all indices

      void initializePopulatedIndicesToAll();

      void setPopulatedIndicesToValue(double val);

      void plotOccupiedUV();
      void plotWt();

      void renormalize();

      void debugPrint(bool print);

      bool debugPrint_;

    private:

      // Reset arrays for dft of a different size

      void resize();

      bool estimateErrInMeanFromData_;

      // The number of points that enter into each point of the UV
      // grid

    public:

      std::vector<unsigned> nPt_;
      double wtSumTotal_;

      double uuSum_;
      double vvSum_;
      double uvSum_;

      // The sum of weights for each point in the UV grid

      std::vector<double> wtSum_;

      // The sum of squared weights for each point in the UV grid

      std::vector<double> wt2Sum_;

      // The weighted mean of the complex visibilities

      fftw_complex* mean_;

      // Overload this from the base class, for reasons explained
      // within the .cc file

      void computeInverseTransform(fftw_plan* invPlan=0);
      void computeInverseTransformNoRenorm(fftw_plan* invPlan=0);

      void copyFromStore();
      void copyToStore();

      void getEstimatedSynthesizedBeam(Angle& fwhmMaj, Angle& fwhmMin, Angle& posAngle);

    public:
      
      // A vector that stores the list of indices in the mean_ array
      // that actually contain data

      std::vector<unsigned> populatedIndices_;

      // For convenience, we also store the UV coordinate that
      // corresponds to each populated index

      std::vector<double>   populatedU_;
      std::vector<double>   populatedV_;

      // The weighted error in the mean of the complex visibilities

      fftw_complex* errorInMean_;
      fftw_complex* store_;

      // True when errors have been calculated

      bool errorInMeanIsValid_;

    }; // End class UvDataGridder

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_UVDATAGRIDDER_H
