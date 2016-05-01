// $Id: $

#ifndef GCP_UTIL_TRANS_H
#define GCP_UTIL_TRANS_H

/**
 * @file Trans.h
 * 
 * Tagged: Thu Jan 16 16:58:45 PST 2014
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
namespace gcp {
  namespace util {

    class Trans {
    public:

      /**
       * Constructor.
       */
      Trans();

      /**
       * Destructor.
       */
      virtual ~Trans();

      float dx_;
      float dy_;
      unsigned nx_;
      unsigned ny_;
      unsigned ndata_;
      float xmins_;
      float ymins_;
      float xmaxs_;
      float ymaxs_;
      float* zdata_;
      bool reverseX_;
      bool reverseY_;

      float valNearestToPoint(float x, float y);
      float convolveAroundPoint(float x, float y, float* data);
      void rectify(float& min, float& max);
      void printStats(float x1, float x2, float y1, float y2);

    }; // End class Trans

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_TRANS_H
