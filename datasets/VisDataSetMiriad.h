// $Id: VisDataSetMiriad.h,v 1.3 2012/03/15 01:19:32 eml Exp $

#ifndef GCP_DATASETUVFS_VISDATASETMIRIAD_H
#define GCP_DATASETUVFS_VISDATASETMIRIAD_H

/**
 * @file VisDataSet.h
 * 
 * Tagged: Wed Jun 16 10:24:00 PDT 2010
 * 
 * @version: $Revision: 1.3 $, $Date: 2012/03/15 01:19:32 $
 * 
 * @author tcsh: Erik Leitch
 */
#include "gcp/datasets/VisDataSet.h"

#include "gcp/fftutil/ObsInfo.h"

#include "gcp/util/MiriadIo.h"

#include <map>
#include <vector>

namespace gcp {
  namespace datasets {

    class VisDataSetMiriad : public gcp::datasets::VisDataSet {
    public:

      /**
       * Constructor.
       */
      VisDataSetMiriad(gcp::util::ThreadPool* pool=0);

      /**
       * Destructor.
       */
      virtual ~VisDataSetMiriad();

      // Initialize this data set from a file.

      void getGroup(unsigned iGroup, gcp::util::ObsInfo::Vis& vis);

      void openFileReader(std::string fileName);
      void closeFileReader();

      void initializeAntennaInformation(std::string fileName);
  

      void updateFrequencyInformation();
      void updateObservationInformation();
      void updateVisibilityInformation();

    private:

      gcp::util::MiriadIo reader_;

      virtual void initializeFileReader(std::string fileName);

    }; // End class VisDataSetMiriad

  } // End namespace datasets
} // End namespace gcp



#endif // End #ifndef GCP_DATASETUVFS_VISDATASETMIRIAD_H
