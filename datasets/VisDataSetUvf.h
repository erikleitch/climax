// $Id: VisDataSetUvf.h,v 1.3 2012/03/15 01:19:32 eml Exp $

#ifndef GCP_DATASETUVFS_VISDATASETUVF_H
#define GCP_DATASETUVFS_VISDATASETUVF_H

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

#include "gcp/util/ThreadPool.h"

#include <map>
#include <vector>

namespace gcp {
  namespace datasets {

    class VisDataSetUvf : public gcp::datasets::VisDataSet {
    public:

      /**
       * Constructor.
       */
      VisDataSetUvf(gcp::util::ThreadPool* pool=0);

      /**
       * Destructor.
       */
      virtual ~VisDataSetUvf();

      // Inherited interface from DataSet

      void writeCompositeModelToFile(std::string fileName, double sigma);

    private:

      gcp::util::ThreadPool* pool_;
      gcp::util::FitsUvfReader reader_;
      std::vector<double> ifOffsetHz_;
      std::vector<double> bwHz_;

      void getGroup(unsigned iGroup, gcp::util::ObsInfo::Vis& vis);

      void openFileReader(std::string fileName);
      void closeFileReader();
      void updateObservationInformation(gcp::util::FitsUvfReader& reader);
      void updateObservationInformation(gcp::datasets::VisDataSet* dataset);

      void initializeAntennaInformation(std::string fileName);
      void initializeFrequencyInformation(std::string fileName);

      void updateFrequencyInformation();
      void updateObservationInformation();
      void updateVisibilityInformation();

    public:

      void addDataSet(std::string file);

    }; // End class VisDataSetUvf

  } // End namespace datasets
} // End namespace gcp



#endif // End #ifndef GCP_DATASETUVFS_VISDATASETUVF_H
