// $Id: OsInfo.h,v 1.2 2012/05/29 16:20:29 eml Exp $

#ifndef GCP_UTIL_OSINFO_H
#define GCP_UTIL_OSINFO_H

/**
 * @file OsInfo.h
 * 
 * Tagged: Fri Jun 11 11:01:39 PDT 2010
 * 
 * @version: $Revision: 1.2 $, $Date: 2012/05/29 16:20:29 $
 * 
 * @author tcsh: Erik Leitch
 */
namespace gcp {
  namespace util {

    class OsInfo {
    public:

      /**
       * Constructor.
       */
      OsInfo();

      /**
       * Destructor.
       */
      virtual ~OsInfo();

      static bool isBigEndian();
      static bool isLittleEndian();

      static void printBits(unsigned int iVal);

      static int getNumberOfCpus(void);

    private:

    }; // End class OsInfo

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_OSINFO_H
