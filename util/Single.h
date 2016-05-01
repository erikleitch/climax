// $Id: $

#ifndef GCP_UTIL_SINGLE_H
#define GCP_UTIL_SINGLE_H

/**
 * @file Single.h
 * 
 * Tagged: Wed Feb  4 10:23:31 PST 2015
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author username: Command not found.
 */
namespace gcp {
  namespace util {

    class Single {
    public:

      static bool createInstance();

      static void printHello();

    private:

      static int count_;
      //      static bool init_;

      /**
       * Constructor.
       */
      Single();

      /**
       * Destructor.
       */
      virtual ~Single();

    }; // End class Single

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_SINGLE_H
