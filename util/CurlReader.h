// $Id: CurlReader.h,v 1.1 2009/08/17 23:36:37 eml Exp $

#ifndef GCP_UTIL_CURLREADER_H
#define GCP_UTIL_CURLREADER_H

/**
 * @file CurlReader.h
 * 
 * Tagged: Fri Aug 14 16:23:24 PDT 2009
 * 
 * @version: $Revision: 1.1 $, $Date: 2009/08/17 23:36:37 $
 * 
 * @author tcsh: Erik Leitch
 */
#include <sstream>

namespace gcp {
  namespace util {

    class CurlReader {
    public:

      /**
       * Constructor.
       */
      CurlReader();

      /**
       * Destructor.
       */
      virtual ~CurlReader();

      std::string getUrl(std::string url);

    protected:

      enum {
	OPTION_FALSE = 0,
	OPTION_TRUE  = 1
      };

      // A stream containing the last read information returned from
      // the URL fetch

      std::ostringstream lastRead_;

      // A function that will be called to handle data from a call to
      // 'perform'

      static size_t handleData(void* buffer, size_t size, size_t nmemb, void* userp);

    }; // End class CurlReader

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_CURLREADER_H
