#include "gcp/util/CurlReader.h"
#include "gcp/util/Exception.h"

#include <curl/curl.h>

using namespace std;
using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
CurlReader::CurlReader() {}

/**.......................................................................
 * Destructor.
 */
CurlReader::~CurlReader() {}


std::string CurlReader::getUrl(std::string url)
{
  lastRead_.str("");

  std::ostringstream testStr;

  // global libcURL init

  curl_global_init( CURL_GLOBAL_ALL ) ;

  // create a context, sometimes known as a handle.
  // Think of it as a lookup table, or a source of config data.

  CURL* ctx = curl_easy_init() ;
  
  if(ctx == NULL) {
    ThrowError("Unable to initialize cURL interface");
  }
  
  // BEGIN: configure the handle:
  
  // target url:

  curl_easy_setopt(ctx, CURLOPT_USERAGENT, "curl 7.21.4 (universal-apple-darwin11.0) libcurl/7.21.4 OpenSSL/0.9.8y zlib/1.2.5");

  curl_easy_setopt( ctx , CURLOPT_URL,  url.c_str() ) ;

  // no progress bar:

  curl_easy_setopt( ctx , CURLOPT_NOPROGRESS , OPTION_TRUE) ;

  // (sending response headers to stderr)

  // This next line causes a seg fault when trying to pass data to my user function
  //  curl_easy_setopt( ctx , CURLOPT_WRITEHEADER , stderr) ;

  // response content: same choices as headers
  // send it to stdout to prove that libcURL differentiates the two

  curl_easy_setopt( ctx , CURLOPT_WRITEFUNCTION, handleData);

  curl_easy_setopt( ctx , CURLOPT_WRITEDATA, (void*)this);

  // END: configure the handle 
	
  // action!

  const CURLcode rc = curl_easy_perform( ctx ) ;

  if(rc != CURLE_OK) {

    std::cerr << "Error from cURL: " << curl_easy_strerror( rc ) << std::endl ;
    std::cerr << "Error from cURL: " << std::endl;
    
  } else {

    // get some info about the xfer:

    double statDouble ;
    long statLong ;
    char* statString = NULL ;
    
    // known as CURLINFO_RESPONSE_CODE in later curl versions

    if(CURLE_OK == curl_easy_getinfo(ctx, CURLINFO_HTTP_CODE, &statLong)){
      COUT("Response code:  " << statLong);
    }
    
    if(CURLE_OK == curl_easy_getinfo(ctx, CURLINFO_CONTENT_TYPE, &statString)){
      COUT("Content type:   " << statString);
    }
    
    if(CURLE_OK == curl_easy_getinfo(ctx, CURLINFO_SIZE_DOWNLOAD, &statDouble)){
      COUT("Download size:  " << statDouble << " bytes");
    }
    
    if(CURLE_OK == curl_easy_getinfo(ctx, CURLINFO_SPEED_DOWNLOAD, &statDouble)){
      COUT("Download speed: " << statDouble << " bytes/sec");
    }
    
  }

  // cleanup

  curl_easy_cleanup(ctx);
  curl_global_cleanup();

  return lastRead_.str();
}

size_t CurlReader::handleData(void* buffer, size_t size, size_t nmemb, void* userp)
{
  CurlReader* reader = (CurlReader*)userp;
  reader->lastRead_ << (char*)buffer;
  return size*nmemb;
}
