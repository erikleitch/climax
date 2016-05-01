// $Id: emacs_macros,v 1.1.1.1.4.1 2006/01/25 08:01:42 spt Exp $

#ifndef GCP_UTIL_SORT_H
#define GCP_UTIL_SORT_H

/**
 * @file Sort.h
 * 
 * Tagged: Wed Feb 14 09:35:58 NZDT 2007
 * 
 * @version: $Revision: 1.1.1.1.4.1 $, $Date: 2006/01/25 08:01:42 $
 * 
 * @author Erik Leitch
 */

#include <iostream>
#include <vector>
#include <string>

namespace gcp {
  namespace util {

    class Sort {
    public:

      /**
       * Constructor.
       */
      Sort();

      /**
       * Destructor.
       */
      virtual ~Sort();

      static std::vector<std::string> sort(std::vector<std::string>& entries);

    }; // End class Sort

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_SORT_H
