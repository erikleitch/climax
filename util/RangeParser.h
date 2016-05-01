// $Id: RangeParser.h,v 1.1 2013/11/19 23:54:10 eml Exp $

#ifndef GCP_UTIL_RANGEPARSER_H
#define GCP_UTIL_RANGEPARSER_H

/**
 * @file RangeParser.h
 * 
 * Tagged: Thu Nov  7 13:13:38 PST 2013
 * 
 * @version: $Revision: 1.1 $, $Date: 2013/11/19 23:54:10 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/String.h"

#include <vector>

namespace gcp {
  namespace util {

    class RangeParser {
    public:

      std::vector<unsigned> extractIndexRange(gcp::util::String& antStr, unsigned lowestValid=0, unsigned highestValid=0,
					      unsigned baseIndex=0, bool actualIndex=true);

      unsigned parseIndexExpression(gcp::util::String& str, 
				    unsigned baseIndex,   unsigned actualIndex, 
				    unsigned lowestValid, unsigned highestValid);

      void parseIndexOperands(gcp::util::String& str, unsigned& op1, unsigned& op2, std::string op,
			      unsigned baseIndex,   unsigned actualIndex, 
			      unsigned lowestValid, unsigned highestValid);

      void addIndex(std::vector<unsigned>& indices, unsigned index, unsigned lowestValid, unsigned highestValid);

      unsigned firstEvenIndex(unsigned lowestValid, unsigned highestValid);
      unsigned firstOddIndex(unsigned lowestValid, unsigned highestValid);

    }; // End class RangeParser

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_RANGEPARSER_H
