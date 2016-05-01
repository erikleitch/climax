#ifndef GCP_UTIL_RANGE_H
#define GCP_UTIL_RANGE_H

/**
 * @file Range.h
 * 
 * Tagged: Fri Sep 17 15:48:20 PDT 2004
 * 
 * @author Erik Leitch
 */
#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/String.h"

#include <string>
#include <vector>

namespace gcp {
  namespace util {
    
    template<class Type>
      class Range;

    template <class Type>
      std::ostream& operator<<(std::ostream& os, 
			       const Range<Type>& range);

    // Class for managing a range

    template<class Type>
      class Range {
    public:
      
      /**
       * Constructor.
       */
      Range();

      /**
       * Constructor.
       */
      Range(Type start, Type stop);
      
      /**
       * Destructor.
       */
      virtual ~Range();
      
      void setStart(Type start);
      void setStop(Type stop);

      Type start();
      Type stop();

      /** 
       * Declare a Range printing method
       */
      friend std::ostream& operator << <>
	(std::ostream& os, const Range<Type>& range);

      /**
       * Declare an operator for incrementing this object
       */
      Range& operator+=(Type incr);

      /**
       * Declare an operator for multiplying this object
       */
      Range& operator*=(Type mult);

      std::vector<unsigned> extractIndexRange(String& indStr, unsigned lowestValid, unsigned highestValid, 
					      unsigned baseIndex, bool actualIndex);


      static void addIndex(std::vector<unsigned>& indices, unsigned index, 
		    unsigned lowestValid, unsigned highestValid);

      static unsigned parseIndexExpression(String& str, 
					   unsigned baseIndex,   unsigned actualIndex, 
					   unsigned lowestValid, unsigned highestValid);

      static void parseIndexOperands(String& str, unsigned& op1, unsigned& op2, std::string op,
				     unsigned baseIndex,   unsigned actualIndex, 
				     unsigned lowestValid, unsigned highestValid);

      static unsigned firstOddIndex(unsigned lowestValid, unsigned highestValid);
      static unsigned firstEvenIndex(unsigned lowestValid, unsigned highestValid);

    private:
      
      Type start_;
      bool startInit_;
      Type stop_;
      bool stopInit_;
      
    }; // End class Range

    /**
     * Constructor
     */
    template<class Type>
      Range<Type>::Range() {
      startInit_ = false;
      stopInit_   = false;
    };

    // Constructor with initialization

    template<class Type>
      Range<Type>::Range(Type start, Type stop) {
      startInit_ = false;
      stopInit_   = false;

      setStart(start);
      setStop(stop);
    };

    /**
     * Destructor
     */
    template<class Type>
      Range<Type>::~Range() {};

    /**
     * Set the start of the range.
     */
    template<class Type>
      void Range<Type>::setStart(Type start) {
      start_ = start;
      startInit_ = true;
    }

    /**
     * Set the end of the range.
     */
    template<class Type>
      void Range<Type>::setStop(Type stop) {
      stop_ = stop;
      stopInit_ = true;
    }

    /**
     * Get the start of the range.
     */
    template<class Type>
      Type Range<Type>::start() {
      if(!startInit_) {
	ThrowError("Start value has not been initialized");
      }
      return start_;
    }

    /**
     * Get the end of the range.
     */
    template<class Type>
      Type Range<Type>::stop() {
      if(!startInit_) {
	ThrowError("End value has not been initialized");
      }
      return stop_;
    }

    /**
     * Print out a range
     */
    template <class Type>
      std::ostream& operator<<(std::ostream& os, 
			       const Range<Type>& range)
    {
      os << "(" << range.start_ << "-" << range.stop_ << ")";
      return os;
    }
    

    /**
     * Increment this object
     */
    template <class Type>
      Range<Type>& Range<Type>::operator+=(Type incr) 
    {
      start_ += incr;
      stop_  += incr;

      return *this;
    }

    /**
     * Multiply this object
     */
    template <class Type>
      Range<Type>& Range<Type>::operator*=(Type mult) 
    {
      start_ *= mult;
      stop_  *= mult;
      return *this;
    }

    /**.......................................................................
     * Parse an index range of the form [1,2,3] or [3-5] or 6, into a
     * vector of requested indices.
     */
    template <class Type>
      std::vector<unsigned> Range<Type>::
      extractIndexRange(String& indStr, unsigned lowestValid, unsigned highestValid, unsigned baseIndex, 
			bool actualIndex)
    {
      std::vector<unsigned> indices;
  
      //------------------------------------------------------------
      // If the string contains a "[", then this is a list or range of
      // indices
      //------------------------------------------------------------
  
      if(indStr.contains("[")) {

	// Is this a range? (ie, [1-15] )

	if(indStr.contains("-")) {
      
	  // Was an increment specified?
      
	  unsigned iStart, iStop, incr=1;

	  String startStr;
	  String stopStr;
      
	  if(indStr.contains(";")) {
	    startStr = indStr.findNextInstanceOf("[", true, "-", true, false);
	    stopStr  = indStr.findNextInstanceOf("-", true, ";", true, false);
	
	    String incrStr  = indStr.findNextInstanceOf(";", true, "]", true, false);
	    incr = incrStr.toInt();
	  } else {
	    startStr = indStr.findNextInstanceOf("[", true, "-", true, false);
	    stopStr  = indStr.findNextInstanceOf("-", true, "]", true, false);
	  }
      
	  iStart = parseIndexExpression(startStr, baseIndex, actualIndex, lowestValid, highestValid);
	  iStop  = parseIndexExpression(stopStr,  baseIndex, actualIndex, lowestValid, highestValid);

	  if(iStart > iStop) {
	    ThrowColorError("Invalid range: " << indStr 
			    << " (Your start index must be less than your stop index)", "red");
	  }
      
	  for(unsigned iInd=iStart; iInd <= iStop; iInd += incr) {
	    addIndex(indices, iInd, lowestValid, highestValid);
	  }

	  // Else a list of indices? (ie, [1,2,3] )
      
	} else if(indStr.contains(",")) {

	  COUT("indStr = " << indStr);

	  String ind;

	  ind = indStr.findNextInstanceOf("[", true, ",", true, false);

	  COUT("ind is now = " << ind);

	  if(!ind.isEmpty()) {
	    addIndex(indices, parseIndexExpression(ind, baseIndex, actualIndex, lowestValid, highestValid), 
		     lowestValid, highestValid);
	  } else {
	    ThrowColorError("Invalid list: '" << indStr 
			    << "' (A list of indices should be of the form: [1,2,3])", "red");
	  }

	  do {
	    ind = indStr.findNextInstanceOf(",", true, ",", true, false);

	    COUT("ind is now = " << ind << " empty = " << ind.isEmpty());

	    if(!ind.isEmpty()) {
	      addIndex(indices, parseIndexExpression(ind, baseIndex, actualIndex, lowestValid, highestValid), 
		       lowestValid, highestValid);
	    }

	  } while(!ind.isEmpty());

	  ind = indStr.findNextInstanceOf(",", true, "]", true, false);

	  COUT("ind is now = " << ind);

	  if(!ind.isEmpty()) {
	    addIndex(indices, parseIndexExpression(ind, baseIndex, actualIndex, lowestValid, highestValid), 
		     lowestValid, highestValid);
	  } else {
	    ThrowColorError("Invalid list(2): " << indStr 
			    << " (A list of indices should be of the form: [1,2,3])", "red");
	  }

	  // Else a single index was specified

	} else {
	  String ind = indStr.findNextInstanceOf("[", true, "]", false);
	  addIndex(indices, parseIndexExpression(ind, baseIndex, actualIndex, lowestValid, highestValid), 
		   lowestValid, highestValid);
	}

	//------------------------------------------------------------
	// Else this was not a list.  Check for allowed regexp (*)
	//------------------------------------------------------------

      } else {

	if(indStr[0] == '*') {

	  for(unsigned iInd=lowestValid; iInd <= highestValid; iInd++) {
	    addIndex(indices, iInd, lowestValid, highestValid);
	  }

	  // Else this was a single number

	} else {
	  addIndex(indices, indStr.toInt(), lowestValid, highestValid);
	}
      }

      return indices;
    }


    /**.......................................................................
     * Add an index to a vector of requested indices
     */
    template<class Type>
      void Range<Type>::addIndex(std::vector<unsigned>& indices, unsigned index, 
				 unsigned lowestValid, unsigned highestValid)
    {
      if(index < lowestValid || index > highestValid) {
	ThrowColorError("Invalid index: " 
			<< index << " (should be " << lowestValid << "-" << highestValid << ")", "red");
      }

      indices.push_back(index);
    }
    
    /**.......................................................................
     * Parse an index expression of the form N+M*i
     */
    template<class Type>
      unsigned Range<Type>::parseIndexExpression(String& str, 
						 unsigned baseIndex,   unsigned actualIndex, 
						 unsigned lowestValid, unsigned highestValid)
    {
      unsigned op1, op2;
      unsigned index = baseIndex;

      if(str.contains("+")) {
	parseIndexOperands(str, op1, op2, "+", baseIndex, actualIndex, lowestValid, highestValid);
	index = op1+op2;
      } else if(str.contains("*")) {
	parseIndexOperands(str, op1, op2, "*", baseIndex, actualIndex, lowestValid, highestValid);
	index = op1*op2;
      } else if(str.contains("ANY")) {

	if(baseIndex == 0 && actualIndex) {
	  ThrowColorError("No base index has been specified for an implicit rule: " << str, "red");
	}

	if(!actualIndex) {
	  index = lowestValid;
	}

      } else if(str.contains("EVEN")) {

	if(actualIndex) {
	  if(baseIndex % 2 != 0) {
	    ThrowColorError("An implicit rule requires an even index (" << str << "), but you have specified an odd index: " << baseIndex, "red");
	  } 
	} else {
	  index = firstEvenIndex(lowestValid, highestValid);
	}

      } else if(str.contains("ODD")) {

	if(actualIndex) {
	  if(baseIndex % 2 == 0) {
	    ThrowColorError("An implicit rule requires an odd index (" << str << "), but you have specified an even index: " << baseIndex, "red");
	  }
	} else {
	  index = firstOddIndex(lowestValid, highestValid);
	}
      } else {
	index = str.toInt();
      }

      return index;
    };

    /**.......................................................................
     * Return the value of two operands in an expression like: 'i op j'
     */
    template<class Type>
      void Range<Type>::parseIndexOperands(String& str, unsigned& op1, unsigned& op2, std::string op,
					   unsigned baseIndex,   unsigned actualIndex, 
					   unsigned lowestValid, unsigned highestValid)

    {
      String op1Str = str.findNextInstanceOf(" ", false, op, true, false);
      String op2Str = str.findNextInstanceOf(op,  true,  op, false, false);

      op1 = parseIndexExpression(op1Str, baseIndex, actualIndex, lowestValid, highestValid);
      op2 = parseIndexExpression(op2Str, baseIndex, actualIndex, lowestValid, highestValid);
    };

    template<class Type>
      unsigned Range<Type>::firstEvenIndex(unsigned lowestValid, unsigned highestValid)
    {
      if(lowestValid % 2 == 0)
	return lowestValid;
      
      if(lowestValid+1 > highestValid) {
	
	ThrowColorError("No valid even index can be constructed from the range: [" << lowestValid << "-" << highestValid << "]", "red");
	
	// We will never get here -- just to avoid compiler warnings.
	
	return 0;
	
      } else {
	return lowestValid+1;
      }
    };
    
    template<class Type>
      unsigned Range<Type>::firstOddIndex(unsigned lowestValid, unsigned highestValid)
    {
      if(lowestValid % 2 != 0)
	return lowestValid;
      
      if(lowestValid+1 > highestValid) {
	
	ThrowColorError("No valid odd index can be constructed from the range: [" << lowestValid << "-" << highestValid << "]", "red");
	
	// We will never get here -- just to avoid compiler warnings.
	
	return 0;
	
      } else {
	return lowestValid+1;
      }
    };

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_RANGE_H
