// $Id: emacs_macros,v 1.1.1.1.4.1 2006/01/25 08:01:42 spt Exp $

#ifndef GCP_UTIL_CIRCBUF_H
#define GCP_UTIL_CIRCBUF_H

/**
 * @file CircBuf.h
 * 
 * Tagged: Tue Jan 23 17:04:06 NZDT 2007
 * 
 * @version: $Revision: 1.1.1.1.4.1 $, $Date: 2006/01/25 08:01:42 $
 * 
 * @author Erik Leitch
 */

#include <iostream>
#include <valarray>

#include "gcp/util/Exception.h"

namespace gcp {
  namespace util {

    // Forward declarations needed before we can declare the output operators below

    template <class type>
      class CircBuf;

    template <class type>
      class Node;

    // Output Operators for classes defined in this file

    template<class type>
      std::ostream& operator<<(std::ostream& os, const CircBuf<type>& obj);

    template<class type>
      std::ostream& operator<<(std::ostream& os, const Node<type>& obj);

    //=======================================================================
    // A utility class to encapsulate a single node in the circular buffer
    //=======================================================================

    template<class type>
      class Node {
      public:
      
      unsigned id_; 
      type node_;
      
      Node* next_;
      Node* prev_;
      
      Node() {
	id_    = 0;
	node_  = 0;
	next_  = 0;
	prev_  = 0;
      }
      
      virtual ~Node() {};
    };

    //=======================================================================
    // Circular buffer class methods.  
    //
    // This class is intended to be used as storage for a sequential
    // stream of objects where at any given point in time we are
    // interested in an n-element stretch of them.  
    //
    // An example is storage for on-the-fly computation of an n-point
    // FFT of timestream data.  As new samples are clocked in, the
    // oldest ones are simply overwritten, so that a call to copy
    // returns an array n-elements long containing the latest data.
    // 
    //=======================================================================

    template <class type>
      class CircBuf {
      public:

      // Constructor.

      CircBuf(unsigned n=1);

      // Copy Constructor.

      CircBuf(const CircBuf<type>& objToBeCopied);

      // Copy Constructor.

      CircBuf(CircBuf<type>& objToBeCopied);

      // Const Assignment Operator.

      void operator=(const CircBuf<type>& objToBeAssigned);

      // Assignment Operator.

      void operator=(CircBuf<type>& objToBeAssigned);

      // Destructor.

      virtual ~CircBuf();

      // Push a new object onto the head of the queue.  If the buffer
      // is full, this will displace the oldest one.
	
      void push(type);
	
      // Return the newest/oldest element in the buffer

      type newest();
      type oldest();

      // Return the length of the buffer

      unsigned size();

      // Return the number of objects currently in the buffer (could
      // be < size())

      unsigned nInBuffer();

      // resize the buffer

      void resize(unsigned n);

      // (re-) initialize the buffer

      void initialize(unsigned n);

      // Unsafe method to copy the current queue to an array

      void copy(type* ptr);

      // Safe method to copy the current queue to an array

      std::valarray<type> copy();

      private:

      unsigned nTotal_;
      unsigned nInBuf_;

      Node<type>* head_;
      Node<type>* tail_;

      std::valarray<Node<type> > nodes_;

    }; // End template class CircBuf

    //-----------------------------------------------------------------------
    // Implementation
    //-----------------------------------------------------------------------

    // Constructors

    template <class type>
      CircBuf<type>::CircBuf(unsigned n)
      {
	initialize(n);
      }

    // Copy Constructor.

    template <class type>
      CircBuf<type>::CircBuf(const CircBuf<type>& objToBeCopied)
      {
	*this = (CircBuf<type>&)objToBeCopied;
      }

    // Copy Constructor.

    template <class type>
      CircBuf<type>::CircBuf(CircBuf<type>& objToBeCopied)
      {
	*this = objToBeCopied;
      };

    
    // Const Assignment Operator.

    template <class type>
      void CircBuf<type>::operator=(const CircBuf<type>& objToBeAssigned)
      {
	*this = (CircBuf<type>&)objToBeAssigned;
      };
    
    // Assignment Operator.

    template <class type>    
      void CircBuf<type>::operator=(CircBuf<type>& objToBeAssigned)
      {
	std::cout << "Calling default assignment operator for class: CircBuf" << std::endl;
      };
    
    // Output Operator.

    template <class type>
      std::ostream& operator<<(std::ostream& os, CircBuf<type>& obj)
      {
	COUT("Default output operator called");
      }

    // Destructor.

    template <class type>
      CircBuf<type>::~CircBuf() {};

    // Push a new object onto the head of the queue.  If the queue is
    // full, this will displace the oldest one.
    
    template <class type>
      void CircBuf<type>::push(type obj) 
      {
	tail_->node_ = obj;
	tail_ = tail_->next_;

	if(tail_ == head_) {
	  head_ = head_->next_;
	} else {
	  nInBuf_++;
	}
      }
    
    // Return the oldest (head) sample

    template <class type>
      type CircBuf<type>::oldest() 
      {
	return head_->node_;
      }

    template <class type>
      type CircBuf<type>::newest() 
      {
	return tail_->prev_->node_;
      }

    // Return the length of the buffer
    
    template<class type>
      unsigned CircBuf<type>::size() 
      {
	return nTotal_;
      }
    
    // Return the number of objects currently in the buffer (could
    // be < size())
    
    template<class type>
      unsigned CircBuf<type>::nInBuffer() 
      {
	return nInBuf_;
      }
    
    // resize the buffer
    
    template<class type>
      void CircBuf<type>::resize(unsigned n)
      {
	initialize(n);
      }

    template<class type>
      void CircBuf<type>::initialize(unsigned n)
      {
	nInBuf_ = 0;

	// Make the queue larger than the number of nodes needed by
	// one.  This is so that tests like for(iter = head; iter !=
	// tail; iter++) will do the correct thing (otherwise, this
	// loop would stop one short)

	if(nodes_.size() != n+1)
	  nodes_.resize(n+1);

	nTotal_ = n;
	
	// And re-link the list
	
	for(unsigned iNode=0; iNode <= n; iNode++) {
	  nodes_[iNode].id_   = iNode;
	  nodes_[iNode].node_ = 0;
	  nodes_[iNode].prev_ = &nodes_[iNode==0 ? n : iNode-1];
	  nodes_[iNode].next_ = &nodes_[iNode==n ? 0 : iNode+1];
	}
	
	// And set the queue ptr pointing to the head of the list
	
	head_ = &nodes_[0];
	
	// And set the queue ptr pointing to the tail of the list.  This
	// should always point to the next unused slot in the array.
	// Since our arry is 1 element longer than we need for
	// computation, this will start pointing to the 0th element, and
	// eventually will always point to element nMaxNodes-1
	
	tail_ = &nodes_[0];
      }

    // Unsafe method to copy the current queue to an array

    template<class type>
      void CircBuf<type>::copy(type* ptr)
      {
	unsigned i=0;

	for(Node<type>* iter = head_; iter != tail_; iter = iter->next_, i++) 
	  *(ptr+i) = iter->node_;
      }

    // Safe method to copy the current queue to an array

    template<class type>
      std::valarray<type> CircBuf<type>::copy()
      {
	std::valarray<type> arr(nInBuf_);
	copy(&arr[0]);
	return arr;
      }
    
    // Print out a buffer

    template <class type>
      std::ostream& operator<<(std::ostream& os, 
			       const CircBuf<type>& buf)
      {
	COUT("Calling default output operator for CircBuf");
      }

    // Print out a buffer

    template <class type>
      std::ostream& operator<<(std::ostream& os, 
			       const Node<type>& node)
      {
	COUT(node.node_ << " " << node.id_);
      }


  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_CIRCBUF_H
