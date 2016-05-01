#ifndef GCP_UTIL_BITMASK_H
#define GCP_UTIL_BITMASK_H

/**
 * @file BitMask.h
 * 
 * Tagged: Tue Mar  6 16:39:02 PST 2012
 * 
 * @version: $Revision: 1.1 $, $Date: 2012/03/15 01:19:32 $
 * 
 * @author username: Erik Leitch
 */
#include <valarray>
#include <sstream>

namespace gcp {
  namespace util {

    class BitMask {
    public:

      // Constructor.

      BitMask();

      // Destructor.

      virtual ~BitMask();

      // Resize this bit mask to be at least as large as the requested
      // number of bits

      void resize(unsigned nBit);

      // Set up this object to manage an externally provided byte
      // array

      void setTo(std::valarray<unsigned char>* byteVec);

      // Set the requested bit high

      void setBitHigh(unsigned iBit);
      void setAllBitsHigh();

      // Set the requested bit low

      void setBitLow(unsigned iBit);
      void setAllBitsLow();

      // Return if the requested bit is high or low

      bool bitIsHigh(unsigned iBit);
      bool bitIsLow(unsigned iBit);

      bool firstNBitsAreHigh(unsigned nBit);

      // Return if all bits are high

      bool allBitAreHigh();

      // Return the requested bit value

      unsigned bitVal(unsigned iBit);

      // Print operator

      friend std::ostream& operator<<(std::ostream& os, BitMask& bitMask);

      // Assignment operators: bm1 = bm2

      void operator=(const BitMask& bitMask);
      void operator=(BitMask& bitMask);

      // Operator for bm1 |= bm2

      void operator|=(BitMask& bitMask);
      void operator|=(const BitMask& bitMask);
      BitMask operator|(BitMask& bitMask);
      BitMask operator|(const BitMask& bitMask);

      // Operator for bm1 &= bm2

      void operator&=(BitMask& bitMask);
      void operator&=(const BitMask& bitMask);
      BitMask operator&(BitMask& bitMask);
      BitMask operator&(const BitMask& bitMask);

      void operator=(std::string bitMaskStr);

    public:

      void getBitAndByte(unsigned iBit);

      unsigned iByte_;
      unsigned iBitInByte_;
      unsigned byteVal_;
      bool internalMemory_;

      std::valarray<unsigned char>* byteVecPtr_;

    }; // End class BitMask

    std::ostream& operator<<(std::ostream& os, BitMask& bitMask);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_BITMASK_H
