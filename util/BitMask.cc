#include "gcp/util/BitMask.h"
#include "gcp/util/Exception.h"

#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
BitMask::BitMask() 
{
  byteVecPtr_ = 0;
}

/**.......................................................................
 * Destructor.
 */
BitMask::~BitMask() 
{
  if(byteVecPtr_ && internalMemory_) {
    delete byteVecPtr_;
    byteVecPtr_ = 0;
  }
}

/**.......................................................................
 * Resize this object to accomodate the specified number of bits
 */
void BitMask::resize(unsigned nBit)
{
  unsigned nByte = nBit / 8 + (nBit % 8 == 0 ? 0 : 1);

  // If the pointer already points to something, resize that object

  if(byteVecPtr_ && nByte != byteVecPtr_->size()) {
    byteVecPtr_->resize(nByte);
  } else {
    byteVecPtr_ = new std::valarray<unsigned char>(nByte);
    internalMemory_ = true;
  }
}

/**.......................................................................
 * Resize this object to accomodate the specified number of bits
 */
void BitMask::setTo(std::valarray<unsigned char>* byteVec)
{
  if(byteVecPtr_ && internalMemory_) {
    delete byteVecPtr_;
  }

  byteVecPtr_ = byteVec;
  internalMemory_ = false;
}

/**.......................................................................
 * Set the specified bit high
 */
void BitMask::setBitHigh(unsigned iBit)
{
  getBitAndByte(iBit);
  (*byteVecPtr_)[iByte_] |= (1 << iBitInByte_);
}

/**.......................................................................
 * Return true if the specified bit is high
 */
bool BitMask::bitIsHigh(unsigned iBit)
{
  return bitVal(iBit) == 1;
}

/**.......................................................................
 * Set the specified bit low
 */
void BitMask::setBitLow(unsigned iBit)
{
  getBitAndByte(iBit);
  (*byteVecPtr_)[iByte_] &= ~(1 << iBitInByte_);
}

/**.......................................................................
 * Return true if the specified bit is low
 */
bool BitMask::bitIsLow(unsigned iBit)
{
  return bitVal(iBit) == 0;
}

/**.......................................................................
 * Convert from bit number to the byte and bit number in that byte
 */
void BitMask::getBitAndByte(unsigned iBit)
{
  iByte_      = iBit / 8;
  iBitInByte_ = iBit % 8;

  if(byteVecPtr_ == 0) {
    ThrowError("Invalid bit number (" << iBit << ") for this object: no byte vector has been assigned");
  }
  
  if(iByte_ > byteVecPtr_->size()-1) {
    ThrowError("Invalid bit number (" << iBit << ") for this object: valid range is 0-" << byteVecPtr_->size()*8-1);
  }

  byteVal_    = (*byteVecPtr_)[iByte_];
}

/**.......................................................................
 * Return the specified bit value
 */
unsigned BitMask::bitVal(unsigned iBit)
{
  getBitAndByte(iBit);
  return (byteVal_ >> iBitInByte_) & 0x1;
}

/**.......................................................................
 * Set all bits high
 */
void BitMask::setAllBitsHigh()
{
  for(unsigned iByte=0; iByte < byteVecPtr_->size(); iByte++) {
    (*byteVecPtr_)[iByte] = 0xFF;
  }
}

/**.......................................................................
 * Set all bits low
 */
void BitMask::setAllBitsLow()
{
  for(unsigned iByte=0; iByte < byteVecPtr_->size(); iByte++) {
    (*byteVecPtr_)[iByte] = 0x00;
  }
}

/**.......................................................................
 * Friend operator for printing this object
 */
std::ostream& gcp::util::operator<<(std::ostream& os, BitMask& bitMask)
{
  // How many digits in the max bit number this object contains?

  unsigned nDig = (unsigned)log10(bitMask.byteVecPtr_->size() * 8)+1;

  for(unsigned iByte=0; iByte < bitMask.byteVecPtr_->size(); iByte++) {
    os << std::right << setw(nDig) << setfill(' ') << iByte*8 << setw(7) << setfill(' ') << ' ';
  }

  // Print a newline

  os << std::endl;

  for(unsigned iByte=0; iByte < bitMask.byteVecPtr_->size(); iByte++) {

    // For each byte, print blanks for the header

    os << std::right << setw(nDig-1) << setfill(' ') << ' ';

    // And print each bit of this byte

    unsigned char byte = (*bitMask.byteVecPtr_)[iByte];
    for(unsigned iBit=0; iBit < 8; iBit++) {
      os << setw(1) << ((byte >> iBit) & 0x1);
    }

  }

  return os;
}

/**.......................................................................
 * Define an operator for ORing two bitmasks together.
 */
void BitMask::operator|=(const BitMask& bm)
{
  return operator|=((BitMask&)bm);
}

/**.......................................................................
 * Define an operator for ORing two bitmasks together.
 */
void BitMask::operator|=(BitMask& bm)
{
  unsigned s1 = byteVecPtr_->size();
  unsigned s2 = bm.byteVecPtr_->size();
  unsigned nByte = s1 < s2 ? s1 : s2;

  for(unsigned iByte=0; iByte < nByte; iByte++) {
    (*byteVecPtr_)[iByte] |= (*bm.byteVecPtr_)[iByte];
  }
}

/**.......................................................................
 * Define an operator for ORing two bitmasks together.
 */
BitMask BitMask::operator|(const BitMask& bm)
{
  return operator|((BitMask&)bm);
}

BitMask BitMask::operator|(BitMask& bm)
{
  BitMask bmRet = *this;
  bmRet |= bm;
  return bmRet;
}

/**.......................................................................
 * Define an operator for ANDing two bitmasks together.
 */
void BitMask::operator&=(const BitMask& bm)
{
  return operator&=((BitMask&)bm);
}

/**.......................................................................
 * Define an operator for ANDing two bitmasks together.
 */
void BitMask::operator&=(BitMask& bm)
{
  unsigned s1 = byteVecPtr_->size();
  unsigned s2 = bm.byteVecPtr_->size();
  unsigned nByte = s1 < s2 ? s1 : s2;

  for(unsigned iByte=0; iByte < nByte; iByte++) {
    (*byteVecPtr_)[iByte] &= (*bm.byteVecPtr_)[iByte];
  }
}

/**.......................................................................
 * Define an operator for ANDing two bitmasks together.
 */
BitMask BitMask::operator&(const BitMask& bm)
{
  return operator&((BitMask&)bm);
}

BitMask BitMask::operator&(BitMask& bm)
{
  BitMask bmRet = *this;
  bmRet &= bm;
  return bmRet;
}

/**.......................................................................
 * Assignment operators
 */
void BitMask::operator=(const BitMask& bm)
{
  operator=((BitMask&) bm);
}

void BitMask::operator=(BitMask& bm)
{
  byteVecPtr_ = bm.byteVecPtr_;

}

void BitMask::operator=(std::string bitMaskStr)
{
  // If we are not currently managing anything, set up for internal
  // management of a bitmask of the length of the passed string

  if(byteVecPtr_ == 0) {
    resize(bitMaskStr.size());
  }

  // Iterate over the smaller of the string size or the number of bits
  // we are managing

  unsigned nBit    = byteVecPtr_->size() * 8;;
  unsigned nBitStr = bitMaskStr.size();
  unsigned nBitMin = nBit < nBitStr ? nBit : nBitStr;

  for(unsigned iBit=0; iBit < nBitMin; iBit++) {
    if(bitMaskStr[iBit] == '0') {
      setBitLow(iBit);
    } else {
      setBitHigh(iBit);
    }
  }
}

bool BitMask::firstNBitsAreHigh(unsigned nBit)
{
  unsigned nByte = nBit / 8;
  unsigned nBits = nBit % 8;

  if(byteVecPtr_ == 0) {
    ThrowError("No byte vector has been assigned");
  }
  
  if(nByte > byteVecPtr_->size()) {
    ThrowError("Invalid number of bits (" << nBit << ") for this object: valid range is 0-" << byteVecPtr_->size()*8-1);
  }

  unsigned char byteVal;
  for(unsigned iByte=0; iByte < nByte; iByte++) {
    byteVal = (*byteVecPtr_)[iByte];
    if(byteVal != 0xFF) {
      return false;
    }
  }

  if(nBits > 0) {
    byteVal = (*byteVecPtr_)[nByte];
    for(unsigned iBit=0; iBit < nBits; iBit++) {
      if(!((byteVal >> iBit) & 0x1)) {
	return false;
      }

    }
  }

  return true;
}
