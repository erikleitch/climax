#include <iostream>
#include <iomanip>

#include <cmath>
#include <vector>

#include "gcp/program/Program.h"
#include "gcp/util/Angle.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Declination.h"
#include "gcp/util/Debug.h"
#include "gcp/util/OsInfo.h"
#include "gcp/util/Timer.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "nint",  "10", "i", "Number of integer elements"},
  { "ibit",  "11", "i", "Bit number of integer elements"},
  { "set",   "t",  "b", "true to set high, false to set low"},
  { "print", "t",  "b", "true to print"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

class BitMask 
{
public:

  BitMask(unsigned int* ptr, unsigned length);
  virtual ~BitMask();

  void setHigh(unsigned iBit);
  void setLow(unsigned iBit);
  unsigned bitVal(unsigned iBit);
  void iBitToIntBit(unsigned iBit, unsigned& iInt, unsigned& iIntBit);

private:

  unsigned int* ptr_;
  unsigned length_;
};

BitMask::BitMask(unsigned int* ptr, unsigned length)
{
  ptr_    = ptr;
  length_ = length;
}

BitMask::~BitMask()
{
}

/**.......................................................................
 * Retrieve this bit val
 */
unsigned BitMask::bitVal(unsigned iBit)
{
  unsigned iInt, iIntBit;
  iBitToIntBit(iBit, iInt, iIntBit);

  // Get the value of this bit

  return (ptr_[iInt]>>iIntBit) & 0x1;
}

/**.......................................................................
 * Set this bit low
 */
void BitMask::setLow(unsigned iBit)
{
  unsigned iInt, iIntBit;
  iBitToIntBit(iBit, iInt, iIntBit);

  // And set it low

  ptr_[iInt] &= ~(1<<iIntBit);
}

/**.......................................................................
 * Set this bit high
 */
void BitMask::setHigh(unsigned iBit)
{
  unsigned iInt, iIntBit;
  iBitToIntBit(iBit, iInt, iIntBit);

  // And set it high
  
  ptr_[iInt] |= (1<<iIntBit);
}
 
void BitMask::iBitToIntBit(unsigned iBit, unsigned& iInt, unsigned& iIntBit)
{
  // Find which int this bit is part of

  iInt    = iBit / 32;
  
  // And which bit of this int it corresponds to

  iIntBit = iBit % 32;
}

int Program::main()
{
#if 1
  char* cptr=0;
  COUT("Size type = " << sizeof(cptr));
  COUT("Size type = " << sizeof(size_t*));
#else
  unsigned nInt  = Program::getIntegerParameter("nint");

  //  checkRatio(nInt);

  unsigned iBit = Program::getIntegerParameter("ibit");
  
  std::vector<unsigned int> intVec(nInt);

  BitMask bitMask(&intVec[0], intVec.size());

  // As a test, set high always.  If requested to set low, it should
  // be properly undone below

  bitMask.setHigh(iBit);

  if(Program::getBooleanParameter("set")) {
    bitMask.setHigh(iBit);
  } else {
    bitMask.setLow(iBit);
  }

  if(Program::getBooleanParameter("print")) {
    COUT("Bit val is now: " << bitMask.bitVal(iBit));
  }

  for(unsigned iBit=0; iBit < intVec.size(); iBit++) {
    OsInfo::printBits(intVec[iBit]);
  }
#endif

  return 0;
}

void checkRatio(unsigned nInt)
{
  unsigned nByte = 4 * nInt;

  std::vector<unsigned char> charArr(nByte);
  std::vector<unsigned int> intArr(nInt);

  std::vector<unsigned char> charArr2(nByte);
  std::vector<unsigned int> intArr2(nInt);

  COUT("Here 0");
  Timer timer;

  timer.start();

  unsigned int* intPtr1 = &intArr[0];
  unsigned int* intPtr2 = &intArr2[0];
  for(unsigned i=0; i < nInt; i++) {
    *(intPtr1+i) |= *(intPtr2+i);
  }

  timer.stop();

  COUT("Here 1");

  double delta1 = timer.deltaInSeconds();
 
  timer.start();

  unsigned char* charPtr1 = &charArr[0];
  unsigned char* charPtr2 = &charArr2[0];
  for(unsigned i=0; i < nByte; i++) {
    *(charPtr1+i) |= *(charPtr2+i);
  }

  timer.stop();

  COUT("Here 2");  

  double delta2 = timer.deltaInSeconds();

  COUT("Here 3");
  COUT("Ratio of char to int elapsed time is: " << delta2/delta1);
}
