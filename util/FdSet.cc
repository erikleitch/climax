#define __FILEPATH__ "util/common/FdSet.cc"

#include <sys/select.h>

#include "gcp/util/Debug.h"
#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/FdSet.h"

using namespace std;
using namespace gcp::util;

#define CHECK_FD(fd) \
  {\
    LogStream errStr;\
    if(fd < 0) {\
      errStr.appendMessage(true, "File descriptor is < 0");\
/*      errStr.report();*/\
      return;\
    }\
  }

/**.......................................................................
 * Constructor.
 */
FdSet::FdSet() 
{
  clear();
};

/**.......................................................................
 * Destructor.
 */
FdSet::~FdSet() {};

/**.......................................................................
 * Private method to zero the set of descriptors to be watched for
 * readability.
 */
void FdSet::zeroReadFdSet()
{
  // Zero the set of descriptors to be watched for input
  
  FD_ZERO(&readFdSetSave_);
}

/**.......................................................................
 * Private method to zero the set of descriptors to be watched for
 * readability.
 */
void FdSet::zeroWriteFdSet()
{
  // Zero the set of descriptors to be watched for writability
  
  FD_ZERO(&writeFdSetSave_);
}

/**.......................................................................
 * Private method to zero the set of descriptors to be watched for
 * exceptions
 */
void FdSet::zeroExceptionFdSet()
{
  // Zero the set of descriptors to be watched for exceptions
  
  FD_ZERO(&exceptionFdSetSave_);
}

/**.......................................................................
 * Private method to register a file descriptor to be watched for
 * input.
 */
void FdSet::registerReadFd(int fd)
{
  // Register this descriptor to be watched for output

  CHECK_FD(fd);

  FD_SET(fd, &readFdSetSave_);
  
  if(fd + 1 > fdSetSize_)
    fdSetSize_ = fd + 1;
}

/**.......................................................................
 * Private method to register a file descriptor to be watched for
 * output.
 */
void FdSet::registerWriteFd(int fd)
{
  // And register this descriptor to be watched for input

  CHECK_FD(fd);
  
  FD_SET(fd, &writeFdSetSave_);
  
  if(fd + 1 > fdSetSize_)
    fdSetSize_ = fd + 1;
}

/**.......................................................................
 * Private method to register a file descriptor to be watched for
 * output.
 */
void FdSet::registerExceptionFd(int fd)
{
  // And register this descriptor to be watched for input

  CHECK_FD(fd);
  
  FD_SET(fd, &exceptionFdSetSave_);
  
  if(fd + 1 > fdSetSize_)
    fdSetSize_ = fd + 1;
}

/**.......................................................................
 * Remove a file descriptor from the mask of descriptors to be watched
 * for readability.
 */
void FdSet::clearFromReadFdSet(int fd) 
{
  CHECK_FD(fd);

  FD_CLR(fd, &readFdSetSave_); 
}

/**.......................................................................
 * Remove a file descriptor from the mask of descriptors to be watched
 * for writeability.
 */
void FdSet::clearFromWriteFdSet(int fd)
{
  CHECK_FD(fd);

  FD_CLR(fd, &writeFdSetSave_);
}

/**.......................................................................
 * Remove a file descriptor from the mask of descriptors to be watched
 * for exceptionability.
 */
void FdSet::clearFromExceptionFdSet(int fd)
{
  CHECK_FD(fd);

  FD_CLR(fd, &exceptionFdSetSave_);
}

/**.......................................................................
 * Return a pointer to the set of read file descriptors
 */
fd_set* FdSet::readFdSet()
{
  // Copy the contents of our saved fd_set

  readFdSet_ = readFdSetSave_;

  // And return the copy

  return &readFdSet_;
}

/**.......................................................................
 * Return a pointer to the set of write file descriptors
 */
fd_set* FdSet::writeFdSet()
{
  // Copy the contents of our saved fd_set

  writeFdSet_ = writeFdSetSave_;

  // And return the copy

  return &writeFdSet_;
}

/**.......................................................................
 * Return a pointer to the set of exception file descriptors
 */
fd_set* FdSet::exceptionFdSet()
{
  // Copy the contents of our saved fd_set

  exceptionFdSet_ = exceptionFdSetSave_;

  // And return the copy

  return &exceptionFdSet_;
}

/**.......................................................................
 * Return the size of the largest file descriptor in the set.
 */
int FdSet::size() 
{
  return fdSetSize_;
}

/**.......................................................................
 * Return true if the file descriptor is set in the read set.
 */
bool FdSet::isSetInRead(int fd)
{
  return fd < 0 ? false : FD_ISSET(fd, &readFdSet_);
}

/**.......................................................................
 * Return true if the file descriptor is set in the write set.
 */
bool FdSet::isSetInWrite(int fd) 
{
  return fd < 0 ? false : FD_ISSET(fd, &writeFdSet_);
}

/**.......................................................................
 * Return true if the file descriptor is set in the exception set.
 */
bool FdSet::isSetInException(int fd) 
{
  return fd < 0 ? false : FD_ISSET(fd, &exceptionFdSet_);
}

void FdSet::print()
{
#if MAC_OSX == 0
    cout << "readFdSetSave_ is: " << endl;
    for(unsigned ibit=0; ibit < __FD_SETSIZE/__NFDBITS; ibit++)
      cout << __FDS_BITS(&readFdSetSave_)[ibit] << " ";
    cout << endl;

    cout << "readFdSet_ is: " << endl;
    for(unsigned ibit=0; ibit < __FD_SETSIZE/__NFDBITS; ibit++)
      cout << __FDS_BITS(&readFdSet_)[ibit] << " ";
    cout << endl;
#endif
}

/**.......................................................................
 * Zero our read and write fd sets.
 */
void FdSet::clear() 
{
  fdSetSize_ = 0;
  zeroReadFdSet();
  zeroWriteFdSet();
  zeroExceptionFdSet();
}

/**.......................................................................
 * Clear an fd from both read and write fd sets.
 */
void FdSet::clear(int fd) 
{
  clearFromReadFdSet(fd);
  clearFromWriteFdSet(fd);
  clearFromExceptionFdSet(fd);
}
