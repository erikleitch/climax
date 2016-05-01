#include "gcp/util/Debug.h"
#include "gcp/util/Monitor.h"
#include "gcp/util/RegDescription.h"
#include "gcp/util/RegParser.h"

#include "gcp/matlab/MexHandler.h"

#include "mex.h"
#include "matrix.h"

#include <iostream.h>
#include <math.h>

using namespace std;
using namespace gcp::util;
using namespace gcp::matlab;

/**.......................................................................
 * Entry point from the matlab environment
 */
void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  gcp::util::Logger::installStdoutPrintFn(MexHandler::stdoutPrintFn);
  gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);
  
  gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
  gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
  gcp::util::ErrHandler::installLogFn(MexHandler::logFn);
  
  unsigned n3=40, n2=30, n1=20, n0=10;
  std::vector<unsigned int> cdims(4);
  std::vector<unsigned int> mdims(4);

  cdims[0] = n0;
  cdims[1] = n1;
  cdims[2] = n2;
  cdims[3] = n3;

  mdims[0] = cdims[3];
  mdims[1] = cdims[2];
  mdims[2] = cdims[1];
  mdims[3] = cdims[0];

  plhs[0] = MexHandler::createMatlabArray(4, (int*)&mdims[0], DataType::UINT);
  
  unsigned int* uintPtr = (unsigned int*)mxGetData(plhs[0]);
  
  std::vector<unsigned> cinds;
  std::vector<unsigned> minds;

  MexHandler::getIndicesMatlab(minds, mdims);

  for(unsigned index=0; index < minds.size(); index++) {
    uintPtr[index] = minds[index];
  }

  return;
}

