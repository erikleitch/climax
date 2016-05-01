/**.......................................................................
 * MATLAB Mex file for writing out UVF files
 *
 * Use like:
 *
 * d=gcpMatReadArc({'array.frame.record','corr.band0.usb[0][0]',
 *                  'antenna*.tracker.actual double', 'antenna*.tracker.source string'},
 *                  '06-jan-2005:15','06-jan-2005:16',
 *                  '/data/gcpdaq/arc','/home/gcpdaq/carma_unstable/gcp/array/conf/cal');
 */
#include <iostream>
#include <math.h>

#include "gcp/util/Debug.h"
#include "gcp/util/HourAngle.h"
#include "gcp/util/Frequency.h"

#include "gcp/matlab/MexHandler.h"
#include "gcp/matlab/UvfWriter.h"

#include "mex.h"
#include "matrix.h"

using namespace gcp::util;
using namespace gcp::matlab;

void writeUvfFile(int nrhs, const mxArray* prhs[]);
void writeFakeUvfFile(const mxArray* prhs[]);

void mexFunction(int nlhs, mxArray      *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  gcp::util::Logger::installStdoutPrintFn(MexHandler::stdoutPrintFn);
  gcp::util::Logger::installStderrPrintFn(MexHandler::stderrPrintFn);
  
  gcp::util::ErrHandler::installThrowFn(MexHandler::throwFn);
  gcp::util::ErrHandler::installReportFn(MexHandler::reportFn);
  gcp::util::ErrHandler::installLogFn(MexHandler::logFn);

  if(nrhs < 8) {
    std::cerr << "Wrong number of arguments. " << std::endl
	      << "Should be: fileName data[NsxNbxNfx3]"
	      << " JD[Nsx2] srcname coord[3] freq[Nf]"
	      << " XYZ[3xnTelescope] uvw[NsxNbx3]" << std::endl;
    return;
  }

  writeUvfFile(nrhs, prhs);
}

void writeUvfFile(int nrhs, const mxArray* prhs[])
{
  UvfWriter uvf;

  // optional arguments

  if (nrhs>8) {
    uvf.setNumberOfTelescopes(prhs[8]);
  }

  if (nrhs>10) {
    uvf.setBaselines(prhs[10]);
  }
  if (nrhs>11) {
    uvf.setCoordsToJ2000(prhs[11]);
  }
  if (nrhs>12) {
    uvf.setFirstTelescopeNum(prhs[12]);
  }

  // Set the filename

  uvf.setFileName(prhs[0]);
  uvf.setUvfData(prhs[1]);
  uvf.setDate(prhs[2]);

  if (nrhs>9) {
    uvf.setDeltaIfFrequencies(prhs[9]);
  }

  uvf.setSourceName(prhs[3]);

  uvf.setCoord(prhs[4]);
  uvf.setFreq(prhs[5]);

  uvf.setXyz(prhs[6]);
  uvf.setUvw(prhs[7]);

  uvf.writeUvfFile();
}

void writeFakeUvfFile(const mxArray* prhs[])
{
  UvfWriter uvf;

  // Set the filename

  uvf.setFileName(prhs[0]);
  uvf.setUvfData(prhs[1]);
  uvf.setDate(prhs[2]);

  uvf.setSourceName(prhs[3]);

  uvf.setCoord(prhs[4]);
  uvf.setFreq(prhs[5]);

  uvf.setXyz(prhs[6]);
  uvf.setUvw(prhs[7]);

  uvf.writeFakeUvfFile();
}

