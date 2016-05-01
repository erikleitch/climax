#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/util/Debug.h"
#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/IoLock.h"
#include "gcp/util/Exception.h"
#include "gcp/util/FitsBinTableReader.h"
#include "gcp/util/FitsUvfReader.h"
#include "gcp/util/OsInfo.h"
#include "gcp/util/Timer.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/fftutil/Dft2d.h"
#include "gcp/fftutil/Image.h"

#include "gcp/datasets/VisDataSetUvf.h"

#include "gcp/models/GaussianClusterModel.h"
#include "gcp/models/PtSrcModel.h"

using namespace std;
using namespace gcp::datasets;
using namespace gcp::models;
using namespace gcp::program;
using namespace gcp::util;

KeyTabEntry Program::keywords[] = {
  { "dev",      "/xs",            "s", "Pgplot device"},
  { "file",     "/Users/eml/projects/climax/climaxTestSuite/A1914.uvf",               "s", "FITS file to read in"},
  { "file1",     "gen4.uvf",               "s", "FITS file to read in"},
  { "file2",     "gen6.uvf",               "s", "FITS file to read in"},
  { "zeropad",  "f",              "b", "Zeropad the array?"},
  { "phase",    "0.0",            "d", "Phase (degrees)"},
  { "period",   "64",             "d", "Period (pixels)"},
  { "freq",     "30",             "d", "freq (GHz)"},
  { "perc",     "0.98",           "d", "percent correlation"},
  { "nthread",  "1",              "i", "number of threads to use"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

void testCoadd(std::string& file);
void testReader(std::string& file);
void testChisqCalcPtSrc(std::string& file, double perc);
void testSimulator(std::string& file, double perc);
void testReplacePtSrc(std::string& file, double perc);
void testReplaceGauss(std::string& file, double perc);
void testReplaceGaussMultiThread(std::string& file, double perc, unsigned nthread);
void computeChisqPtSrc(Angle& xOff, Angle& yOff, PtSrcModel& ptsrc, VisDataSet& vds);
void testChisqGauss(VisDataSet& vds);
void testUvCat(std::string& file1, std::string& file2);
void testInit(std::string file);

int Program::main()
{
  Debug::setLevel(Debug::DEBUGNONE);

  std::string file = Program::getStringParameter("file");
  double perc      = Program::getDoubleParameter("perc");
  unsigned nthread = Program::getIntegerParameter("nthread");
  
  //testChisqCalc(file, perc);

  std::string file1 = Program::getStringParameter("file1");
  std::string file2 = Program::getStringParameter("file2");
  //  testUvCat(file1, file2);
  testInit(file);
  //testCoadd(file);
  //  testChisqCalcPtSrc(file, perc);
  //  testSimulator(file, perc);
  //testReplacePtSrc(file, perc);
  //  testReplaceGauss(file, perc);
  //testReplaceGaussMultiThread(file, perc, nthread);
}

/**.......................................................................
 * Test the chi-squared calculation
 */
void testUvCat(std::string& file1, std::string& file2)
{
  VisDataSetUvf vds1(0);
  VisDataSetUvf vds2(0);

  vds1.setName("vds1");
  vds2.setName("vds2");

  std::ostringstream os;

  os.str("");
  os << file1 << ", shift=0,0 deg";
  vds1.setParameter("file", os.str());
  vds1.setParameter("uvmax", "2000");
  vds1.setParameter("interactive", "true");
  vds1.setParameter("store", "true");
  vds1.initializeCommonParameters();

  os.str("");
  os << file2 << ", shift=0,0 deg";
  vds2.setParameter("file", os.str());
  vds2.setParameter("uvmax", "2000");
  vds2.setParameter("interactive", "true");
  vds2.setParameter("store", "true");
  vds2.initializeCommonParameters();

  vds1.storeDataInternallyOnReadin(true);
  vds1.loadData(false);

  vds2.storeDataInternallyOnReadin(true);
  vds2.loadData(false);

  ObsInfo& obs1 = vds1.obs_;
  ObsInfo& obs2 = vds2.obs_;

  std::vector<ObsInfo::Vis> vis;
  vis.resize(obs1.nGroup_ + obs2.nGroup_);

  double jdlast = obs1.visibilities_[obs1.visibilities_.size()-1].jd_;

  for(unsigned iVis=0; iVis < obs1.nGroup_; iVis++)
    vis[iVis] = obs1.visibilities_[iVis];

  double jdfirst = obs2.visibilities_[0].jd_;
  for(unsigned iVis=0; iVis < obs2.nGroup_; iVis++) {
    double djd = obs2.visibilities_[iVis].jd_ - jdfirst;
    obs2.visibilities_[iVis].jd_ = jdlast + djd;
    vis[obs1.nGroup_ + iVis] = obs2.visibilities_[iVis];
  }
  
  obs1.setNumberOfGroups(obs1.nGroup_ + obs2.nGroup_);
  obs1.visibilities_ = vis;

  FitsIoHandler fitsio;
  fitsio.setFirstTelescopeNum(1);
  fitsio.writeUvfFile("test.uvf", obs1);
}

/**.......................................................................
 * Test the chi-squared calculation
 */
void testCoadd(std::string& file)
{
  VisDataSetUvf vds1(0);
  VisDataSetUvf vds2(0);
  VisDataSetUvf vds3(0);

  vds1.setName("vds1");
  vds2.setName("vds2");
  vds3.setName("vds3");

  std::ostringstream os;

  os.str("");
  os << file << ", shift=0,0 deg";
  vds1.setParameter("file", os.str());
  vds1.setParameter("uvmax", "2000");
  vds1.setParameter("interactive", "true");
  vds1.setParameter("store", "true");
  vds1.initializeCommonParameters();

  os.str("");
  //  os << file << ", shift=0,0 deg";
  os << file << ", shift=0.1,0.1 deg";
  vds2.setParameter("file", os.str());
  vds2.setParameter("uvmax", "2000");
  vds2.setParameter("interactive", "true");
  vds2.initializeCommonParameters();

  os.str("");
  //  os << file << ", shift=0,0 deg";
  os << file << ", shift=0.1,-0.1 deg";
  vds3.setParameter("file", os.str());
  vds3.setParameter("uvmax", "2000");
  vds3.setParameter("interactive", "true");
  vds3.initializeCommonParameters();

  vds1.storeDataInternallyOnReadin(true);
  vds1.loadData(false);
  vds1.writeUvfFile("test2.uvf");

  vds2.loadData(false);
  vds3.loadData(false);

  PgUtil::open("/xs");
  PgUtil::setInteractive(true);

  //  vds1 += vds2;

  vds1.display();

  PgUtil::setZmin(0.0);
  PgUtil::setZmax(0.0);

  vds1 += vds2;
  vds1 += vds3;

  vds1.displayPrimaryBeams();

  vds1.display();
}

void testChisqGauss(VisDataSet& vds)
{
  // Define a model component
  
  GaussianClusterModel cluster;
  Temperature tmax;
  Angle sigma;
  Angle xoff;
  Angle yoff;
  double arcminval = 3.0;
  ChisqVariate chisq;

  // Compute chisq

  Timer timer;

  //------------------------------------------------------------
  // First iteration with actual model
  //------------------------------------------------------------

  timer.start();

  tmax.setMicroK(5000);

  Flux fnorm;
  fnorm.setJy(1.0);
  //  cluster.setNormalization(fnorm);

  sigma.setArcMinutes(0.001);
  cluster.setSigma(sigma);

  xoff.setArcMinutes(60.0/512 * 5);
  yoff.setArcMinutes(0.0);

  cluster.setXOffset(xoff);
  cluster.setYOffset(yoff);

  vds.addModel(cluster);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte());
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");

  //------------------------------------------------------------
  // First iteration with -X offset instead
  //------------------------------------------------------------

  timer.start();

  tmax.setMicroK(5000);
  //  cluster.setNormalization(tmax);
  sigma.setArcMinutes(2.0);
  cluster.setSigma(sigma);

  xoff.setArcMinutes(-arcminval);
  yoff.setArcMinutes(arcminval);
  cluster.setXOffset(xoff);
  cluster.setYOffset(yoff);

  vds.addModel(cluster);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte());
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");

  //------------------------------------------------------------
  // Second iteration with -Y offset instead
  //------------------------------------------------------------

  timer.start();

  tmax.setMicroK(5000);
  //  cluster.setNormalization(tmax);
  sigma.setArcMinutes(2.0);
  cluster.setSigma(sigma);

  xoff.setArcMinutes(arcminval);
  yoff.setArcMinutes(-arcminval);
  cluster.setXOffset(xoff);
  cluster.setYOffset(yoff);

  vds.addModel(cluster);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte());
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");

  //------------------------------------------------------------
  // Third iteration with -X and -Y offset instead
  //------------------------------------------------------------

  timer.start();

  tmax.setMicroK(5000);
  //  cluster.setNormalization(tmax);
  sigma.setArcMinutes(2.0);
  cluster.setSigma(sigma);

  xoff.setArcMinutes(-arcminval);
  yoff.setArcMinutes(-arcminval);
  cluster.setXOffset(xoff);
  cluster.setYOffset(yoff);

  vds.addModel(cluster);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte());
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");
}

/**.......................................................................
 * Test the chi-squared calculation
 */
void testChisqCalcPtSrc(std::string& file, double perc)
{
  VisDataSetUvf vds(0);

  // Determine what we can from the UVF file

  vds.initializeFromFile(file);

  // Now iterate through the UVF file, reading each group to determine
  // how many baselines actually belong to each baseline group, and
  // what the maximum UV radius encountered is.

  vds.countData(file);

  // Here's where we actually read the data

  Image image;
  image.xAxis().setAngularSize(Angle(Angle::Degrees(), 4.0));
  image.yAxis().setAngularSize(Angle(Angle::Degrees(), 4.0));
  image.xAxis().setNpix(256);
  image.yAxis().setNpix(256);

  vds.loadData(file, image);

  // Now calculate primary beams

  vds.computePrimaryBeams();

  vds.display();

  // Define a model component

  PtSrcModel ptsrc;

  Flux flux;
  flux.setJy(1.0);
  ptsrc.setFlux(flux);

  Angle xoff, yoff;
  xoff.setArcMinutes(4*60.0/256 * 10);
  yoff.setArcMinutes(0);

  ptsrc.setOffset(xoff, yoff);
  ptsrc.setSpectralIndex(0.0);

  ChisqVariate chisq;

  // Compute chisq

  Timer timer;

  COUT("Here 0 xoff = " << xoff);

  //------------------------------------------------------------
  // Next iteration with something else
  //------------------------------------------------------------

  timer.start();

  xoff.setArcMinutes(4*60.0/256 * 6);
  yoff.setArcMinutes(0);

  ptsrc.setOffset(xoff, yoff);

  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff);
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");

  COUT("Here 4");

  //------------------------------------------------------------
  // Next iteration with something else
  //------------------------------------------------------------

  timer.start();

  xoff.setArcMinutes(4*60.0/256 * 7);
  yoff.setArcMinutes(0);

  ptsrc.setOffset(xoff, yoff);

  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff);
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");

  COUT("Here 4");


  //------------------------------------------------------------
  // Next iteration with something else
  //------------------------------------------------------------

  timer.start();

  xoff.setArcMinutes(4*60.0/256 * 8);
  yoff.setArcMinutes(0);

  ptsrc.setOffset(xoff, yoff);

  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff);
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");

  COUT("Here 4");

  //------------------------------------------------------------
  // Next iteration with something else
  //------------------------------------------------------------

  timer.start();

  xoff.setArcMinutes(4*60.0/256 * 9);
  yoff.setArcMinutes(0);

  ptsrc.setOffset(xoff, yoff);

  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff);
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");

  COUT("Here 4");

  //------------------------------------------------------------
  // Next iteration with something else
  //------------------------------------------------------------

  timer.start();

  xoff.setArcMinutes(4*60.0/256 * 10);
  yoff.setArcMinutes(0);

  ptsrc.setOffset(xoff, yoff);

  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff);
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");

  COUT("Here 4");

  //------------------------------------------------------------
  // Next iteration with something else
  //------------------------------------------------------------

  timer.start();

  xoff.setArcMinutes(4*60.0/256 * 11);
  yoff.setArcMinutes(0);

  ptsrc.setOffset(xoff, yoff);

  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff);
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");

  COUT("Here 4");

  //------------------------------------------------------------
  // Next iteration with something else
  //------------------------------------------------------------

  timer.start();

  xoff.setArcMinutes(4*60.0/256 * 12);
  yoff.setArcMinutes(0);

  ptsrc.setOffset(xoff, yoff);

  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff);
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");

  COUT("Here 4");

  //------------------------------------------------------------
  // Next iteration with something else
  //------------------------------------------------------------

  timer.start();

  xoff.setArcMinutes(4*60.0/256 * 13);
  yoff.setArcMinutes(0);

  ptsrc.setOffset(xoff, yoff);

  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff);
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");

  //------------------------------------------------------------
  // Next iteration with something else
  //------------------------------------------------------------

  timer.start();

  xoff.setArcMinutes(4*60.0/256 * 14);
  yoff.setArcMinutes(0);

  ptsrc.setOffset(xoff, yoff);

  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff);
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");

  //------------------------------------------------------------
  // Next iteration with something else
  //------------------------------------------------------------

  timer.start();

  xoff.setArcMinutes(4*60.0/256 * 15);
  yoff.setArcMinutes(0);

  ptsrc.setOffset(xoff, yoff);

  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff);
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");

  //------------------------------------------------------------
  // Next iteration with something else
  //------------------------------------------------------------

  timer.start();

  xoff.setArcMinutes(4*60.0/256 * 16);
  yoff.setArcMinutes(0);

  ptsrc.setOffset(xoff, yoff);

  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff);
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");
}

void testSimulator(std::string& file, double perc)
{
  VisDataSetUvf vds;

  vds.setupForSimulation(true);

  // Determine what we can from the UVF file

  vds.initializeFromFile(file);

  // Now iterate through the UVF file, reading each group to determine
  // how many baselines actually belong to each baseline group, and
  // what the maximum UV radius encountered is.

  vds.countData(file);

  // Here's where we actually read the data

  COUT("Loading data:");
  vds.loadData(file, perc);

#if 1
  COUT("About to write UVF file");
  vds.writeUvfFile("uvfTest.uvf");
#endif

  // Test the simulation routine

  Lla lla;
  lla.setCoordSystem(COORD_GEOCENTRIC);
  lla.latitude_ = Angle(Angle::Degrees(), "-20:00:00");
  vds.getObs().setArrayLocation(lla);

  HourAngle start, stop, delta;

  start.setHours(-3);
  stop.setHours(3);
  delta.setHours(0.1);

  vds.getObs().setObsHa(start, stop, delta);

  COUT("Here 0");

  Image image1;
  image1.setNpix(512);
  image1.setAngularSize(Angle(Angle::Degrees(), 1.0));
  image1.setUnits(Unit::UNITS_JY);
  Angle off(Angle::Degrees(), 0.1);
  image1.createDeltaFunctionImage(2, off, off);

  Image image2;
  image2.setNpix(512);
  image2.setAngularSize(Angle(Angle::Degrees(), 1.0));
  image2.setUnits(Unit::UNITS_JY);
  Angle sigma(Angle::ArcMinutes(), 0.35);
  Angle off2(Angle::Degrees(), -0.1);
  image2.createGaussianImage(2, sigma, off2, off2);

  Image image3;
  image3.setNpix(512);
  image3.setAngularSize(Angle(Angle::Degrees(), 1.0));
  image3.setUnits(Unit::UNITS_JY);
  image3.createDeltaFunctionImage(1);

  Image image4;
  image4.setNpix(512);
  image4.setAngularSize(Angle(Angle::Degrees(), 1.0));
  image4.setUnits(Unit::UNITS_K);
  sigma.setArcMinutes(2.0);
  image4.createGaussianImage(-0.5, sigma);

  COUT("About to add image");

  //  vds.addImage(image1);
  //  vds.addImage(image2);
  //  vds.addImage(image3);
  vds.addImage(image4);

  // Now calculate primary beams

  COUT("Calculating primary beams:");
  vds.computePrimaryBeams();

  COUT("About to observe");
  vds.observe();

  COUT("About to write UVF file");
  vds.writeUvfFile("uvfTest.uvf");
}

void testReplacePtSrc(std::string& file, double perc)
{
  //------------------------------------------------------------
  // Test creating images of different sizes
  //------------------------------------------------------------

  //------------------------------------------------------------
  // Create an image from a point source model
  //------------------------------------------------------------

  Angle xOff, yOff;
  //xOff.setArcMinutes(30*0.117188);
  //  yOff.setArcMinutes(30*0.117188);
  xOff.setArcMinutes(0);
  yOff.setArcMinutes(0);

  Frequency freq;
  freq.setGHz(30.0);

  // Define a model component

  PtSrcModel ptsrc;
  ptsrc.setSpectralIndex(0.0);
  ptsrc.setGHz(30.0);
  ptsrc.setOffset(xOff, yOff);
  ptsrc.setJy(1.0);

  Image image;
  image.setNpix(1024);
  image.xAxis().setAngularSize(Angle(Angle::Degrees(), 2.0));
  image.yAxis().setAngularSize(Angle(Angle::Degrees(), 2.0));

  ptsrc.fillSzImage(image, freq);

  image.display();

  //------------------------------------------------------------
  // Now we have an image of the offset point source
  //------------------------------------------------------------

  VisDataSetUvf vds;

  vds.setupForSimulation(true);

  // Determine what we can from the UVF file

  vds.initializeFromFile(file);

  // Now iterate through the UVF file, reading each group to determine
  // how many baselines actually belong to each baseline group, and
  // what the maximum UV radius encountered is.

  vds.countData(file);

  // Here's where we actually read the data

  vds.loadData(file, image);

  //------------------------------------------------------------
  // Now replace internally stored visibilities with simulated ones
  //------------------------------------------------------------

  COUT("About to add image");

  vds.addImage(image);

  // Now calculate primary beams

  vds.computePrimaryBeams();

  Flux noiseRms;
  noiseRms.setMilliJy(0.1);

#if 0
  vds.getObs().setFixedNoiseRms(noiseRms);
  vds.getObs().setNoiseType(ObsInfo::NOISE_FIXED);
  vds.getObs().seed(1);
#endif
  vds.observe();

#if 1
  vds.writeUvfFile("uvfTest.uvf");
  return;
#endif

  // Now that we have internally replaced visibilities with simulated
  // vis, 'load' the data again as if reading from an external file,
  // to grid the simulated visibilities

  vds.loadData(file, image);

  // Now re-calculate primary beams to the size that the gridded data
  // will support

  vds.computePrimaryBeams();

  //------------------------------------------------------------
  // Test computing chisq for different point source models
  //------------------------------------------------------------

  ptsrc.setJy(1.0);
  xOff.setArcMinutes(10*0.117188);
  yOff.setArcMinutes(0.0);
  computeChisqPtSrc(xOff, yOff, ptsrc, vds);

  ptsrc.setJy(0.9);
  computeChisqPtSrc(xOff, yOff, ptsrc, vds);

  ptsrc.setJy(0.7);
  computeChisqPtSrc(xOff, yOff, ptsrc, vds);

  ptsrc.setJy(0.5);
  computeChisqPtSrc(xOff, yOff, ptsrc, vds);

  ptsrc.setJy(0.2);
  computeChisqPtSrc(xOff, yOff, ptsrc, vds);

  ptsrc.setJy(1.1);
  computeChisqPtSrc(xOff, yOff, ptsrc, vds);

  ptsrc.setJy(1.0);
  yOff.setArcMinutes(0.0);
  computeChisqPtSrc(xOff, yOff, ptsrc, vds);

  ptsrc.setJy(0.5);
  computeChisqPtSrc(xOff, yOff, ptsrc, vds);

  xOff.setArcMinutes(1.5);
  computeChisqPtSrc(xOff, yOff, ptsrc, vds);

  xOff.setArcMinutes(1.17188);
  computeChisqPtSrc(xOff, yOff, ptsrc, vds);
}

void computeChisqPtSrc(Angle& xOff, Angle& yOff, PtSrcModel& ptsrc, VisDataSet& vds)
{
  Timer timer;
  timer.start();
  ptsrc.setOffset(xOff, yOff);
  vds.addModel(ptsrc);
  gcp::util::ChisqVariate chisq = vds.computeChisq();
  timer.stop();

  COUT("xOff = " << xOff << " yOff = " << yOff << " Chisq = " << chisq << " pte = " << chisq.pte());
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");
}

 
void testReplaceGauss(std::string& file, double perc)
{
  COUT("Here 0");
  VisDataSetUvf vds;
  vds.setupForSimulation(true);

  // Determine what we can from the UVF file

  COUT("Here 1");
  vds.initializeFromFile(file);

  // Now iterate through the UVF file, reading each group to determine
  // how many baselines actually belong to each baseline group, and
  // what the maximum UV radius encountered is.

  COUT("Here 2");
  vds.countData(file);

  // Here's where we actually read the data

  COUT("Here 3");
  vds.loadData(file, perc);

  // Now replace internally stored visibilities with simulated ones

  Image image;
  image.setNpix(512);
  image.setAngularSize(Angle(Angle::Degrees(), 1.0));
  image.setUnits(Unit::UNITS_UK);
  Angle sigma;

  sigma.setArcMinutes(2);
  image.createGaussianImage(-5000, sigma);

  vds.addImage(image);

  // Now calculate primary beams

  COUT("Abotu to compute primaryu beams");
  vds.computePrimaryBeams();

  Flux noiseRms;
  noiseRms.setMilliJy(10.0);

#if 0
  vds.getObs().setFixedNoiseRms(noiseRms);
  vds.getObs().setNoiseType(ObsInfo::NOISE_FIXED);
  vds.getObs().seed(1);
#endif
  vds.observe();

#if 0
  vds.writeUvfFile("uvfTest.uvf");
  return;
#endif

  // Now that we have internally replaced visibilities with simulated
  // vis, 'load' the data again as if reading from an external file

  vds.loadData(file, perc);

  // Now re-calculate primary beams to the size that the gridded data
  // will support

  vds.computePrimaryBeams();

  // Define a model component

  GaussianClusterModel cluster;
  Temperature tmax;

  // Compute chisq

  Timer timer;

  tmax.setMicroK(-5000);
  //  cluster.setNormalization(tmax);
  sigma.setArcMinutes(2);
  cluster.setSigma(sigma);

  timer.start();
  vds.addModel(cluster);
  timer.stop();

  ChisqVariate chisq = vds.computeChisq();



  COUT("Chisq = " << chisq << " pte = " << chisq.pte());
  COUT("Add model required: " << timer.deltaInSeconds() << " seconds");

  tmax.setMicroK(-5000);
  //  cluster.setNormalization(tmax);
  sigma.setArcMinutes(1);
  cluster.setSigma(sigma);

  timer.start();
  vds.addModel(cluster);
  timer.stop();

  chisq = vds.computeChisq();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte());
  COUT("Add model required: " << timer.deltaInSeconds() << " seconds");

  tmax.setMicroK(-8000);
  //  cluster.setNormalization(tmax);
  sigma.setArcMinutes(2);
  cluster.setSigma(sigma);

  timer.start();
  vds.addModel(cluster);
  timer.stop();

  chisq = vds.computeChisq();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte());
  COUT("Add model required: " << timer.deltaInSeconds() << " seconds");
}

void testReplaceGaussMultiThread(std::string& file, double perc, unsigned nthread)
{
  ThreadPool pool(nthread);
  pool.spawn();

  VisDataSetUvf vds(&pool);
  vds.setupForSimulation(true);

  // Determine what we can from the UVF file

  CTOUT("Here 0");
  vds.initializeFromFile(file);

  // Now iterate through the UVF file, reading each group to determine
  // how many baselines actually belong to each baseline group, and
  // what the maximum UV radius encountered is.

  COUT("Here 1");
  vds.countData(file);

  // Here's where we actually read the data

  COUT("Here 2");
  vds.loadData(file, perc);

  // Now replace internally stored visibilities with simulated ones

  unsigned padfac = 2.0;
  Image image;
  image.setNpix(padfac * 512);
  image.setAngularSize(Angle(Angle::Degrees(), padfac * 1.9));
  image.setUnits(Unit::UNITS_UK);

  GaussianClusterModel gauss, gauss2;


  Angle xoff, yoff, off;
  Temperature tnorm;
  Angle majSig, minSig, rotAng;

  xoff.setArcMinutes(2);
  yoff.setArcMinutes(3);
  //xoff.setArcMinutes(0);
  //yoff.setArcMinutes(0.0);

  gauss.setXOffset(xoff);
  gauss.setYOffset(yoff);

  Flux fnorm;
  fnorm.setJy(0.2);

  tnorm.setK(1);
  //  gauss.setNormalization(tnorm);
  //  gauss.setNormalization(fnorm);

  minSig.setArcMinutes(2);
  majSig.setArcMinutes(2);

  gauss.setMinSigma(minSig);
  gauss.setMajSigma(majSig);

  rotAng.setDegrees(0);
  gauss.setRotationAngle(rotAng);

  Frequency freq;
  freq.setGHz(30.0);

  gauss.fillSzImage(image, freq);
  image.display();

  Image image2 = image;

  xoff.setArcMinutes(0);
  yoff.setArcMinutes(0.0);

#if 0
  gauss2 = gauss;
  gauss2.setXOffset(xoff);
  gauss2.setYOffset(yoff);
  
  gauss2.fillSzImage(image2, freq);
  image2.display();
#endif

  unsigned imind;

  for(unsigned i=0; i < padfac*512; i++) {
    for(unsigned j=0; j < padfac*512; j++) {

      imind = j * padfac*512 + i;

      if(i < (padfac-1)*(padfac*512)/(2*padfac) || i >= (padfac+1)*(padfac*512)/(2*padfac)) {
	image.data_[imind] = 0.0;
      }

      if(j < (padfac-1)*(padfac*512)/(2*padfac) || j >= (padfac+1)*(padfac*512)/(2*padfac)) {
	image.data_[imind] = 0.0;
      }

    }
  }

  image.display();

#if 0
  // Stuff a point-source model in instead

  PtSrcModel ptsrc;
  Flux flux;
  flux.setJy(1.0);
  ptsrc.setFlux(flux);
  xoff.setArcMinutes(1.17188);
  yoff.setArcMinutes(0.0);
  ptsrc.setOffset(xoff, yoff);
  ptsrc.setSpectralIndex(0.0);
  Frequency freq;
  freq.setGHz(30.0);
  ptsrc.fillSzImage(image, freq);
#endif

  image.display();

  COUT("Here 3");
  vds.addImage(image);

  // Now calculate primary beams

  COUT("Here 4");
  vds.computePrimaryBeams();

  Flux noiseRms;
  noiseRms.setMilliJy(25);

#if 0
  vds.getObs().setFixedNoiseRms(noiseRms);
  vds.getObs().setNoiseType(ObsInfo::NOISE_FIXED);
  vds.getObs().seed(1);
#endif
  COUT("Here 5");
  vds.observe();

#if 1
  vds.writeUvfFile("uvfTest.uvf");
  return;
#endif

  // Now that we have internally replaced visibilities with simulated
  // vis, 'load' the data again as if reading from an external file

  vds.loadData(file, perc);

  // Now re-calculate primary beams to the size that the gridded data
  // will support

  vds.computePrimaryBeams();

  testChisqGauss(vds);
}

#if 0
  // Define a model component

  GaussianClusterModel cluster;
  Temperature tmax;

  // Compute chisq

  Timer timer;

  timer.start();

  tmax.setMicroK(5000);
//  cluster.setNormalization(tmax);
  sigma.setArcMinutes(2);
  cluster.setSigma(sigma);

  off.setArcMinutes(5);
  cluster.setXOffset(off);
  cluster.setYOffset(off);
  vds.addModel(cluster);
  ChisqVariate chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte());
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");

  timer.start();

  tmax.setMicroK(5000);
//  cluster.setNormalization(tmax);
  sigma.setArcMinutes(2);
  cluster.setSigma(sigma);

  xoff.setArcMinutes(-5);
  yoff.setArcMinutes(5);
  cluster.setXOffset(xoff);
  cluster.setYOffset(yoff);
  vds.clearModel();
  vds.addModel(cluster);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte());
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");

  timer.start();

  tmax.setMicroK(5000);
//  cluster.setNormalization(tmax);
  sigma.setArcMinutes(2);
  cluster.setSigma(sigma);

  xoff.setArcMinutes(-5);
  yoff.setArcMinutes(-5);
  cluster.setXOffset(xoff);
  cluster.setYOffset(yoff);

  vds.clearModel();
  vds.addModel(cluster);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte());
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");

  timer.start();

  tmax.setMicroK(5000);
//  cluster.setNormalization(tmax);
  sigma.setArcMinutes(2);
  cluster.setSigma(sigma);

  xoff.setArcMinutes(-5);
  yoff.setArcMinutes(5);
  cluster.setXOffset(xoff);
  cluster.setYOffset(yoff);

  vds.clearModel();
  vds.addModel(cluster);
  chisq = vds.computeChisq();

  timer.stop();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte());
  COUT("Computation required: " << timer.deltaInSeconds() << " seconds");
}
#endif

void testReader(std::string& file)
{
  gcp::util::FitsUvfReader reader(file);

  COUT("RA  = " << reader.obsra());
  COUT("DEC = " << reader.obsdec());

  std::vector<double> freqs = reader.freqs();
  double freq0 = freqs[0];

  gcp::util::FitsBinTableReader tableReader(file);

  bool stop=false;
  while(!stop) {
    try {
      tableReader.getNextTableInfo();
      COUT("Found table: " << tableReader.tableName());

      for(int iCol=tableReader.nCol()-1; iCol >= 0; iCol--) {
	COUT("column = " << tableReader.colName(iCol));

	if(tableReader.colName(iCol) == "IF FREQ") {
	  COUT(tableReader.colDataType(iCol));
	  std::vector<double> vals = tableReader.getDoubleData(iCol);
	  for(unsigned i=0; i < vals.size(); i++) {
	    COUT(vals[i] + freq0);
	  }
	}

	if(tableReader.colName(iCol) == "STATION") {
	  std::vector<string> vals = tableReader.getStringData(iCol);
	  for(unsigned i=0; i < vals.size(); i++) {
	    COUT(vals[i]);
	  }
	}

	if(tableReader.colName(iCol) == "ANNAME") {
	  std::vector<string> vals = tableReader.getStringData(iCol);
	  for(unsigned i=0; i < vals.size(); i++) {
	    COUT(vals[i]);
	  }
	}

	if(tableReader.colName(iCol) == "MNTSTA") {
	  std::vector<string> vals = tableReader.getStringData(iCol);
	  for(unsigned i=0; i < vals.size(); i++) {
	    COUT(vals[i]);
	  }
	}

      }

    } catch(...) {
      stop = true;
    }
  }

  FitsUvfReader::Vis vis;
  reader.readData(0, vis);
  COUT(vis);
  //  for(unsigned iGroup=0; iGroup < reader.nGroup(); iGroup++) {
  //    reader.readData(iGroup, vis);
  //  }

  std::vector<double> ifs = reader.ifs();
  COUT("Number of IFs: " << reader.nIf());
}

/**.......................................................................
 * Test calculating Chi-squared against simulated data
 */
void testSimChisq(std::string& file, double perc)
{
  // We want to create a VisDataSet object, then install a model to
  // generate simulated visibilities.  Then compute visibilities on
  // the same grid as the 

  VisDataSetUvf vds;

  vds.setupForSimulation(true);

  // Determine what we can from the UVF file

  vds.initializeFromFile(file);

  // Now iterate through the UVF file, reading each group to determine
  // how many baselines actually belong to each baseline group, and
  // what the maximum UV radius encountered is.

  vds.countData(file);

  // Here's where we actually read the data

  COUT("Loading data:");
  vds.loadData(file, perc);
}

void testInit(std::string file)
{
  VisDataSetUvf uvf;
  uvf.setParameter("file", file);
  uvf.loadDataMultiple(false);
}
