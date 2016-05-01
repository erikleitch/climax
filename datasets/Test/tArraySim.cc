#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/models/GaussianClusterModel.h"
#include "gcp/models/PtSrcModel.h"

#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/IoLock.h"
#include "gcp/util/Exception.h"
#include "gcp/util/FitsBinTableReader.h"
#include "gcp/util/FitsUvfReader.h"

#include "gcp/pgutil/PgUtil.h"

#include "gcp/fftutil/Dft2d.h"
#include "gcp/fftutil/Image.h"
#include "gcp/fftutil/ObsInfo.h"

#include "gcp/datasets/VisDataSetUvf.h"

using namespace std;
using namespace gcp::datasets;
using namespace gcp::models;
using namespace gcp::program;
using namespace gcp::util;

KeyTabEntry Program::keywords[] = {
  { "dev",      "/xs",            "s", "Pgplot device"},
  { "file",     "",               "s", "FITS file to write out"},
  { "zeropad",  "f",              "b", "Zeropad the array?"},
  { "phase",    "0.0",            "d", "Phase (degrees)"},
  { "period",   "64",             "d", "Period (pixels)"},
  { "freq",     "30",             "d", "freq (GHz)"},
  { "perc",     "0.98",           "d", "percent correlation"},
  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

void testReader(std::string& file);
void setupObservation(VisDataSetUvf& vds);
void addAnt(ObsInfo& obs, Length& east, Length& north, Length& up, unsigned index);
void observeModel(VisDataSetUvf& vds);
void calcChisqSim(VisDataSetUvf& vds);
void calcChisqFile(std::string fileName);
void calcChisqModel(VisDataSetUvf& vds);
void fillImageWithModel(Image& image);

void calcChisqModelPtSrc(VisDataSetUvf& vds);
void fillImageWithModelPtSrc(Image& image);

void calcChisqModelGauss(VisDataSetUvf& vds);
void fillImageWithModelGauss(Image& image);

/**.......................................................................
 * Main -- Generate a pure simulated UVF file
 */
int Program::main()
{
  VisDataSetUvf vds;
  setupObservation(vds);
  observeModel(vds);

  vds.writeUvfFile(Program::getStringParameter("file"));

  return 0;

  calcChisqSim(vds);
  calcChisqFile(Program::getStringParameter("file"));

  return 0;
}

void calcChisqFile(std::string fileName)
{
  VisDataSetUvf vds;

  Image image;
  image.setNpix(256);
  image.setAngularSize(Angle(Angle::Degrees(), 4.0));

  vds.initializeFromFile(fileName);
  vds.countData(fileName);
  vds.loadData(fileName, image);

  calcChisqModel(vds);
}


void calcChisqSim(VisDataSetUvf& vds)
{
  // Now regrid the fake data, as if reading it in

  vds.setupForSimulation(true);

  Image image;
  image.setNpix(256);
  image.setAngularSize(Angle(Angle::Degrees(), 4.0));

  vds.loadData("", image);

  calcChisqModel(vds);
}

void calcChisqModel(VisDataSetUvf& vds)
{
  calcChisqModelPtSrc(vds);
  //calcChisqModelGauss(vds);
}

void calcChisqModelGauss(VisDataSetUvf& vds)
{
  GaussianClusterModel model;
  Angle xoff, yoff, off;
  
  xoff.setArcMinutes(4*60.0/256 * 0);
  yoff.setArcMinutes(4*60.0/256 * 0);

  model.setXOffset(xoff);
  model.setYOffset(yoff);

  Angle minSig, majSig;
  minSig.setArcMinutes(1.0);
  majSig.setArcMinutes(2.0);

  model.setMinSigma(minSig);
  model.setMajSigma(majSig);

  Angle rot;
  rot.setDegrees(30);
  model.setRotationAngle(rot);

#if 0
  Flux fnorm;
  fnorm.setJy(0.2);
  //  model.setNormalization(fnorm);
#else
  Temperature tnorm;
  tnorm.setMicroK(-800);
  //  model.setNormalization(tnorm);
#endif

  ChisqVariate chisq;

  vds.addModel(model);
  chisq = vds.computeChisq();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff << " yoff = " << yoff);

  minSig.setArcMinutes(1.0);
  model.setMinSigma(minSig);
  vds.addModel(model);
  chisq = vds.computeChisq();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff << " yoff = " << yoff);

  minSig.setArcMinutes(3.0);
  model.setMinSigma(minSig);
  vds.addModel(model);
  chisq = vds.computeChisq();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff << " yoff = " << yoff);

  minSig.setArcMinutes(1.9);
  model.setMinSigma(minSig);
  vds.addModel(model);
  chisq = vds.computeChisq();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff << " yoff = " << yoff);

  minSig.setArcMinutes(2.01);
  model.setMinSigma(minSig);
  vds.addModel(model);
  chisq = vds.computeChisq();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff << " yoff = " << yoff);

  yoff.setArcMinutes(-4*60.0/256 * 5);
  model.setYOffset(yoff);
  minSig.setArcMinutes(2.01);
  model.setMinSigma(minSig);
  vds.addModel(model);
  chisq = vds.computeChisq();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff << " yoff = " << yoff);
}

void calcChisqModelPtSrc(VisDataSetUvf& vds)
{
  COUT("Done loading data");

  PtSrcModel ptsrc;

  Flux flux;
  flux.setJy(0.2);
  ptsrc.setFlux(flux);
  ptsrc.setSpectralIndex(0.0);

  Angle xoff, yoff;
  ChisqVariate chisq;

  yoff.setArcMinutes(4*60.0/256 * 5);

  xoff.setArcMinutes(4*60.0/256 * 8);
  ptsrc.setOffset(xoff, yoff);
  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff);

  xoff.setArcMinutes(4*60.0/256 * 9);
  ptsrc.setOffset(xoff, yoff);
  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff);

  xoff.setArcMinutes(4*60.0/256 * 10);
  ptsrc.setOffset(xoff, yoff);
  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff);

  xoff.setArcMinutes(4*60.0/256 * 11);
  ptsrc.setOffset(xoff, yoff);
  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff);

  xoff.setArcMinutes(4*60.0/256 * 20);
  ptsrc.setOffset(xoff, yoff);
  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff);

  xoff.setArcMinutes( 4*60.0/256 * 10);
  yoff.setArcMinutes(-4*60.0/256 * 5);
  ptsrc.setOffset(xoff, yoff);
  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff  << " " << yoff);

  xoff.setArcMinutes(-4*60.0/256 * 10);
  yoff.setArcMinutes(-4*60.0/256 * 5);
  ptsrc.setOffset(xoff, yoff);
  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff  << " " << yoff);

  xoff.setArcMinutes(-4*60.0/256 * 10);
  yoff.setArcMinutes( 4*60.0/256 * 5);
  ptsrc.setOffset(xoff, yoff);
  vds.addModel(ptsrc);
  chisq = vds.computeChisq();

  COUT("Chisq = " << chisq << " pte = " << chisq.pte() << " " << xoff << " " << yoff);
}

void setupObservation(VisDataSetUvf& vds)
{
  ObsInfo& obs = vds.getObs();

  unsigned nAnt = 8;
  
  Lla lla;

  Angle lng, lat;
  Length alt;

  lng.setDegrees(-100);
  lat.setDegrees(35);
  alt.setFeet(4000);

  lla.longitude_ = lng;
  lla.latitude_  = lat;
  lla.altitude_  = alt;

  obs.setArrayLocation(lla);
  obs.setNumberOfAntennas(11);
  obs.setAntennaType(Antenna::ANT_SZA);

  Length east, north, up;
  
  up.setMeters(0.0);

  east.setMeters(0.0);
  north.setMeters(0.0);
  addAnt(obs, east, north, up, 0);

  east.setMeters(0.0);
  north.setMeters(3.5);
  addAnt(obs, east, north, up, 1);

  east.setMeters(0.0);
  north.setMeters(-3.55);
  addAnt(obs, east, north, up, 2);

  east.setMeters(3.52);
  north.setMeters(0.0);
  addAnt(obs, east, north, up, 3);

  east.setMeters(-3.51);
  north.setMeters(0.0);
  addAnt(obs, east, north, up, 4);

  east.setMeters(3.54);
  north.setMeters(3.5);
  addAnt(obs, east, north, up, 5);

  east.setMeters(3.53);
  north.setMeters(-3.52);
  addAnt(obs, east, north, up, 6);

  east.setMeters(-3.51);
  north.setMeters(3.51);
  addAnt(obs, east, north, up, 7);

  east.setMeters(-3.57);
  north.setMeters(-3.58);
  addAnt(obs, east, north, up, 8);

  east.setMeters(-2*3.57);
  north.setMeters(0.0);
  addAnt(obs, east, north, up, 9);

  east.setMeters(2*3.57);
  north.setMeters(0.0);
  addAnt(obs, east, north, up, 10);

  obs.plotAntennas();

  obs.setSourceName("A1316");

  HourAngle ra;
  //  ra.setHours(12.0 - 0.005);
  ra.setHours(12.0);

  Declination dec;
  //  dec.setDegrees(35 - 2.0/60);
  dec.setDegrees(35);

  obs.setObsRa(ra);
  obs.setObsDec(dec);
  obs.setObsEquinox(2000);

  vds.setRa(ra);
  vds.setDec(dec);

  HourAngle start, stop, delta;
  start.setHours(-6);
  stop.setHours(6);
  delta.setHours(0.1);

  obs.setObsHa(start, stop, delta);
  obs.setTelescopeName("SZA");
  obs.setInstrumentName("CLIMAX");

  unsigned nFreq = 16;
  std::vector<Frequency> freqs(nFreq);

  for(unsigned iFreq=0; iFreq < nFreq; iFreq++) {
    freqs[iFreq].setGHz(27 + iFreq*0.5);
  }

  std::vector<Frequency> bws(nFreq);

  for(unsigned iFreq=0; iFreq < nFreq; iFreq++) {
    bws[iFreq].setMHz(500);
  }

  obs.setFrequencyInformation(freqs, bws);
  obs.setNumberOfStokesParameters(1);
}


void addAnt(ObsInfo& obs, Length& east, Length& north, Length& up, unsigned index)
{
  LengthTriplet enu;

  enu.setCoordSystem(COORD_ENU);

  enu.east_  = east;
  enu.north_ = north;
  enu.up_    = up;

  obs.setAntennaLocation(enu, index);
}

void observeModel(VisDataSetUvf& vds) 
{
  Image image;
  image.setNpix(256);
  image.setAngularSize(Angle(Angle::Degrees(), 4.0));
  image.setRaDec(vds.ra_, vds.dec_);

  fillImageWithModel(image);

  image.display();
  
  // First we need to determine unique baseline groupings, as on read-in
  
  vds.determineUniqueBaselineGroupings(vds.getObs());
  
  // Only now can we add the images
  
  COUT("here 0");
  vds.addImage(image);
  
  // And calculate primary beams
  
  COUT("here 1");
  vds.computePrimaryBeams();
  
  // Now we can observe
  
  Flux noiseRms;
  noiseRms.setJy(0.005);

#if 0
  vds.getObs().setFixedNoiseRms(noiseRms);
  vds.getObs().setNoiseType(ObsInfo::NOISE_FIXED);
  vds.getObs().seed(1);
#endif
  vds.observe(vds.getObs());
}

void fillImageWithModel(Image& image)
{
  //  fillImageWithModelPtSrc(image);
  fillImageWithModelGauss(image);
}

void fillImageWithModelPtSrc(Image& image)
{
  PtSrcModel ptsrc;
  Angle xoff, yoff, off;
  
  xoff.setArcMinutes(4*60.0/256 * 5);
  yoff.setArcMinutes(4*60.0/256 * 0);

  ptsrc.setXOffset(xoff);
  ptsrc.setYOffset(yoff);

  Flux fnorm;
  fnorm.setJy(1.0);

  ptsrc.setFlux(fnorm);
  ptsrc.setSpectralIndex(0.0);

  Frequency freq;
  freq.setGHz(30.0);
  ptsrc.setFrequency(freq);

  ptsrc.fillSzImage(image, freq);
}

void fillImageWithModelGauss(Image& image)
{
  GaussianClusterModel model;
  Angle xoff, yoff, off;

  HourAngle modRa;
  modRa.setHours(12.0);

  Declination modDec;
  modDec.setDegrees(35);

  model.setRa(modRa);
  model.setDec(modDec);

  xoff.setArcMinutes(4*60.0/256 * 0);
  yoff.setArcMinutes(4*60.0/256 * 0);

  model.setXOffset(xoff);
  model.setYOffset(yoff);

  Angle minSig, majSig;
  minSig.setArcMinutes(1.0);
  majSig.setArcMinutes(2.0);

  model.setMinSigma(minSig);
  model.setMajSigma(majSig);

  Angle rot;
  rot.setDegrees(30);
  model.setRotationAngle(rot);

#if 0
  Flux fnorm;
  fnorm.setJy(0.2);
  //  model.setNormalization(fnorm);
#else
  Temperature tnorm;
  tnorm.setMicroK(-800);
  //  model.setNormalization(tnorm);
#endif

  Frequency freq;
  freq.setGHz(30.0);

  model.fillSzImage(image, freq);
}
