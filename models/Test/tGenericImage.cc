#include <iostream>
#include <math.h>

#include "gcp/program/Program.h"

#include "gcp/models/GenericImage.h"

#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"
#include "gcp/util/IoLock.h"
#include "gcp/util/Exception.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;
using namespace gcp::models;

KeyTabEntry Program::keywords[] = {
  { "dev",      "/xs",            "s", "Pgplot device"},
  { "file",     "/Users/eml/Desktop/Downloads/skv204978076416.fits",               "s", "FITS file to read in"},
  { "xdeg",       "1.5",         "s", "Size of the X-axis, in degrees"},
  { "ydeg",       "1.5",         "s", "Size of the Y-axis, in degrees"},

  { END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS,END_OF_KEYWORDS},
};

void Program::initializeUsage() {};

int Program::main()
{
  std::string file = Program::getStringParameter("file");
  std::string xdeg = Program::getStringParameter("xdeg");
  std::string ydeg = Program::getStringParameter("ydeg");

  GenericImage gi;

  gi.setParameter("file",  file);
  gi.setParameter("xsize", xdeg, "degrees");
  gi.setParameter("ysize", ydeg, "degrees");

  gi.getVar("norm")->setVal(-8000.0, "muK");
  gi.getVar("norm")->setUnits("muK");

  gi.initializeImage();
  gi.difmapDisplay();

  // Now create an image and fill it 

  Image image = gi.image_;

  Frequency freq;
  freq.setGHz(30);

  gi.fillImage(DataSetType::DATASET_2D, image, &freq);

  image.display();

  return 0;
}
