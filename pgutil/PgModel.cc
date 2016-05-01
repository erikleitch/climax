#include "gcp/pgutil/PgModel.h"
#include "gcp/pgutil/PgUtil.h"

#include "cpgplot.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
PgModel::PgModel() 
{
  drawCenter_ = true;
  fill_       = false;
}

/**.......................................................................
 * Destructor.
 */
PgModel::~PgModel() {}

PgModel::PgModel(PgModel& mod) {
  *this = mod;
}

PgModel::PgModel(const PgModel& mod) {
  *this = mod;
}

void PgModel::operator=(const PgModel& mod) {
  *this = (PgModel&)mod;
}

void PgModel::operator=(PgModel& mod) {

  drawCenter_ = mod.drawCenter_;
  fill_ = mod.fill_;
  type_ = mod.type_;
  xMid_ = mod.xMid_;
  yMid_ = mod.yMid_;

  peak_  = mod.peak_;
  xRad1_ = mod.xRad1_;
  yRad1_ = mod.yRad1_;

  xRad2_ = mod.xRad2_;
  yRad2_ = mod.yRad2_;
    
  rad1_ = mod.rad1_;
  rad2_ = mod.rad2_;

  radMin_ = mod.radMin_;
  radMax_ = mod.radMax_;

  rot_   = mod.rot_;
  angle_ = mod.angle_;

  x_ = mod.x_;
  y_ = mod.y_;

  xMin_ = mod.xMin_;
  xMax_ = mod.xMax_;
  yMin_ = mod.yMin_;
  yMax_ = mod.yMax_;
}

void PgModel::rectify() {

  float dx = xRad1_ - xMid_;
  float dy = yRad1_ - yMid_;

  rad1_ = sqrt(dx*dx + dy*dy);

  // Get the rotation angle from the dot product of the first vector with the x-axis

  if(dy < 0)
    rot_.setRadians(-acos(dx/rad1_));
  else
    rot_.setRadians(+acos(dx/rad1_));

  COUT("Rotation angle is now: " << rot_);

  if(!hasTwoAxes()) {
    rad2_  = rad1_;
    xRad2_ = xRad1_;
    yRad2_ = yRad1_;
  }

  dx = xRad2_ - xMid_;
  dy = yRad2_ - yMid_;
    
  rad2_ = sqrt(dx*dx + dy*dy);
    
  // Force the second axis perpendicular to the first one
    
  xRad2_ = -rad2_ * sin(rot_.radians()) + xMid_;
  yRad2_ =  rad2_ * cos(rot_.radians()) + yMid_;

  if(rad1_ < rad2_) {
    radMin_ = rad1_;
    radMax_ = rad2_;

    dx = xRad2_ - xMid_;
    COUT("Setting angle to yRad = " << yRad2_ << " rad = " << rad2_);
    if(dx < 0)
      angle_.setRadians(-acos((yRad2_ - yMid_)/rad2_));
    else
      angle_.setRadians(-acos((yRad2_ - yMid_)/rad2_));
      
  } else {
    radMin_ = rad2_;
    radMax_ = rad1_;

    dx = xRad1_ - xMid_;
    COUT("Setting angle to yRad = " << yRad1_ << " rad = " << rad1_);
    if(dx < 0)
      angle_.setRadians(-acos((yRad1_ - yMid_)/rad1_));
    else
      angle_.setRadians(acos((yRad1_ - yMid_)/rad1_));
  }
    
  precalculateShape();
}

void PgModel::precalculateShape() {
  x_.resize(500);
  y_.resize(500);

  double dtheta = 2*M_PI/(x_.size()-1);
  double x,y, theta, xrat, yrat;

  double cr = cos(rot_.radians());
  double sr = sin(rot_.radians());

  for(unsigned i=0; i < x_.size(); i++) {
    theta = dtheta * i;
    xrat = cos(theta);
    yrat = sin(theta);
    x = xrat * rad1_;
    y = yrat * rad2_;

    x_[i] = x * cr - y * sr + xMid_;
    y_[i] = x * sr + y * cr + yMid_;
  }
}

void PgModel::draw() 
{
  int ci;
  cpgqci(&ci);
  int lw;
  cpgqlw(&lw);
  float ch;
  cpgqch(&ch);

  if(!(type_ == TYPE_BOX || type_ == TYPE_BEAM)) {
    cpgsci(0);
    cpgsch(1.2*ch);
    cpgslw(10*lw);
    render();
  }

  cpgsci(getColor());
  cpgslw(lw);
  cpgsch(ch);
  render();

  cpgsci(ci);
};

void PgModel::render() 
{
  float x1,x2,y1,y2;
  cpgqwin(&x1,&x2,&y1,&y2);

  int ci;
  cpgqci(&ci);

  if(type_ == TYPE_DELTA) {

    if(!(xMid_ > x1 && xMid_ < x2 && yMid_ > y1 && yMid_ < y2))
      return;

    cpgpt(1, &xMid_, &yMid_, 12);

  } else if(type_ == TYPE_BOX)  {

    if(!(xMin_ > x1 && xMax_ < x2 && yMin_ > y1 && yMax_ < y2))
      return;

    cpgmove(xMin_, yMin_);
    cpgdraw(xMin_, yMax_);
    cpgdraw(xMax_, yMax_);
    cpgdraw(xMax_, yMin_);
    cpgdraw(xMin_, yMin_);

  } else {

    if(!(xMid_ > x1 && xMid_ < x2 && yMid_ > y1 && yMid_ < y2))
      return;

    float dx = xRad1_ - xMid_;
    float dy = yRad1_ - yMid_;

    if(drawCenter_) {
      cpgmove(xMid_, yMid_);
      cpgdraw(xMid_ + dx, yMid_ + dy);
      cpgmove(xMid_, yMid_);
      cpgdraw(xMid_ - dx, yMid_ - dy);
      
      dx = xRad2_ - xMid_;
      dy = yRad2_ - yMid_;
      
      cpgmove(xMid_, yMid_);
      cpgdraw(xMid_ + dx, yMid_ + dy);
      cpgmove(xMid_, yMid_);
      cpgdraw(xMid_ - dx, yMid_ - dy);
    }

    if(fill_) {
      int oldCi;
      cpgqci(&oldCi);

      cpgsci(15);
      cpgsfs(1);
      cpgpoly(x_.size(), &x_[0], &y_[0]);
      cpgsci(oldCi);
    }

    cpgline(x_.size(), &x_[0], &y_[0]);
  }
};

bool PgModel::hasTwoAxes() {
  return type_ == TYPE_ELBETA || type_ == TYPE_GAUSS;
};

int PgModel::getColor() {
  switch (type_) {
  case TYPE_ARNAUD:
    return 5;
    break;
  case TYPE_GAUSS:
    return 8;
    break;
  case TYPE_BETA:
    return 7;
    break;
  case TYPE_ELBETA:
    return 10;
    break;
  case TYPE_DELTA:
    return 6;
  case TYPE_BEAM:
    return 5;
    break;
  case TYPE_BOX:
    return 10;
    break;
  default:
    return 11;
    break;
  }
};

/**.......................................................................
 * Use current PgUtil callbacks to print model parameters
 */
void PgModel::print(std::string& unit, unsigned i, Trans& trans)
{
  ostringstream os;
  os << "mod" << i;

  std::string name;
  name = os.str();

  switch (type_) {
  case TYPE_DELTA:

    COUT("");
    COUT("addmodel type=ptsrc name=" << name << ";");
    COUT(name << ".Sradio = " << peak_ << " " << unit << ";");
    COUT(name << formatPosition("xoff", xMid_));
    COUT(name << formatPosition("yoff", yMid_));

    break;

  case TYPE_ARNAUD:

    COUT("");
    COUT("addmodel type=arnaudmodel name=" << name << ";");
    COUT(name << ".Sradio = " << peak_ << " " << unit << ";");
    COUT(name << formatPosition("xoff", xMid_));
    COUT(name << formatPosition("yoff", yMid_));
    COUT(name << formatPosition("thetaCore", rad1_));

    break;

  case TYPE_BETA:

    COUT("");
    COUT("addmodel type=betamodel name=" << name << ";");
    COUT(name << ".Sradio = " << peak_ << " " << unit << ";");
    COUT(name << formatPosition("xoff", xMid_));
    COUT(name << formatPosition("yoff", yMid_));
    COUT(name << formatPosition("thetaCore", rad1_));
    COUT(name << ".beta = 0.77;");

    break;

  case TYPE_ELBETA:

    COUT("");
    COUT("addmodel type=betamodel name=" << name << ";");
    COUT(name << ".Sradio = " << peak_ << " " << unit << ";");
    COUT(name << formatPosition("xoff", xMid_));
    COUT(name << formatPosition("yoff", yMid_));
    COUT(name << formatPosition("thetaCore", radMax_));
    COUT(name << ".axialRatio = " << radMin_/radMax_ << ";");
    COUT(name << ".rotang = " << angle_.degrees() << " deg;");
    COUT(name << ".beta = 0.77;");

    break;

  case TYPE_GAUSS:

    COUT("");
    COUT("addmodel type=gauss2d name=" << name << ";");
    COUT(name << ".Sradio = " << peak_ << " " << unit << ";");
    COUT(name << formatPosition("xoff", xMid_));
    COUT(name << formatPosition("yoff", yMid_));
    COUT(name << formatPosition("majsigma", radMax_));
    COUT(name << ".axialRatio = " << radMin_/radMax_ << ";");
    COUT(name << ".rotang = " << angle_.degrees() << " deg;");

    break;

  default:
    break;
  }
}

std::string PgModel::formatPosition(std::string varName, double pos)
{
  string xstr, ystr;
  std::ostringstream os;

  if(PgUtil::unitCallback_) {
    (*PgUtil::unitCallback_)(pos, 0.0, xstr, ystr, PgUtil::unitCallbackArgs_);
    os << "." << varName << " = " << xstr << ";";
  }

  return os.str();
}

PgModel::Type PgModel::keyToType(char key)
{
  switch (key) {
  case 'P':
    return TYPE_DELTA;
    break;
  case 'B':
    return TYPE_BETA;
    break;
  case 'E':
    return TYPE_ELBETA;
    break;
  case 'A':
    return TYPE_ARNAUD;
    break;
  case 'G':
    return TYPE_GAUSS;
    break;
  default:
    ThrowError("Unrecognized type key: " << key);
    return TYPE_NONE;
    break;
  }
}
