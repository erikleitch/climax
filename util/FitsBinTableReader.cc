#include "gcp/util/Exception.h"
#include "gcp/util/FitsBinTableReader.h"
#include "gcp/util/String.h"

#include<iostream>
#include<string>

#include<string.h>

#include "fitsio.h"

using namespace std;
using namespace gcp::util;

#define ASSIGN_VALS(type, fn)			\
  {						\
    std::vector<type> vals = fn(colNo);		\
    dvals.resize(vals.size());			\
    for(unsigned i=0; i < vals.size(); i++)	\
      dvals[i] = (double)vals[i];		\
  }

void FitsBinTableReader::initialize()
{
  tType_    = 0;
  tForm_    = 0;
  tUnit_    = 0;
  nKeyword_ = 0;
  nRow_     = 0;
  nCol_     = 0;
}

/**.......................................................................
 * Constructor.
 */
FitsBinTableReader::FitsBinTableReader(fitsfile* fptr)
{
  initialize();
  fptr_ = fptr;
}

/**.......................................................................
 * Constructor.
 */
FitsBinTableReader::FitsBinTableReader(fitsfile* fptr, std::string extension)
{
  initialize();
  fptr_ = fptr;
  getTableInfo(extension);
}

/**.......................................................................
 * Constructor.
 */
FitsBinTableReader::FitsBinTableReader(std::string file) 
{
  initialize();
  open(file);
}

/**.......................................................................
 * Constructor.
 */
FitsBinTableReader::FitsBinTableReader(std::string file, std::string extension) 
{
  initialize();
  open(file);
  getTableInfo(extension);
}

/**.......................................................................
 * Constructor.
 */
FitsBinTableReader::FitsBinTableReader(std::string file, int extnum)
{
  initialize();
  open(file);
  getTableInfoByExtnum(extnum);
}

/**.......................................................................
 * Constructor.
 */
void FitsBinTableReader::setTo(fitsfile* fptr, std::string extension)
{
  initialize();
  fptr_ = fptr;
  getTableInfo(extension);
}


void FitsBinTableReader::getTableInfoByExtnum(int extnum)
{
  int hdutype;
  if(ffmrhd(fptr_, 1, &hdutype, &status_) > 0)
    ThrowFitsError("Problem moving to next HDU");

  getTableInfo(hdutype);
}

void FitsBinTableReader::getNextTableInfo()
{
  int hdutype = moveToNextTable();
  getTableInfo(hdutype);
}

void FitsBinTableReader::getTableInfo(std::string tabspec)
{
  int hdutype = moveToTable(tabspec);
  getTableInfo(hdutype);
}

void FitsBinTableReader::getTableInfo(int hdutype)
{
  int morekeys, keynum;

  if(ffghsp(fptr_, &nKeyword_, &morekeys, &status_) > 0)
    ThrowFitsError("Unable to get header information");

  char* tmp1 = (char*)malloc(nKeyword_ * 20);
  char* tmp2 = (char*)malloc(nKeyword_ * 20);
  char* tmp3 = (char*)malloc(nKeyword_ * 20);

  tType_ = (char**)malloc(nKeyword_ * sizeof(char*));
  tForm_ = (char**)malloc(nKeyword_ * sizeof(char*));
  tUnit_ = (char**)malloc(nKeyword_ * sizeof(char*));

  for(unsigned i=0; i < nKeyword_; i++) {
    tType_[i] = tmp1 + 20*i;
    tForm_[i] = tmp2 + 20*i;
    tUnit_[i] = tmp3 + 20*i;
  }

  long pcount;

  if(hdutype == 2) {
    if(ffghbn(fptr_, nKeyword_, &nRow_, &nCol_, tType_, 
	      tForm_, tUnit_, &binName_[0], &pcount, &status_) > 0)
      ThrowFitsError("Error getting header information");
  } else {
    long rowLen;
    long tbcol;
    if(ffghtb(fptr_, nKeyword_, &rowLen, &nRow_, &nCol_, tType_, &tbcol,
	      tForm_, tUnit_, &binName_[0], &status_) > 0)
      ThrowFitsError("Error getting header information");
  }

  //------------------------------------------------------------
  // Now get the type, width and dimensions of each column
  //------------------------------------------------------------

  long naxes[FITS_BIN_MAX_DIM];
  int naxis;

  columns_.resize(nCol_);
  for(unsigned iCol=0; iCol < nCol_; iCol++) {
    columns_[iCol].name_ = tType_[iCol];
    columns_[iCol].unit_ = tUnit_[iCol];

    if(ffgtcl(fptr_, iCol+1, &columns_[iCol].typecode_, &columns_[iCol].repeat_,
	      &columns_[iCol].width_, &status_))
      ThrowFitsError("Error getting column information");

    if(ffgtdm(fptr_, iCol+1, FITS_BIN_MAX_DIM, &naxis,
	      naxes, &status_))
      ThrowFitsError("Error getting column information");

    columns_[iCol].naxes_.resize(naxis);

    for(unsigned i=0; i < naxis; i++)
      columns_[iCol].naxes_[i] = naxes[i];
  }
}

/**.......................................................................
 * Destructor.
 */
FitsBinTableReader::~FitsBinTableReader() 
{
  if(tType_) {
    if(tType_[0]) {
      free(tType_[0]);
    }
    free(tType_);
    tType_ = 0;
  }

  if(tForm_) {
    if(tForm_[0]) {
      free(tForm_[0]);
    }
    free(tForm_);
    tForm_ = 0;
  }

  if(tUnit_) {
    if(tUnit_[0]) {
      free(tUnit_[0]);
    }
    free(tUnit_);
    tUnit_ = 0;
  }


}

/**.......................................................................
 * Read data for the given column
 */
void FitsBinTableReader::getData(unsigned colNo, void* array, 
				 void* nulval, int* anynul)
{
  Column& col = columns_[colNo];
  if(fits_read_col(fptr_, col.typecode_, colNo+1, 1, 1, col.repeat_*nRow_, 
		   nulval, array, anynul, &status_))
    ThrowFitsError("Problem reading data");
}

std::vector<bool> FitsBinTableReader::getBoolData(unsigned colNo)
{
  std::vector<unsigned char> cVec = getUcharData(colNo);
  std::vector<bool> bVec;

  bVec.resize(cVec.size());

  for(unsigned i=0; i < cVec.size(); i++)
    bVec[i] = (bool)cVec[i];

  return bVec;
}
std::vector<unsigned char> FitsBinTableReader::getUcharData(unsigned colNo)
{
  std::vector<unsigned char> vals(colRepeat(colNo)*nRow());
  int anynul=0;
  unsigned char nulval=0;
  getData(colNo, (void*)&vals[0], &nulval, &anynul);

  return vals;
}
std::vector<char> FitsBinTableReader::getCharData(unsigned colNo)
{
  std::vector<char> vals(colRepeat(colNo)*nRow());
  int anynul=0;
  char nulval=0;
  getData(colNo, (void*)&vals[0], &nulval, &anynul);

  return vals;
}
std::vector<unsigned short> FitsBinTableReader::getUshortData(unsigned colNo)
{
  std::vector<unsigned short> vals(colRepeat(colNo)*nRow());
  int anynul=0;
  unsigned short nulval=0;
  getData(colNo, (void*)&vals[0], &nulval, &anynul);

  return vals;
}
std::vector<short> FitsBinTableReader::getShortData(unsigned colNo)
{
  std::vector<short> vals(colRepeat(colNo)*nRow());
  int anynul=0;
  short nulval=0;
  getData(colNo, (void*)&vals[0], &nulval, &anynul);

  return vals;
}
std::vector<unsigned int> FitsBinTableReader::getUintData(unsigned colNo)
{
  std::vector<unsigned int> vals(colRepeat(colNo)*nRow());
  int anynul=0;
  unsigned int nulval=0;
  getData(colNo, (void*)&vals[0], &nulval, &anynul);

  return vals;
}
std::vector<int> FitsBinTableReader::getIntData(unsigned colNo)
{
  std::vector<int> vals(colRepeat(colNo)*nRow());
  int anynul=0;
  int nulval=0;
  getData(colNo, (void*)&vals[0], &nulval, &anynul);

  return vals;
}
std::vector<unsigned long> FitsBinTableReader::getUlongData(unsigned colNo)
{
  std::vector<unsigned long> vals(colRepeat(colNo)*nRow());
  int anynul=0;
  unsigned long nulval=0;
  getData(colNo, (void*)&vals[0], &nulval, &anynul);

  return vals;
}
std::vector<long> FitsBinTableReader::getLongData(unsigned colNo)
{
  std::vector<long> vals(colRepeat(colNo)*nRow());
  int anynul=0;
  long nulval=0;
  getData(colNo, (void*)&vals[0], &nulval, &anynul);

  return vals;
}
std::vector<float> FitsBinTableReader::getFloatData(unsigned colNo)
{
  std::vector<float> vals(colRepeat(colNo)*nRow());
  int anynul=0;
  float nulval=0;
  getData(colNo, (void*)&vals[0], &nulval, &anynul);

  return vals;
}

std::string FitsBinTableReader::getUnit(std::string name)
{
  for(unsigned iCol=0; iCol < columns_.size(); iCol++) {
    if(columns_[iCol].name_ == name) {
      return columns_[iCol].unit_;
    }
  }
  return "?";
}

std::vector<double> FitsBinTableReader::getDoubleData(unsigned colNo)
{
  std::vector<double> vals(colRepeat(colNo)*nRow());
  int anynul=0;
  double nulval=0;
  getData(colNo, (void*)&vals[0], &nulval, &anynul);

  return vals;
}

/**.......................................................................
 * Read string data
 */
std::vector<std::string> FitsBinTableReader::getStringData(unsigned colNo)
{
  std::vector<std::string> vals;

  Column& col = columns_[colNo];

  char* arr[nRow()];
  for(unsigned i=0; i < nRow(); i++)
    arr[i] = (char*)malloc(col.repeat_+1);

  char nulval = '\0';
  int anynul;
  char** cptr = arr;

  if(ffgcvs(fptr_, colNo+1, 1, 1, nRow_, 
	    &nulval, cptr, &anynul, &status_))
    ThrowFitsError("Problem reading data");

  for(unsigned i=0; i < nRow(); i++) {
    vals.push_back(arr[i]);
    free(arr[i]);
  }

  return vals;
}

std::vector<double> FitsBinTableReader::getData(std::string colName)
{
  return getData(getColNo(colName));
}

std::vector<double> FitsBinTableReader::getData(unsigned colNo)
{
  std::vector<double> dvals;

  switch(colDataType(colNo)) {
  case DataType::DOUBLE:
    ASSIGN_VALS(double, getDoubleData);
    break;
  case DataType::FLOAT:
    ASSIGN_VALS(float, getFloatData);
    break;
  case DataType::UINT:
    ASSIGN_VALS(unsigned int, getUintData);
    break;
  case DataType::INT:
    ASSIGN_VALS(int, getIntData);
    break;
  case DataType::ULONG:
    ASSIGN_VALS(unsigned long, getUlongData);
    break;
  case DataType::LONG:
    ASSIGN_VALS(long, getLongData);
    break;
  case DataType::USHORT:
    ASSIGN_VALS(unsigned short, getUshortData);
    break;
  case DataType::SHORT:
    ASSIGN_VALS(short, getShortData);
    break;
  case DataType::UCHAR:
    ASSIGN_VALS(unsigned char, getUcharData);
    break;
  case DataType::CHAR:
    ASSIGN_VALS(char, getCharData);
    break;
  case DataType::BOOL:
    ASSIGN_VALS(bool, getBoolData);
    break;
  case DataType::STRING:
    //    ASSIGN_VALS(short, getStringData);
    break;
  default:
    ThrowError("No return function for datatype: " << colDataType(colNo));
    break;
  }

  return dvals;
}

FitsBinTableReader::Column* FitsBinTableReader::getColumn(std::string name)
{
  for(unsigned iCol=0; iCol < columns_.size(); iCol++) {
    if(columns_[iCol].name_ == name) {
      return &columns_[iCol];
    }
  }

  ThrowError("No column named: " << name);

  return 0;
}

unsigned FitsBinTableReader::getColNo(std::string name)
{
  for(unsigned iCol=0; iCol < columns_.size(); iCol++) {
    if(columns_[iCol].name_ == name) {
      return iCol;
    }
  }

  ThrowError("No column named: " << name);
  return 0;
}

std::vector<long> FitsBinTableReader::getColNaxes(std::string name)
{
  Column* col = getColumn(name);
  return col->naxes_;
}

/**.......................................................................
 * Write the contents of this object to an ostream
 */
ostream& 
gcp::util::operator<<(ostream& os, FitsBinTableReader::Column& col)
{
  os << "Name   = " << col.name_     << std::endl;
  os << "Unit   = " << col.unit_     << std::endl;
  os << "Type   = " << FitsBinTableReader::colDataType((int)col.typecode_) << std::endl;
  os << "Repeat = " << col.repeat_   << std::endl;
  os << "Width  = " << col.width_    << std::endl;
  os << "Naxis  = " << col.naxes_.size()    << std::endl;

  for(unsigned i=0; i < col.naxes_.size(); i++) {
    os << "dim[" << i << "] = " << col.naxes_[i] << std::endl;
  }

  return os;
}

DataType::Type FitsBinTableReader::colDataType(unsigned colNo)
{
  return colDataType((int)columns_[colNo].typecode_);
}

DataType::Type FitsBinTableReader::colDataType(int typecode)
{
  switch (typecode) {
  case TBYTE:
    return DataType::UCHAR;
    break;
  case TLOGICAL:
    return DataType::BOOL;
    break;
  case TSTRING:
    return DataType::STRING;
    break;
  case TUSHORT:
    return DataType::USHORT;
    break;
  case TSHORT:
    return DataType::SHORT;
    break;
  case TUINT:
    return DataType::UINT;
    break;
  case TINT:
    return DataType::INT;
    break;
  case TULONG:
    return DataType::ULONG;
    break;
  case TLONG:
    return DataType::LONG;
    break;
  case TFLOAT:
    return DataType::FLOAT;
    break;
  case TDOUBLE:
    return DataType::DOUBLE;
    break;
  case TCOMPLEX:
    return DataType::COMPLEX_FLOAT;
    break;
  case TDBLCOMPLEX:
    return DataType::COMPLEX_DOUBLE;
    break;
  case TBIT:
    return DataType::BIT;
  default:
    ThrowError("Unrecognized type: " << typecode);
    return DataType::NONE;
    break;
  }
}

void FitsBinTableReader::printColumns()
{
  for(unsigned iCol=0; iCol < nCol(); iCol++) {
    COUT("Column " << iCol << " " << std::endl << columns_[iCol]);
  }
}

std::vector<bool> FitsBinTableReader::getBoolData(std::string colName)
{
  return getBoolData(getColNo(colName));
}
std::vector<unsigned char> FitsBinTableReader::getUcharData(std::string colName)
{
  return getUcharData(getColNo(colName));
}
std::vector<char> FitsBinTableReader::getCharData(std::string colName)
{
  return getCharData(getColNo(colName));
}
std::vector<unsigned short> FitsBinTableReader::getUshortData(std::string colName)
{
  return getUshortData(getColNo(colName));
}
std::vector<short> FitsBinTableReader::getShortData(std::string colName)
{
  return getShortData(getColNo(colName));
}
std::vector<unsigned int> FitsBinTableReader::getUintData(std::string colName)
{
  return getUintData(getColNo(colName));
}
std::vector<int> FitsBinTableReader::getIntData(std::string colName)
{
  return getIntData(getColNo(colName));
}
std::vector<unsigned long> FitsBinTableReader::getUlongData(std::string colName)
{
  return getUlongData(getColNo(colName));
}
std::vector<long> FitsBinTableReader::getLongData(std::string colName)
{
  return getLongData(getColNo(colName));
}
std::vector<float> FitsBinTableReader::getFloatData(std::string colName)
{
  return getFloatData(getColNo(colName));
}
std::vector<double> FitsBinTableReader::getDoubleData(std::string colName)
{
  return getDoubleData(getColNo(colName));
}
std::vector<std::string> FitsBinTableReader::getStringData(std::string colName)
{
  return getStringData(getColNo(colName));
}

FitsReader::Axis FitsBinTableReader::getAxis(std::string valContains)
{
  FitsReader::Axis axis;
  char keyname[80], keyval[80], comment[80];
  bool reading = false;

  for(unsigned iKey=0; iKey < nKeyword_; iKey++) {
    if(ffgkyn(fptr_, iKey, keyname, keyval, comment, &status_) > 0)
      ThrowFitsError("Error getting keyword");

    String keyNameStr(keyname);
    String keyValStr(keyval);

    if(reading && keyNameStr.contains("TYP")) {
      reading = false;
      break;
    }

    if(keyValStr.contains(valContains) && keyNameStr.contains("TYP")) {
      axis.n_           = 0;
      axis.type_        = keyval;
      axis.typeComment_ = comment;
      axis.isPresent_   = true;
      reading = true;
    }

    if(reading && keyNameStr.contains("CRVL"))
      axis.refVal_ = keyValStr.toDouble();

    if(reading && keyNameStr.contains("CRPX"))
      axis.refPix_ = keyValStr.toDouble();

    if(reading && keyNameStr.contains("CDLT"))
      axis.delta_  = keyValStr.toDouble();

    if(reading && keyNameStr.contains("CUNI"))
      axis.unit_   = keyValStr.str();
  }

  if(!axis.isPresent_)
    ThrowError("No axis matching '" << valContains << "' is present");

  return axis;
}
