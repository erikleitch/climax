#define __FILEPATH__ "util/DataType.cc"

#include "gcp/util/Complex.h"
#include "gcp/util/Debug.h"
#include "gcp/util/DataType.h"
#include "gcp/util/Exception.h"
#include "gcp/util/LogStream.h"

using namespace std;

using namespace gcp::util;

void DataType::initialize() {
  type_ = NONE;
  
  ptr_     = 0;
  isArray_ = false;
  allowReassignment_ = true;
  hasValue_ = false;
}

void DataType::allowAssignmentFromDifferentType(bool allow)
{
  allowReassignment_ = allow;
}


DataType::DataType()
{
  initialize();
}

void DataType::setType(Type type)
{
  type_ = type;
}

DataType::DataType(bool b)
{
  initialize();

  data_.b = b;
  type_ = BOOL;
  hasValue_ = true;
}

DataType::DataType(unsigned char uc)
{
  initialize();

  data_.uc = uc;
  type_ = UCHAR;
  hasValue_ = true;
}

DataType::DataType(char c)
{
  initialize();

  data_.c = c;
  type_ = CHAR;
  hasValue_ = true;
}

DataType::DataType(unsigned short us)
{
  initialize();

  data_.us = us;
  type_ = USHORT;
  hasValue_ = true;
}

DataType::DataType(short s)
{
  initialize();

  data_.s = s;
  type_ = SHORT;
  hasValue_ = true;
}

DataType::DataType(unsigned int ui)
{
  initialize();

  data_.ui = ui;
  type_ = UINT;
  hasValue_ = true;
}

DataType::DataType(int i)
{
  initialize();

  data_.i = i;
  type_ = INT;
  hasValue_ = true;
}

DataType::DataType(unsigned long ul)
{
  initialize();

  data_.ul = ul;
  type_ = ULONG;
  hasValue_ = true;
}

DataType::DataType(long l)
{
  initialize();

  data_.l = l;
  type_ = LONG;
  hasValue_ = true;
}

DataType::DataType(float f)
{
  initialize();

  data_.f = f;
  type_ = FLOAT;
  hasValue_ = true;
}

DataType::DataType(double d)
{
  initialize();

  data_.d = d;
  type_ = DOUBLE;
  hasValue_ = true;
}

DataType::DataType(Complex<float> cf)
{
  initialize();

  data_.cf = *cf.data();
  type_ = COMPLEX_FLOAT;
  hasValue_ = true;
}

DataType::DataType(std::string str)
{
  initialize();

  data_.str = str;
  type_ = STRING;
  hasValue_ = true;
}

DataType::~DataType() {};

bool& DataType::hasValue()
{
  return hasValue_;
}

/**.......................................................................
 * Return the size, in bytes, of the requested type
 */
unsigned DataType::sizeOf(Type type)
{
  switch(type) {
  case BOOL:
    return sizeof(bool);
    break;
  case UCHAR:
    return sizeof(unsigned char);
    break;
  case CHAR:
    return sizeof(char);
    break;
  case USHORT:
    return sizeof(unsigned short);
    break;
  case SHORT:
    return sizeof(short);
    break;
  case UINT:
    return sizeof(unsigned int);
    break;
  case INT:
    return sizeof(int);
    break;
  case ULONG:
    return sizeof(unsigned long);
    break;
  case LONG:
    return sizeof(long);
    break;
  case FLOAT:
    return sizeof(float);
    break;
  case DOUBLE:
    return sizeof(double);
    break;
  case COMPLEX_FLOAT:
    return sizeof(Complex<float>::Data);
    break;
  default:
    return 0;
    break;
  }
}

/**.......................................................................
 * Assignment operators
 */

// Data type assignment operators

void DataType::operator=(bool b)
{
  data_.b = b;
  type_ = BOOL;
  ptr_ = 0;
  isArray_ = false;
  hasValue_ = true;
}

void DataType::operator=(unsigned char uc)
{
  data_.uc = uc;
  type_ = UCHAR;
  ptr_ = 0;
  isArray_ = false;
  hasValue_ = true;
}

void DataType::operator=(char c)
{
  data_.c = c;
  type_ = CHAR;
  ptr_ = 0;
  isArray_ = false;
  hasValue_ = true;
}

void DataType::operator=(unsigned short us)
{
  data_.us = us;
  type_ = USHORT;
  ptr_ = 0;
  isArray_ = false;
  hasValue_ = true;
}

void DataType::operator=(short s)
{
  data_.s = s;
  type_ = SHORT;
  ptr_ = 0;
  isArray_ = false;
  hasValue_ = true;
}

void DataType::operator=(unsigned int ui)
{
  data_.ui = ui;
  type_ = UINT;
  ptr_ = 0;
  isArray_ = false;
  hasValue_ = true;
}

void DataType::operator=(int i)
{
  data_.i = i;
  type_ = INT;
  ptr_ = 0;
  isArray_ = false;
  hasValue_ = true;
}

void DataType::operator=(unsigned long ul)
{
  data_.ul = ul;
  type_ = ULONG;
  ptr_ = 0;
  isArray_ = false;
  hasValue_ = true;
}

void DataType::operator=(long l)
{
  data_.l = l;
  type_ = LONG;
  ptr_ = 0;
  isArray_ = false;
  hasValue_ = true;
}

void DataType::operator=(float f)
{
  data_.f = f;
  type_ = FLOAT;
  ptr_ = 0;
  isArray_ = false;
  hasValue_ = true;
}

void DataType::operator=(double d)
{
  data_.d = d;
  type_ = DOUBLE;
  ptr_ = 0;
  isArray_ = false;
  hasValue_ = true;
}

void DataType::operator=(Complex<float> cf)
{
  data_.cf = *cf.data();
  type_ = COMPLEX_FLOAT;
  ptr_ = 0;
  isArray_ = false;
  hasValue_ = true;
}

void DataType::operator=(std::string str)
{
  data_.str = str;
  type_ = STRING;
  ptr_ = 0;
  isArray_ = false;
  hasValue_ = true;
}

// Pointer assignment operators

void DataType::operator=(bool* b)
{
  data_.b = *b;
  type_ = BOOL;
  ptr_ = (void*)b;
  isArray_ = true;
  hasValue_ = true;
}

void DataType::operator=(unsigned char* uc)
{
  data_.uc = *uc;
  type_ = UCHAR;
  ptr_ = (void*)uc;
  isArray_ = true;
  hasValue_ = true;
}

void DataType::operator=(char* c)
{
  data_.c = *c;
  type_ = CHAR;
  ptr_ = (void*)c;
  isArray_ = true;
  hasValue_ = true;
}

void DataType::operator=(unsigned short* us)
{
  data_.us = *us;
  type_ = USHORT;
  ptr_ = (void*)us;
  isArray_ = true;
  hasValue_ = true;
}

void DataType::operator=(short* s)
{
  data_.s = *s;
  type_ = SHORT;
  ptr_ = (void*)s;
  isArray_ = true;
  hasValue_ = true;
}

void DataType::operator=(unsigned int* ui)
{
  data_.ui = *ui;
  type_ = UINT;
  ptr_ = (void*)ui;
  isArray_ = true;
  hasValue_ = true;
}

void DataType::operator=(int* i)
{
  data_.i = *i;
  type_ = INT;
  ptr_ = (void*)i;
  isArray_ = true;
  hasValue_ = true;
}

void DataType::operator=(unsigned long* ul)
{
  data_.ul = *ul;
  type_ = ULONG;
  ptr_ = (void*)ul;
  isArray_ = true;
  hasValue_ = true;
}

void DataType::operator=(long* l)
{
  data_.l = *l;
  type_ = LONG;
  ptr_ = (void*)l;
  isArray_ = true;
  hasValue_ = true;
}

void DataType::operator=(float* f)
{
  data_.f = *f;
  type_ = FLOAT;
  ptr_ = (void*)f;
  isArray_ = true;
  hasValue_ = true;
}

void DataType::operator=(double* d)
{
  data_.d = *d;
  type_ = DOUBLE;
  ptr_ = (void*)d;
  isArray_ = true;
  hasValue_ = true;
}

void DataType::operator=(Complex<float>* cf)
{
  data_.cf = *(cf->data());
  type_ = COMPLEX_FLOAT;
  ptr_ = (void*)cf;
  isArray_ = true;
  hasValue_ = true;
}

void DataType::operator=(std::string* str)
{
  data_.str = *str;
  type_ = STRING;
  ptr_ = 0;
  isArray_ = false;
  hasValue_ = true;
}

// Generic assignment operators

void DataType::operator=(const DataType& dataType)
{
  data_    = dataType.data_;
  type_    = dataType.type_;
  ptr_     = dataType.ptr_;
  isArray_ = dataType.isArray_;
  hasValue_ = dataType.hasValue_;
}

void DataType::operator=(DataType& dataType)
{
  data_    = dataType.data_;
  type_    = dataType.type_;
  ptr_     = dataType.ptr_;
  isArray_ = dataType.isArray_;
  hasValue_ = dataType.hasValue_;
}

// Subtraction operators

void DataType::operator-=(DataType& dataType)
{
  checkType(dataType);

  switch (type_) {
  case UCHAR:
    data_.uc = (unsigned char)abs((int)data_.uc - (int)dataType.data_.uc);
    break;
  case CHAR:
    data_.c -= dataType.data_.c;
    break;
  case BOOL:
    data_.b = (bool)abs((int)data_.b - (int)dataType.data_.b);
    break;
  case USHORT:
    data_.us = (unsigned short)abs((int)data_.us - (int)dataType.data_.us);
    break;
  case SHORT:
    data_.s -= dataType.data_.s;
    break;
  case UINT:
    data_.ui = (unsigned int)abs((int)data_.ui - (int)dataType.data_.ui);
    break;
  case INT:
    data_.i -= dataType.data_.i;
    break;
  case ULONG:
    data_.ul = (unsigned long)labs((long)data_.ul - (long)dataType.data_.ul);
    break;
  case LONG:
    data_.l -= dataType.data_.l;
    break;
  case FLOAT:
    data_.f -= dataType.data_.f;
    break;
  case DOUBLE:
    data_.d -= dataType.data_.d;
    break;
  case COMPLEX_FLOAT:
    {
      Complex<float> cf1(data_.cf);
      Complex<float> cf2(dataType.data_.cf);
      Complex<float> diff = cf1 - cf2;
      data_.cf = *diff.data();
      data_.cf = *diff.data();
    }
    break;
  default:
    ThrowError("Unrecognized data type");
    break;
  }

  
}

bool DataType::operator==(DataType& dataType)
{
  checkType(dataType);
  switch (type_) {
  case UCHAR:
    DBPRINT(true, Debug::DEBUG15, "data_.uc = " << data_.uc
	    << "dataType.data_.uc = " << dataType.data_.uc);

    return data_.uc == dataType.data_.uc;
    break;
  case CHAR:
    return data_.c == dataType.data_.c;
    break;
  case BOOL:
    return data_.b == dataType.data_.b;
    break;
  case USHORT:
    return data_.us == dataType.data_.us;
    break;
  case SHORT:
    return data_.s == dataType.data_.s;
    break;
  case UINT:
    return data_.ui == dataType.data_.ui;
    break;
  case INT:
    return data_.i == dataType.data_.i;
    break;
  case ULONG:
    return data_.ul == dataType.data_.ul;
    break;
  case LONG:
    return data_.l == dataType.data_.l;
    break;
  case FLOAT:
    return data_.f == dataType.data_.f;
    break;
  case DOUBLE:
    return data_.d == dataType.data_.d;
    break;
  case COMPLEX_FLOAT:
    {
      Complex<float> cf1(data_.cf);
      Complex<float> cf2(dataType.data_.cf);
      return cf1 == cf2;
    }
    break;
  default:
    return false;
    break;
  }
}

bool DataType::operator>(DataType& dataType)
{
  checkType(dataType);
  switch (type_) {
  case UCHAR:
    return data_.uc > dataType.data_.uc;
    break;
  case CHAR:
    return data_.c > dataType.data_.c;
    break;
  case BOOL:
    return data_.b > dataType.data_.b;
    break;
  case USHORT:
    return data_.us > dataType.data_.us;
    break;
  case SHORT:
    return data_.s > dataType.data_.s;
    break;
  case UINT:
    return data_.ui > dataType.data_.ui;
    break;
  case INT:
    return data_.i > dataType.data_.i;
    break;
  case ULONG:
    return data_.ul > dataType.data_.ul;
    break;
  case LONG:
    return data_.l > dataType.data_.l;
    break;
  case FLOAT:
    return data_.f > dataType.data_.f;
    break;
  case DOUBLE:
    return data_.d > dataType.data_.d;
    break;
  case COMPLEX_FLOAT:
    {
      Complex<float> cf1(data_.cf);
      Complex<float> cf2(dataType.data_.cf);
      return cf1 > cf2;
    }
    break;
  default:
    return false;
    break;
  }
}

bool DataType::operator>=(DataType& dataType)
{
  checkType(dataType);
  switch (type_) {
  case UCHAR:
    return data_.uc >= dataType.data_.uc;
    break;
  case CHAR:
    return data_.c >= dataType.data_.c;
    break;
  case BOOL:
    return data_.b >= dataType.data_.b;
    break;
  case USHORT:
    return data_.us >= dataType.data_.us;
    break;
  case SHORT:
    return data_.s >= dataType.data_.s;
    break;
  case UINT:
    return data_.ui >= dataType.data_.ui;
    break;
  case INT:
    return data_.i >= dataType.data_.i;
    break;
  case ULONG:
    return data_.ul >= dataType.data_.ul;
    break;
  case LONG:
    return data_.l >= dataType.data_.l;
    break;
  case FLOAT:
    return data_.f >= dataType.data_.f;
    break;
  case DOUBLE:
    return data_.d >= dataType.data_.d;
    break;
  case COMPLEX_FLOAT:
    {
      Complex<float> cf1(data_.cf);
      Complex<float> cf2(dataType.data_.cf);
      return cf1 >= cf2;
    }
    break;
  default:
    return false;
    break;
  }
}

bool DataType::operator<(DataType& dataType)
{
  checkType(dataType);
  switch (type_) {
  case UCHAR:
    return data_.uc < dataType.data_.uc;
    break;
  case CHAR:
    return data_.c < dataType.data_.c;
    break;
  case BOOL:
    return data_.b < dataType.data_.b;
    break;
  case USHORT:
    return data_.us < dataType.data_.us;
    break;
  case SHORT:
    return data_.s < dataType.data_.s;
    break;
  case UINT:
    return data_.ui < dataType.data_.ui;
    break;
  case INT:
    return data_.i < dataType.data_.i;
    break;
  case ULONG:
    return data_.ul < dataType.data_.ul;
    break;
  case LONG:
    return data_.l < dataType.data_.l;
    break;
  case FLOAT:
    return data_.f < dataType.data_.f;
    break;
  case DOUBLE:
    return data_.d < dataType.data_.d;
    break;
  case COMPLEX_FLOAT:
    {
      Complex<float> cf1(data_.cf);
      Complex<float> cf2(dataType.data_.cf);
      return cf1 < cf2;
    }
    break;
  default:
    return false;
    break;
  }
}

bool DataType::operator<=(DataType& dataType)
{
  checkType(dataType);
  switch (type_) {
  case UCHAR:
    return data_.uc <= dataType.data_.uc;
    break;
  case CHAR:
    return data_.c <= dataType.data_.c;
    break;
  case BOOL:
    return data_.b <= dataType.data_.b;
    break;
  case USHORT:
    return data_.us <= dataType.data_.us;
    break;
  case SHORT:
    return data_.s <= dataType.data_.s;
    break;
  case UINT:
    return data_.ui <= dataType.data_.ui;
    break;
  case INT:
    return data_.i <= dataType.data_.i;
    break;
  case ULONG:
    return data_.ul <= dataType.data_.ul;
    break;
  case LONG:
    return data_.l <= dataType.data_.l;
    break;
  case FLOAT:
    return data_.f <= dataType.data_.f;
    break;
  case DOUBLE:
    return data_.d <= dataType.data_.d;
    break;
  case COMPLEX_FLOAT:
    {
      Complex<float> cf1(data_.cf);
      Complex<float> cf2(dataType.data_.cf);
      return cf1 <= cf2;
    }
    break;
  default:
    return false;
    break;
  }
}

ostream& 
gcp::util::operator<<(ostream& os, DataType::Type type)
{
  switch (type) {
  case DataType::UCHAR:
    os << "UCHAR";
    break;
  case DataType::CHAR:
    os << "CHAR";
    break;
  case DataType::BOOL:
    os << "BOOL";
    break;
  case DataType::USHORT:
    os << "USHORT";
    break;
  case DataType::SHORT:
    os << "SHORT";
    break;
  case DataType::UINT:
    os << "UINT";
    break;
  case DataType::INT:
    os << "INT";
    break;
  case DataType::ULONG:
    os << "ULONG";
    break;
  case DataType::LONG:
    os << "LONG";
    break;
  case DataType::FLOAT:
    os << "FLOAT";
    break;
  case DataType::DOUBLE:
    os << "DOUBLE";
    break;
  case DataType::COMPLEX_FLOAT:
    os << "COMPLEX_FLOAT";
    break;
  case DataType::STRING:
    os << "STRING";
    break;
  case DataType::BIT:
    os << "BIT";
    break;
  default:
    os << "UNKNOWN(" << static_cast<unsigned int>(type) << ")";
    ThrowError("Unknown data type");
    break;
  }
  return os;
}

 
void DataType::checkType(DataType& dataType)
{
  if(dataType.type_ != type_)
    ThrowError("Data types don't match");
}
void* DataType::data()
{
  switch (type_) {
  case UCHAR:
    return (void*)&data_.uc;
    break;
  case CHAR:
    return (void*)&data_.c;
    break;
  case BOOL:
    return (void*)&data_.b;
    break;
  case USHORT:
    return (void*)&data_.us;
    break;
  case SHORT:
    return (void*)&data_.s;
    break;
  case UINT:
    return (void*)&data_.ui;
    break;
 case INT:
    return (void*)&data_.i;
    break;
  case ULONG:
    return (void*)&data_.ul;
    break;
  case LONG:
    return (void*)&data_.l;
    break;
  case FLOAT:
    return (void*)&data_.f;
    break;
  case DOUBLE:
    return (void*)&data_.d;
    break;
  case COMPLEX_FLOAT:
    return (void*)&data_.cf;
    break;
  default:
    ThrowError("No data type specified");
    return 0;
    break;
  }
}

/**.......................................................................
 * If this data type object is storing an array, increment to the next
 * element in the array, else do nothing.
 */
void DataType::operator+=(unsigned n)
{
  for(unsigned i=0; i < n; i++) {
    ++(*this);
  }
}

/**.......................................................................                                          
 * If this data type object is storing an array, increment to the next                                              
 * element in the array, else do nothing.                                                                           
 */
void DataType::operator++()
{
  // Do nothing if we are not harbouring an array                                                                   

  if(!isArray_)
    return;

  // Else re-cast from the void ptr and use our native assignmemnt                                                  
  // operator to store the next value in the array                                                                  

  else {
    switch (type_) {
    case UCHAR:
      *this = ((unsigned char*)(ptr_))+1;
      break;
    case CHAR:
      *this = ((char*)(ptr_))+1;
      break;
    case BOOL:
      *this = ((bool*)(ptr_))+1;
      break;
    case USHORT:
      *this = ((unsigned short*)(ptr_))+1;
      break;
    case SHORT:
     *this = ((short*)(ptr_))+1;
      break;
    case UINT:
      *this = ((unsigned int*)(ptr_))+1;
      break;
    case INT:
      *this = ((int*)(ptr_))+1;
      break;
    case ULONG:
      *this = ((unsigned long*)(ptr_))+1;
      break;
    case LONG:
      *this = ((long*)(ptr_))+1;
      break;
    case FLOAT:
      *this = ((float*)(ptr_))+1;
      break;
    case DOUBLE:
      *this = ((double*)(ptr_))+1;
      break;
    case COMPLEX_FLOAT:
      *this = ((Complex<float>*)(ptr_))+1;
      break;
    default:
      ThrowError("No data type specified");
      break;
    }
  }
}

std::ostream&
gcp::util::operator<<(std::ostream& os, DataType& dataType)
{
  switch (dataType.type_) {
  case DataType::UCHAR:
    os << (int)dataType.data_.uc;
    break;
  case DataType::CHAR:
    os << dataType.data_.c;
    break;
  case DataType::BOOL:
    os << dataType.data_.b;
    break;
  case DataType::USHORT:
    os << dataType.data_.us;
    break;
  case DataType::SHORT:
    os << dataType.data_.s;
    break;
  case DataType::UINT:
    os << dataType.data_.ui;
    break;
  case DataType::INT:
    os << dataType.data_.i;
    break;
  case DataType::ULONG:
    os << dataType.data_.ul;
    break;
  case DataType::LONG:
    os << dataType.data_.l;
    break;
  case DataType::FLOAT:
    os << dataType.data_.f;
    break;
  case DataType::DOUBLE:
    os << dataType.data_.d;
    break;
  case DataType::STRING:
    os << dataType.data_.str;
    break;
  case DataType::COMPLEX_FLOAT:
    {
      Complex<float> cf(dataType.data_.cf);
      os << cf;
    }
    break;
  default:
    os << "Unrecognized type";
    break;
  }

  return os;
}

double DataType::getValAsDouble()
{
  switch (type_) {
  case UCHAR:
    return (double)data_.uc;
    break;
  case CHAR:
    return (double)data_.c;
    break;
  case BOOL:
    return (double)data_.b;
    break;
  case USHORT:
    return (double)data_.us;
    break;
  case SHORT:
    return (double)data_.s;
    break;
  case UINT:
    return (double)data_.ui;
    break;
  case INT:
    return (double)data_.i;
    break;
  case ULONG:
    return (double)data_.ul;
    break;
  case LONG:
    return (double)data_.l;
    break;
  case FLOAT:
    return (double)data_.f;
    break;
  case DOUBLE:
    return (double)data_.d;
    break;
  default:
    ThrowError("Invalid Data type conversion");
    return 0.0;
    break;
  }
}

void DataType::convertToAbs()
{
  switch (type_) {
  case UCHAR:
    break;
  case CHAR:
    data_.c = (char)abs((int)data_.c);
    break;
  case BOOL:
    break;
  case USHORT:
    break;
  case SHORT:
    data_.s = (short)abs((int)data_.s);
    break;
  case UINT:
    break;
  case INT:
    data_.i = (int)abs((int)data_.i);
    break;
  case ULONG:
    break;
  case LONG:
    data_.l = (long)labs(data_.l);
    break;
  case FLOAT:
    data_.f = fabsf(data_.f);
    break;
  case DOUBLE:
    data_.d = fabs(data_.d);
    break;
  case COMPLEX_FLOAT:
    break;
  default:
    break;
  }
}

DataType::Type DataType::typeOf(bool* obj)
{
  return DataType::BOOL;
}
DataType::Type DataType::typeOf(unsigned char* obj)
{
  return DataType::UCHAR;
}
DataType::Type DataType::typeOf(char* obj)
{
  return DataType::CHAR;
}
DataType::Type DataType::typeOf(unsigned short* obj)
{
  return DataType::USHORT;
}
DataType::Type DataType::typeOf(short* obj)
{
  return DataType::SHORT;
}
DataType::Type DataType::typeOf(unsigned int* obj)
{
  return DataType::UINT;
}
DataType::Type DataType::typeOf(int* obj)
{
  return DataType::INT;
}
DataType::Type DataType::typeOf(unsigned long* obj)
{
  return DataType::ULONG;
}
DataType::Type DataType::typeOf(long* obj)
{
  return DataType::LONG;
}
DataType::Type DataType::typeOf(float* obj)
{
  return DataType::FLOAT;
}
DataType::Type DataType::typeOf(double* obj)
{
  return DataType::DOUBLE;
}
DataType::Type DataType::typeOf(gcp::util::Complex<float>* obj)
{
  return DataType::COMPLEX_FLOAT;
}
DataType::Type DataType::typeOf(gcp::util::Complex<float>::Data* obj)
{
  return DataType::COMPLEX_FLOAT;
}

DataType::Type DataType::typeOf(std::string* str)
{
  return DataType::STRING;
}

bool DataType::getBoolVal()
{
  if(type_ != DataType::BOOL) {
    ThrowColorError("This datatype is not a Bool type", "red");
  }

  return data_.b;
}

unsigned char DataType::getUcharVal()
{
  if(type_ != DataType::UCHAR) {
    ThrowColorError("This datatype is not a Uchar type", "red");
  }

  return data_.uc;
}

char DataType::getCharVal()
{
  if(type_ != DataType::CHAR) {
    ThrowColorError("This datatype is not a Char type", "red");
  }

  return data_.c;
}

unsigned short DataType::getUshortVal()
{
  if(type_ != DataType::USHORT) {
    ThrowColorError("This datatype is not a Ushort type", "red");
  }

  return data_.us;
}

short DataType::getShortVal()
{
  if(type_ != DataType::SHORT) {
    ThrowColorError("This datatype is not a Short type", "red");
  }

  return data_.s;
}

unsigned int DataType::getUintVal()
{
  if(type_ != DataType::UINT) {
    ThrowColorError("This datatype is not a Uint type", "red");
  }

  return data_.ui;
}

int DataType::getIntVal()
{
  if(type_ != DataType::INT) {
    ThrowColorError("This datatype is not a Int type", "red");
  }

  return data_.i;
}

unsigned long DataType::getUlongVal()
{
  if(type_ != DataType::ULONG) {
    ThrowColorError("This datatype is not a Ulong type", "red");
  }

  return data_.ul;
}

long DataType::getLongVal()
{
  if(type_ != DataType::LONG) {
    ThrowColorError("This datatype is not a Long type", "red");
  }

  return data_.l;
}

float DataType::getFloatVal()
{
  if(type_ != DataType::FLOAT) {
    ThrowColorError("This datatype is not a Float type", "red");
  }

  return data_.f;
}

double DataType::getDoubleVal()
{
  if(type_ != DataType::DOUBLE) {
    ThrowColorError("This datatype is not a Double type", "red");
  }

  return data_.d;
}

std::string DataType::getStringVal()
{
  if(type_ != DataType::STRING) {
    ThrowColorError("This datatype is not a String type", "red");
  }

  return data_.str;
}
