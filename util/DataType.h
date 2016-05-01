// $Id: DataType.h,v 1.3 2012/05/09 21:17:48 eml Exp $

#ifndef GCP_UTIL_DATATYPE_H
#define GCP_UTIL_DATATYPE_H

/**
 * @file DataType.h
 * 
 * Tagged: Tue Jun 22 22:32:16 UTC 2004
 * 
 * @version: $Revision: 1.3 $, $Date: 2012/05/09 21:17:48 $
 * 
 * @author 
 */
#include "gcp/util/Complex.h"

#include <ostream>

namespace gcp {
  namespace util {
    
    /**
     * Enumerate data types
     */
    class DataType {
    public:
      
      enum Type {
	NONE          = 0x0,
	UNKNOWN       = 0x0,
	UCHAR         = 0x1,
	CHAR          = 0x2,
	BOOL          = 0x4,
	USHORT        = 0x8,
	SHORT         = 0x10,
	UINT          = 0x20,
	INT           = 0x40,
	ULONG         = 0x80,
	LONG          = 0x100,
	FLOAT         = 0x200,
	DOUBLE        = 0x400,
	COMPLEX_FLOAT = 0x800,

        // Not handled yet in this object, but included for
        // compatibility with FITS format types

        COMPLEX_DOUBLE= 0x1000,
        STRING        = 0x2000,
        DATE          = 0x4000,
        BIT           = 0x8000
      };

      
      DataType();
      DataType(bool b);
      DataType(unsigned char uc);
      DataType(char ch);
      DataType(unsigned short us);
      DataType(short s);
      DataType(unsigned int ui);
      DataType(int i);
      DataType(unsigned long ul);
      DataType(long l);
      DataType(float f);
      DataType(double d);
      DataType(gcp::util::Complex<float> cf);
      DataType(std::string str);

      void initialize();

      // If true, then assignment operators will allow reassignment of
      // this object to the assigned type.  If false, an attempt to
      // assign from a different type will throw.

      void allowAssignmentFromDifferentType(bool allow);

      // Set the type of this object

      void setType(Type type);

      bool& hasValue();

      virtual ~DataType();
      
      /**
       * Return the size, in bytes, of the requested type
       */
      static unsigned sizeOf(Type type);
      
      double getValAsDouble();

      /**
       * Accessor methods
       */
      bool getBoolVal();
      unsigned char getUcharVal();
      char getCharVal();
      unsigned short getUshortVal();
      short getShortVal();
      unsigned int getUintVal();
      int getIntVal();
      unsigned long getUlongVal();
      long getLongVal();
      float getFloatVal();
      double getDoubleVal();
      std::string getStringVal();

      /**
       * Assignment operators
       */
      void operator=(bool b);
      void operator=(unsigned char uc);
      void operator=(char ch);
      void operator=(unsigned short us);
      void operator=(short s);
      void operator=(unsigned int ui);
      void operator=(int i);
      void operator=(unsigned long ul);
      void operator=(long l);
      void operator=(float f);
      void operator=(double d);
      void operator=(gcp::util::Complex<float> cf);
      void operator=(std::string str);

      /**
       * Assignment operators for pointers
       */
      void operator=(bool* b);
      void operator=(unsigned char* uc);
      void operator=(char* ch);
      void operator=(unsigned short* us);
      void operator=(short* s);
      void operator=(unsigned int* ui);
      void operator=(int* i);
      void operator=(unsigned long* ul);
      void operator=(long* l);
      void operator=(float* f);
      void operator=(double* d);
      void operator=(gcp::util::Complex<float>* cf);
      void operator=(std::string* str);

      void operator=(DataType& dataType);
      void operator=(const DataType& dataType);
      
      void operator-=(DataType& dataType);

      bool operator==(DataType& dataType);
      bool operator>(DataType& dataType);
      bool operator>=(DataType& dataType);
      bool operator<(DataType& dataType);
      bool operator<=(DataType& dataType);

      void operator++();
      void operator+=(unsigned n);

      void convertToAbs();

      friend std::ostream& operator<<(std::ostream& os, DataType& dataType);
      friend std::ostream& operator<<(std::ostream& os, DataType::Type type);

      /**                                                                                                           
       * Return a void ptr to the data for this data type
       */
      void* data();

      // The actual data in this DataType will be stored as a union

      struct {	
	bool b;
	unsigned char uc;
	char c;
	unsigned short us;
	short s;
	unsigned int ui;
	int i;
	unsigned long ul;
	long l;
	float f;
	double d;
	gcp::util::Complex<float>::Data cf;
	std::string str;
      } data_;
      
      Type type_;
      bool allowReassignment_;
      bool hasValue_;

      void checkType(DataType& dataType);

      // True if we are indexing an array of data

      bool isArray_;

      // A void ptr to the array

      void* ptr_;

      // Overloaded utility methods

      static Type typeOf(bool* obj);
      static Type typeOf(unsigned char* obj);
      static Type typeOf(char* obj);
      static Type typeOf(unsigned short* obj);
      static Type typeOf(short* obj);
      static Type typeOf(unsigned int* obj);
      static Type typeOf(int* obj);
      static Type typeOf(unsigned long* obj);
      static Type typeOf(long* obj);
      static Type typeOf(float* obj);
      static Type typeOf(double* obj);
      static Type typeOf(gcp::util::Complex<float>* obj);
      static Type typeOf(gcp::util::Complex<float>::Data* obj);
      static Type typeOf(std::string* str);

    }; // End class DataType
    
  } // End namespace util
} // End namespace gcp




#endif // End #ifndef GCP_UTIL_DATATYPE_H
