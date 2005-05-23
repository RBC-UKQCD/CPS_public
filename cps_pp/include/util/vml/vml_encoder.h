#ifndef _PAB_VML_ENCODER_H_
#define _PAB_VML_ENCODER_H_
/* Peter Boyle 2004 Virtual markup language*/
#include <util/vml/vml.h>
#include <stdio.h>
/*
 * Try to abstract out references, inheritance, arrays, structs etc..
 */

CPS_START_NAMESPACE

class GenericEncoder
{
 public:
  /*
   * Interface to represent the basic types in XDR routines
   */
  virtual bool_t Void(void)=0;
  virtual bool_t Short ( VML *vmls, char *name, short &val ) =0;
  virtual bool_t UnsignedShort ( VML *vmls, char *name, unsigned short &val ) =0;
  virtual bool_t Int ( VML *vmls, char *name, int &val ) =0;
  virtual bool_t UnsignedInt ( VML *vmls, char *name, unsigned int &val ) =0;
  virtual bool_t Long ( VML *vmls, char *name, long &val ) =0;
  virtual bool_t UnsignedLong ( VML *vmls, char *name, unsigned long &val ) =0;
  virtual bool_t LongLong ( VML *vmls, char *name, long long &val ) =0;
  virtual bool_t UnsignedLongLong ( VML *vmls, char *name, unsigned long long &val ) =0;
  virtual bool_t Int8 ( VML *vmls, char *name, int8_t &val ) =0;
  virtual bool_t UnsignedInt8 ( VML *vmls, char *name, uint8_t &val ) =0;
  virtual bool_t Int16 ( VML *vmls, char *name, int16_t &val ) =0;
  virtual bool_t UnsignedInt16 ( VML *vmls, char *name, uint16_t &val ) =0;
  virtual bool_t Int32 ( VML *vmls, char *name, int32_t &val ) =0;
  virtual bool_t UnsignedInt32 ( VML *vmls, char *name, uint32_t &val ) =0;
  virtual bool_t Int64 ( VML *vmls, char *name, int64_t &val ) =0;
  virtual bool_t UnsignedInt64 ( VML *vmls, char *name, uint64_t &val ) =0;
  virtual bool_t Char ( VML *vmls, char *name, char &val ) =0;
  virtual bool_t UnsignedChar ( VML *vmls, char *name, unsigned char &val ) =0;
  virtual bool_t Bool ( VML *vmls, char *name, bool_t &val ) =0;
  virtual bool_t Enum ( VML *vmls, char *ename, char *name, char *&value ) =0;
  virtual bool_t Double(VML *vmls, char *name, double &d)=0;
  virtual bool_t Float (VML *vmls, char *name, float  &f)=0;

  virtual bool_t Array     ( VML *vmls, char *type, char *name, 
			     char * &vals, int &nvals,
			     int &sizeofone, vmlproc_t do_one, int DoAlloc = 1 )=0;

  virtual bool_t Bytes( VML *vmls, char *name, char *&vals, int &length ) =0;

  virtual bool_t String (VML *vmls, char *name, char *&str )=0;

  virtual bool_t Reference ( VML *vmls, char *type, char *name, 
			     char *&ref, vmlproc_t do_ref, int sizeofone  )=0;
#define NOSUP printf("Not implemented"); exit(-1);
  virtual bool_t Inherits  ( VML *vmls, char *type, char *name, 
			     char *&ref, vmlproc_t do_ref  ){NOSUP};

  
  virtual bool_t StructBegin( VML *vmls, char *type, char *instance )=0;
  virtual bool_t StructEnd  ( VML *vmls, char *type, char *instance )=0;
  virtual bool_t ClassBegin ( VML *vmls, char *name, char *instance ){NOSUP};
  virtual bool_t ClassEnd   ( VML *vmls, char *name, char *instance ){NOSUP};

};

class TextEncoder : public GenericEncoder 
{
  virtual bool_t Void(void);
  virtual bool_t Short ( VML *vmls, char *name, short &val ) ;
  virtual bool_t UnsignedShort ( VML *vmls, char *name, unsigned short &val ) ;
  virtual bool_t Int ( VML *vmls, char *name, int &val ) ;
  virtual bool_t UnsignedInt ( VML *vmls, char *name, unsigned int &val ) ;
  virtual bool_t Long ( VML *vmls, char *name, long &val ) ;
  virtual bool_t UnsignedLong ( VML *vmls, char *name, unsigned long &val ) ;
  virtual bool_t LongLong ( VML *vmls, char *name, long long &val ) ;
  virtual bool_t UnsignedLongLong ( VML *vmls, char *name, unsigned long long &val ) ;
  virtual bool_t Int8 ( VML *vmls, char *name, int8_t &val ) ;
  virtual bool_t UnsignedInt8 ( VML *vmls, char *name, uint8_t &val ) ;
  virtual bool_t Int16 ( VML *vmls, char *name, int16_t &val ) ;
  virtual bool_t UnsignedInt16 ( VML *vmls, char *name, uint16_t &val ) ;
  virtual bool_t Int32 ( VML *vmls, char *name, int32_t &val ) ;
  virtual bool_t UnsignedInt32 ( VML *vmls, char *name, uint32_t &val ) ;
  virtual bool_t Int64 ( VML *vmls, char *name, int64_t &val ) ;
  virtual bool_t UnsignedInt64 ( VML *vmls, char *name, uint64_t &val ) ;
  virtual bool_t Char ( VML *vmls, char *name, char &val ) ;
  virtual bool_t UnsignedChar ( VML *vmls, char *name, unsigned char &val ) ;
  virtual bool_t Enum ( VML *vmls, char *ename, char *name, char *&value);
  virtual bool_t Bool ( VML *vmls, char *name, bool_t &val ) ;
  virtual bool_t Double(VML *vmls, char *name, double &d);
  virtual bool_t Float (VML *vmls, char *name, float  &f);



  virtual bool_t Array     ( VML *vmls, char *type, char *name, 
			     char * &vals, int &nvals,
			     int &sizeofone, vmlproc_t do_one ,
			     int DoAlloc = 1);

  virtual bool_t Bytes( VML *vmls, char *name, char *&vals, int &length ) ;

  virtual bool_t String (VML *vmls, char *name, char *&str );

  virtual bool_t Reference ( VML *vmls, char *type, char *name, 
			     char *&ref, vmlproc_t do_ref , int sizeofone );
#define NOSUP printf("Not implemented"); exit(-1);
  virtual bool_t Inherits  ( VML *vmls, char *type, char *name, 
		     char *&ref, vmlproc_t do_ref  ){NOSUP};

  
  virtual bool_t StructBegin( VML *vmls, char *type, char *instance );
  virtual bool_t StructEnd  ( VML *vmls, char *type, char *instance );
  virtual bool_t ClassBegin ( VML *vmls, char *name, char *instance );
  virtual bool_t ClassEnd   ( VML *vmls, char *name, char *instance);

};

class XmlEncoder : public GenericEncoder
{
  virtual bool_t Void(void);
  virtual bool_t Short ( VML *vmls, char *name, short &val ) ;
  virtual bool_t UnsignedShort ( VML *vmls, char *name, unsigned short &val ) ;
  virtual bool_t Int ( VML *vmls, char *name, int &val ) ;
  virtual bool_t UnsignedInt ( VML *vmls, char *name, unsigned int &val ) ;
  virtual bool_t Long ( VML *vmls, char *name, long &val ) ;
  virtual bool_t UnsignedLong ( VML *vmls, char *name, unsigned long &val ) ;
  virtual bool_t LongLong ( VML *vmls, char *name, long long &val ) ;
  virtual bool_t UnsignedLongLong ( VML *vmls, char *name, unsigned long long &val ) ;
  virtual bool_t Int8 ( VML *vmls, char *name, int8_t &val ) ;
  virtual bool_t UnsignedInt8 ( VML *vmls, char *name, uint8_t &val ) ;
  virtual bool_t Int16 ( VML *vmls, char *name, int16_t &val ) ;
  virtual bool_t UnsignedInt16 ( VML *vmls, char *name, uint16_t &val ) ;
  virtual bool_t Int32 ( VML *vmls, char *name, int32_t &val ) ;
  virtual bool_t UnsignedInt32 ( VML *vmls, char *name, uint32_t &val ) ;
  virtual bool_t Int64 ( VML *vmls, char *name, int64_t &val ) ;
  virtual bool_t UnsignedInt64 ( VML *vmls, char *name, uint64_t &val ) ;
  virtual bool_t Char ( VML *vmls, char *name, char &val ) ;
  virtual bool_t UnsignedChar ( VML *vmls, char *name, unsigned char &val ) ;
  virtual bool_t Bool ( VML *vmls, char *name, bool_t &val ) ;
  virtual bool_t Enum ( VML *vmls, char *ename, char *name, char *&value);
  virtual bool_t Double(VML *vmls, char *name, double &d);
  virtual bool_t Float (VML *vmls, char *name, float  &f);

  virtual bool_t Array     ( VML *vmls, char *type, char *name, 
			     char * &vals, int &nvals,
			     int &sizeofone, vmlproc_t do_one , 
			     int DoAlloc = 1);

  virtual bool_t Bytes( VML *vmls, char *name, char *&vals, int &length ) ;

  virtual bool_t String (VML *vmls, char *name, char *&str );

  virtual bool_t Reference ( VML *vmls, char *type, char *name, 
			     char *&ref, vmlproc_t do_ref , int sizeofone );
#define NOSUP printf("Not implemented"); exit(-1);
  virtual bool_t Inherits  ( VML *vmls, char *type, char *name, 
		     char *&ref, vmlproc_t do_ref  ){NOSUP};

  virtual bool_t StructBegin( VML *vmls, char *name, char *instance );
  virtual bool_t StructEnd  ( VML *vmls, char *name, char *instance );
  virtual bool_t ClassBegin ( VML *vmls, char *name, char *instance );
  virtual bool_t ClassEnd   ( VML *vmls, char *name, char *instance );

};

class XdrEncoder : public GenericEncoder
{
  virtual bool_t Void(void);
  virtual bool_t Short ( VML *vmls, char *name, short &val ) ;
  virtual bool_t UnsignedShort ( VML *vmls, char *name, unsigned short &val ) ;
  virtual bool_t Int ( VML *vmls, char *name, int &val ) ;
  virtual bool_t UnsignedInt ( VML *vmls, char *name, unsigned int &val ) ;
  virtual bool_t Long ( VML *vmls, char *name, long &val ) ;
  virtual bool_t UnsignedLong ( VML *vmls, char *name, unsigned long &val ) ;
  virtual bool_t LongLong ( VML *vmls, char *name, long long &val ) ;
  virtual bool_t UnsignedLongLong ( VML *vmls, char *name, unsigned long long &val ) ;
  virtual bool_t Int8 ( VML *vmls, char *name, int8_t &val ) ;
  virtual bool_t UnsignedInt8 ( VML *vmls, char *name, uint8_t &val ) ;
  virtual bool_t Int16 ( VML *vmls, char *name, int16_t &val ) ;
  virtual bool_t UnsignedInt16 ( VML *vmls, char *name, uint16_t &val ) ;
  virtual bool_t Int32 ( VML *vmls, char *name, int32_t &val ) ;
  virtual bool_t UnsignedInt32 ( VML *vmls, char *name, uint32_t &val ) ;
  virtual bool_t Int64 ( VML *vmls, char *name, int64_t &val ) ;
  virtual bool_t UnsignedInt64 ( VML *vmls, char *name, uint64_t &val ) ;
  virtual bool_t Char ( VML *vmls, char *name, char &val ) ;
  virtual bool_t UnsignedChar ( VML *vmls, char *name, unsigned char &val ) ;
  virtual bool_t Bool ( VML *vmls, char *name, bool_t &val ) ;
  virtual bool_t Double(VML *vmls, char *name, double &d);
  virtual bool_t Float (VML *vmls, char *name, float  &f);

  virtual bool_t Array     ( VML *vmls, char *type, char *name, 
			     char * &vals, int &nvals,
			     int &sizeofone, vmlproc_t do_one , 
			     int DoAlloc = 1);

  virtual bool_t Bytes( VML *vmls, char *name, char *&vals, int &length ) ;

  virtual bool_t String (VML *vmls, char *name, char *&str );

  virtual bool_t Reference ( VML *vmls, char *type, char *name, 
			     char *&ref, vmlproc_t do_ref , int sizeofone );

  virtual bool_t Inherits  ( VML *vmls, char *type, char *name, 
		     char *&ref, vmlproc_t do_ref  ){NOSUP};

  virtual bool_t StructBegin( VML *vmls, char *name, char *instance ){NOSUP};
  virtual bool_t StructEnd  ( VML *vmls, char *name, char *instance ){NOSUP};
  virtual bool_t ClassBegin ( VML *vmls, char *name, char *instance ){NOSUP};
  virtual bool_t ClassEnd   ( VML *vmls, char *name, char *instance ){NOSUP};

};
CPS_END_NAMESPACE
#endif
