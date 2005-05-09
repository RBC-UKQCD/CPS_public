#include <config.h>
#include <util/vml/vml_encoder.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
#include <string.h>

#if 0
  /* Peter Boyle 2005
   * Interface to represent the basic types in XDR routines
   */

bool_t XdrEncoder::Void(void)
{
  return true;
}


bool_t XdrEncoder::Enum ( VML *vmls,char *ename, char *name, char *&value)
{
  return xdr_enum(vmls->xmls,(enum_t *)value);
}
bool_t XdrEncoder::Short ( VML *vmls, char *name, short &val ) 
{
  return xdr_short(vmls->xmls,&val);
}
bool_t XdrEncoder::UnsignedShort (VML *vmls, char *name, unsigned short &val) 
{
  return xdr_u_short(vmls->xmls,&val);
}
bool_t XdrEncoder::Int ( VML *vmls, char *name, int &val ) 
{
  return xdr_int(vmls->xmls,&val);
}
bool_t XdrEncoder::UnsignedInt ( VML *vmls, char *name, unsigned int &val ) 
{
  return xdr_u_int(vmls->xmls,&val);
}
bool_t XdrEncoder::Long ( VML *vmls, char *name, long &val ) 
{
  return xdr_long(vmls->xmls,&val);
}
bool_t XdrEncoder::UnsignedLong ( VML *vmls, char *name, unsigned long &val )
{
  return xdr_u_long(vmls->xmls,&val);
}
bool_t XdrEncoder::LongLong ( VML *vmls, char *name, long long &val ) 
{
  return xdr_longlong_t(vmls->xmls,&val);
}
bool_t XdrEncoder::UnsignedLongLong ( VML *vmls, char *name, 
				       unsigned long long &val ) 
{
  return xdr_u_longlong_t(vmls->xmls,&val);
}
bool_t XdrEncoder::Int8 ( VML *vmls, char *name, int8_t &val ) 
{
  return xdr_int8_t(vmls->xmls,&val);
}
bool_t XdrEncoder::UnsignedInt8 ( VML *vmls, char *name, uint8_t &val ) 
{
  return xdr_u_int8_t(vmls->xmls,&val);
}
bool_t XdrEncoder::Int16 ( VML *vmls, char *name, int16_t &val ) 
{
  return xdr_int16_t(vmls->xmls,&val);
}
bool_t XdrEncoder::UnsignedInt16 ( VML *vmls, char *name, uint16_t &val ) 
{
  return xdr_u_int16_t(vmls->xmls,&val);
}
bool_t XdrEncoder::Int32 ( VML *vmls, char *name, int32_t &val ) 
{
  return xdr_int32_t(vmls->xmls,&val);
}
bool_t XdrEncoder::UnsignedInt32 ( VML *vmls, char *name, uint32_t &val ) 
{
  return xdr_u_int32_t(vmls->xmls,&val);
}
bool_t XdrEncoder::Int64 ( VML *vmls, char *name, int64_t &val ) 
{
  return xdr_int64_t(vmls->xmls,&val);
}
bool_t XdrEncoder::UnsignedInt64 ( VML *vmls, char *name, uint64_t &val ) 
{
  return xdr_u_int64_t(vmls->xmls,&val);
}
bool_t XdrEncoder::Char ( VML *vmls, char *name, char &val ) 
{
  return xdr_char(vmls->xmls,&val);
}
bool_t XdrEncoder::UnsignedChar ( VML *vmls, char *name, unsigned char &val ) 
{
  return xdr_u_char(vmls->xmls,&val);
}
bool_t XdrEncoder::Bool ( VML *vmls, char *name, bool_t &val ) 
{
  return xdr_bool_t(vmls->xmls,&val);
}
bool_t XdrEncoder::Double(VML *vmls, char *name, double &d)
{
  return xdr_double(vmls->xmls,&val);
}
bool_t XdrEncoder::Float (VML *vmls, char *name, float  &f)
{
  return xdr_float(vmls->xmls,&val);
}
bool_t XdrEncoder::Array     ( VML *vmls, char *type, char *name, 
				char * &vals, int &nvals,
				int &sizeofone, vmlproc_t do_one, int DoAlloc )
{

  unsigned int size=nvals;
  if ( !xdr_u_int(vmls->xdrs,&size) ) {
    return false;
  }
  
  if ( vmls->x_op == VML_DECODE ) {
    if ( ( vals == NULL ) || ( nvals < size ) ) {
      vals = malloc(size*sizeofone);
    }
  }

  for ( int i=0;i<size;i++ ) { 
    if ( ! do_one (vmls->xdrs, "elem",(void *)((unsigned long)vals+i*sizeofone)) ) { 
      return false;
    }
  }
  return true;
}

bool_t XdrEncoder::Bytes( VML *vmls, char *name, char *&vals, 
			   int &length) 
{
  return xdr_bytes(vmls->xdrs,&vals,(unsigned int *)&length,length);
}

bool_t XdrEncoder::String (VML *vmls, char *name, char *&str )
{
  return xdr_string(vmls->xdrs,&str,0);
}

bool_t XdrEncoder::Reference ( VML *vmls, char *type, char *name, 
				char *&ref, vmlproc_t do_ref,int sizeofone  )
{
  return Array(vmls,"reference",name,ref,length,sizeofone,do_ref);
}

bool_t XdrEncoder::StructBegin( VML *vmls, char *type, char *instance )
{
  return true;
}
bool_t XdrEncoder::StructEnd  ( VML *vmls, char *type, char *instance )
{
  return true;
}

bool_t XdrEncoder::ClassBegin( VML *vmls, char *type, char *instance )
{
  return true;
}
bool_t XdrEncoder::ClassEnd  ( VML *vmls, char *type, char *instance )
{
  return true;
}

#endif

