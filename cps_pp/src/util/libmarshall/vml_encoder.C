#ifndef NO_CPS
#include <config.h>
#include <util/vml/vml_encoder.h>
CPS_START_NAMESPACE
USING_NAMESPACE_CPS
#else
#include <rpc/types.h>
#include <rpc/rpc.h>
#include <vml.h>
#include <vml_encoder.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
#include <string.h>
  /* Peter Boyle 2004
   * Interface to represent the basic types in XDR routines
   */
#define DEB(A) fprintf(stderr,A);
bool_t TextEncoder::Void(void)
{
  return true;
}

#define SimpleEncode(val,name,fmt) \
  char format[256]; \
  char line[256]; \
  sprintf(format,fmt,name);\
  if ( vmls->x_op == VML_ENCODE ) { \
    sprintf(line,format,val);\
    if(!vmls->Puts(line)) return false;\
    return true;\
  } else if ( vmls->x_op == VML_DECODE ) { \
    if (!vmls->Gets((char *)line,256)) { DEB("Gets"); return false;}\
    if ( sscanf(line,format,&val) == 1 ) { \
      return true;\
    } else { \
      DEB("scanf"); DEB(name); DEB(line); return false; \
    } \
  }  \
  return true;

#define SimpleEncodeDouble(val,name,fmt1,fmt2) \
  char format[256]; \
  char line[256]; \
  if ( vmls->x_op == VML_ENCODE ) { \
    sprintf(format,fmt1,name);\
    sprintf(line,format,val);\
    if(!vmls->Puts(line)) return false;\
    return true;\
  } else if ( vmls->x_op == VML_DECODE ) { \
    sprintf(format,fmt2,name);\
    if (!vmls->Gets((char *)line,256)) { DEB("Gets"); return false;}\
    if ( sscanf(line,format,&val) == 1 ) { \
      return true;\
    } else { \
      DEB("scanf"); DEB(name); DEB(line); return false; \
    } \
  }  \
  return true;

bool_t TextEncoder::Enum ( VML *vmls,char *ename, char *name, char *&value)
{
  char format[256]; 
  char line[256]; 
  sprintf(format,"%s %s = %%s\n",ename,name);
  if ( vmls->x_op == VML_ENCODE ) { 
    sprintf(line,format,value);
    if(!vmls->Puts(line)) return false;
    return true;
  } else if ( vmls->x_op == VML_DECODE ) { 
    if (!vmls->Gets((char *)line,256)) return false;
    value = (char *) malloc(strlen(line));

    if ( sscanf(line,format,value) == 1 ) { 
      return true;
    } else { 
      value = NULL; DEB("enum"); DEB(name);DEB(line);
      return false; 
    } 
  }  
  return true;
}
bool_t TextEncoder::Short ( VML *vmls, char *name, short &val ) 
{
  SimpleEncode (val,name,"short %s = %%d\n");
}
bool_t TextEncoder::UnsignedShort (VML *vmls, char *name, unsigned short &val) 
{
  SimpleEncode (val,name,"unsigned short %s = %%d\n");
}
bool_t TextEncoder::Int ( VML *vmls, char *name, int &val ) 
{
  SimpleEncode (val,name,"int %s = %%d\n");
}
bool_t TextEncoder::UnsignedInt ( VML *vmls, char *name, unsigned int &val ) 
{
  SimpleEncode (val,name,"unsigned int %s = %%d\n");
}
bool_t TextEncoder::Long ( VML *vmls, char *name, long &val ) 
{
  SimpleEncode (val,name,"long %s = %%ld\n");
}
bool_t TextEncoder::UnsignedLong ( VML *vmls, char *name, unsigned long &val )
{
  SimpleEncode (val,name,"unsigned long %s = 0x%%lx\n");
}
bool_t TextEncoder::LongLong ( VML *vmls, char *name, long long &val ) 
{
  SimpleEncode (val,name,"long long %s = %%lld\n");
}
bool_t TextEncoder::UnsignedLongLong ( VML *vmls, char *name, 
				       unsigned long long &val ) 
{
  SimpleEncode (val,name,"unsigned long long %s = %%llu\n");
}
bool_t TextEncoder::Int8 ( VML *vmls, char *name, int8_t &val ) 
{
  SimpleEncode (val,name,"int8 %s = %%d\n");
}
bool_t TextEncoder::UnsignedInt8 ( VML *vmls, char *name, uint8_t &val ) 
{
  SimpleEncode (val,name,"uint8 %s = %%u\n");
}
bool_t TextEncoder::Int16 ( VML *vmls, char *name, int16_t &val ) 
{
  SimpleEncode (val,name,"int16 %s = %%d\n");
}
bool_t TextEncoder::UnsignedInt16 ( VML *vmls, char *name, uint16_t &val ) 
{
  SimpleEncode (val,name,"uint16 %s = %%u\n");
}
bool_t TextEncoder::Int32 ( VML *vmls, char *name, int32_t &val ) 
{
  SimpleEncode (val,name,"int32 %s = %%d\n");
}
bool_t TextEncoder::UnsignedInt32 ( VML *vmls, char *name, uint32_t &val ) 
{
  SimpleEncode (val,name,"int32 %s = %%x\n");
}
bool_t TextEncoder::Int64 ( VML *vmls, char *name, int64_t &val ) 
{
  SimpleEncode (val,name,"int64 %s = %%ll\n");
}
bool_t TextEncoder::UnsignedInt64 ( VML *vmls, char *name, uint64_t &val ) 
{
  SimpleEncode (val,name,"uint64 %s = %%llu\n");
}
bool_t TextEncoder::Char ( VML *vmls, char *name, char &val ) 
{
  SimpleEncode (val,name,"char %s = %%c\n");
}
bool_t TextEncoder::UnsignedChar ( VML *vmls, char *name, unsigned char &val ) 
{
  SimpleEncode (val,name,"unsigned char %s = %%d\n");
}
bool_t TextEncoder::Bool ( VML *vmls, char *name, bool_t &val ) 
{
  SimpleEncode (val,name,"bool %s = %%d\n");
}
bool_t TextEncoder::Double(VML *vmls, char *name, double &d)
{
  SimpleEncodeDouble (d,name,"double %s = %%24.16le\n","double %s = %%le\n");
}
bool_t TextEncoder::Float (VML *vmls, char *name, float  &f)
{
  SimpleEncode (f,name,"float %s = %%e\n");
}
bool_t TextEncoder::Array     ( VML *vmls, char *type, char *name, 
				char * &vals, int &nvals,
				int &sizeofone, vmlproc_t do_one, int DoAlloc )
{
  char tmp[256];
  char line[256];
  if ( vmls->x_op == VML_ENCODE ) { 

    sprintf(line,"%s %s[%d] = { \n",type,name,nvals);
    if(!vmls->Puts(line)) return false;
    for(int i=0;i<nvals;i++) { 
      sprintf(tmp,"%s[%d]",name,i);
      do_one(vmls,tmp,(void *)((unsigned long)vals+i*sizeofone));
    }
    if(!vmls->Puts("}\n")) return false;
    return true;
  } else if ( vmls->x_op == VML_DECODE ) { 
    if (!vmls->Gets((char *)line,256)) return false;
    sprintf(tmp,"%s %s[%%d] = { \n",type,name);
    sscanf(line,tmp,&nvals);

    if ( DoAlloc ) { 
      if ( nvals )
	vals = (char *)malloc(sizeofone*nvals);
      else 
	vals = NULL;
    }
    
    for(int i=0;i<nvals;i++) { 
      if ( vals == NULL ) { DEB("NULL array\n"); return false;};
      sprintf(tmp,"%s[%d]",name,i);
      bool_t ret = do_one(vmls,tmp,(void *)((unsigned long)vals+i*sizeofone));
	  if (ret == false) return false;
    }
    /*Skip the closing bracket...*/
    if (!vmls->Gets((char *)line,256)) { DEB("Array close\n"); return false;}
    return true;

  }  
  return true;
}

bool_t TextEncoder::Bytes( VML *vmls, char *name, char *&vals, 
			   int &length) 
{
  int sizeofone = 1;
  return Array(vmls,"byte",name,vals,length,sizeofone,(vmlproc_t)vml_u_char);
}

bool_t TextEncoder::String (VML *vmls, char *name, char *&str )
{
  int sizeofone = 1;
  int length = 0;
  char line[1024];
  char fmt[128];

  if ( vmls->x_op == VML_ENCODE ) { 

   char * temp_str = str;
    if ( temp_str == NULL ) temp_str = "";
    sprintf(line,"string %s = \"%s\"\n",name,temp_str);
    if(!vmls->Puts(line)) return false;

  } else { 

    if (!vmls->Gets((char *)line,1024)) {DEB("Gets failed");return false;}
    str = (char *) malloc(strlen(line));
    sprintf(fmt,"string %s = \"%%s\"\n",name);
    sscanf(line,fmt,str);
    if ( str[strlen(str)-1] == '\"' ) { 
      str[strlen(str)-1] = '\0';
    }
  }
  return true;
}

bool_t TextEncoder::Reference ( VML *vmls, char *type, char *name, 
				char *&ref, vmlproc_t do_ref,int sizeofone  )
{
  int length   = 1;
  return Array(vmls,"reference",name,ref,length,sizeofone,do_ref);
}

bool_t TextEncoder::StructBegin( VML *vmls, char *type, char *instance )
{
  char line[256];
  if ( vmls->x_op == VML_ENCODE ) { 
    sprintf(line,"struct %s %s = {\n",type,instance);
    if(!vmls->Puts(line)) return false;
    return true;
  } else if ( vmls->x_op == VML_DECODE ) { 
    if (!vmls->Gets((char *)line,256)) return false;
    return true;
  }  
  return true;
}
bool_t TextEncoder::StructEnd  ( VML *vmls, char *type, char *instance )
{
  char line[256];
  if ( vmls->x_op == VML_ENCODE ) { 
    if (!vmls->Puts("}\n"))return false;
    return true;
  } else if ( vmls->x_op == VML_DECODE ) { 
    if (!vmls->Gets((char *)line,256)) return false;
    return true;
  }  
  return true;
}

bool_t TextEncoder::ClassBegin( VML *vmls, char *type, char *instance )
{
  char line[256];
  if ( vmls->x_op == VML_ENCODE ) { 
    sprintf(line,"class %s %s = {\n",type,instance);
    if(!vmls->Puts(line)) return false;
    return true;
  } else if ( vmls->x_op == VML_DECODE ) { 
    if (!vmls->Gets((char *)line,256)) return false;
    return true;
  }  
  return true;
}
bool_t TextEncoder::ClassEnd  ( VML *vmls, char *type, char *instance )
{
  char line[256];
  if ( vmls->x_op == VML_ENCODE ) { 
    if (!vmls->Puts("}\n"))return false;
    return true;
  } else if ( vmls->x_op == VML_DECODE ) { 
    if (!vmls->Gets((char *)line,256)) return false;
    return true;
  }  
  return true;
}

#ifndef NO_CPS
CPS_END_NAMESPACE
#endif
