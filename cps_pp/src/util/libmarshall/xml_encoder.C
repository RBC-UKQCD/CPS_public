#ifndef NO_CPS
#include <config.h>
#include <util/vml/vml_encoder.h>
#else
#include <vml_encoder.h>
#endif

#include <stdlib.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
#include <string.h>

#ifndef NO_CPS
CPS_START_NAMESPACE
USING_NAMESPACE_CPS
#endif
  /*
   * Interface to represent the basic types in XDR routines
   */
bool_t XmlEncoder::Void(void)
{
  return true;
}

#define SimpleEncode(val,name,xmltype,fmt) \
  char format[256]; \
  char line[256]; \
  sprintf(format,fmt,name,xmltype,xmltype);\
  if ( vmls->x_op == VML_ENCODE ) { \
    sprintf(line,format,val);\
    vmls->Puts(line);\
    return true;\
  } else if ( vmls->x_op == VML_DECODE ) { \
    vmls->Gets((char *)line,256);\
    if ( sscanf(line,format,&val) == 1 ) { \
      return true;\
    } else { \
      return false; \
    } \
  }  \
  return true;

#define DoIt(A,xmltype) SimpleEncode(val,name,#xmltype,\
       "<item><name>%s</name><value><%s>%%" #A "</%s></value></item>\n");
bool_t XmlEncoder::Short ( VML *vmls, char *name, short &val ) 
{
  DoIt (d,i2);
}
bool_t XmlEncoder::UnsignedShort (VML *vmls, char *name, unsigned short &val) 
{
  DoIt (u,ui2);
}
bool_t XmlEncoder::Enum ( VML *vmls, char *ename, char *name, char *&value ) 
{

  char format[512]; 
  char line[512]; 

  sprintf(format,"<item><name>%s</name><value><%s>%%s </%s></value></item>\n",name,ename,ename);

  if ( vmls->x_op == VML_ENCODE ) {

    sprintf(line,format, value); 
    vmls->Puts(line);
    return true;
  } else if ( vmls->x_op == VML_DECODE ) { 

    if ( !value ) value = (char *)malloc(256);

    vmls->Gets((char *)line,256);
    fprintf(stderr,"Line %s\n",line);
    fprintf(stderr,"Format %s\n",format);

    if ( sscanf(line,format,value) == 1 ) { 
      return true;
    } else { 
      fprintf(stderr,"Enum XML Bummer...\n");
      return false; 
    } 
  }  
  return true;
}
bool_t XmlEncoder::Int ( VML *vmls, char *name, int &val ) 
{
  DoIt (d,int);
}
bool_t XmlEncoder::UnsignedInt ( VML *vmls, char *name, unsigned int &val ) 
{
  DoIt (u,ui4);
}
bool_t XmlEncoder::Long ( VML *vmls, char *name, long &val ) 
{
  DoIt (ld,i4);
}
bool_t XmlEncoder::UnsignedLong ( VML *vmls, char *name, unsigned long &val )
{
  DoIt (lu,ui4);
}
bool_t XmlEncoder::LongLong ( VML *vmls, char *name, long long &val ) 
{
  DoIt (lld,i8);
}
bool_t XmlEncoder::UnsignedLongLong ( VML *vmls, char *name, 
				       unsigned long long &val ) 
{
  DoIt (llu,ui8);
}
bool_t XmlEncoder::Int8 ( VML *vmls, char *name, int8_t &val ) 
{
  DoIt (d,i1);
}
bool_t XmlEncoder::UnsignedInt8 ( VML *vmls, char *name, uint8_t &val ) 
{
  DoIt (u,ui1);
}
bool_t XmlEncoder::Int16 ( VML *vmls, char *name, int16_t &val ) 
{
  DoIt (d,i2);
}
bool_t XmlEncoder::UnsignedInt16 ( VML *vmls, char *name, uint16_t &val ) 
{
  DoIt (u,ui2);
}
bool_t XmlEncoder::Int32 ( VML *vmls, char *name, int32_t &val ) 
{
  DoIt (ld,i4);
}
bool_t XmlEncoder::UnsignedInt32 ( VML *vmls, char *name, uint32_t &val ) 
{
  DoIt (lu,ui4);
}
bool_t XmlEncoder::Int64 ( VML *vmls, char *name, int64_t &val ) 
{
  DoIt (lld,i8);
}
bool_t XmlEncoder::UnsignedInt64 ( VML *vmls, char *name, uint64_t &val ) 
{
  DoIt (llu,ui8);
}
bool_t XmlEncoder::Char ( VML *vmls, char *name, char &val ) 
{
  DoIt (c,char);
}
bool_t XmlEncoder::UnsignedChar ( VML *vmls, char *name, unsigned char &val ) 
{
  DoIt (u,ui1);
}
bool_t XmlEncoder::Bool ( VML *vmls, char *name, bool_t &val ) 
{
    DoIt (u,boolean);
}
bool_t XmlEncoder::Double(VML *vmls, char *name, double &val)
{
    DoIt (le,double);
}
bool_t XmlEncoder::Float (VML *vmls, char *name, float  &val)
{
    DoIt (e,float);
}
bool_t XmlEncoder::Array     ( VML *vmls, char *type, char *name, 
				char * &vals, int &nvals,
				int &sizeofone, vmlproc_t do_one, int DoAlloc )
{
  char tmp[256];
  char line[256];
  if ( vmls->x_op == VML_ENCODE ) { 
    sprintf(line,"<%s><name>%s</name><size>%d</size><data>\n",type,name,nvals);
    vmls->Puts(line);
    for(int i=0;i<nvals;i++) { 
      sprintf(tmp,"%s[%d]",name,i);
      do_one(vmls,tmp,(void *)((unsigned long)vals+i*sizeofone));
    }
    sprintf(line,"</data></%s>\n",type);
    vmls->Puts(line);
    return true;
  } else if ( vmls->x_op == VML_DECODE ) { 
    vmls->Gets((char *)line,256);
    sprintf(tmp,"<%s><name>%s</name><size>%%d</size><data>\n",type,name);
    sscanf(line,tmp,&nvals);
    if ( DoAlloc ) { 
      if ( nvals )
	vals = (char *)malloc(sizeofone*nvals);
      else 
	vals = NULL;
    }
    
    for(int i=0;i<nvals;i++) { 
      if ( vals == NULL ) return false;
      sprintf(tmp,"%s[%d]",name,i);
      do_one(vmls,tmp,(void *)((unsigned long)vals+i*sizeofone));
    }
    /*Skip the closing...*/
    vmls->Gets((char *)line,256);
    return true;

  }  
  return true;
}

bool_t XmlEncoder::Bytes( VML *vmls, char *name, char *&vals, 
			   int &length) 
{
  int sizeofone = 1;
  return Array(vmls,"bytes",name,vals,length,sizeofone,(vmlproc_t)vml_u_char);
}

bool_t XmlEncoder::String (VML *vmls, char *name, char *&value )
{

  char format[512]; 
  char line[512]; 

  sprintf(format,"<item><name>%s</name><value><string>%%s </string></value></item>\n",name);

  if ( vmls->x_op == VML_ENCODE ) {

    if( !value ) value = "";
    sprintf(line,format, value); 
    vmls->Puts(line);
    return true;
  } else if ( vmls->x_op == VML_DECODE ) { 

    if ( !value ) value = (char *)malloc(256);

    vmls->Gets((char *)line,256);

    if ( sscanf(line,format,value) == 1 ) { 
      return true;
    } else { 
      return false; 
    } 
  }  
  return true;
}

bool_t XmlEncoder::Reference ( VML *vmls, char *type, char *name, 
			       char *&ref, vmlproc_t do_ref ,int sizeofone )
{
  int length   = 1;
  return Array(vmls,"reference",name,ref,length,sizeofone,do_ref);
}

bool_t XmlEncoder::StructBegin( VML *vmls, char *type, char *instance )
{
  char line[256];
  if ( vmls->x_op == VML_ENCODE ) { 
    sprintf(line,"<struct><name>%s</name><instance>%s</instance><data>\n",type,instance);
    vmls->Puts(line);
    return true;
  } else if ( vmls->x_op == VML_DECODE ) { 
    vmls->Gets((char *)line,256);
    return true;
  }  
  return true;
}
bool_t XmlEncoder::StructEnd  ( VML *vmls, char *type, char *instance )
{
  char line[256];
  if ( vmls->x_op == VML_ENCODE ) { 
    vmls->Puts("</data></struct>\n");
    return true;
  } else if ( vmls->x_op == VML_DECODE ) { 
    vmls->Gets((char *)line,256);
    return true;
  }  
  return true;
}

bool_t XmlEncoder::ClassBegin( VML *vmls, char *type, char *instance )
{
  char line[256];
  if ( vmls->x_op == VML_ENCODE ) { 
    sprintf(line,"<class><name>%s</name><instance>%s</instance><data>\n",type,instance);
    vmls->Puts(line);
    return true;
  } else if ( vmls->x_op == VML_DECODE ) { 
    vmls->Gets((char *)line,256);
    return true;
  }  
  return true;
}
bool_t XmlEncoder::ClassEnd  ( VML *vmls, char *type, char *instance )
{
  char line[256];
  if ( vmls->x_op == VML_ENCODE ) { 
    vmls->Puts("</data></class>\n");
    return true;
  } else if ( vmls->x_op == VML_DECODE ) { 
    vmls->Gets((char *)line,256);
    return true;
  }  
  return true;
}
#ifndef NO_CPS
CPS_END_NAMESPACE
#endif
