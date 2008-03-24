/*
 * Peter Boyle, 2004
 * vml.c, Generic VML routines implementation.
 */

#ifndef NO_CPS
#include <config.h>
#endif
#include <stdio.h>
#include <limits.h>
#include <string.h>
#define DEB(A) fprintf(stderr,A)

#ifndef NO_CPS
#include <util/qcdio.h>
#include <util/vml/types.h>
#include <util/vml/vml.h>
#include <util/vml/vml_encoder.h>
CPS_START_NAMESPACE
USING_NAMESPACE_CPS
#else
#include <rpc/types.h>
#include <rpc/rpc.h>
#include <vml.h>
#include <vml_encoder.h>
#endif

static GenericEncoder *VmlEncode = new TextEncoder;


extern "C" { 

/*
 * constants specific to the vml "protocol"
 */
#define LASTUNSIGNED	((u_int) 0-1)

/*
 * for unit alignment
 */
static const char vml_zero[4] = {0, 0, 0, 0};

/*
 * Free a data structure using VML
 * Not a filter, but a convenient utility nonetheless
 */
void
vml_free (vmlproc_t proc, char *name, char *objp)
{
  VML x;

  x.x_op = VML_FREE;
  (*proc) (&x, name, objp);
}


/*
 * VML nothing
 */
bool_t
vml_void (void)
{
  return VML_TRUE;
}

/*
 * VML integers
 */
bool_t
vml_int (VML *vmls, char *name,int *ip)
{
  return VmlEncode->Int(vmls,name,*ip);
}

/*
 * VML unsigned integers
 */
bool_t
vml_u_int (VML *vmls, char *name,  u_int *up)
{
  return VmlEncode->UnsignedInt(vmls,name,*up);
}

/*
 * VML long integers
 * same as vml_u_long - open coded to save a proc call!
 */
bool_t
vml_long (VML *vmls, char *name, long *lp)
{
  return VmlEncode->Long(vmls,name,*lp);
}

/*
 * VML unsigned long integers
 * same as vml_long - open coded to save a proc call!
 */
bool_t
vml_u_long (VML *vmls, char *name, u_long *ulp)
{
  return VmlEncode->UnsignedLong(vmls,name,*ulp);
}

/*
 * VML hyper integers
 * same as vml_u_hyper - open coded to save a proc call!
 */
bool_t
vml_hyper (VML *vmls, char *name, quad_t *llp)
{
  bool_t ret;
  long long lll = *llp;
  ret = VmlEncode->LongLong(vmls,name,lll);
  *llp = lll;
  return ret;
}


/*
 * VML hyper integers
 * same as vml_hyper - open coded to save a proc call!
 */
bool_t
vml_u_hyper (VML *vmls, char *name, u_quad_t *ullp)
{
  bool_t ret;
  unsigned long long lll = *ullp;
  ret = VmlEncode->UnsignedLongLong(vmls,name,lll);
  *ullp = lll;
  return ret;
}

bool_t
vml_longlong_t (VML *vmls, char *name, quad_t *llp)
{
  return vml_hyper (vmls, name, llp);
}

bool_t
vml_u_longlong_t (VML *vmls, char *name, u_quad_t *ullp)
{
  return vml_u_hyper (vmls, name, ullp);
}

/*
 * VML short integers
 */
bool_t
vml_short (VML *vmls, char *name, short *sp)
{
  return VmlEncode->Short(vmls,name,*sp);
}

/*
 * VML unsigned short integers
 */
bool_t
vml_u_short (VML *vmls, char *name, u_short *usp)
{
  return VmlEncode->UnsignedShort(vmls,name,*usp);
}


/*
 * VML a char
 */
bool_t
vml_char (VML *vmls, char *name, char *cp)
{
  return VmlEncode->Char(vmls,name,*cp);
}

/*
 * VML an unsigned char
 */
bool_t
vml_u_char (VML *vmls, char *name, u_char *cp)
{
  return VmlEncode->UnsignedChar(vmls,name,*cp);
}

/*
 * VML booleans
 */
bool_t
vml_bool (VML *vmls, char *name, bool_t *bp)
{
  return VmlEncode->Bool(vmls,name,*bp);
}
char *vml_enum_string ( enum_t *val, struct vml_enum_map *map )
{
  for ( int i=0;map[i].name!=NULL; i++ ) {
    if ( map[i].val == *val ) {
      return map[i].name;
    }
  }
  fprintf(stderr,"logic bomb in vml_enum_string\n");
  // exit(-1);
  return NULL;
}
enum_t *vml_enum_val  ( char *string, struct vml_enum_map * map)
{
  for ( int i=0;map[i].name!=NULL; i++ ) {
    if ( !strcmp(string,map[i].name) ) {
      return &map[i].val;
    }
  }
  fprintf(stderr,"logic bomb in vml_enum_val in \"%s\" \n",string); exit(-1);
}
/*
 * VML enumerations
 */
bool_t
vml_enum (VML *vmls, char *name, enum_t *ep,struct vml_enum_map *list)
{
  char val_buf[256];
  char * val_str = val_buf;
  /*Fix me... error codes?*/
  if ( vmls->x_op == VML_ENCODE ) { 
    val_str = vml_enum_string(ep,list);
	if(!val_str) {fprintf(stderr,"%d not in %s\n",ep,list->enum_name);exit(-1);}
  }
  
  VmlEncode->Enum(vmls,list[0].enum_name,name,val_str);
  if ( vmls->x_op == VML_DECODE ) { 
    *ep = *(vml_enum_val(val_str,list));
  }
  return VML_TRUE;
}

/*
 * VML opaque data
 * Allows the specification of a fixed size sequence of opaque bytes.
 * cp points to the opaque object and cnt gives the byte length.
 */
#define BYTES_PER_VML_UNIT 4
bool_t
vml_opaque (VML *vmls, char *name, caddr_t cp, u_int cnt)
{
  int icnt = cnt;
  char *tmp = (char *)cp;
  VmlEncode->Bytes(vmls,name,tmp,icnt);
  cnt = icnt;
  return VML_TRUE;
}

/*
 * VML counted bytes
 * *cpp is a pointer to the bytes, *sizep is the count.
 * If *cpp is NULL maxsize bytes are allocated
 */
bool_t
vml_bytes (
  VML *vmls,
  char *name,
  char **cpp,
  u_int *sizep,
  u_int maxsize)
{
  int icnt = *sizep;
  VmlEncode->Bytes(vmls,name,*cpp,icnt);
  *sizep = icnt;
  return VML_TRUE;
}

/*
 * Implemented here due to commonality of the object.
 */
bool_t
vml_netobj (
	    VML *vmls,
	    char *name,
	    struct netobj *np)
{

  return vml_bytes (vmls, name, &np->n_bytes, &np->n_len, MAX_NETOBJ_SZ);
}

/*
 * VML a discriminated union
 * Support routine for discriminated unions.
 * You create an array of vmldiscrim structures, terminated with
 * an entry with a null procedure pointer.  The routine gets
 * the discriminant value and then searches the array of vmldiscrims
 * looking for that value.  It calls the procedure given in the vmldiscrim
 * to handle the discriminant.  If there is no specific routine a default
 * routine may be called.
 * If there is no specific or default routine an error is returned.
 */
bool_t
vml_union (
	   VML *vmls,
	   char *name,
     enum_t *dscmp,	/* enum to decide which arm to work on */
     char *unp,		/* the union itself */
     const struct vml_discrim *choices,	/* [value, vml proc] for each arm */
     vmlproc_t dfault)		/* default vml routine */
{
  enum_t dscm;
  long tmp;

  if ( dscmp ) tmp = *dscmp;
  /*
   * we deal with the discriminator;  it's an enum
   */
  if (!vml_long (vmls, name, &tmp))
    {
      return VML_FALSE;
    }
  dscm = tmp;
  *dscmp = tmp;

  /*
   * search choices for a value that matches the discriminator.
   * if we find one, execute the vml routine for that value.
   */
  for (; choices->proc != NULL_vmlproc_t; choices++)
    {
      if (choices->value == dscm)
	return (*(choices->proc)) (vmls, name, unp, LASTUNSIGNED);
    }

  /*
   * no match - execute the default vml routine if there is one
   */
  return ((dfault == NULL_vmlproc_t) ? VML_FALSE :
	  (*dfault) (vmls, name, unp, LASTUNSIGNED));
}


/*
 * Non-portable vml primitives.
 * Care should be taken when moving these routines to new architectures.
 */


/*
 * VML null terminated ASCII strings
 * vml_string deals with "C strings" - arrays of bytes that are
 * terminated by a NULL character.  The parameter cpp references a
 * pointer to storage; If the pointer is null, then the necessary
 * storage is allocated.  The last parameter is the max allowed length
 * of the string as specified by a protocol.
 */
bool_t
vml_string (
	    VML *vmls,
	    char *name,
	    char **cpp,
	    u_int maxsize)
{
  return VmlEncode->String(vmls,name,*cpp);
}

/*
 * Wrapper for vml_string that can be called directly from
 * routines like clnt_call
 */
bool_t
vml_wrapstring (
		VML *vmls,
		char * name,
		char **cpp)
{
  if (vml_string (vmls, name, cpp, LASTUNSIGNED))
    {
      return VML_TRUE;
    }
  return VML_FALSE;
}

bool_t
vml_float(VML *vmls,
	  char *name,
	  float *fp)
{
  return VmlEncode->Float(vmls,name,*fp);
}


bool_t
vml_double(VML *vmls,
	   char *name,
	   double *dp)
{
  return VmlEncode->Double(vmls,name,*dp);
}

/* VML 64bit integers */
bool_t
vml_int64_t (VML *vmls, char *name, int64_t *ip)
{
  return VmlEncode->Int64(vmls,name,*ip);
}

/* VML 64bit unsigned integers */
bool_t
vml_uint64_t (VML *vmls, char *name, uint64_t *uip)
{
  return VmlEncode->UnsignedInt64(vmls,name,*uip);
}

/* VML 32bit integers */
bool_t
vml_int32_t (VML *vmls, char *name, int32_t *lp)
{
  return VmlEncode->Int32(vmls,name,*lp);
}

/* VML 32bit unsigned integers */
bool_t
vml_uint32_t (VML *vmls, char *name, uint32_t *ulp)
{
  return VmlEncode->UnsignedInt32(vmls,name,*ulp);
}

/* VML 16bit integers */
bool_t
vml_int16_t (VML *vmls, char *name, int16_t *ip)
{
  return VmlEncode->Int16(vmls,name,*ip);
}

/* VML 16bit unsigned integers */
bool_t
vml_uint16_t (VML *vmls, char *name ,uint16_t *uip)
{
  return VmlEncode->UnsignedInt16(vmls,name,*uip);
}

/* VML 8bit integers */
bool_t
vml_int8_t (VML *vmls, char *name, int8_t *ip)
{
  return VmlEncode->Int8(vmls,name,*ip);
}

/* VML 8bit unsigned integers */
bool_t
vml_uint8_t (VML *vmls, char *name, uint8_t *uip)
{
  return VmlEncode->UnsignedInt8(vmls,name,*uip);
}



bool_t
vml_array (
	   VML *vmls,
	   char *name,
	   caddr_t *addrp,		/* array pointer */
	   u_int *sizep,		/* number of elements */
	   u_int maxsize,		/* max numberof elements */
	   u_int elsize,		/* size in bytes of each element */
	   vmlproc_t elproc)	/* vml routine to handle each element */
{
  int nvals  = *sizep;
  int sz = elsize;
  char *tmp = (char *)(* addrp);
  bool_t ret = VmlEncode->Array(vmls,(char *)"Array",
				name,tmp,nvals,sz,elproc);
  *addrp = tmp;
  *sizep = nvals;
  return ret;
}
bool_t
vml_vector (VML *vmls,
	    char *name,
	    char *basep,
	    u_int nelem,
	    u_int elemsize,
	    vmlproc_t vml_elem)
{
  int nvals  = nelem;
  int sz = elemsize;

  /* For VECTOR
   * Disable allocate, and never pass on a free request
   */
  if ( vmls->x_op == VML_FREE ) return true;

  bool_t ret = VmlEncode->Array(vmls,(char *)"Vector",
				name,basep,nvals,
				sz,
				vml_elem,0);
  return ret;
}
bool_t
vml_reference (
	       VML *vmls,
	       char *name,
	       caddr_t *pp, /* the pointer to work on */
	       u_int size,  /* size of the object pointed to */
	       vmlproc_t proc) /* vml routine to handle the object */
{
  char *tmp = (char *)*pp;
  char *type = "Reference";
  bool_t ret = VmlEncode->Reference(vmls,type,name,tmp,proc,size);
  *pp = tmp;
  return ret;
}
bool_t
vml_pointer (
	     VML *vmls,
	     char *name,
	     char **objpp,
	     u_int obj_size,
	     vmlproc_t vml_obj)
{

  bool_t more_data;

  more_data = (*objpp != NULL);
  if (!vml_bool (vmls, name, &more_data))
    {
      return VML_FALSE;
    }
  if (!more_data)
    {
      *objpp = NULL;
      return VML_TRUE;
    }
  return vml_reference (vmls, name, objpp, obj_size, vml_obj);
}


void vml_struct_begin (VML *vmls, char *type, char *instance)
{
  VmlEncode->StructBegin(vmls,type,instance);
}
void vml_struct_end   (VML *vmls, char *type, char *instance)
{
  VmlEncode->StructEnd(vmls,type,instance);
}
void vml_class_begin (VML *vmls, char *type, char *instance)
{
  VmlEncode->ClassBegin(vmls,type,instance);
}
void vml_class_end   (VML *vmls, char *type, char *instance)
{
  VmlEncode->ClassEnd(vmls,type,instance);
}

}

bool_t VML::Puts(char *string)
{
  switch ( StreamType ) { 
  case VML_STDIO:
  case VML_DESCRIPTOR:
#ifndef NO_CPS
    Fprintf(x_fp,string);
#else
    fprintf(x_fp,string);
    fflush(x_fp);
#endif
    return true;
    break;
  case VML_MEM:
    int len = strlen(string);
    if ( x_handy - len  > 0 ) { 
      strcpy((char *)x_private,string);
      x_private = x_private + len;
      x_handy   = x_handy - len;
      return true;
    } 
    return false;
    break;
  }
  return false;
}
char *VML::Gets(char *string,int n)
{
  char *ret;
  switch ( StreamType ) { 
  case VML_STDIO:
  case VML_DESCRIPTOR:
    ret = fgets(string, n, x_fp);
    if ( ret == NULL ) DEB("bad fgets\n");
    return ret;
    break;
  case VML_MEM:
    char *cp=(char *)x_private;
    int count=0;
    unsigned long x_private_new;
    while( (*cp!='\0') && (*cp!='\n') && count < n-2 ) { 
      string[count++] = *cp++;
    }
    string[count++] = '\n';
    string[count++] = '\0';
    if (*cp=='\n') cp++;
    
    x_private_new = (unsigned long)cp;
    x_handy= x_handy - (x_private_new - (unsigned long)x_private);
    x_private = (char *)x_private_new;
    if ( x_handy >= 0 ) return string;
    else return NULL;
    break;
  }
  DEB("Bad StreamType in Gets\n");
  return NULL;
}
bool_t VML::Create(char *buf,int length, enum vml_op op)
{
  x_op = op;
  x_private = x_base = buf;
  x_handy = length;
  StreamType = VML_MEM;
  return true;
}
bool_t VML::Create(int fd,enum vml_op op)
{
  
  FILE *fp = NULL;
  if ( op == VML_ENCODE ) fp = fdopen(fd,"w");
  if ( op == VML_DECODE ) fp = fdopen(fd,"r");
  if ( fp == NULL ) return VML_FALSE;
  return Create(fp,op);
}
bool_t VML::Create(char *file, enum vml_op op)
{
  FILE *fp = NULL;
#ifndef NO_CPS
  if ( op == VML_ENCODE ) { 
    fp  = Fopen(file,"w");
  } else if ( op == VML_DECODE ) { 
    fp  = fopen(file,"r");
  }
#else
  if ( op == VML_ENCODE ) { 
    fp  = fopen(file,"w");
  } else if ( op == VML_DECODE ) { 
    fp  = fopen(file,"r");
  }
#endif
  if ( fp == NULL ) return VML_FALSE;
  return Create(fp,op);

}
bool_t VML::Create(FILE *file, enum vml_op op)
{
  x_fp = file;
  x_op = op;
  StreamType = VML_STDIO;
  if ( file ) return true;
  return false;
}

void VML::Destroy(void) 
{
  if ( StreamType == VML_MEM ) return;
#ifndef NO_CPS
  else if (x_fp) Fclose(x_fp);
#else
  else if (x_fp) fclose(x_fp);
#endif
}
bool_t vmlmem_create (VML *__vmls, char * __addr,
				int __size, enum vml_op __xop)
{
  return __vmls->Create(__addr,__size,__xop);
}
bool_t vmlstdio_create (VML *__vmls, FILE *fp, enum vml_op __xop)
{ 
  return __vmls->Create(fp,__xop);
}
bool_t vmlfd_create (VML *__vmls, int fd, enum vml_op __xop)
{
  return __vmls->Create(fd,__xop);
}

bool_t vmlfile_create (VML *__vmls, char *file, enum vml_op __xop)
{
  return __vmls->Create(file,__xop);
}
void vml_markup_type(enum vml_markup type)
{
  switch (type) {
  case VML_TEXT:
     VmlEncode = new TextEncoder;
     break;
  case VML_XML:
     VmlEncode = new XmlEncoder;
     break;
  default:
    exit(-1);
  };

  return ;
}
void vml_destroy (VML *__vmls)
{
  __vmls->Destroy();
}
#ifndef NO_CPS
CPS_END_NAMESPACE
#endif
