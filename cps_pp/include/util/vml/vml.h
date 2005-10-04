/* Peter Boyle 2004.
 * vml.h, Serialisation routines to connect Sunrpc/XDR to VML
 */

#ifndef _VML_H
#define _VML_H 1

#include <config.h>
#include <util/vml/vml_config.h>
#include <util/vml/types.h>
#include <stdio.h>
#include <stdlib.h>

CPS_START_NAMESPACE
struct VML;
typedef struct VML VML;



#ifdef __cplusplus 
extern "C" { 
#endif
#ifdef AIX
typedef unsigned long long u_quad_t;
typedef long long quad_t;
#endif
/*
 *
 * Each data type provides a single procedure which takes two arguments:
 *
 *      bool_t
 *      vmlproc(vml, name, argresp)
 *              VML *vmls;
 *              char *name;
 *              <type> *argresp;
 *
 */

/*PAB Mar 2005... add mapping of enum's to the text string..*/
struct vml_enum_map { 
  char *enum_name;
  char *name;
  enum_t val;
};
char *vml_enum_string ( enum_t *val, struct vml_enum_map * );
enum_t *vml_enum_val  ( char *string, struct vml_enum_map * );

/*
 * Vml operations.  VML_ENCODE causes the type to be encoded into the
 * stream.  VML_DECODE causes the type to be extracted from the stream.
 * VML_FREE can be used to release the space allocated by an VML_DECODE
 * request.
 */
enum vml_op {
  VML_ENCODE = 0,
  VML_DECODE = 1,
  VML_FREE = 2
};

/*
 * The VML handle.
 * Contains operation which is being applied to the stream,
 * an operations vector for the particular implementation (e.g. see vml_mem.c),
 * and two private fields for the use of the particular implementation.
 */
  enum vml_type { 
    VML_MEM,
    VML_DESCRIPTOR,
    VML_STDIO
  };

 
/*
 * A vmlproc_t exists for each data type which is to be encoded or decoded.
 * adds the member name.
 *
 * The second argument to the vmlproc_t is a pointer to an opaque pointer.
 * The opaque pointer generally points to a structure of the data type
 * to be decoded.  If this pointer is 0, then the type routines should
 * allocate dynamic storage of the appropriate size and return it.
 * bool_t       (*vmlproc_t)(VML *, caddr_t *);
 */
typedef bool_t (*vmlproc_t)  (VML *, char *name ,void *,...);


/*
 * Support struct for discriminated unions.
 * You create an array of vmldiscrim structures, terminated with
 * a entry with a null procedure pointer.  The vml_union routine gets
 * the discriminant value and then searches the array of structures
 * for a matching value.  If a match is found the associated vml routine
 * is called to handle that part of the union.  If there is
 * no match, then a default routine may be called.
 * If there is no match and no default routine it is an error.
 */
#define NULL_vmlproc_t ((vmlproc_t)0)
struct vml_discrim
{
  int value;
  vmlproc_t proc;
};

/*
 * These are the "generic" vml routines.
 * None of these can have const applied because it's not possible to
 * know whether the call is a read or a write to the passed parameter
 * also, the VML structure is always updated by some of these calls.
 */
extern bool_t vml_void (void);
extern bool_t vml_short (VML *__vmls, char *name,short *__sp);
extern bool_t vml_u_short (VML *__vmls, char *name, u_short *__usp);
extern bool_t vml_int (VML *__vmls, char *name, int *__ip);
extern bool_t vml_u_int (VML *__vmls, char *name, u_int *__up);
extern bool_t vml_long (VML *__vmls, char *name, long *__lp);
extern bool_t vml_u_long (VML *__vmls, char *name, u_long *__ulp);
extern bool_t vml_hyper (VML *__vmls, char *name, quad_t *__llp);
extern bool_t vml_u_hyper (VML *__vmls, char *name, u_quad_t *__ullp);
extern bool_t vml_longlong_t (VML *__vmls, char *name, quad_t *__llp);
extern bool_t vml_u_longlong_t (VML *__vmls, char *name, u_quad_t *__ullp);
extern bool_t vml_int8_t (VML *__vmls, char *name, int8_t *__ip);
extern bool_t vml_uint8_t (VML *__vmls, char *name, uint8_t *__up);
extern bool_t vml_int16_t (VML *__vmls, char *name, int16_t *__ip);
extern bool_t vml_uint16_t (VML *__vmls, char *name, uint16_t *__up);
extern bool_t vml_int32_t (VML *__vmls, char *name, int32_t *__ip);
extern bool_t vml_uint32_t (VML *__vmls, char *name, uint32_t *__up);
extern bool_t vml_int64_t (VML *__vmls, char *name, int64_t *__ip);
extern bool_t vml_uint64_t (VML *__vmls, char *name, uint64_t *__up);
extern bool_t vml_bool (VML *__vmls, char *name, bool_t *__bp);
extern bool_t vml_enum (VML *__vmls, char *name, enum_t *__ep,struct vml_enum_map *list);
extern bool_t vml_array (VML * _vmls, char *name, caddr_t *__addrp, u_int *__sizep,
			      u_int __maxsize, u_int __elsize,
			      vmlproc_t __elproc);
extern bool_t vml_bytes (VML *__vmls, char *name, char **__cpp, u_int *__sizep,
			      u_int __maxsize);
extern bool_t vml_opaque (VML *__vmls, char *name, caddr_t __cp, u_int __cnt);
extern bool_t vml_string (VML *__vmls, char *name, char **__cpp, u_int __maxsize);
extern bool_t vml_union (VML *__vmls, char *name, enum_t *__dscmp, char *__unp,
			      const struct vml_discrim *__choices,
			      vmlproc_t dfault);
extern bool_t vml_char (VML *__vmls, char *name, char *__cp);
extern bool_t vml_u_char (VML *__vmls, char *name, u_char *__cp);
extern bool_t vml_vector (VML *__vmls, char *name, char *__basep, u_int __nelem,
			       u_int __elemsize, vmlproc_t __vml_elem);
extern bool_t vml_float (VML *__vmls, char *name, float *__fp);
extern bool_t vml_double (VML *__vmls, char *name, double *__dp);
extern bool_t vml_reference (VML *__vmls, char *name, caddr_t *__xpp, u_int __size,
				  vmlproc_t __proc);
extern bool_t vml_pointer (VML *__vmls, char *name, char **__objpp,
				u_int __obj_size, vmlproc_t __vml_obj);
extern bool_t vml_wrapstring (VML *__vmls, char *name, char **__cpp);
extern u_long vml_sizeof (vmlproc_t, void *);

/*
 * Common opaque bytes objects used by many rpc protocols;
 * declared here due to commonality.
 */
extern bool_t vml_netobj (VML *__vmls, char *name, struct netobj *__np);

/*
 * These are the public routines for the various implementations of
 * vml streams.
 */

/* VML I/O types */
bool_t vmlmem_create (VML *__vmls, char * __addr,
		      int __size, enum vml_op __xop);
bool_t vmlstdio_create (VML *__vmls, FILE *fp, enum vml_op __xop);
bool_t vmlfd_create (VML *__vmls, int fd, enum vml_op __xop);
bool_t vmlfile_create (VML *__vmls, char *file, enum vml_op __xop);

void vml_destroy (VML *__vmls);


/* VML markup language type selection*/
enum vml_markup  { VML_TEXT, VML_XML, VML_XDR };
extern void vml_markup_type(enum vml_markup type);

/* VML pseudo records for tcp */
extern void vmlrec_create (VML *__vmls, u_int __sendsize,
				u_int __recvsize, caddr_t __tcp_handle,
				int (*__readit) (char *, char *, int),
				int (*__writeit) (char *, char *, int));

/* make end of vml record */
extern bool_t vmlrec_endofrecord (VML *__vmls, bool_t __sendnow);

/* move to beginning of next record */
extern bool_t vmlrec_skiprecord (VML *__vmls);

/* true if no more input */
extern bool_t vmlrec_eof (VML *__vmls);

/* free memory buffers for vml */
extern void vml_free (vmlproc_t __proc, char *name, char *__objp);

extern void vml_struct_begin (VML *vmls, char *type, char *instance);
extern void vml_struct_end   (VML *vmls, char *type, char *instance);
extern void vml_class_begin (VML *vmls, char *type, char *instance);
extern void vml_class_end   (VML *vmls, char *type, char *instance);

#ifdef __cplusplus 
}
#endif
struct VML
  {
    public :
    enum vml_op x_op;		/* operation; fast additional param */

    /*
     * VML type dependent data
     */
    enum vml_type StreamType;
    FILE *x_fp;
    int   x_descriptor;

    caddr_t x_private;		/* pointer to private data */
    caddr_t x_base;		/* private used for position info */
    int x_handy;		/* extra private word */
    /*
     * Methods available 
     */
    int Puts(char *str);
    char *Gets(char *str,int sz);
    bool_t Create(char * bf,int length, vml_op op);
    bool_t Create(int fd, vml_op op);
    bool_t Create(char *file, vml_op op);
    bool_t Create(FILE *file, vml_op op);
    void Destroy(void);
  };

CPS_END_NAMESPACE
#endif /* rpc/vml.h */
