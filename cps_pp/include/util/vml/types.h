#ifndef _VML_TYPES_H
#define _VML_TYPES_H
#include <conf.h>
#define VML_TRUE ((long) 1)
#define VML_FALSE ((long) 0)
#ifndef TRUE 
#define TRUE VML_TRUE
#define FALSE VML_FALSE
#endif
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#ifndef HAVE_BOOL_T
typedef int bool_t;
#endif

#ifndef HAVE_ENUM_T
typedef int enum_t;
#endif

#ifndef HAVE_INT8_T
typedef char int8_t;
#endif

#ifndef HAVE_INT16_T
typedef short int16_t;
#endif

#ifndef HAVE_INT32_T
typedef int int32_t;
#endif

#ifndef HAVE_INT64_T
typedef long long int64_t;
#endif

#ifndef HAVE_UINT8_T
typedef unsigned char uint8_t;
#endif

#ifndef HAVE_UINT16_T
typedef unsigned short uint16_t;
#endif

#ifndef HAVE_UINT32_T
typedef unsigned int uint32_t;
#endif
#ifndef HAVE_UINT64_T
typedef unsigned long long uint64_t;
#endif
#define MAX_NETOBJ_SZ 1024
#ifndef HAVE_UINT_T
typedef unsigned int u_int;
#endif
#ifndef HAVE_U_QUAD_T
typedef unsigned long long u_quad_t;
#endif
#ifndef HAVE_QUAD_T
typedef long long int quad_t;
#endif
#ifndef HAVE_NETOBJ
struct netobj {
        u_int   n_len;
        char    *n_bytes;
};
typedef struct netobj netobj;
#endif
#if 0
#include <rpc/types.h>
#include <rpc/rpc.h>
#endif
#endif
