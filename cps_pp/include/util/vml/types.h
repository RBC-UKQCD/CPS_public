#ifndef _VML_TYPES_H
#define _VML_TYPES_H
#define VML_TRUE ((long) 1)
#define VML_FALSE ((long) 0)
typedef int bool_t;
typedef int enum_t;
#if TARGET == QCDOC
typedef char int8_t;
typedef short int16_t;
typedef int int32_t;
typedef long long int64_t;
typedef unsigned char uint8_t;
typedef unsigned short uint16_t;
typedef unsigned int long uint32_t;
typedef unsigned long long uint64_t;
#define MAX_NETOBJ_SZ 1024
typedef unsigned int u_int;
struct netobj {
        u_int   n_len;
        char    *n_bytes;
};
typedef struct netobj netobj;
#else
#include <rpc/types.h>
#include <rpc/rpc.h>
#endif
#endif
