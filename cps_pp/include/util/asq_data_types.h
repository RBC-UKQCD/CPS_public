#ifndef ASQ_DATA_TYPES_H
#define ASQ_DATA_TYPES_H

#include <qcdocos/scu_dir_arg.h>
#ifdef SCIDAC
#ifdef ASQD_SINGLE
typedef float Float;
typedef float IFloat;
#define FNAME(pre,suf) pre##_F_##suf
#else
typedef double Float;
typedef double IFloat;
#define FNAME(pre,suf) pre##_D_##suf
#endif
#else
#include <precision.h>
#endif

#endif
