#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  BG/L operating system interface declarations and includes
*/

/*******************************************************************
*
* bgl_sys.h
*
* BG/L operating system interface declarations and includes
*
* The minimum amount of definitions needed to do simple maipulations
* directly on the BG/L hardware.
*
* Several of the definitions and functions below were copied from 
* (in alphabetical order) G. Almasi, J. Castanos and M. Giampapa
*
*******************************************************************/

#ifndef INCLUDED_BGL_SYS_H
#define INCLUDED_BGL_SYS_H

CPS_END_NAMESPACE


CPS_START_NAMESPACE

//------------------------------------------------------------------
// #define memory addresses
//------------------------------------------------------------------


//------------------------------------------------------------------
// #define sizes
//------------------------------------------------------------------
#define BGL_L1_ALIGNSIZE  32
#define BGL_QUAD_ALIGNSIZE  16
#define BGL_DWORD_ALIGNSIZE 8


//------------------------------------------------------------------
// # define the basic Quad operations macros
//------------------------------------------------------------------
#define STR(s) #s

#define QuadLoad(v,f)  asm volatile( "lfpdx " STR(f) ",0,%0" :: "r" (v) : "fr" STR(f) )

#define QuadStore(v,f) asm volatile( "stfpdx " STR(f) ",0,%0" :: "r" (v) : "memory" )

#define QuadLoadU(v,f,u)  asm volatile( "lfpdux " STR(f) ",%0,  %1" \
                                        :                           \
                                        : "r" (v),                  \
                                          "r" (u)                   \
                                        : "fr" STR(f) )             \

#define QuadLoadNU(v,f,u)  asm volatile( "lfpdx " STR(f) ",%0,  %1" \
                                        :                           \
                                        : "r" (v),                  \
                                          "r" (u)                   \
                                        : "fr" STR(f) )             \

#define QuadStoreU(v,f,u) asm volatile( "stfpdux " STR(f) ",%0, %1" \
                                        :                           \
                                        : "r" (v),                  \
                                          "r" (u)                   \
                                        : "memory" )                \

#define QuadStoreNU(v,f,u) asm volatile( "stfpdx " STR(f) ",%0, %1" \
                                        :                           \
                                        : "r" (v),                  \
                                          "r" (u)                   \
                                        : "memory" )                \

#define QuadMove(src,dst,f)  asm volatile( "lfpdx  " STR(f) ",0,%0;" \
                                           "stfpdx " STR(f) ",0,%1;" \
                                           :                         \
                                           : "r" (src),              \
                                             "r" (dst)               \
                                           : "memory",               \
                                             "fr" STR(f) )

#define DCBT(src) asm volatile( "dcbt 0, %0 " \
                                :             \
                                : "r" (src) ) \
 
#define DCBTST(src) asm volatile( "dcbtst 0, %0 " \
                                  :             \
                                  : "r" (src) ) \
 
#define DCBZ(src) asm volatile( "dcbz 0, %0 " \
                                :             \
                                : "r" (src) ) \
 

#define DCBI(src) asm volatile( "dcbi 0, %0 " \
                                :             \
                                : "r" (src) ) \
 

#define DCBF(src) asm volatile( "dcbf 0, %0 " \
                                :             \
                                : "r" (src) ) \
 


//------------------------------------------------------------------
// Type definitions
//------------------------------------------------------------------

//------------------------------------------------------------------
// The definition of a quad word = 128 bits -- 16 Byte alligned
//------------------------------------------------------------------
typedef struct _BGLQuad_tag
{
  unsigned w0;
  unsigned w1;
  unsigned w2;
  unsigned w3;
} BGLQuad __attribute__((aligned(BGL_QUAD_ALIGNSIZE)));



//------------------------------------------------------------------
// Inline functions
//------------------------------------------------------------------

//------------------------------------------------------------------
// Ensure that order-of-execution for memory accesses is respected.
//------------------------------------------------------------------
static inline void BGLCPS_Msync()
{
  asm volatile("sync" ::: "memory");
}





#endif


CPS_END_NAMESPACE
