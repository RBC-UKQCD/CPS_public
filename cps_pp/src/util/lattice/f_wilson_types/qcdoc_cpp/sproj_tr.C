#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  gamma matrix projection and spin trace routines.

  Used by derivatives of the FwilsonTypes class.
  
  $Id: sproj_tr.C,v 1.1 2004-11-20 00:06:48 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-11-20 00:06:48 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_wilson_types/qcdoc_cpp/sproj_tr.C,v 1.1 2004-11-20 00:06:48 chulwoo Exp $
//  $Id: sproj_tr.C,v 1.1 2004-11-20 00:06:48 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_wilson_types/qcdoc_cpp/sproj_tr.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
///------------------------------------------------------------------
//
// sproj_tr.C
//
// These functions return a color matrix in f constructed from
// the spinors v, w using: 
// f_(i,j) = Tr_spin[ (1 +/- gamma_mu) v_i w^dag_j ]
//
// num_blk is the number of spinors v, w. The routines 
// accumulate the sum over spinors in f.
//
// stride is the number of IFloats between spinors (a stride = 0
// means that the spinors are consecutive in memory)
// 
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/sproj_tr.h>
CPS_START_NAMESPACE

#ifndef _GCC_
#ifdef _XLC_
#include <builtins.h>
#define	  __builtin_prefetch(A) __dcbt(A)
#else
#define	  __builtin_prefetch(A)
#endif
#endif
 

enum{MATRIX_SIZE = 18};
#define SPROJ_START() \
  int row, col, blk ; \
  IFloat *tf, *tv, *tw ; \ 
  register IFloat f00r=0.,f00i=0.,f01r=0.,f01i=0.,f02r=0.,f02i=0.;\
  register IFloat f10r=0.,f10i=0.,f11r=0.,f11i=0.,f12r=0.,f12i=0.;\
  register IFloat f20r=0.,f20i=0.,f21r=0.,f21i=0.,f22r=0.,f22i=0.;\
  IFloat v1_0,v1_1,v1_2,v1_3,v1_4,v1_5;\
  IFloat w1_0,w1_1,w1_2,w1_3,w1_4,w1_5;\
  IFloat *vnext = v+24 + v_stride;\
  IFloat *wnext = w+24 + w_stride;\
  for (blk=0; blk<num_blk; blk++) {\
    tv = v ;\
    tw = w ;

#define SPROJ_COMP() \
	f00r += v1_0*w1_0+v1_1*w1_1;\
	f00i += v1_1*w1_0-v1_0*w1_1;\
	f01r += v1_0*w1_2+v1_1*w1_3;\
	f01i += v1_1*w1_2-v1_0*w1_3;\
	f02r += v1_0*w1_4+v1_1*w1_5;\
	f02i += v1_1*w1_4-v1_0*w1_5;\
	  __builtin_prefetch(vnext); vnext +=4;\
	  __builtin_prefetch(wnext); wnext +=4;\
	f10r += v1_2*w1_0+v1_3*w1_1;\
	f10i += v1_3*w1_0-v1_2*w1_1;\
	f11r += v1_2*w1_2+v1_3*w1_3;\
	f11i += v1_3*w1_2-v1_2*w1_3;\
	f12r += v1_2*w1_4+v1_3*w1_5;\
	f12i += v1_3*w1_4-v1_2*w1_5;\
	  __builtin_prefetch(vnext); vnext +=4;\
	  __builtin_prefetch(wnext); wnext +=4;\
	f20r += v1_4*w1_0+v1_5*w1_1;\
	f20i += v1_5*w1_0-v1_4*w1_1;\
	f21r += v1_4*w1_2+v1_5*w1_3;\
	f21i += v1_5*w1_2-v1_4*w1_3;\
	f22r += v1_4*w1_4+v1_5*w1_5;\
	f22i += v1_5*w1_4-v1_4*w1_5;\
	  __builtin_prefetch(vnext); vnext +=4;\
	  __builtin_prefetch(wnext); wnext +=4;

#define SPROJ_END() \ 
    v += 24 + v_stride ;\
    w += 24 + w_stride ;\
    vnext +=  v_stride ;\
    wnext +=  w_stride ;\
  }\
  tf = f;\
  *tf = f00r;tf++; *tf = f00i;tf++;\
  *tf = f01r;tf++; *tf = f01i;tf++;\
  *tf = f02r;tf++; *tf = f02i;tf++;\
  *tf = f10r;tf++; *tf = f10i;tf++;\
  *tf = f11r;tf++; *tf = f11i;tf++;\
  *tf = f12r;tf++; *tf = f12i;tf++;\
  *tf = f20r;tf++; *tf = f20i;tf++;\
  *tf = f21r;tf++; *tf = f21i;tf++;\
  *tf = f22r;tf++; *tf = f22i;tf++;\
  return ;

extern "C" {
//------------------------------------------------------------------
// sproj with (1 + gamma_0)
//------------------------------------------------------------------
void sprojTrXp(IFloat *f, IFloat *v, IFloat *w, int num_blk,
               int v_stride, int w_stride)
{
	SPROJ_START()

	v1_0 = (*(tv   ) - *(tv+19));
	v1_1 = (*(tv+ 1) + *(tv+18));
	v1_2 = (*(tv+ 2) - *(tv+21));
	v1_3 = (*(tv+ 3) + *(tv+20));
	v1_4 = (*(tv+ 4) - *(tv+23));
	v1_5 = (*(tv+ 5) + *(tv+22));
	w1_0 = (*(tw   ) - *(tw+19));
	w1_1 = (*(tw+ 1) + *(tw+18));
	w1_2 = (*(tw+ 2) - *(tw+21));
	w1_3 = (*(tw+ 3) + *(tw+20));
	w1_4 = (*(tw+ 4) - *(tw+23));
	w1_5 = (*(tw+ 5) + *(tw+22));

	SPROJ_COMP()

	v1_0 = (*(tv+ 6) - *(tv+13));
	v1_1 = (*(tv+ 7) + *(tv+12));
	v1_2 = (*(tv+ 8) - *(tv+15));
	v1_3 = (*(tv+ 9) + *(tv+14));
	v1_4 = (*(tv+10) - *(tv+17));
	v1_5 = (*(tv+11) + *(tv+16));
	w1_0 = (*(tw+ 6) - *(tw+13));
	w1_1 = (*(tw+ 7) + *(tw+12));
	w1_2 = (*(tw+ 8) - *(tw+15));
	w1_3 = (*(tw+ 9) + *(tw+14));
	w1_4 = (*(tw+10) - *(tw+17));
	w1_5 = (*(tw+11) + *(tw+16));

	SPROJ_COMP()
	SPROJ_END()
}


//------------------------------------------------------------------
// sproj with (1 - gamma_0)
//------------------------------------------------------------------
void sprojTrXm(IFloat *f, IFloat *v, IFloat *w, int num_blk,
               int v_stride, int w_stride)
{
	SPROJ_START()
	v1_0 = (*(tv   ) + *(tv+19));
	v1_1 = (*(tv+ 1) - *(tv+18));
	v1_2 = (*(tv+ 2) + *(tv+21));
	v1_3 = (*(tv+ 3) - *(tv+20));
	v1_4 = (*(tv+ 4) + *(tv+23));
	v1_5 = (*(tv+ 5) - *(tv+22));
	w1_0 = (*(tw   ) + *(tw+19));
	w1_1 = (*(tw+ 1) - *(tw+18));
	w1_2 = (*(tw+ 2) + *(tw+21));
	w1_3 = (*(tw+ 3) - *(tw+20));
	w1_4 = (*(tw+ 4) + *(tw+23));
	w1_5 = (*(tw+ 5) - *(tw+22));
	SPROJ_COMP()
	v1_0 = (*(tv+ 6) + *(tv+13));
	v1_1 = (*(tv+ 7) - *(tv+12));
	v1_2 = (*(tv+ 8) + *(tv+15));
	v1_3 = (*(tv+ 9) - *(tv+14));
	v1_4 = (*(tv+10) + *(tv+17));
	v1_5 = (*(tv+11) - *(tv+16));
	w1_0 = (*(tw+ 6) + *(tw+13));
	w1_1 = (*(tw+ 7) - *(tw+12));
	w1_2 = (*(tw+ 8) + *(tw+15));
	w1_3 = (*(tw+ 9) - *(tw+14));
	w1_4 = (*(tw+10) + *(tw+17));
	w1_5 = (*(tw+11) - *(tw+16));
	SPROJ_COMP()
	SPROJ_END()
}


//------------------------------------------------------------------
// sproj with (1 + gamma_1)
//------------------------------------------------------------------
void sprojTrYp(IFloat *f, IFloat *v, IFloat *w, int num_blk,
               int v_stride, int w_stride)
{
	SPROJ_START()
	v1_0 = (*(tv   ) - *(tv+18));
	v1_1 = (*(tv+ 1) - *(tv+19));
	v1_2 = (*(tv+ 2) - *(tv+20));
	v1_3 = (*(tv+ 3) - *(tv+21));
	v1_4 = (*(tv+ 4) - *(tv+22));
	v1_5 = (*(tv+ 5) - *(tv+23));
	w1_0 = (*(tw   ) - *(tw+18));
	w1_1 = (*(tw+ 1) - *(tw+19));
	w1_2 = (*(tw+ 2) - *(tw+20));
	w1_3 = (*(tw+ 3) - *(tw+21));
	w1_4 = (*(tw+ 4) - *(tw+22));
	w1_5 = (*(tw+ 5) - *(tw+23));
	SPROJ_COMP()
	v1_0 = (*(tv+ 6) + *(tv+12));
	v1_1 = (*(tv+ 7) + *(tv+13));
	v1_2 = (*(tv+ 8) + *(tv+14));
	v1_3 = (*(tv+ 9) + *(tv+15));
	v1_4 = (*(tv+10) + *(tv+16));
	v1_5 = (*(tv+11) + *(tv+17));
	w1_0 = (*(tw+ 6) + *(tw+12));
	w1_1 = (*(tw+ 7) + *(tw+13));
	w1_2 = (*(tw+ 8) + *(tw+14));
	w1_3 = (*(tw+ 9) + *(tw+15));
	w1_4 = (*(tw+10) + *(tw+16));
	w1_5 = (*(tw+11) + *(tw+17));
	SPROJ_COMP()
	SPROJ_END()
}

//------------------------------------------------------------------
// sproj with (1 - gamma_1)
//------------------------------------------------------------------
void sprojTrYm(IFloat *f, IFloat *v, IFloat *w, int num_blk,
               int v_stride, int w_stride)
{
	SPROJ_START()
	v1_0 = (*(tv   ) + *(tv+18));
	v1_1 = (*(tv+ 1) + *(tv+19));
	v1_2 = (*(tv+ 2) + *(tv+20));
	v1_3 = (*(tv+ 3) + *(tv+21));
	v1_4 = (*(tv+ 4) + *(tv+22));
	v1_5 = (*(tv+ 5) + *(tv+23));
	w1_0 = (*(tw   ) + *(tw+18));
	w1_1 = (*(tw+ 1) + *(tw+19));
	w1_2 = (*(tw+ 2) + *(tw+20));
	w1_3 = (*(tw+ 3) + *(tw+21));
	w1_4 = (*(tw+ 4) + *(tw+22));
	w1_5 = (*(tw+ 5) + *(tw+23));
	SPROJ_COMP()
	v1_0 = (*(tv+ 6) - *(tv+12));
	v1_1 = (*(tv+ 7) - *(tv+13));
	v1_2 = (*(tv+ 8) - *(tv+14));
	v1_3 = (*(tv+ 9) - *(tv+15));
	v1_4 = (*(tv+10) - *(tv+16));
	v1_5 = (*(tv+11) - *(tv+17));
	w1_0 = (*(tw+ 6) - *(tw+12));
	w1_1 = (*(tw+ 7) - *(tw+13));
	w1_2 = (*(tw+ 8) - *(tw+14));
	w1_3 = (*(tw+ 9) - *(tw+15));
	w1_4 = (*(tw+10) - *(tw+16));
	w1_5 = (*(tw+11) - *(tw+17));
	SPROJ_COMP()
	SPROJ_END()
}


//------------------------------------------------------------------
// sproj with (1 + gamma_2)
//------------------------------------------------------------------
void sprojTrZp(IFloat *f, IFloat *v, IFloat *w, int num_blk,
               int v_stride, int w_stride)
{
	SPROJ_START()
	v1_0 = (*(tv   ) - *(tv+13));
	v1_1 = (*(tv+ 1) + *(tv+12));
	v1_2 = (*(tv+ 2) - *(tv+15));
	v1_3 = (*(tv+ 3) + *(tv+14));
	v1_4 = (*(tv+ 4) - *(tv+17));
	v1_5 = (*(tv+ 5) + *(tv+16));
	w1_0 = (*(tw   ) - *(tw+13));
	w1_1 = (*(tw+ 1) + *(tw+12));
	w1_2 = (*(tw+ 2) - *(tw+15));
	w1_3 = (*(tw+ 3) + *(tw+14));
	w1_4 = (*(tw+ 4) - *(tw+17));
	w1_5 = (*(tw+ 5) + *(tw+16));
	SPROJ_COMP()
	v1_0 = (*(tv+ 6) + *(tv+19));
	v1_1 = (*(tv+ 7) - *(tv+18));
	v1_2 = (*(tv+ 8) + *(tv+21));
	v1_3 = (*(tv+ 9) - *(tv+20));
	v1_4 = (*(tv+10) + *(tv+23));
	v1_5 = (*(tv+11) - *(tv+22));
	w1_0 = (*(tw+ 6) + *(tw+19));
	w1_1 = (*(tw+ 7) - *(tw+18));
	w1_2 = (*(tw+ 8) + *(tw+21));
	w1_3 = (*(tw+ 9) - *(tw+20));
	w1_4 = (*(tw+10) + *(tw+23));
	w1_5 = (*(tw+11) - *(tw+22));
	SPROJ_COMP()
	SPROJ_END()
}


//------------------------------------------------------------------
// sproj with (1 - gamma_2)
//------------------------------------------------------------------
void sprojTrZm(IFloat *f, IFloat *v, IFloat *w, int num_blk,
               int v_stride, int w_stride)
{
	SPROJ_START()
	v1_0 = (*(tv   ) + *(tv+13));
	v1_1 = (*(tv+ 1) - *(tv+12));
	v1_2 = (*(tv+ 2) + *(tv+15));
	v1_3 = (*(tv+ 3) - *(tv+14));
	v1_4 = (*(tv+ 4) + *(tv+17));
	v1_5 = (*(tv+ 5) - *(tv+16));
	w1_0 = (*(tw   ) + *(tw+13));
	w1_1 = (*(tw+ 1) - *(tw+12));
	w1_2 = (*(tw+ 2) + *(tw+15));
	w1_3 = (*(tw+ 3) - *(tw+14));
	w1_4 = (*(tw+ 4) + *(tw+17));
	w1_5 = (*(tw+ 5) - *(tw+16));
	SPROJ_COMP()
	v1_0 = (*(tv+ 6) - *(tv+19));
	v1_1 = (*(tv+ 7) + *(tv+18));
	v1_2 = (*(tv+ 8) - *(tv+21));
	v1_3 = (*(tv+ 9) + *(tv+20));
	v1_4 = (*(tv+10) - *(tv+23));
	v1_5 = (*(tv+11) + *(tv+22));
	w1_0 = (*(tw+ 6) - *(tw+19));
	w1_1 = (*(tw+ 7) + *(tw+18));
	w1_2 = (*(tw+ 8) - *(tw+21));
	w1_3 = (*(tw+ 9) + *(tw+20));
	w1_4 = (*(tw+10) - *(tw+23));
	w1_5 = (*(tw+11) + *(tw+22));
	SPROJ_COMP()
	SPROJ_END()
}


//------------------------------------------------------------------
// sproj with (1 + gamma_3)
//------------------------------------------------------------------
void sprojTrTp(IFloat *f, IFloat *v, IFloat *w, int num_blk,
               int v_stride, int w_stride)
{
	SPROJ_START()
	v1_0 = (*(tv   ) + *(tv+12));
	v1_1 = (*(tv+ 1) + *(tv+13));
	v1_2 = (*(tv+ 2) + *(tv+14));
	v1_3 = (*(tv+ 3) + *(tv+15));
	v1_4 = (*(tv+ 4) + *(tv+16));
	v1_5 = (*(tv+ 5) + *(tv+17));
	w1_0 = (*(tw   ) + *(tw+12));
	w1_1 = (*(tw+ 1) + *(tw+13));
	w1_2 = (*(tw+ 2) + *(tw+14));
	w1_3 = (*(tw+ 3) + *(tw+15));
	w1_4 = (*(tw+ 4) + *(tw+16));
	w1_5 = (*(tw+ 5) + *(tw+17));
	SPROJ_COMP()
	v1_0 = (*(tv+ 6) + *(tv+18));
	v1_1 = (*(tv+ 7) + *(tv+19));
	v1_2 = (*(tv+ 8) + *(tv+20));
	v1_3 = (*(tv+ 9) + *(tv+21));
	v1_4 = (*(tv+10) + *(tv+22));
	v1_5 = (*(tv+11) + *(tv+23));
	w1_0 = (*(tw+ 6) + *(tw+18));
	w1_1 = (*(tw+ 7) + *(tw+19));
	w1_2 = (*(tw+ 8) + *(tw+20));
	w1_3 = (*(tw+ 9) + *(tw+21));
	w1_4 = (*(tw+10) + *(tw+22));
	w1_5 = (*(tw+11) + *(tw+23));
	SPROJ_COMP()
	SPROJ_END()
}


//------------------------------------------------------------------
// sproj with (1 - gamma_3)
//------------------------------------------------------------------
void sprojTrTm(IFloat *f, IFloat *v, IFloat *w, int num_blk,
               int v_stride, int w_stride)
#if 1
{
	SPROJ_START()
	v1_0 = (*(tv   ) - *(tv+12));
	v1_1 = (*(tv+ 1) - *(tv+13));
	v1_2 = (*(tv+ 2) - *(tv+14));
	v1_3 = (*(tv+ 3) - *(tv+15));
	v1_4 = (*(tv+ 4) - *(tv+16));
	v1_5 = (*(tv+ 5) - *(tv+17));
	w1_0 = (*(tw   ) - *(tw+12));
	w1_1 = (*(tw+ 1) - *(tw+13));
	w1_2 = (*(tw+ 2) - *(tw+14));
	w1_3 = (*(tw+ 3) - *(tw+15));
	w1_4 = (*(tw+ 4) - *(tw+16));
	w1_5 = (*(tw+ 5) - *(tw+17));
	SPROJ_COMP()
	v1_0 = (*(tv+ 6) - *(tv+18));
	v1_1 = (*(tv+ 7) - *(tv+19));
	v1_2 = (*(tv+ 8) - *(tv+20));
	v1_3 = (*(tv+ 9) - *(tv+21));
	v1_4 = (*(tv+10) - *(tv+22));
	v1_5 = (*(tv+11) - *(tv+23));
	w1_0 = (*(tw+ 6) - *(tw+18));
	w1_1 = (*(tw+ 7) - *(tw+19));
	w1_2 = (*(tw+ 8) - *(tw+20));
	w1_3 = (*(tw+ 9) - *(tw+21));
	w1_4 = (*(tw+10) - *(tw+22));
	w1_5 = (*(tw+11) - *(tw+23));
	SPROJ_COMP()
	SPROJ_END()
}
#else
{
  int row, col, blk ;

  IFloat *tf, *tv, *tw ;

  tf = f ;

  for (row=0; row<MATRIX_SIZE; row++) 
    *tf++ = 0.0 ;

  for (blk=0; blk<num_blk; blk++) {
    tf = f ;
    tv = v ;
    tw = w ;
    for (row=0; row<3; row++) {
      for (col=0; col<3; col++) {
        *tf++ +=   (*(tv   ) - *(tv+12)) * *(tw   )
                 + (*(tv+ 1) - *(tv+13)) * *(tw+ 1)
                 + (*(tv+ 6) - *(tv+18)) * *(tw+ 6)
                 + (*(tv+ 7) - *(tv+19)) * *(tw+ 7)
                 + (*(tv+12) - *(tv   )) * *(tw+12)
                 + (*(tv+13) - *(tv+ 1)) * *(tw+13)
                 + (*(tv+18) - *(tv+ 6)) * *(tw+18)
                 + (*(tv+19) - *(tv+ 7)) * *(tw+19) ;

        *tf++ +=   (*(tv+ 1) - *(tv+13)) * *(tw   )
                 - (*(tv   ) - *(tv+12)) * *(tw+ 1)
                 + (*(tv+ 7) - *(tv+19)) * *(tw+ 6)
                 - (*(tv+ 6) - *(tv+18)) * *(tw+ 7)
                 + (*(tv+13) - *(tv+ 1)) * *(tw+12)
                 - (*(tv+12) - *(tv   )) * *(tw+13)
                 + (*(tv+19) - *(tv+ 7)) * *(tw+18)
                 - (*(tv+18) - *(tv+ 6)) * *(tw+19) ;

        tw += 2 ;
      }
      tw = w ;
      tv += 2 ;
    }
    v += 24 + v_stride ;
    w += 24 + w_stride ;
  }
  return ;
}
#endif
}

CPS_END_NAMESPACE
