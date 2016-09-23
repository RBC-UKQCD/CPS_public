#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  gamma matrix projection and spin trace routines.

  Used by derivatives of the FwilsonTypes class.
  
  $Id: sproj_tr.C,v 1.4 2004/08/18 11:58:04 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/08/18 11:58:04 $
//  $Header: /space/cvs/cps/cps++/src/util/lattice/f_wilson_types/noarch/sproj_tr.C,v 1.4 2004/08/18 11:58:04 zs Exp $
//  $Id: sproj_tr.C,v 1.4 2004/08/18 11:58:04 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /space/cvs/cps/cps++/src/util/lattice/f_wilson_types/noarch/sproj_tr.C,v $
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

enum{MATRIX_SIZE = 18};


//------------------------------------------------------------------
// sproj with (1 + gamma_0)
//------------------------------------------------------------------
void sprojTrXp(IFloat *f, IFloat *v, IFloat *w, int num_blk,
               int v_stride, int w_stride)
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
        *tf++ +=   (*(tv   ) - *(tv+19)) * *(tw   )
                 + (*(tv+ 1) + *(tv+18)) * *(tw+ 1)
                 + (*(tv+ 6) - *(tv+13)) * *(tw+ 6)
                 + (*(tv+ 7) + *(tv+12)) * *(tw+ 7)
                 + (*(tv+12) + *(tv+ 7)) * *(tw+12)
                 + (*(tv+13) - *(tv+ 6)) * *(tw+13)
                 + (*(tv+18) + *(tv+ 1)) * *(tw+18)
                 + (*(tv+19) - *(tv   )) * *(tw+19) ;

        *tf++ +=   (*(tv+ 1) + *(tv+18)) * *(tw   )
                 - (*(tv   ) - *(tv+19)) * *(tw+ 1)
                 + (*(tv+ 7) + *(tv+12)) * *(tw+ 6)
                 - (*(tv+ 6) - *(tv+13)) * *(tw+ 7)
                 + (*(tv+13) - *(tv+ 6)) * *(tw+12)
                 - (*(tv+12) + *(tv+ 7)) * *(tw+13)
                 + (*(tv+19) - *(tv   )) * *(tw+18)
                 - (*(tv+18) + *(tv+ 1)) * *(tw+19) ;

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


//------------------------------------------------------------------
// sproj with (1 - gamma_0)
//------------------------------------------------------------------
void sprojTrXm(IFloat *f, IFloat *v, IFloat *w, int num_blk,
               int v_stride, int w_stride)
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
        *tf++ +=   (*(tv   ) + *(tv+19)) * *(tw   )
                 + (*(tv+ 1) - *(tv+18)) * *(tw+ 1)
                 + (*(tv+ 6) + *(tv+13)) * *(tw+ 6)
                 + (*(tv+ 7) - *(tv+12)) * *(tw+ 7)
                 + (*(tv+12) - *(tv+ 7)) * *(tw+12)
                 + (*(tv+13) + *(tv+ 6)) * *(tw+13)
                 + (*(tv+18) - *(tv+ 1)) * *(tw+18)
                 + (*(tv+19) + *(tv   )) * *(tw+19) ;

        *tf++ +=   (*(tv+ 1) - *(tv+18)) * *(tw   )
                 - (*(tv   ) + *(tv+19)) * *(tw+ 1)
                 + (*(tv+ 7) - *(tv+12)) * *(tw+ 6)
                 - (*(tv+ 6) + *(tv+13)) * *(tw+ 7)
                 + (*(tv+13) + *(tv+ 6)) * *(tw+12)
                 - (*(tv+12) - *(tv+ 7)) * *(tw+13)
                 + (*(tv+19) + *(tv   )) * *(tw+18)
                 - (*(tv+18) - *(tv+ 1)) * *(tw+19) ;

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


//------------------------------------------------------------------
// sproj with (1 + gamma_1)
//------------------------------------------------------------------
void sprojTrYp(IFloat *f, IFloat *v, IFloat *w, int num_blk,
               int v_stride, int w_stride)
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
        *tf++ +=   (*(tv   ) - *(tv+18)) * *(tw   )
                 + (*(tv+ 1) - *(tv+19)) * *(tw+ 1)
                 + (*(tv+ 6) + *(tv+12)) * *(tw+ 6)
                 + (*(tv+ 7) + *(tv+13)) * *(tw+ 7)
                 + (*(tv+12) + *(tv+ 6)) * *(tw+12)
                 + (*(tv+13) + *(tv+ 7)) * *(tw+13)
                 + (*(tv+18) - *(tv   )) * *(tw+18)
                 + (*(tv+19) - *(tv+ 1)) * *(tw+19) ;

        *tf++ +=   (*(tv+ 1) - *(tv+19)) * *(tw   )
                 - (*(tv   ) - *(tv+18)) * *(tw+ 1)
                 + (*(tv+ 7) + *(tv+13)) * *(tw+ 6)
                 - (*(tv+ 6) + *(tv+12)) * *(tw+ 7)
                 + (*(tv+13) + *(tv+ 7)) * *(tw+12)
                 - (*(tv+12) + *(tv+ 6)) * *(tw+13)
                 + (*(tv+19) - *(tv+ 1)) * *(tw+18)
                 - (*(tv+18) - *(tv   )) * *(tw+19) ;

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


//------------------------------------------------------------------
// sproj with (1 - gamma_1)
//------------------------------------------------------------------
void sprojTrYm(IFloat *f, IFloat *v, IFloat *w, int num_blk,
               int v_stride, int w_stride)
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
        *tf++ +=   (*(tv   ) + *(tv+18)) * *(tw   )
                 + (*(tv+ 1) + *(tv+19)) * *(tw+ 1)
                 + (*(tv+ 6) - *(tv+12)) * *(tw+ 6)
                 + (*(tv+ 7) - *(tv+13)) * *(tw+ 7)
                 + (*(tv+12) - *(tv+ 6)) * *(tw+12)
                 + (*(tv+13) - *(tv+ 7)) * *(tw+13)
                 + (*(tv+18) + *(tv   )) * *(tw+18)
                 + (*(tv+19) + *(tv+ 1)) * *(tw+19) ;

        *tf++ +=   (*(tv+ 1) + *(tv+19)) * *(tw   )
                 - (*(tv   ) + *(tv+18)) * *(tw+ 1)
                 + (*(tv+ 7) - *(tv+13)) * *(tw+ 6)
                 - (*(tv+ 6) - *(tv+12)) * *(tw+ 7)
                 + (*(tv+13) - *(tv+ 7)) * *(tw+12)
                 - (*(tv+12) - *(tv+ 6)) * *(tw+13)
                 + (*(tv+19) + *(tv+ 1)) * *(tw+18)
                 - (*(tv+18) + *(tv   )) * *(tw+19) ;

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


//------------------------------------------------------------------
// sproj with (1 + gamma_2)
//------------------------------------------------------------------
void sprojTrZp(IFloat *f, IFloat *v, IFloat *w, int num_blk,
               int v_stride, int w_stride)
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
        *tf++ +=   (*(tv   ) - *(tv+13)) * *(tw   )
                 + (*(tv+ 1) + *(tv+12)) * *(tw+ 1)
                 + (*(tv+ 6) + *(tv+19)) * *(tw+ 6)
                 + (*(tv+ 7) - *(tv+18)) * *(tw+ 7)
                 + (*(tv+12) + *(tv+ 1)) * *(tw+12)
                 + (*(tv+13) - *(tv   )) * *(tw+13)
                 + (*(tv+18) - *(tv+ 7)) * *(tw+18)
                 + (*(tv+19) + *(tv+ 6)) * *(tw+19) ;

        *tf++ +=   (*(tv+ 1) + *(tv+12)) * *(tw   )
                 - (*(tv   ) - *(tv+13)) * *(tw+ 1)
                 + (*(tv+ 7) - *(tv+18)) * *(tw+ 6)
                 - (*(tv+ 6) + *(tv+19)) * *(tw+ 7)
                 + (*(tv+13) - *(tv   )) * *(tw+12)
                 - (*(tv+12) + *(tv+ 1)) * *(tw+13)
                 + (*(tv+19) + *(tv+ 6)) * *(tw+18)
                 - (*(tv+18) - *(tv+ 7)) * *(tw+19) ;

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


//------------------------------------------------------------------
// sproj with (1 - gamma_2)
//------------------------------------------------------------------
void sprojTrZm(IFloat *f, IFloat *v, IFloat *w, int num_blk,
               int v_stride, int w_stride)
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
        *tf++ +=   (*(tv   ) + *(tv+13)) * *(tw   )
                 + (*(tv+ 1) - *(tv+12)) * *(tw+ 1)
                 + (*(tv+ 6) - *(tv+19)) * *(tw+ 6)
                 + (*(tv+ 7) + *(tv+18)) * *(tw+ 7)
                 + (*(tv+12) - *(tv+ 1)) * *(tw+12)
                 + (*(tv+13) + *(tv   )) * *(tw+13)
                 + (*(tv+18) + *(tv+ 7)) * *(tw+18)
                 + (*(tv+19) - *(tv+ 6)) * *(tw+19) ;

        *tf++ +=   (*(tv+ 1) - *(tv+12)) * *(tw   )
                 - (*(tv   ) + *(tv+13)) * *(tw+ 1)
                 + (*(tv+ 7) + *(tv+18)) * *(tw+ 6)
                 - (*(tv+ 6) - *(tv+19)) * *(tw+ 7)
                 + (*(tv+13) + *(tv   )) * *(tw+12)
                 - (*(tv+12) - *(tv+ 1)) * *(tw+13)
                 + (*(tv+19) - *(tv+ 6)) * *(tw+18)
                 - (*(tv+18) + *(tv+ 7)) * *(tw+19) ;

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


//------------------------------------------------------------------
// sproj with (1 + gamma_3)
//------------------------------------------------------------------
void sprojTrTp(IFloat *f, IFloat *v, IFloat *w, int num_blk,
               int v_stride, int w_stride)
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
        *tf++ +=   (*(tv   ) + *(tv+12)) * *(tw   )
                 + (*(tv+ 1) + *(tv+13)) * *(tw+ 1)
                 + (*(tv+ 6) + *(tv+18)) * *(tw+ 6)
                 + (*(tv+ 7) + *(tv+19)) * *(tw+ 7)
                 + (*(tv+12) + *(tv   )) * *(tw+12)
                 + (*(tv+13) + *(tv+ 1)) * *(tw+13)
                 + (*(tv+18) + *(tv+ 6)) * *(tw+18)
                 + (*(tv+19) + *(tv+ 7)) * *(tw+19) ;

        *tf++ +=   (*(tv+ 1) + *(tv+13)) * *(tw   )
                 - (*(tv   ) + *(tv+12)) * *(tw+ 1)
                 + (*(tv+ 7) + *(tv+19)) * *(tw+ 6)
                 - (*(tv+ 6) + *(tv+18)) * *(tw+ 7)
                 + (*(tv+13) + *(tv+ 1)) * *(tw+12)
                 - (*(tv+12) + *(tv   )) * *(tw+13)
                 + (*(tv+19) + *(tv+ 7)) * *(tw+18)
                 - (*(tv+18) + *(tv+ 6)) * *(tw+19) ;

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


//------------------------------------------------------------------
// sproj with (1 - gamma_3)
//------------------------------------------------------------------
void sprojTrTm(IFloat *f, IFloat *v, IFloat *w, int num_blk,
               int v_stride, int w_stride)
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

CPS_END_NAMESPACE
