#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_wilson_types/comsrc/sigmaproj_tr.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: sigmaproj_tr.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:36  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:25  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: sigmaproj_tr.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_wilson_types/comsrc/sigmaproj_tr.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
///------------------------------------------------------------------
//
// Sigmaproj_tr.C
//
// These functions return a color matrix in f constructed from
// the spinors v, w using: 
// f_(i,j) = 1/2 Tr_spin[ Sigma_{mu,nu} v_i w^dag_j ]
//
// num_blk is the number of spinors v, w. The routines 
// accumulate the sum over spinors in f.
//
// stride is the number of IFloats between spinors (a stride = 0
// means that the spinors are consecutive in memory)
// 
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/sproj_tr.h>
CPS_START_NAMESPACE

enum{MATRIX_SIZE = 18};


//------------------------------------------------------------------
// Sigmaproj with Sigma_{0,1}
//------------------------------------------------------------------
void SigmaprojTrXY(IFloat *f, IFloat *v, IFloat *w, int num_blk,
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
        *tf++ +=   (*(tv+ 1)) * *(tw   )
                 - (*(tv   )) * *(tw+ 1)
                 - (*(tv+ 7)) * *(tw+ 6)
                 + (*(tv+ 6)) * *(tw+ 7)
                 + (*(tv+13)) * *(tw+12)
                 - (*(tv+12)) * *(tw+13)
                 - (*(tv+19)) * *(tw+18)
                 + (*(tv+18)) * *(tw+19) ;

        *tf++ += - (*(tv   )) * *(tw   )
                 - (*(tv+ 1)) * *(tw+ 1)
                 + (*(tv+ 6)) * *(tw+ 6)
                 + (*(tv+ 7)) * *(tw+ 7)
                 - (*(tv+12)) * *(tw+12)
                 - (*(tv+13)) * *(tw+13)
                 + (*(tv+18)) * *(tw+18)
                 + (*(tv+19)) * *(tw+19) ;

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
// Sigmaproj with Sigma_{1,0}
//------------------------------------------------------------------
void SigmaprojTrYX(IFloat *f, IFloat *v, IFloat *w, int num_blk,
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
        *tf++ += - (*(tv+ 1)) * *(tw   )
                 + (*(tv   )) * *(tw+ 1)
                 + (*(tv+ 7)) * *(tw+ 6)
                 - (*(tv+ 6)) * *(tw+ 7)
                 - (*(tv+13)) * *(tw+12)
                 + (*(tv+12)) * *(tw+13)
                 + (*(tv+19)) * *(tw+18)
                 - (*(tv+18)) * *(tw+19) ;

        *tf++ +=   (*(tv   )) * *(tw   )
                 + (*(tv+ 1)) * *(tw+ 1)
                 - (*(tv+ 6)) * *(tw+ 6)
                 - (*(tv+ 7)) * *(tw+ 7)
                 + (*(tv+12)) * *(tw+12)
                 + (*(tv+13)) * *(tw+13)
                 - (*(tv+18)) * *(tw+18)
                 - (*(tv+19)) * *(tw+19) ;

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
// Sigmaproj with Sigma_{0,2}
//------------------------------------------------------------------
void SigmaprojTrXZ(IFloat *f, IFloat *v, IFloat *w, int num_blk,
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

	*tf++ += - (*(tv+ 6)) * *(tw   )
                 - (*(tv+ 7)) * *(tw+ 1)
                 + (*(tv   )) * *(tw+ 6)
                 + (*(tv+ 1)) * *(tw+ 7)
                 - (*(tv+18)) * *(tw+12)
                 - (*(tv+19)) * *(tw+13)
                 + (*(tv+12)) * *(tw+18)
                 + (*(tv+13)) * *(tw+19) ;

	*tf++ += -  (*(tv+ 7)) * *(tw   )
                 + (*(tv+ 6)) * *(tw+ 1)
                 + (*(tv+ 1)) * *(tw+ 6)
                 - (*(tv   )) * *(tw+ 7)
                 - (*(tv+19)) * *(tw+12)
                 + (*(tv+18)) * *(tw+13)
                 + (*(tv+13)) * *(tw+18)
                 - (*(tv+12)) * *(tw+19) ;

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
// Sigmaproj with Sigma_{2,0}
//------------------------------------------------------------------
void SigmaprojTrZX(IFloat *f, IFloat *v, IFloat *w, int num_blk,
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

        *tf++ +=   (*(tv+ 6)) * *(tw   )
                 + (*(tv+ 7)) * *(tw+ 1)
                 - (*(tv   )) * *(tw+ 6)
                 - (*(tv+ 1)) * *(tw+ 7)
                 + (*(tv+18)) * *(tw+12)
                 + (*(tv+19)) * *(tw+13)
                 - (*(tv+12)) * *(tw+18)
                 - (*(tv+13)) * *(tw+19) ;


        *tf++ +=   (*(tv+ 7)) * *(tw   )
                 - (*(tv+ 6)) * *(tw+ 1)
                 - (*(tv+ 1)) * *(tw+ 6)
                 + (*(tv   )) * *(tw+ 7)
                 + (*(tv+19)) * *(tw+12)
                 - (*(tv+18)) * *(tw+13)
                 - (*(tv+13)) * *(tw+18)
                 + (*(tv+12)) * *(tw+19) ;

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
// Sigmaproj with Sigma_{0,3}
//------------------------------------------------------------------
void SigmaprojTrXT(IFloat *f, IFloat *v, IFloat *w, int num_blk,
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

	*tf++ += - (*(tv+ 7)) * *(tw   )
                 + (*(tv+ 6)) * *(tw+ 1)
                 - (*(tv+ 1)) * *(tw+ 6)
                 + (*(tv   )) * *(tw+ 7)
                 + (*(tv+19)) * *(tw+12)
                 - (*(tv+18)) * *(tw+13)
                 + (*(tv+13)) * *(tw+18)
                 - (*(tv+12)) * *(tw+19) ;
	
 	*tf++ +=   (*(tv+ 6)) * *(tw   )
                 + (*(tv+ 7)) * *(tw+ 1)
                 + (*(tv   )) * *(tw+ 6)
                 + (*(tv+ 1)) * *(tw+ 7)
                 - (*(tv+18)) * *(tw+12)
                 - (*(tv+19)) * *(tw+13)
                 - (*(tv+12)) * *(tw+18)
                 - (*(tv+13)) * *(tw+19) ;
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
// Sigmaproj with Sigma_{3,0}
//------------------------------------------------------------------
void SigmaprojTrTX(IFloat *f, IFloat *v, IFloat *w, int num_blk,
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

        *tf++ +=   (*(tv+ 7)) * *(tw   )
                 - (*(tv+ 6)) * *(tw+ 1)
                 + (*(tv+ 1)) * *(tw+ 6)
                 - (*(tv   )) * *(tw+ 7)
                 - (*(tv+19)) * *(tw+12)
                 + (*(tv+18)) * *(tw+13)
                 - (*(tv+13)) * *(tw+18)
                 + (*(tv+12)) * *(tw+19) ;

        *tf++ += - (*(tv+ 6)) * *(tw   )
                 - (*(tv+ 7)) * *(tw+ 1)
                 - (*(tv   )) * *(tw+ 6)
                 - (*(tv+ 1)) * *(tw+ 7)
                 + (*(tv+18)) * *(tw+12)
                 + (*(tv+19)) * *(tw+13)
                 + (*(tv+12)) * *(tw+18)
                 + (*(tv+13)) * *(tw+19) ;


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
// Sigmaproj with Sigma_{1,2}
//------------------------------------------------------------------
void SigmaprojTrYZ(IFloat *f, IFloat *v, IFloat *w, int num_blk,
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

	*tf++ +=   (*(tv+ 7)) * *(tw   )
                 - (*(tv+ 6)) * *(tw+ 1)
                 + (*(tv+ 1)) * *(tw+ 6)
                 - (*(tv   )) * *(tw+ 7)
                 + (*(tv+19)) * *(tw+12)
                 - (*(tv+18)) * *(tw+13)
                 + (*(tv+13)) * *(tw+18)
                 - (*(tv+12)) * *(tw+19) ;
	
 	*tf++ += - (*(tv+ 6)) * *(tw   )
                 - (*(tv+ 7)) * *(tw+ 1)
                 - (*(tv   )) * *(tw+ 6)
                 - (*(tv+ 1)) * *(tw+ 7)
                 - (*(tv+18)) * *(tw+12)
                 - (*(tv+19)) * *(tw+13)
                 - (*(tv+12)) * *(tw+18)
                 - (*(tv+13)) * *(tw+19) ;
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
// Sigmaproj with Sigma_{2,1}
//------------------------------------------------------------------
void SigmaprojTrZY(IFloat *f, IFloat *v, IFloat *w, int num_blk,
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

        *tf++ += - (*(tv+ 7)) * *(tw   )
                 + (*(tv+ 6)) * *(tw+ 1)
                 - (*(tv+ 1)) * *(tw+ 6)
                 + (*(tv   )) * *(tw+ 7)
                 - (*(tv+19)) * *(tw+12)
                 + (*(tv+18)) * *(tw+13)
                 - (*(tv+13)) * *(tw+18)
                 + (*(tv+12)) * *(tw+19) ;

        *tf++ +=   (*(tv+ 6)) * *(tw   )
                 + (*(tv+ 7)) * *(tw+ 1)
                 + (*(tv   )) * *(tw+ 6)
                 + (*(tv+ 1)) * *(tw+ 7)
                 + (*(tv+18)) * *(tw+12)
                 + (*(tv+19)) * *(tw+13)
                 + (*(tv+12)) * *(tw+18)
                 + (*(tv+13)) * *(tw+19) ;


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
// Sigmaproj with Sigma_{1,3}
//------------------------------------------------------------------
void SigmaprojTrYT(IFloat *f, IFloat *v, IFloat *w, int num_blk,
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

	*tf++ += - (*(tv+ 6)) * *(tw   )
                 - (*(tv+ 7)) * *(tw+ 1)
                 + (*(tv   )) * *(tw+ 6)
                 + (*(tv+ 1)) * *(tw+ 7)
                 + (*(tv+18)) * *(tw+12)
                 + (*(tv+19)) * *(tw+13)
                 - (*(tv+12)) * *(tw+18)
                 - (*(tv+13)) * *(tw+19) ;

	*tf++ += -  (*(tv+ 7)) * *(tw   )
                 + (*(tv+ 6)) * *(tw+ 1)
                 + (*(tv+ 1)) * *(tw+ 6)
                 - (*(tv   )) * *(tw+ 7)
                 + (*(tv+19)) * *(tw+12)
                 - (*(tv+18)) * *(tw+13)
                 - (*(tv+13)) * *(tw+18)
                 + (*(tv+12)) * *(tw+19) ;

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
// Sigmaproj with Sigma_{3,1}
//------------------------------------------------------------------
void SigmaprojTrTY(IFloat *f, IFloat *v, IFloat *w, int num_blk,
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

        *tf++ +=   (*(tv+ 6)) * *(tw   )
                 + (*(tv+ 7)) * *(tw+ 1)
                 - (*(tv   )) * *(tw+ 6)
                 - (*(tv+ 1)) * *(tw+ 7)
                 - (*(tv+18)) * *(tw+12)
                 - (*(tv+19)) * *(tw+13)
                 + (*(tv+12)) * *(tw+18)
                 + (*(tv+13)) * *(tw+19) ;


        *tf++ +=   (*(tv+ 7)) * *(tw   )
                 - (*(tv+ 6)) * *(tw+ 1)
                 - (*(tv+ 1)) * *(tw+ 6)
                 + (*(tv   )) * *(tw+ 7)
                 - (*(tv+19)) * *(tw+12)
                 + (*(tv+18)) * *(tw+13)
                 + (*(tv+13)) * *(tw+18)
                 - (*(tv+12)) * *(tw+19) ;

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
// Sigmaproj with Sigma_{2,3}
//------------------------------------------------------------------
void SigmaprojTrZT(IFloat *f, IFloat *v, IFloat *w, int num_blk,
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
        *tf++ += - (*(tv+ 1)) * *(tw   )
                 + (*(tv   )) * *(tw+ 1)
                 + (*(tv+ 7)) * *(tw+ 6)
                 - (*(tv+ 6)) * *(tw+ 7)
                 + (*(tv+13)) * *(tw+12)
                 - (*(tv+12)) * *(tw+13)
                 - (*(tv+19)) * *(tw+18)
                 + (*(tv+18)) * *(tw+19) ;

        *tf++ +=   (*(tv   )) * *(tw   )
                 + (*(tv+ 1)) * *(tw+ 1)
                 - (*(tv+ 6)) * *(tw+ 6)
                 - (*(tv+ 7)) * *(tw+ 7)
                 - (*(tv+12)) * *(tw+12)
                 - (*(tv+13)) * *(tw+13)
                 + (*(tv+18)) * *(tw+18)
                 + (*(tv+19)) * *(tw+19) ;

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
// Sigmaproj with Sigma_{3,2}
//------------------------------------------------------------------
void SigmaprojTrTZ(IFloat *f, IFloat *v, IFloat *w, int num_blk,
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
        *tf++ +=   (*(tv+ 1)) * *(tw   )
                 - (*(tv   )) * *(tw+ 1)
                 - (*(tv+ 7)) * *(tw+ 6)
                 + (*(tv+ 6)) * *(tw+ 7)
                 - (*(tv+13)) * *(tw+12)
                 + (*(tv+12)) * *(tw+13)
                 + (*(tv+19)) * *(tw+18)
                 - (*(tv+18)) * *(tw+19) ;

        *tf++ += - (*(tv   )) * *(tw   )
                 - (*(tv+ 1)) * *(tw+ 1)
                 + (*(tv+ 6)) * *(tw+ 6)
                 + (*(tv+ 7)) * *(tw+ 7)
                 + (*(tv+12)) * *(tw+12)
                 + (*(tv+13)) * *(tw+13)
                 - (*(tv+18)) * *(tw+18)
                 - (*(tv+19)) * *(tw+19) ;

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
