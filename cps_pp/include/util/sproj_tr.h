#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Declaration of spin projection and trace routines.

  $Id: sproj_tr.h,v 1.2 2003-07-24 16:53:53 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/sproj_tr.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: sproj_tr.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:31  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:19  anj
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
//  $RCSfile: sproj_tr.h,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/sproj_tr.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// sproj_tr.h
//
// These functions return a color matrix in f constructed from
// the spinors v, w using: 
//
// f_(i,j) = Tr_spin[ (1 +/- gamma_mu) v_i w^dag_j ]
// for the sprojTr* functions and
//
// f_(i,j) = 1/2 Tr_spin[ Sigma_{mu,nu} v_i w^dag_j ]
// for the SigmaprojTr* functions.
//
// num_blk is the number of spinors v, w. The routines 
// accumulate the sum over spinors in f.
//
// stride is the number of IFloats between spinors (a stride = 0
// means that the spinors are consecutive in memory)
//
//------------------------------------------------------------------

#ifndef INCLUDED_SPROJ_TR_H
#define INCLUDED_SPROJ_TR_H

CPS_END_NAMESPACE
#include <util/data_types.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// sprojTr* functions
//------------------------------------------------------------------

void sprojTrXm(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with (1 - gamma_0)

void sprojTrYm(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with (1 - gamma_1)

void sprojTrZm(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with (1 - gamma_2)

void sprojTrTm(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with (1 - gamma_3)

void sprojTrXp(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with (1 + gamma_0)

void sprojTrYp(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with (1 + gamma_1)

void sprojTrZp(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with (1 + gamma_2)

void sprojTrTp(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with (1 + gamma_3)



//------------------------------------------------------------------
// Sigmaproj* functions
//------------------------------------------------------------------

void SigmaprojTrXY(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with Sigma_{0,1}

void SigmaprojTrYX(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with Sigma_{1,0}

void SigmaprojTrXZ(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with Sigma_{0,2}

void SigmaprojTrZX(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with Sigma_{2,0}

void SigmaprojTrXT(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with Sigma_{0,3}

void SigmaprojTrTX(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with Sigma_{3,0}

void SigmaprojTrYZ(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with Sigma_{1,2}

void SigmaprojTrZY(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with Sigma_{2,1}

void SigmaprojTrYT(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with Sigma_{1,3}

void SigmaprojTrTY(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with Sigma_{3,1}

void SigmaprojTrZT(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with Sigma_{2,3}

void SigmaprojTrTZ(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with Sigma_{3,2}

void SigmaprojTrXX(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with Sigma_{0,0}

void SigmaprojTrYY(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with Sigma_{1,1}

void SigmaprojTrZZ(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with Sigma_{2,2}

void SigmaprojTrTT(IFloat *f, IFloat *v_plus_mu, IFloat *w, int num_blk,
               int v_plus_mu_stride, int w_stride) ;
//!< Projection with Sigma_{3,3}

#endif

CPS_END_NAMESPACE
