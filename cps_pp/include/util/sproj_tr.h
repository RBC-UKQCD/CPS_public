#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Declaration of spin projection and trace routines.

  $Id: sproj_tr.h,v 1.5 2004-11-21 21:40:00 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-11-21 21:40:00 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/sproj_tr.h,v 1.5 2004-11-21 21:40:00 chulwoo Exp $
//  $Id: sproj_tr.h,v 1.5 2004-11-21 21:40:00 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: sproj_tr.h,v $
//  $Revision: 1.5 $
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
extern "C"{
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

}



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
