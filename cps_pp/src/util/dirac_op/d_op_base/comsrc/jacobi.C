#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOp class Ritz eigensolver methods.

  $Id: jacobi.C,v 1.4 2004/08/18 11:57:48 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/08/18 11:57:48 $
//  $Header: /space/cvs/cps/cps++/src/util/dirac_op/d_op_base/comsrc/jacobi.C,v 1.4 2004/08/18 11:57:48 zs Exp $
//  $Id: jacobi.C,v 1.4 2004/08/18 11:57:48 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_base/comsrc/jacobi.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/* JACOB */

/* This subroutine contains a "single node" Jacobi routine */
/* to be used with the Ritz functional eigenvialue/vector finder. */


/*  Psi		Eigenvectors			(Modify) */
/*  N_eig	Eigenvalue number 		(Read) */
/*  lambda	Diagonals / Eigenvalues		(Modify) */
/*  off_diag	Upper triang off-diag matrix elems	(Modify) */
/*  Toler	Tolerance for off-diag elems	(Read) */
/*  N_max	Maximal number of Jacobi iters	(Read) */
/*  N_Count	Number of Jacobi iters		(Write) */

CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <math.h>
CPS_START_NAMESPACE


int DiracOp::Jacobi(Vector **psi, int N_eig, Float *lambda, Complex *off_diag, 
		    Float Toler, int N_max)

{ /* Local Variables */
  Complex ctmp1;
  Complex ctmp2;
  Complex v12;
  Complex v21;

  Float v11;
  Float dd;
  Float ftmp;
  Float diff_l;
  Float theta;
  Float t;
  Float acc;
  Float c;
  Float s;
  Float al1;
  Float al2;
  int k;
  int i;
  int j;
  int ij;
  int m;
  int mi;
  int mj;
  int i_rot;
  int n_count = 0;
  Float tol_sq = Toler * Toler;
  char *fname = "Jacobi(V*,F,...)";

// Flash the LED and then turn it off
//------------------------------------------------------------------
  VRB.LedFlash(cname,fname,3);
  VRB.LedOff(cname,fname);
  VRB.Func(cname,fname);

// Set the node checkerboard size of the fermion field
//------------------------------------------------------------------
  int f_size = RitzLatSize();
    
// Set the node checkerboard size of the fermion field
//------------------------------------------------------------------
  Vector *psi_t1 = (Vector *) smalloc(f_size * sizeof(Float));
  if(psi_t1 == 0)
    ERR.Pointer(cname,fname, "psi_t1");
  VRB.Smalloc(cname,fname, "psi_t1", psi_t1, f_size * sizeof(Float));

  Vector *psi_t2 = (Vector *) smalloc(f_size * sizeof(Float));
  if(psi_t2 == 0)
    ERR.Pointer(cname,fname, "psi_t2");
  VRB.Smalloc(cname,fname, "psi_t2", psi_t2, f_size * sizeof(Float));

  for(k = 0; k <= N_max; k++)
  {
    i_rot = 0;
    ij = 0;
    for(j = 1; j < N_eig; j++)
      for(i = 0; i < j; i++)
      {
	dd = norm(off_diag[ij]);
	ftmp = fabs(lambda[i] * lambda[j] * tol_sq);

	if (dd > ftmp)
	{
	  /* Make a rotation to set off-diagonal part to zero */
	  i_rot++;
	  dd = sqrt(dd);
	  acc = 100.0 * dd;
	  diff_l = lambda[j] - lambda[i];
	  ftmp = fabs(diff_l);

	  if ((ftmp+acc) == ftmp)
	  {
	    t = dd / diff_l;
	  }
	  else
	  {
	    theta = 0.5 * diff_l / dd;
	    t = sqrt(1.0 + theta*theta);
	    ftmp = fabs(theta);
	    t = 1.0 / (ftmp+t);
	    if (theta < 0.0) t = -t;
	  }

	  if (diff_l >= 0)
	  {
	    c = sqrt(1.0 + t*t);
	    c = 1.0 / c;
	    s = - t * c;
	  }
	  else
	  {
	    s = sqrt(1.0 + t*t);
	    s = 1.0 / s;
	    c = t * s;
	  }

	  ftmp = c * c;
	  al1 = ftmp * lambda[i];
	  al2 = ftmp * lambda[j];
	  ftmp = s * s;
	  al1 += ftmp * lambda[j];
	  al2 += ftmp * lambda[i];
	  ftmp = 2.0 * dd * s * c;
	  al1 += ftmp;
	  al2 -= ftmp;
	  lambda[i] = al1;
	  lambda[j] = al2;
	  v11 = c;
	  ftmp = s / dd;
	  v12 = ftmp * off_diag[ij];
	  v21 = -conj(v12);
	  off_diag[ij] = 0.0;

	  /* Now rotate the eigenvectors */
	  psi_t1->CopyVec(psi[i], f_size);
	  psi_t1->VecTimesEquFloat(v11, f_size);
	  psi_t1->CTimesV1PlusV2(-v21, psi[j], psi_t1, f_size);

	  psi_t2->CopyVec(psi[j], f_size);
	  psi_t2->VecTimesEquFloat(v11, f_size);
	  psi_t2->CTimesV1PlusV2(-v12, psi[i], psi_t2, f_size);

	  psi[i]->CopyVec(psi_t1, f_size);
	  psi[j]->CopyVec(psi_t2, f_size);

	  /* Rotate the other matrix elements */
	  for(m = 0; m < N_eig; m++)
	  {
	    if( m != i && m != j )
	    {
	      if( m < i )
	      {
		mi = i * (i-1) / 2 + m;
		mj = j * (j-1) / 2 + m;
		ctmp1 = off_diag[mi] * v11 - off_diag[mj] * v21;
		ctmp2 = off_diag[mj] * v11 - off_diag[mi] * v12;
		off_diag[mi] = ctmp1;
		off_diag[mj] = ctmp2;
	      }
	      else if( m < j)
	      {
		mi = m * (m-1) / 2 + i;
		mj = j * (j-1) / 2 + m;
		ctmp1 = conj(off_diag[mi]) * v11 - off_diag[mj] * v21;
		ctmp2 = off_diag[mj] * v11 - conj(off_diag[mi]) * v12;
		off_diag[mi] = conj(ctmp1);
		off_diag[mj] = ctmp2;
	      }
	      else
	      {
		mi = m * (m-1) / 2 + i;
		mj = m * (m-1) / 2 + j;
		ctmp1 = conj(off_diag[mi]) * v11 - conj(off_diag[mj]) * v21;
		ctmp2 = conj(off_diag[mj]) * v11 - conj(off_diag[mi]) * v12;
		off_diag[mi] = conj(ctmp1);
		off_diag[mj] = conj(ctmp2);
	      }
	    }
	  }
	}

	ij++;
      }

    if( i_rot == 0 )
    {
      n_count = k;
      VRB.Result(cname, fname, "Jacobi converged after %d iters\n", k);

      /* Sort the eigenvalues */
      for(j = 1; j < N_eig; j++)
	for(i = 0; i < j; i++)
	{
	  if( fabs(lambda[j]) < fabs(lambda[i]) )
	  {
	    ftmp = lambda[i];
	    lambda[i] = lambda[j];
	    lambda[j] = ftmp;

	    psi_t1->CopyVec(psi[i], f_size);
	    psi_t2->CopyVec(psi[j], f_size);
	    psi[j]->CopyVec(psi_t1, f_size);
	    psi[i]->CopyVec(psi_t2, f_size);
	  }
	}

      VRB.Sfree(cname,fname, "psi_t2", psi_t2);
      sfree(psi_t2);
      VRB.Sfree(cname,fname, "psi_t1", psi_t1);
      sfree(psi_t1);

// Flash the LED and then turn it on
//------------------------------------------------------------------
      VRB.FuncEnd(cname,fname);
      VRB.LedFlash(cname,fname,2);
      VRB.LedOn(cname,fname);

      return n_count;
    }
  }

  VRB.Sfree(cname,fname, "psi_t2", psi_t2);
  sfree(psi_t2);
  VRB.Sfree(cname,fname, "psi_t1", psi_t1);
  sfree(psi_t1);

// Flash the LED and then turn it on
//------------------------------------------------------------------
  VRB.FuncEnd(cname,fname);
  VRB.LedFlash(cname,fname,2);
  VRB.LedOn(cname,fname);

  n_count = k;
  ERR.General(cname,fname, "too many Jacobi iterations");
  return n_count;
}

CPS_END_NAMESPACE
