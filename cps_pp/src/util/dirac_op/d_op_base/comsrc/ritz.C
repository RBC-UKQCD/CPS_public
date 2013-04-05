#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOp class Ritz eigensolver methods.

  $Id: ritz.C,v 1.12 2013-04-05 17:46:30 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:46:30 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_base/comsrc/ritz.C,v 1.12 2013-04-05 17:46:30 chulwoo Exp $
//  $Id: ritz.C,v 1.12 2013-04-05 17:46:30 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.12 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_base/comsrc/ritz.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// RITZ

// This subroutine minimizes the Ritz functional with a CG based
// algorithm to find the n-th lowest eigenvalue of a hermitian A

// lambda_n = min_z <z|Az>/<z|z>

// In the subspace orthogonal to the (n-1) lower eigenvectors of A

// This version includes the "early" termination criterion of
// Kalkreuter and Simma.

// Algorithm:

//  Apsi[0] :=  A . Psi[0] ;
//  mu[0]    :=  < Psi[0] | A . Psi[0] > ; 	Initial Ritz value
// Note: we assume  < Psi[0] |Psi[0] > = 1
//  p[1] = g[0]   :=  ( A - mu[0] ) Psi[0] ;	Initial direction
//  IF |g[0]| <= max(RsdR_a, RsdR_r |mu[0]|) THEN RETURN;	Converged?
//  FOR k FROM 1 TO MaxCG DO			CG iterations
//      s1 = (mu[k-1] + < p[k] | A p[k] >/p[k]|^2) / 2;
//      s2 = (mu[k-1] - < p[k] | A p[k] >/p[k]|^2) / 2;
//      s3 = |g[k-1]|^2 / |p[k]|;
//      Compute a and cos(theta), sin(theta)
//      (see DESY internal report, September 1994, by Bunk, Jansen,
//       Luescher and Simma)
//      lambda = mu[k] = mu[k-1] - 2 a sin^2(theta);
//      Psi[k] = cos(theta) Psi[k-1] + sin(theta) p[k]/|p[k]|;
//      Apsi[k] = cos(theta) Apsi[k-1] + sin(theta) A p[k]/|p[k]|;
//      g[k] = Apsi[k] - mu[k] Psi[k]
//      IF |g[k]| <= max(RsdR_a, RsdR_r |mu[k]|) 
//           && |del_lam| < Rsdlam |lambda| THEN RETURN;	Converged?
//      b[k+1] := cos(theta) |g[k]|^2 / |g[k-1]|^2;
//      p[k+1] := g[k] + b[k+1] (p[k] - Psi[k] < Psi[k] | p[k] >);	New direction

/* Arguments: */

/*  Psi_all	Eigenvectors			(Modify) */
/*  N_eig	Eigenvalue number 		(Read) */
/*  lambda	N_eig-th Eigenvalue		(Write) */
/*  RsdR_a	(absolute) residue		(Read) */
/*  RsdR_r	(relative) residue		(Read) */
/*  Rsdlam	relative accuracy of lambda	(Read) */
/*  Cutl_zero	Zero cut-off for lambda		(Read) */
/*  n_renorm	Renormalize every n_renorm iter.	(Read) */
/*  Ncb		Number of sublattices		(Read) */
/*  N_Count	Number of CG iteration		(Write) */
/*  Kalk_Sim	Use Kalkreuter-Simma criterion  (Read) */
/*  N_min	Minimal CG iterations		(Read) */
/*  N_max	Maximal CG iterations		(Read) */
/*  Cv_fact	"Convergence factor" required	(Read) */
/*  ProjApsiP	flag for projecting A.psi	(Read) */

 /* Local Variables: */

/*  psi			New eigenvector */
/*  p			Direction vector */
/*  Apsi		Temporary for  A.psi */
/*  Ap			Temporary for  A.p, and other */
/*  mu			Ritz functional value */
/*  g2			| g[k] |**2 */
/*  p2			| p[k] |**2 */
/*  k			CG iteration counter */
/*  b                   beta[k+1] */
/*  and some others... */

/* Global Variables: */

/*  MaxCG       Maximum number of CG iterations allowed */

/* Subroutines: */

/*  A		Apply matrix A to vector */


CPS_END_NAMESPACE
#include <stdio.h>
#include <stdlib.h>
CPS_START_NAMESPACE


CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <comms/glb.h>
#include <util/verbose.h>
#include <util/error.h>
#include <math.h>
CPS_START_NAMESPACE

//! Vector orthogonalisation
/*!
  Implementation of the Gram-Schmidt method to orthogonalise vectors.
  \param psi An array of vectors to be orthogonalised.
  \param Npsi The number of vectors to be orthogonalised.
  \param vec An array of vectors against which each vector in \a psi is to
  be orthogonalised.
  \pre The vectors in the array \a vec should have unit norm.
  \param Nvec The number of vectors in the array \a vec.
  \param f_size The length, in floating point numbers, of the vectors.
  \post Each vector in the array \a psi is orthogonal to all of the vectors
  in the array \a vec.
*/
void DiracOp::GramSchm(Vector **psi, int Npsi, Vector **vec, int Nvec, int f_size) 
{
  Complex xp;

  for(int s = 0; s < Npsi; ++s)
  {
    for(int i = 0; i < Nvec; ++i)
    {
      xp = vec[i]->CompDotProductGlbSum(psi[s], f_size);

      /* psi[s] = psi[s] - <vec[i],psi[s]> vec[i] */
      psi[s]->CTimesV1PlusV2(-xp, vec[i], psi[s], f_size);
    }
  }
}
// same as above, but orthogonalize only one vector
void DiracOp::GramSchm(Vector *psi, Vector **vec, int Nvec, int f_size) 
{
  Complex xp;

  for(int i = 0; i < Nvec; ++i)
    {
      xp = vec[i]->CompDotProductGlbSum(psi, f_size);
      /* psi = psi - <vec[i],psi> vec[i] */
      psi->CTimesV1PlusV2(-xp, vec[i], psi, f_size);
    }
}

/*!
  The eigensolver implemented here finds the \e n th lowest eigenvalue
  \e lambda_n of a of a hermitian matrix \a A by using a Conjugate Gradient
  based method to minimize the Ritz functional in the subspace orthogonal
  to the \e (n-1) lower eigenvectors of \e A.
*/

int DiracOp::Ritz(Vector **psi_all, int N_eig, Float &lambda, 
		  Float RsdR_a, Float RsdR_r, Float Rsdlam, 
		  Float Cutl_zero, int n_renorm, int Kalk_Sim,
		  int N_min, int N_max, Float Cv_fact,
		  int MaxCG, int ProjApsiP)

{ /* Local Variables */

  int n_count;
  Complex xp;
  Float mu;  // Double
  Float p2;  // Double
  Float g2;  // Double
  Float a;  // Double
  Float d;  // Double
  Float s1;  // Double
  Float s2;  // Double
  Float s3;  // Double
  Float g2_0;
  Float del_lam;
  Float rsd;
  Float rsda_sq;
  Float rsdr_sq;
  Float ct;
  Float st;
  Float b;
  Float acc;
  int k;
  int nn;
  char *fname = "Ritz(V*,F,...)";

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
  Vector *psi = (Vector *) smalloc(cname,fname, "psi", f_size * sizeof(Float));

  Vector *p = (Vector *) smalloc(cname,fname, "p", f_size * sizeof(Float));

  Vector *Apsi = (Vector *) smalloc(cname,fname, "Apsi", f_size * sizeof(Float));
  Vector *Ap = (Vector *) smalloc(cname,fname, "Ap", f_size * sizeof(Float));

  VRB.Input(cname,fname,"RsdR_a=%e RsdR_r=%e Rsdlam=%e Cv_fact=%e\n", 
      RsdR_a,RsdR_r,Rsdlam,Cv_fact);

  //TIZB acc = 2.0e-6;
  acc = 2.0e-24;
  acc *= acc;
  rsda_sq = RsdR_a * RsdR_a;
  rsdr_sq = RsdR_r * RsdR_r;
  rsda_sq = (rsda_sq > acc) ? rsda_sq : acc;


  /*  Make Psi_all(N_eig-1) orthogonal to the previous eigenvectors  */
  nn = N_eig - 1;

  /*  First copy it into a the local psi  */
  psi->CopyVec(psi_all[nn], f_size);

  GramSchm(&psi, 1, psi_all, nn, f_size);

  /* Normalize */
  d = sqrt(psi->NormSqGlbSum(f_size));
  ct = 1.0 / d;
  psi->VecTimesEquFloat(ct, f_size);

  /* Now we can start */
  /*  Apsi[0]   :=  A . Psi[0]  */
  RitzMat(Apsi, psi);

  xp = psi->CompDotProductGlbSum(Apsi, f_size);
  mu = psi->ReDotProductGlbSum(Apsi, f_size);

  /*  Project to orthogonal subspace, if wanted  */
  /** Should not be necessary, following Kalkreuter-Simma **/
  if (ProjApsiP == 1)
    GramSchm(&Apsi, 1, psi_all, nn, f_size);

  /*  mu  := < Psi[0] | A Psi[0] >  */
  mu = psi->ReDotProductGlbSum(Apsi, f_size);

  /*  p[0] = g[0]   :=  ( A - mu[0] ) Psi[0]  =  Apsi - mu psi */
  lambda = mu;
  p->FTimesV1PlusV2(-lambda, psi, Apsi, f_size);
  

  /*  g2 = p2 = |g[0]|^2 = |p[0]|^2 */
  g2 = p->NormSqGlbSum(f_size);
  p2 = g2;
  g2_0 = g2;
  g2_0 *= Cv_fact;
  g2_0 = (g2_0 > acc) ? g2_0 : acc;

  if (ProjApsiP == 1)
    VRB.Result(cname,fname,"nn = %d, lambda=%g, g2=%g, g2_0=%g\n",
	    nn, (IFloat)mu, (IFloat)g2, (IFloat)g2_0);

  /*  IF |g[0]| <= min(RsdR_a, RsdR_r |mu[0]|) THEN RETURN; */
  rsd = rsdr_sq * lambda * lambda;
  rsd = (rsda_sq > rsd) ? rsda_sq : rsd;
  if ( Kalk_Sim == 0 && N_min <= 0 && g2 <= rsd )
  {
    /*  Copy psi back into psi_all(N_eig-1)  */
    psi_all[nn]->CopyVec(psi, f_size);

    sfree(cname,fname, "Ap", Ap);
    sfree(cname,fname, "Apsi", Apsi);
    sfree(cname,fname, "p", p);
    sfree(cname,fname, "psi", psi);

    n_count = 0;
    return n_count;
  }

  del_lam = mu;

  /*  FOR k FROM 1 TO MaxCG DO */
  for(k = 1; k <= MaxCG; ++k)
  {
    if (k % 100 == 0){
      VRB.Result(cname,fname,"nn = %d, iter=%d, lambda=%e, del_lam=%e\n", 
          nn, k, (IFloat)mu, del_lam);
      VRB.Result(cname,fname, "g2/g2_0=%e/%e=Cv_fact=%e > %e\n",
(IFloat)g2, (IFloat)g2_0, (IFloat)g2/g2_0, (IFloat)Cv_fact);
    }
//   if (k % 100 == 0)
//      VRB.Result(cname,fname,"nn = %d, lambda=%g, g2=%g, g2_0=%g\n",
//		 nn, (IFloat)mu, (IFloat)g2, (IFloat)g2_0);

    /*  Ap = A * p  */
    RitzMat(Ap, p);

    /*  Project to orthogonal subspace, if wanted  */
    /** Should not be necessary, following Kalkreuter-Simma **/
    if (ProjApsiP == 1)
      GramSchm(&Ap, 1, psi_all, nn, f_size);

    /*  d = < p | A p >  */
    d = p->ReDotProductGlbSum(Ap, f_size);

    d = d / p2;
    s1 = 0.5 * (mu+d);
    s2 = 0.5 * (mu-d);
    p2 = sqrt(p2);
    p2 = 1.0 / p2;
    s3 = g2 * p2;
    a = (s2 > 0.0) ? s2 : -s2;
    if( a >= s3 )
    {
      d = s3 / s2;
      d = sqrt(1.0 + d*d);
      a = a * d;
    }
    else
    {
      d = s2 / s3;
      d = sqrt(1.0 + d*d);
      a = s3 * d;
    }

    s2 /= a;			/* Now s2 is cos(delta) */
    s3 /= a;			/* Now s3 is sin(delta) */
    if( s2 > 0 )
    {
      s2 = 0.5 * (1.0+s2);
      d = sqrt(s2);
      d = - d;			/* Now d is sin(theta) */
      s2 = -0.5 * s3 / d;	/* Now s2 is cos(theta) */
    }
    else
    {
      s2 = 0.5 * (1.0-s2);
      s2 = sqrt(s2);		/* Now s2 is cos(theta) */
      d = -0.5 * s3 / s2;	/* Now d is sin(theta) */
    }

    /* mu[k] = mu[k-1] - 2 a d^2 */
    s1 = 2.0 * a * d * d;
    mu -= s1;
    lambda = mu;
    del_lam = s1;

    st = d*p2;		        /* Now st is sin(theta)/|p| */
    ct = s2;		        /* Now ct is cos(theta) */

    /*  Psi[k] = ct Psi[k-1] + st p[k-1] */
    /*  Apsi[k] = ct Apsi[k-1] + st Ap */
    psi->VecTimesEquFloat(ct, f_size);
    psi->FTimesV1PlusV2(st, p, psi, f_size);
    Apsi->VecTimesEquFloat(ct, f_size);
    Apsi->FTimesV1PlusV2(st, Ap, Apsi, f_size);

    /*  Ap = g[k] = Apsi[k] - mu[k] Psi[k] */
    Ap->FTimesV1PlusV2(-lambda, psi, Apsi, f_size);

    /*  g2  =  |g[k]|**2 = |Ap|**2 */
    s1 = g2;			/* Now s1 is |g[k-1]|^2 */
    g2 = Ap->NormSqGlbSum(f_size);

    /*+  */
    /*  IF |g[k]| <= min(RsdR_a, RsdR_r |mu[k]|)
	   && |del_lam| <= Rsdlam*|lambda|  THEN RETURN; */
    rsd = rsdr_sq * lambda * lambda;
    rsd = (rsda_sq > rsd) ? rsda_sq : rsd;
    st = Rsdlam * fabs(mu);	/* old value of st is no longer needed */

/*    if ( TO_REAL(g2) <= rsd && del_lam <= st ) */
    if ( (Kalk_Sim == 0 && k >= N_min && (fabs(mu) < Cutl_zero ||
	       (g2 <=  rsd && del_lam <= st ) ) ) ||
	 (Kalk_Sim == 1 && ( del_lam <= st || fabs(mu) < Cutl_zero ||
	       (k >= N_min && (g2 < g2_0 || k >= N_max ) ) ) ) )
    {
      VRB.Result(cname, fname, "Converged at iter %d, lambda = %g, del_lam = %g\n",
		 k, (IFloat)lambda, (IFloat)del_lam);
      VRB.Result(cname, fname, "  rsd = %g, g2 = %g, g2_0 = %g\n",
		 (IFloat)rsd, (IFloat)g2, (IFloat)g2_0);
      n_count = k;

      /* Renormalize and recompute lambda */
      /*  Project to orthogonal subspace  */
      GramSchm(&psi, 1, psi_all, nn, f_size);

      d = sqrt(psi->NormSqGlbSum(f_size));
      ct = 1.0 / d;
      psi->VecTimesEquFloat(ct, f_size);
      d -= 1.0;
      VRB.Result(cname, fname, "Deviation at convergence: %g\n", (IFloat)d);

      /*  Apsi  :=  A . Psi  */
      RitzMat(Apsi, psi);

      /*  mu  := < Psi | A Psi >  */
      s1 = mu;
      mu = psi->ReDotProductGlbSum(Apsi, f_size);
      lambda = mu;
      VRB.Result(cname, fname, "Mu-s at convergence: old %g vs. nn %g\n", 
		 (IFloat)s1, (IFloat)mu);

      /*  Copy psi back into psi_all(N_eig-1)  */
      psi_all[nn]->CopyVec(psi, f_size);

      VRB.Sfree(cname,fname, "Ap", Ap);
      sfree(Ap);
      VRB.Sfree(cname,fname, "Apsi", Apsi);
      sfree(Apsi);
      VRB.Sfree(cname,fname, "p", p);
      sfree(p);
      VRB.Sfree(cname,fname, "psi", psi);
      sfree(psi);

// Flash the LED and then turn it on
//------------------------------------------------------------------
      VRB.FuncEnd(cname,fname);
      VRB.LedFlash(cname,fname,2);
      VRB.LedOn(cname,fname);
      
      return k;
    }

    /* b = beta[k] = cos(theta) |g[k]|^2 / |g[k-1]|^2 */
    b = ct * g2 / s1;
    ct *= 0.05 * b;
    d = sqrt(g2);
    ct /= (p2*d);
    if( ct > 1.0 )
    {
      /* Restart: p[k] = g[k] = Ap */
      VRB.Result(cname, fname, "Restart at iter %d since beta = %g\n", 
		 k, (IFloat)b);
      p->CopyVec(Ap, f_size);
    }
    else
    {
      /* xp = < Psi[k] | p[k-1] > */
      xp = psi->CompDotProductGlbSum(p, f_size);

      /* p[k] = g[k] + b (p[k-1] - xp psi[k]) */
      p->CTimesV1PlusV2(-xp, psi, p, f_size);
      p->FTimesV1PlusV2(b, p, Ap, f_size);
    }

    if( k%n_renorm == 0 )
    {
      /* Renormalize, and re-orthogonalize */
      /*  Project to orthogonal subspace  */
      GramSchm(&psi, 1, psi_all, nn, f_size);

      /* Normalize */
      d = sqrt(psi->NormSqGlbSum(f_size));
      ct = 1.0 / d;
      psi->VecTimesEquFloat(ct, f_size);
      d -= 1.0;

      if (ProjApsiP == 1)
	VRB.Result(cname,fname,"Deviation at iter %d: %g\n", k, (IFloat)d);

      /*  Apsi  :=  A . Psi  */
      RitzMat(Apsi, psi);

      /*  Project to orthogonal subspace, if wanted  */
      /** Should not be necessary, following Kalkreuter-Simma **/
      if (ProjApsiP == 1)
	GramSchm(&Apsi, 1, psi_all, nn, f_size);

      /*  mu  := < Psi | A Psi >  */
      s1 = mu;
      mu = psi->ReDotProductGlbSum(Apsi, f_size);
      if (ProjApsiP == 1)
	VRB.Result(cname,fname,"Mu-s at iter %d: old %g vs. new %g\n", 
		   k, (IFloat)s1, (IFloat)mu);

      /*  g[k] = Ap = ( A - mu ) Psi  */
      lambda = mu;
      Ap->FTimesV1PlusV2(-lambda, psi, Apsi, f_size);

      /*  g2  =  |g[k]|**2 = |Ap|**2 */
      g2 = Ap->NormSqGlbSum(f_size);

      VRB.Result(cname,fname,"g2 at iter %d: %g\n", k, (IFloat)g2);

      /*  Project p[k] to orthogonal subspace  */
      GramSchm(&p, 1, psi_all, nn, f_size);

      /*  Make p[k] orthogonal to Psi[k]  */
      GramSchm(&p, 1, &psi, 1, f_size);

      /*  Make < g[k] | p[k] > = |g[k]|**2: p[k] = p_old[k] + xp g[k],
	  xp = (g2 - < g | p_old >) / g2; g[k] = Ap */
      xp = g2 - Ap->CompDotProductGlbSum(p, f_size);
      ct = 1.0 / g2;
      xp *= ct;
      p->CTimesV1PlusV2(xp, Ap, p, f_size);
    }
    else if (ProjApsiP == 1 && nn > 0)
    {
      /*  Project psi and p to orthogonal subspace  */
      GramSchm(&psi, 1, psi_all, nn, f_size);
      GramSchm(&p, 1, psi_all, nn, f_size);
    }

    /*  p2  =  |p[k]|**2 */
    p2 = p->NormSqGlbSum(f_size);

    //TIZB
#if 1
    if (ProjApsiP == 1)
      VRB.Result(cname,fname,"At iter %d, lambda = %g, del_lam = %g, g2 = %g\n",
	     k, (IFloat)lambda, (IFloat)del_lam, (IFloat)g2);
#endif
  }

  /*  Copy psi back into psi_all(N_eig-1)  */
  psi_all[nn]->CopyVec(psi, f_size);

  n_count = MaxCG;
  VRB.Sfree(cname,fname, "Ap", Ap);
  sfree(Ap);
  VRB.Sfree(cname,fname, "Apsi", Apsi);
  sfree(Apsi);
  VRB.Sfree(cname,fname, "p", p);
  sfree(p);
  VRB.Sfree(cname,fname, "psi", psi);
  sfree(psi);

  ERR.General(cname,fname, "too many CG/Ritz iterations: %d\n",n_count);
  return n_count;
}

CPS_END_NAMESPACE
