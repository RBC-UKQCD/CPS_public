#include <config.h>
CPS_START_NAMESPACE
 /*! \file
  \brief  Definition of DiracOpWilsonTypes class eigensolver methods.

  $Id: eigen_wilson.C,v 1.3 2004-04-27 03:51:20 cwj Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: cwj $
//  $Date: 2004-04-27 03:51:20 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson_types/eigen_wilson.C,v 1.3 2004-04-27 03:51:20 cwj Exp $
//  $Id: eigen_wilson.C,v 1.3 2004-04-27 03:51:20 cwj Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: eigen_wilson.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson_types/eigen_wilson.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*  Compute the spectrol flow using the Ritz functional minimization */
/*  routine, if desired in the Kalkreuter-Simma algorithm */

/*  Compute low lying eigenvalues of the hermitian (but not pos. def.) */
/*  matrix RitzEigMat using the Ritz functional minimization routine, */
/*  if desired in the Kalkreuter-Simma algorithm */

/*  psi		Eigenvectors			(Modify) */
/*  lambda_H	Eigenvalues			(Write) */
/*  valid_eig	flag for valid eigenvalues	(Write) */
/*  eig_arg	Argument structure		(Read) */
/*  NCG_tot	Total number of CG iter		(Return value) */

/*  Parameters in eig_arg */
/*  N_eig	Eigenvalue number		(Read) */
/*  N_min	Minimal CG iterations		(Read) */
/*  N_max	Maximal CG iterations		(Read) */
/*  Kalk_Sim	Use Kalkreuter-Simma criterion	(Read) */
/*  n_renorm	Renormalize every n_renorm iter.	(Read) */
/*  N_KS_max	Max number of Kalkreuter-Simma iterations	(Read) */
/*  RsdR_a	(absolute) residue		(Read) */
/*  RsdR_r	(relative) residue		(Read) */
/*  Rsdlam	relative accuracy of lambda	(Read) */
/*  Cv_fact	"Convergence factor" required	(Read) */
/*  n_KS	total number of Kalkreuter-Simma iterations	(Write) */
/*  ProjApsiP	flag for projecting A.psi	(Read) */

/*  n_valid	number of valid eigenvalues	(Write) */

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


int DiracOpWilsonTypes::RitzEig(Vector **psi, Float lambda_H[], int valid_eig[], EigArg *eig_arg)
{
  int i;
  int n_jacob;
  
  // Initialize
  Float RsdR_a = eig_arg->RsdR_a;
  Float RsdR_r = eig_arg->RsdR_r;
  Float Rsdlam = eig_arg->Rsdlam;
  Float Cv_fact = eig_arg->Cv_fact;
  int Kalk_Sim = eig_arg->Kalk_Sim;
  int N_eig = eig_arg->N_eig;
  int N_min = eig_arg->N_min;
  int N_max = eig_arg->N_max;
  int N_KS_max = eig_arg->N_KS_max;
  int n_renorm = eig_arg->n_renorm;
  int MaxCG = eig_arg->MaxCG;
  int ProjApsiP = eig_arg->ProjApsiP;

  // Local vars which are for extension purposes
  Float lambda_t;
  Float dummy;
  Float del_lamb;
  int n;
  int j;
  int ij;
  int NCG_tot;
  int n_KS;
  int n_valid;

  char *fname = "RitzEig(V**,V**,I)";

// Flash the LED and then turn it off
//------------------------------------------------------------------
  VRB.LedFlash(cname,fname,3);
  VRB.LedOff(cname,fname);
  VRB.Func(cname,fname);


// Print out input parameters
//------------------------------------------------------------------
  VRB.Input(cname,fname,
	    "mass = %g\n",IFloat(eig_arg->mass));
  VRB.Input(cname,fname,
	    "N_eig = %d\n",N_eig);


//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------
  if (N_eig == 1)
    Kalk_Sim = 0;

  if (Kalk_Sim == 0)
    N_KS_max = 0;

  // Determine machine accuracy
  Float acc = 2.0e-6;
  acc *= acc;
  Float rsdl_sq = Rsdlam * Rsdlam;
  rsdl_sq = (rsdl_sq > acc) ? rsdl_sq : acc;
  Float rsdl_zero = 10.0 * rsdl_sq;


// Set the node checkerboard size of the fermion field
//------------------------------------------------------------------
  int f_size = RitzLatSize();
    
// Set the node checkerboard size of the fermion field
//------------------------------------------------------------------
  Float *lambda = (Float *) smalloc(N_eig * sizeof(Float));
  if(lambda == 0)
    ERR.Pointer(cname,fname, "lambda");
  VRB.Smalloc(cname,fname, "lambda", lambda, N_eig * sizeof(Float));

  Float *lambda_old = (Float *) smalloc(N_eig * sizeof(Float));
  if(lambda_old == 0)
    ERR.Pointer(cname,fname, "lambda_old");
  VRB.Smalloc(cname,fname, "lambda_old", lambda_old, N_eig * sizeof(Float));

  Complex *off_diag;
  if (N_eig > 1)
  {
    off_diag = (Complex *) smalloc(N_eig*(N_eig-1)/2 * sizeof(Complex));
    if(off_diag == 0)
      ERR.Pointer(cname,fname, "off_diag");
    VRB.Smalloc(cname,fname, "off_diag", off_diag, N_eig*(N_eig-1)/2 * sizeof(Complex));

    for(int i=0; i < N_eig*(N_eig-1)/2; ++i)
      off_diag[i] = 0.0;
  }

  NCG_tot = 0;

  /* Note: psi will be normalized in Ritz! */

  for(n=0; n < N_eig; ++n)
  {
    lambda_old[n] = 1.0;
    valid_eig[n] = 0;
  }

  del_lamb = 1.0;

  if (Kalk_Sim == 1)
  {
    /* while (del_lamb > Rsdlam) */
    for(n_KS = 1; n_KS <= N_KS_max && del_lamb > Rsdlam; n_KS++)
    {
      ij = 0;
      for(n = 1, i = 0; n <= N_eig; n++, i++)
      {
	NCG_tot += Ritz(psi, n, lambda_t, RsdR_a, RsdR_r, Rsdlam, rsdl_sq,
			n_renorm, Kalk_Sim, N_min, N_max, Cv_fact, MaxCG,
			ProjApsiP);

	lambda[i] = lambda_t;
	VRB.Debug(cname,fname,"yes KS: lambda[%d] = %g\n",i,(IFloat)lambda[i]);

	Vector *tmp = (Vector *) smalloc(f_size * sizeof(Float));
	if(tmp == 0)
	  ERR.Pointer(cname,fname, "tmp");
	VRB.Smalloc(cname,fname, "tmp", tmp, f_size * sizeof(Float));

	RitzMat(tmp, psi[i]);

	for(j = 0; j < i; j++)
	  off_diag[ij++] = psi[j]->CompDotProductGlbSum(tmp, f_size);

	VRB.Sfree(cname,fname, "tmp", tmp);
	sfree(tmp);
      }

      n_jacob = Jacobi(psi, N_eig, lambda, off_diag, RsdR_a, 50);

      for(i = 0; i < N_eig; i++)
	VRB.Result(cname,fname,"KS=%d: lambda[%d] = %g\n",n_KS,i,(IFloat)(lambda[i]));

      lambda_old[0] -= lambda[0];
      lambda_t = fabs(lambda_old[0]);
      if (lambda_t > rsdl_sq)
      {
	del_lamb = fabs(lambda_t / lambda[0]);
      }
      else
	del_lamb = 0.0;

      lambda_old[0] = lambda[0];
      for(i = 1; i < N_eig; i++)
      {
	lambda_old[i] -= lambda[i];
	lambda_t = fabs(lambda_old[i]);
	if (lambda_t > rsdl_sq)
	{
	  lambda_t = fabs(lambda_t / lambda[i]);
	  if (lambda_t > del_lamb)
	    del_lamb = lambda_t;
	}
	lambda_old[i] = lambda[i];
      }
    }

    if (n_KS >= N_KS_max)
      VRB.Result(cname,fname,"NONConversion_warning: n_KS = %d  N_KS_max = %d",
		 n_KS,N_KS_max);

    /* We have just found the eigenv of  H^2 */
    /* Diagonalize again for  H.v = lambda.v */
    Vector *tmp = (Vector *) smalloc(f_size * sizeof(Float));
    if(tmp == 0)
      ERR.Pointer(cname,fname, "tmp");
    VRB.Smalloc(cname,fname, "tmp", tmp, f_size * sizeof(Float));

    ij = 0;
    for(i = 0; i < N_eig; i++)
    {
      RitzEigMat(tmp, psi[i]);

      lambda_H[i] = psi[i]->ReDotProductGlbSum(tmp, f_size);
      lambda_old[i] = fabs(lambda_H[i]);

      for(j = 0; j < i; j++)
	off_diag[ij++] = psi[j]->CompDotProductGlbSum(tmp, f_size);
    }

    VRB.Sfree(cname,fname, "tmp", tmp);
    sfree(tmp);
  }
  else	/* Kalk_Sim == NO */
  {
    ij = 0;
    for(n = 1, i = 0; n <= N_eig; n++, i++)
    {
      NCG_tot += Ritz(psi, n, lambda_t, RsdR_a, RsdR_r, Rsdlam, rsdl_sq,
		      n_renorm, 0, 0, 0, Cv_fact, MaxCG, ProjApsiP);

      lambda[i] = lambda_t;
      VRB.Debug(cname,fname,"no KS: lambda[%d] = %g\n",i,(IFloat)lambda[i]);

      /* We have just found the eigenv of  H^2 */
      /* Diagonalize again for  H.v = lambda.v */
      Vector *tmp = (Vector *) smalloc(f_size * sizeof(Float));
      if(tmp == 0)
	ERR.Pointer(cname,fname, "tmp");
      VRB.Smalloc(cname,fname, "tmp", tmp, f_size * sizeof(Float));

      RitzEigMat(tmp, psi[i]);
	
      lambda_H[i] = psi[i]->ReDotProductGlbSum(tmp, f_size);
      lambda_old[i] = fabs(lambda_H[i]);

      for(j = 0; j < i; j++)
	off_diag[ij++] = psi[j]->CompDotProductGlbSum(tmp, f_size);

      VRB.Sfree(cname,fname, "tmp", tmp);
      sfree(tmp);
    }
  }


  // Check for a common failure mode. If all goes well, the eigenvalues
  // should be in increasing order out of Jacobi. Sometimes, the last one
  // is not since it didn't converge well enough, or there are nearby 
  // eigenvalues. The last one then will often have too small of an eigenvalue
  /// and be out of order. Toss it out if this happens.
  if (N_eig > 2)
  {
    lambda_t = lambda_old[N_eig-1] - lambda_old[N_eig-2] + Rsdlam;
    dummy = lambda_old[N_eig-2] - lambda_old[N_eig-3] + Rsdlam;
    if (lambda_t < 0.0 && dummy > 0.0)
    {
      n = N_eig - 1;
    }
    else
      n = N_eig;
  }
  else
    n = N_eig;

  dummy = sqrt(Rsdlam);

  if (N_eig > 1)
  {
    n_jacob = Jacobi(psi, n, lambda_H, off_diag, RsdR_a, 50);

    /* Label eigenvalues okay, or not */
    n_valid = 0;
    for(n = 0; n < N_eig; n++)
    {
      valid_eig[n] = 0;
      lambda_t = lambda_H[n]*lambda_H[n];
      if( lambda_t < rsdl_zero )
      {
	for(j = n_valid; j < N_eig; j++)
	  if( lambda[j] < rsdl_zero )
	  {
	    valid_eig[n] = 1;
	    n_valid++;
	    break;
	  }
      }
      else
      {
	for(j = n_valid; j < N_eig; j++)
	{
	  del_lamb = fabs(lambda[j] - lambda_t);
	  if( del_lamb < dummy*lambda[j] )
	  {
	    valid_eig[n] = 1;
	    n_valid++;
	    break;
	  }
	}
      }
    }
    
    VRB.Result(cname,fname,"Final Jacobi: n_jacob=%d n_valid=%d NCG=%d\n",
	       n_jacob,n_valid,NCG_tot); 
    for(n = 0; n < N_eig; n++)
      VRB.Result(cname,fname,"lambda_old[%d] = %g\n",n,(IFloat)lambda_old[n]);
    for(n = 0; n < N_eig; n++)
      VRB.Result(cname,fname,"lambda_H[%d] = %g\n",n,(IFloat)lambda_H[n]);
    for(n = 0; n < N_eig; n++)
      VRB.Result(cname,fname,"valid_eig[%d] = %d\n",n,valid_eig[n]);

    VRB.Sfree(cname,fname, "off_diag", off_diag);
    sfree(off_diag);
  }
  else	/* N_eig = 1 */
  {
    del_lamb = fabs(lambda_H[0] * lambda_H[0] - lambda[0]);
    if( del_lamb < dummy*lambda[0] )
    {
      valid_eig[0] = 1;
      n_valid = 1;
    }
    else
    {
      valid_eig[0] = 0;
      n_valid = 0;
    }

    VRB.Result(cname,fname,"Final Result: n_valid=%d  NCG=%g\n",n_valid,NCG_tot); 
    for(n = 0; n < N_eig; n++)
      VRB.Result(cname,fname,"lambda[%d] = %g\n",n,(IFloat)lambda[n]);
    for(n = 0; n < N_eig; n++)
      VRB.Result(cname,fname,"lambda_H[%d] = %g\n",n,(IFloat)lambda_H[n]);
    for(n = 0; n < N_eig; n++)
      VRB.Result(cname,fname,"valid_eig[%d] = %d\n",n,valid_eig[n]);
  }

  VRB.Sfree(cname,fname, "lambda_old", lambda_old);
  sfree(lambda_old);
  VRB.Sfree(cname,fname, "lambda", lambda);
  sfree(lambda);

// Flash the LED and then turn it on
//------------------------------------------------------------------
  VRB.FuncEnd(cname,fname);
  VRB.LedFlash(cname,fname,2);
  VRB.LedOn(cname,fname);
      
  return NCG_tot;
}

CPS_END_NAMESPACE
