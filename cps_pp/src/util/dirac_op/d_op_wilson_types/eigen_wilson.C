#include <config.h>

CPS_START_NAMESPACE
 /*! \file
  \brief  Definition of DiracOpWilsonTypes class eigensolver methods.

*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
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
#include <util/qcdio.h>
#include <math.h>
CPS_START_NAMESPACE
void MatHermElements( DiracOpWilsonTypes * dirac_op,
		      Vector ** psi, int n_vec, int f_size, 
 		      Float * diag, Complex * off_diag );

inline void PrintDot(char *fname, char *vname, Vector *tmp, int f_size){
//      printf("%s: %s=%e\n",fname,vname,tmp->ReDotProductGlbSum(tmp, f_size));
}

/*
  A normalised gram-shcmidt orthoganalisation the input vectors
  should be normalised
*/

inline Float Norm( Vector* psi, int f_size )
{
  Float norm( sqrt(psi->NormSqGlbSum(f_size)) );
//  printf("Norm=%e\n",norm);
  psi->VecTimesEquFloat(1.0/norm,f_size);
  return norm;
}

Float GramSchmNorm( Vector* psi, 
                    Vector* vec, 
                    int f_size  )
{
  const Complex xp(vec->CompDotProductGlbSum(psi, f_size));
  psi->CTimesV1PlusV2  (-xp, vec, psi, f_size);
  return Norm(psi,f_size); 
}



int DiracOpWilsonTypes::RitzEig(Vector **psi, Float lambda_H[], int valid_eig[], EigArg *eig_arg)
{
  char *fname = "RitzEig(V**,V**,I)";
  VRB.Func(cname,fname);


  // Initialise constants
  const Float RsdR_a  ( eig_arg->RsdR_a    );
  const Float RsdR_r  ( eig_arg->RsdR_r    );
  const Float Rsdlam  ( eig_arg->Rsdlam    );
  const Float Cv_fact ( eig_arg->Cv_fact   );
  const int N_eig     ( eig_arg->N_eig     );
  const int N_min     ( eig_arg->N_min     );
  const int N_max     ( eig_arg->N_max     );
  const int n_renorm  ( eig_arg->n_renorm  );
  const int MaxCG     ( eig_arg->MaxCG     );
  const int ProjApsiP ( eig_arg->ProjApsiP );

  int Kalk_Sim ( eig_arg->Kalk_Sim  );
  int N_KS_max ( eig_arg->N_KS_max  );

 
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
  // can't do a kal-sim if I only have one eigenvector
  if (N_eig == 1)    Kalk_Sim = 0;

  if (Kalk_Sim == 0)    N_KS_max = 0;

  // Determine machine accuracy
  Float acc = 2.0e-6;
  acc *= acc;
  Float rsdl_sq = Rsdlam * Rsdlam;
  rsdl_sq = (rsdl_sq > acc) ? rsdl_sq : acc;
  Float rsdl_zero = 10.0 * rsdl_sq;


// Set the node checkerboard size of the fermion field
//------------------------------------------------------------------
  const int f_size = RitzLatSize();

  // Local vars 
  Float lambda_t, del_lamb;
  int j, n, ij, NCG_tot, n_KS;
  int i;
  int n_jacob;

  
  Float * lambda = lambda_H;

   /*
    we are going to pass lambda^2 back to the calling function
    the calling function has lambda_H[N_eig...2N_eig-1] reserved so
    we let lambda_old point to lambda_H[N_eig], and lambda^2 will 
    be put there automatically
  */

  Float *lambda_old = & lambda_H[N_eig]; 
   
// Set the node checkerboard size of the fermion field
//------------------------------------------------------------------
  Complex *off_diag (0x0);
  if (N_eig > 1)
  {
    off_diag = (Complex *) smalloc(cname,fname, "off_diag", N_eig*(N_eig-1)/2 * sizeof(Complex));

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
              /* 
                 call the ritz. This call will return a negative number
                 in the case of an error. In this case the only error is
                 hitting the max_cg, in which case it will return -1;
              */
	      int ncnt =  Ritz(psi, n, lambda_t, RsdR_a, RsdR_r, Rsdlam, 0.0 /*rsdl_sq*/,
			       n_renorm, Kalk_Sim, N_min, N_max, Cv_fact, MaxCG,
			       ProjApsiP);

              // pass the error "up"
              if ( ncnt < 0 ) { if ( off_diag !=0x0 ) { sfree(off_diag); } ; return ncnt; }

	      NCG_tot += ncnt;
              VRB.Result(cname,fname ,"PCNTINFO KS=%d n_eig=%d ncnt=%d\n", n_KS, n, ncnt);

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

          VRB.Result(cname, fname, "NCG_tot=%d\n", NCG_tot);

	  for(i = 0; i < N_eig; i++)
	    VRB.Result(cname,fname,"KS=%d: lambda[%d] = %g\n",n_KS,i,(IFloat)(lambda[i]));

	  /*
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
	  */

	  //calculate the maximum of the relative changes of the eigenvalues
          
	  if ( n_KS == 1 ) 
            {
              for(i = 0; i < N_eig; i++) lambda_old[i] = lambda[i];
            }
          else 
            {
              del_lamb = 0.0;
              for(i = 0; i < N_eig; i++)
                {
                  lambda_t   = fabs((lambda_old[i]-lambda[i])/lambda[i]);
                  if   (lambda_t > del_lamb)
                    del_lamb = lambda_t;
                  lambda_old[i]=lambda[i];
                }
            }
	}
      

      // check to see if maxed out on KS steps.

      if (n_KS >= N_KS_max) {
          /* 
             want to be able to recover from this in alg_eig (by
             moving onto the next mass, so just return a negative
             value for the iteration number and cope with the problem
             higer up. Choose -2 so as to not clash with the error code
             from the ritz call
          */
	VRB.Warn(cname,fname,"NONConversion_warning: n_KS = %d  N_KS_max = %d\n",
		 n_KS,N_KS_max);
	
	if ( off_diag !=0x0 ) { sfree(off_diag); }
	return -2;
      }

    /* We have just found the eigenv of  H^2 */
    /* Diagonalize again for  H.v = lambda.v */

      /*
    Vector *tmp = (Vector *) smalloc(cname,fname, "tmp", tmp, f_size * sizeof(Float));

    ij = 0;
    for(i = 0; i < N_eig; i++)
    {
      RitzEigMat(tmp, psi[i]);

      lambda_H[i] = psi[i]->ReDotProductGlbSum(tmp, f_size);
      lambda_old[i] = fabs(lambda_H[i]);

      for(j = 0; j < i; j++)
	off_diag[ij++] = psi[j]->CompDotProductGlbSum(tmp, f_size);
    }

    sfree(cname,fname, "tmp", tmp);

      */


  }
  else	/* Kalk_Sim == NO */
  {
    ij = 0;
    for(n = 1, i = 0; n <= N_eig; n++, i++)
    {
      NCG_tot += Ritz(psi, n, lambda_t, RsdR_a, RsdR_r, Rsdlam, 0.0 /*rsdl_sq*/,
		      n_renorm, 0, 0, 0, Cv_fact, MaxCG, ProjApsiP);

//      printf("lambda_old[%d]=%e\n",i,lambda_t);
      lambda_old[i] = lambda_t;
      VRB.Debug(cname,fname,"no KS: lambda[%d] = %g\n",i,(IFloat)lambda_old[i]);

      /* We have just found the eigenv of  H^2 */
      /* Diagonalize again for  H.v = lambda.v */
      /*
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
      */
    }
  }



  /*
    We have just found the eigenv of  H^2 
    Diagonalize again for  H.v = lambda.v 
    (1) by the projection:   [1  \pm  H / sqrt(lambda^2)] .  v

    Currently select only one eigen vector from degenerated eigenvectors
    by the projection. 
  */

  /*
    calculates \psi^{\dagger} D_H \psi and puts the diagonal in lambda_H
    and the off_diagonal in off_diag
  */
  MatHermElements(this, psi, N_eig, f_size, lambda_H, off_diag); 


  /* 
     start off considering all the eigenvectors
     as possibly valid
  */
  for (i=0;i<N_eig;i++) { valid_eig[i] = 1; }



 /*
    then iteratively remove the eigenvectors that
    are invalid ( only count the elements of the
    jacobi matrix between valid eigenvectors ) 
  */
  
  bool looping(true);

  while ( looping )
    {
      looping = false;
      for (i=0;i<N_eig;i++)
        {
          if ( valid_eig[i] )
            {
              int xy(0),x,y;
              const Float rEig2( 1/lambda_old[i] );
              Float sum        (lambda_H[i]*lambda_H[i] * rEig2);
              // loop over all off-diagonal element and pick up all those
              // of the correct row
              for (x=0;x<N_eig;x++)
                {
                  for (y=0;y<x;y++)
                    {
                      if ( ( x==i && valid_eig[y] ) || ( y==i && valid_eig[x] ) )
                        {
                            Float value( abs(off_diag[xy]) );
                          sum+=value*value*rEig2;
                        }
                      xy++;
                    }
                }
              // sum should now have the sum of the squares of the elements
              // of row i of the matrix. If this deviates from 1 then we have 
              // some missing eigenvectors
              if ( sum < 0.98 )
                {
                  valid_eig[i] = 0;
                  // tagged a new invalid eigenvalue -> need
                  // a new iteration to take this into account
                  looping = true;
                }
            } // if valid_eig
        } // loop over eignevalue
    } // while ( looping )

  /*
    dump these out to the eigenvector output file so
    we can check for degeneracy
  */

  {
    int ij = 0;
    int i, j;
    FILE* fp(Fopen(eig_arg->fname,"a"));
    Fprintf(fp,"Eig2 before projection/jacobi\n");
//    printf("Eig2 before projection/jacobi\n");
    for(i = 0; i < N_eig; i++)
      {
        Fprintf(fp,"%i %e %i\n",i,lambda_old[i],valid_eig[i]);
      }
    for(i = 0; i < N_eig; i++) 
      {
        Fprintf(fp,"Jacobi matrix: %i %i %e\n",i,i,(float)lambda_H[i]);
        for(j = 0; j < i; j++)
          {
            Fprintf(fp,"Jacobi matrix: %i %i %e %e\n"
                    ,i
                    ,j
                    ,off_diag[ij].real()
                    ,off_diag[ij].imag()   );
            ij++;
          }
      }
    Fclose(fp);
  }
    
  int N_dp(0);
  for(i=0;i<N_eig;i++)
    {
      if ( valid_eig[i] == 1 ) N_dp++;
      else { break; }
    }
  // invalidate the rest of the eigenvectors
  VRB.Flow(cname,fname,"N_dp=%d\n",N_dp);
  for (i++;i<N_eig;i++) { valid_eig[i] = 0; }
  
  //allocate tmp vector
  Vector *tmp = (Vector *) smalloc(f_size * sizeof(Float));
  if(tmp == 0) ERR.Pointer("",fname, "tmp");
  Vector *tmp2 = (Vector *) smalloc(f_size * sizeof(Float));
  if(tmp2 == 0) ERR.Pointer("",fname, "tmp2");
  
  while ( N_dp < N_eig )
    {
      valid_eig[N_dp] = 1;
      
      const Float lam(sqrt(fabs(lambda_old[N_dp])));
      VRB.Flow(cname,fname,"lam=%e\n",lam);
            
      // tmp=  (HemiteMatrix) .  psi[N_dp] 
      RitzEigMat(tmp,psi[N_dp]);
      PrintDot(fname,"psi[N_sp]",psi[N_dp],f_size);
      PrintDot(fname,"tmp",tmp,f_size);

      const Float proj(tmp->ReDotProductGlbSum(psi[N_dp], f_size));
      VRB.Flow(cname,fname,"proj=%e\n",proj);
      
      Float sign(1);
      if ( proj <  0 ) { sign = -1; }

            
      tmp2->FTimesV1PlusV2( sign/lam, tmp, psi[N_dp], f_size );
      PrintDot(fname,"tmp2",tmp2,f_size);
      Rcomplex snorm(0,0),onorm(0,0);
      for (i=0;i<N_dp;i++)
        {
          snorm+=tmp2->CompDotProductGlbSum(psi[i],f_size);
        }
      tmp->FTimesV1PlusV2( -sign/lam, tmp, psi[N_dp], f_size );
      PrintDot(fname,"tmp",tmp,f_size);
      for (i=0;i<N_dp;i++)
        {
          onorm+=tmp->CompDotProductGlbSum(psi[i],f_size);
        }
      if ( norm(onorm) > norm(snorm) )
        {
          psi[N_dp]->CopyVec( tmp2, f_size );
        }
      else
        {
          psi[N_dp]->CopyVec( tmp, f_size );
        }
      Norm(psi[N_dp],f_size);
      PrintDot(fname,"psi[N_sp]",psi[N_dp],f_size);
      N_dp++;
    }
//  printf("GramSchm");
  for (i=0;i<N_eig;i++)
    {
      for(j=i+1;j<N_eig;j++)
        {
          GramSchmNorm(psi[j],psi[i],f_size);
          PrintDot(fname,"psi[i]",psi[i],f_size);
          PrintDot(fname,"psi[j]",psi[j],f_size);
        }
    }

  

  sfree(tmp);
  sfree(tmp2);
  
  /*
    Diagonalize again for  H.v = lambda.v 
    (2)  by the Jacobi diagonalization
    do Jacobi diagonalization only using 0..N_dp vector 
  */
  int min_dim(0);
  for(i=0;i<N_eig;i++){ if(valid_eig[i]) min_dim++; }
  
  /*
    re-calculates \psi^{\dagger} D_H \psi and puts the diagonal in lambda_H
    and the off_diagonal in off_diag
  */
  
  MatHermElements(this, psi, min_dim, f_size, lambda_H, off_diag); 
  
  {  
    int ij = 0;
    int i, j;
    FILE* fp(Fopen(eig_arg->fname,"a"));
    Fprintf(fp,"Jacobi matrix after projection\n");
//    printf("Jacobi matrix after projection\n");
    for(i = 0; i < N_eig; i++) 
      {
        Fprintf(fp,"Jacobi matrix: %i %i %e\n",i,i,(float)lambda_H[i]);
        for(j = 0; j < i; j++)
          {
            Fprintf(fp,"Jacobi matrix: %i %i %e %e\n"
                    ,i
                    ,j
                    ,off_diag[ij].real()
                    ,off_diag[ij].imag()   );
            ij++;
          }
      }
    Fclose(fp);
  }
  
  n_jacob = Jacobi(psi, min_dim, lambda_H, off_diag, 1.0e-5, 50);

  if(N_eig>1) sfree(off_diag);


  // end of calculation, dump some valuable information.
  
  VRB.Result(cname,fname,"Final Jacobi: n_jacob=%d  NCG=%d\n", 
	     n_jacob,NCG_tot); 
  
  for(n = 0; n < N_eig; n++)
    VRB.Result(cname,fname,"lambda_HH[%d] = %g\n",n,(float)lambda_old[n]);
  for(n = 0; n < N_eig; n++)
    VRB.Result(cname,fname,"lambda_H[%d] = %g\n",n,(float)lambda_H[n]);
  

  return NCG_tot;




#if 0
  // ################ OLD VERSION OF CODE ############
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

    VRB.Result(cname,fname,"Final Result: n_valid=%d  NCG=%d\n",n_valid,NCG_tot); 
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

#endif  // OLD CODE END
}



//-------------------------------------------------------------------
// Calculate the Matrix elements for the hermitian dirac operator.
// only the real part of the diagonal(the imaginary part should be 0.0
// and half of the complex off-diagonal elements are calculated.
//-------------------------------------------------------------------
void MatHermElements( DiracOpWilsonTypes * dirac_op,
		      Vector ** psi, int n_vec, int f_size, 
 		      Float * diag, Complex * off_diag ){

  char * fname  = "MatHermElements";
  Vector *tmp = (Vector *) smalloc("DiracOpWilsonTypes",fname,"tmp",f_size * sizeof(Float));
//  if(tmp == 0) ERR.Pointer("",fname, "tmp");
  
  int ij = 0;
  int i, j;

  for(i = 0; i < n_vec; i++) {
      dirac_op->RitzEigMat(tmp, psi[i]);
      diag[i] = psi[i]->ReDotProductGlbSum(tmp, f_size);
#if 0
      Float *tmp_p = (Float*)tmp;
      for(j = 0; j < 10; j++)
        printf("tmp[%d]=%e\n",j,*(tmp_p+j));
      tmp_p = (Float*)psi[i];
      for(j = 0; j < 10; j++)
        printf("psi[%d][%d]=%e\n",i,j,*(tmp_p+j));
      printf("%s: diag[%d]=%e\n",fname,i,diag[i]);
#endif
      PrintDot(fname,"psi[i]",psi[i],f_size);
      PrintDot(fname,"tmp",tmp,f_size);
      for(j = 0; j < i; j++){
	off_diag[ij++] = psi[j]->CompDotProductGlbSum(tmp, f_size);
#if 0
        tmp_p = (Float*) off_diag+ij;
        printf("%s: off_diag[%d]=%e %e\n",fname,ij-1,*tmp_p, *(tmp_p+1));
#endif
      }
  }

  sfree(tmp);
}



CPS_END_NAMESPACE
