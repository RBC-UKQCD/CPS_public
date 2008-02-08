#include <config.h>
#include <stdio.h>
CPS_START_NAMESPACE
 /*! \file
   \brief  Definition of DiracOp class BiCGstab(n) solver method.
 */

CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <comms/glb.h>
#include <util/verbose.h>
#include <util/error.h>
#include <math.h>

#define PROFILE
#ifdef PROFILE
#include <util/time_cps.h>
#endif

CPS_START_NAMESPACE

/*!
  Solves \f$ M f_{out} = f_{in} \f$ for \f$ f_{out}\f$
  where \a M is the (possibly odd-even preconditioned) fermionic matrix.

  \param psi The solution vector
  \param chi The source vector
  \param chi_norm The square norm of the source vector
  \param n Which type of BiCGstab(n) are we doing
  \return The number of iterations performed.
 */
int DiracOpWilsonTypes::BiCGstab(Vector *psi, Vector *chi, 
				 Float chi_norm, int n, Float *true_res)
{
  char fname[30];
  sprintf(fname,"BiCGstab(V*,V*,F,%d,F*)",n);
  VRB.Func(cname,fname);
  
// Flash the LED and then turn it off
//------------------------------------------------------------------
  VRB.LedFlash(cname,fname,3);
  VRB.LedOff(cname,fname);
  VRB.Func(cname,fname);


// Print out input parameters
//------------------------------------------------------------------
  VRB.Input(cname,fname,
	    "max_num_iter = %d\n",dirac_arg->max_num_iter);
  VRB.Input(cname,fname,
	    "mass = %e\n",IFloat(dirac_arg->mass));
  VRB.Input(cname,fname,
	    "src_norm_sq = %e\n",IFloat(chi_norm));
  VRB.Input(cname,fname,
	    "n = %d\n", n);

//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------

  int k, f_size;

  if (n <= 0)
    ERR.General(cname, fname, "Invalid n parameter: %d\n",n);
  const int n_max=10;
  if (n > n_max)
    ERR.General(cname, fname, "n(%d)>n_max(%d)\n",n,n_max);

  if(lat.Fclass() == F_CLASS_CLOVER)
    f_size = lat.FsiteSize()*GJP.VolNodeSites() / 2;
  else
    f_size = lat.FsiteSize()*GJP.VolNodeSites() / (lat.FchkbEvl()+1);

  Vector **r = (Vector**)smalloc((n+1)*sizeof(Vector*),"r",fname,cname);
  Vector **u = (Vector**)smalloc((n+1)*sizeof(Vector*),"u",fname,cname);
  for (int i=0; i<=n; i++) {
    r[i] = (Vector*)smalloc(f_size*sizeof(Float),"r[i]",cname,fname);
    u[i] = (Vector*)smalloc(f_size*sizeof(Float),"u[i]",cname,fname);
  }
  
  Vector *r0 = (Vector*)smalloc(f_size*sizeof(Float),"r0",fname,cname);
  
  Complex rho0(1.0,0.0), rho1;
  Complex alpha(0.0,0.0), omega(1.0,0.0), beta;

  Float sigma[n_max+1];
  Complex gamma[n_max+1], gamma_prime[n_max+1], gamma_prime_prime[n_max+1];
  Complex tau[n_max+1][n_max+1];
    
  Float rsd_sq;
  Float rsdcg_sq;
  int converged = 0;

  // If source isn't supplied, calculate it
  if (chi_norm == 0.0) {
    chi_norm = chi->NormSqNode(f_size);
    DiracOpGlbSum(&chi_norm);
  }

  // If source norm = 0, solution must be 0
  if (chi_norm == 0.0) {
    psi->VecZero(f_size);
    return 0;
  }

  rsdcg_sq = dirac_arg->stop_rsd * dirac_arg->stop_rsd;
  rsd_sq = chi_norm*rsdcg_sq;
  VRB.Flow(cname,fname, "stp_cnd =%e\n", IFloat(rsdcg_sq));

  // Calculate initial residual vector
  MatPc(r[0],psi);
  r[0] -> FTimesV1MinusV2(1.0,chi,r[0],f_size);
  sigma[0] = r[0] -> NormSqNode(f_size);
  DiracOpGlbSum(&sigma[0]);

  VRB.Flow(cname,fname,"|res[0]|^2 = %e\n", sigma[0]);

  // Set r0
  r0 -> CopyVec(r[0], f_size);

  u[0] -> VecZero(f_size);

#ifdef PROFILE
  struct timeval start;
  struct timeval end;
  CGflops = 0;
  gettimeofday(&start,NULL);
#endif

  // for k=0 until MaxCG do
  // if |psi[k+1] - psi[k]| <= RsdCG |psi[k+1]| then return
  for (k=1; k<=dirac_arg->max_num_iter && !converged; k+=n) {

    rho0 *= -omega;
    
    // BiCG part
    for (int j=0; j<n; j++) {
      rho1 = r0 -> CompDotProductGlbSum(r[j], f_size);
      beta = alpha * rho1 / rho0;
      rho0 = rho1;
      for (int i=0; i<=j; i++) {
	u[i] -> CTimesV1PlusV2(-beta,u[i],r[i],f_size);
      }
      MatPc(u[j+1],u[j]);
      alpha = rho0 / r0->CompDotProductGlbSum(u[j+1], f_size);
      for (int i=0; i<=j; i++) {
	r[i] -> CTimesV1PlusV2(-alpha, u[i+1], r[i], f_size);
      }
      MatPc(r[j+1], r[j]);
      psi -> CTimesV1PlusV2(alpha, u[0], psi, f_size);
    }

    // MR part
    for (int j=1; j<=n; j++) {
      for (int i=1; i<j; i++) {
	tau[i][j] = r[i] -> CompDotProductGlbSum(r[j],f_size) / sigma[i];//
	r[j] -> CTimesV1PlusV2(-tau[i][j], r[i], r[j], f_size);
      }
      sigma[j] = r[j] -> NormSqNode(f_size);
      DiracOpGlbSum(&sigma[j]);
      gamma_prime[j] = (r[j] -> CompDotProductGlbSum(r[0], f_size))/sigma[j];
    }
    
    gamma[n] = gamma_prime[n];
    omega = gamma[n];
    for (int j=n-1; j>0; j--) {
      gamma[j] = gamma_prime[j];
      for (int i=j+1; i<=n; i++) gamma[j] -= tau[j][i]*gamma[i];
    }

    for (int j=1; j<n; j++) {
      gamma_prime_prime[j] = gamma[j+1];
      for (int i=j+1; i<n; i++) gamma_prime_prime[j] += tau[j][i]*gamma[i+1];
    }

    u[0] -> CTimesV1PlusV2(-gamma[n], u[n], u[0], f_size);
    psi -> CTimesV1PlusV2(gamma[1], r[0], psi, f_size);
    r[0] -> CTimesV1PlusV2(-gamma_prime[n], r[n], r[0], f_size);    

    for (int j=1; j<n; j++) {
      u[0] -> CTimesV1PlusV2(-gamma[j], u[j], u[0], f_size);
      psi -> CTimesV1PlusV2(gamma_prime_prime[j], r[j], psi, f_size);
      r[0] -> CTimesV1PlusV2(-gamma_prime[j], r[j], r[0], f_size);
    }

    sigma[0] = r[0] -> NormSqNode(f_size);
    DiracOpGlbSum(&sigma[0]);

    if (sigma[0] < rsd_sq) converged = 1;

    VRB.Flow(cname,fname,
	     "|res[%d]|^2 = %e\n", k-1+n, sigma[0]);

  }
  
#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(cname,fname,CGflops,&start,&end); 
#endif

  // It has not reached stp_cnd: Issue a warning
  if (k >= dirac_arg->max_num_iter) 
    VRB.Warn(cname,fname,
	     "CG reached max iterations = %d. |res|^2 = %e\n",
	     k, sigma[0]);

  if (converged) {
    VRB.Result(cname,fname,
	       "True |res| / |src| = %e, iter = %d\n", sqrt(sigma[0]), k-1);
    if (true_res != 0) *true_res = sqrt(sigma[0]);
  }

  sfree(r0, "r0", fname, cname);
  for (int i=0; i<n+1; i++) {
    sfree(u[i], "u[i]", fname, cname);
    sfree(r[i], "r[i]", fname, cname);
  }
  sfree(u, "u", fname, cname);
  sfree(r, "r", fname, cname);

  return k;
}

CPS_END_NAMESPACE
