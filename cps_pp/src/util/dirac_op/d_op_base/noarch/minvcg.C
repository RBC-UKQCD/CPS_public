#include <stdio.h>
#include <config.h>
CPS_START_NAMESPACE
 /*! \file
   \brief  Definition of DiracOp class multishift CG solver method.

   $Id: minvcg.C,v 1.8 2004-12-21 19:02:39 chulwoo Exp $
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
CPS_START_NAMESPACE

/*!
  Solves \f$ (M^\dagger M + shift) f_{out} = f_{in} \f$ for \f$ f_{out}\f$
  for a given number of shifts,
  where \a M is the (possibly odd-even preconditioned) fermionic matrix.

  \param psi The solution vectors
  \param chi The source vector
  \param chi_norm The square norm of the source vector
  \param shift The shifts
  \param The mass parameters  for each shift
  \param Nmass The number of shifts
  \param isz The smallest shift
  \param RsdCG The target residuals for each shift
  \param type The type of multimass inverter.
   -  If type is ::MULTI, then regular multishift inversion is performed with
      each solution stored separately.
   -  If type is ::SINGLE, then each solution is multiplied by an
      amount in parameter \a alpha and summed to a single solution vector.
  \param alpha The contribution of each shifted solution to the total
  solution vector if \a type is SINGLE.

  \return The number of iterations performed.
 */
int DiracOp::MInvCG(Vector *psi, Vector *chi, Float chi_norm, Float *mass, 
		    int Nmass, int isz, Float *RsdCG,
		    MultiShiftSolveType type, Float *alpha)
{
  char *fname = "MInvCG(V*,V**,...)";
  VRB.Func(cname,fname);
  
// Flash the LED and then turn it off
//------------------------------------------------------------------
  VRB.LedFlash(cname,fname,3);
  VRB.LedOff(cname,fname);
  VRB.Func(cname,fname);


// Print out input parameters
//------------------------------------------------------------------
  VRB.Input(cname,fname,
	    "number of shifts = %d\n",Nmass);
  VRB.Input(cname,fname,
	    "smallest shift stop_rsd = %e\n",IFloat(RsdCG[isz]));
  VRB.Input(cname,fname,
	    "max_num_iter = %d\n",dirac_arg->max_num_iter);
  VRB.Input(cname,fname,
	    "mass = %e\n",IFloat(dirac_arg->mass));
  VRB.Input(cname,fname,
	    "src_norm_sq = %e\n",IFloat(chi_norm));

//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------

  int iz, k, s;
  int n_vec;

  if(lat.Fclass() == F_CLASS_CLOVER)
      n_vec = GJP.VolNodeSites() / 2;
  else
      n_vec = GJP.VolNodeSites() / (lat.FchkbEvl()+1);
  const int f_size = n_vec * lat.FsiteSize();
  
  
  
  Vector *r = (Vector *)smalloc(f_size * sizeof(Float),
				cname,fname,"r");
  
  Vector *Ap = (Vector *)smalloc(f_size * sizeof(Float),
				 cname,fname, "Ap");
  
  Vector **p = (Vector **)smalloc(Nmass * sizeof(Vector*),
				  cname,fname, "p");
  for (s=0; s<Nmass; s++) 
      *(p+s) = (Vector*)smalloc(f_size * sizeof(Float),
				cname,fname, "p[i]");
  
  
  int convP;
  int *convsP = (int*)smalloc(Nmass*sizeof(int));
  
  Float a=0, as, b, bp, b_tmp;
  Float *bs = (Float*)smalloc(Nmass * sizeof(Float));
  Float **z = (Float**)smalloc(2 * sizeof(Float*));
  for (s=0; s<2; s++) *(z+s) = (Float*)smalloc(Nmass * sizeof(Float));
  Float css, ztmp;
  
  Float c, cs, d, cp;
  Float *dot=0;

  Float *rsd_sq = (Float*)smalloc(Nmass*sizeof(Float)); 
  Float *rsdcg_sq = (Float*)smalloc(Nmass*sizeof(Float)); 

  // If source norm = 0, solution must be 0
  if (chi_norm == 0.0) {
    if (type == SINGLE) psi->VecTimesEquFloat(0.0,f_size);
    else for (k=0; k<Nmass; k++) (psi + n_vec*k)->VecTimesEquFloat(0.0,f_size);
    return 0;
  }
  
  r-> CopyVec(chi,f_size);
  for (s=0; s<Nmass; s++) p[s] -> CopyVec(chi,f_size);
  cp = chi_norm;

  for (s=0; s<Nmass; s++) {
    rsdcg_sq[s] = RsdCG[s]*RsdCG[s];
    rsd_sq[s] = cp*rsdcg_sq[s];
  }

  /*  d = <p, A.p>  */
  if (mass[isz] > 0) {
    MatPcDagMatPc(Ap,p[isz],dot);
    Ap -> FTimesV1PlusV2(mass[isz],p[isz],Ap, f_size);
    d = p[isz] -> ReDotProductGlbSum(Ap, f_size);
  } else {
    MatPcDagMatPc(Ap,p[isz],&d);
    DiracOpGlbSum(&d);
  }
  IFloat *Ap_tmp = (IFloat *)Ap;
  VRB.Flow(cname,fname,"Ap= %e pAp =%e\n",*Ap_tmp,d);

  b = -cp/d;

  z[0][isz] = 1.0;
  z[1][isz] = 1.0;
  bs[isz] = b;
  iz = 1;
  
  for (s=0; s<Nmass; s++) {
    if (s==isz) continue;
    z[1-iz][s] = 1.0;
    z[iz][s] = 1.0 / ( 1.0 - b*(mass[s] - mass[isz]) );
    bs[s] = b*z[iz][s];
  }

  // r[1] += b[0] A.p[0]
  r -> FTimesV1PlusV2(b,Ap,r,f_size);
  // c = |r[1]|^2
  c = r -> NormSqGlbSum(f_size);
  VRB.Flow(cname,fname,"|r[1]^2 =%e\n",c);
  
  // Psi[1] -= b[0] p[0] =- b[0] chi;
  if (type == SINGLE) {
    for (s=0; s<Nmass; s++) {
      b_tmp = bs[s] * alpha[s];
      psi -> FTimesV1PlusV2(-b_tmp,chi,psi,f_size);
    }
  } else {
    for (s=0; s<Nmass; s++) (psi + n_vec*s)-> FTimesV1PlusV2(-bs[s],chi,psi+n_vec*s,f_size);  
  }

  // Check the convergance of the first solution
  for (s=0; s<Nmass; s++) convsP[s] = 0;
  
  convP = (c < rsd_sq[isz]) ? 1 : 0;
  
  // a[k+1] = |r[k]**2/ |r[k-1]|**2
  a = c/cp;
  
  // for k=1 until MaxCG do
  // if |psi[k+1] - psi[k]| <= RsdCG |psi[k+1]| then return
  for (k=1; k<=dirac_arg->max_num_iter && !convP; k++) {
    // a[k+1] = |r[k]**2/ |r[k-1]|**2
    a = c/cp;
  VRB.Flow(cname,fname,"a =%e\n",a);

    // p[k+1] = r[k+1] + a[k+1] p[k]
    //   Compute the shifted as
    //   ps[k+1] = zs[k+1] r[k+1] + a[k+1] ps[k]
    for (s=0; s<Nmass; s++) {
      if (convsP[s]) continue;
      if (s==isz) {
	p[s] -> FTimesV1PlusV2(a,p[s],r,f_size);
      } else {
	as = a * (z[iz][s] * bs[s]) / (z[1-iz][s] * b);
	p[s] -> VecTimesEquFloat(as,f_size);	
	p[s] -> FTimesV1PlusV2(z[iz][s],r,p[s],f_size);	
      }
    IFloat *Ap_tmp = (IFloat *)p[s];
  VRB.Flow(cname,fname,"isz = %d as = %e p[%d] =%e\n",isz, as, s,*Ap_tmp);
    }
    
    // cp = |r[k]**2
    cp = c;
    
    // b[k] = |r[k]**2 / <p[k], Ap[k]>
    if (mass[isz] > 0) {
      MatPcDagMatPc(Ap,p[isz],dot);
      Ap -> FTimesV1PlusV2(mass[isz],p[isz],Ap, f_size);
      d = p[isz] -> ReDotProductGlbSum(Ap, f_size);
    } else {
      MatPcDagMatPc(Ap,p[isz],&d);
      DiracOpGlbSum(&d);
    }
    IFloat *Ap_tmp = (IFloat *)Ap;
  VRB.Flow(cname,fname,"Ap =%e  |b[%d]^2 =%e\n",*Ap_tmp,k,d);

    bp = b;
    b = -cp/d;
    
    //Compute the shifted bs and z
    bs[isz] = b;
    iz = 1 - iz;
    for (s=0; s<Nmass; s++) {
      if (s==isz || convsP[s]) continue;      
      ztmp = z[1-iz][s]*z[iz][s]*bp / 
	( b*a*(z[iz][s]-z[1-iz][s]) + z[iz][s]*bp*(1-b*(mass[s] - mass[isz])));
      bs[s] = b*ztmp / z[1-iz][s];
      z[iz][s] = ztmp;
    }
    
    // r[k+1] += b[k] A.p[k]
    r -> FTimesV1PlusV2(b,Ap,r,f_size);
    // c = |r[k]|**2
    c = r-> NormSqGlbSum(f_size);
    
    // p[k+1] = r[k+1] + a[k+1] p[k]
    //   Compute the shifted as
    //   ps[k+1] = zs[k+1] r[k+1] + a[k+1] ps[k]
    // Psi[k+1] -= b[k] p[k]

    if (type == SINGLE)
      for (s=0; s<Nmass; s++) {
	if (convsP[s]) continue;
	psi->FTimesV1PlusV2(-bs[s]*alpha[s],p[s],psi,f_size);
      }
    else
      for (s=0; s<Nmass; s++) {
	if (convsP[s]) continue;
	(psi + n_vec*s)->FTimesV1PlusV2(-bs[s],p[s],psi+n_vec*s,f_size);
      }
    
    // if |psi[k+1] -psi[k]| <= rsdCG |psi[k+1]| then return
    // or if |r[k+1]| <= RsdCG |chi| then return
    for (s=0; s<Nmass; s++) {
      if (convsP[s]) continue;
      // convergance methods
      
      //check norm of shifted residuals
      css = c * z[iz][s] * z[iz][s];
      convsP[s] = (css < rsd_sq[s]) ? 1 : 0;

      if (convsP[s]){
	RsdCG[s] = sqrt(css);
	VRB.Result(cname,fname,"%d shift converged, iter = %d, res = %e\n",s,k,css);
      }
    }    
    
    convP = convsP[isz];
    // if zero solution has converged, exit unless other solutions have not
    if (convP) for (s=0; s<Nmass; s++) if (!convsP[s]) convP = 0;
  }
  
  // free arrays and vectors
  for (s=0; s<Nmass; s++)
      sfree(*(p+s), cname, fname, "p[s]");
    
  sfree(p, cname, fname, "p");
  sfree(Ap, cname, fname, "Ap");
  sfree(r, cname, fname, "r");
  sfree(bs);
  sfree(*(z+1));
  sfree(*z);
  sfree(z);
  sfree(convsP);
  sfree(rsdcg_sq);
  sfree(rsd_sq);
  
  if (k >= dirac_arg->max_num_iter)   // It has not reached stp_cnd: Issue a warning
    VRB.Warn(cname,fname,"CG reached max iterations = %d. |res|^2 = %e\n",k, css);
  return k;
}

CPS_END_NAMESPACE
