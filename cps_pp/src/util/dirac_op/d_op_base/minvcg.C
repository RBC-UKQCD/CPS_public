#include <config.h>
CPS_START_NAMESPACE
 /*! \file
   \brief  Definition of DiracOpBase class multishift CG solver method.
   
 */

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

//! Multishift CG invertor used in RHMC.
int DiracOp::MInvCG(Vector **psi, Vector *chi, Float chi_norm, Float *mass, int Nmass, 
		    int isz, Float *RsdCG, Vector **EigVec, int NEig)
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
	    "src_norm_sq = %e\n",IFloat(chi_norm*chi_norm));

//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------

  int DeflateReorthFreq = 10000000;
  int n_count;
  int iz, k, s;
  int convP;
  int *convsP = (int*)smalloc(Nmass*sizeof(int));
  
  int f_size;
  if(lat.Fclass() == F_CLASS_CLOVER) {
    f_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  } else {
    f_size = GJP.VolNodeSites() * lat.FsiteSize() / (lat.FchkbEvl()+1);
  }
  
  Vector **p = (Vector **)smalloc(Nmass * sizeof(Vector*));
  if(p == 0) ERR.Pointer(cname,fname, "p");
  VRB.Smalloc(cname,fname, "psi", psi, f_size * sizeof(Vector*));
  for (s=0; s<Nmass; s++) {
    *(p+s) = (Vector*)smalloc(f_size * sizeof(Vector));
    if(*(p+s) == 0) ERR.Pointer(cname,fname, "p[i]");
    VRB.Smalloc(cname,fname,"p[s]", p[s], sizeof(Vector));
  }
  
  Vector *r = (Vector *)smalloc(f_size * sizeof(Vector));
  if(r == 0) ERR.Pointer(cname,fname, "r");
  VRB.Smalloc(cname,fname,"r", r, sizeof(Vector));
  
  Vector *Ap = (Vector *)smalloc(f_size * sizeof(Vector));
  if(Ap == 0) ERR.Pointer(cname,fname, "Ap");
  VRB.Smalloc(cname,fname,"Ap", Ap, sizeof(Vector));
  
  Vector *ltmp = (Vector *)smalloc(f_size * sizeof(Vector));
  if(ltmp == 0) ERR.Pointer(cname,fname, "ltmp");
  VRB.Smalloc(cname,fname,"ltmp", ltmp, sizeof(Vector));
  
  Float a, as, b, bp;
  Float *bs = (Float*)smalloc(Nmass * sizeof(Float));
  Float **z = (Float**)smalloc(2 * sizeof(float*));
  for (s=0; s<2; s++) *(z+s) = (Float*)smalloc(Nmass * sizeof(Float));
  Float css, ztmp;
  
  Float c, cs, d, cp;
  
  Float *rsd_sq = (Float*)smalloc(Nmass*sizeof(Float)); 
  Float *rsdcg_sq = (Float*)smalloc(Nmass*sizeof(Float)); 
  
  Float *dot=0;
  
  // If source norm = 0, solution must be 0
  if (chi_norm == 0.0) {
    for (k=0; k<Nmass; k++) psi[k]->VecTimesEquFloat(0.0,f_size);
    return 0;
  }
  
#ifdef MINV_CHRONO
  if (psi[isz] -> NormSqGlbSum(f_size) == 0.0) {
    r-> CopyVec(chi,f_size);
    for (s=0; s<Nmass; s++) p[s] -> CopyVec(chi,f_size);    
    cp = chi_norm * chi_norm;
    //printf("MInv: zero start\n");
  } else {
    MatPcDagMatPc(Ap,psi[isz],dot);
    Ap ->FTimesV1PlusV2(mass[isz],psi[isz],Ap,f_size);
    r -> FTimesV1MinusV2(1.0,chi,Ap,f_size);
    for (s=0; s<Nmass; s++) p[s] -> CopyVec(r,f_size);
    cp = r -> NormSqGlbSum(f_size);
    //printf("MInv: non-zero start\n");
  }
#else
  // zero psi[0] if not using chronological invertor
  for (s=0; s<Nmass; s++) psi[s] -> VecTimesEquFloat(0.0, f_size);
  r-> CopyVec(chi,f_size);
  for (s=0; s<Nmass; s++) p[s] -> CopyVec(chi,f_size);
  cp = chi_norm*chi_norm;
#endif
  
  for (s=0; s<Nmass; s++) {
    rsdcg_sq[s] = RsdCG[s]*RsdCG[s];
    rsd_sq[s] = cp*rsdcg_sq[s];
  }

  ltmp-> CopyVec(p[isz],f_size);
  MatPcDagMatPc(Ap,ltmp,dot);
  
  // Project out eigen vectors (not necessary for staggered?)
  GramSchm(&Ap, 1, EigVec, NEig, f_size);
  
  /*  d = <p, A.p>  */
  Ap -> FTimesV1PlusV2(mass[isz],p[isz],Ap,f_size);
  d = p[isz] -> ReDotProductGlbSum(Ap, f_size);
  
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
  
  // Psi[1] -= b[0] p[0] =- b[0] chi;
  for (s=0; s<Nmass; s++) psi[s]-> FTimesV1PlusV2(-bs[s],chi,psi[s],f_size);
  
  // c = |r[1]|^2
  
  c = r -> NormSqGlbSum(f_size);
  
  // Check the convergance of the first solution
  for (s=0; s<Nmass; s++) convsP[s] = 0;
  
  css = c; // was double-float conversion
  convP = (css < rsd_sq[isz]) ? 1 : 0;
  
  // for k=1 until MaxCG do
  // if |psi[k+1] - psi[k]| <= RsdCG |psi[k+1]| then return
  for (k=1; k<=dirac_arg->max_num_iter && !convP; k++) {
    
    // a[k+1] = |r[k]**2/ |r[k-1]|**2
    a = c/cp;
    
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
    }
    
    // Project out eigenvectors
    
    if (k% DeflateReorthFreq ==0)
      {
	ltmp->CopyVec(p[isz],f_size);
	GramSchm(&ltmp, 1, EigVec, NEig, f_size);
	p[isz]->CopyVec(ltmp,f_size);
      }
    
    // cp = |r[k]**2
    cp = c;
    
    // b[k] = |r[k]**2 / <p[k], Ap[k]>
    //   First compute d = <p,Ap>
    //   Ap = A. p
    
    ltmp -> CopyVec(p[isz],f_size);
    MatPcDagMatPc(Ap,ltmp,dot);
    
    if (k% DeflateReorthFreq ==0) GramSchm(&Ap,1,EigVec,NEig,f_size);
    
    Ap -> FTimesV1PlusV2(mass[isz],p[isz],Ap,f_size);
    d = p[isz] -> ReDotProductGlbSum(Ap, f_size);
    
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
    
    // Project out eigenvectors
    if (k% DeflateReorthFreq ==0) GramSchm(&r, 1, EigVec, NEig, f_size);
    
    // Psi[k+1] -= b[k] p[k]
    for (s=0; s<Nmass; s++) {
      if (convsP[s]) continue;
      psi[s]->FTimesV1PlusV2(-bs[s],p[s],psi[s],f_size);
    }
    
    // c = |r[k]|**2
    c = r-> NormSqGlbSum(f_size);
    
    // if |psi[k+1] -psi[k]| <= rsdCG |psi[k+1]| then return
    // or if |r[k+1]| <= RsdCG |chi| then return
    for (s=0; s<Nmass; s++) {
      if (convsP[s]) continue;
      // convergance methods
      
#ifndef RELATIVE_ERROR
      //check norm of shifted residuals
      css = c * z[iz][s] * z[iz][s];
      convsP[s] = (css < rsd_sq[s]) ? 1 : 0;
      
#else
      //check relative error of solution
      cs = 0;
      d = 0;
      cs = p[s]-> NormSqGlbSum(p[s], f_size);
      d = psi[s]-> NormSqGlbSum(psi[s], f_size);
      cs*= bs[s]*bs[s];
      d *= rsdcg_sq[s];
      convsP[s] = (cs < d) ? 1 : 0;      
#endif
      //      if (convsP[s]) printf("j = %d, mass = %e converged at %d iterations\n", s, mass[s], k);
    }    
    
    convP = convsP[isz];
    n_count = k;
  }
  
  // Project out eigenvectors
  GramSchm(psi, Nmass, EigVec, NEig, f_size);
  
  // free arrays and vectors
  for (s=0; s<Nmass; s++) {
    VRB.Sfree(cname,fname,"p[s]",*(p+s));
    sfree(*(p+s));
  }
  
  VRB.Sfree(cname,fname,"p",p);
  sfree(p);
  VRB.Sfree(cname,fname,"Ap",Ap);
  sfree(Ap);
  VRB.Sfree(cname,fname,"r",r);
  sfree(r);
  VRB.Sfree(cname,fname,"ltmp",ltmp);
  sfree(ltmp);
  sfree(bs);
  sfree(*(z+1));
  sfree(*z);
  sfree(z);
  sfree(convsP);
  sfree(rsdcg_sq);
  sfree(rsd_sq);
  
  if (n_count == dirac_arg->max_num_iter)   // It has not reached stp_cnd: Issue a warning
    VRB.Warn(cname,fname,"CG reached max iterations = %d. |res|^2 = %e\n",n_count, css);
  return k;
}

CPS_END_NAMESPACE
