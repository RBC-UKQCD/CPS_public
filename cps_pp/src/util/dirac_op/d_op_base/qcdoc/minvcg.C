#define PROFILE

#include <config.h>
CPS_START_NAMESPACE
 /*! \file
   \brief  Definition of DiracOpBase class multishift CG solver method.

   $Id: minvcg.C,v 1.4 2004-06-17 16:21:13 zs Exp $
 */

CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <qalloc.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <comms/glb.h>
#include <util/verbose.h>
#include <util/error.h>
#include <math.h>
#include <util/time.h>
#include <util/enum.h>
CPS_START_NAMESPACE

extern "C" { 
  void vaxpy(Float *scale,Vector *mult,Vector *add, int ncvec);
  void vaxpy_norm(Float *scale,Vector *mult,Vector *add, int ncvec, Float *norm);
  void vaxpy_vxdot(Float *scale,Vector *mult,Vector *add, int ncvec, Float *norm);  
}

//! Multishift CG invertor used in RHMC.
int DiracOp::MInvCG(Vector **psi_slow, Vector *chi, Float chi_norm, Float *mass, int Nmass, 
		    int isz, Float *RsdCG, MultiShiftSolveType type, Float *alpha_slow)
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
  int f_size;
  if(lat.Fclass() == F_CLASS_CLOVER) {
    f_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  } else {
    f_size = GJP.VolNodeSites() * lat.FsiteSize() / (lat.FchkbEvl()+1);
  }
  
  Vector *r = (Vector *)qalloc(QCOMMS|QFAST,f_size * sizeof(Float));
  if(r == 0) r = (Vector *)qalloc(QCOMMS,f_size * sizeof(Float));
  if(r == 0) ERR.Pointer(cname,fname, "r");
  VRB.Smalloc(cname,fname,"r", r, sizeof(Float));
  
  Vector *Ap = (Vector *)qalloc(QCOMMS|QFAST,f_size * sizeof(Float));
  if(Ap == 0) Ap = (Vector *)qalloc(QCOMMS,f_size * sizeof(Float));
  if(Ap == 0) ERR.Pointer(cname,fname, "Ap");
  VRB.Smalloc(cname,fname,"Ap", Ap, sizeof(Float));
  
  Vector **p = (Vector **)qalloc(QCOMMS|QFAST,Nmass * sizeof(Vector*));
  if(p == 0) p = (Vector **)qalloc(QCOMMS,Nmass * sizeof(Vector*));
  if(p == 0) ERR.Pointer(cname,fname, "p");
  VRB.Smalloc(cname,fname, "p", p, f_size * sizeof(Vector*));

  Vector **psi = (Vector **) qalloc(QCOMMS|QFAST,Nmass * sizeof(Vector*));
  if(psi == 0) psi = (Vector **)qalloc(QCOMMS,Nmass * sizeof(Vector*));
  if(psi == 0) ERR.Pointer(cname,fname, "psi");
  VRB.Smalloc(cname,fname, "psi", psi, f_size * sizeof(Vector*));

  if (type == SINGLE) {
    *psi = (Vector*)qalloc(QCOMMS|QFAST,f_size * sizeof(Float));
    if(*psi == 0) *psi = (Vector*)qalloc(QCOMMS,f_size * sizeof(Float));
    if(*psi == 0) ERR.Pointer(cname,fname, "psi[0]");
    VRB.Smalloc(cname,fname,"psi[0]", psi[0], sizeof(Float));
    psi[0] -> CopyVec(psi_slow[0],f_size);
  }

  for (s=0; s<Nmass; s++) {
    *(p+s) = (Vector*)qalloc(QCOMMS|QFAST,f_size * sizeof(Float));
    if(*(p+s) == 0) *(p+s) = (Vector*)qalloc(QCOMMS,f_size * sizeof(Float));
    if(*(p+s) == 0) ERR.Pointer(cname,fname, "p[s]");
    VRB.Smalloc(cname,fname,"p[s]", p[s], sizeof(Float));
    if (type == MULTI) {
      *(psi+s) = (Vector*)qalloc(QCOMMS|QFAST,f_size * sizeof(Float));
      if(*(psi+s) == 0) *(psi+s) = (Vector*)qalloc(QCOMMS,f_size * sizeof(Float));
      if(*(psi+s) == 0) ERR.Pointer(cname,fname, "psi[s]");
      VRB.Smalloc(cname,fname,"psi[s]", psi[s], sizeof(Float));
      psi[s] -> CopyVec(psi_slow[s],f_size);
    }
  }
  
  int convP, converged;
  int *convsP = (int*)qalloc(QCOMMS|QFAST,Nmass*sizeof(int));
  if(convsP == 0) convsP = (int*)qalloc(QCOMMS,Nmass*sizeof(int));
  
  Float a, as, b, bp;
  Float *bs = (Float*)qalloc(QCOMMS|QFAST,Nmass * sizeof(Float));
  if(bs == 0) bs = (Float*)qalloc(QCOMMS,Nmass * sizeof(Float));
  Float **z = (Float**)qalloc(QCOMMS|QFAST,2 * sizeof(Float*));
  if(z == 0) z = (Float**)qalloc(QCOMMS,2 * sizeof(Float*));
  for (s=0; s<2; s++) {
    *(z+s) = (Float*)qalloc(QCOMMS|QFAST,Nmass * sizeof(Float));
    if (*(z+s) == 0) *(z+s) = (Float*)qalloc(QCOMMS,Nmass * sizeof(Float));
  }
  Float css, ztmp;
  
  Float c, cp, d;
  
  Float *rsd_sq = (Float*)qalloc(QCOMMS|QFAST,Nmass*sizeof(Float)); 
  if(rsd_sq == 0) rsd_sq = (Float*)qalloc(QCOMMS,Nmass*sizeof(Float)); 
  Float *rsdcg_sq = (Float*)qalloc(QCOMMS|QFAST,Nmass*sizeof(Float)); 
  if(rsdcg_sq == 0) rsdcg_sq = (Float*)qalloc(QCOMMS,Nmass*sizeof(Float)); 
  
  Float *at = (Float*)qalloc(QCOMMS|QFAST,Nmass*sizeof(Float));
  if (at == 0) at = (Float*)qalloc(QCOMMS,Nmass*sizeof(Float));

  Float *alpha, *b_tmp;
  if (type == SINGLE) {
    alpha = (Float*)qalloc(QCOMMS|QFAST,Nmass*sizeof(Float));
    if (alpha==0) alpha = (Float*)qalloc(QCOMMS|QFAST,Nmass*sizeof(Float));
    for (s=0; s<Nmass; s++) alpha[s] = alpha_slow[s];
  }

  b_tmp = (Float*)qalloc(QCOMMS|QFAST,sizeof(Float));
  if (b_tmp==0) b_tmp = (Float*)qalloc(QCOMMS,sizeof(Float));

  // If source norm = 0, solution must be 0
  if (chi_norm == 0.0) {
    if (type == SINGLE) psi[0]->VecTimesEquFloat(0.0,f_size);
    else for (k=0; k<Nmass; k++) psi[k]->VecTimesEquFloat(0.0,f_size);
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
    MatPcDagMatPc(Ap,p[isz]);
    vaxpy_vxdot(mass+isz,p[isz],Ap,f_size/6,&d);
  } else {
    MatPcDagMatPc(Ap,p[isz],&d);
  }
  DiracOpGlbSum(&d);
 
  b = -cp/d;
  
  z[0][isz] = 1.0;
  z[1][isz] = 1.0;
  bs[isz] = -b;
  iz = 1;
  
  for (s=0; s<Nmass; s++) {
    if (s==isz) continue;
    z[1-iz][s] = 1.0;
    z[iz][s] = 1.0 / ( 1.0 - b*(mass[s] - mass[isz]));
    bs[s] = -b*z[iz][s];
  }

  // r[1] += b[0] A.p[0]
  // c = |r[1]|^2
  vaxpy_norm(&b,Ap,r,f_size/6,&c);
  DiracOpGlbSum(&c);

  // Psi[1] -= b[0] p[0] =- b[0] chi;
  if (type == SINGLE) {
    *b_tmp = 0;
    for (s=0; s<Nmass; s++) *b_tmp += alpha[s]*bs[s];
    vaxpy(b_tmp,chi,psi[0],f_size/6);
  } else {
    for (s=0; s<Nmass; s++) vaxpy(bs+s,chi,psi[s],f_size/6);
  }

  // Check the convergance of the first solution
  for (s=0; s<Nmass; s++) {
    convsP[s] = 0;
    at[s] = (Float)1.0;
  }
  
  convP = (c < rsd_sq[isz]) ? 1 : 0;
  int nflops;
  if(lat.Fclass() == F_CLASS_ASQTAD) {
    nflops = (1170+12*Nmass)*GJP.VolNodeSites();  
  } else if (lat.Fclass() == F_CLASS_STAG) {
    nflops = (594+12*Nmass)*GJP.VolNodeSites();  
  }
  int Nmass_loop = Nmass;

  // for k=1 until MaxCG do
  // if |psi[k+1] - psi[k]| <= RsdCG |psi[k+1]| then return
  for (k=1; k<=dirac_arg->max_num_iter && !convP; k++) {
    timeval start,end;
    // a[k+1] = |r[k]**2/ |r[k-1]|**2
    a = c/cp;
#ifdef PROFILE
    gettimeofday(&start,NULL);
#endif
    // p[k+1] = r[k+1] + a[k+1] p[k]
    //   Compute the shifted as
    //   ps[k+1] = zs[k+1] r[k+1] + a[k+1] ps[k]
    for (s=0; s<Nmass_loop; s++) {
      as = -a * (z[iz][s] * bs[s]) / (z[1-iz][s] * b);
      at[s] *= as;
      as = z[iz][s]/at[s];
      vaxpy(&as,r,p[s],f_size/6);
    }

    // cp = |r[k]**2
    cp = c;
    
    // b[k] = |r[k]**2 / <p[k], Ap[k]>    
    if (mass[isz] > 0) {
      MatPcDagMatPc(Ap,p[isz]);
      vaxpy_vxdot(mass+isz,p[isz],Ap,f_size/6,&d);
    } else {
      MatPcDagMatPc(Ap,p[isz],&d);
    }
    DiracOpGlbSum(&d);
    
    bp = b;
    b = -cp/(d*at[isz]*at[isz]);

    //Compute the shifted bs and z
    bs[isz] = -b;
    iz = 1 - iz;
    for (s=1; s<Nmass_loop; s++) {
      //if (convsP[s]) continue;      
      ztmp = z[1-iz][s]*z[iz][s]*bp / 
	( b*a*(z[iz][s]-z[1-iz][s]) + z[iz][s]*bp*(1-b*(mass[s] - mass[isz])));
      bs[s] = -b*ztmp / z[1-iz][s];
      z[iz][s] = ztmp;
    }
    
    // r[k+1] += b[k] A.p[k]
    // c = |r[k]|**2
    *b_tmp = b*at[isz];
    vaxpy_norm(b_tmp,Ap,r,f_size/6,&c);
    DiracOpGlbSum(&c);

    // Psi[k+1] -= b[k] p[k]
    if (type == SINGLE)
      for (s=0; s<Nmass_loop; s++) {
	*b_tmp = bs[s]*alpha[s]*at[s];
	vaxpy(b_tmp,p[s],psi[0],f_size/6);
      }
    else 
      for (s=0; s<Nmass_loop; s++) {
	*b_tmp = bs[s]*at[s];
	vaxpy(b_tmp,p[s],psi[s],f_size/6);
      }    

#ifdef PROFILE
    gettimeofday(&end,NULL);
    printf("DiracOpAsqtad::MInvCG:: ");
    print_flops(nflops,&start,&end); 
#endif

    // if |psi[k+1] -psi[k]| <= rsdCG |psi[k+1]| then return
    // or if |r[k+1]| <= RsdCG |chi| then return
    converged = 0;
    for (s=0; s<Nmass_loop; s++) {
      //if (convsP[s]) continue;
      //check norm of shifted residuals
      css = c * z[iz][s] * z[iz][s];
      convsP[s] = (css < rsd_sq[s]) ? 1 : 0;
      if (convsP[s]) {
	converged++;
	nflops -= 12*GJP.VolNodeSites();
	VRB.Result(cname,fname,"%d shift converged, iter = %d, res = %e\n",s,k,css);
      }
    }
    Nmass_loop -= converged;
    
    convP = convsP[isz];
    // if zero solution has converged, exit unless other solutions have not
    //if (convP) for (s=0; s<Nmass; s++) if (!convsP[s]) convP = 0;
  }

  // free arrays and vectors
  for (s=0; s<Nmass; s++) {
    if (type == MULTI || s==0) {
      // Copy solution vectors from fast memory to DDR
      psi_slow[s] -> CopyVec(psi[s],f_size);
      VRB.Sfree(cname,fname,"psi[s]",*(psi+s));
      qfree(*(psi+s));
    }
    VRB.Sfree(cname,fname,"p[s]",*(p+s));
    qfree(*(p+s));
  }
  
  VRB.Sfree(cname,fname,"p",p);
  qfree(p);
  VRB.Sfree(cname,fname,"psi",psi);
  qfree(psi);
  VRB.Sfree(cname,fname,"Ap",Ap);
  qfree(Ap);
  VRB.Sfree(cname,fname,"r",r);
  qfree(r);
  qfree(bs);
  qfree(b_tmp);
  qfree(*(z+1));
  qfree(*z);
  qfree(z);
  qfree(convsP);
  qfree(rsdcg_sq);
  qfree(rsd_sq);
  qfree(at);
  if (type == SINGLE) qfree(alpha);
  
  if (k >= dirac_arg->max_num_iter)   // It has not reached stp_cnd: Issue a warning
    VRB.Warn(cname,fname,"CG reached max iterations = %d. |res|^2 = %e\n",k, css);
  return k;
}

CPS_END_NAMESPACE
