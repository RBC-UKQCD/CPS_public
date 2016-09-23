#define PROFILE

#include <stdio.h>
#include <config.h>
CPS_START_NAMESPACE
 /*! \file
   \brief  Definition of DiracOp class multishift CG solver method.

  $Id: minvcg.C,v 1.6 2012/03/26 13:50:11 chulwoo Exp $
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

#ifdef PROFILE
#include <util/time_cps.h>
#endif

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
  \param isz The relative location of the smallest shift
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
int DiracOp::MInvCG(Vector **psi, Vector *chi, Float chi_norm, Float *mass, 
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
  
  if( (lat.Fclass() != F_CLASS_DWF) && (GJP.Snodes()>1) )
    ERR.General(cname,fname,"Fermion class type inconsistent with spread-out S dimension\n");


// Print out input parameters
//------------------------------------------------------------------
  VRB.Result(cname,fname,
	    "number of shifts = %d\n",Nmass);
  VRB.Result(cname,fname,
	    "smallest shift stop_rsd = %e\n",IFloat(RsdCG[0]));
  VRB.Result(cname,fname,
	    "max_num_iter = %d\n",dirac_arg->max_num_iter);
  VRB.Result(cname,fname,
	    "mass = %e\n",IFloat(dirac_arg->mass));
  VRB.Result(cname,fname,
	    "src_norm_sq = %e\n",IFloat(chi_norm));

//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------

  int iz, k, s;
  int f_size;

  if(lat.Fclass() == F_CLASS_CLOVER)
    f_size = lat.FsiteSize()*GJP.VolNodeSites() / 2;
  else
    f_size = lat.FsiteSize()*GJP.VolNodeSites() / (lat.FchkbEvl()+1);

  Vector *r = (Vector*)smalloc(f_size*sizeof(Float),cname,fname,"r");
  
  Vector *Ap = (Vector*)smalloc(f_size*sizeof(Float),cname,fname,"Ap");
  
  Vector **p = (Vector**)smalloc(Nmass*sizeof(Vector*),cname,fname,"p");

  for (s=0; s<Nmass; s++) 
    *(p+s) = (Vector*)smalloc(f_size * sizeof(Float),
			      cname,fname, "p[i]");
    
  int convP;
  int *convsP = (int*)smalloc(Nmass*sizeof(int));
  int *converged = (int*)smalloc(Nmass*sizeof(int));
  
  Float a=0, as, b, bp, b_tmp;
  Float *bs = (Float*)smalloc(Nmass * sizeof(Float));
  Float **z = (Float**)smalloc(2 * sizeof(Float*));
  for (s=0; s<2; s++) *(z+s) = (Float*)smalloc(Nmass * sizeof(Float));
  Float css, ztmp;
  
  Float c, cs, d, cp;
  Float *dot=0;

  Float *rsd_sq = (Float*)smalloc(Nmass*sizeof(Float)); 
  Float *rsdcg_sq = (Float*)smalloc(Nmass*sizeof(Float)); 

  if (type == MULTI) {
    for (int i=0; i<Nmass; i++)
      psi[i] -> VecZero(f_size);
  }
  
  if (type == SINGLE || type == MULTI) {
    r-> CopyVec(chi,f_size);
    cp = chi_norm;
  } else if (type == GENERAL) {
    MatPcDagMatPc(r,psi[0]);
    r -> FTimesV1MinusV2(1.0,chi,r,f_size);
    cp = r -> NormSqGlbSum(f_size);
  }

  for (s=0; s<Nmass; s++) p[s] -> CopyVec(r,f_size);

  for (s=0; s<Nmass; s++) {
    rsdcg_sq[s] = RsdCG[s]*RsdCG[s];
    rsd_sq[s] = cp*rsdcg_sq[s];
    converged[s] = 0;
  }

  /*  d = <p, A.p>  */
  if (mass[0] > 0) {
    MatPcDagMatPc(Ap,p[0]);
    Ap -> FTimesV1PlusV2(mass[0],p[0],Ap, f_size);
    d = p[0] -> ReDotProductGlbSum(Ap, f_size);
  } else {
    MatPcDagMatPc(Ap,p[0],&d);
    DiracOpGlbSum(&d);
  }
  IFloat *Ap_tmp = (IFloat *)Ap;
  VRB.Flow(cname,fname,"Ap= %e pAp =%e\n",*Ap_tmp,d);

  b = -cp/d;

  z[0][0] = 1.0;
  z[1][0] = 1.0;
  bs[0] = b;
  iz = 1;
  
  for (s=0; s<Nmass; s++) {
    if (s==0) continue;
    z[1-iz][s] = 1.0;
    z[iz][s] = 1.0 / ( 1.0 - b*(mass[s] - mass[0]) );
    bs[s] = b*z[iz][s];
  }

  // r[1] += b[0] A.p[0]
  r -> FTimesV1PlusV2(b,Ap,r,f_size);
  // c = |r[1]|^2
  c = r -> NormSqGlbSum(f_size);
  VRB.Flow(cname,fname,"|r[1]|^2 =%e\n", c);

  // Psi[1] -= b[0] p[0] =- b[0] chi;
  if (type == SINGLE) {
    for (s=0; s<Nmass; s++) {
      b_tmp = bs[s] * alpha[s];
      psi[0] -> FTimesV1PlusV2(-b_tmp,chi,psi[0],f_size);
    }
  } else {
    for (s=0; s<Nmass; s++) 
      psi[s]-> FTimesV1PlusV2(-bs[s],chi,psi[s],f_size);  
  }

  // Check the convergance of the first solution
  for (s=0; s<Nmass; s++) convsP[s] = 0;
  
  convP = (c < rsd_sq[0]) ? 1 : 0;
  
  // a[k+1] = |r[k]**2/ |r[k-1]|**2
  a = c/cp;
  
#ifdef PROFILE
  struct timeval start;
  struct timeval end;
  CGflops = 0;
  gettimeofday(&start,NULL);
#endif

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
      if (s==0) {
	p[s] -> FTimesV1PlusV2(a,p[s],r,f_size);
	CGflops += f_size*2;
      } else {
	as = a * (z[iz][s] * bs[s]) / (z[1-iz][s] * b);
	p[s] -> VecTimesEquFloat(as,f_size);	
	p[s] -> FTimesV1PlusV2(z[iz][s],r,p[s],f_size);	
	CGflops += f_size*3;
      }
    }
    
    // cp = |r[k]**2
    cp = c;
    
    // b[k] = |r[k]**2 / <p[k], Ap[k]>
    if (mass[0] > 0) {
      MatPcDagMatPc(Ap,p[0],dot);
      Ap -> FTimesV1PlusV2(mass[0],p[0],Ap, f_size);
      d = p[0] -> ReDotProductGlbSum(Ap, f_size);
      CGflops += f_size*4;
    } else {
      MatPcDagMatPc(Ap,p[0],&d);
      DiracOpGlbSum(&d);
    }

    bp = b;
    b = -cp/d;
    
    //Compute the shifted bs and z
    bs[0] = b;
    iz = 1 - iz;
    for (s=0; s<Nmass; s++) {
      if (s==0 || convsP[s]) continue;      
      ztmp = z[1-iz][s]*z[iz][s]*bp / 
	( b*a*(z[iz][s]-z[1-iz][s]) + z[iz][s]*bp*(1-b*(mass[s] - mass[0])));
      bs[s] = b*ztmp / z[1-iz][s];
      z[iz][s] = ztmp;
    }
    
    // r[k+1] += b[k] A.p[k]
    r -> FTimesV1PlusV2(b,Ap,r,f_size);
    // c = |r[k]|**2
    c = r-> NormSqGlbSum(f_size);
    CGflops += f_size*4;    

    // p[k+1] = r[k+1] + a[k+1] p[k]
    //   Compute the shifted as
    //   ps[k+1] = zs[k+1] r[k+1] + a[k+1] ps[k]
    // Psi[k+1] -= b[k] p[k]

    if (type == SINGLE)
      for (s=0; s<Nmass; s++) {
	if (convsP[s]) continue;
	VRB.Result(cname,fname,"bs[%d]=%g psi[%d]=%g\n",s,bs[s],s,psi[s]);
	psi[0]->FTimesV1PlusV2(-bs[s]*alpha[s],p[s],psi[0],f_size);
	CGflops += f_size*2;
      }
    else
      for (s=0; s<Nmass; s++) {
	if (convsP[s]) continue;
	VRB.Result(cname,fname,"bs[%d]=%g psi[%d]=%g\n",s,bs[s],s,psi[s]);
	psi[s]->FTimesV1PlusV2(-bs[s],p[s],psi[s],f_size);
	CGflops += f_size*2;
      }
    
    // if |psi[k+1] -psi[k]| <= rsdCG |psi[k+1]| then return
    // or if |r[k+1]| <= RsdCG |chi| then return
    for (s=0; s<Nmass; s++) {
      if (convsP[s]) continue;
      //check norm of shifted residuals
      css = c * z[iz][s] * z[iz][s];
      convsP[s] = (css < rsd_sq[s]) ? 1 : 0;
      if (convsP[s]){
	RsdCG[s] = css;
	converged[s] = k;
      }
    }    
    
    convP = convsP[0];
    // if zero solution has converged, exit unless other solutions have not
    if (convP) for (s=0; s<Nmass; s++) if (!convsP[s]) convP = 0;
  }
  
#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(cname,fname,CGflops,&start,&end); 
#endif

  if (k >= dirac_arg->max_num_iter)   // It has not reached stp_cnd: Issue a warning
    VRB.Warn(cname,fname,"CG reached max iterations = %d. |res|^2 = %e\n",k, css);

  for (s=Nmass-1; s>=0; s--) {
    if (convsP[s]) {
      VRB.Result(cname,fname,"%d shift converged, iter = %d, res^2 = %e\n",
		 s+isz,converged[s],RsdCG[s]);
      RsdCG[s] = sqrt(RsdCG[s]);
    } else {
      RsdCG[s] = c*z[iz][s]*z[iz][s];
      VRB.Result(cname,fname,
		 "%d shift did not converge, iter = %d, res^2 = %e\n",
		 s+isz,k,RsdCG[s]);
      RsdCG[s] = sqrt(RsdCG[s]);
    }
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
  sfree(converged);
  sfree(convsP);
  sfree(rsdcg_sq);
  sfree(rsd_sq);
  
  return k;
}

CPS_END_NAMESPACE
