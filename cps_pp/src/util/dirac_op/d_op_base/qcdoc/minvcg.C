#undef PROFILE

#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOpBase class multishift CG solver method.

  $Id: minvcg.C,v 1.23 2008/02/08 18:35:07 chulwoo Exp $
*/

CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include <util/smalloc.h>
#include <qalloc.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <comms/glb.h>
#include <util/verbose.h>
#include <util/error.h>
#include <math.h>
#include <util/checksum.h>
#ifdef HAVE_STRINGS_H
#include <strings.h> // for bzero()
#endif

#define PROFILE

#ifdef PROFILE
#include <util/time_cps.h>
#endif

CPS_START_NAMESPACE


extern "C" { 
  void vaxpy(Float *scale,Vector *mult,Vector *add, int ncvec);
  void vaxpy_norm(Float *scale,Vector *mult,Vector *add, int ncvec, Float *norm);
  void vaxpy_vxdot(Float *scale,Vector *mult,Vector *add, int ncvec, Float *norm);  
}

//! Multishift CG invertor used in RHMC.
int DiracOp::MInvCG(Vector **psi_slow, Vector *chi, Float chi_norm, Float *mass, 
		    int Nmass, int isz, Float *RsdCG, MultiShiftSolveType type,
		    Float *alpha)
{
  char *fname = "MInvCG(V*,V**,...)";
  VRB.Func(cname,fname);
  IFloat *Ap_tmp = (IFloat *)chi;
  VRB.Flow(cname,fname,"chi= %e\n",*Ap_tmp);
  
  // Flash the LED and then turn it off
  //------------------------------------------------------------------
  VRB.LedFlash(cname,fname,3);
  VRB.LedOff(cname,fname);

  if( (lat.Fclass() != F_CLASS_DWF) && (GJP.Snodes()>1) )
    ERR.General(cname,fname,"Fermion class type inconsistent with spread-out S dimension\n");

  // Print out input parameters
  //------------------------------------------------------------------
  VRB.Input(cname,fname,
	    "number of shifts = %d\n",Nmass);
  VRB.Input(cname,fname,
	    "smallest shift stop_rsd = %e\n",IFloat(RsdCG[0]));
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
  if(lat.Fclass() == F_CLASS_CLOVER)
    f_size = lat.FsiteSize()*GJP.VolNodeSites() / 2;
  else
    f_size = lat.FsiteSize()*GJP.VolNodeSites() / (lat.FchkbEvl()+1);

  unsigned long loc_csum = local_checksum((Float *)chi,f_size);
  CSM.SaveCsum(CSUM_EVL_SRC,loc_csum);
  CSM.Clear(CSUM_GLB_LOC);
  CSM.Clear(CSUM_GLB_SUM);

  Vector *r = (Vector *)fmalloc(f_size * sizeof(Float),
				cname,fname, "r");
   
  Vector *Ap = (Vector *)fmalloc(f_size * sizeof(Float),
				 cname,fname, "Ap");
    
  Vector **p = (Vector **)smalloc(Nmass * sizeof(Vector*),
				  cname,fname, "p");

  Vector **psi = (Vector **)smalloc(Nmass * sizeof(Vector*),
				     cname,fname, "psi");
  
  for (s=0; s<Nmass; s++) {
    p[s] = (Vector*)fmalloc(f_size * sizeof(Float),
			    cname,fname, "p[s]");

    if (type == MULTI || s==0) {
      psi[s] = (Vector*)fmalloc(f_size * sizeof(Float),
				cname,fname, "psi[s]");
      if (type == MULTI) bzero((char *)psi[s],sizeof(Float)*f_size);
      else psi[s] -> CopyVec(psi_slow[s],f_size);
    }

  }
  
  int convP;
  int *convsP = (int*)smalloc(cname,fname,"convsP",Nmass*sizeof(int));
  int *converged = (int*)smalloc(cname,fname,"converged",Nmass*sizeof(int));

  Float a, as, b, bp;
  Float *bs = (Float*)smalloc(cname,fname,"bs",Nmass * sizeof(Float));


  Float *z[2];
  for (s=0; s<2; s++) {
    z[s] = (Float*)smalloc(Nmass * sizeof(Float), cname,fname,"z[s]");
  }

  Float css, ztmp;  
  Float c, cp, d;
  
  Float *rsd_sq = (Float*)smalloc(Nmass*sizeof(Float), cname,fname, "rsd_sq");

  Float *rsdcg_sq = (Float*)smalloc(Nmass*sizeof(Float), cname,fname, "rsdcg_sq");
  
  Float *at = (Float*)smalloc(Nmass*sizeof(Float),  cname,fname, "at");

  // code for reproducing CG
  Vector *sol_store;
  int test_num = 0;
  unsigned int *d_store;
  int test_freq = GJP.CGreprodFreq();

  if (test_freq && (CGcount % test_freq == 0) ) {
    test_num = 1;
  VRB.Result(cname,fname,"reproducibility testing interval reachaed\n");
// Allocate space for storing solution
//------------------------------------------------------------------
    sol_store = (Vector *) smalloc(cname,fname,"sol_store", f_size* sizeof(Float));

    sol_store->CopyVec(chi, f_size);


// Allocate space for storing d
//------------------------------------------------------------------
    d_store = (unsigned int *) smalloc(cname,fname, "d_store", (dirac_arg->max_num_iter-1) * sizeof(unsigned int));
  
    for ( int n = 0; n < dirac_arg->max_num_iter-1; n++ )  d_store[n] = 0;

  }   // --- end code for setting up the reproducing CG


  unsigned long x_loc_csum ; 
#ifdef PROFILE
  struct timeval start;
  struct timeval end;
#endif

  for ( int test = 0; test < test_num+1; test++ ) {
    if( test_num != 0 ) 
      {
	VRB.Result(cname,fname,"Repro CG test number = %d\n",test);
      }
    if (test == 1) chi-> CopyVec(sol_store, f_size);

  Float b_tmp;

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
    VRB.Flow(cname,fname,"rsq_sd[%d]= %0.16e\n",s,rsd_sq[s]);
    converged[s] = 0;
  }

  Ap_tmp = (IFloat *)p[0];
  VRB.Flow(cname,fname,"mass[%d]=%e p[%d]= %e\n",0,mass[0],0,*Ap_tmp);
  /*  d = <p, A.p>  */
  if (mass[0] > 0) {
    MatPcDagMatPc(Ap,p[0]);
    Ap_tmp = (IFloat *)Ap;
    VRB.Flow(cname,fname,"Ap= %e\n",*Ap_tmp);
    vaxpy_vxdot(mass,p[0],Ap,f_size/6,&d);
    VRB.Flow(cname,fname,"Ap= %e\n",*Ap_tmp);
  } else {
    MatPcDagMatPc(Ap,p[0],&d);
  }
  x_loc_csum = local_checksum((Float *)Ap,f_size);

  DiracOpGlbSum(&d);
  Ap_tmp = (IFloat *)Ap;
  VRB.Flow(cname,fname,"Ap= %e pAp =%e\n",*Ap_tmp,d);
 
  b = -cp/d;
  
  *(z[0]) = 1.0;
  *(z[1]) = 1.0;
  bs[0] = -b;
  iz = 1;
  
  for (s=0; s<Nmass; s++) {
    if (s==0) continue;
    z[1-iz][s] = 1.0;
    z[iz][s] = 1.0 / ( 1.0 - b*(mass[s] - mass[0]));
    bs[s] = -b*z[iz][s];
  }

  // r[1] += b[0] A.p[0]
  // c = |r[1]|^2
  vaxpy_norm(&b,Ap,r,f_size/6,&c);
  DiracOpGlbSum(&c);
  VRB.Flow(cname,fname, "|res[%d]|^2 = %0.16e\n", 0, IFloat(c));

  // Psi[1] -= b[0] p[0] =- b[0] chi;
  if (type == SINGLE) {
    b_tmp = 0;
    for (s=0; s<Nmass; s++) b_tmp += alpha[s]*bs[s];
    vaxpy(&b_tmp,chi,psi[0],f_size/6);
  } else {
    for (s=0; s<Nmass; s++) vaxpy(bs+s,chi,psi[s],f_size/6);
  }

  // Check the convergance of the first solution
  for (s=0; s<Nmass; s++) {
    convsP[s] = 0;
    at[s] = (Float)1.0;
  }
  
  convP = (c < rsd_sq[0]) ? 1 : 0;

#ifdef PROFILE
  CGflops = 0;
  gettimeofday(&start,NULL);
#endif

  
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
      as = -a * (z[iz][s] * bs[s]) / (z[1-iz][s] * b);
      at[s] *= as;
      as = z[iz][s]/at[s];
      vaxpy(&as,r,p[s],f_size/6);
      CGflops += f_size*2;
    }

    // cp = |r[k]**2
    cp = c;
    
    // b[k] = |r[k]**2 / <p[k], Ap[k]>    
    if (mass[0] > 0) {
      MatPcDagMatPc(Ap,p[0]);
      vaxpy_vxdot(mass,p[0],Ap,f_size/6,&d);
      CGflops += f_size*4;
    } else {
      MatPcDagMatPc(Ap,p[0],&d);
    }
//    x_loc_csum = x_loc_csum ^ local_checksum((Float *)Ap,f_size);
    DiracOpGlbSum(&d);


    if (test_num != 0 ) {
      unsigned int mmp_checksum = local_checksum((Float *)Ap,f_size);
      //    VRB.Result(cname,fname,"DEBUG iter %d checksum = %p\n",k,mmp_checksum);
      // fprintf(stderr, "NODE (%d %d %d %d %d)FAILS checksum[%d] %p\n",
      //        GJP.XnodeCoor(),GJP.YnodeCoor(),GJP.ZnodeCoor(),GJP.TnodeCoor(),GJP.SnodeCoor(),k,mmp_checksum );

      /* Check reproducibility */
      if ( test == 0) d_store[ k ] = mmp_checksum;
      else if ( mmp_checksum != d_store[ k ] ){
        fprintf(stderr, "NODE (%d %d %d %d %d)FAILS TO REPRODUCE\n",
        GJP.XnodeCoor(),GJP.YnodeCoor(),GJP.ZnodeCoor(),GJP.TnodeCoor(),GJP.SnodeCoor());
        fprintf(stderr,"mmp =%p mmp_store = %p\n",mmp_checksum,d_store[k]);
// Temporary hack to exit immediately
        Float *null_p = NULL; *null_p = 0.;
        InterruptExit(-1, "NODE FAILS TO REPRODUCE");
      }
    } // end of repro CG test
    
    bp = b;
    b = -cp/(d*at[0]*at[0]);

    //Compute the shifted bs and z
    bs[0] = -b;
    iz = 1 - iz;
    for (s=0; s<Nmass; s++) {
      if (s==0 || convsP[s]) continue;      
      ztmp = z[1-iz][s]*z[iz][s]*bp / 
	( b*a*(z[iz][s]-z[1-iz][s]) + z[iz][s]*bp*(1-b*(mass[s] - mass[0])));
      bs[s] = -b*ztmp / z[1-iz][s];
      z[iz][s] = ztmp;
    }
    
    // r[k+1] += b[k] A.p[k]
    // c = |r[k]|**2
    b_tmp = b*at[0];
    vaxpy_norm(&b_tmp,Ap,r,f_size/6,&c);
    CGflops += f_size*4;
    DiracOpGlbSum(&c);
    VRB.Flow(cname,fname, "|res[%d]|^2 = %0.16e\n", k, IFloat(c));

    // Psi[k+1] -= b[k] p[k]
    if (type == SINGLE)
      for (s=0; s<Nmass; s++) {
	if (convsP[s]) continue;
	b_tmp = bs[s]*alpha[s]*at[s];
	vaxpy(&b_tmp,p[s],psi[0],f_size/6);
	CGflops += f_size*2;
      }
    else 
      for (s=0; s<Nmass; s++) {
	if (convsP[s]) continue;
	b_tmp = bs[s]*at[s];
	vaxpy(&b_tmp,p[s],psi[s],f_size/6);
	CGflops += f_size*2;
      }    

    // if |psi[k+1] -psi[k]| <= rsdCG |psi[k+1]| then return
    // or if |r[k+1]| <= RsdCG |chi| then return
    for (s=0; s<Nmass; s++) {
      if (convsP[s]) continue;
      //check norm of shifted residuals
      css = c * z[iz][s] * z[iz][s];
      convsP[s] = (css < rsd_sq[s]) ? 1 : 0;
      if (convsP[s]) {
	RsdCG[s] = css;
	converged[s] = k;
      }
    }
    
    convP = convsP[0];
    // if zero solution has converged, exit unless other solutions have not
    if (convP) for (s=0; s<Nmass; s++) if (!convsP[s]) convP = 0;
  }

  } // end of the loop over the repro test

#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(cname,fname,CGflops,&start,&end); 
#endif

  if (k >= dirac_arg->max_num_iter)  // It has not reached stp_cnd: Issue a warning
    ERR.General(cname,fname,"CG reached max iterations = %d. |res| = %e\n",k, css);
    VRB.Result(cname,fname,"iterations = %d. |res| = %e\n",k, css);


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

  loc_csum = 0x0;
  // free arrays and vectors
  for (s=0; s<Nmass; s++) {

    if (type == MULTI || s==0) {
      // Copy solution vectors from fast memory to DDR
      psi_slow[s] -> CopyVec(psi[s],f_size);
      ffree(cname,fname,"psi[s]",psi[s]);
      loc_csum = loc_csum ^ local_checksum((Float *)psi_slow[s],f_size);
    }
    ffree(cname,fname,"p[s]",p[s]);
  }
  CSM.SaveCsum(CSUM_EVL_SOL,loc_csum);
  CSM.SaveCsum(CSUM_MMP_SUM,x_loc_csum);
  CSM.SaveCsumSum(CSUM_GLB_LOC);
  CSM.SaveCsumSum(CSUM_GLB_SUM);
  
  sfree(cname,fname,"p",p);
  sfree(cname,fname,"psi",psi);
  ffree(cname,fname,"Ap",Ap);
  ffree(cname,fname,"r",r);
  sfree(cname,fname,"bs",bs);
  sfree(cname,fname,"z[1]",z[1]);
  sfree(cname,fname,"z[0]",z[0]);
  sfree(cname,fname,"converged",converged);
  sfree(cname,fname,"convsP",convsP);
  sfree(cname,fname,"rsdcg_sq",rsdcg_sq);
  sfree(cname,fname,"rsd_sq",rsd_sq);
  sfree(cname,fname,"at",at);
  
  CGcount++;
  return k;
}

CPS_END_NAMESPACE
