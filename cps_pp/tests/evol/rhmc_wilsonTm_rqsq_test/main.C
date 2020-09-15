#include<config.h>

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#include<util/lattice.h>
#include<util/lattice/fbfm.h>
#include<util/random.h>
#include<util/time_cps.h>

#include<alg/alg_hmc.h>
#include<alg/bfm_arg.h>
#include<alg/common_arg.h>
#include<alg/hmc_arg.h>
#include<alg/hmd_arg.h>

#include<alg/alg_int.h>
#include<alg/int_arg.h>
#include<alg/alg_wline.h>
#include<alg/no_arg.h>
#include<alg/do_arg.h>
#include<alg/alg_plaq.h>
#include<alg/alg_pbp.h>
#include<alg/alg_tcharge.h>
#include<alg/alg_smear.h>
#include<alg/ape_smear_arg.h>
#include<alg/pbp_arg.h>
#include<alg/alg_remez.h>
#include<alg/alg_w_spect.h>
#include<alg/array_arg.h>
#include<alg/alg_fix_gauge.h>

#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/qcdio.h>
#include<util/WriteLatticePar.h>
#include<util/ReadLatticePar.h>
#include<util/qioarg.h>
#include<util/dirac_op.h>

#include <util/lat_cont.h>

#include <chroma.h>
#include <omp.h>
#include <pthread.h>
//-------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

const char *cname = "";

HmcArg hmc_arg;

ActionGaugeArg gauge_arg;
ActionRationalQuotientArg rat_quo_arg;
ActionRationalQuotientArg rat_quo_arg_2;

// define all integrators; vml file arguments will not change
IntABArg ab1_arg;
IntABArg ab2_arg;

EvoArg evo_arg;
DoArg do_arg;
PbpArg pbp_arg;
NoArg no_arg;

void checkpoint(int traj);

#define decode_vml(arg_name)\
  printf("Decoding %s\n",#arg_name);\
  do{                                       \
        if ( ! arg_name.Decode(#arg_name".vml", #arg_name) )            \
            ERR.General(cname, fname, "Bad " #arg_name ".vml.\n");      \
    } while(0)

void decode_vml_all(void)
{
    char *fname = "decode_vml_all()";

    decode_vml(do_arg);
    decode_vml(hmc_arg);
    decode_vml(evo_arg);
    decode_vml(gauge_arg);
    decode_vml(rat_quo_arg);
    decode_vml(rat_quo_arg_2);
    decode_vml(ab1_arg);
    decode_vml(ab2_arg);
    decode_vml(pbp_arg);
}

void truncate_it(CommonArg *common_arg, const char stem[], int traj);
void measure_plaq(CommonArg &common_arg);
void measure_pbp(CommonArg &common_arg, int traj);
void run_hmc(CommonArg &common_arg, int traj, AlgIntAB &int_ab);


void setup(int argc, char *argv[])
{
    const char *fname = "setup()";

    Start(&argc, &argv);

    if(argc < 2) {
        ERR.General(cname, fname, "Must provide VML directory.\n");
    }

    if(chdir(argv[1]) != 0) {
        ERR.General(cname, fname, "Changing directory to %s failed.\n", argv[1]);
    }

    decode_vml_all();

    if(chdir(evo_arg.work_directory) != 0) {
        ERR.General(cname, fname, "Changing directory to %s failed.\n", evo_arg.work_directory);
    }
    VRB.Result(cname, fname, "Reading VML files successfully.\n");

    GJP.Initialize(do_arg);
    //LRG.Initialize();

    VRB.Result(cname, fname, "VRB.Level(%d)\n", do_arg.verbose_level);
    VRB.Level(do_arg.verbose_level);
}

Float* rand_4d_canonical_fermion(Lattice &lat){
  long f_size = (long)24 * GJP.VolNodeSites();
  Float *v1 = (Float *)pmalloc(sizeof(Float) * f_size);
  printf("Making random gaussian 4d vector\n");
  lat.RandGaussVector((Vector*)v1, 0.5, 2, CANONICAL, FOUR_D);
  printf("Finished making random gaussian vector\n");
  return v1;
}

int MInvCG_CPS(Vector **psi, Vector *chi, Float chi_norm, Float *mass, 
	       int Nmass, int isz, Float *RsdCG,
	       MultiShiftSolveType type, Float *alpha, DiracOp &dop, Lattice &lat, CgArg &cg_arg)
{
  char *fname = "MInvCG(V*,V**,...) [Duplicate of d_op_base/noarch version]";
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
	    "max_num_iter = %d\n",cg_arg.max_num_iter);
  VRB.Result(cname,fname,
	    "mass = %e\n",IFloat(cg_arg.mass));
  VRB.Result(cname,fname,
	    "src_norm_sq = %e\n",IFloat(chi_norm));

//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------

  int iz, k, s;
  size_t f_size;

  if(lat.Fclass() == F_CLASS_CLOVER)
    f_size = lat.FsiteSize()*GJP.VolNodeSites() / 2;
  else
    f_size = lat.FsiteSize()*GJP.VolNodeSites() / (lat.FchkbEvl()+1);

  if(GJP.Gparity()) f_size*=2;

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
    dop.MatPcDagMatPc(r,psi[0]);
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
    dop.MatPcDagMatPc(Ap,p[0]); //Ap = Mpc^dag Mpc p
    Ap -> FTimesV1PlusV2(mass[0],p[0],Ap, f_size); //Ap = mass[0]*p + Ap
    d = p[0] -> ReDotProductGlbSum(Ap, f_size); //d = p.Ap =  p. (Mpc^dag Mpc + mass[0])p
  } else { //CK: Why do we do something different for negative mass??
    dop.MatPcDagMatPc(Ap,p[0],&d);
    glb_sum(&d);
  }
  IFloat *Ap_tmp = (IFloat *)Ap;
  VRB.Flow(cname,fname,"Ap= %e pAp =%e\n",*Ap_tmp,d);

  b = -cp/d;

  VRB.Flow(cname,fname,"b = -cp/d = -%e/%e =%e\n",cp,d,b);

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
  
  // for k=1 until MaxCG do
  // if |psi[k+1] - psi[k]| <= RsdCG |psi[k+1]| then return
  for (k=1; k<=cg_arg.max_num_iter && !convP; k++) {
    // a[k+1] = |r[k]**2/ |r[k-1]|**2
    a = c/cp;
    VRB.Flow(cname,fname,"a =%e, |r[%d]]^2 = %e\n",a,k,c);

    // p[k+1] = r[k+1] + a[k+1] p[k]
    //   Compute the shifted as
    //   ps[k+1] = zs[k+1] r[k+1] + a[k+1] ps[k]
    for (s=0; s<Nmass; s++) {
      if (convsP[s]) continue;
      if (s==0) {
	p[s] -> FTimesV1PlusV2(a,p[s],r,f_size);
      } else {
	as = a * (z[iz][s] * bs[s]) / (z[1-iz][s] * b);
	p[s] -> VecTimesEquFloat(as,f_size);	
	p[s] -> FTimesV1PlusV2(z[iz][s],r,p[s],f_size);	
      }
    }
    
    // cp = |r[k]**2
    cp = c;
    
    // b[k] = |r[k]**2 / <p[k], Ap[k]>
    if (mass[0] > 0) {
      dop.MatPcDagMatPc(Ap,p[0],dot);
      Ap -> FTimesV1PlusV2(mass[0],p[0],Ap, f_size);
      d = p[0] -> ReDotProductGlbSum(Ap, f_size);
    } else {
      dop.MatPcDagMatPc(Ap,p[0],&d);
      glb_sum(&d);
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

    // p[k+1] = r[k+1] + a[k+1] p[k]
    //   Compute the shifted as
    //   ps[k+1] = zs[k+1] r[k+1] + a[k+1] ps[k]
    // Psi[k+1] -= b[k] p[k]

    if (type == SINGLE)
      for (s=0; s<Nmass; s++) {
	if (convsP[s]){ 
	Float *tmp_p = (Float *)psi[0];
	VRB.Result(cname,fname,"bs[%d]=%g psi[%d]=%g\n",s,bs[s],s,*tmp_p);
 	continue;}
	psi[0]->FTimesV1PlusV2(-bs[s]*alpha[s],p[s],psi[0],f_size);
      }
    else
      for (s=0; s<Nmass; s++) {
	if (convsP[s]){ 
 	continue;}
	psi[s]->FTimesV1PlusV2(-bs[s],p[s],psi[s],f_size);
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
  
  if (k >= cg_arg.max_num_iter)   // It has not reached stp_cnd: Issue a warning
    VRB.Warn(cname,fname,"CG reached max iterations = %d. |res|^2 = %e\n",k, css);

  for (s=Nmass-1; s>=0; s--) {
    if (convsP[s]) {
      Float res2 = RsdCG[s];
      RsdCG[s] = sqrt(RsdCG[s]);
      VRB.Result(cname,fname,"%d shift converged, iter = %d, res^2 = %e, |res| = %e\n",
		 s+isz,converged[s],res2,RsdCG[s]);

      //CK: add calculation of true residual
      Vector* mmp = (Vector*)pmalloc(f_size*sizeof(Float));
      dop.MatPcDagMatPc(mmp, psi[s]);
      mmp -> FTimesV1PlusV2(mass[s],psi[s],mmp, f_size);
      
      Float src_norm_sq = chi->NormSqNode(f_size);
      glb_sum(&src_norm_sq);

      Vector* res = (Vector*)pmalloc(f_size*sizeof(Float));
      res->CopyVec(chi, f_size);
      res->VecMinusEquVec(mmp, f_size);
      Float res_norm_sq_cur = res->NormSqNode(f_size);
      glb_sum(&res_norm_sq_cur);
      Float tmp = res_norm_sq_cur / src_norm_sq;
      tmp = sqrt(tmp);
      VRB.Result(cname,fname,
		 "True %d: |res| / |src| = %e\n", s,IFloat(tmp));
      pfree(mmp);
      pfree(res);

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

int InvCgShift_CPS(Vector *out, 
		   Vector *in, 
		   Float src_norm_sq, 
		   Float *true_res,
		   Float *shift,
		   DiracOp &dop,
		   Lattice& lat, CgArg &cg_arg){
  int itr;                       // Current number of CG iterations
  int max_itr;                       // Max number of CG iterations
  Float stp_cnd;                   // Stop if residual^2 <= stp_cnd
  Float res_norm_sq_prv;          // The previous step |residual|^2
  Float res_norm_sq_cur;           // The current step |residual|^2 
  Float a;
  Float b;
  Float d;
  int i, ic, icb;
  const char *fname = "InvCgShift(V*,V*,F,F*) [Duplicate of d_op_base/noarch version]";
  IFloat *temp;


// Print out input parameters
//------------------------------------------------------------------
  VRB.Input(cname,fname,
	    "stop_rsd = %e\n",IFloat(cg_arg.stop_rsd));
  VRB.Input(cname,fname,
	    "max_num_iter = %d\n",cg_arg.max_num_iter);
  VRB.Input(cname,fname,
	    "mass = %e\n",IFloat(cg_arg.mass));
  VRB.Input(cname,fname,
	    "src_norm_sq = %e\n",IFloat(src_norm_sq));


//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------

// Set the source vector pointer
//------------------------------------------------------------------
  Vector *src = in;

// Set the solution vector pointer
//------------------------------------------------------------------
  Vector *sol = out;

// Set the node checkerboard size of the fermion field
//------------------------------------------------------------------

  size_t f_size_cb;

  if(lat.Fclass() == F_CLASS_CLOVER) {
    f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  } else {
    f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / (lat.FchkbEvl()+1);
  }
    
  if(GJP.Gparity()) f_size_cb*=2;

  Vector *res = (Vector *) smalloc(f_size_cb * sizeof(Float));
  Vector *dir = (Vector *) smalloc(f_size_cb * sizeof(Float));
  Vector *mmp = (Vector *) smalloc(f_size_cb * sizeof(Float));

// If src_norm_sq is not provided calculate it
//------------------------------------------------------------------
  if(src_norm_sq == 0){
    src_norm_sq = src->NormSqNode(f_size_cb); //CK: in G-parity situation we want the norm^2 of the whole 2-flavour double-wrapped source
    glb_sum(&src_norm_sq);
  }
  VRB.Flow(cname,fname,"src_norm_sq=%e\n",src_norm_sq);

// Calculate stopping condition
//------------------------------------------------------------------
  stp_cnd = src_norm_sq * cg_arg.stop_rsd * cg_arg.stop_rsd;
  VRB.Flow(cname,fname, 
	   "stp_cnd =%e\n", IFloat(stp_cnd));

// Make IFloat pointers out of Vector pointers
//------------------------------------------------------------------
  IFloat *f_sol = (IFloat *) sol; 
  IFloat *f_dir = (IFloat *) dir; 
  IFloat *f_res = (IFloat *) res; 
  IFloat *f_mmp = (IFloat *) mmp; 

//------------------------------------------------------------------
// Initial step:
// res = src - MatPcDagMatPc * sol
// dir = res
// if( |res|^2 <= stp_cnd ){ 
//   n_count = 0
//   free memory
//   return
// }
//------------------------------------------------------------------
  Float *in_f =  (Float *) sol;
  // Mmp = MatPcDagMatPc * sol
  dop.MatPcDagMatPc(mmp, sol);
  if (shift){
    mmp -> FTimesV1PlusV2(*shift,sol,mmp, f_size_cb);
  }
  //print_vec( mmp, "mmp");

  // res = src
  res->CopyVec(src, f_size_cb);
  //print_vec( res, "res");

  // res -= mmp
  res->VecMinusEquVec(mmp, f_size_cb);
  //print_vec( res, "res");

  // dir = res
  dir->CopyVec(res, f_size_cb);  
  //print_vec( dir, "dir");

  // res_norm_sq_cur = res * res
  res_norm_sq_cur = res->NormSqNode(f_size_cb);
  //printf("res_norm_sq_cur=%e\n",res_norm_sq_cur);
  glb_sum(&res_norm_sq_cur);

  // if( |res|^2 <= stp_cnd ) we are done
  VRB.Flow(cname,fname,
  	   "|res[0]|^2 = %e\n", IFloat(res_norm_sq_cur));
  itr = 0;
  max_itr = 9999;
  if(res_norm_sq_cur <= stp_cnd) max_itr = 0;
  //printf("max_itr=%d\n",max_itr);


//------------------------------------------------------------------
// Loop over CG iterations
//------------------------------------------------------------------

  for(i=0; i < max_itr; i++){
    itr = itr + 1;
    res_norm_sq_prv = res_norm_sq_cur;

    // mmp = MatPcDagMatPc * dir
    // d = <dir, MatPcDagMatPc*dir>

    dop.MatPcDagMatPc(mmp, dir, &d);
    if (shift){
      mmp -> FTimesV1PlusV2(*shift,dir,mmp, f_size_cb);
      Float dir_sq = dir -> NormSqNode(f_size_cb);
      glb_sum(&dir_sq);
      d += (*shift) * dir_sq;
    }
    //printf("d=%e\n",d);
    //print_vec( mmp, "mmp");
  
    glb_sum(&d);
    VRB.Flow(cname,fname, "d = %e\n", IFloat(d));

    // If d = 0 we are done
    if(d == 0.0) {
      VRB.Warn(cname,fname,"d(%e) = 0.0!!\n",d);
      //	exit(5);
      break;
      //??? or should we give a warning or error? Yes we should, really.
    }

    a = res_norm_sq_prv / d;
    VRB.Flow(cname,fname, "a = %e\n", IFloat(a));

    // Set circular buffer
    //    setCbufCntrlReg(4, CBUF_MODE4);

    // sol = a * dir + sol;
    sol->FTimesV1PlusV2(a, dir, sol, f_size_cb);
    //print_vec( sol, "sol");

    // res = - a * (MatPcDagMatPc * dir) + res;
    res->FTimesV1PlusV2(-a, mmp, res, f_size_cb);
    //print_vec( res, "res");

    // res_norm_sq_cur = res * res
    res_norm_sq_cur = res->NormSqNode(f_size_cb);
    glb_sum(&res_norm_sq_cur);

    // if( |res|^2 <= stp_cnd ) we are done
    VRB.Flow(cname,fname,
	     "|res[%d]|^2 = %e\n", itr, IFloat(res_norm_sq_cur));
    if(res_norm_sq_cur <= stp_cnd) break;

    b = res_norm_sq_cur / res_norm_sq_prv;
    VRB.Flow(cname,fname, "b = %e\n", IFloat(b));

    // dir = b * dir + res;
    dir->FTimesV1PlusV2(b, dir, res, f_size_cb);
    //print_vec( dir, "dir");
  }

  // It has not reached stp_cnd: Issue a warning
  if(itr == cg_arg.max_num_iter - 1){
    VRB.Warn(cname,fname,
	      "CG reached max iterations = %d. |res|^2 = %e\n",
	     itr+1, IFloat(res_norm_sq_cur) );
  }

//------------------------------------------------------------------
// Done. Finish up and return
//------------------------------------------------------------------
  // Calculate and set true residual: 
  // true_res = |src - MatPcDagMatPc * sol| / |src|
  dop.MatPcDagMatPc(mmp, sol);
  if (shift){
    mmp -> FTimesV1PlusV2(*shift,sol,mmp, f_size_cb);
  }
  res->CopyVec(src, f_size_cb);
  res->VecMinusEquVec(mmp, f_size_cb);
  res_norm_sq_cur = res->NormSqNode(f_size_cb);
  glb_sum(&res_norm_sq_cur);
  Float tmp = res_norm_sq_cur / src_norm_sq;
  tmp = sqrt(tmp);
  if(true_res != 0){
    *true_res = tmp;
  }
  VRB.Flow(cname,fname,
	     "True |res| / |src| = %e, iter = %d\n", IFloat(tmp), itr+1);
  

  // Free memory
  VRB.Sfree(cname,fname, "mmp", mmp);
  sfree(mmp);
  VRB.Sfree(cname,fname, "dir", dir);
  sfree(dir);
  VRB.Debug("b ============\n");
  VRB.Sfree(cname,fname, "res", res);
  sfree(res);

  VRB.Debug("a ============\n");

  // Return number of iterations
  return itr+1;
}

int InvCg_CPS(Vector *out, 
	      Vector *in, 
	      Float src_norm_sq, 
	      Float *true_res,
	      DiracOp &dop,
	      Lattice& lat, CgArg &cg_arg){
  return InvCgShift_CPS(out,in,src_norm_sq,true_res,NULL,dop,lat,cg_arg);
}



int MinvCGtest(){
  //Test the BFM multi-mass solver

  VRB.Level(VERBOSE_FLOW_LEVEL);

  GnoneFwilsonTm* lattice = new GnoneFwilsonTm;

  Float targ_resid = 1e-8;

  CgArg cg_arg;
  cg_arg.mass = -1.8;
  cg_arg.epsilon = 0.5;
  cg_arg.max_num_iter = 10000;
  cg_arg.stop_rsd = targ_resid;
  cg_arg.true_rsd = targ_resid;
  cg_arg.Inverter = CG;
  
  DiracOpWilsonTm dop(*lattice, (Vector*)0, (Vector*)0, &cg_arg, CNV_FRM_NO);

  //technically we need a WILSON ordered fermion, but as it is random it doesn't matter
  Float* in = rand_4d_canonical_fermion(*lattice);

  size_t f_size_cb =  GJP.VolNodeSites() * 24/2;
  Float src_norm_sq = ((Vector*)in)->NormSqNode(f_size_cb);
  glb_sum(&src_norm_sq);
  
  Vector **out_1 = (Vector**)pmalloc(3 * sizeof(Vector*));
  Vector **out_2 = (Vector**)pmalloc(3 * sizeof(Vector*));
  
  int nmass =3;

  //Float mass[3] = {-1.8, -1.0, -0.4};
  Float mass[3] = {0.1, 0.2, 0.3};

  for(int j=0;j<3;j++){
    out_1[j] = (Vector*)pmalloc(f_size_cb*sizeof(Float));
    out_2[j] = (Vector*)pmalloc(f_size_cb*sizeof(Float));

    for(int i=0;i<f_size_cb;i++){
      ((Float*)out_1[j])[i] = 0.0; ((Float*)out_2[j])[i] = 0.0;
    }
  }

  Float true_rsd[3] = {targ_resid, targ_resid, targ_resid};

  MInvCG_CPS((Vector**)out_2, (Vector*)in, src_norm_sq, &mass[0], nmass, 0, &true_rsd[0], MULTI, NULL, dop, *lattice, cg_arg); 

  for(int i=0;i<nmass;i++) true_rsd[i] = targ_resid;

  bool fail(false);

#if 0
  //do a quick test to ensure the usual CPS CG inverter (which actually uses bfm here but has been tested and is correct) and the multimass get the same answer
  dop.InvCg((Vector*)out_1[0], (Vector*)in, 0.0, &true_rsd[0]); //This works
  //InvCg_CPS((Vector*)out_1[0], (Vector*)in, 0.0, &true_rsd[0], dop, *lattice, cg_arg); //this also works

  Float* o10f = (Float*)out_1[0];
  Float* o20f = (Float*)out_2[0];

  for(int i=0;i<f_size_cb;i++){
    if( fabs(o10f[i]-o20f[i])>1e-08 ){
      printf("InvCGtest fail %d: %f %f, ratio %f\n",i,o10f[i],o20f[i],o10f[i]/o20f[i]);
      fail=true;
    }
  }
  if(fail){
    printf("Failed InvCg test\n"); exit(-1);
  }else printf("Passed InvCg test\n");
  //test worked
#endif


  Float kappa = 1.0 / (2.0 * (cg_arg.mass + 4.0));

  //In InCg BFM solution vector has to be multiplied by 0.25/kappa^2 when converting to cps
  //   sol = (MdagM)^-1 src
  //Thus (MdagM)^-1 BFM is 4*kappa^2 larger than (MdagM)^-1 CPS
  //Thus MdagM BFM is 0.25/kappa^2 smaller than MdagM CPS
  //Here we invert  MdagM + shift
  //We therefore need to convert the shift to  shift*0.25/kappa^2  such that the normalization matches that of the MdagM operator

  //Residuals
  // true_res = |src - (MatPcDagMatPc + shift) * sol| / |src|
  //BFM calculation;  // true_res = |src - (MatPcDagMatPc + shift)_BFM * sol_BFM| / |src| = |src - (MatPcDagMatPc + shift)_CPS*0.25/kappa^2 * sol_BFM| / |src| = |src - (MatPcDagMatPc + shift)_CPS * sol_CPS| / |src|
  //Hence residuals should agree between BFM and CPS
  

  for(int i=0;i<nmass;i++){
    // mass[i] *= 0.25/kappa/kappa;  //0.25/kappa/kappa;
    true_rsd[i] = targ_resid;
  }






  dop.MInvCG((Vector**)out_1, (Vector*)in, src_norm_sq, &mass[0], nmass, 0, &true_rsd[0], MULTI, NULL); //bfm version
  
  fail=false;
  for(int j=0;j<nmass;j++){
    Float* o1f = (Float*)out_1[j];
    Float* o2f = (Float*)out_2[j];

    for(int i=0;i<f_size_cb;i++){
      if( fabs(o1f[i]-o2f[i])>1e-08 ){
	printf("MInvCGtest %d fail %d: %f %f, ratio %f\n",j,i,o1f[i],o2f[i],o1f[i]/o2f[i]);
	fail=true;
      }
    }

    if(!fail) printf("Passed MInvCg test for shift %d\n",j);
  }
  if(fail){
    printf("Failed MInvCg test\n"); exit(-1);
  }else printf("Passed MInvCg test\n");
    
  pfree(in);
  for(int j=0;j<3;j++){
    pfree(out_1[j]);
    pfree(out_2[j]);
  }
  pfree(out_1);
  pfree(out_2);
  delete lattice;

  return 0;
}



int main(int argc, char *argv[])
{
    const char *fname = "main()";

    setup(argc, argv);

    //VRB.ElapsedTime("CPS", fname);

    // OpenMP test
    VRB.Result( "CPS", "main", "omp_get_num_threads[1] -> %d", omp_get_num_threads() );
    
    #pragma omp parallel
    {
      if ( UniqueID() == 0 && omp_get_thread_num() == 0 ) {
        VRB.Result( "CPS", "main", "omp_get_num_threads[2] -> %d", omp_get_num_threads() );
      }
    }

    VRB.Result( "CPS", "main", "omp_get_num_threads[3] -> %d", omp_get_num_threads() );
    
    return MinvCGtest();

    //////////////////////////////////////////////////////////////////////
    // creating actions and integrators
    AlgMomentum mom;
    AlgActionGauge gauge(mom, gauge_arg);
    AlgActionRationalQuotient rat_quo(mom, rat_quo_arg);
    AlgActionRationalQuotient rat_quo_2(mom, rat_quo_arg_2);

    IntABArg sum_arg;
    sum_arg.A_steps = 1;
    sum_arg.B_steps = 1;
    sum_arg.level = EMBEDDED_INTEGRATOR;
    AlgIntSum sum(rat_quo, rat_quo_2, sum_arg);

    //!< Construct numerical integrators
    AlgIntAB &ab1 = AlgIntAB::Create(mom, gauge,       ab1_arg);
    AlgIntAB &ab2 = AlgIntAB::Create(ab1, sum,         ab2_arg);
    //////////////////////////////////////////////////////////////////////

    int traj = evo_arg.traj_start;
    for(int conf = 0; conf< evo_arg.gauge_configurations; ++conf) {
        CommonArg common_arg_plaq;
        CommonArg common_arg_pbp;
        CommonArg common_arg_hmc;

        truncate_it(&common_arg_plaq , evo_arg.plaquette_stem      , traj);
        truncate_it(&common_arg_pbp  , evo_arg.pbp_stem            , traj);
        truncate_it(&common_arg_hmc  , evo_arg.evo_stem            , traj);

        // Inner trajectory loop
        for(int i = 0; i < evo_arg.gauge_unload_period; ++i, ++traj) {

            measure_plaq(common_arg_plaq);
            measure_pbp(common_arg_pbp, traj);

            //VRB.ElapsedTime("CPS", "main[2]");
            VRB.Result( "CPS", "main", "omp_get_num_threads[4] -> %d", omp_get_num_threads() );
            run_hmc(common_arg_hmc, traj, ab2);
        }//End of inter-cfg sweep

        checkpoint(traj);

    } //End config loop

    AlgIntAB::Destroy(ab2);
    AlgIntAB::Destroy(ab1);

    End();

    VRB.Result(cname, fname, "Program ended normally.\n");
    
    return 0;
}

#undef encode_vml
#define encode_vml(arg_name, traj) do{                                  \
        char vml_file[256];                                             \
        sprintf(vml_file, #arg_name".%d", traj);                        \
        if( !arg_name.Encode(vml_file, #arg_name) ){                    \
            ERR.General(cname, fname, #arg_name " encoding failed.\n"); \
        }                                                               \
    }while(0)

void checkpoint(int traj)
{
    const char *fname="checkpoint()";

    char lat_file[256];
    char rng_file[256];

    Float time = -dclock();

    // Save this config to disk
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);

    sprintf(lat_file,"%s.%d",evo_arg.gauge_file_stem,traj);
    QioArg wt_arg(lat_file,0.001);

    wt_arg.ConcurIONumber=evo_arg.io_concurrency;
    WriteLatticeParallel wl;
    wl.setHeader(evo_arg.ensemble_id,evo_arg.ensemble_label,traj);
    wl.write(lat,wt_arg);

    if(!wl.good())
        ERR.General(cname,fname,"Failed write lattice %s",lat_file);

    LatticeFactory::Destroy();

    // Save the RNG's
    sprintf(rng_file,"%s.%d",evo_arg.rng_file_stem,traj);
    if ( !LRG.Write(rng_file) )
        ERR.General(cname,fname,"Failed write RNG file %s",rng_file);

    // Update the parameter files for restart
    do_arg.start_seed_filename = rng_file;
    do_arg.start_seed_kind = START_SEED_FILE;
    do_arg.start_conf_filename = lat_file;
    do_arg.start_conf_kind = START_CONF_FILE;
    evo_arg.traj_start     = traj;

    encode_vml(hmc_arg, traj);
    encode_vml(gauge_arg, traj);
    encode_vml(rat_quo_arg, traj);
    encode_vml(rat_quo_arg_2, traj);
    encode_vml(ab1_arg, traj);
    encode_vml(ab2_arg, traj);
    encode_vml(pbp_arg, traj);
    encode_vml(do_arg, traj);
    encode_vml(evo_arg, traj);

    time += dclock();
    print_flops("","checkpoint()",0,time);
}

void truncate_it(CommonArg *common_arg, const char stem[], int traj)
{
    char fnbuf[1024];
    sprintf(fnbuf, "%s.%d", stem, traj);
    FILE *truncate_it = Fopen(fnbuf, "w");
    Fclose(truncate_it);
    common_arg->set_filename(fnbuf);
}

void measure_plaq(CommonArg &common_arg)
{
    const char *fname = "measure_plaq()";

    Float dtime = -dclock();

    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON);
    AlgPlaq plaq(lat, &common_arg, &no_arg);
    plaq.run();
    LatticeFactory::Destroy();

    dtime += dclock();
    print_flops("AlgPlaq", "run()", 0, dtime);	
}

void measure_pbp(CommonArg &common_arg, int traj)
{
  return; 


    const char *fname = "measure_pbp()";

    // fix pbp_arg
    pbp_arg.src_u_s = 0;
    pbp_arg.src_l_s = 0;
    pbp_arg.snk_u_s = 0;
    pbp_arg.snk_l_s = 0;

    const int g_int = evo_arg.gauge_unload_period;
    if (traj % g_int == 0 && evo_arg.measure_pbp) {
        Float dtime = -dclock();

        LRGState rng_state;
        rng_state.GetStates();

        Lattice &lat = LatticeFactory::Create(F_CLASS_WILSON_TM, G_CLASS_NONE);
        VRB.Result( "cps", "measure_pbp", "LatticeFactory::Create(F_CLASS_BFM, G_CLASS_NONE)" );

        AlgPbp pbp(lat, &common_arg, &pbp_arg);

        for(int pbp_counter = 0; pbp_counter < evo_arg.measure_pbp; pbp_counter++) {
            pbp.run();
        }
        LatticeFactory::Destroy();
        rng_state.SetStates();

        dtime += dclock();
        print_flops("AlgPbp", "run()", 0, dtime);	
    }
}

void run_hmc(CommonArg &common_arg, int traj, AlgIntAB &int_ab)
{
    const char *fname = "run_hmc()";

    Float dtime = -dclock();

    if ( (evo_arg.reproduce_interval > 0) &&
         (traj % evo_arg.reproduce_interval) == 0 ) {
        VRB.Result(cname,fname,"Running traj %d with reproduction\n",traj);
        hmc_arg.reproduce = REPRODUCE_YES;
    } else {
        VRB.Result(cname,fname,"Running traj %d without reproduction\n",traj);
        hmc_arg.reproduce = REPRODUCE_NO;	
    }
    
    AlgHmc hmc(int_ab, common_arg, hmc_arg);
    hmc.run();

    dtime += dclock();
    print_flops("AlgHmc", "run()", 0, dtime);
}
