//#include <config.h>
//#include <util/time_cps.h>
#include <sys/time.h>
#include <qalloc.h>
#include "asq_data_types.h"
#include "asqtad_int.h"

//CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOp class CG solver methods.

*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-05-10 05:51:23 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_asqtad/qcdoc/asqtad_cg.C,v 1.10 2012-05-10 05:51:23 chulwoo Exp $
//  $Id: asqtad_cg.C,v 1.10 2012-05-10 05:51:23 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: asqtad_cg.C,v $
//  $Revision: 1.10 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_asqtad/qcdoc/asqtad_cg.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// inv_cg.C
//
//------------------------------------------------------------------
#include <math.h>
#include <stdio.h>

#undef PROFILE


#ifdef PROFILE
#include <time.h>
#include <sys/time_cps.h>
void report_flops(int flops, struct timeval *start,struct timeval *end);
#endif

#ifdef  PARALLEL
//Uncomment the following line to activate reproducibility test
#endif


#ifdef ASQD_SINGLE
extern "C"   void asq_vaxmy_cpp(Float *scale,PTvector *mult,PTvector *sub,int ncvec);
#define asq_vaxmy(A,B,C,D) asq_vaxmy_cpp(A,B,C,D)
extern "C"   void asq_vaxmy_vxdot_s(Float *scale, PTvector *mult, PTvector *sub, 
               int ncvec, Float *norm);
#define asq_vaxmy_vxdot(A,B,C,D,E) asq_vaxmy_vxdot_s(A,B,C,D,E)
extern "C"   void asq_vaxpy3_s(Float *res,Float *scale,Float *mult,Float *add, int ncvec);
#define asq_vaxpy3(A,B,C,D,E) asq_vaxpy3_s(A,B,C,D,E)
extern "C"   void asq_vaxpy3_norm_s(Float *res,Float *scale,Float *mult,
               Float *add, int ncvec,Float *norm);
#define asq_vaxpy3_norm(A,B,C,D,E,F) asq_vaxpy3_norm_s(A,B,C,D,E,F)
#else
extern "C"   void asq_vaxmy(Float *scale,PTvector *mult,PTvector *sub,int ncvec);
extern "C"   void asq_vaxmy_vxdot(Float *scale, PTvector *mult, PTvector *sub, 
               int ncvec, Float *norm);
extern "C"   void asq_vaxpy3(Float *res,Float *scale,Float *mult,Float *add, 
              int ncvec);
extern "C"   void asq_vaxpy3_norm(Float *res,Float *scale,Float *mult,
               Float *add, int ncvec,Float *norm);
#endif

Float asq_print_flops(char *cname, char *fname, unsigned long long nflops, struct timeval *start, struct timeval
*end){
        int sec = end->tv_sec - start->tv_sec;
        int usec = end->tv_usec - start->tv_usec;
        Float time = sec + 1.e-6*usec;
        printf("%s::%s: %e flops /%e seconds = %e MFlops\n",
        cname,fname,(Float)nflops,time,(Float)nflops/(time*1.e6));
        return nflops/time;
}


// The granularity used in the interleaving
#define GRAN 12

#undef PROFILE
void AsqD::Dslash(Float *out, Float *in ){
 unsigned long long nflops = 1146*vol;
 static int called=0;
  printf("f_size_cb=%d\n",f_size_cb);
#ifdef PROFILE
  struct timeval start,end;
  gettimeofday(&start,NULL);
#endif
  dirac(out+f_size_cb, in, node_odd, 0);
  dirac(out, in+f_size_cb, 1-node_odd, 0);
#ifdef PROFILE
  gettimeofday(&end,NULL);
  printf("DiracOpAsqtad::MatPcDagMatPc:: ");
  print_flops(nflops,&start,&end);
#endif
}

#undef PROFILE
void AsqD::MdagM(Float *mass_sq, Float *out, Float *in, int m_odd,
 Float *dot_prd){
 unsigned long long nflops = 1146*vol;
 static int called=0;
 if (!called){
// printf("mass_sq=%e\n",mass_sq);
  called=1;
 }
 int odd = (node_odd+m_odd)%2;
#ifdef PROFILE
  printf("in=%p frm_tmp=%p out=%p\n",in,frm_tmp,out);
  struct timeval start,end;
  gettimeofday(&start,NULL);
#endif
  dirac((Float *)frm_tmp, (Float *)in, odd, 0);
#if 0
    Float *tmp_p = frm_tmp;
    for(int j = 0;j<vol*3;j++)
      if(fabs(*(tmp_p+j))>1e-5) fprintf(stderr,"%d %d %d %d frm_tmp[%d]=%e\n",
           CoorT(),CoorX(),CoorY(),CoorZ(),j,*(tmp_p+j));
#endif
  dirac((Float *)out, (Float *)frm_tmp, 1-odd, 0);

#ifdef PROFILE
  gettimeofday(&end,NULL);
  asq_print_flops(cname,"MdagM",nflops,&start,&end);
  gettimeofday(&start,NULL);
#endif


  if( dot_prd !=0 ){
	asq_vaxmy_vxdot(mass_sq,in,out,f_size_cb/6,dot_prd);
//    CGflOps +=f_size_cb*4;
    nflops =f_size_cb*4;
  } else {
//    CGflops +=f_size_cb*2;
    nflops =f_size_cb*2;
    asq_vaxmy(mass_sq,in,out,f_size_cb/6);
  }
#ifdef PROFILE
  gettimeofday(&end,NULL);
  asq_print_flops(cname,"MdagM",nflops,&start,&end);
#endif

}

//------------------------------------------------------------------
/*!
  Solves \f$ M^\dagger M out = in \f$ for \a out using the Conjugate
  Gradient method, where \a M is the
  fermion matrix, possibly odd-even preconditioned, possibly a single parity
  of the odd-even preconditioned fermion matrix.
  The residual used for the stopping
  criterion is \f$ |f_{in} - M^\dagger M f_{out}| / |f_{in}| \f$.

  \param out The initial guess of solution vector.
  \param in The source vector
  \param src_norm_sq The square norm of the source vector. If this is set to
  zero it will be calculated inside this method.
  \param true_res Whether or not to report the true residual. This will
  point to the true residual  if it initially points to something non-zero.
  \return The number of solver iterations.
  \post \a f_out contains the solution vector.
  \post true_res The true residual, if this was non-zero to start with.
*/
//------------------------------------------------------------------
//static gauge_agg *uc_l_save[2];
//static gauge_agg *uc_nl_save[2];
//static Float *dir_save;
int AsqD::cg_called = 0;
#define PROFILE
int AsqD::InvCg( InvArg *inv_arg, Float *out, 
		   Float *in, 
//		   Float src_norm_sq, 
		   Float *true_res,
		   int odd ){
  int itr;                       // Current number of CG iterations
  int max_itr;                       // Max number of CG iterations
  Float stp_cnd;                   // Stop if residual^2 <= stp_cnd
  Float res_norm_sq_prv;          // The previous step |residual|^2
  Float res_norm_sq_cur;           // The current step |residual|^2 
  Float a;
  Float b;
  Float d;
  int i;
  char *fname = "InvCg(V*,V*,F,F*)";

//  Float *tmp_f;

#if 0
// Print out input parameters
//------------------------------------------------------------------
  VRB.Input(cname,fname,
	    "stop_rsd = %e\n",Float(dirac_arg->stop_rsd));
  VRB.Input(cname,fname,
	    "max_num_iter = %d\n",dirac_arg->max_num_iter);
  VRB.Input(cname,fname,
	    "mass = %e\n",Float(dirac_arg->mass));
  VRB.Input(cname,fname,
	    "src_norm_sq = %e\n",Float(src_norm_sq));
#endif


//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------

// Set the source vector pointer
//------------------------------------------------------------------
  PTvector *src = in;
//  Float *src_tmp = (Float *)src;

// Set the solution vector pointer
//------------------------------------------------------------------
  PTvector *sol = out;

// Allocate memory for the solution/residual field.
//------------------------------------------------------------------
  Float *res = (Float *) qalloc(QCOMMS|QFAST,f_size_cb * sizeof(Float));
  if(res == 0){
    res = (Float *) qalloc(QCOMMS,f_size_cb * sizeof(Float));
//    printf("res=%p\n",res);
  }
  if(res == 0) PointerErr(cname,fname, "res");

// Allocate memory for the direction vector dir.
//------------------------------------------------------------------
  PTvector *dir;
  if(vol >1024) dir=0;
    else dir = (PTvector *) qalloc(QCOMMS|QFAST,f_size_cb * sizeof(Float));
  if(dir == 0){
    dir = (PTvector *) qalloc(QCOMMS,f_size_cb * sizeof(Float));
  }
  if(dir == 0)
    PointerErr(cname,fname, "dir");

// Allocate mem. for the result vector of matrix multiplication mmp.
//------------------------------------------------------------------
  PTvector *mmp;
  if(vol >1024) mmp=0;
    else mmp = (PTvector *) qalloc(QCOMMS|QFAST,f_size_cb * sizeof(Float));
  if(mmp == 0){
    mmp = (PTvector *) qalloc(QCOMMS,f_size_cb * sizeof(Float));
  }
 // printf("dir=%p mmp=%p\n",dir,mmp);
  if(mmp == 0)
    PointerErr(cname,fname, "mmp");

// If src_norm_sq is not provided calculate it
//------------------------------------------------------------------
  Float  src_norm_sq = NormSqNode(src,f_size_cb);
    Sum(&src_norm_sq);

// Calculate stopping condition
//------------------------------------------------------------------
  stp_cnd = src_norm_sq * inv_arg->stop_rsd * inv_arg->stop_rsd;
  //printf("AsqD::InvCg(): stp_cnd =%e\n", Float(stp_cnd));
  Float mass_sq = 4*inv_arg->mass*inv_arg->mass;
  //  printf("mass_sq=%e\n",mass_sq);
  
#ifdef PROFILE
  struct timeval start;
  struct timeval end;
    unsigned long long CGflops    = 0;
    unsigned long long nflops = 0;
    unsigned long long nflops_tmp;
    gettimeofday(&start,NULL);
#endif

  inv_arg->final_iter = 0;

  for(int restart=0;restart<=inv_arg->restart;restart++){

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
    // Mmp = MatPcDagMatPc * sol
    MdagM(&mass_sq,mmp, sol,odd);
  
    // res = src
    CopyVec(res,src, f_size_cb);
  
    // res -= mmp
    VecMinusEquVec(res,mmp, f_size_cb);
  
    // dir = res
    CopyVec(dir,res, f_size_cb);  
  
    // res_norm_sq_cur = res * res
    res_norm_sq_cur = NormSqNode(dir,f_size_cb);
  
    Sum(&res_norm_sq_cur);
  
    // if( |res|^2 <= stp_cnd ) we are done
//    printf("AsqD::InvCg: |res[0]|^2 = %e\n", Float(res_norm_sq_cur));
    itr = 0;
    max_itr = inv_arg->niter-1;
    if(res_norm_sq_cur <= stp_cnd) max_itr = 0;

//------------------------------------------------------------------
// Loop over CG iterations
//------------------------------------------------------------------
//  Gint::SynchMachine();
    for(i=0; i < max_itr; i++){
      itr++;
      res_norm_sq_prv = res_norm_sq_cur;
  
  #if 0
      Float *tmp_p = dir;
      for(int j = 0;j<vol*3;j++)
        if(fabs(*(tmp_p+j))>1e-5) fprintf(stderr,"%d %d %d %d dir[%d]=%e\n",
             CoorT(),CoorX(),CoorY(),CoorZ(),j,*(tmp_p+j));
  #endif
  
      // mmp = MatPcDagMatPc * dir
      // d = <dir, MatPcDagMatPc*dir>
      MdagM(&mass_sq,mmp, dir, odd,&d);
  
      CGflops += 1146*vol+f_size_cb*4;
  
      Sum(&d);
      CGflops +=f_size_cb*2;
  
      // If d = 0 we are done
      if(d == 0.0) break;
      //??? or should we give a warning or error? Yes we should, really.
  
      a = res_norm_sq_prv / d;
      // sol = a * dir + sol;
   
      asq_vaxpy3(sol,&a,dir,sol,f_size_cb/6);
      CGflops +=f_size_cb*2;
   
   
      // res = - a * (MatPcDagMatPc * dir) + res;
      // res_norm_sq_cur = res * res
  
     // a *= -1.0;
     Float ma = -a;
      asq_vaxpy3_norm(res,&ma,mmp,res,f_size_cb/6,&res_norm_sq_cur);
      Sum(&res_norm_sq_cur);
      CGflops +=f_size_cb*6;
  
      b = res_norm_sq_cur / res_norm_sq_prv;
  
      // dir = b * dir + res;
      asq_vaxpy3(dir,&b,dir,res,f_size_cb/6);
      CGflops +=f_size_cb*2;
  
  
      // if( |res|^2 <= stp_cnd ) we are done
  //    printf("%s::%s:|res[%d]|^2 = %e\n", cname,fname,itr, Float(res_norm_sq_cur));
      if(res_norm_sq_cur <= stp_cnd) break;
  
    }

    if(res_norm_sq_cur > stp_cnd) 
      printf("CG reached max iterations = %d, restart=%d, |res|^2 = %e\n", itr+1, restart, Float(res_norm_sq_cur) );
  
    inv_arg->final_iter += itr;
  }
  // It has not reached stp_cnd: Issue a warning
  if(res_norm_sq_cur > stp_cnd) 
    printf("CG reached max iterations = %d, restart=%d, |res|^2 = %e\n", inv_arg->final_iter , inv_arg->restart, Float(res_norm_sq_cur) );

//------------------------------------------------------------------
// Done. Finish up and return
//------------------------------------------------------------------
  // Calculate and set true residual: 
  // true_res = |src - MatPcDagMatPc * sol| / |src|

  MdagM(&mass_sq, mmp, sol,odd);
  CGflops += 1146*vol+f_size_cb*2;
  CopyVec(dir,src, f_size_cb);
  VecMinusEquVec(dir,mmp, f_size_cb);
  res_norm_sq_cur = NormSqNode(dir,f_size_cb);
  Sum(&res_norm_sq_cur);
  Float tmp = res_norm_sq_cur / src_norm_sq;
  tmp = sqrt(tmp);
  if(true_res != 0){
    *true_res = tmp;
  }
//  printf("%s::%s: True |res| / |src| = %e, iter = %d\n", cname,fname, Float(tmp), itr+1);
#ifdef PROFILE
  gettimeofday(&end,NULL);
//  asq_print_flops(cname,fname,CGflops,&start,&end); 
  inv_arg->final_flop = (double)CGflops;
  int sec = end.tv_sec - start.tv_sec;
  int usec = end.tv_usec - start.tv_usec;
  inv_arg->final_sec = sec + 1.e-6*usec;
#endif
  inv_arg->final_rsq = res_norm_sq_cur;
//  inv_arg->final_iter = itr;

  // Free memory
  qfree(mmp);
  qfree(dir);
  qfree(res);

  // Return number of iterations
//  printf("AsqD::InvCG\n");
  cg_called++;
  return inv_arg->final_iter;

}



//CPS_END_NAMESPACE
