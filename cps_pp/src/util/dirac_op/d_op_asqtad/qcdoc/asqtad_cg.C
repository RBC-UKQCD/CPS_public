//#include <config.h>
//#include <util/time.h>
#include <qalloc.h>
#include <util/asqtad_int.h>

//CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOp class CG solver methods.

*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2005-03-07 00:22:22 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_asqtad/qcdoc/asqtad_cg.C,v 1.3 2005-03-07 00:22:22 chulwoo Exp $
//  $Id: asqtad_cg.C,v 1.3 2005-03-07 00:22:22 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: asqtad_cg.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_asqtad/qcdoc/asqtad_cg.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// inv_cg.C
//
//------------------------------------------------------------------
#if 0
CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <comms/nga_reg.h>
#include <comms/cbuf.h>
#include <qcdocos/gint.h>
CPS_START_NAMESPACE
#endif
#include <math.h>
#include <stdio.h>

#undef PROFILE


#ifdef PROFILE
#include <time.h>
#include <sys/time.h>
void report_flops(int flops, struct timeval *start,struct timeval *end);
#endif

#ifdef  PARALLEL
//Uncomment the following line to activate reproducibility test
#undef REPRODUCE_TEST
//#undef REPRODUCE_TEST
#endif

#ifdef  REPRODUCE_TEST
CPS_END_NAMESPACE
#include <comms/sysfunc.h>
CPS_START_NAMESPACE
#endif

extern "C" { 
  void asqd_vaxmy(Float *scale,vector *mult,vector *sub,int ncvec);
  void asqd_vaxmy_vxdot(Float *scale, vector *mult, vector *sub, int
ncvec, Float *norm);
  void invcg_r_norm(Float *resa, Float *scale, Float *mult, Float *add, 
		      int ncvec, Float *norm);
  void invcg_xp_update(Float *out1, Float *out2, Float *A, Float *B, 
		       Float *mult, Float *add, int size);
  void asqd_vaxpy3(Float *res,Float *scale,Float *mult,Float *add, int
ncvec);
  void asqd_vaxpy3_norm(Float *res,Float *scale,Float *mult,Float *add,
int ncvec,Float *norm);

}

// The granularity used in the interleaving
#define GRAN 12

void AsqD::MdagM(Float *mass_sq, Float *out, Float *in, 
 Float *dot_prd){
 long nflops = 1146*vol;
 static int called=0;
 if (!called){
 printf("mass_sq=%e\n",mass_sq);
  called=1;
 }
#undef PROFILE
#ifdef PROFILE
  struct timeval start,end;
  gettimeofday(&start,NULL);
#endif
  dirac((Float *)frm_tmp, (Float *)in, 0, 0);
  dirac((Float *)out, (Float *)frm_tmp, 1, 0);
//  CGflops += nflops;


  if( dot_prd !=0 ){
	asqd_vaxmy_vxdot(mass_sq,in,out,f_size_cb/6,dot_prd);
//    CGflOps +=f_size_cb*4;
    nflops +=f_size_cb*4;
  } else {
//    CGflops +=f_size_cb*2;
    nflops +=f_size_cb*2;
    asqd_vaxmy(mass_sq,in,out,f_size_cb/6);
  }

#ifdef PROFILE
  gettimeofday(&end,NULL);
  printf("DiracOpAsqtad::MatPcDagMatPc:: ");
  print_flops(nflops,&start,&end);
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
int AsqD::InvCg( InvArg *inv_arg, Float *out, 
		   Float *in, 
//		   Float src_norm_sq, 
		   Float *true_res){
  int f_size_cb;     // Node checkerboard size of the fermion field
  int itr;                       // Current number of CG iterations
  int max_itr;                       // Max number of CG iterations
  Float stp_cnd;                   // Stop if residual^2 <= stp_cnd
  Float res_norm_sq_prv;          // The previous step |residual|^2
  Float res_norm_sq_cur;           // The current step |residual|^2 
  Float a;
  Float b;
  Float d;
  int i, j;
  char *fname = "InvCg(V*,V*,F,F*)";

  printf("AsqD::InvCG\n");
  Float *tmp_f;
  tmp_f = (Float *)in;
  printf("in[0]=%0.14e\n",*tmp_f);
  tmp_f = (Float *)fat;
  printf("Fat[0]=%0.14e\n",*tmp_f);
  tmp_f = (Float *)naik;
  printf("Naik[0]=%0.14e\n",*tmp_f);
// Flash the LED and then turn it off
// Flash the LED and then turn it off
// Flash the LED and then turn it off
//------------------------------------------------------------------
//  VRB.LedFlash(cname,fname,3);
//  VRB.LedOff(cname,fname);
//  VRB.Func(cname,fname);


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
  vector *src = in;
//  Float *src_tmp = (Float *)src;

// Set the solution vector pointer
//------------------------------------------------------------------
  vector *sol = out;

// Set the node checkerboard size of the fermion field
//------------------------------------------------------------------

    f_size_cb = vol * 6/2;
  
  if (f_size_cb % GRAN != 0) 
    fprintf(stderr,"%s::%s: Field length %d is not a multiple of granularity %d\n", cname,fname,GRAN, f_size_cb);

// Allocate memory for the solution/residual field.
//------------------------------------------------------------------
  Float *X = (Float *) qalloc(QCOMMS|QFAST,2*f_size_cb * sizeof(Float));
  if(X == 0){
    X = (Float *) qalloc(QCOMMS,2*f_size_cb * sizeof(Float));
    printf("X=%p\n",X);
  }
  if(X == 0) PointerErr(cname,fname, "X");
//  VRB.Smalloc(cname,fname, "X", X, 2*f_size_cb * sizeof(Float));
//

// Allocate memory for the direction vector dir.
//------------------------------------------------------------------
  vector *dir;
  if(vol >1024) dir=0;
    else dir = (vector *) qalloc(QCOMMS|QFAST,f_size_cb * sizeof(Float));
  if(dir == 0){
    dir = (vector *) qalloc(QCOMMS,f_size_cb * sizeof(Float));
  }
  if(dir == 0)
    PointerErr(cname,fname, "dir");
//  VRB.Smalloc(cname,fname, "dir", dir, f_size_cb * sizeof(Float));

// Allocate mem. for the result vector of matrix multiplication mmp.
//------------------------------------------------------------------
  vector *mmp;
  if(vol >1024) mmp=0;
    else mmp = (vector *) qalloc(QCOMMS|QFAST,f_size_cb * sizeof(Float));
  if(mmp == 0){
    mmp = (vector *) qalloc(QCOMMS,f_size_cb * sizeof(Float));
  }
  printf("dir=%p mmp=%p\n",dir,mmp);
  if(mmp == 0)
    PointerErr(cname,fname, "mmp");
//  VRB.Smalloc(cname,fname, "mmp", mmp, f_size_cb * sizeof(Float));

// If src_norm_sq is not provided calculate it
//------------------------------------------------------------------
//  if(src_norm_sq == 0){
  Float  src_norm_sq = NormSqNode(src,f_size_cb);
    Sum(&src_norm_sq);
//  }

// Calculate stopping condition
//------------------------------------------------------------------
  stp_cnd = src_norm_sq * inv_arg->stop_rsd * inv_arg->stop_rsd;
//  VRB.Flow(cname,fname, "stp_cnd =%e\n", Float(stp_cnd));

#ifdef REPRODUCE_TEST 
  
// Allocate space for storing solution
//------------------------------------------------------------------
  vector *sol_store = (vector *) smalloc(f_size_cb * sizeof(Float));
  if(sol_store == 0) PointerErr(cname,fname, "sol_store");
//  VRB.Smalloc(cname,fname, "sol_store", sol_store, f_size_cb * sizeof(Float));

// Allocate space for storing d
//------------------------------------------------------------------
  Float *d_store = (Float *) smalloc( dirac_arg->max_num_iter-1 * sizeof(Float));

  if(d_store == 0) PointerErr(cname,fname, "d_store");
//  VRB.Smalloc(cname,fname, "d_store", d_store, dirac_arg->max_num_iter-1 * sizeof(Float));

  for ( int n = 0; n < dirac_arg->max_num_iter-1; n++ )  d_store[n] = 0;

  CopyVec(sol_store,sol, f_size_cb);

  for ( int test = 0; test < 2; test++ ) {
    if (test == 1) sol-> CopyVec(sol_store, f_size_cb);
    
#endif

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
  Float mass_sq = 4*inv_arg->mass*inv_arg->mass;
  printf("mass_sq=%e\n",mass_sq);
  MdagM(&mass_sq,mmp, sol);

  // res = src
  CopyVec(dir,src, f_size_cb);

  // res -= mmp
  VecMinusEquVec(dir,mmp, f_size_cb);

  // dir = res
  //dir->CopyVec(res, f_size_cb);  

  Float *Fsol = (Float*)sol;
  Float *Fdir = (Float*)dir;
  Float *Fmmp = (Float*)mmp;
  Float *Xptr;

  // Interleave solution and residual
  Xptr = X;
  for (j=0; j<f_size_cb/GRAN;j++) {
    for (i=0; i<GRAN; i++) *Xptr++ = *(Fsol+j*GRAN+i);
    for (i=0; i<GRAN; i++) *Xptr++ = *(Fdir+j*GRAN+i);
  }

  // res_norm_sq_cur = res * res
  res_norm_sq_cur = NormSqNode(dir,f_size_cb);

  Sum(&res_norm_sq_cur);

  // if( |res|^2 <= stp_cnd ) we are done
//  VRB.Flow(cname,fname, "|res[0]|^2 = %e\n", Float(res_norm_sq_cur));
  itr = 0;
  max_itr = inv_arg->niter-1;
  if(res_norm_sq_cur <= stp_cnd) max_itr = 0;


#ifdef PROFILE
  struct timeval start;
  struct timeval end;
  struct timeval linalg_tmp;
  struct timeval linalg_start;
  struct timeval linalg_end;

  CGflops    = 0;
  int nflops = 0;
  int nflops_tmp;
  gettimeofday(&start,NULL);

#endif

//------------------------------------------------------------------
// Loop over CG iterations
//------------------------------------------------------------------
//  Gint::SynchMachine();

  for(i=0; i < max_itr; i++){
//    timeval start,end;
    itr++;
    res_norm_sq_prv = res_norm_sq_cur;

    // mmp = MatPcDagMatPc * dir
    // d = <dir, MatPcDagMatPc*dir>
    MdagM(&mass_sq,mmp, dir, &d);

#ifdef REPRODUCE_TEST 

    /* Check reproducibility */
    if ( test == 0) d_store[ i ] = d;
    else if ( d != d_store[ i ] )
      InterruptExit(-1, "NODE FAILS TO REPRODUCE");
    /* End of Check */

#endif
#ifdef PROFILE
//    gettimeofday(&linalg_tmp,NULL);
//    nflops_tmp = 0;
#endif
  
    Sum(&d);

    // If d = 0 we are done
    if(d == 0.0) break;
    //??? or should we give a warning or error? Yes we should, really.

    a = -res_norm_sq_prv / d;

    // res = - a * (MatPcDagMatPc * dir) + res;
    // res_norm_sq_cur = res * res

    invcg_r_norm(X+GRAN, &a, Fmmp, X+GRAN, f_size_cb/GRAN, &res_norm_sq_cur);
    Sum(&res_norm_sq_cur);
#ifdef PROFILE
//    nflops_tmp +=f_size_cb*4;
#endif
//    CGflops+=f_size_cb*4;

    a = -a;
    b = res_norm_sq_cur / res_norm_sq_prv;

    // sol = a * dir + sol;
    //sol->FTimesV1PlusV2(a, dir, sol, f_size_cb);
    // dir = b * dir + res;
    invcg_xp_update(X, Fdir, &a, &b, Fdir, X, f_size_cb/GRAN);


#ifdef PROFILE
//    linalg_start = linalg_tmp;
//    gettimeofday(&linalg_end,NULL);
//    nflops =nflops_tmp+f_size_cb*4;
#endif
//    CGflops+=f_size_cb*4;

    // if( |res|^2 <= stp_cnd ) we are done
//    VRB.Flow(cname,fname, "|res[%d]|^2 = %e\n", itr, Float(res_norm_sq_cur));
//    printf("%s::%s:|res[%d]|^2 = %e\n", cname,fname,itr, Float(res_norm_sq_cur));
    if(res_norm_sq_cur <= stp_cnd) break;

  }

#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(cname,fname,CGflops,&start,&end); 
#endif

  // It has not reached stp_cnd: Issue a warning
  if(itr == inv_arg->niter - 1){
 //   VRB.Warn(cname,fname, "CG reached max iterations = %d. |res|^2 = %e\n", itr+1, Float(res_norm_sq_cur) );
  }

//------------------------------------------------------------------
// Done. Finish up and return
//------------------------------------------------------------------
  // Calculate and set true residual: 
  // true_res = |src - MatPcDagMatPc * sol| / |src|
  Xptr = X-GRAN;
  for (j=0; j<f_size_cb; j++) {
    if (j%GRAN==0) Xptr += GRAN;
    *(Fsol++) = *(Xptr++);
  }

  MdagM(&mass_sq, mmp, sol);
  CopyVec(dir,src, f_size_cb);
  VecMinusEquVec(dir,mmp, f_size_cb);
  res_norm_sq_cur = NormSqNode(dir,f_size_cb);
  Sum(&res_norm_sq_cur);
  Float tmp = res_norm_sq_cur / src_norm_sq;
  tmp = sqrt(tmp);
  if(true_res != 0){
    *true_res = tmp;
  }
  tmp_f = (Float *)sol;
  printf("out[0]=%0.14e\n",*tmp_f);
  tmp_f = (Float *)mmp;
  printf("mmp[0]=%0.14e\n",*tmp_f);
  printf("%s::%s: True |res| / |src| = %e, iter = %d\n", cname,fname, Float(tmp), itr+1);

#ifdef REPRODUCE_TEST 
  }
  VRB.Sfree(cname, fname,"d_store", d_store);
  sfree(d_store);
  VRB.Sfree(cname, fname,"sol_store", sol_store);
  sfree(sol_store);

#endif

  // Free memory
//  VRB.Sfree(cname,fname, "mmp", mmp);
  qfree(mmp);
//  VRB.Sfree(cname,fname, "dir", dir);
  qfree(dir);
//  VRB.Debug("b ============\n");
//  VRB.Sfree(cname,fname, "X", X);
  Free(X);

//  VRB.Debug("a ============\n");

// Flash the LED and then turn it on
//------------------------------------------------------------------
//  VRB.FuncEnd(cname,fname);
//  VRB.LedFlash(cname,fname,2);
//  VRB.LedOn(cname,fname);

  // Return number of iterations
  printf("AsqD::InvCG\n");
  return itr+1;

}

#if 0

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
  \return The number of solver iterations.
  \post \a f_out contains the solution vector.
*/
//------------------------------------------------------------------
int DiracOp::InvCg(vector *out, vector *in, Float src_norm_sq)
{ return InvCg(out, in, src_norm_sq, 0); }


//------------------------------------------------------------------
/*!
  Solves \f$ M^\dagger M out = in \f$ for \a out, where \a M is the
  (possibly odd-even preconditioned) fermionic matrix, using the Conjugate
  Gradient method,
  The residual used for the stopping criterion  is
  \f$ |f_{in} - M^\dagger M f_{out}| / |f_{in}| \f$.

  \param out The initial guess of solution vector.
  \param in The source vector
  \param true_res Whether or not to report the true residual. This will
  point to the true residual  if it initially points to something non-zero.
  \return The number of solver iterations.
  \post \a f_out contains the solution vector.
  \post true_res The true residual, if this was non-zero to start with.
*/
// Same as original but with src_norm_sq=0.0
//------------------------------------------------------------------
int DiracOp::InvCg(vector *out, vector *in, Float *true_res)
{ return InvCg(out, in, 0.0, true_res); }


//------------------------------------------------------------------
// Same as original but with src_norm_sq=0.0, true_res=0
/*!
  Solves \f$ M^\dagger M out = in \f$ for \a out using the Conjugate
  Gradient method, where \a M is the
  fermion matrix, possibly odd-even preconditioned, possibly a single parity
  of the odd-even preconditioned fermion matrix.
  The residual used for the stopping
  criterion is \f$ |f_{in} - M^\dagger M f_{out}| / |f_{in}| \f$.
  
  \param out The initial guess of solution vector.
  \param in The source vector
  \return The number of solver iterations.
  \post \a f_out contains the solution vector.
*/
//------------------------------------------------------------------
int DiracOp::InvCg(vector *out, vector *in)
{ return InvCg(out, in, 0.0, 0); }


//------------------------------------------------------------------
// Same as original but with 
// in=f_in, out=f_out
/*!
  Solves \f$ M^\dagger M out = in \f$ for \a out using the Conjugate
  Gradient method, where \a M is the
  fermion matrix, possibly odd-even preconditioned, possibly a single parity
  of the odd-even preconditioned fermion matrix.
  The residual used for the stopping
  criterion is \f$ |f_{in} - M^\dagger M f_{out}| / |f_{in}| \f$.

  The initial guess of solution vector is the vector \a f_field_out passed
  as a constructor argument.
  The source vector is the vector \a f_field_in passed
  as a constructor argument.

  \param src_norm_sq The square norm of the source vector. If this is set to
  zero it will be calculated inside this method.
  \param true_res Whether or not to report the true residual. This will
  point to the true residual  if it initially points to something non-zero.
  \return The number of solver iterations.
  \post \a f_field_out contains the solution vector.
  \post true_res The true residual, if this was non-zero to start with.
*/
//------------------------------------------------------------------
int DiracOp::InvCg(Float src_norm_sq, Float *true_res)
{ return InvCg(f_out, f_in, src_norm_sq, true_res); }


//------------------------------------------------------------------
// Same as original but with 
// in=f_in, out=f_out, true_res=0
/*!
  Solves \f$ M^\dagger M out = in \f$ for \a out using the Conjugate
  Gradient method, where \a M is the
  fermion matrix, possibly odd-even preconditioned, possibly a single parity
  of the odd-even preconditioned fermion matrix.
  The residual used for the stopping
  criterion is \f$ |f_{in} - M^\dagger M f_{out}| / |f_{in}| \f$.

  The initial guess of solution vector is the vector \a f_field_out passed
  as a constructor argument.
  The source vector is the vector \a f_field_in passed
  as a constructor argument.

  \param src_norm_sq The square norm of the source vector. If this is set to
  zero it will be calculated inside this method.
  \return The number of solver iterations.
  \post \a f_field_out contains the solution vector.
*/
//------------------------------------------------------------------
int DiracOp::InvCg(Float src_norm_sq)
{ return InvCg(f_out, f_in, src_norm_sq, 0); }


//------------------------------------------------------------------
// Same as original but with 
// in=f_in, out=f_out, src_norm_sq=0.0
/*!
  Solves \f$ M^\dagger M out = in \f$ for \a out using the Conjugate
  Gradient method, where \a M is the
  fermion matrix, possibly odd-even preconditioned, possibly a single parity
  of the odd-even preconditioned fermion matrix.
  The residual used for the stopping
  criterion is \f$ |f_{in} - M^\dagger M f_{out}| / |f_{in}| \f$.

  The initial guess of solution vector is the vector \a f_field_out passed
  as a constructor argument.
  The source vector is the vector \a f_field_in passed
  as a constructor argument.

  \param true_res Whether or not to report the true residual. This will
  point to the true residual  if it initially points to something non-zero.
  \return The number of solver iterations.
  \post \a f_field_out contains the solution vector.
  \post true_res The true residual, if this was non-zero to start with.
*/
//------------------------------------------------------------------
int DiracOp::InvCg(Float *true_res)
{ return InvCg(f_out, f_in, 0.0, true_res); }


//------------------------------------------------------------------
// Same as original but with 
// in=f_in, out=f_out, src_norm_sq=0.0, true_res=0
/*!
  Solves \f$ M^\dagger M out = in \f$ for \a out using the Conjugate
  Gradient method, where \a M is the
  fermion matrix, possibly odd-even preconditioned, possibly a single parity
  of the odd-even preconditioned fermion matrix.
  The residual used for the stopping
  criterion is \f$ |f_{in} - M^\dagger M f_{out}| / |f_{in}| \f$.

  The initial guess of solution vector is the vector \a f_field_out passed
  as a constructor argument.
  The source vector is the vector \a f_field_in passed
  as a constructor argument.

  \return The number of solver iterations.
  \post \a f_field_out contains the solution vector.
*/
//------------------------------------------------------------------
int DiracOp::InvCg(void)
{ return InvCg(f_out, f_in, 0.0, 0); }
#endif

#ifdef PROFILE
#include <stdio.h>
void report_flops(int flops, struct timeval *start,struct timeval *end)
{

  double t;
  double mflops;

  t = ( end->tv_usec - start->tv_usec )*1.E-6;
  t+= ( end->tv_sec - start->tv_sec );

  mflops = (flops * 1.E-6) / t;
  printf("\t%ld:%ld -> %ld:%ld\n",
	 start->tv_sec,start->tv_usec,
	 end->tv_sec,end->tv_usec
	 );
  printf("\t%d flops %le seconds %lf Mflop/s\n",flops,t,mflops);
}
#endif

//CPS_END_NAMESPACE
