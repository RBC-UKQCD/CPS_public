#include <config.h>
#include <qalloc.h>
#include <util/time.h>

CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOp class CG solver methods.

*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-07-09 04:15:17 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_base/qcdoc/inv_cg.C,v 1.7 2004-07-09 04:15:17 chulwoo Exp $
//  $Id: inv_cg.C,v 1.7 2004-07-09 04:15:17 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: inv_cg.C,v $
//  $Revision: 1.7 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_base/qcdoc/inv_cg.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// inv_cg.C
//
//------------------------------------------------------------------
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
#include <math.h>
#include <stdio.h>
#include <qcdocos/gint.h>
CPS_START_NAMESPACE

#undef PROFILE

// PAB generic Dop flops reporting
//int DiracOp::CGflops;

#ifdef PROFILE
#include <time.h>
#include <sys/time.h>
void report_flops(int flops, struct timeval *start,struct timeval *end);
#endif

#ifdef  PARALLEL
//Uncomment the following line to activate reproducibility test
//#define REPRODUCE_TEST
#undef REPRODUCE_TEST
#endif

#ifdef  REPRODUCE_TEST
CPS_END_NAMESPACE
#include <comms/sysfunc.h>
CPS_START_NAMESPACE
#endif

#if 0
extern "C" { 
  void vaxpy3(Vector *res,Float *scale,Vector *mult,Vector *add, int ncvec);
  void vaxpy3_norm(Vector *res,Float *scale,Vector *mult,Vector *add, int ncvec,Float *norm);
}
#endif
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
int DiracOp::InvCg(Vector *out, 
		   Vector *in, 
		   Float src_norm_sq, 
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
  int i;
  char *fname = "InvCg(V*,V*,F,F*)";

// Flash the LED and then turn it off
//------------------------------------------------------------------
  VRB.LedFlash(cname,fname,3);
  VRB.LedOff(cname,fname);
  VRB.Func(cname,fname);


// Print out input parameters
//------------------------------------------------------------------
  VRB.Input(cname,fname,
	    "stop_rsd = %e\n",IFloat(dirac_arg->stop_rsd));
  VRB.Input(cname,fname,
	    "max_num_iter = %d\n",dirac_arg->max_num_iter);
  VRB.Input(cname,fname,
	    "mass = %e\n",IFloat(dirac_arg->mass));
  VRB.Input(cname,fname,
	    "src_norm_sq = %e\n",IFloat(src_norm_sq));


//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------

// Set the source vector pointer
//------------------------------------------------------------------
  Vector *src = in;
  IFloat *src_tmp = (IFloat *)src;

// Set the solution vector pointer
//------------------------------------------------------------------
  Vector *sol = out;

// Set the node checkerboard size of the fermion field
//------------------------------------------------------------------

  if(lat.Fclass() == F_CLASS_CLOVER) {
    f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  } else {
    f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / (lat.FchkbEvl()+1);
  }
    
// Allocate memory for the residual vector res.
//------------------------------------------------------------------
  Vector *res = (Vector *) qalloc(QCOMMS|QFAST,f_size_cb * sizeof(Float));
  if(res == 0){
    res = (Vector *) qalloc(QCOMMS,f_size_cb * sizeof(Float));
    printf("res=%p\n",res);
  }
  if(res == 0)
    ERR.Pointer(cname,fname, "res");
  VRB.Smalloc(cname,fname, "res", res, f_size_cb * sizeof(Float));

// Allocate memory for the direction vector dir.
//------------------------------------------------------------------
  Vector *dir;
  if(GJP.VolNodeSites() >4096) dir=0;
    else dir = (Vector *) qalloc(QCOMMS|QFAST,f_size_cb * sizeof(Float));
  if(dir == 0){
    dir = (Vector *) qalloc(QCOMMS,f_size_cb * sizeof(Float));
  }
  if(dir == 0)
    ERR.Pointer(cname,fname, "dir");
  VRB.Smalloc(cname,fname, "dir", dir, f_size_cb * sizeof(Float));

// Allocate mem. for the result vector of matrix multiplication mmp.
//------------------------------------------------------------------
  Vector *mmp;
  if(GJP.VolNodeSites() >4096) mmp=0;
    else mmp = (Vector *) qalloc(QCOMMS|QFAST,f_size_cb * sizeof(Float));
  if(mmp == 0){
    mmp = (Vector *) qalloc(QCOMMS,f_size_cb * sizeof(Float));
  }
  if(mmp == 0)
    ERR.Pointer(cname,fname, "mmp");
  VRB.Smalloc(cname,fname, "mmp", mmp, f_size_cb * sizeof(Float));

// If src_norm_sq is not provided calculate it
//------------------------------------------------------------------
  if(src_norm_sq == 0){
    src_norm_sq = src->NormSqNode(f_size_cb);
    DiracOpGlbSum(&src_norm_sq);
  }

// Calculate stopping condition
//------------------------------------------------------------------
  stp_cnd = src_norm_sq * dirac_arg->stop_rsd * dirac_arg->stop_rsd;
  VRB.Flow(cname,fname, "stp_cnd =%e\n", IFloat(stp_cnd));

#ifdef REPRODUCE_TEST 
  
// Allocate space for storing solution
//------------------------------------------------------------------
  Vector *sol_store = (Vector *) smalloc(f_size_cb * sizeof(Float));
  if(sol_store == 0) ERR.Pointer(cname,fname, "sol_store");
  VRB.Smalloc(cname,fname, "sol_store", sol_store, f_size_cb * sizeof(Float));

// Allocate space for storing d
//------------------------------------------------------------------
  Float *d_store = (Float *) smalloc( dirac_arg->max_num_iter-1 * sizeof(Float));

  if(d_store == 0) ERR.Pointer(cname,fname, "d_store");
  VRB.Smalloc(cname,fname, "d_store", d_store, dirac_arg->max_num_iter-1 * sizeof(Float));

  for ( int n = 0; n < dirac_arg->max_num_iter-1; n++ )  d_store[n] = 0;

  sol_store->CopyVec(sol, f_size_cb);

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
  MatPcDagMatPc(mmp, sol);

  // res = src
  res->CopyVec(src, f_size_cb);

  // res -= mmp
  res->VecMinusEquVec(mmp, f_size_cb);

  // dir = res
  dir->CopyVec(res, f_size_cb);  

  // res_norm_sq_cur = res * res
  res_norm_sq_cur = res->NormSqNode(f_size_cb);

  DiracOpGlbSum(&res_norm_sq_cur);

  // if( |res|^2 <= stp_cnd ) we are done
  VRB.Flow(cname,fname,
  	   "|res[0]|^2 = %e\n", IFloat(res_norm_sq_cur));
  itr = 0;
  max_itr = dirac_arg->max_num_iter-1;
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
    timeval start,end;
    itr++;
    res_norm_sq_prv = res_norm_sq_cur;

    // mmp = MatPcDagMatPc * dir
    // d = <dir, MatPcDagMatPc*dir>
    MatPcDagMatPc(mmp, dir, &d);

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
  
    DiracOpGlbSum(&d);

    // If d = 0 we are done
    if(d == 0.0) {
      break;
      //??? or should we give a warning or error? Yes we should, really.
    }

    a = res_norm_sq_prv / d;

    // sol = a * dir + sol;
    //sol->FTimesV1PlusV2(a, dir, sol, f_size_cb);

    vaxpy3(sol,&a,dir,sol,f_size_cb/6);


#ifdef PROFILE
//    nflops_tmp +=f_size_cb*2;
#endif
    CGflops+=f_size_cb*2;

    //PAB. Should be able to use "vaxpy_norm here"
    // res = - a * (MatPcDagMatPc * dir) + res;
    // res_norm_sq_cur = res * res

   // a *= -1.0;
  double ma = -a;
    vaxpy3_norm(res,&ma,mmp,res,f_size_cb/6,&res_norm_sq_cur);
    DiracOpGlbSum(&res_norm_sq_cur);
#ifdef PROFILE
//    nflops_tmp +=f_size_cb*4;
#endif
    CGflops+=f_size_cb*4;


    // if( |res|^2 <= stp_cnd ) we are done
    VRB.Flow(cname,fname, "|res[%d]|^2 = %e\n", itr, IFloat(res_norm_sq_cur));
    if(res_norm_sq_cur <= stp_cnd) break;

    b = res_norm_sq_cur / res_norm_sq_prv;

    // dir = b * dir + res;
    vaxpy3(dir,&b,dir,res,f_size_cb/6);
#ifdef PROFILE
//    linalg_start = linalg_tmp;
//    gettimeofday(&linalg_end,NULL);
//    nflops =nflops_tmp+f_size_cb*2;
#endif
    CGflops+=f_size_cb*2;
  }

#ifdef PROFILE
    gettimeofday(&end,NULL);

    report_flops(CGflops,&start,&end); 
#endif

  // It has not reached stp_cnd: Issue a warning
  if(itr == dirac_arg->max_num_iter - 1){
    VRB.Warn(cname,fname, "CG reached max iterations = %d. |res|^2 = %e\n",
	     itr+1, IFloat(res_norm_sq_cur) );
  }

//------------------------------------------------------------------
// Done. Finish up and return
//------------------------------------------------------------------
  // Calculate and set true residual: 
  // true_res = |src - MatPcDagMatPc * sol| / |src|
  MatPcDagMatPc(mmp, sol);
  res->CopyVec(src, f_size_cb);
  res->VecMinusEquVec(mmp, f_size_cb);
  res_norm_sq_cur = res->NormSqNode(f_size_cb);
  DiracOpGlbSum(&res_norm_sq_cur);
  Float tmp = res_norm_sq_cur / src_norm_sq;
  tmp = sqrt(tmp);
  if(true_res != 0){
    *true_res = tmp;
  }
  VRB.Result(cname,fname,
	     "True |res| / |src| = %e, iter = %d\n", IFloat(tmp), itr+1);

#ifdef REPRODUCE_TEST 
  }
  VRB.Sfree(cname, fname,"d_store", d_store);
  sfree(d_store);
  VRB.Sfree(cname, fname,"sol_store", sol_store);
  sfree(sol_store);

#endif

  // Free memory
  VRB.Sfree(cname,fname, "mmp", mmp);
  qfree(mmp);
  VRB.Sfree(cname,fname, "dir", dir);
  qfree(dir);
  VRB.Debug("b ============\n");
  VRB.Sfree(cname,fname, "res", res);
  sfree(res);

  VRB.Debug("a ============\n");

// Flash the LED and then turn it on
//------------------------------------------------------------------
  VRB.FuncEnd(cname,fname);
  VRB.LedFlash(cname,fname,2);
  VRB.LedOn(cname,fname);

  // Return number of iterations
  return itr+1;

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
  \return The number of solver iterations.
  \post \a f_out contains the solution vector.
*/
//------------------------------------------------------------------
int DiracOp::InvCg(Vector *out, Vector *in, Float src_norm_sq)
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
int DiracOp::InvCg(Vector *out, Vector *in, Float *true_res)
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
int DiracOp::InvCg(Vector *out, Vector *in)
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

CPS_END_NAMESPACE
