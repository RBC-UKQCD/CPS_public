#include <config.h>
#ifdef USE_BLAS
#include <stdio.h>
#include <stdlib.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOp class CG solver methods.

  $Id: inv_cg_blas.C,v 1.3 2013-04-05 17:51:13 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:51:13 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_base/noarch/inv_cg_blas.C,v 1.3 2013-04-05 17:51:13 chulwoo Exp $
//  $Id: inv_cg_blas.C,v 1.3 2013-04-05 17:51:13 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: inv_cg_blas.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_base/noarch/inv_cg_blas.C,v $
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
#include <util/time_cps.h>
#include <util/qblas_extend.h>
//#include <comms/nga_reg.h>
#include <comms/cbuf.h>
#include <math.h>
#if TARGET == BGL
#include <sys/bgl/bgl_sys_all.h>
#endif
CPS_START_NAMESPACE

#ifdef  PARALLEL
//Uncomment the following line to activate reproducibility test
//#define REPRODUCE_TEST
#endif
#undef REPRODUCE_TEST

#define PROFILE

#ifdef  REPRODUCE_TEST
CPS_END_NAMESPACE
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE
#endif

static int bgl_cg_count = 0;

//------------------------------------------------------------------
// Circular buffer zero wait state access setting
//------------------------------------------------------------------
const unsigned CBUF_MODE4 = 0xcb18c1ff;

static int f_size_cb;     // Node checkerboard size of the fermion field

static inline void print_vec( Vector *vec, char *name){
  Float temp_f = vec->NormSqNode(f_size_cb);
  Float *temp_p = (Float *)vec;
  glb_sum(&temp_f);
  VRB.Flow("","print_vec()", "%s = %e %e \n", name,IFloat(temp_f),*temp_p);
}

inline       IFloat* IFl( Vector* x ){ return reinterpret_cast<IFloat*>(x); }
inline const IFloat* IFl( const Vector* x ){ 
  return reinterpret_cast<const IFloat*>(x); 
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
int DiracOp::InvCg(Vector *out, 
		   Vector *in, 
		   Float src_norm_sq, 
		   Float *true_res){
  char *fname = "InvCg(V*,V*,F,F*)";
  VRB.Func(cname,fname);

  int itr;              // Current number of CG iterations
  int max_itr;          // Max number of CG iterations
  Float stp_cnd;        // Stop if residual^2 <= stp_cnd
  Float res_norm_sq_prv;// The previous step |residual|^2
  Float res_norm_sq_cur;// The current step |residual|^2 
  Float a;
  Float b;
  Float d;
  int i, ic, icb;
  
  IFloat *temp;
  
  // Print out input parameters
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
  Vector *src = in;

  // Set the solution vector pointer
  Vector *sol = out;

  // Set the node checkerboard size of the fermion field
  if(lat.Fclass() == F_CLASS_CLOVER) {
    f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  } else {
    f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / (lat.FchkbEvl()+1);
  }
  
  // Allocate memory for the residual vector res.
  Vector *res = (Vector *) smalloc(f_size_cb * sizeof(Float));
  if(res == 0)
    ERR.Pointer(cname,fname, "res");
  VRB.Smalloc(cname,fname, "res", res, f_size_cb * sizeof(Float));
  
  // Allocate memory for the direction vector dir.
  //------------------------------------------------------------------
  Vector *dir = (Vector *) smalloc(f_size_cb * sizeof(Float));
  if(dir == 0)
    ERR.Pointer(cname,fname, "dir");
  VRB.Smalloc(cname,fname, "dir", dir, f_size_cb * sizeof(Float));
  
  // Allocate mem. for the result vector of matrix multiplication mmp.
  //------------------------------------------------------------------
  Vector *mmp = (Vector *) smalloc(f_size_cb * sizeof(Float));
  if(mmp == 0)
    ERR.Pointer(cname,fname, "mmp");
  VRB.Smalloc(cname,fname, "mmp", mmp, f_size_cb * sizeof(Float));
  
  // If src_norm_sq is not provided calculate it
  //------------------------------------------------------------------
  if(src_norm_sq == 0){
    src_norm_sq = cblas_ddot(f_size_cb,IFl(src));
    DiracOpGlbSum(&src_norm_sq);
  }
  VRB.Flow(cname,fname,"src_norm_sq=%e\n",src_norm_sq);

  // Calculate stopping condition
  //------------------------------------------------------------------
  stp_cnd = src_norm_sq * dirac_arg->stop_rsd * dirac_arg->stop_rsd;
  VRB.Flow(cname,fname,"stp_cnd =%e\n",IFloat(stp_cnd));
  
#ifdef REPRODUCE_TEST 
  
  // Allocate space for storing solution
  //------------------------------------------------------------------
  Vector *sol_store = (Vector *) smalloc(f_size_cb * sizeof(Float));
  if(sol_store == 0)
    ERR.Pointer(cname,fname, "sol_store");
  VRB.Smalloc(cname,fname, "sol_store", sol_store, f_size_cb * sizeof(Float));
  
  // Allocate space for storing d
  //------------------------------------------------------------------
  Float *d_store = (Float *) smalloc( (dirac_arg->max_num_iter-1) * sizeof(Float));
  if(d_store == 0)
    ERR.Pointer(cname,fname, "d_store");
  VRB.Smalloc(cname,fname, "d_store", d_store, (dirac_arg->max_num_iter-1) * sizeof(Float));
  
  for ( int n = 0; n < dirac_arg->max_num_iter-1; n++ )  d_store[n] = 0;
  
  cblas_dcopy(f_size_cb,IFl(sol),IFl(sol_store));
  
  for ( int test = 0; test < 2; test++ ) {
    if (test == 1) {
      cblas_dcopy(f_size_cb,IFl(sol_store),IFl(sol));
    }
  
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
    print_vec( mmp, "mmp");
    
    // res = src
    cblas_dcopy(f_size_cb,IFl(src),IFl(res));
    print_vec( res, "res");
    
    // res -= mmp
    cblas_daxpy(f_size_cb,-1,IFl(mmp),IFl(res));

    print_vec( res, "res");
    
    // dir = res
    cblas_dcopy(f_size_cb,IFl(res),IFl(dir));
    print_vec( dir, "dir");
    
    // res_norm_sq_cur = res * res
    //res_norm_sq_cur = res->NormSqNode(f_size_cb); //!BLAS
    res_norm_sq_cur = cblas_ddot(f_size_cb,IFl(res));

    //printf("res_norm_sq_cur=%e\n",res_norm_sq_cur);
    DiracOpGlbSum(&res_norm_sq_cur);
    
    // if( |res|^2 <= stp_cnd ) we are done
    VRB.Flow(cname,fname,
             "|res[0]|^2 = %e\n", IFloat(res_norm_sq_cur));
    itr = 0;
    max_itr = dirac_arg->max_num_iter-1;
    if(res_norm_sq_cur <= stp_cnd) max_itr = 0;
        
    //------------------------------------------------------------------
    // Loop over CG iterations
    //------------------------------------------------------------------
    double cg_time;
#ifdef PROFILE
    struct timeval start;
    struct timeval end;
    CGflops    = 0;
    gettimeofday(&start,NULL);
    cg_time = -dclock();
#if TARGET == BGL
    unsigned long long start_time = rts_get_timebase();
#endif
#endif
    
    for(i=0; i < max_itr; i++){
      
      itr = itr + 1; 
      res_norm_sq_prv = res_norm_sq_cur;
      
      // mmp = MatPcDagMatPc * dir
      // d = <dir, MatPcDagMatPc*dir>
      
      MatPcDagMatPc(mmp, dir, &d);
      
      print_vec( mmp, "mmp");
      
#ifdef REPRODUCE_TEST 
      
      /* Check reproducibility */
      if ( test == 0) d_store[ i ] = d;
      else if ( d != d_store[ i ] )
        InterruptExit(-1, "NODE FAILS TO REPRODUCE");
      /* End of Check */
      
#endif
      
      DiracOpGlbSum(&d);
      VRB.Flow(cname,fname, "d = %e\n", IFloat(d));
      
      // If d = 0 we are done
      if(d == 0.0) {
        ERR.General(cname,fname,"d(%e) = 0.0!!\n",d);
	exit(5);
        break;
        //??? or should we give a warning or error? Yes we should, really.
      }
      
      a = res_norm_sq_prv / d;
      VRB.Flow(cname,fname, "a = %e\n", IFloat(a));
      
      // sol = a * dir + sol;
      cblas_daxpy(f_size_cb,a,IFl(dir),IFl(sol));

      print_vec( sol, "sol");
      
      // res = - a * (MatPcDagMatPc * dir) + res;
      cblas_daxpy(f_size_cb,-a,IFl(mmp),IFl(res));
      print_vec( res, "res");
      
      // res_norm_sq_cur = res * res
      res_norm_sq_cur = cblas_ddot(f_size_cb,IFl(res));
      
      DiracOpGlbSum(&res_norm_sq_cur);
      
      // if( |res|^2 <= stp_cnd ) we are done
      VRB.Flow(cname,fname,
               "|res[%d]|^2 = %e\n", itr, IFloat(res_norm_sq_cur));

      if(res_norm_sq_cur <= stp_cnd) break;
      
      b = res_norm_sq_cur / res_norm_sq_prv;
      VRB.Flow(cname,fname, "b = %e\n", IFloat(b));
      
      // dir = b * dir + res;
      // slow!
      cblas_dscal(f_size_cb,b,IFl(dir));
      cblas_daxpy(f_size_cb,1,IFl(res),IFl(dir));

      print_vec( dir, "dir");
      CGflops+=f_size_cb*8;
      
    }

#ifdef PROFILE
    gettimeofday(&end,NULL);
    cg_time += dclock();
#if TARGET == BGL
    unsigned long long stop_time = rts_get_timebase();
    int inv_time = stop_time - start_time;
    Float perf = inv_time;
    perf = perf / (itr+1);
    perf = perf / GJP.VolNodeSites();
    perf = GJP.SnodeSites()* 678.0 * 100.0 / perf;
    
    //----------------------------------------------------------------
    // Performance reporting
    //----------------------------------------------------------------
    if(!UniqueID() && bgl_cg_count%1 == 0){
      printf("INVERTER TIME IN PCYCLES = %llu\n", stop_time - start_time);
      printf("INVERTER PERFORMANCE     = %3.1f\%\n", perf);
    }
    bgl_cg_count++;
#endif
    unsigned long long flops_per_site = CGflops;
    flops_per_site /= (GJP.VolNodeSites()*(itr+1));
//    print_flops(cname,fname,CGflops,&start,&end);
    print_flops(cname,fname,CGflops,cg_time);
    VRB.Result(cname,fname,"flops_per_site=%llu\n",flops_per_site);
#endif
    
    // It has not reached stp_cnd: Issue a warning
    if(itr == dirac_arg->max_num_iter - 1){
      VRB.Warn(cname,fname,
               "CG reached max iterations = %d. |res|^2 = %e\n",
               itr+1, IFloat(res_norm_sq_cur) );
    }
    
    //------------------------------------------------------------------
    // Done. Finish up and return
    //------------------------------------------------------------------
    // Calculate and set true residual: 
    // true_res = |src - MatPcDagMatPc * sol| / |src|
    MatPcDagMatPc(mmp, sol);

    cblas_dcopy(f_size_cb,IFl(src),IFl(res));
    cblas_daxpy(f_size_cb,-1,IFl(mmp),IFl(res));
    res_norm_sq_cur = cblas_ddot(f_size_cb,IFl(res));
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
  sfree(mmp);
  VRB.Sfree(cname,fname, "dir", dir);
  sfree(dir);
  VRB.Debug("b ============\n");
  VRB.Sfree(cname,fname, "res", res);
  sfree(res);
  VRB.Debug("a ============\n");
  
  VRB.FuncEnd(cname,fname);
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


CPS_END_NAMESPACE
#endif
