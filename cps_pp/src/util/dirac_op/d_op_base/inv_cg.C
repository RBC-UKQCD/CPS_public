#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_base/inv_cg.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: inv_cg.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.9  2002/03/11 22:27:02  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.6.2.1  2002/03/08 16:36:36  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.6  2001/08/16 10:50:15  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.4  2001/07/03 17:01:00  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.3  2001/06/21 15:40:16  anj
//  Updated the _TARTAN ifdefs, using PARALLEL instead (where appropriate).Anj
//
//  Revision 1.2  2001/06/19 18:12:36  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:05  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: inv_cg.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_base/inv_cg.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// inv_cg.C
//
//------------------------------------------------------------------
CPS_END_NAMESPACE
#include<config.h>
#include<util/dirac_op.h>
#include<util/lattice.h>
#include<util/smalloc.h>
#include<util/vector.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<comms/nga_reg.h>
#include<comms/cbuf.h>
#include <math.h>
CPS_START_NAMESPACE

#ifdef  PARALLEL
//Uncomment the following line to activate reproducibility test
//#define REPRODUCE_TEST
#endif

#ifdef  REPRODUCE_TEST
CPS_END_NAMESPACE
#include <sysfunc.h>
CPS_START_NAMESPACE
#endif

//------------------------------------------------------------------
// Circular buffer zero wait state access setting
//------------------------------------------------------------------
const unsigned CBUF_MODE4 = 0xcb18c1ff;

//------------------------------------------------------------------
//
// The Conjugate Gradient inverter.
// It inverts the matrix MatPcDagMatPc.
// This is the preconditioned (if relevant) matrix that
// is used in the HMC evolution.
//
// source is *in, initial guess and solution is *out.
//
// src_norm_sq is the norm of the source squared (dot product 
// of source * source). If src_norm_sq = 0 it is calculated
// inside InvCg
//
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
//
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
  int i, ic, icb;
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
    src_norm_sq = src->NormSqNode(f_size_cb);
    DiracOpGlbSum(&src_norm_sq);
  }

// Calculate stopping condition
//------------------------------------------------------------------
  stp_cnd = src_norm_sq * dirac_arg->stop_rsd * dirac_arg->stop_rsd;
  VRB.Flow(cname,fname, 
	   "stp_cnd =%e\n", IFloat(stp_cnd));

// Make IFloat pointers out of Vector pointers
//------------------------------------------------------------------
  IFloat *f_sol = (IFloat *) sol; 
  IFloat *f_dir = (IFloat *) dir; 
  IFloat *f_res = (IFloat *) res; 
  IFloat *f_mmp = (IFloat *) mmp; 

// Calculate the cram buffers size (must divide f_size_cb exactly)
//------------------------------------------------------------------
  int cram_buf_size = CRAM_SCRATCH_SIZE / 2;
  for(i=0; i< CRAM_SCRATCH_SIZE / 2; i++){
    cram_buf_size = cram_buf_size - i;
    if(f_size_cb % cram_buf_size == 0) break;
  }
  int cram_buf_size_sof = cram_buf_size * sizeof(Float);
  int cram_blocks = f_size_cb / cram_buf_size;

// Set pointers to two cram buffers
//------------------------------------------------------------------
#ifdef _TARTAN
  IFloat *cram_a = (IFloat *) CRAM_SCRATCH_ADDR;
  IFloat *cram_b = (IFloat *) (CRAM_SCRATCH_ADDR + cram_buf_size);
#else
  IFloat cram_a[CRAM_SCRATCH_SIZE/2];
  IFloat cram_b[CRAM_SCRATCH_SIZE/2];
#endif



#ifdef REPRODUCE_TEST 
  
// Allocate space for storing solution
//------------------------------------------------------------------
  Vector *sol_store = (Vector *) smalloc(f_size_cb * sizeof(Float));
  if(sol_store == 0)
    ERR.Pointer(cname,fname, "sol_store");
  VRB.Smalloc(cname,fname, "sol_store", sol_store, f_size_cb * sizeof(Float));

// Allocate space for storing d
//------------------------------------------------------------------
  Float *d_store = (Float *) smalloc( dirac_arg->max_num_iter-1 * sizeof(Float));

  if(d_store == 0)
    ERR.Pointer(cname,fname, "d_store");
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


//------------------------------------------------------------------
// Loop over CG iterations
//------------------------------------------------------------------
  for(i=0; i < max_itr; i++){

    itr = itr + 1;
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
  
    DiracOpGlbSum(&d);

    // If d = 0 we are done
    if(d == 0.0) {
      break;
      //??? or should we give a warning or error
    }

    a = res_norm_sq_prv / d;

    // Set circular buffer
    setCbufCntrlReg(4, CBUF_MODE4);

    // sol = a * dir + sol;
    // sol->FTimesV1PlusV2(a, dir, sol, f_size_cb);
    for(icb = 0; icb < cram_blocks; icb++){
      ic = icb * cram_buf_size;
      moveMem(cram_a, f_sol+ic+BANK4_BASE, cram_buf_size_sof);
      moveMem(cram_b, f_dir+ic+BANK4_BASE+BANK_SIZE, cram_buf_size_sof);
      fTimesV1PlusV2(f_sol+ic, a, cram_b, cram_a, cram_buf_size);
    }

    // res = - a * (MatPcDagMatPc * dir) + res;
    // res->FTimesV1PlusV2(-a, mmp, res, f_size_cb);
    for(icb = 0; icb < cram_blocks; icb++){
      ic = icb * cram_buf_size;
      moveMem(cram_a, f_res+ic+BANK4_BASE, cram_buf_size_sof);
      moveMem(cram_b, f_mmp+ic+BANK4_BASE+BANK_SIZE, cram_buf_size_sof);
      fTimesV1PlusV2(f_res+ic, -a, cram_b, cram_a, cram_buf_size);
    }

    // res_norm_sq_cur = res * res
    res_norm_sq_cur = res->NormSqNode(f_size_cb);
    DiracOpGlbSum(&res_norm_sq_cur);

    // if( |res|^2 <= stp_cnd ) we are done
    VRB.Flow(cname,fname,
	     "|res[%d]|^2 = %e\n", itr, IFloat(res_norm_sq_cur));
    if(res_norm_sq_cur <= stp_cnd) break;

    b = res_norm_sq_cur / res_norm_sq_prv;

    // dir = b * dir + res;
    // dir->FTimesV1PlusV2(b, dir, res, f_size_cb);
    for(icb = 0; icb < cram_blocks; icb++){
      ic = icb * cram_buf_size;
      moveMem(cram_a, f_res+ic+BANK4_BASE, cram_buf_size_sof);
      moveMem(cram_b, f_dir+ic+BANK4_BASE+BANK_SIZE, cram_buf_size_sof);
      fTimesV1PlusV2(f_dir+ic, b, cram_b, cram_a, cram_buf_size);
    }


  }

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
  sfree(mmp);
  VRB.Sfree(cname,fname, "dir", dir);
  sfree(dir);
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
// Same as original but with true_res=0
//------------------------------------------------------------------
int DiracOp::InvCg(Vector *out, Vector *in, Float src_norm_sq)
{ return InvCg(out, in, src_norm_sq, 0); }


//------------------------------------------------------------------
// Same as original but with src_norm_sq=0.0
//------------------------------------------------------------------
int DiracOp::InvCg(Vector *out, Vector *in, Float *true_res)
{ return InvCg(out, in, 0.0, true_res); }


//------------------------------------------------------------------
// Same as original but with src_norm_sq=0.0, true_res=0
//------------------------------------------------------------------
int DiracOp::InvCg(Vector *out, Vector *in)
{ return InvCg(out, in, 0.0, 0); }


//------------------------------------------------------------------
// Same as original but with 
// in=f_in, out=f_out
//------------------------------------------------------------------
int DiracOp::InvCg(Float src_norm_sq, Float *true_res)
{ return InvCg(f_out, f_in, src_norm_sq, true_res); }


//------------------------------------------------------------------
// Same as original but with 
// in=f_in, out=f_out, true_res=0
//------------------------------------------------------------------
int DiracOp::InvCg(Float src_norm_sq)
{ return InvCg(f_out, f_in, src_norm_sq, 0); }


//------------------------------------------------------------------
// Same as original but with 
// in=f_in, out=f_out, src_norm_sq=0.0
//------------------------------------------------------------------
int DiracOp::InvCg(Float *true_res)
{ return InvCg(f_out, f_in, 0.0, true_res); }


//------------------------------------------------------------------
// Same as original but with 
// in=f_in, out=f_out, src_norm_sq=0.0, true_res=0
//------------------------------------------------------------------
int DiracOp::InvCg(void)
{ return InvCg(f_out, f_in, 0.0, 0); }

CPS_END_NAMESPACE
