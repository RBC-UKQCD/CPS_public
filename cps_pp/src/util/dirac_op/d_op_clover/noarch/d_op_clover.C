#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_clover/noarch/d_op_clover.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: d_op_clover.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.7  2002/03/11 22:27:03  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.4.2.1  2002/03/08 16:36:37  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.4  2001/08/16 10:50:16  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:12:37  anj
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
//  $RCSfile: d_op_clover.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_clover/noarch/d_op_clover.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

//----------------------------------------------------------------------
//
// d_op_clover.C
//
// DiracOpClover is derived from the DiracOp base class. 
// DiracOpClover is the front end for a library that contains
// all Dirac operators associated with Clover fermions.
//
//----------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/dirac_op.h>
CPS_START_NAMESPACE


CPS_END_NAMESPACE
#include<util/lattice.h>
#include<util/vector.h>
#include<util/enum.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/wilson.h>
#include<util/clover.h>
#include<mem/p2v.h>
CPS_START_NAMESPACE

//----------------------------------------------------------------------
// Constructor
//----------------------------------------------------------------------
DiracOpClover::DiracOpClover(Lattice & latt,
			     Vector *f_field_out,
			     Vector *f_field_in,
			     CgArg *arg,
			     CnvFrmType cnv_frm_flg) 
  : DiracOpWilsonTypes(latt, f_field_out, f_field_in, arg, cnv_frm_flg)
{
  cname = "DiracOpClover";
  char *fname = "DiracOpClover(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func(cname,fname);

  // Do the necessary conversions
  //--------------------------------------------------------------------
  if(cnv_frm == CNV_FRM_YES)
    lat.Convert(WILSON, f_out, f_in);
  else if(cnv_frm == CNV_FRM_NO)
    lat.Convert(WILSON);

  // Initialize parameters
  //--------------------------------------------------------------------
  DiracArg(dirac_arg);


  // Initialize the pointer to the initialized Clover structure
  // (the structure has been initialized by the Lattice::Fclover
  // constructor.
  //--------------------------------------------------------------------
  clover_lib_arg = (Clover *)(lat.FdiracOpInitPtr());


  // Calculate the clover matrices for the ODD checkerboard
  //--------------------------------------------------------------------
  CloverMatChkb(CHKB_ODD, 0);


  // Calculate the inversed clover matrices for the EVEN checkerboard
  //--------------------------------------------------------------------
  CloverMatChkb(CHKB_EVEN, 1);

  // Copy optimized code into its execution place (CRAM)
  //----------------------------------------------------------------
  p2vCloverLib();
}


//----------------------------------------------------------------------
// Destructor
//----------------------------------------------------------------------
DiracOpClover::~DiracOpClover() {
  char *fname = "~DiracOpClover()";
  VRB.Func(cname,fname);

  //--------------------------------------------------------------------
  // Do the necessary conversions
  //--------------------------------------------------------------------
  if(cnv_frm == CNV_FRM_YES)
    lat.Convert(CANONICAL, f_out, f_in);
  else if(cnv_frm == CNV_FRM_NO)
    lat.Convert(CANONICAL);

}


//----------------------------------------------------------------------
// DiracArg(CgArg *arg):
// It sets the dirac_arg pointer to arg and initializes
// kappa.
//----------------------------------------------------------------------
// Modified in by Ping for anisotropic lattices and clover improvement
//----------------------------------------------------------------------
void DiracOpClover::DiracArg(CgArg *arg){
  dirac_arg = arg;

  // The following formulae holds only for isotropic lattices
  //
  //  kappa = 1.0 / (2.0 * (dirac_arg->mass + 4.0));
  //  omega = omega_xi = kappa * clover_lib_arg->clover_coef; 


  // kappa = 1.0 / [2((dirac_arg->mass + vel_t) xi + 3 vel_s)]
  kappa = 0.5 / ((dirac_arg->mass + GJP.XiVXi()) * GJP.XiBare() + 
                 3*GJP.XiV());
 
  // omega = Csw_s / [2((dirac_arg->mass + 1.0) xi + 3 vel)]
  omega = kappa * GJP.CloverCoeff();
 
  // kappa = vel / [2((dirac_arg->mass + 1.0) xi + 3 vel)]
  kappa *= GJP.XiV();  
 
  // omega_xi = Csw_t vel^2 / [2 xi ((dirac_arg->mass + 1.0) xi + 3 vel)]
  omega_xi = kappa * GJP.CloverCoeffXi() * GJP.XiV();
  if (omega_xi != 0.0) 
    omega_xi /=  (GJP.XiBare() * GJP.XiVXi() * GJP.XiVXi());
}


 
//------------------------------------------------------------------
// Dslash(Vector *out, Vector *in, ChkbType cb, DagType dag) :
// Not implemented.
//------------------------------------------------------------------
void DiracOpClover::Dslash(Vector *out, 
                           Vector *in, 
                           ChkbType cb, 
                           DagType dag) {
  ERR.NotImplemented(cname,"Dslash");
}


//----------------------------------------------------------------------
// MatPcDagMatPc :
// MatPcDagMatPc is the Hermitian matrix M^dag M, where M is
// the even/odd preconditioned Dirac Operator matrix.        
// MatPcDagMatPc connects only odd-->odd sites.
// The in, out fields are defined on the odd checkerboard.
// If dot_prd is not 0 then the dot product (on node)
// <out, in> = <MatPcDagMatPc*in, in> = <MatPc*in, MatPc*in>
//----------------------------------------------------------------------
void DiracOpClover::MatPcDagMatPc(Vector *out, 
				  Vector *in, 
				  Float *dot_prd) 
{
  char *fname = "MatPcDagMatPc(V*,V*,F*)";
  VRB.Func(cname,fname);

  // use clover_lib_arg->frm_buf1 as local buffer
  //--------------------------------------------------------------------
  Vector *frm_buf1 = (Vector *)(clover_lib_arg->frm_buf1);  

  // frm_buf1 = MatPc in
  //--------------------------------------------------------------------
  MatPc(frm_buf1, in);

  // *dot_prd = <frm_buf1, frm_buf1>
  //--------------------------------------------------------------------
  if (dot_prd) {
    int vec_size = lat.FsiteSize()*GJP.VolNodeSites()/2;  
    *dot_prd = frm_buf1->NormSqNode(vec_size);    
  }

  // out = MatPcDag MatPc in
  //--------------------------------------------------------------------
  MatPcDag(out, frm_buf1);
}

//----------------------------------------------------------------------
// MatPc(Vector *out, Vector *in) :
// MatPc is the preconditioned fermion matrix.  
// MatPc connects only odd-->odd sites.
// The in, out fields are defined on the odd checkerboard.
//   out = (Aoo - kappa*kappa*Doe Aee^inv Deo) in
// where Aoo and Aee are hermitian matrices.
//----------------------------------------------------------------------
// Local buffer:
//    clover_lib_arg->frm_buf0
//----------------------------------------------------------------------
void DiracOpClover::MatPcDagOrNot(Vector *out, 
				  const Vector *in, int dag) const 
{
  char *fname = "MatPcDagOrNot";
  VRB.Func(cname,fname);

  // use clover_lib_arg->frm_buf0 as local buffer
  //--------------------------------------------------------------------
  Vector *frm_buf0 = (Vector *)(clover_lib_arg->frm_buf0);  
  const int half_sites = GJP.VolNodeSites()/2;

  // frm_buf0 = Doe Aee^inv Deo in
  //--------------------------------------------------------------------
  {
    Wilson* wilson_p = clover_lib_arg->wilson_p;
    wilson_dslash((IFloat *)frm_buf0, (IFloat *)gauge_field, (IFloat *)in, 
		  1, dag, wilson_p); 
    clover_mat_mlt((IFloat *)out, 
		   (const IFloat *)(lat.Aux0Ptr()), 
		   (const IFloat *)frm_buf0, 
		   half_sites);    
    wilson_dslash((IFloat *)frm_buf0, (IFloat *)gauge_field, (IFloat *)out, 
		  0, dag, wilson_p); 
  }
  
  // out = (Aoo - kappa*kappa Doe Aee^inv Deo) in
  //--------------------------------------------------------------------
  {
    clover_mat_mlt((IFloat *)out, 
		   (const IFloat *)(lat.Aux1Ptr()), 
		   (const IFloat *)in, 
		   half_sites);    
    fTimesV1PlusV2((IFloat *)out, 
		   -kappa*kappa, 
		   (const IFloat *)frm_buf0, 
		   (const IFloat *)out,                // out here OK
		   lat.FsiteSize()*half_sites);    
  }
  
}



//----------------------------------------------------------------------
// int MatInv(Vector *out, Vector *in, Float *true_res);
//----------------------------------------------------------------------
// Purpose:
//      calculates out where A * out = in 
// Arguments:
//    A:   the fermion matrix (Dirac operator) with no preconditioning.
//    in:  the fermion field source vector, defined on the whole lattice.
//         unchanged on return.
//    out: the initial guess, defined on the whole lattice.
//         on return is the solution.
//    true_res:  if not 0,
//         on return holds the value of the true residual of the CG.
//         *true_res = |src - MatPcDagMatPc * sol| / |src|
//    return:
//         the total number of CG iterations.
// Algorithm:
//    i.  The preconditioned matrix is inverted with
//        i.1: on the odd checkerboard:  the conjugate gradient 
//        i.2: on the even checkerboard: explicit inversion of 
//                                       those local and hermitian
//                                       clover matrices.
//    ii. steps: (no extra storage required)
//        1. out_even0 = Aee^inv in_even
//        2. out_odd  = (MatPcDagMatPc)^inv MatPcDag 
//                      (kappa Doe out_even0 + in_odd)
//        3. out_even  = Aee^inv (kappa Deo out_odd + in_even)
//----------------------------------------------------------------------
int DiracOpClover::MatInv(Vector *out, Vector *in, Float *true_res,
			  PreserveType prs_in) 
{
  char *fname = "MatInv(V*,V*,F*)";
  VRB.Func(cname,fname);

  // Initializations
  //--------------------------------------------------------------------
  int half_sites = GJP.VolNodeSites()/2;    
  int vec_size = lat.FsiteSize() * half_sites;    
  IFloat *in_even = (IFloat *)in + vec_size;
  IFloat *out_even = (IFloat *)out+ vec_size;  
  IFloat *A_even = (IFloat *)lat.Aux0Ptr();

  // out_even = out_even0 = Aee^inv even_in
  //--------------------------------------------------------------------
  clover_mat_mlt(out_even, A_even, in_even, half_sites);

  // use clover_lib_arg->frm_buf1 as local buffer
  //--------------------------------------------------------------------
  Vector *frm_buf1 = (Vector *)(clover_lib_arg->frm_buf1);  
  Wilson* wilson_p = clover_lib_arg->wilson_p;

  // frm_buf1 = kappe Doe out_even0 + in_odd
  //--------------------------------------------------------------------
  wilson_dslash((IFloat *)frm_buf1, (IFloat *)gauge_field,  
		(IFloat *)out_even, 0, 0, 
		wilson_p);                   // frm_buf1 = Doe out_even0
  fTimesV1PlusV2((IFloat *)frm_buf1, 
		 kappa, 
		 (const IFloat *)frm_buf1, 
		 (const IFloat *)in,     
		 vec_size);


  // out_even = MatPcDag [kappe * Doe out_even0 + in_odd]
  //--------------------------------------------------------------------
  MatPcDag((Vector *)out_even, frm_buf1);


  // out_odd = (MatPcDagMatPc)^inv MatPcDag 
  //           [kappe * Doe out_even0 + in_odd]           done!
  //--------------------------------------------------------------------
  int iter = InvCg(out, (Vector *)out_even, true_res);


  // frm_buf1 = kappa Deo out_odd + in_even
  //--------------------------------------------------------------------
  wilson_dslash((IFloat *)out_even, (IFloat *)gauge_field, 
		(IFloat *)out, 1, 0, 
		wilson_p);          // out_even = Deo out_odd
  fTimesV1PlusV2((IFloat *)frm_buf1, 
		 kappa, 
		 (const IFloat *)out_even, 
		 (const IFloat *)in_even,     
		 vec_size);

  // out_even = Aee^inv (kappa Deo out_odd + in_even)     done!
  //--------------------------------------------------------------------
  clover_mat_mlt((IFloat *)out_even, A_even, 
		 (const IFloat*)frm_buf1, half_sites);

  VRB.FuncEnd(cname,fname);
  return iter;
}


//------------------------------------------------------------------
// Overloaded function is same as original 
// but true_res=0.
//------------------------------------------------------------------
int DiracOpClover::MatInv(Vector *out, Vector *in, PreserveType prs_in)
{ return MatInv(out, in, 0, prs_in); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in and out = f_out.
//------------------------------------------------------------------
int DiracOpClover::MatInv(Float *true_res, PreserveType prs_in)
{ return MatInv(f_out, f_in, true_res, prs_in); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in, out = f_out, true_res=0.
//------------------------------------------------------------------
int DiracOpClover::MatInv(PreserveType prs_in)
{ return MatInv(f_out, f_in, 0, prs_in); }


//----------------------------------------------------------------------
// int MatEvlInv(Vector *out, Vector *in, Float *true_res);
//----------------------------------------------------------------------
// Purpose:
//      calculates out where A * out = in 
// Arguments:
//    A:   the preconditioned fermion matrix (Dirac operator) that
//         appears in the HMC evolution.
//    in:  the fermion field source vector, defined on the whole lattice.
//         unchanged on return.
//    out: the initial guess, defined on the whole lattice.
//         on return is the solution.
//    true_res:  if not 0,
//         on return holds the value of the true residual of the CG.
//         *true_res = |src - MatPcDagMatPc * sol| / |src|
//    return:
//         the total number of CG iterations.
// Algorithm:
//   The preconditioned matrix is inverted with
//     1: on the odd checkerboard:  the conjugate gradient 
//     2: on the even checkerboard: explicit inversion of 
//                                 those local and hermitian clover matrices.
//----------------------------------------------------------------------
int DiracOpClover::MatEvlInv(Vector *out, Vector *in, Float *true_res) 
{
  char *fname = "MatEvlInv(V*,V*,F*)";
  VRB.Func(cname,fname);

  // out_even = Aee^inv Aee^inv in_even                      done!
  //--------------------------------------------------------------------
  {
    int half_sites = GJP.VolNodeSites()/2;    
    int vec_size = lat.FsiteSize() * half_sites;    
    IFloat *in_even = (IFloat *)in + vec_size;
    IFloat *out_even = (IFloat *)out + vec_size;  
    IFloat *A_even = (IFloat *)lat.Aux0Ptr();
    IFloat *frm_buf0 = (IFloat *)(clover_lib_arg->frm_buf0);  
    clover_mat_mlt(frm_buf0, A_even, in_even, half_sites);
    clover_mat_mlt(out_even, A_even, frm_buf0, half_sites);
  }

  // out_odd = (MatPcDagMatPc)^inv in_odd                    done!
  //--------------------------------------------------------------------
  return InvCg(out,in,true_res);
}

//------------------------------------------------------------------
// Overloaded function is same as original
// but true_res=0.
//------------------------------------------------------------------
int DiracOpClover::MatEvlInv(Vector *out, Vector *in)
{ return MatEvlInv(out, in, 0); }


//------------------------------------------------------------------
// Overloaded function is same as original
// but in = f_in and out = f_out.
//------------------------------------------------------------------
int DiracOpClover::MatEvlInv(Float *true_res)
{ return MatEvlInv(f_out, f_in, true_res); }


//------------------------------------------------------------------
// Overloaded function is same as original
// but in = f_in, out = f_out, true_res=0.
//------------------------------------------------------------------
int DiracOpClover::MatEvlInv(void)
{ return MatEvlInv(f_out, f_in, 0); }


//------------------------------------------------------------------
// Mat(Vector *out, Vector *in) :
// Mat is the unpreconditioned fermion matrix.  
// Mat works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpClover::Mat(Vector *out, Vector *in) {  
  ERR.NotImplemented(cname,"Mat");
}


//------------------------------------------------------------------
// MatDag(Vector *out, Vector *in) :
// MatDag is the unpreconditioned fermion matrix.  
// MatDag works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpClover::MatDag(Vector *out, Vector *in) {
  ERR.NotImplemented(cname,"MatDag");
}


//------------------------------------------------------------------
// MatHerm(Vector *out, Vector *in) :
// MatHerm is gamma5*Mat.
// MatHerm works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpClover::MatHerm(Vector *out, Vector *in) {
  char *fname = "MatHerm(V*,V*)";
  VRB.Func(cname,fname);

  int temp_size = GJP.VolNodeSites() * lat.FsiteSize();
  Vector *temp = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp == 0) 
    ERR.Pointer(cname, fname, "temp");
  VRB.Smalloc(cname,fname, "temp", temp, temp_size);

  Mat(temp, in);
  MultGamma(out, temp, 15, GJP.VolNodeSites());
  
  VRB.Sfree(cname, fname, "temp", temp);
  sfree(temp);
}


//------------------------------------------------------------------
// Lingling changed the prvious defintion of the fermion fields by
// some constant factors. chi is the solution to MatPcInv on the 
// odd sites.  The user passes two full size
// CANONICAL fermion vectors with conversion enabled to the
// constructor.  Using chi, the function fills these vectors;
// the result may be used to compute the HMD fermion force.
//------------------------------------------------------------------

void DiracOpClover::CalcHmdForceVecs(Vector *chi)
{
  char *fname = "CalcHmdForceVecs(V*)" ;
  VRB.Func(cname,fname) ;

  if (f_out == 0)
    ERR.Pointer(cname, fname, "f_out") ;

  if (f_in == 0)
    ERR.Pointer(cname, fname, "f_in") ;

//------------------------------------------------------------------
// f_out stores (chi,rho), f_in stores (psi,sigma)
//------------------------------------------------------------------

  Vector *chi_new, *rho, *psi, *sigma ;

  int half_sites = GJP.VolNodeSites() / 2 ;
  int f_size_cb = lat.FsiteSize() * half_sites;

  chi_new = f_out ;

  chi_new->CopyVec(chi, f_size_cb) ;

  psi = f_in ;

  MatPc(psi,chi) ;

  rho = (Vector *)((Float *)f_out + f_size_cb) ;

  Vector *frm_buf0 = (Vector *)(clover_lib_arg->frm_buf0);

  wilson_dslash((IFloat *)frm_buf0, (IFloat *)gauge_field, (IFloat *)chi, 
		CHKB_ODD, DAG_NO, clover_lib_arg-> wilson_p) ;

  clover_mat_mlt((IFloat *)rho, 
		 (const IFloat *)(lat.Aux0Ptr()), 
		 (const IFloat *)frm_buf0, 
		 half_sites);
  
  rho->VecTimesEquFloat(kappa,f_size_cb) ;
  
  sigma = (Vector *)((Float *)f_in + f_size_cb) ;
  
  wilson_dslash((IFloat *)frm_buf0, (IFloat *)gauge_field, (IFloat *)psi, 
		CHKB_ODD, DAG_YES, clover_lib_arg-> wilson_p) ;

  clover_mat_mlt((IFloat *)sigma, 
		   (const IFloat *)(lat.Aux0Ptr()), 
		   (const IFloat *)frm_buf0, 
		   half_sites);

  sigma->VecTimesEquFloat(kappa,f_size_cb) ; 

  return ;
}

CPS_END_NAMESPACE
