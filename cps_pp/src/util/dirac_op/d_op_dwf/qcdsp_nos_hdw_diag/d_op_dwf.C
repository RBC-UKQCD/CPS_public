#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:06 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdsp_nos_hdw_diag/d_op_dwf.C,v 1.2 2004-06-04 21:14:06 chulwoo Exp $
//  $Id: d_op_dwf.C,v 1.2 2004-06-04 21:14:06 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdsp_nos_hdw_diag/d_op_dwf.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// d_op_dwf.C
//
// DiracOpDwf is derived from the DiracOp base class. 
// DiracOpDwf is the front end for a library that contains
// all Dirac operators associated with Dwf fermions.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/dirac_op.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/wilson.h>
#include<util/dwf.h>
#include<mem/p2v.h>
#include<comms/glb.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
DiracOpDwf::DiracOpDwf(Lattice & latt,
			     Vector *f_field_out,
			     Vector *f_field_in,
			     CgArg *arg,
			     CnvFrmType cnv_frm_flg) :
			     DiracOpWilsonTypes(latt, 
						f_field_out,
						f_field_in, 
						arg,
						cnv_frm_flg)
{
  cname = "DiracOpDwf";
  char *fname = "DiracOpDwf(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func(cname,fname);


  //----------------------------------------------------------------
  // Do the necessary conversions
  //----------------------------------------------------------------
  if(cnv_frm == CNV_FRM_YES)
    lat.Convert(WILSON, f_out, f_in);
  else if(cnv_frm == CNV_FRM_NO)
    lat.Convert(WILSON);

  //----------------------------------------------------------------
  // Initialize parameters
  //----------------------------------------------------------------
  DiracArg(dirac_arg);

  //----------------------------------------------------------------
  // Initialize the pointer to the initialized Dwf structure
  // (the structure has been initialized by the Lattice::Fdwf
  // constructor.
  //----------------------------------------------------------------
  dwf_lib_arg = lat.FdiracOpInitPtr();

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // call to dwf_init but is done here again in case the user
  // has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the Fdwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Copy optimized code into its execution place (CRAM)
  //----------------------------------------------------------------
  p2vWilsonLib();

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
DiracOpDwf::~DiracOpDwf() {
  char *fname = "~DiracOpDwf()";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Do the necessary conversions
  //----------------------------------------------------------------
  if(cnv_frm == CNV_FRM_YES)
    lat.Convert(CANONICAL, f_out, f_in);
  else if(cnv_frm == CNV_FRM_NO)
    lat.Convert(CANONICAL);

}


//------------------------------------------------------------------
// DiracArg(CgArg *arg):
// It sets the dirac_arg pointer to arg and initializes
// mass.
//------------------------------------------------------------------
void DiracOpDwf::DiracArg(CgArg *arg){
  dirac_arg = arg;
  mass = dirac_arg->mass;
}


//------------------------------------------------------------------
// MatPcDagMatPc :
// MatPcDagMatPc is the fermion matrix that appears in the HMC 
// evolution. It is a Hermitian matrix.
// The in, out fields are defined on the checkerboard lattice.
// If dot_prd is not 0 then the dot product (on node)
// <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.
//------------------------------------------------------------------
void DiracOpDwf::MatPcDagMatPc(Vector *out, 
			       Vector *in, 
			       Float *dot_prd){

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  dwf_mdagm(out, 
	    gauge_field, 
	    in, 
	    dot_prd,
	    mass,
	    (Dwf *) dwf_lib_arg);
}


//------------------------------------------------------------------
// Dslash(Vector *out, Vector *in, ChkbType cb, DagType dag) :
// Dslash is the derivative part of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
// cb is not used.
//------------------------------------------------------------------
void DiracOpDwf::Dslash(Vector *out, 
			Vector *in, 
			ChkbType cb, 
			DagType dag) {

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  dwf_dslash(out, 
	     gauge_field, 
	     in, 
	     mass,
	     cb,
	     dag,
	     (Dwf *) dwf_lib_arg);
}

//------------------------------------------------------------------
// MatPc(Vector *out, Vector *in) :
// MatPc is the fermion matrix.  
// The in, out fields are defined on the checkerboard lattice.
//------------------------------------------------------------------
void DiracOpDwf::MatPc(Vector *out, Vector *in) {  

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  dwf_m(out, 
	gauge_field, 
	in, 
	mass,
	(Dwf *) dwf_lib_arg);
}

//------------------------------------------------------------------
// MatPcDag(Vector *out, Vector *in) :
// MatPcDag is the dagger of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice.
//------------------------------------------------------------------
void DiracOpDwf::MatPcDag(Vector *out, Vector *in) {

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  dwf_mdag(out, 
	   gauge_field, 
	   in, 
	   mass,
	   (Dwf *) dwf_lib_arg);
}


//------------------------------------------------------------------
// int MatInv(Vector *out, Vector *in, 
//            Float *true_res, PreserveType prs_in);
// The inverse of the unconditioned Dirac Operator 
// using Conjugate gradient.
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// prs_in is used to specify if the source
// in should be preserved or not. If not the memory usage
// is less by half the size of a fermion vector.
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int DiracOpDwf::MatInv(Vector *out, 
		       Vector *in, 
		       Float *true_res,
		       PreserveType prs_in) {
  char *fname = "MatInv(V*,V*,F*)";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  Vector *temp2;
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  // check out if converted
  //for (int ii = 0; ii < 2 * temp_size; ii++) {
  //  VRB.Result(cname, fname, "in[%d] = %e\n", ii, 
  //  *((IFloat *)in + ii));
  //  VRB.Result(cname, fname, "out[%d] = %e\n", ii, 
  //  *((IFloat *)out + ii));
  //}

  Vector *temp = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp == 0) ERR.Pointer(cname, fname, "temp");
  VRB.Smalloc(cname,fname, "temp", temp, temp_size * sizeof(Float));

  if(prs_in == PRESERVE_YES){
    temp2 = (Vector *) smalloc(temp_size * sizeof(Float));
    if (temp2 == 0) ERR.Pointer(cname, fname, "temp2");
    VRB.Smalloc(cname,fname, "temp2", temp2, temp_size * sizeof(Float));
  }

  // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (IFloat *) in + temp_size );

  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (IFloat *) out + temp_size );

  Dslash(temp, even_in, CHKB_EVEN, DAG_NO);

  fTimesV1PlusV2((IFloat *)temp, (IFloat) dwf_arg->dwf_kappa, (IFloat *)temp,
    (IFloat *)in, temp_size);

  // save source
  if(prs_in == PRESERVE_YES){
    moveMem((IFloat *)temp2, (IFloat *)in, 
		temp_size * sizeof(IFloat) / sizeof(char));
  }

  MatPcDag(in, temp);

  int iter = InvCg(out,in,true_res);

  // restore source
  if(prs_in == PRESERVE_YES){
    moveMem((IFloat *)in, (IFloat *)temp2, 
		temp_size * sizeof(IFloat) / sizeof(char));
  }

  Dslash(temp, out, CHKB_ODD, DAG_NO);

  fTimesV1PlusV2((IFloat *)even_out, (IFloat) dwf_arg->dwf_kappa, (IFloat *)temp,
    (IFloat *) even_in, temp_size);

  VRB.Sfree(cname, fname, "temp", temp);
  sfree(temp);

  if(prs_in == PRESERVE_YES){
    VRB.Sfree(cname, fname, "temp2", temp2);
    sfree(temp2);
  }

  return iter;
}


//------------------------------------------------------------------
// Overloaded function is same as original 
// but true_res=0.
//------------------------------------------------------------------
int DiracOpDwf::MatInv(Vector *out, Vector *in, PreserveType prs_in)
{ return MatInv(out, in, 0, prs_in); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in and out = f_out.
//------------------------------------------------------------------
int DiracOpDwf::MatInv(Float *true_res, PreserveType prs_in)
{ return MatInv(f_out, f_in, true_res, prs_in); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in, out = f_out, true_res=0.
//------------------------------------------------------------------
int DiracOpDwf::MatInv(PreserveType prs_in)
{ return MatInv(f_out, f_in, 0, prs_in); }


//------------------------------------------------------------------
// Mat(Vector *out, Vector *in) :
// Mat is the unpreconditioned fermion matrix.  
// Mat works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpDwf::Mat(Vector *out, Vector *in) {  
  char *fname = "Mat(V*,V*)";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (IFloat *) in + temp_size );
  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (IFloat *) out + temp_size );

  Dslash(out, even_in, CHKB_EVEN, DAG_NO);

  fTimesV1PlusV2((IFloat *)out, -((Dwf *)dwf_lib_arg)->dwf_kappa, 
		 (IFloat *)out, (IFloat *)in, temp_size);

  Dslash(even_out, in, CHKB_ODD, DAG_NO);
  
  fTimesV1PlusV2((IFloat *)even_out, -((Dwf *)dwf_lib_arg)->dwf_kappa, 
		 (IFloat *)even_out, (IFloat *)even_in, temp_size);

}


//------------------------------------------------------------------
// MatDag(Vector *out, Vector *in) :
// MatDag is the unpreconditioned fermion matrix.  
// MatDag works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpDwf::MatDag(Vector *out, Vector *in) {
  char *fname = "MatDag(V*,V*)";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (IFloat *) in + temp_size );
  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (IFloat *) out + temp_size );

  Dslash(out, even_in, CHKB_EVEN, DAG_YES);

  fTimesV1PlusV2((IFloat *)out, -((Dwf *) dwf_lib_arg)->dwf_kappa, 
		 (IFloat *)out, (IFloat *)in, temp_size);

  Dslash(even_out, in, CHKB_ODD, DAG_YES);
  
  fTimesV1PlusV2((IFloat *)even_out, -((Dwf *) dwf_lib_arg)->dwf_kappa,
		 (IFloat *)even_out, (IFloat *)even_in, temp_size);

}


//------------------------------------------------------------------
// MatHerm(Vector *out, Vector *in) :
// MatHerm is gamma5*R*Mat.
// MatHerm works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpDwf::MatHerm(Vector *out, Vector *in) {
  char *fname = "MatHerm(V*,V*)";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize();
  Vector *temp = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp == 0) 
    ERR.Pointer(cname, fname, "temp");
  VRB.Smalloc(cname,fname, "temp", temp, temp_size);

  Mat(out, in);
  lat.Freflex(temp, out);
  MultGamma(out, temp, 15, GJP.VolNodeSites()*GJP.SnodeSites());
  
  VRB.Sfree(cname, fname, "temp", temp);
  sfree(temp);

}


//------------------------------------------------------------------
// GRF
// chi is the solution to MatPcInv.  The user passes two full size
// CANONICAL fermion vectors with conversion enabled to the
// constructor.  Using chi, the function fills these vectors;
// the result may be used to compute the HMD fermion force.
//------------------------------------------------------------------

void DiracOpDwf::CalcHmdForceVecs(Vector *chi)
{
  char *fname = "CalcHmdForceVecs(V*)" ;
  VRB.Func(cname,fname) ;

  if (f_out == 0)
    ERR.Pointer(cname, fname, "f_out") ;

  if (f_in == 0)
    ERR.Pointer(cname, fname, "f_in") ;

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------

//------------------------------------------------------------------
// f_out stores (chi,rho), f_in stores (psi,sigma)
//------------------------------------------------------------------

  Vector *chi_new, *rho, *psi, *sigma ;

  int f_size_cb = 12 * GJP.VolNodeSites() * GJP.SnodeSites() ;

  chi_new = f_out ;

  chi_new->CopyVec(chi, f_size_cb) ;

  psi = f_in ;

  MatPc(psi,chi) ;

  {
    Float kappa = ((Dwf *)dwf_lib_arg)->dwf_kappa ;
    psi->VecTimesEquFloat(-kappa*kappa,f_size_cb) ;
  }

  rho = (Vector *)((Float *)f_out + f_size_cb) ;

  Dslash(rho, chi, CHKB_ODD, DAG_NO) ;

  sigma = (Vector *)((Float *)f_in + f_size_cb) ;

  Dslash(sigma, psi, CHKB_ODD, DAG_YES) ;

  return ;
}

//------------------------------------------------------------------
// DiracOpGlbSum(Float *): 
// The global sum used by InvCg. If s_nodes = 1
// it is the usual global sum. If s_nodes > 1 it
// is the 5-dimensional globals sum glb_sum_five.
//------------------------------------------------------------------
void DiracOpDwf::DiracOpGlbSum(Float *float_p) {
  if(GJP.Snodes() == 1) {
    glb_sum(float_p);
  }
  else {
    glb_sum_five(float_p);
  }
}

CPS_END_NAMESPACE
