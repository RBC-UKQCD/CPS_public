#include <config.h>
#include <stdio.h>
#include <alg/no_arg.h>
#include <alg/common_arg.h>
#include <alg/alg_plaq.h>

CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOpMobius class methods.

  $Id: d_op_mobius.C,v 1.6 2013-06-07 19:26:34 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-06-07 19:26:34 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/d_op_mobius.C,v 1.6 2013-06-07 19:26:34 chulwoo Exp $
//  $Id: d_op_mobius.C,v 1.6 2013-06-07 19:26:34 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: d_op_mobius.C,v $
//  $Revision: 1.6 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/d_op_mobius.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// d_op_mobius.C
//
// DiracOpMobius is derived from the DiracOp base class. 
// DiracOpMobius is the front end for a library that contains
// all Dirac operators associated with Dwf fermions.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/wilson.h>
#include <util/time_cps.h>
#include <util/dwf.h>
#include <util/mobius.h>
//#include <mem/p2v.h>
#include <comms/glb.h>

//#define USE_BLAS
#ifdef USE_BLAS
#include "noarch/blas-subs.h"
#endif

CPS_START_NAMESPACE



//------------------------------------------------------------------
/*!
  Only one instance of this class is allowed to be in existence at
  any time.
  \param latt The lattice on which this Dirac operator is defined
  \param f_field_out A (pointer to) a spin-colour field (optionally). 
  \param f_field_in A (pointer to) a spin-colour field (optionally).
  \param arg Parameters for the solver.
  \param cnv_frm_flag Whether the lattice fields should be converted to
  to a new storage order appropriate for the type of fermion action.
  If this is ::CNV_FRM_NO, then just the gauge field is converted.
  If this is ::CNV_FRM_YES, then the fields \a f_field_out and \a f_field_in
  are also converted: This assumes they are initially in the same order as
  the gauge field.
 */
//------------------------------------------------------------------
static  Matrix *new_gauge_field;
static  Matrix *old_gauge_field;
DiracOpMobius::DiracOpMobius(Lattice & latt,
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
  cname = "DiracOpMobius";
  char *fname = "DiracOpMobius(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func(cname,fname);


  //----------------------------------------------------------------
  // Do the necessary conversions
  //----------------------------------------------------------------
#undef PROFILE
#ifdef PROFILE
  Float time = -dclock();
#endif
  if(cnv_frm == CNV_FRM_YES)
    lat.Convert( DWF_4D_EOPREC_EE, f_out, f_in);
  else if(cnv_frm == CNV_FRM_NO)
    lat.Convert(WILSON);
#ifdef PROFILE
  time += dclock();
  print_flops("lattice","Convert()",0,time);
#endif

  //----------------------------------------------------------------
  // Initialize parameters
  //----------------------------------------------------------------
  DiracArg(dirac_arg);

  //----------------------------------------------------------------
  // Initialize the pointer to the initialized Mobius structure
  // (the structure has been initialized by the Lattice::Fmobius
  // constructor.
  //----------------------------------------------------------------
  mobius_lib_arg = lat.FdiracOpInitPtr();

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // call to mobius_init but is done here again in case the user
  // has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the Fmobius object.
  //----------------------------------------------------------------
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  mobius_arg->ls = GJP.SnodeSites();
  mobius_arg->mobius_kappa_b = 
    1.0 / ( 2 * (GJP.Mobius_b() *(4 - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
  mobius_arg->mobius_kappa_c = 
    1.0 / ( 2 * (GJP.Mobius_c() *(4 - GJP.DwfHeight()) - GJP.DwfA5Inv()) );

}


//------------------------------------------------------------------
/*!
  If the storage order of any fields was changed by the constructor
  then they are changed to the canonical order by the destructor.
*/
//------------------------------------------------------------------
DiracOpMobius::~DiracOpMobius() {
  char *fname = "~DiracOpMobius()";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Do the necessary conversions
  //----------------------------------------------------------------
#undef PROFILE
#ifdef PROFILE
  Float time = -dclock();
#endif
  if(cnv_frm == CNV_FRM_YES)
    lat.Convert(CANONICAL, f_out, f_in);
  else if(cnv_frm == CNV_FRM_NO)
    lat.Convert(CANONICAL);
#ifdef PROFILE
  time += dclock();
  print_flops("lattice","Convert()",0,time);
#endif

}


//------------------------------------------------------------------
// DiracArg(CgArg *arg):
// It sets the dirac_arg pointer to arg and initializes
// mass.
//------------------------------------------------------------------
void DiracOpMobius::DiracArg(CgArg *arg){
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
void DiracOpMobius::MatPcDagMatPc(Vector *out, 
				  Vector *in, 
				  Float *dot_prd){

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  mobius_arg->ls = GJP.SnodeSites();
  mobius_arg->mobius_kappa_b = 
    1.0 / ( 2 * (GJP.Mobius_b() *(4 - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
  mobius_arg->mobius_kappa_c = 
    1.0 / ( 2 * (GJP.Mobius_c() *(4 - GJP.DwfHeight()) - GJP.DwfA5Inv()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  mobius_mdagm(out, 
	       gauge_field, 
	       in, 
	       dot_prd,
	       mass,
	       (Dwf *) mobius_lib_arg);
}

//------------------------------------------------------------------
// MatPcDagMatPc :
// MatPcDagMatPc is the fermion matrix that appears in the HMC 
// evolution. It is a Hermitian matrix.
// The in, out fields are defined on the checkerboard lattice.
// If dot_prd is not 0 then the dot product (on node)
// <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.
//
// When dirac_arg->eigen_shift is non-zero, it shift the spectrum of matrix:
//    MatPcDagMatPc = H^2  ->  (H-shift)(H-shift)
// where H = Gamma_5 MatPc
//
// For other fermions, one could also implement similar shifts.
// For wilson, H = gamma_5 MatPC .
//------------------------------------------------------------------
void DiracOpMobius::MatPcDagMatPcShift(Vector *out, 
				       Vector *in, 
				       Float *dot_prd){

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------

  const Float shift = dirac_arg->eigen_shift;
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  mobius_arg->ls = GJP.SnodeSites();
  mobius_arg->mobius_kappa_b = 
    1.0 / ( 2 * (GJP.Mobius_b() *(4 - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
  mobius_arg->mobius_kappa_c = 
    1.0 / ( 2 * (GJP.Mobius_c() *(4 - GJP.DwfHeight()) - GJP.DwfA5Inv()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  
  // we still check if shift is really needed
  if( shift == 0.0 )  {
    mobius_mdagm(out, 
		 gauge_field, 
		 in, 
		 dot_prd,
		 mass,
		 (Dwf *) mobius_lib_arg);
  } else {
    
    //mobius_arg->eigen_shift = dirac_arg->eigen_shift;
    
    mobius_mdagm_shift(out, 
		       gauge_field, 
		       in, 
		       dot_prd,
		       mass,
		       (Dwf *) mobius_lib_arg,
		       shift);
  }
  
}


//------------------------------------------------------------------
// Dslash(Vector *out, Vector *in, ChkbType cb, DagType dag) :
// Dslash is the derivative part of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
// cb is not used.
//------------------------------------------------------------------
void DiracOpMobius::Dslash(Vector *out, 
			   Vector *in, 
			   ChkbType cb, 
			   DagType dag) {
  
  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  mobius_arg->ls = GJP.SnodeSites();
  mobius_arg->mobius_kappa_b = 
    1.0 / ( 2 * (GJP.Mobius_b() *(4 - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
  mobius_arg->mobius_kappa_c = 
    1.0 / ( 2 * (GJP.Mobius_c() *(4 - GJP.DwfHeight()) - GJP.DwfA5Inv()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  mobius_dslash(out, 
		gauge_field, 
		in, 
		mass,
		cb,
		dag,
		(Dwf *) mobius_lib_arg);
}

//------------------------------------------------------------------
/*!
  The preconditioned matrix connects sites of even parity.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/
//------------------------------------------------------------------
void DiracOpMobius::MatPc(Vector *out, Vector *in) {  

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  mobius_arg->ls = GJP.SnodeSites();
  mobius_arg->mobius_kappa_b = 
    1.0 / ( 2 * (GJP.Mobius_b() *(4 - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
  mobius_arg->mobius_kappa_c = 
    1.0 / ( 2 * (GJP.Mobius_c() *(4 - GJP.DwfHeight()) - GJP.DwfA5Inv()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  mobius_m(out, 
	   gauge_field, 
	   in, 
	   mass,
	   (Dwf *) mobius_lib_arg);
}

//------------------------------------------------------------------
/*!
  The preconditioned matrix connects sites of even parity.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/
//------------------------------------------------------------------
void DiracOpMobius::MatPcDag(Vector *out, Vector *in) {

  char *fname = "MatPcDag(*V,*V)";
  VRB.Func(fname,cname);
  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  mobius_arg->ls = GJP.SnodeSites();
  mobius_arg->mobius_kappa_b = 
    1.0 / ( 2 * (GJP.Mobius_b() *(4 - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
  mobius_arg->mobius_kappa_c = 
    1.0 / ( 2 * (GJP.Mobius_c() *(4 - GJP.DwfHeight()) - GJP.DwfA5Inv()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  IFloat *tmp = (IFloat *)in;
  mobius_mdag(out, 
	      gauge_field, 
	      in, 
	      mass,
	      (Dwf *) mobius_lib_arg);
  tmp = (IFloat *)out;
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
int DiracOpMobius::MatInv(Vector *out, 
			  Vector *in, 
			  Float *true_res,
			  PreserveType prs_in) {
  char *fname = "MatInv(V*,V*,F*)";
  VRB.Func(cname,fname);
  
  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
  //printf("KAPPA_B %g\n",((Dwf*)mobius_lib_arg)->mobius_kappa_b); exit(0);
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  mobius_arg->ls = GJP.SnodeSites();
  mobius_arg->mobius_kappa_b = 
    1.0 / ( 2 * (GJP.Mobius_b() *(4 - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
  mobius_arg->mobius_kappa_c = 
    1.0 / ( 2 * (GJP.Mobius_c() *(4 - GJP.DwfHeight()) - GJP.DwfA5Inv()) );
  Float minus_kappa_b = -mobius_arg->mobius_kappa_b;
  Float kappa_b = - minus_kappa_b;
  Float norm;

  //printf("KAPPA_B %g\n",kappa_b); exit(0);

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  Vector *temp2;
  Vector *temp3;
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  Vector *temp = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp == 0) ERR.Pointer(cname, fname, "temp");
  VRB.Smalloc(cname,fname, "temp", temp, temp_size * sizeof(Float));

  temp2 = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp2 == 0) ERR.Pointer(cname, fname, "temp2");
  VRB.Smalloc(cname,fname, "temp2", temp2, temp_size * sizeof(Float));

  // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (IFloat *) in + temp_size );

  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (IFloat *) out + temp_size );

  // prepare source
  // mult by Dminus to compare with Hantao
#if 1
  temp3 = (Vector *) smalloc(2*temp_size * sizeof(Float));
  if (temp3 == 0) ERR.Pointer(cname, fname, "temp3");
  VRB.Smalloc(cname,fname, "temp3", temp3, temp_size * sizeof(Float)); 
  Dminus(temp3,in);
  moveFloat((IFloat *)in, (IFloat *)temp3, 2*temp_size);
  VRB.Sfree(cname, fname, "temp3", temp3);
  sfree(temp3);
#endif

  mobius_m5inv(temp, even_in, mass, DAG_NO, mobius_arg);  
  mobius_dslash_4(temp2, gauge_field, temp, CHKB_ODD, DAG_NO, mobius_arg, mass);
  fTimesV1PlusV2((IFloat *)temp, kappa_b, (IFloat *)temp2,
		 (IFloat *)in, temp_size);

  // save source
  if(prs_in == PRESERVE_YES){
    moveMem((IFloat *)temp2, (IFloat *)in, 
	    temp_size * sizeof(IFloat) / sizeof(char));
  }

  int iter;
  switch (dirac_arg->Inverter) {
  case CG:
    MatPcDag(in, temp);
    iter = InvCg(out,in,true_res);
    break;
  case BICGSTAB:
    iter = BiCGstab(out,temp,0.0,dirac_arg->bicgstab_n,true_res);
  case LOWMODEAPPROX :
    MatPcDag(in, temp);
    iter = InvLowModeApprox(out,in, dirac_arg->fname_eigen, dirac_arg->neig, true_res );
    break;
  case CG_LOWMODE_DEFL : 
    MatPcDag(in, temp);
    InvLowModeApprox(out,in, dirac_arg->fname_eigen, dirac_arg->neig, true_res );   
    iter = InvCg(out,in,true_res);
    break;
  default:
    ERR.General(cname,fname,"InverterType %d not implemented\n",
		dirac_arg->Inverter);
  }


  // check solution
  //norm = out->NormSqGlbSum(temp_size);
  //printf("Norm out %.14e\n",norm);
  //norm = in->NormSqGlbSum(temp_size);
  //printf("Norm in %.14e\n",norm);
  //MatPcDagMatPc(temp,out);  
  //norm = temp->NormSqGlbSum(temp_size);
  //printf("Norm MatPcDagMatPc*out %.14e\n",norm);
  //exit(0);

  // restore source
  if(prs_in == PRESERVE_YES){
    moveMem((IFloat *)in, (IFloat *)temp2, 
	    temp_size * sizeof(IFloat) / sizeof(char));
  }

  // TIZB check below carefully !

  // construct even site solution
  // psi_e = M5inv . in_e - kappa*M5inv*D_Weo*psi_o

  //TIZB mobius_dslash_4(temp, gauge_field, out, CHKB_ODD, DAG_NO, mobius_arg);
  mobius_dslash_4(temp, gauge_field, out, CHKB_EVEN, DAG_NO, mobius_arg, mass);
  mobius_m5inv(even_out, temp, mass, DAG_NO, mobius_arg);
  mobius_m5inv(temp, even_in, mass, DAG_NO, mobius_arg);
  fTimesV1PlusV2((IFloat *)even_out, kappa_b, (IFloat *)even_out,
		 (IFloat *)temp, temp_size);
  
  VRB.Sfree(cname, fname, "temp2", temp2);
  sfree(temp2);
  VRB.Sfree(cname, fname, "temp", temp);
  sfree(temp);

 
  return iter;
}


//------------------------------------------------------------------
// Overloaded function is same as original 
// but true_res=0.
//------------------------------------------------------------------
int DiracOpMobius::MatInv(Vector *out, Vector *in, PreserveType prs_in)
{ return MatInv(out, in, 0, prs_in); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in and out = f_out.
//------------------------------------------------------------------
int DiracOpMobius::MatInv(Float *true_res, PreserveType prs_in)
{ return MatInv(f_out, f_in, true_res, prs_in); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in, out = f_out, true_res=0.
//------------------------------------------------------------------
int DiracOpMobius::MatInv(PreserveType prs_in)
{ return MatInv(f_out, f_in, 0, prs_in); }


//------------------------------------------------------------------
// Mat(Vector *out, Vector *in) :
// Mat is the unpreconditioned fermion matrix.  
// Mat works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpMobius::Mat(Vector *out, Vector *in) {  
  char *fname = "Mat(V*,V*)";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  mobius_arg->ls = GJP.SnodeSites();
  mobius_arg->mobius_kappa_b = 
    1.0 / ( 2 * (GJP.Mobius_b() *(4 - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
  mobius_arg->mobius_kappa_c = 
    1.0 / ( 2 * (GJP.Mobius_c() *(4 - GJP.DwfHeight()) - GJP.DwfA5Inv()) );
  Float kappa = mobius_arg->mobius_kappa_b;
  Float minus_kappa = -kappa;
  Float kappa_ratio = mobius_arg->mobius_kappa_b/ mobius_arg->mobius_kappa_c;

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (IFloat *) in + temp_size );
  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (IFloat *) out + temp_size );
  // temp
  Vector *frm_tmp2 = (Vector *) mobius_arg->frm_tmp2;

  //odd part
  //mobius_dslash_4(out, gauge_field, even_in, CHKB_EVEN, DAG_NO, mobius_arg, mass);
  mobius_dslash_4(out, gauge_field, even_in, CHKB_ODD, DAG_NO, mobius_arg, mass);
  out->VecTimesEquFloat(minus_kappa, temp_size); 
  // intialize to zero since using the "plus-equal version"
  for(int i=0;i<temp_size;i++){
    *((IFloat*)frm_tmp2+i)=0.0;
  }
  mobius_dslash_5_plus(frm_tmp2, in, mass, 0, mobius_arg);
  fTimesV1PlusV2((IFloat*)frm_tmp2, kappa_ratio, (IFloat*)frm_tmp2, (IFloat *)in, temp_size);
  out->VecAddEquVec(frm_tmp2, temp_size); 

  //even part
  mobius_dslash_4(even_out, gauge_field, in, CHKB_EVEN, DAG_NO, mobius_arg, mass);
  even_out->VecTimesEquFloat(minus_kappa, temp_size); 
  // intialize to zero since using the "plus-equal version"
  for(int i=0;i<temp_size;i++){
    *((IFloat*)frm_tmp2+i)=0.0;
  }
  mobius_dslash_5_plus(frm_tmp2, even_in, mass, 0, mobius_arg);
  fTimesV1PlusV2((IFloat*)frm_tmp2, kappa_ratio, (IFloat*)frm_tmp2, (IFloat *)even_in, temp_size);
  even_out->VecAddEquVec(frm_tmp2, temp_size);

}


void DiracOpMobius::Dminus(Vector *out, Vector *in) {  
  char *fname = "Dminus(V*,V*)";
  VRB.Func(cname,fname);

  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  Float kappa_c_inv_div2 = 0.5*( 2 * (GJP.Mobius_c() *(4 - GJP.DwfHeight()) - GJP.DwfA5Inv()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  // points to the odd part of fermion source 
  Vector *odd_in = (Vector *) ( (IFloat *) in + temp_size );
  // points to the odd part of fermion solution
  Vector *odd_out = (Vector *) ( (IFloat *) out + temp_size );

  mobius_dminus(out, gauge_field, odd_in, CHKB_ODD, DAG_NO, mobius_arg);
  mobius_dminus(odd_out, gauge_field, in, CHKB_EVEN, DAG_NO, mobius_arg);
  // out = (c*D_W-1)*in
  fTimesV1PlusV2((IFloat*)out, kappa_c_inv_div2, (IFloat*)in, (IFloat *)out, 2*temp_size);

}


//------------------------------------------------------------------
// MatDag(Vector *out, Vector *in) :
// MatDag is the unpreconditioned fermion matrix.  
// MatDag works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpMobius::MatDag(Vector *out, Vector *in) {
  char *fname = "MatDag(V*,V*)";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.MobiusA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  mobius_arg->ls = GJP.SnodeSites();
  mobius_arg->mobius_kappa_b = 
    1.0 / ( 2 * (GJP.Mobius_b() *(4 - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
  mobius_arg->mobius_kappa_c = 
    1.0 / ( 2 * (GJP.Mobius_c() *(4 - GJP.DwfHeight()) - GJP.DwfA5Inv()) );
  Float kappa = mobius_arg->mobius_kappa_b;
  Float minus_kappa = -kappa;
  Float kappa_ratio = mobius_arg->mobius_kappa_b/mobius_arg->mobius_kappa_c;

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (IFloat *) in + temp_size );
  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (IFloat *) out + temp_size );
  // temp
  Vector *frm_tmp2 = (Vector *) mobius_arg->frm_tmp2;

  //odd part
  mobius_dslash_4(out, gauge_field, even_in, CHKB_ODD, DAG_YES, mobius_arg, mass);
  //mobius_dslash_4(out, gauge_field, even_in, CHKB_EVEN, DAG_YES, mobius_arg, mass);
  out->VecTimesEquFloat(kappa, temp_size); 
  // intialize to zero since using the "plus-equal version"
  for(int i=0;i<temp_size;i++){
    *((IFloat*)frm_tmp2+i)=0.0;
  }
  mobius_dslash_5_plus(frm_tmp2, in, mass, DAG_YES, mobius_arg);
  fTimesV1PlusV2((IFloat*)frm_tmp2, kappa_ratio, (IFloat*)frm_tmp2, (IFloat *)in, temp_size);
  out->VecAddEquVec(frm_tmp2, temp_size); 

  //even part
  //mobius_dslash_4(even_out, gauge_field, in, CHKB_ODD, DAG_YES, mobius_arg, mass);
  mobius_dslash_4(even_out, gauge_field, in, CHKB_EVEN, DAG_YES, mobius_arg, mass);
  even_out->VecTimesEquFloat(kappa, temp_size); 
  // intialize to zero since using the "plus-equal version"
  for(int i=0;i<temp_size;i++){
    *((IFloat*)frm_tmp2+i)=0.0;
  }
  mobius_dslash_5_plus(frm_tmp2, even_in, mass, DAG_YES, mobius_arg);
  fTimesV1PlusV2((IFloat*)frm_tmp2, kappa_ratio, (IFloat*)frm_tmp2, (IFloat *)even_in, temp_size);
  even_out->VecAddEquVec(frm_tmp2, temp_size); 

}


//------------------------------------------------------------------
// MatHerm(Vector *out, Vector *in) :
// MatHerm is gamma5*R*Mat.
// MatHerm works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpMobius::MatHerm(Vector *out, Vector *in) {
  char *fname = "MatHerm(V*,V*)";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  mobius_arg->ls = GJP.SnodeSites();
  mobius_arg->mobius_kappa_b = 
    1.0 / ( 2 * (GJP.Mobius_b() *(4 - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
  mobius_arg->mobius_kappa_c = 
    1.0 / ( 2 * (GJP.Mobius_c() *(4 - GJP.DwfHeight()) - GJP.DwfA5Inv()) );

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
/*!
  \pre This method is to be used when the instance of this object has been
  created with \a f_field_out and \a f_field_in pointers to spin-colour
  vectors defined over the whole lattice. 

  \param chi A spin-colour vector defined on odd parity lattice sites.

  \post The vector \a f_field_out is \f$ (1+D)\chi \f$

  and the vector \a f_field_in is \f$ (D^\dagger-\kappa^2 M)\chi \f$

  where \e M is the odd-even preconditioned fermion matrix connecting odd to
  odd parity sites and \e D is the hopping term connecting odd to
  even parity sites. Recall that \a chi is defined on odd sites only.
  The new vectors are in odd-even order.
*/
//------------------------------------------------------------------

void DiracOpMobius::CalcHmdForceVecs(Vector *chi)
{
  char *fname = "CalcHmdForceVecs(V*)" ;
  VRB.Func(cname,fname) ;

  if (f_out == 0)
    ERR.Pointer(cname, fname, "f_out") ;

  if (f_in == 0)
    ERR.Pointer(cname, fname, "f_in") ;

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // and DiracOpMobius constructors but is done again in case the
  // user has modified the GJP.MobiusA5Inv(), GJP.MobiusHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpMobius object.
  //----------------------------------------------------------------
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  mobius_arg->ls = GJP.SnodeSites();
  mobius_arg->mobius_kappa_b = 
    1.0 / ( 2 * (GJP.Mobius_b() *(4 - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
  mobius_arg->mobius_kappa_c = 
    1.0 / ( 2 * (GJP.Mobius_c() *(4 - GJP.DwfHeight()) - GJP.DwfA5Inv()) );

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

//  fprintf(stderr,"psi=%p chi=%p rho=%p sigma=%p\n",psi,chi,rho,sigma);
  MatPc(psi,chi) ;
//  fprintf(stderr,"MatPc\n");

  {
    Float kappa = ((Dwf *)mobius_lib_arg)->mobius_kappa_b ;
    psi->VecTimesEquFloat(-kappa*kappa,f_size_cb) ;
  }

  rho = (Vector *)((Float *)f_out + f_size_cb) ;

  Dslash(rho, chi, CHKB_ODD, DAG_NO) ;
//  fprintf(stderr,"Dslash\n");

  sigma = (Vector *)((Float *)f_in + f_size_cb) ;

  Dslash(sigma, psi, CHKB_ODD, DAG_YES) ;
//  fprintf(stderr,"Dslash\n");

  return ;
}

//------------------------------------------------------------------
// DiracOpGlbSum(Float *): 
// The global sum used by InvCg. If s_nodes = 1
// it is the usual global sum. If s_nodes > 1 it
// is the 5-dimensional globals sum glb_sum_five.
//------------------------------------------------------------------
void DiracOpMobius::DiracOpGlbSum(Float *float_p) {
//  if(GJP.Snodes() == 1) {
//    glb_sum(float_p);
//  }
//  else {
    glb_sum_five(float_p);
//  }
}



#ifndef USE_BLAS
#define MOVE_FLOAT( pa, pb, n )  moveFloat(pa, pb, n) 
#define VEC_TIMESEQU_FLOAT(py, fact, n ) vecTimesEquFloat( py, fact, n)
#define AXPY(n, fact, px, py)  fTimesV1PlusV2(py, fact, px, py, n)
#else
#define MOVE_FLOAT( pa, pb, n )  cblas_dcopy(n, pb, 1, pa, 1)
#define VEC_TIMESEQU_FLOAT(py, fact, n ) cblas_dscal( n,  fact, py,1 )
#define AXPY(n, fact, px, py)  cblas_daxpy(n, fact, px,1,py,1)
#endif


//------------------------------------------------------------------
// RitzMat(Vector *out, Vector *in) :
// RitzMat is the base operator used in in Ritz.
// RitzMat works on the full or half lattice.
// The in, out fields are defined on the full or half lattice.
//------------------------------------------------------------------
void DiracOpMobius::RitzMat(Vector *out, Vector *in) {
  char *fname = "RitzMat(V*,V*)";
  VRB.Func(cname,fname);

  //printf("single ritzmat %d\n",dirac_arg->RitzMatOper);
  switch(dirac_arg->RitzMatOper)
    {
      // Now always call MatPcDagMatPcShift even for MATPCDAG_MATPC case
      //
      //case MATPCDAG_MATPC:
      //MatPcDagMatPc(out, in);
      //break;
    case MATPCDAG_MATPC:
    case MATPCDAG_MATPC_SHIFT:
      MatPcDagMatPcShift(out, in);
      break;
    case MATPC_HERM :
      MatPcHerm(out, in);
      break;

    case MAT_HERM:
    case MATDAG_MAT:
      MatDagMat(out, in);
      break;
    case NEG_MATPCDAG_MATPC:
      MatPcDagMatPc(out, in);
      out -> VecNegative(out, RitzLatSize());
      break;    
      
    case NEG_MATDAG_MAT:
      MatDagMat(out, in);
      out->VecNegative(out, RitzLatSize());
      break;
      
    default:
      ERR.General(cname,fname,"RitzMatOper %d not implemented\n",
		  dirac_arg->RitzMatOper);
    }

#if 0
  //debug
  const int size = RitzLatSize(); // this is number of Float, f_size .
  Complex deb = out->CompDotProductGlbSum(in, size);
  Float d_deb = in->NormSqGlbSum(size);
  printf("single Ritz %e %e %e\n", deb.real()/d_deb, deb.imag()/d_deb,d_deb);
  //debug
#endif
}



// PolynomialAccerelation
//
//  Q = [ -2 RitzMat + (alpha + beta) ] / [ alpha - beta ]
//
//  Output:  out =  T_n(Q) in
//
//  T_0 = 1,    T_1 = Q
//   T_{n+1}(Q) =  2 Q T_n(Q)  - T_{n-1}(Q)
//  
// Calling virtual RitzMat(V*,V*)
//
//   alpha = param[0]^2
//   beta = ( param[1] + fabs(eigend_shift) )^2
//
void DiracOpMobius::RitzMat(Vector *out, Vector *in,
			 MatrixPolynomialArg* cheby_arg) {
  char *fname = "RitzMat(V*,V*,MatrixPolyArg*)";
  VRB.Func(cname,fname);

  double time_start=dclock();
    
  //debug const
  const Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  const Float shift = dirac_arg->eigen_shift;

  const int Npol = cheby_arg-> Npol;
  const int size = RitzLatSize(); // this is number of Float, f_size .
  //
  // Q = 2 / (alpha-beta)  RitzMat  -   (alpha+beta)/(alpha-beta)
  // 2 Q =   c1  (  c0 Ddag D  -   1 )
  //  c1 = 2 (alpha+beta)/(alpha-beta),   c0 =  2 / (alpha+beta)
  //
  const Float alpha = pow( cheby_arg-> params.params_val[0], 2);
  const Float beta  = pow( cheby_arg-> params.params_val[1] + fabs(shift), 2);
  //printf("alpha=%e beta=%e\n", alpha,beta);
  
  const Float c1 =   2.0*(alpha+beta)/(alpha-beta);
  const Float c0 =   2.0/(alpha+beta);

  Vector *tmp  = (Vector*)cheby_arg->tmp1;
  Vector *tmp2 = (Vector*)cheby_arg->tmp2;
  

  //  tmp2 =  T_0 v = v = in
  //tmp2 -> CopyVec(in, size);
  MOVE_FLOAT( (Float*)tmp2, (Float*)in, size );
  //  tmp =  T_1 v = Q v = Q in
  //  QV = 0.5* (2Q)V = 0.5 c1 ( c0 Ddag D - 1)

  RitzMat(tmp, in);
  
#if 0
  tmp->VecTimesEquFloat(c0, size);
  tmp->VecMinusEquVec(in,size);
  tmp->VecTimesEquFloat(0.5*c1, size);
  //  tmp =  0.5 c1 ( c0 DdagD in - in )
#else
  VEC_TIMESEQU_FLOAT((Float*)tmp,0.5*c1*c0, size);
  AXPY( size, -0.5*c1, (Float*)in, (Float*)tmp );
#endif
  
  // debug
  //out->CopyVec(tmp,size);
  //printf("cheby %f %f\n", alpha,beta);
  
  // loop over
  for(int i=2; i<=Npol; ++i){
#if 0
    // out = 2 Q tmp
    RitzMat(out, tmp);

    out->VecTimesEquFloat(c0, size);
    out->VecMinusEquVec(tmp,size);
    out->VecTimesEquFloat(c1, size);

    // out = out - tmp2
    out->VecMinusEquVec(tmp2, size);
#else
    // out = c1 (  c0 DagD tmp - tmp) -tmp2

    RitzMat(out, tmp);

    VEC_TIMESEQU_FLOAT((Float*)out, c1*c0, size);
    AXPY( size , - c1, (Float*)tmp, (Float*)out);
    AXPY( size,  -1.0, (Float*)tmp2, (Float*)out);
#endif

    
    if( i!=Npol) {
#if 0
      // tmp2 = tmp
      tmp2->CopyVec(tmp, size);
      // tmp = out
      tmp->CopyVec(out, size);
#else
      // tmp2 = tmp
      Vector* swap_tmp2 = tmp2;
      tmp2 = tmp;
      tmp = swap_tmp2;
      // tmp = out
      tmp->CopyVec(out, size);
#endif
    }
  }
 
  
  if(!UniqueID()) printf("mpoly tot = %e\n", dclock() - time_start);
  
#if 0
  //debug
  Complex deb = out->CompDotProductGlbSum(in, size);
  printf("debug %e %e %e\n", deb.real(), deb.imag(),in->NormSqGlbSum(size));
  //debug
#endif
  
}


//!! N.B. This overwrites contents of  mobius_arg->frm_tmp2
void DiracOpMobius::MatPcHerm(Vector *out, Vector *in) {
  char *fname = "MatPcHerm(V*,V*)";
  VRB.Func(cname,fname);
  Dwf *mobius_arg = (Dwf *) mobius_lib_arg;
  Vector* vtmp = (Vector*)(mobius_arg->frm_tmp1);
  
  MatPc(vtmp,in);
  ReflectAndMultGamma5( out, vtmp,  mobius_arg->vol_4d/2, mobius_arg->ls);
  
}

// specific to dwf 
void ReflectAndMultGamma5( Vector *out, const Vector *in,  int nodevol, int ls)
{
  char *fname = "MultGamma5(V*,V*,i)";
  VRB.Func("",fname);
  for(int s=0; s< ls; ++s) { 
    IFloat *p = (IFloat *)out + 24*nodevol*s;
    IFloat *q = (IFloat *)in + 24*nodevol*(ls-1-s);
    for(int n = 0; n < nodevol; ++n)
      {
	int i;
	for(i = 0; i < 12; ++i)
	  *p++ = *q++;
	
	for(i = 0; i < 12; ++i)
	  *p++ = - *q++;
      }
  }

}

void HermicianDWF_ee( Vector* vtmp, Vector* evec, Float mass, Lattice* lattice, Vector* Apsi )
{
	CgArg cg_arg;
	cg_arg.mass = mass;
	cg_arg.RitzMatOper = MATPC_HERM; // could be MATPCDAG_MATPC;
	DiracOpDwf dop( *lattice, 0, 0, &cg_arg, CNV_FRM_NO );

	dop. MatPc(Apsi, evec);
	ReflectAndMultGamma5( vtmp, Apsi,  
			      GJP.VolNodeSites()/2, GJP.SnodeSites() );
}

CPS_END_NAMESPACE
