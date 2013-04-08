#include<stdio.h>
#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOpStag class methods.

  $Id: d_op_stag.C,v 1.12 2013-04-08 20:50:00 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//  $Author: chulwoo $
//  $Date: 2013-04-08 20:50:00 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_stag/noarch/d_op_stag.C,v 1.12 2013-04-08 20:50:00 chulwoo Exp $
//  $Id: d_op_stag.C,v 1.12 2013-04-08 20:50:00 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: d_op_stag.C,v $
//  $Revision: 1.12 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_stag/noarch/d_op_stag.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// d_op_stag.C
//
// DiracOpStag is derived from the DiracOpStagTypes class.
// DiracOpStag is the front end for a library that contains
// all Dirac operators associated with Staggered fermions.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/vector.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/stag.h>
#include <comms/cbuf.h>
#include <comms/glb.h>
#include <math.h>
CPS_START_NAMESPACE



const unsigned CBUF_MODE1 = 0xcb911548;
const unsigned CBUF_MODE2 = 0xcca52112;
const unsigned CBUF_MODE3 = 0xc98c6106;
const unsigned CBUF_MODE4 = 0xcca52112;


//------------------------------------------------------------------
/*!
  Only one instance of this class is allowed to be in existence at
  any time.
  \param latt The lattice on which this Dirac operator is defined
  \param f_field_out A (pointer to) a spin-colour field (optionally). 
  \param f_field_in A (pointer to) a spin-colour field (optionally).
  \param arg Parameters for the solver.
  \param convert Whether the lattice fields should be converted to
  to a new storage order appropriate for the type of fermion action.
  If this is ::CNV_FRM_NO, then just the gauge field is converted.
  If this is ::CNV_FRM_YES, then the fields \a f_field_out and \a f_field_in
  are also converted: This assumes they are initially in the same order as
  the gauge field.
 */
//------------------------------------------------------------------
DiracOpStag::DiracOpStag(Lattice & latt,
			 Vector *f_field_out,
			 Vector *f_field_in,
			 CgArg *arg,
			 CnvFrmType convert) :
			 DiracOpStagTypes(latt, 
					  f_field_out,
					  f_field_in, 
					  arg,
					  convert)
{
  cname = "DiracOpStag";
  char *fname = "DiracOpStag(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Do the necessary conversions
  //----------------------------------------------------------------
  if(cnv_frm == CNV_FRM_YES)
    lat.Convert(STAG, f_out, f_in);
  else if(cnv_frm == CNV_FRM_NO)
    lat.Convert(STAG);

  //----------------------------------------------------------------
  // Set the node checkerboard size of the fermion field
  //----------------------------------------------------------------
  f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  //----------------------------------------------------------------
  // Allocate memory for the temporary fermion vector frm_tmp.
  //----------------------------------------------------------------
  frm_tmp = (Vector *) smalloc(f_size_cb * sizeof(Float));
  if(frm_tmp == 0)
    ERR.Pointer(cname,fname, "frm_tmp");
  VRB.Smalloc(cname,fname, "frm_tmp", 
	      frm_tmp, f_size_cb * sizeof(Float));

  //----------------------------------------------------------------
  // Initialize parameters
  //----------------------------------------------------------------
  DiracArg(dirac_arg);

}


//------------------------------------------------------------------
/*!
  If the storage order of any fields was changed by the constructor
  then they are changed to the canonical order by the destructor.
*/
//------------------------------------------------------------------
DiracOpStag::~DiracOpStag() {
  char *fname = "~DiracOpStag()";
  VRB.Func(cname,fname);

  if(cnv_frm == CNV_FRM_YES)
    lat.Convert(CANONICAL, f_out, f_in);
  else if(cnv_frm == CNV_FRM_NO)
    lat.Convert(CANONICAL);

  //----------------------------------------------------------------
  // Free memory
  //----------------------------------------------------------------
  VRB.Sfree(cname,fname, "frm_tmp", frm_tmp);
  sfree(frm_tmp);
}


//------------------------------------------------------------------
// DiracArg(CgArg *arg):
// It sets the dirac_arg pointer to arg and initializes
// mass_sq = 4 * mass^2.
//------------------------------------------------------------------
  void DiracOpStag::DiracArg(CgArg *arg){
    dirac_arg = arg;

    // Added for anisotropic lattices
    //------------------------------------------------------------------
    mass_rs = dirac_arg -> mass * GJP.XiBare()/GJP.XiV();
    mass_sq = 4 * mass_rs * mass_rs;
    // End modification

  }


//------------------------------------------------------------------
// MatPcDagMatPc :
// MatPcDagMatPc is the fermion matrix that appears in the HMC 
// evolution. It is a Hermitian matrix where M is
// the even/odd preconditioned Dirac Operator matrix.        
// MatPcDagMatPc connects only even-->even sites.
// The in, out fields are defined on the even checkerboard.
// If dot_prd is not 0 then the dot product (on node)
// <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.
//------------------------------------------------------------------
void DiracOpStag::MatPcDagMatPc(Vector *out, 
				Vector *in, 
				Float *dot_prd){
//  setCbufCntrlReg(1, CBUF_MODE1);
//  setCbufCntrlReg(2, CBUF_MODE2);
//  setCbufCntrlReg(3, CBUF_MODE3);
//  setCbufCntrlReg(4, CBUF_MODE4);
  stag_dirac(frm_tmp, in, 0, 0);
  stag_dirac(out, frm_tmp, 1, 0);
  out->FTimesV1MinusV2(mass_sq, in, out, f_size_cb);

  if( dot_prd !=0 ){
    *dot_prd = dotProduct((IFloat *) in, (IFloat *) out, f_size_cb);
  }
}


//------------------------------------------------------------------
// Dslash(Vector *out, Vector *in, ChkbType cb, DagType dag) :
// Dslash is the derivative part of the fermion matrix. 
// Dslash conects only odd-->even or even-->odd sites.
// The in, out fields are defined on a checkerboard.
// cb refers to the checkerboard of the in field.
//------------------------------------------------------------------
void DiracOpStag::Dslash(Vector *out, 
			 Vector *in, 
			 ChkbType cb, 
			 DagType dag) {
			        
  setCbufCntrlReg(1, CBUF_MODE1);
  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(3, CBUF_MODE3);
  setCbufCntrlReg(4, CBUF_MODE4);

  stag_dirac(out, 
	in, 
	int(cb),
	int(dag));
}


//------------------------------------------------------------------
// Dslash(Vector *out, Vector *in, ChkbType cb, DagType dag) :
// Dslash is the derivative part of the fermion matrix. 
// Dslash conects only odd-->even or even-->odd sites.
// The in, out fields are defined on a checkerboard.
// cb refers to the checkerboard of the in field.
// dir_flag is flag which takes value 0 when all direction contribute to D,
// 1 - when only the special anisotropic direction contributes to D,
// 2 - when all  except the special anisotropic direction.
  /*!
    A vector defined on lattice sites of single parity is multiplied by
    all or some of the derivative part of the fermion matrix; on an
    anisotropic lattice
    \param out The resulting vector.
    \param in The vector to be multiplied.
    \param cb The parity on which the vector \a in is defined.
    \param dag Whether to multiply by the hermitian conjugate matrix:
    Should be set to 1 to multiply by the hermitian conjugate, 0 otherwise.
    \param dir_flag Which parts of the derivative matrix to use:
    When 0, the derivative matrix is the sum over derivatives in all
    directions.
    When 1, the derivative matrix is the  derivative in the anisotropic
    direction.
    When 2, the derivative matrix is the sum over derivatives in all
    directions but the anisotropic direction.
  */
//-------------------------------------------------------------------
void DiracOpStag::Dslash(Vector *out, 
			 Vector *in, 
			 ChkbType cb, 
			 DagType dag,
			 int dir_flag) {

  setCbufCntrlReg(1, CBUF_MODE1);
  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(3, CBUF_MODE3);
  setCbufCntrlReg(4, CBUF_MODE4);

  stag_dirac(out, 
	in, 
	int(cb),
	int(dag),
	dir_flag);
}


//------------------------------------------------------------------
// int MatInv(Vector *out, Vector *in, 
//            Float *true_res, PreserveType prs_in);
// The inverse of the Dirac Operator (D+m)
// using Conjugate gradient.
// Assume: the vector in contains both even and odd src.
//               even part is the 1st part.
// Return: the vector out contains both even and odd solutions.
//       the even solution is the 1st part.
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// prs_in is not used. The source in is always preserved.
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int DiracOpStag::MatInv(Vector *out, 
			Vector *in, 
			Float *true_res,
			PreserveType prs_in) {
  char *fname = "MatInv(V*,V*,F*)";
  VRB.Func(cname,fname);

  setCbufCntrlReg(1, CBUF_MODE1);
  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(3, CBUF_MODE3);
  setCbufCntrlReg(4, CBUF_MODE4);

  IFloat *k_e = (IFloat *)in;
  IFloat *k_o = k_e+f_size_cb;

  Vector *tmp = (Vector *) smalloc(f_size_cb * sizeof(Float));
  if(tmp == 0)
    ERR.Pointer(cname,fname, "tmp");
  VRB.Smalloc(cname,fname, "tmp", 
	      tmp, f_size_cb * sizeof(Float));

  // tmp = (2m - D)k
  stag_dirac(tmp, (Vector *)k_o, 1, 0);
  fTimesV1MinusV2((IFloat *)tmp, 2.*mass_rs, k_e,
  	(IFloat *)tmp, f_size_cb);


  int iter;
  switch (dirac_arg->Inverter) {
  case CG:
    iter = InvCg(out,tmp,true_res);
    break;
  case LOWMODEAPPROX :
    iter = InvLowModeApprox(out,tmp, dirac_arg->fname_eigen, dirac_arg->neig, true_res );
    break;
  case CG_LOWMODE_DEFL :
    InvLowModeApprox(out,tmp, dirac_arg->fname_eigen, dirac_arg->neig, true_res );   
    iter = InvCg(out,tmp,true_res);
    break;
  default:
    ERR.General(cname,fname,"InverterType %d not implemented\n",
                dirac_arg->Inverter);
  }

  // calculate odd solution
  IFloat *x_e = (IFloat *)out;
  IFloat *x_o = x_e+f_size_cb;
  moveMem(x_o, k_o, f_size_cb*sizeof(Float) / sizeof(char));
  stag_dirac(tmp, (Vector *)x_e, 0, 0);
  vecMinusEquVec(x_o, (IFloat *)tmp, f_size_cb);
  vecTimesEquFloat(x_o, 0.5/mass_rs, f_size_cb);

  sfree(tmp);

  return iter;
}


//------------------------------------------------------------------
// Overloaded function is same as original 
// but true_res=0.
//------------------------------------------------------------------
int DiracOpStag::MatInv(Vector *out, Vector *in, PreserveType prs_in)
{ return MatInv(out, in, 0, prs_in); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in and out = f_out.
//------------------------------------------------------------------
int DiracOpStag::MatInv(Float *true_res, PreserveType prs_in)
{ return MatInv(f_out, f_in, true_res, prs_in); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in, out = f_out, true_res=0.
//------------------------------------------------------------------
int DiracOpStag::MatInv(PreserveType prs_in)
{ return MatInv(f_out, f_in, 0, prs_in); }


//------------------------------------------------------------------
// RitzEigMat(Vector *out, Vector *in) :
// RitzEigMat is the base operator used in in RitzEig.
// RitzEigMat works on the full or half lattice.
// The in, out fields are defined on the full or half lattice.
//------------------------------------------------------------------
void DiracOpStag::RitzEigMat(Vector *out, Vector *in) {
  ERR.NotImplemented(cname,"RitzEigMat");
}

CPS_END_NAMESPACE
