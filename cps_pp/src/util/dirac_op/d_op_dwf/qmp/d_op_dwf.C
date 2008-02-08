#include <config.h>
#include <stdio.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOpDwf class methods.

  $Id: d_op_dwf.C,v 1.3 2008-02-08 18:35:07 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008-02-08 18:35:07 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qmp/d_op_dwf.C,v 1.3 2008-02-08 18:35:07 chulwoo Exp $
//  $Id: d_op_dwf.C,v 1.3 2008-02-08 18:35:07 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: d_op_dwf.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qmp/d_op_dwf.C,v $
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
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/wilson.h>
#include <util/time_cps.h>
#include <util/dwf.h>
#include <mem/p2v.h>
#include <comms/glb.h>
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
}


//------------------------------------------------------------------
/*!
  If the storage order of any fields was changed by the constructor
  then they are changed to the canonical order by the destructor.
*/
//------------------------------------------------------------------
DiracOpDwf::~DiracOpDwf() {

}


//------------------------------------------------------------------
// DiracArg(CgArg *arg):
// It sets the dirac_arg pointer to arg and initializes
// mass.
//------------------------------------------------------------------
void DiracOpDwf::DiracArg(CgArg *arg){
   ERR.NotImplemented(cname, "DiracArg");
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
   ERR.NotImplemented(cname, "MatPcDagMatPc");
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
   ERR.NotImplemented(cname, "Dslash");
}

//------------------------------------------------------------------
/*!
  The preconditioned matrix connects sites of odd parity.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/
//------------------------------------------------------------------
void DiracOpDwf::MatPc(Vector *out, Vector *in) {  
   ERR.NotImplemented(cname, "MatPc");
}

//------------------------------------------------------------------
/*!
  The preconditioned matrix connects sites of odd parity.
  \param out The resulting vector.
  \param in The vector to be multiplied.
*/
//------------------------------------------------------------------
void DiracOpDwf::MatPcDag(Vector *out, Vector *in) {

   ERR.NotImplemented(cname, "MatPcDag");
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
   ERR.NotImplemented(cname, "MatInv");
}


//------------------------------------------------------------------
// Overloaded function is same as original 
// but true_res=0.
//------------------------------------------------------------------
int DiracOpDwf::MatInv(Vector *out, Vector *in, PreserveType prs_in)
{    ERR.NotImplemented(cname, "MatInv");}


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in and out = f_out.
//------------------------------------------------------------------
int DiracOpDwf::MatInv(Float *true_res, PreserveType prs_in)
{    ERR.NotImplemented(cname, "MatInv");}


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in, out = f_out, true_res=0.
//------------------------------------------------------------------
int DiracOpDwf::MatInv(PreserveType prs_in)
{    ERR.NotImplemented(cname, "MatInv");}


//------------------------------------------------------------------
// Mat(Vector *out, Vector *in) :
// Mat is the unpreconditioned fermion matrix.  
// Mat works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpDwf::Mat(Vector *out, Vector *in) {  
   ERR.NotImplemented(cname, "Mat");
}


//------------------------------------------------------------------
// MatDag(Vector *out, Vector *in) :
// MatDag is the unpreconditioned fermion matrix.  
// MatDag works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpDwf::MatDag(Vector *out, Vector *in) {
   ERR.NotImplemented(cname, "MatDag");
}


//------------------------------------------------------------------
// MatHerm(Vector *out, Vector *in) :
// MatHerm is gamma5*R*Mat.
// MatHerm works on the full lattice.
// The in, out fields are defined on the full lattice.
//------------------------------------------------------------------
void DiracOpDwf::MatHerm(Vector *out, Vector *in) {
   ERR.NotImplemented(cname, "MatHerm");
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

void DiracOpDwf::CalcHmdForceVecs(Vector *chi)
{
  ERR.NotImplemented(cname, "CalcHmdForceVecs");
}

//------------------------------------------------------------------
// DiracOpGlbSum(Float *): 
// The global sum used by InvCg. If s_nodes = 1
// it is the usual global sum. If s_nodes > 1 it
// is the 5-dimensional globals sum glb_sum_five.
//------------------------------------------------------------------
void DiracOpDwf::DiracOpGlbSum(Float *float_p) {
  ERR.NotImplemented(cname, "DiracOpGlbSum");
}


CPS_END_NAMESPACE
