#include<config.h>
//#include<qalloc.h>
CPS_START_NAMESPACE
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
#include<util/dirac_op.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/stag.h>
#include<comms/cbuf.h>
#include<comms/glb.h>
#include<comms/scu.h>
//#include<mem/p2v.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
DiracOpStag::DiracOpStag(Lattice & latt,
			 Vector *f_field_out,
			 Vector *f_field_in,
			 CgArg *arg,
			 CnvFrmType cnv_frm_flg) :
			 DiracOpStagTypes(latt, 
					  f_field_out,
					  f_field_in, 
					  arg,
					  cnv_frm_flg)
{
 cname = "DiracOpStag";
}

//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
DiracOpStag::~DiracOpStag() {
}


//------------------------------------------------------------------
// DiracArg(CgArg *arg):
// It sets the dirac_arg pointer to arg and initializes
// mass_sq = 4 * mass^2.
//------------------------------------------------------------------
void DiracOpStag::DiracArg(CgArg *arg){
   ERR.NotImplemented(cname,"DiracArg");
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
  ERR.NotImplemented(cname,"MatPcDagMatPc");
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
  ERR.NotImplemented(cname,"MatInv");
}


//------------------------------------------------------------------
// Dslash(Vector *out, Vector *in, ChkbType cb, DagType dag, int dir_flag) :
// Dslash is the derivative part of the fermion matrix. 
// Dslash conects only odd-->even or even-->odd sites.
// The in, out fields are defined on a checkerboard.
// cb refers to the checkerboard of the in field.
// dir_flag is flag which takes value 0 when all direction contribute to D,
// 1 - when only the special anisotropic direction contributes to D,
// 2 - when all  except the special anisotropic direction. 
//------------------------------------------------------------------
void DiracOpStag::Dslash(Vector *out, 
			 Vector *in, 
			 ChkbType cb, 
			 DagType dag,
			 int dir_flag) {
  ERR.NotImplemented(cname,"MatInv");
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
   ERR.NotImplemented(cname,"MatInv");
}


//------------------------------------------------------------------
// Overloaded function is same as original 
// but true_res=0.
//------------------------------------------------------------------
int DiracOpStag::MatInv(Vector *out, Vector *in, PreserveType prs_in)
{   ERR.NotImplemented(cname,"MatInv"); }


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in and out = f_out.
//------------------------------------------------------------------
int DiracOpStag::MatInv(Float *true_res, PreserveType prs_in)
{   ERR.NotImplemented(cname,"MatInv");}


//------------------------------------------------------------------
// Overloaded function is same as original 
// but in = f_in, out = f_out, true_res=0.
//------------------------------------------------------------------
int DiracOpStag::MatInv(PreserveType prs_in)
{   ERR.NotImplemented(cname,"MatInv"); }

//------------------------------------------------------------------
// RitzMat(Vector *out, Vector *in) :
// RitzMat is the base operator used in in Ritz.
// RitzMat works on the full or half lattice.
// The in, out fields are defined on the full or half lattice.
//------------------------------------------------------------------
void DiracOpStag::RitzMat(Vector *out, Vector *in) {
  ERR.NotImplemented(cname,"RitzMat");
}

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
