#include<stdio.h>
#include<config.h>
#include<math.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of Fasqtad class.

*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_asqtad/f_asqtad.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_asqtad.C
//
// Fasqtad is derived from FstagTypes and is relevant to
// staggered fermions
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/dirac_op.h>
#include <util/asqtad.h>
#include <util/vector.h>
#include <util/gjp.h>
CPS_START_NAMESPACE


#if TARGET == QCDOC
void  set_pt (Fasqtad *lat);
#endif


//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
Fasqtad::Fasqtad()
{
  cname = "Fasqtad";
  char *fname = "Fasqtad()";
  VRB.Func(cname,fname);

  asqtad_dirac_init(GaugeField());
#if TARGET == QCDOC
  set_pt (this);
#endif

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
Fasqtad::~Fasqtad()
{
  char *fname = "~Fasqtad()";
  VRB.Func(cname,fname);
  asqtad_destroy_dirac_buf();
}


//------------------------------------------------------------------
// Returns the type of fermion class.
//------------------------------------------------------------------
FclassType Fasqtad::Fclass() const{
  return F_CLASS_ASQTAD;
}



//------------------------------------------------------------------
// int FmatEvlInv(Vector *f_out, Vector *f_in, 
//                CgArg *cg_arg, 
//                Float *true_res,
//		  CnvFrmType cnv_frm = CNV_FRM_YES):
// It calculates f_out where A * f_out = f_in and
// A is the fermion matrix that appears in the HMC 
// evolution ([Dirac^dag Dirac]). The inversion is done
// with the conjugate gradient. cg_arg is the structure
// that contains all the control parameters, f_in is the
// fermion field source vector, f_out should be set to be
// the initial guess and on return is the solution.
// f_in and f_out are defined on a checkerboard.
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int Fasqtad::FmatEvlInv(Vector *f_out, Vector *f_in, 
		      CgArg *cg_arg, 
		      Float *true_res,
		      CnvFrmType cnv_frm)
{
  int iter;
  char *fname = "FmatEvlInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  DiracOpAsqtad asqtad(*this, f_out, f_in, cg_arg, cnv_frm);
  iter = asqtad.InvCg(true_res);

  asqtad.Dslash(f_tmp, f_out, CHKB_EVEN, DAG_NO);

  // Return the number of iterations
  return iter;
}




//------------------------------------------------------------------
// int FmatEvlMInv(Vector **f_out, Vector *f_in, 
//                Float shift[], int Nshift, 
//                CgArg *cg_arg, Float *true_res,
//		  CnvFrmType cnv_frm = CNV_FRM_YES):
// It calculates f_out where (A + shift)* f_out = f_in and
// A is the fermion matrix that appears in the HMC 
// evolution ([Dirac^dag Dirac]) and shift is a real shift of the 
// fermion matrix, with Nshift such shifts. The inversion is done 
// with the multishift conjugate gradient. cg_arg is the structure
// that contains all the control parameters, f_in is the
// fermion field source vector, f_out is the array of solution 
// vectors, f_in and f_out are defined on a checkerboard.
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int Fasqtad::FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
			 int Nshift, int isz, CgArg *cg_arg,
			 CnvFrmType cnv_frm, MultiShiftSolveType type, 
			 Float *alpha, Vector** f_out_d)
{
  char *fname = "FmatMInv(V**, V*, .....)";
  VRB.Func(cname,fname);

  Float dot = f_in -> NormSqGlbSum(e_vsize);

  Float RsdCG[Nshift];
  for (int s=0; s<Nshift; s++) RsdCG[s] = cg_arg->stop_rsd;

  //Fake the constructor
  DiracOpAsqtad asqtad(*this, f_out[0], f_in, cg_arg, cnv_frm);

  int iter = asqtad.MInvCG(f_out,f_in,dot,shift,Nshift,isz,RsdCG,type,alpha);

  if (type == MULTI && f_out_d != 0)
    for (int s=0; s<Nshift; s++)
      asqtad.Dslash(f_out_d[s],f_out[s],CHKB_EVEN,DAG_NO);

  return iter;

}

//------------------------------------------------------------------
// int FmatInv(Vector *f_out, Vector *f_in, 
//             CgArg *cg_arg, 
//             Float *true_res,
//             CnvFrmType cnv_frm = CNV_FRM_YES,
//             PreserveType prs_f_in = PRESERVE_YES):
// It calculates f_out where A * f_out = f_in and
// A is the fermion matrix (Dirac operator). The inversion
// is done with the conjugate gradient. cg_arg is the 
// structure that contains all the control parameters, f_in 
// is the fermion field source vector, f_out should be set 
// to be the initial guess and on return is the solution.
// f_in and f_out are defined on the whole lattice.
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// cnv_frm is used to specify if f_in should be converted 
// from canonical to fermion order and f_out from fermion 
// to canonical. 
// prs_f_in is used to specify if the source
// f_in should be preserved or not. If not the memory usage
// is less by half the size of a fermion vector.
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int Fasqtad::FmatInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm,
		   PreserveType prs_f_in)
{
  int iter;
  char *fname = "FmatInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  DiracOpAsqtad stag(*this, f_out, f_in, cg_arg, cnv_frm);
  
  iter = stag.MatInv(true_res, prs_f_in);
  
  // Return the number of iterations
  return iter;
}



//------------------------------------------------------------------
// int FeigSolv(Vector **f_eigenv, Float *lambda, int valid_eig[],
//              EigArg *eig_arg, 
//              CnvFrmType cnv_frm = CNV_FRM_YES):
//------------------------------------------------------------------
int Fasqtad::FeigSolv(Vector **f_eigenv, Float lambda[], 
		    Float chirality[], int valid_eig[],
		    Float **hsum,
		    EigArg *eig_arg, 
		    CnvFrmType cnv_frm)
{
  int iter;
  char *fname = "FeigSolv(EigArg*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = eig_arg -> mass;
  cg_arg.RitzMatOper = eig_arg->RitzMatOper;
  int N_eig = eig_arg->N_eig;

  // Call constructor and solve for eigenvectors.
  // Use null pointers to fake out constructor.
  Vector *v1 = (Vector *)0;
  Vector *v2 = (Vector *)0;

  DiracOpAsqtad stag(*this, v1, v2, &cg_arg, CNV_FRM_NO);
  
  iter = stag.RitzEig(f_eigenv, lambda, valid_eig, eig_arg);
  
  // Modified for anisotropic lattices
  Float factor = GJP.XiV()/GJP.XiBare();
  // Chirality is trivial
  for(int i=0; i < N_eig; ++i) {
    chirality[i] = 1.0;
    lambda[i] *= factor;
  }
  // End modification

#if 0
  // !!! THIS DOES NOT WORK YET !!!
  // Slice-sum the eigenvector density to make a 1D vector
  if (eig_arg->print_hsum)
    for(i=0; i < N_eig; ++i)
      slice_sum_sq(hsum[i], f_eigenv[i], eig_arg->hsum_dir);
#endif

  // Return the number of iterations
  return iter;
}

//------------------------------------------------------------------
// SetPhi(Vector *phi, Vector *frm_e, Vector *frm_o, Float mass):
// It sets the pseudofermion field phi from frm_e, frm_o.
//------------------------------------------------------------------
void Fasqtad::SetPhi(Vector *phi, Vector *frm_e, Vector *frm_o, 
		   Float mass){
  char *fname = "SetPhi(V*,V*,V*,F)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = mass;
//  static int called=0;

  DiracOpAsqtad stag(*this, phi, frm_o, &cg_arg, CNV_FRM_NO);
  stag.Dslash(phi, frm_o, CHKB_ODD, DAG_NO);

  // Modified for anisotropic lattices
  //------------------------------------------------------------------
  fTimesV1MinusV2((IFloat *)phi, 2.*mass*GJP.XiBare()/GJP.XiV(), 
	(IFloat *)frm_e, (IFloat *)phi, e_vsize);

}



//------------------------------------------------------------------
// Float BhamiltonNode(Vector *boson, Float mass):
// The boson Hamiltonian of the node sublattice.
//------------------------------------------------------------------
Float Fasqtad::BhamiltonNode(Vector *boson, Float mass){
  char *fname = "BhamiltonNode(V*,F)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = mass;

  DiracOpAsqtad stag(*this, f_tmp, boson, &cg_arg, CNV_FRM_NO);
  Float ham;
  stag.MatPcDagMatPc(f_tmp, boson, &ham);

  return ham;
}



//------------------------------------------------------------------
// Fdslash(Vector *f_out, Vector *f_in, CgArg *cg_arg, CnvFrmType cnv_frm,
//                    int dir_flag) :
// Fdslash is the derivative part of the fermion matrix. 
// Fdslash calculates both odd-->even and even-->odd sites.
// dir_flag is flag which takes value 0 when all direction contribute to D,
// 1 - when only the special anisotropic direction contributes to D,
// 2 - when all  except the special anisotropic direction. 
//------------------------------------------------------------------
void Fasqtad::Fdslash(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
		    CnvFrmType cnv_frm, int dir_flag)
{
  int offset;
  char *fname = "Fdslash(V*,V*,CgArg*,CnvFrmType,int)";
  VRB.Func(cname,fname);
  
  DiracOpAsqtad stag(*this, f_out, f_in, cg_arg, cnv_frm);
  offset = GJP.VolNodeSites() * FsiteSize() / (2 * VECT_LEN);
  
  stag.Dslash(f_out, f_in+offset, CHKB_ODD, DAG_NO, dir_flag);
  stag.Dslash(f_out+offset, f_in, CHKB_EVEN, DAG_NO, dir_flag);

}

CPS_END_NAMESPACE
