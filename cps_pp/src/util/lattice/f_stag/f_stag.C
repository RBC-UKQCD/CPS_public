#include<config.h>
#include<util/qcdio.h>
#include<stdlib.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of Fstag class.

  $Id: f_stag.C,v 1.26 2013-06-07 19:26:34 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_stag/f_stag.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_stag.C
//
// Fstag is derived from FstagTypes and is relevant to
// staggered fermions
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/dirac_op.h>
#include <util/stag.h>
#include <util/gjp.h>
//#include <comms/nga_reg.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <comms/cbuf.h>
#include <math.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Initialize static variables.
//------------------------------------------------------------------



//  CRAM temp buffer
#ifdef _TARTAN

static Matrix *mp0 = (Matrix *)CRAM_SCRATCH_ADDR;	// ihdot
static Matrix *mp1 = mp0 + 1;
static Vector *vp0 = (Vector *) (CRAM_SCRATCH_ADDR+2*18*sizeof(IFloat));
static Vector *vp1 = vp0 + 1;
static Vector *vp2 = vp1 + 1;
static Vector *vp3 = vp2 + 1;

#else

static Matrix mt0;
static Matrix mt1;
static Matrix *mp0 = &mt0;		// ihdot
static Matrix *mp1 = &mt1;
static Vector vt0;
static Vector vt2;
static Vector vt3;
static Vector *vp0 = &vt0;
static Vector *vp2 = &vt2;
static Vector *vp3 = &vt3;

#endif 





//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
Fstag::Fstag()
{
  cname = "Fstag";
  char *fname = "Fstag()";
  VRB.Func(cname,fname);
  stag_dirac_init(GaugeField());
  
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
Fstag::~Fstag()
{
  char *fname = "~Fstag()";
  VRB.Func(cname,fname);
  stag_destroy_dirac_buf();
}


//------------------------------------------------------------------
// FclassType Fclass():
// It returns the type of fermion class.
//------------------------------------------------------------------
FclassType Fstag::Fclass() const{
  return F_CLASS_STAG;
}



//-------------------------------------------------------
//  // get V_mu(x) = U_mu(x)^dag * X(x+mu)
//  get V_mu(x) = U_mu(x) * X(x+mu)
//  end up with BANK3
//  begin with BANK4_BASE+BANK_SIZE
//-------------------------------------------------------
void Fstag::
getUDagX(Vector& v, const Vector *cvp, int *x, int mu) const
{
    Vector v_tmp1;
    Matrix *uoff = GaugeField()+GsiteOffset(x)+mu;

    setCbufCntrlReg(3, CBUF_MODE3);
    setCbufCntrlReg(4, CBUF_MODE4);


    //----------------------------------------
    // mp1 = U(x)
    //----------------------------------------
    // mp1->Dagger((IFloat *)uoff+BANK4_BASE+BANK_SIZE);
    moveMem(mp1, (IFloat *)uoff,
    	MATRIX_SIZE*sizeof(IFloat));

    // Modified for anisotropic lattices
    //------------------------------------------------------------------
    if (mu == GJP.XiDir()) 
      vecTimesEquFloat((IFloat *)mp1, GJP.XiVXi()/GJP.XiV(), MATRIX_SIZE);
    // End modification

    //----------------------------------------
    //  choose the right phase
    //----------------------------------------
    int eta_u = 0;
    for(int u = 0; u < mu; ++u) {
        eta_u += x[u];
    }
    eta_u = (eta_u & 1)? -1 : 1;

    const Vector *p = &v_tmp1;
    if(x[mu] == node_sites[mu]-1) {	// x+mu off node
	x[mu] = 0;
	getPlusData((IFloat *)&v_tmp1, (IFloat *)(cvp+FsiteOffsetChkb(x)),
	    VECT_LEN, mu);
        x[mu] = node_sites[mu]-1;
        if(bc[mu]) eta_u = -eta_u;

    } else { // x+mu on node
        x[mu]++;
	p = cvp+FsiteOffsetChkb(x);
	x[mu]--;
    }

    //----------------------------------------
    // copy *p to CRAM
    //----------------------------------------
    moveMem(vp0, (IFloat *)p, VECT_LEN*sizeof(IFloat));

    if(eta_u == -1) {
        vecNegative((IFloat *)mp1, (IFloat *)mp1, MATRIX_SIZE);
    }

    uDotXEqual((IFloat *)&v, (IFloat *)mp1, (IFloat *)vp0);
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
int Fstag::FmatEvlInv(Vector *f_out, Vector *f_in, 
		      CgArg *cg_arg, 
		      Float *true_res,
		      CnvFrmType cnv_frm)
{
  int iter;
  char *fname = "FmatEvlInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  printf("f_out = %e, f_in = %e\n", f_out->NormSqGlbSum(e_vsize), f_in->NormSqGlbSum(e_vsize));
  DiracOpStag stag(*this, f_out, f_in, cg_arg, cnv_frm);
  iter = stag.InvCg(&(cg_arg->true_rsd));
  if (true_res) *true_res = cg_arg ->true_rsd;

  stag.Dslash(f_tmp, f_out, CHKB_EVEN, DAG_NO);

  // Return the number of iterations
  return iter;
}

//------------------------------------------------------------------
// int FmatEvlMInv(Vector *f_out, Vector *f_in, 
//                Float shift[], int Nshift, 
//                CgArg **cg_arg, Float *true_res,
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
int Fstag::FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
		       int Nshift, int isz, CgArg **cg_arg,
		       CnvFrmType cnv_frm, MultiShiftSolveType type, 
		       Float *alpha, Vector **f_out_d)
{
  char *fname = "FmatMInv(V*, V*, .....)";
  VRB.Func(cname,fname);

  Float dot = f_in -> NormSqGlbSum4D(e_vsize);

  Float *RsdCG = new Float[Nshift];
  for (int s=0; s<Nshift; s++) RsdCG[s] = cg_arg[s]->stop_rsd;

  Float trueMass;
  massRenormalise(&(cg_arg[0]->mass), &trueMass, Nshift, shift, RENORM_FORWARDS);

  //Fake the constructor
  Vector *tmp;
  DiracOpStag stag(*this, tmp, f_in, cg_arg[0], cnv_frm);
  int iter = stag.MInvCG(f_out,f_in,dot,shift,Nshift,isz,RsdCG,type,alpha);  

  if (type == MULTI && f_out_d != 0)
    for (int s=0; s<Nshift; s++)
      stag.Dslash(f_out_d[s], f_out[s], CHKB_EVEN, DAG_NO);
  for (int s=0; s<Nshift; s++) cg_arg[s]->true_rsd = RsdCG[s];

  massRenormalise(&(cg_arg[0]->mass), &trueMass, Nshift, shift, RENORM_BACKWARDS);

  delete[] RsdCG;
  return iter;

}

//------------------------------------------------------------------
// Lattice class api to the chronological inverter
//------------------------------------------------------------------
void Fstag::FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
			 Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm)
{

  char *fname = "FminResExt(V*, V*, V**, V**, int, CgArg *, CnvFrmType)";
  VRB.Func(cname,fname);
  
  DiracOpStag stag(*this, sol, source, cg_arg, cnv_frm);
  stag.MinResExt(sol,source,sol_old,vm,degree);
  
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
int Fstag::FmatInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm,
		   PreserveType prs_f_in)
{
  int iter;
  char *fname = "FmatInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  DiracOpStag stag(*this, f_out, f_in, cg_arg, cnv_frm);
  
  iter = stag.MatInv(true_res, prs_f_in);
  
  // Return the number of iterations
  return iter;
}


//------------------------------------------------------------------
// int FeigSolv(Vector **f_eigenv, Float *lambda, int *valid_eig,
//              EigArg *eig_arg, 
//              CnvFrmType cnv_frm = CNV_FRM_YES):
//------------------------------------------------------------------
int Fstag::FeigSolv(Vector **f_eigenv, Float *lambda, 
		    Float *chirality, int *valid_eig,
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
#if 0
  // IS THIS NECESSARY ???
  if(cnv_frm == CNV_FRM_YES)
    for(int i=0; i < N_eig; ++i)
      Fconvert(f_eigenv[i], WILSON, StrOrd());
#endif

  // Call constructor and solve for eigenvectors.
  // Use null pointers to fake out constructor.
  Vector *v1 = (Vector *)0;
  Vector *v2 = (Vector *)0;

  DiracOpStag stag(*this, v1, v2, &cg_arg, CNV_FRM_NO);
 
  stag.RitzMat(f_eigenv[0],f_eigenv[1]);
  //for(int i=0;i<GJP.VolNodeSites()*3;i++){
   //  printf("VEC %d %g\n",i,*((Float*)f_eigenv[0]+i));
  //} 
  //exit(0);
  iter = stag.RitzEig(f_eigenv, lambda, valid_eig, eig_arg);
  
#if 0
  // IS THIS NECESSARY ???
  if(cnv_frm == CNV_FRM_YES)
    for(int i=0; i < N_eig; ++i)
      Fconvert(f_eigenv[i], CANONICAL, StrOrd());
#endif

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
// Solve  A * f_eigenv = lambda * f_eigenv where
// A is the fermion matrix (Dirac operator). The solution
// is done with the Lanczos algorithm. eig_arg is the
// structure that contains all the control parameters, f_eigenv
// is the fermion field eigenvectors, lambda are the
// returned eigenvalues.
// f_eigenv is defined on the whole lattice.
//------------------------------------------------------------------
int Fstag::FeigSolv(Vector **f_eigenv, Float *lambda,
                    LanczosArg *eig_arg,
                    CnvFrmType cnv_frm)
{
  int iter;
  char *fname = "FeigSolv(V*,F*,LanczosArg*,CnvFrmType)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = eig_arg->mass;
  //cg_arg.RitzMatOper = eig_arg->RitzMatOper;
  //int N_eig = eig_arg->N_eig;
  int nk = eig_arg->nk_lanczos_vectors;
  int np = eig_arg->np_lanczos_vectors;
  int maxiters = eig_arg->maxiters;
  Float stopres = eig_arg->stop_residual;
  MatrixPolynomialArg* cheby_arg = &(eig_arg->matpoly_arg);

  if(cnv_frm == CNV_FRM_YES) // convert only nk, not (nk+np)
    for(int i=0; i < nk; ++i)
      Fconvert(f_eigenv[i], STAG, StrOrd());

  // Call constructor and solve for eigenvectors.
  // Use null pointers to fake out constructor.
  Vector *v1 = (Vector *)0;
  Vector *v2 = (Vector *)0;
  DiracOpStag stag(*this, v1, v2, &cg_arg, CNV_FRM_NO);

  iter = stag.ImpResLanczos(f_eigenv, lambda,  eig_arg);

  if(cnv_frm == CNV_FRM_YES) for(int i=0; i < nk; ++i) // convert only nk, not (nk+np)
                               Fconvert(f_eigenv[i], CANONICAL, StrOrd());

  return iter;
}

//------------------------------------------------------------------
// SetPhi(Vector *phi, Vector *frm_e, Vector *frm_o, Float mass,
//        dag):
// It sets the pseudofermion field phi from frm_e, frm_o.
//------------------------------------------------------------------
Float Fstag::SetPhi(Vector *phi, Vector *frm_e, Vector *frm_o, 
		    Float mass, DagType dag){
  // dag is ignored for staggered

  char *fname = "SetPhi(V*,V*,V*,F)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = mass;

  DiracOpStag stag(*this, phi, frm_o, &cg_arg, CNV_FRM_NO);
  stag.Dslash(phi, frm_o, CHKB_ODD, DAG_NO);

  // Modified for anisotropic lattices
  //------------------------------------------------------------------
  fTimesV1MinusV2((IFloat *)phi, 2.*mass*GJP.XiBare()/GJP.XiV(), 
	(IFloat *)frm_e, (IFloat *)phi, e_vsize);
  // End modification
  
  return 0.0;
}


//------------------------------------------------------------------
// FforceSite(Matrix& force, Vector *frm, int *x, int mu):
// It calculates the fermion force at site x and direction mu.
// frm is the fermion field that resulted from the application
// of the inverter on the pseudofermion field.
//------------------------------------------------------------------
void Fstag::FforceSite(Matrix& force, Vector *frm, int *x, int mu)
{
  char *fname = "FforceSite(M&,V*,i*,i)";
//VRB.Func(cname,fname);

    int x_off = FsiteOffsetChkb(x);

    setCbufCntrlReg(3, CBUF_MODE3);

    //----------------------------------------
    // mp1 and mp2 free
    //
    //  calculate fermion part.
    //
    //----------------------------------------
    if((x[0]+x[1]+x[2]+x[3])%2) { // odd, assume base_odd = 0

        //----------------------------------------
	// U_u(x) X(x+u)
        //----------------------------------------
	getUDagX(*vp2, frm, x, mu);

        //----------------------------------------
	// *vp3 = *vp2;
        //----------------------------------------
	moveMem(vp3, vp2, VECT_LEN*sizeof(IFloat));

	vecMinusEquVec((IFloat *)vp3,
		(IFloat *)&f_tmp[x_off], VECT_LEN);

    } else {	// even
	getUDagX(*vp2, f_tmp, x, mu);

        //----------------------------------------
	// *vp3 = X[x]
        //----------------------------------------
	moveMem(vp3, (IFloat *)&frm[x_off],
	    VECT_LEN*sizeof(IFloat));
	*vp2 += *vp3;
    }

    force.Cross2(*vp2, *vp3);

    mp1->Dagger((IFloat *)&force);
    force.TrLessAntiHermMatrix(*mp1);
}

//!< Routine which allows bosonic staggered formulation.  Not sure why
//!< anyone would want this, but included for completeness.
void Fstag::BforceVector(Vector *in, CgArg *cg_arg) {

  int iter;
  char *fname = "BforceVector(V*)";
  VRB.Func(cname,fname);

  DiracOpStag stag(*this, f_tmp, in, cg_arg, CNV_FRM_NO);
  stag.Dslash(f_tmp, in, CHKB_EVEN, DAG_NO);

}


//------------------------------------------------------------------
// EvolveMomFforce(Matrix *mom, Vector *frm, Float mass, 
//                 Float step_size):
// It evolves the canonical momentum mom by step_size
// using the fermion force. 
//------------------------------------------------------------------
ForceArg Fstag::EvolveMomFforce(Matrix *mom, Vector *frm,
				Float mass, Float dt){
  char *fname = "EvolveMomFforce(M*,V*,F,F)";
  VRB.Func(cname,fname);
 
  setCbufCntrlReg(4, CBUF_MODE4);
  int x[4];
  
  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;

  for(x[0] = 0; x[0] < node_sites[0]; ++x[0]) {
    for(x[1] = 0; x[1] < node_sites[1]; ++x[1]) {
      for(x[2] = 0; x[2] < node_sites[2]; ++x[2]) {
	for(x[3] = 0; x[3] < node_sites[3]; ++x[3]) {
	  
	  Matrix *ihp = mom+GsiteOffset(x);
	  
	  for (int mu = 0; mu < 4; ++mu) {
	    FforceSite(*mp0, frm, x, mu);
	    fTimesV1PlusV2((IFloat *)(ihp+mu), dt,
			   (IFloat *)mp0,
			   (IFloat *)(ihp+mu), 18);
	    Float norm = mp0->norm();
	    Float tmp = sqrt(norm);
	    L1 += tmp;
	    L2 += norm;
	    Linf = (tmp>Linf ? tmp : Linf);
	  }
	}
      }
    }
  }

  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();

  return ForceArg(dt*L1, dt*sqrt(L2), dt*Linf);

}

ForceArg Fstag::EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
		      Float mass, Float step_size) {
  char *fname = "EvolveMomFforce(M*,V*,V*,F,F)";
  ERR.General(cname,fname,"Not Implemented\n");
  return ForceArg(0.0,0.0,0.0);
}

ForceArg Fstag::RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
				 int isz, Float *alpha, Float mass, Float dt,
				 Vector **sol_d, ForceMeasure force_measure) {
  char *fname = "RHMC_EvolveMomFforce";
  char *force_label;

  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;

  int g_size = GJP.VolNodeSites() * GsiteSize();

  Matrix *mom_tmp;

  if (force_measure == FORCE_MEASURE_YES) {
    mom_tmp = (Matrix*)smalloc(g_size*sizeof(Float),cname, fname, "mom_tmp");
    ((Vector*)mom_tmp) -> VecZero(g_size);
    force_label = new char[100];
  } else {
    mom_tmp = mom;
  }

  for (int i=0; i<degree; i++) {
    f_tmp -> CopyVec(sol_d[i], e_vsize);
    ForceArg Fdt = EvolveMomFforce(mom_tmp,sol[i],mass,dt*alpha[i]);
    if (force_measure == FORCE_MEASURE_YES) {
      sprintf(force_label, "Rational, mass = %e, pole = %d:", mass, i+isz);
      Fdt.print(dt,force_label);
    }
  }

  // If measuring the force, need to measure and then sum to mom
  if (force_measure == FORCE_MEASURE_YES) {
    for (int i=0; i<g_size/18; i++) {
      Float norm = (mom_tmp+i)->norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);
    }
    glb_sum(&L1);
    glb_sum(&L2);
    glb_max(&Linf);

    L1 /= 4.0*GJP.VolSites();
    L2 /= 4.0*GJP.VolSites();

    fTimesV1PlusV2((IFloat*)mom, 1.0, (IFloat*)mom_tmp, (IFloat*)mom, g_size);

    delete[] force_label;
    sfree(mom_tmp, cname, fname, "mom_tmp");
  }

  return ForceArg(L1, sqrt(L2), Linf);
}

//------------------------------------------------------------------
// Float BhamiltonNode(Vector *boson, Float mass):
// The boson Hamiltonian of the node sublattice.
//------------------------------------------------------------------
Float Fstag::BhamiltonNode(Vector *boson, Float mass){
  char *fname = "BhamiltonNode(V*,F)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = mass;

  DiracOpStag stag(*this, f_tmp, boson, &cg_arg, CNV_FRM_NO);
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
void Fstag::Fdslash(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
		    CnvFrmType cnv_frm, int dir_flag)
{
  int offset;
  char *fname = "Fdslash(V*,V*,CgArg*,CnvFrmType,int)";
  VRB.Func(cname,fname);
  
  DiracOpStag stag(*this, f_out, f_in, cg_arg, cnv_frm);
  offset = GJP.VolNodeSites() * FsiteSize() / (2 * VECT_LEN);
  
  stag.Dslash(f_out, f_in+offset, CHKB_ODD, DAG_NO, dir_flag);
  stag.Dslash(f_out+offset, f_in, CHKB_EVEN, DAG_NO, dir_flag);

}

//------------------------------------------------------------------
// FdMdmu(Vector *f_out, Vector *f_in, CgArg *cg_arg, CnvFrmType cnv_frm,
//                    int order) :
// FdMdmu is the derivative of the fermion matrix with respect to the
// chemical potential.
// order is the order of the derivative.
//------------------------------------------------------------------
void Fstag::FdMdmu(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
		    CnvFrmType cnv_frm, int order)
{
  int offset;
  char *fname = "FdMdmu(V*,V*,CgArg*,CnvFrmType,int)";
  VRB.Func(cname,fname);
  
  //DiracOpStag stag(*this, f_out, f_in, cg_arg, cnv_frm);
  //offset = GJP.VolNodeSites() * FsiteSize() / (2 * VECT_LEN);
  
  //stag.dMdmu(f_out, f_in+offset, CHKB_ODD, DAG_NO, order);
  //stag.dMdmu(f_out+offset, f_in, CHKB_EVEN, DAG_NO, order);
  ERR.General(cname,fname, "Not Implemented");

}

CPS_END_NAMESPACE
