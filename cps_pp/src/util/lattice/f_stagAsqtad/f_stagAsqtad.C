#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of FstagAsqtad class.

  $Id: f_stagAsqtad.C,v 1.2 2003-10-23 13:38:59 zs Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/dirac_op.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <comms/nga_reg.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <comms/cbuf.h>
CPS_START_NAMESPACE




//------------------------------------------------------------------
// Initialize static variables.
//------------------------------------------------------------------

//  CRAM temp buffer
#ifdef _TARTAN

static Matrix *mp0 = (Matrix *)CRAM_SCRATCH_ADDR;	// ihdot
static Matrix *mp1 = mp0 + 1;
// static Matrix *mp2 = mp1 + 1;
// static Matrix *mp3 = mp2 + 1;
static Vector *vp0 = (Vector *)
    (CRAM_SCRATCH_ADDR+2*MATRIX_SIZE*sizeof(IFloat));
static Vector *vp1 = vp0 + 1;
static Vector *vp2 = vp1 + 1;
static Vector *vp3 = vp2 + 1;

#else

static Matrix mt0;
static Matrix mt1;
static Matrix mt2;
static Matrix mt3;
static Matrix *mp0 = &mt0;		// ihdot
static Matrix *mp1 = &mt1;
// static Matrix *mp2 = &mt2;
// static Matrix *mp3 = &mt3;
static Vector vt0;
//static Vector vt1;
static Vector vt2;
static Vector vt3;
static Vector *vp0 = &vt0;
//static Vector *vp1 = &vt1;
static Vector *vp2 = &vt2;
static Vector *vp3 = &vt3;

#endif 









FstagAsqtad::FstagAsqtad(){
    cname = "FstagAsqtad";
    char *fname = "FstagAsqtad()";
    VRB.Func(cname,fname);
}


FstagAsqtad::~FstagAsqtad(){
    char *fname = "~FstagAsqtad()";
    VRB.Func(cname,fname);
}


FclassType FstagAsqtad::Fclass(){
    return F_CLASS_ASQTAD;
}


int FstagAsqtad::FmatEvlInv(Vector *f_out, Vector *f_in, 
			    CgArg *cg_arg, 
			    Float *true_res,
			    CnvFrmType cnv_frm)
{
    int iter;
    char *fname = "FmatEvlInv(CgArg*,V*,V*,F*,CnvFrmType)";
    VRB.Func(cname,fname);

    DiracOpStag stag(*this, f_out, f_in, cg_arg, cnv_frm);
  
    iter = stag.InvCg(true_res);

    stag.Dslash(f_tmp, f_out, CHKB_EVEN, DAG_NO);

    // Return the number of iterations
    return iter;
}


int FstagAsqtad::FmatEvlInv(Vector *f_out, Vector *f_in, 
			    CgArg *cg_arg, 
			    CnvFrmType cnv_frm)
{ return FmatEvlInv(f_out, f_in, cg_arg, 0, cnv_frm); }


int FstagAsqtad::FmatInv(Vector *f_out, Vector *f_in, 
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


int FstagAsqtad::FmatInv(Vector *f_out, Vector *f_in, 
			 CgArg *cg_arg, 
			 CnvFrmType cnv_frm,
			 PreserveType prs_f_in)
{ return FmatInv(f_out, f_in, cg_arg, 0, cnv_frm, prs_f_in); }


int FstagAsqtad::FeigSolv(Vector **f_eigenv, Float lambda[], 
			  Float chirality[], int valid_eig[],
			  Float **hsum,
			  EigArg *eig_arg, 
			  CnvFrmType cnv_frm)
{
    int iter;
    char *fname = "FeigSolv(EigArg*,V*,F*,CnvFrmType)";
    VRB.Func(cname,fname);
    CgArg cg_arg;
    cg_arg.mass = 0.0;
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

    DiracOpStag wilson(*this, v1, v2, &cg_arg, CNV_FRM_NO);
  
    iter = wilson.RitzEig(f_eigenv, lambda, valid_eig, eig_arg);
  
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

void FstagAsqtad::SetPhi(Vector *phi, Vector *frm_e, Vector *frm_o, 
			 Float mass){
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
}


void FstagAsqtad::FforceSite(Matrix& force, Vector *frm, int *x, int mu){

    char *fname = "FforceSite(M&,V*,i*,i)";
    ERR.NotImplemented(cname,fname);
	
}






Float FstagAsqtad::BhamiltonNode(Vector *boson, Float mass){
    char *fname = "BhamiltonNode(V*,F)";
    VRB.Func(cname,fname);
    CgArg cg_arg;
    cg_arg.mass = mass;

    DiracOpStag stag(*this, f_tmp, boson, &cg_arg, CNV_FRM_NO);
    Float ham;
    stag.MatPcDagMatPc(f_tmp, boson, &ham);

    return ham;
}



void FstagAsqtad::Fdslash(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
			  CnvFrmType cnv_frm, int dir_flag){
    int offset;
    char *fname = "Fdslash(V*,V*,CgArg*,CnvFrmType,int)";
    VRB.Func(cname,fname);
  
    DiracOpStag stag(*this, f_out, f_in, cg_arg, cnv_frm);
    offset = GJP.VolNodeSites() * FsiteSize() / (2 * VECT_LEN);
  
    stag.Dslash(f_out, f_in+offset, CHKB_ODD, DAG_NO, dir_flag);
    stag.Dslash(f_out+offset, f_in, CHKB_EVEN, DAG_NO, dir_flag);

}


CPS_END_NAMESPACE
