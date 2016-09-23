#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of Fclover class.

  $Id: f_clover.C,v 1.21 2007/06/06 16:06:23 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2007/06/06 16:06:23 $
//  $Header: /space/cvs/cps/cps++/src/util/lattice/f_clover/f_clover.C,v 1.21 2007/06/06 16:06:23 chulwoo Exp $
//  $Id: f_clover.C,v 1.21 2007/06/06 16:06:23 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: f_clover.C,v $
//  $Revision: 1.21 $
//  $Source: /space/cvs/cps/cps++/src/util/lattice/f_clover/f_clover.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_clover.C
//
// Fclover is derived from FwilsonTypes and is relevant to
// clover Wilson fermions
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <stdlib.h>
#include <math.h>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/clover.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <comms/scu.h>
#include <util/dense_matrix.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// local constansts
//------------------------------------------------------------------
enum {
  CLOVER_MAT_SIZE = 72,   // COLORS^2 * lat.SpinComponents()^2 / 2;
  HALF_CLOVER_MAT_SIZE = 36
};


//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
Fclover::Fclover()
{
  cname = "Fclover";
  char *fname = "Fclover()";
  VRB.Func(cname,fname);


  //----------------------------------------------------------------
  // Allocate memory for the two clover matrices (one for each
  // checkerboard). Set the pointers to be the aux0_ptr and aux1_ptr
  // pointers of the base class.
  // ???Ping???
  //----------------------------------------------------------------
  int size = CLOVER_MAT_SIZE * GJP.VolNodeSites()*sizeof(Float) / 2;  

  aux0_ptr = (void *) smalloc(size);
  if(aux0_ptr == 0)
    ERR.Pointer(cname,fname, "aux0_ptr");
  VRB.Smalloc(cname,fname, "aux0_ptr",aux0_ptr, size);

  aux1_ptr = (void *) smalloc(size);
  if(aux1_ptr == 0)
    ERR.Pointer(cname,fname, "aux1_ptr");
  VRB.Smalloc(cname,fname, "aux1_ptr",aux1_ptr, size);


  //----------------------------------------------------------------
  // Do initializations before the clover library can be used
  //----------------------------------------------------------------
  static Clover clover_struct;
  f_dirac_op_init_ptr = &clover_struct;
  clover_init(&clover_struct);
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
Fclover::~Fclover()
{
  char *fname = "~Fclover()";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Un-initialize the clover library. Memory is set free here.
  //----------------------------------------------------------------
  clover_end((Clover *) f_dirac_op_init_ptr);

  // Free memory for the clover matrices
  // ???Ping???
  //----------------------------------------------------------------
  VRB.Sfree(cname,fname, "aux1_ptr",aux1_ptr);
  sfree(aux1_ptr);
  VRB.Sfree(cname,fname, "aux0_ptr",aux0_ptr);
  sfree(aux0_ptr);

}


//------------------------------------------------------------------
// FclassType Fclass():
// It returns the type of fermion class.
//------------------------------------------------------------------
FclassType Fclover::Fclass() const{
  return F_CLASS_CLOVER;
}


//------------------------------------------------------------------
// int FchkbEvl() :
// returns 0 => The fermion fields in the evolution
//      are defined on ODD-EVEN checkerboard (whole
//      lattice).
//------------------------------------------------------------------
int Fclover::FchkbEvl() const
{
  return 0;
}


//------------------------------------------------------------------
// int FmatEvlInv(Vector *f_out, Vector *f_in, 
//                CgArg *cg_arg, 
//                Float *true_res,
//		  CnvFrmType cnv_frm = CNV_FRM_YES):
// It calculates f_out where A * f_out = f_in and
// A is the preconditioned fermion matrix that appears
// in the HMC evolution (even/odd  preconditioning 
// of [Dirac^dag Dirac]. The inversion of the odd checkerboard
// piece is done with the conjugate gradient algorithm
// while the inversion of the even checkerboard is done
// using standard explicit hermitian matrix inversion of the
// clover matrix. cg_arg is the structure that contains
// all the control parameters for the CG, f_in is the
// fermion field source vector, f_out should be set to be
// the initial guess and on return is the solution.
// f_in and f_out are defined on the whole latice.
// If true_res !=0 the value of the true residual of the CG and
// is returned in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int Fclover::FmatEvlInv(Vector *f_out, Vector *f_in, 
			CgArg *cg_arg, 
			Float *true_res,
			CnvFrmType cnv_frm)
{
  int iter;
  char *fname = "FmatEvlInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  DiracOpClover clover(*this, f_out, f_in, cg_arg, cnv_frm);
  
//  iter = clover.InvCg(&(cg_arg->true_rsd));
  iter = clover.MatEvlInv(&(cg_arg->true_rsd));
  if (true_res) *true_res = cg_arg ->true_rsd;


  // Return the number of CG iterations
  return iter;
}

int Fclover::FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
			 int Nshift, int isz, CgArg **cg_arg, 
			 CnvFrmType cnv_frm, MultiShiftSolveType type, 
			 Float *alpha, Vector **f_out_d)
{
  char *fname = "MatMInv(Vector **out, Vector *in, Float *shift,...";
  VRB.Func(cname, fname);
  ERR.NotImplemented(cname,fname);
  return -1;
}


void Fclover::FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
			 Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm)
{

  char *fname = "FminResExt(V*, V*, V**, V**, int, CgArg *, CnvFrmType)";
  VRB.Func(cname,fname);
  
  DiracOpClover clover(*this, sol, source, cg_arg, cnv_frm);
  clover.MinResExt(sol,source,sol_old,vm,degree);
  
}

//------------------------------------------------------------------
// int FmatInv(Vector *f_out, Vector *f_in, 
//             CgArg *cg_arg, 
//             Float *true_res,
//             CnvFrmType cnv_frm = CNV_FRM_YES,
//             PreserveType prs_f_in = PRESERVE_YES):
// It calculates f_out where A * f_out = f_in and
// A is the fermion matrix (Dirac operator) with no 
// preconditioning. The preconditioned matrix is inverted 
// and from the result the non-preconditioned f_out is 
// calculated . The inversion of the odd checkerboard
// piece is done with the conjugate gradient algorithm
// while the inversion of the even checkerboard is done
// using standard explicit hermitian matrix inversion of the
// clover matrix. cg_arg is the structure that contains
// all the control parameters for the CG, f_in is the
// fermion field source vector, f_out should be set to be
// the initial guess and on return is the solution.
// f_in and f_out are defined on the whole latice.
// If true_res !=0 the value of the true residual of the CG and
// is returned in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// cnv_frm is used to specify if f_in should be converted 
// from canonical to fermion order and f_out from fermion 
// to canonical. 
// prs_f_in is used to specify if the source
// f_in should be preserved or not. If not the memory usage
// is less by the size of one fermion vector.
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int Fclover::FmatInv(Vector *f_out, Vector *f_in, 
		     CgArg *cg_arg, 
		     Float *true_res,
		     CnvFrmType cnv_frm,
		     PreserveType prs_f_in)
{
  int iter;
  char *fname = "FmatInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  DiracOpClover clover(*this, f_out, f_in, cg_arg, cnv_frm);
  
  iter = clover.MatInv(true_res, prs_f_in);
  
  // Return the number of iterations
  return iter;
}



//------------------------------------------------------------------
// int FeigSolv(Vector **f_eigenv, Float *lambda, int *valid_eig,
//              EigArg *eig_arg, 
//              CnvFrmType cnv_frm = CNV_FRM_YES):
//------------------------------------------------------------------
int Fclover::FeigSolv(Vector **f_eigenv, Float *lambda,
		      Float *chirality, int *valid_eig,
		      Float **hsum,
		      EigArg *eig_arg, 
		      CnvFrmType cnv_frm)
{
  int iter;
  char *fname = "FeigSolv(EigArg*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = eig_arg->mass;
  cg_arg.RitzMatOper = eig_arg->RitzMatOper;

  if(cnv_frm == CNV_FRM_YES)
    for(int i=0; i < eig_arg->N_eig; ++i)
      Fconvert(f_eigenv[i], WILSON, StrOrd());

  // Call constructor and solve for eigenvectors.
  // Use null pointers to fake out constructor.
  Vector *v1 = (Vector *)0;
  Vector *v2 = (Vector *)0;

  DiracOpClover clover(*this, v1, v2, &cg_arg, CNV_FRM_NO);
  
  iter = clover.RitzEig(f_eigenv, lambda, valid_eig, eig_arg);
  
  if(cnv_frm == CNV_FRM_YES)
    for(int i=0; i < eig_arg->N_eig; ++i)
      Fconvert(f_eigenv[i], CANONICAL, StrOrd());

  // Compute chirality
  int f_size = (GJP.VolNodeSites() * FsiteSize());
  Float factor = 4.0 + eig_arg->mass;
  v1 = (Vector *)smalloc(f_size*sizeof(Float));
  if (v1 == 0)
    ERR.Pointer(cname, fname, "v1");
  VRB.Smalloc(cname, fname, "v1", v1, f_size*sizeof(Float));

  for(int i=0; i < eig_arg->N_eig; ++i)
  {
    Gamma5(v1, f_eigenv[i], GJP.VolNodeSites());
    chirality[i] = f_eigenv[i]->ReDotProductGlbSum4D(v1, f_size);
    lambda[i] *= factor;
  }

  VRB.Sfree(cname, fname, "v1", v1);
  sfree(v1);

  // Return the number of iterations
  return iter;
}


//------------------------------------------------------------------
// SetPhi(Vector *phi, Vector *frm1, Vector *frm2, Float mass,
//        DagType dag):
// It sets the pseudofermion field phi from frm1, frm2.
// Modified - now returns the (trivial) value of the action
//------------------------------------------------------------------
Float Fclover::SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
		      Float mass, DagType dag){
  char *fname = "SetPhi(V*,V*,V*,F)";
  VRB.Func(cname,fname);

  CgArg cg_arg;
  cg_arg.mass = mass;

  if (phi == 0) ERR.Pointer(cname,fname,"phi") ;
  if (frm1 == 0) ERR.Pointer(cname,fname,"frm1") ;

  Fconvert(frm1,WILSON,CANONICAL);
  DiracOpClover clover(*this, frm1, frm2, &cg_arg, CNV_FRM_NO);

  const int half_sites = GJP.VolNodeSites()/2;  
  int vec_size = FsiteSize() * half_sites;    
  IFloat *frm1_even = (IFloat *)frm1 + vec_size;
  IFloat *phi_even = (IFloat *)phi + vec_size;  
  IFloat *A_even = (IFloat *)Aux0Ptr();

// phi_odd = (Aoo - kappa*kappa*Doe Aee^inv Deo)^dagger frm1_odd
  if (dag == DAG_YES) clover.MatPcDag(phi, frm1);  
  else clover.MatPc(phi, frm1);  

// Calculate the clover matrices for the EVEN checkerboard
  IFloat *Ap = A_even;
  for (int i = 0; i < GJP.VolNodeSites(); ++i) {
    mat_inv(Ap, Ap, 6, MAT_INV_ALG_LDL_CMPR, 0);    
    Ap += HALF_CLOVER_MAT_SIZE;         
  } 

// phi_even = Aee frm1_even
  clover_mat_mlt((IFloat *)phi_even, 
                   (const IFloat *)(A_even), 
                   (const IFloat *)frm1_even, 
                   half_sites);    

  return FhamiltonNode(frm1, frm1);
  
}


//------------------------------------------------------------------
// EvolveMomFforce(Matrix *mom, Vector *frm, Float mass, 
//                 Float dt):
// It evolves the canonical momentum mom by dt
// using the fermion force.
//------------------------------------------------------------------
ForceArg Fclover::EvolveMomFforce(Matrix *mom, Vector *frm, 
				  Float mass, Float dt) {
  char *fname = "EvolveMomFforce(M*,V*,F,F,F)";
  VRB.Func(cname,fname);
  
  Matrix *gauge = GaugeField() ;
  
  if (Colors() != 3)
    ERR.General(cname,fname,"Wrong nbr of colors.") ;

  if (SpinComponents() != 4)
    ERR.General(cname,fname,"Wrong nbr of spin comp.") ;
  
  if (mom == 0)
    ERR.Pointer(cname,fname,"mom") ;
  
  if (frm == 0)
    ERR.Pointer(cname,fname,"frm") ;
  
//------------------------------------------------------------------
// allocate space for four CANONICAL fermion fields.
//------------------------------------------------------------------
  int f_size = FsiteSize() * GJP.VolNodeSites() ;
  
  char *str_v1 = "v1" ;
  Vector *v1 = (Vector *)smalloc(f_size*sizeof(Float)) ;
  if (v1 == 0)
    ERR.Pointer(cname, fname, str_v1) ;
  VRB.Smalloc(cname, fname, str_v1, v1, f_size*sizeof(Float)) ;
  
  char *str_v2 = "v2" ;
  Vector *v2 = (Vector *)smalloc(f_size*sizeof(Float)) ;
  if (v2 == 0)
    ERR.Pointer(cname, fname, str_v2) ;
  VRB.Smalloc(cname, fname, str_v2, v2, f_size*sizeof(Float)) ;
  
  char *str_v3 = "v3" ;
  Vector *v3 = (Vector *)smalloc(f_size*sizeof(Float)) ;
  if (v3 == 0)
    ERR.Pointer(cname, fname, str_v3) ;
  VRB.Smalloc(cname, fname, str_v3, v3, f_size*sizeof(Float)) ;
  
  char *str_v4 = "v4" ;
  Vector *v4 = (Vector *)smalloc(f_size*sizeof(Float)) ;
  if (v4 == 0)
    ERR.Pointer(cname, fname, str_v4) ;
  VRB.Smalloc(cname, fname, str_v4, v4, f_size*sizeof(Float)) ;

//------------------------------------------------------------------
// allocate space for two CANONICAL fermion field on a site.
//------------------------------------------------------------------

  char *str_site_v1 = "site_v1";
  Float *site_v1 = (Float *)smalloc(FsiteSize()*sizeof(Float));
  if (site_v1 == 0)
    ERR.Pointer(cname, fname, str_site_v1) ;
  VRB.Smalloc(cname, fname, str_site_v1, site_v1,
    FsiteSize()*sizeof(Float)) ;

  char *str_site_v2 = "site_v2";
  Float *site_v2 = (Float *)smalloc(FsiteSize()*sizeof(Float));
  if (site_v2 == 0)
    ERR.Pointer(cname, fname, str_site_v2) ;
  VRB.Smalloc(cname, fname, str_site_v2, site_v2,
    FsiteSize()*sizeof(Float)) ;

  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;

//------------------------------------------------------------------
// Calculate the force vectors for the ODD checkerboard
//------------------------------------------------------------------

  {
    CgArg cg_arg ;
    cg_arg.mass = mass ;

    DiracOpClover clover(*this, v1, v2, &cg_arg, CNV_FRM_YES) ;
    clover.CalcHmdForceVecs(frm) ;
  }

//------------------------------------------------------------------
// Calculate the force vectors for the EVEN checkerboard
//------------------------------------------------------------------

  {
    CgArg cg_arg ;
    cg_arg.mass = mass ;

    DiracOpClover clover(*this, v3, v4, &cg_arg, CNV_FRM_YES) ;

    const int half_sites = GJP.VolNodeSites()/2;  
    int vec_size = FsiteSize() * half_sites;    

    IFloat *frm_even = (IFloat *)frm + vec_size;
    IFloat *v3_even = (IFloat *)v3 + vec_size;
    IFloat *v4_even = (IFloat *)v4 + vec_size;  
    IFloat *A_even = (IFloat *)Aux0Ptr();

// v3_even = frm_even
    ((Vector *) v3_even)->CopyVec((Vector *)frm_even, vec_size) ;

//--------------------------------------------------------------------
// Calculate the clover matrices for the EVEN checkerboard

    IFloat *Ap = A_even;
    for (int i = 0; i < GJP.VolNodeSites(); ++i) {
      mat_inv(Ap, Ap, 6, MAT_INV_ALG_LDL_CMPR, 0);    
      Ap += HALF_CLOVER_MAT_SIZE;         
    } 
    
// v4_even = Aee v3_even
    clover_mat_mlt((IFloat *)v4_even, 
                   (const IFloat *)(A_even), 
                   (const IFloat *)v3_even, 
                   half_sites);    
    
  } 
  int x, y, z, t, lx, ly, lz, lt ;
  
  lx = GJP.XnodeSites() ;
  ly = GJP.YnodeSites() ;
  lz = GJP.ZnodeSites() ;
  lt = GJP.TnodeSites() ;

//------------------------------------------------------------------
// Evolves the canonical momentum mom by dt
// using fermion force without the clover contribution. 
// Start by summing first over direction (mu) and then over site
// to allow SCU transfers to happen face-by-face in the outermost
// loop.
//------------------------------------------------------------------

  int mu ;

  Matrix tmp, f;

  for (mu=0; mu<4; mu++) {
    for (t=0; t<lt; t++)
    for (z=0; z<lz; z++)
    for (y=0; y<ly; y++)
    for (x=0; x<lx; x++) {
      int gauge_offset = x+lx*(y+ly*(z+lz*t)) ;
      int vec_offset = FsiteSize()*gauge_offset ;
      gauge_offset = mu+4*gauge_offset ;		

      Float *v1_plus_mu=NULL ;
      Float *v2_plus_mu=NULL ;
      int vec_plus_mu_offset = FsiteSize() ;
      
      Float kappa = 1.0 / (2.0 * (mass + 4.0));
      Float coeff = 2.0 * dt * kappa;
         
      switch (mu) {
        case 0 :
          vec_plus_mu_offset *= (x+1)%lx+lx*(y+ly*(z+lz*t)) ;
          if ((x+1) == lx) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;                        
            v2_plus_mu = site_v2 ;                        
            if (GJP.XnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
	  break ;
        case 1 :
	  vec_plus_mu_offset *= x+lx*((y+1)%ly+ly*(z+lz*t)) ;
	  if ((y+1) == ly) {
	    getPlusData( (IFloat *)site_v1,
			 (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
	    getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;                        
            v2_plus_mu = site_v2 ;                        
            if (GJP.YnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
          break ;
        case 2 :
          vec_plus_mu_offset *= x+lx*(y+ly*((z+1)%lz+lz*t)) ;
          if ((z+1) == lz) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            if (GJP.ZnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
          break ;
        case 3 :
          vec_plus_mu_offset *= x+lx*(y+ly*(z+lz*((t+1)%lt))) ;
          if ((t+1) == lt) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            if (GJP.TnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
      } // end switch mu

      sproj_tr[mu](   (IFloat *)&tmp,
                      (IFloat *)v1_plus_mu,
                      (IFloat *)v2+vec_offset, 1, 0, 0);

      sproj_tr[mu+4]( (IFloat *)&f,
                      (IFloat *)v2_plus_mu,
                      (IFloat *)v1+vec_offset, 1, 0, 0);

      tmp += f ;
      tmp *=coeff;

      f.DotMEqual(*(gauge+gauge_offset), tmp) ;

      tmp.Dagger(f) ;

      f.TrLessAntiHermMatrix(tmp) ;

      *(mom+gauge_offset) += f ;
      Float norm = f.norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);
    }
              
  }

//------------------------------------------------------------------
// Evolves the canonical momentum mom by dt
// using the clover contribution of fermion force 
//------------------------------------------------------------------

  EvolveMomFforceSupp(mom, v1, v2, v3, v4, mass, dt);
    
//------------------------------------------------------------------
// deallocate space for two CANONICAL fermion fields on a site.
//------------------------------------------------------------------

  VRB.Sfree(cname, fname, str_site_v2, site_v2) ;
  sfree(site_v2) ;

  VRB.Sfree(cname, fname, str_site_v1, site_v1) ;
  sfree(site_v1) ;

//------------------------------------------------------------------
// deallocate space for four CANONICAL fermion fields.
//------------------------------------------------------------------
  VRB.Sfree(cname, fname, str_v4, v4) ;
  sfree(v4) ;

  VRB.Sfree(cname, fname, str_v3, v3) ;
  sfree(v3) ;

  VRB.Sfree(cname, fname, str_v2, v2) ;
  sfree(v2) ;

  VRB.Sfree(cname, fname, str_v1, v1) ;
  sfree(v1) ;

  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();

  return ForceArg(L1, sqrt(L2), Linf);
}

ForceArg Fclover::RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
				   int isz, Float *alpha, Float mass, Float dt,
				   Vector **sol_d, ForceMeasure force_measure) {
  char *fname = "RHMC_EvolveMomFforce";

  ERR.General(cname,fname,"Not implemented\n");
  return ForceArg(0.0,0.0,0.0);
}

//------------------------------------------------------------------
// Float BhamiltonNode(Vector *boson, Float mass):
// The boson Hamiltonian of the node sublattice.
//------------------------------------------------------------------
Float Fclover::BhamiltonNode(Vector *boson, Float mass){
  char *fname = "BhamiltonNode(V*,F)";
  VRB.Func(cname,fname);

  CgArg cg_arg;
  cg_arg.mass = mass;
 
  if (boson == 0)
    ERR.Pointer(cname,fname,"boson");

  const int half_sites = GJP.VolNodeSites()/2;  
  int vec_size = FsiteSize() * half_sites;    
  int f_size = vec_size *2;
   
  Vector *bsn_tmp = (Vector *)
    smalloc(f_size*sizeof(Float));
 
  char *str_tmp = "bsn_tmp" ;
 
  if (bsn_tmp == 0)
    ERR.Pointer(cname,fname,str_tmp) ;
 
  VRB.Smalloc(cname,fname,str_tmp,bsn_tmp,f_size*sizeof(Float));

  DiracOpClover clover(*this, boson, bsn_tmp, &cg_arg, CNV_FRM_NO) ;

  IFloat *boson_even = (IFloat *)boson + vec_size;
  IFloat *bsn_tmp_even = (IFloat *)bsn_tmp + vec_size;  
  IFloat *A_even = (IFloat *)Aux0Ptr();

// bsn_tmp_odd = (Aoo - kappa*kappa*Doe Aee^inv Deo) boson_odd

  clover.MatPc(bsn_tmp, boson);  

//--------------------------------------------------------------------
// Calculate the clover matrices for the EVEN checkerboard

    IFloat *Ap = A_even;
    for (int i = 0; i < GJP.VolNodeSites(); ++i) {
      mat_inv(Ap, Ap, 6, MAT_INV_ALG_LDL_CMPR, 0);    
      Ap += HALF_CLOVER_MAT_SIZE;         
    } 

// bsn_tmp_even = Aee boson_even

  clover_mat_mlt((IFloat *)bsn_tmp_even, 
                   (const IFloat *)(A_even), 
                   (const IFloat *)boson_even, 
                   half_sites);    

  Float ret_val = bsn_tmp->NormSqNode(f_size) ;
 
  VRB.Sfree(cname,fname,str_tmp,bsn_tmp);
 
  sfree(bsn_tmp) ;
 
  return ret_val;
  
}



//------------------------------------------------------------------
// EvolveMomFforceSupp(Matrix *mom, Vector *v1, Vector *v2, 
//                 Float mass, Float dt):
// It evolves the canonical momentum mom by dt
// using the clover contribution of fermion force 
//------------------------------------------------------------------
void Fclover::EvolveMomFforceSupp(Matrix *mom, Vector *v1, Vector *v2,
		Vector *v3, Vector *v4, Float mass, Float dt){
  char *fname = "EvolveMomFforceSupp(M*,V*,F,F,F)";
  VRB.Func(cname,fname);

  Matrix *gauge = GaugeField() ;
 
  char *str_gt = "g_tmp" ;
  Matrix *g_tmp = (Matrix *)smalloc(6*sizeof(Matrix)) ;
  if (g_tmp == 0)
    ERR.Pointer(cname, fname, str_gt) ;
  VRB.Smalloc(cname, fname, str_gt, g_tmp, 6*sizeof(Matrix)) ;

  char *str_g1t = "g1_tmp" ;
  Matrix *g1_tmp = (Matrix *)smalloc(6*sizeof(Matrix)) ;
  if (g1_tmp == 0)
    ERR.Pointer(cname, fname, str_g1t) ;
  VRB.Smalloc(cname, fname, str_g1t, g1_tmp, 6*sizeof(Matrix)) ;

  char *str_g2t = "g2_tmp" ;
  Matrix *g2_tmp = (Matrix *)smalloc(6*sizeof(Matrix)) ;
  if (g2_tmp == 0)
    ERR.Pointer(cname, fname, str_g2t) ;
  VRB.Smalloc(cname, fname, str_g2t, g2_tmp, 6*sizeof(Matrix)) ;

  int mu, nu, loc[4], loc_sites[4];
  int size_Matrix = 2 * Colors() * Colors();

  loc_sites[0] = GJP.XnodeSites() ;
  loc_sites[1] = GJP.YnodeSites() ;
  loc_sites[2] = GJP.ZnodeSites() ;
  loc_sites[3] = GJP.TnodeSites() ;
  
  Matrix tmp, f, f1, f2, f3;

//------------------------------------------------------------------
// Suppose g1_{mu, nu}(m)=1/2 Tr_spin[ Sigma_{mu,nu} (v1 v2^dag + v2 v1^dag) ]
// and g2_{mu, nu}(m)=1/2 Tr_spin[ Sigma_{mu,nu} (v3 v4^dag + v4 v3^dag) ],
// G1* and G2* is used to store the matix for g1 and g2 at different site
  
  Matrix G1, G1_plus_mu, G1_plus_mu_plus_nu, G1_plus_mu_minus_nu,
    G1_plus_nu, G1_minus_nu;
  Matrix G2, G2_plus_mu, G2_plus_mu_plus_nu, G2_plus_mu_minus_nu,
    G2_plus_nu, G2_minus_nu;
	  
  Float kappa = 1.0 / (2.0 * (mass + 4.0));
	
  Float coeff = 0.25 * dt * kappa * GJP.CloverCoeff();
      
// start by summing first over direction (mu) and then over site
// to allow SCU transfers to happen face-by-face in the outermost
// loop.
//------------------------------------------------------------------

  for (mu=0; mu<4; mu++) {   
    for (loc[3]=0; loc[3]<loc_sites[3]; loc[3]++)  // t direction
    for (loc[2]=0; loc[2]<loc_sites[2]; loc[2]++)  // z direction
    for (loc[1]=0; loc[1]<loc_sites[1]; loc[1]++)  // y direction
    for (loc[0]=0; loc[0]<loc_sites[0]; loc[0]++) { // x direction
      int mod = (loc[0]+loc[1]+loc[2]+loc[3]) % 2;
      
      int gauge_base_offset = loc[0]+loc_sites[0]*
				 (loc[1]+loc_sites[1]*
				  (loc[2]+loc_sites[2]*
				   loc[3])) ;
      int vec_offset = FsiteSize() * gauge_base_offset;
      gauge_base_offset *=4;
      
      int gauge_offset = mu+gauge_base_offset ;		
      
      for (nu=0; nu<3; nu++)  { // Three directions other than mu
	int nu_dir = (nu+mu+1)%4;  // nu direction 

// offset needed to get the color matrix in g1 and g2
	int mu_nu = 3*mu+nu;  
	
//------------------------------------------------------------------
// Decide all kinds of offset
// And calculate all the color matrix G* constructed from
// the spinors v1, v2, v3 and v4

	// Site(m) m=(x,y,z,t)
	int loc_tmp[4]; 
	for (int i=0; i<4; i++)
	  loc_tmp[i]=loc[i];
		
	Sigmaproj_tr[mu_nu]( (IFloat *)&f,
			     (IFloat *)v1+vec_offset,
			     (IFloat *)v2+vec_offset, 1, 0, 0);
	
	Sigmaproj_tr[mu_nu]( (IFloat *)&G1,
			     (IFloat *)v2+vec_offset,
			     (IFloat *)v1+vec_offset, 1, 0, 0);

	G1 += f;

	if(mod == 0) {
	  Sigmaproj_tr[mu_nu]( (IFloat *)&f,
			       (IFloat *)v3+vec_offset,
			       (IFloat *)v4+vec_offset, 1, 0, 0);

	  Sigmaproj_tr[mu_nu]( (IFloat *)&G2,
			       (IFloat *)v4+vec_offset,
			       (IFloat *)v3+vec_offset, 1, 0, 0);

	  G2 += f;
	} else {
	  G2.ZeroMatrix();
	}
	
	// Site(m+mu) m=(x,y,z,t)
	loc_tmp[mu] = (loc_tmp[mu]+1)%loc_sites[mu];

	int gauge_plus_mu_offset = loc_tmp[0]+loc_sites[0]*
	                            (loc_tmp[1]+loc_sites[1]*
				       (loc_tmp[2]+loc_sites[2]*
					 loc_tmp[3])) ;
	int vec_plus_mu_offset = FsiteSize()*gauge_plus_mu_offset;
	gauge_plus_mu_offset *=4;
	
	Sigmaproj_tr[mu_nu]( (IFloat *)&f,
			     (IFloat *)v1+vec_plus_mu_offset,
			     (IFloat *)v2+vec_plus_mu_offset, 1, 0, 0);
	
	Sigmaproj_tr[mu_nu]( (IFloat *)&G1_plus_mu,
			     (IFloat *)v2+vec_plus_mu_offset,
			     (IFloat *)v1+vec_plus_mu_offset, 1, 0, 0);

	G1_plus_mu += f;

	if(mod == 1) {
	  Sigmaproj_tr[mu_nu]( (IFloat *)&f,
			       (IFloat *)v3+vec_plus_mu_offset,
			       (IFloat *)v4+vec_plus_mu_offset, 1, 0, 0);
	  
	  Sigmaproj_tr[mu_nu]( (IFloat *)&G2_plus_mu,
			       (IFloat *)v4+vec_plus_mu_offset,
			       (IFloat *)v3+vec_plus_mu_offset, 1, 0, 0);

	  G2_plus_mu += f;

	} else {
	  G2_plus_mu.ZeroMatrix();
	}

	// Site(m+mu+nu) m=(x,y,z,t)
	loc_tmp[nu_dir] = 
	  (loc_tmp[nu_dir]+1)%loc_sites[nu_dir];

	int gauge_plus_mu_plus_nu_offset = loc_tmp[0]+loc_sites[0]*
	                                     (loc_tmp[1]+loc_sites[1]*
					       (loc_tmp[2]+loc_sites[2]*
						loc_tmp[3])) ;
	
	int vec_plus_mu_plus_nu_offset = 
	  FsiteSize()*gauge_plus_mu_plus_nu_offset;
	gauge_plus_mu_plus_nu_offset *=4;
	
	Sigmaproj_tr[mu_nu]( (IFloat *)&f,
			     (IFloat *)v1+vec_plus_mu_plus_nu_offset,
			     (IFloat *)v2+vec_plus_mu_plus_nu_offset, 1, 0, 0);
	
	Sigmaproj_tr[mu_nu]( (IFloat *)&G1_plus_mu_plus_nu,
			     (IFloat *)v2+vec_plus_mu_plus_nu_offset,
			     (IFloat *)v1+vec_plus_mu_plus_nu_offset, 1, 0, 0);

	G1_plus_mu_plus_nu += f;
	  
	if(mod == 0) {
	  Sigmaproj_tr[mu_nu]( (IFloat *)&f,
			       (IFloat *)v3+vec_plus_mu_plus_nu_offset,
			       (IFloat *)v4+vec_plus_mu_plus_nu_offset, 1, 0, 0);
	  
	  Sigmaproj_tr[mu_nu]( (IFloat *)&G2_plus_mu_plus_nu,
			       (IFloat *)v4+vec_plus_mu_plus_nu_offset,
			       (IFloat *)v3+vec_plus_mu_plus_nu_offset, 1, 0, 0);

	  G2_plus_mu_plus_nu += f;

	} else {
	  G2_plus_mu_plus_nu.ZeroMatrix();
	}
	
	// Site(m+mu-nu) m=(x,y,z,t)
	loc_tmp[nu_dir] = 
	  (loc_tmp[nu_dir]-2+2*loc_sites[nu_dir])%loc_sites[nu_dir];

	int gauge_plus_mu_minus_nu_offset = loc_tmp[0]+loc_sites[0]*
	                                      (loc_tmp[1]+loc_sites[1]*
					       (loc_tmp[2]+loc_sites[2]*
						loc_tmp[3])) ;
	int vec_plus_mu_minus_nu_offset = 
	  FsiteSize()*gauge_plus_mu_minus_nu_offset;
	gauge_plus_mu_minus_nu_offset *=4;
	
	Sigmaproj_tr[mu_nu]( (IFloat *)&f,
			     (IFloat *)v1+vec_plus_mu_minus_nu_offset,
			     (IFloat *)v2+vec_plus_mu_minus_nu_offset, 1, 0, 0);
	
	Sigmaproj_tr[mu_nu]( (IFloat *)&G1_plus_mu_minus_nu,
			     (IFloat *)v2+vec_plus_mu_minus_nu_offset,
			     (IFloat *)v1+vec_plus_mu_minus_nu_offset, 1, 0, 0);

	G1_plus_mu_minus_nu += f;
	  
	if(mod == 0) {
	  Sigmaproj_tr[mu_nu]( (IFloat *)&f,
			       (IFloat *)v3+vec_plus_mu_minus_nu_offset,
			       (IFloat *)v4+vec_plus_mu_minus_nu_offset, 1, 0, 0);
	  
	  Sigmaproj_tr[mu_nu]( (IFloat *)&G2_plus_mu_minus_nu,
			       (IFloat *)v4+vec_plus_mu_minus_nu_offset,
			       (IFloat *)v3+vec_plus_mu_minus_nu_offset, 1, 0, 0);

	  G2_plus_mu_minus_nu += f;

	} else {
	  G2_plus_mu_minus_nu.ZeroMatrix();
	}

	// Site(m-nu) m=(x,y,z,t)
	loc_tmp[mu] = 
	  (loc_tmp[mu]-1+loc_sites[mu])%loc_sites[mu];

	int gauge_minus_nu_offset = loc_tmp[0]+loc_sites[0]*
				       (loc_tmp[1]+loc_sites[1]*
					(loc_tmp[2]+loc_sites[2]*
					 loc_tmp[3])) ;
	int vec_minus_nu_offset = FsiteSize()* gauge_minus_nu_offset;
	gauge_minus_nu_offset *=4;
	
	Sigmaproj_tr[mu_nu]( (IFloat *)&f,
			     (IFloat *)v1+vec_minus_nu_offset,
			     (IFloat *)v2+vec_minus_nu_offset, 1, 0, 0);
	
	Sigmaproj_tr[mu_nu]( (IFloat *)&G1_minus_nu,
			     (IFloat *)v2+vec_minus_nu_offset,
			     (IFloat *)v1+vec_minus_nu_offset, 1, 0, 0);

	G1_minus_nu += f;
	  
	if(mod == 1) {
	  Sigmaproj_tr[mu_nu]( (IFloat *)&f,
			       (IFloat *)v3+vec_minus_nu_offset,
			       (IFloat *)v4+vec_minus_nu_offset, 1, 0, 0);
	  
	  Sigmaproj_tr[mu_nu]( (IFloat *)&G2_minus_nu,
			       (IFloat *)v4+vec_minus_nu_offset,
			       (IFloat *)v3+vec_minus_nu_offset, 1, 0, 0);

	  G2_minus_nu += f;

	} else {
	  G2_minus_nu.ZeroMatrix();
	}
	
	// Site(m+nu) m=(x,y,z,t)
	loc_tmp[nu_dir] = 
	  (loc_tmp[nu_dir]+2)%loc_sites[nu_dir];

	int gauge_plus_nu_offset =  loc_tmp[0]+loc_sites[0]*
				       (loc_tmp[1]+loc_sites[1]*
					(loc_tmp[2]+loc_sites[2]*
					 loc_tmp[3])) ;
	int vec_plus_nu_offset = FsiteSize() * gauge_plus_nu_offset;
	gauge_plus_nu_offset *=4;
	
	Sigmaproj_tr[mu_nu]( (IFloat *)&f,
			     (IFloat *)v1+vec_plus_nu_offset,
			     (IFloat *)v2+vec_plus_nu_offset, 1, 0, 0);
	
	Sigmaproj_tr[mu_nu]( (IFloat *)&G1_plus_nu,
			     (IFloat *)v2+vec_plus_nu_offset,
			     (IFloat *)v1+vec_plus_nu_offset, 1, 0, 0);

	G1_plus_nu += f;
	  
	if(mod == 1) {
	  Sigmaproj_tr[mu_nu]( (IFloat *)&f,
			       (IFloat *)v3+vec_plus_nu_offset,
			       (IFloat *)v4+vec_plus_nu_offset, 1, 0, 0);
	  
	  Sigmaproj_tr[mu_nu]( (IFloat *)&G2_plus_nu,
			       (IFloat *)v4+vec_plus_nu_offset,
			       (IFloat *)v3+vec_plus_nu_offset, 1, 0, 0);

	  G2_plus_nu += f;

	} else {
	  G2_plus_nu.ZeroMatrix();
	}
	
//------------------------------------------------------------------
// Get the gauge links and g1 links      
// The link matrix are defined as following
// g_tmp[0] = U_nu^dagger(m)
// g_tmp[1] = U_mu^dagger(m+nu)
// g_tmp[2] = U_nu(m+mu)
// g_tmp[3] = U_nu(m-nu)
// g_tmp[4] = U_mu^daggger(m-nu)
// g_tmp[5] = U_nu^dagger(m+mu-nu)
// g1_tmp[0] = g1_{mu,nu}(m)
// g1_tmp[1] = g1_{mu,nu}(m+mu+nu)
// g1_tmp[2] = g1_{mu,nu}(m+mu-nu)
// g1_tmp[3] = g1_{mu,nu}(m+mu)
// g1_tmp[4] = g1_{mu,nu}(m+nu)
// g1_tmp[5] = g1_{mu,nu}(m-nu)
// if even site, calculate 
// g2_tmp[0] = g2_{mu,nu}(m)
// g2_tmp[1] = g2_{mu,nu}(m+mu+nu)
// g2_tmp[2] = g2_{mu,nu}(m+mu-nu)
// if odd site, calculate
// g2_tmp[3] = g2_{mu,nu}(m+mu)
// g2_tmp[4] = g2_{mu,nu}(m+nu)
// g2_tmp[5] = g2_{mu,nu}(m-nu)
//------------------------------------------------------------------

	// g_tmp[0] = U_nu^dagger(m)
	g_tmp[0].Dagger(*(gauge+gauge_base_offset+nu_dir));  
	
	// g1_tmp[0] = g1_{mu,nu}(m)
	// g2_tmp[0] = g2_{mu,nu}(m) for even site

	g1_tmp[0] = G1;

	if (mod == 0) {
	  g2_tmp[0] = G2;
	}
		
	// g_tmp[1] = U_mu^dagger(m+nu)
	// g1_tmp[4] = g1_{mu,nu}(m+nu)
	// g2_tmp[4] = g2_{mu,nu}(m+nu) for odd site
	// f1 stores =  g1_{mu,nu}(m+mu+nu) (g1_tmp[1])
	// g2_tmp[4] is used to store g2_{mu,nu}(m+mu+nu) for even site
	
	if ((loc[nu_dir]+1) == loc_sites[nu_dir]) {

	  getPlusData( (IFloat *) &f, 
		       (IFloat *) (gauge+gauge_plus_nu_offset+mu), 
		       size_Matrix, nu_dir);
	  g_tmp[1].Dagger(f);
	  
	  getPlusData( (IFloat *) &g1_tmp[4],
		       (IFloat *) &G1_plus_nu,
		       size_Matrix, nu_dir); 

	  getPlusData( (IFloat *) &f1, 
		       (IFloat *) &G1_plus_mu_plus_nu, 
		       size_Matrix, nu_dir);

	  if (mod ==1) {

	    getPlusData( (IFloat *) &g2_tmp[4],
			 (IFloat *) &G2_plus_nu,
			 size_Matrix, nu_dir); 

	  } else {  
	    //g2_tmp[4] is only meaningful when mod=1
	    // here it's used as a temporary storage place

	    getPlusData( (IFloat *) &g2_tmp[4],
                         (IFloat *) &G2_plus_mu_plus_nu, 
			 size_Matrix, nu_dir);	     
	  }
	  
	} else {
	  g_tmp[1].Dagger(*(gauge+gauge_plus_nu_offset+mu));  
	  g1_tmp[4] = G1_plus_nu;
	  f1 = G1_plus_mu_plus_nu;

	  if (mod == 1) {
	    g2_tmp[4] = G2_plus_nu;
	  } else {
	    g2_tmp[4] = G2_plus_mu_plus_nu;
	  }
	  
	}
	
	// g_tmp[3] = U_nu(m-nu)
	// g_tmp[4] = U_mu^daggger(m-nu)
	// g1_tmp[5] = g1_{mu,nu}(m-nu)
	// f2 stores g1_{mu,nu}(m+mu-nu) (g1_tmp[2])
	// g2_tmp[5] = g2_{mu,nu}(m-nu) for odd site
	// g2_tmp[5] stores g2_{mu,nu}(m+mu-nu) (g2_tmp[2]) for even site

	//	if(loc[nu_dir]-1 == 0) {
	if(loc[nu_dir] == 0) {
	  
	  getMinusData( (IFloat *) &g_tmp[3], 
			(IFloat *) (gauge+gauge_minus_nu_offset+nu_dir), 
			size_Matrix, nu_dir);

	  getMinusData( (IFloat *) &f, 
			(IFloat *) (gauge+gauge_minus_nu_offset+mu), 
			size_Matrix, nu_dir);

	  g_tmp[4].Dagger(f);  
	  
	  getMinusData( (IFloat *) &g_tmp[5], 
			(IFloat *)(gauge+gauge_plus_mu_minus_nu_offset+nu_dir), 
			size_Matrix, nu_dir);
	  f.Dagger(g_tmp[5]); // f stores U_nu^dagger(m+mu-nu)
	  
	  getMinusData( (IFloat *) &g1_tmp[5],
			(IFloat *) &G1_minus_nu,
			size_Matrix, nu_dir);

	  getMinusData( (IFloat *) &f2,
			(IFloat *) &G1_plus_mu_minus_nu,
			size_Matrix, nu_dir); 

	  if (mod ==1) {
	    getMinusData( (IFloat *) &g2_tmp[5],
			  (IFloat *) &G2_minus_nu,
			  size_Matrix, nu_dir);
	  } else {
	    getMinusData( (IFloat *) &g2_tmp[5],
			  (IFloat *) &G2_plus_mu_minus_nu,
			  size_Matrix, nu_dir);
	    }

	} else {
	  g_tmp[3] = *(gauge+gauge_minus_nu_offset+nu_dir);
	  g_tmp[4].Dagger(*(gauge+gauge_minus_nu_offset+mu));
	  f.Dagger(*(gauge+gauge_plus_mu_minus_nu_offset+nu_dir));

	  g1_tmp[5] = G1_minus_nu;
	  f2 = G1_plus_mu_minus_nu;

	  if (mod ==1) {
	    g2_tmp[5] = G2_minus_nu;
	  } else {
	    g2_tmp[5] = G2_plus_mu_minus_nu;
	  }
	}	

	// g_tmp[2] = U_nu(m+mu)
	// g_tmp[5] = U_nu^dagger(m+mu-nu)
	// g1_tmp[3] = g1_{mu,nu}(m+mu)
	// g1_tmp[1] = g1_{mu,nu}(m+mu+nu)
	// g1_tmp[2] = g1_{mu,nu}(m+mu-nu)
	// g2_tmp[3] = g2_{mu,nu}(m+mu) for odd site
	// g2_tmp[1] = g2_{mu,nu}(m+mu+nu) for even site
	// g2_tmp[2] = g2_{mu,nu}(m+mu-nu) for even site
		
	if ((loc[mu]+1) == loc_sites[mu]) {

	  getPlusData( (IFloat *) &g_tmp[2],
		       (IFloat *) (gauge+gauge_plus_mu_offset+nu_dir), 
		       size_Matrix, mu) ;

	  getPlusData( (IFloat *) &g_tmp[5],
		       (IFloat *) &f, 
		       size_Matrix, mu) ;

	  getPlusData( (IFloat *) &g1_tmp[3],
		       (IFloat *) &G1_plus_mu,
		       size_Matrix, mu);

	  getPlusData( (IFloat *) &g1_tmp[1],
		       (IFloat *) &f1,
		       size_Matrix, mu);

	  getPlusData( (IFloat *) &g1_tmp[2],
		       (IFloat *) &f2,
		       size_Matrix, mu);

	  if (mod ==0) {
	     getPlusData( (IFloat *) &g2_tmp[1],
			  (IFloat *) &g2_tmp[4],
			  size_Matrix, mu);

	     getPlusData( (IFloat *) &g2_tmp[2],
			  (IFloat *) &g2_tmp[5],
			  size_Matrix, mu);
	  } else {
	    getPlusData( (IFloat *) &g2_tmp[3],
			 (IFloat *) &G2_plus_mu,
			 size_Matrix, mu);
	  }
	} else {

	  g_tmp[2] = *(gauge+gauge_plus_mu_offset+nu_dir);
	  g_tmp[5] = f;

	  g1_tmp[3] = G1_plus_mu;
	  g1_tmp[1] = f1;
	  g1_tmp[2] = f2;

	  if(mod == 0) {
	    g2_tmp[1] = g2_tmp[4];
	    g2_tmp[2] = g2_tmp[5];
	  } else {
	    g2_tmp[3] = G2_plus_mu;
	  }
	}
	
//------------------------------------------------------------------
// Calculate the clover force term contribution
	
	// f1 = U_nu(m+mu)U_mu^dagger(m+nu)
	f1.DotMEqual(g_tmp[2], g_tmp[1]);
	// f2 = U_mu^dagger(m+nu)U_nu^dagger(m)
	f2.DotMEqual(g_tmp[1], g_tmp[0]); 
	// f3 = U_nu(m+mu)U_mu^dagger(m+nu)U_nu^dagger(m)
	f3.DotMEqual(f1, g_tmp[0]);

	if(nu == 0) { // Initialize tmp
	  tmp.DotMEqual(f3, g1_tmp[0]);
	} else {
	  f.DotMEqual(f3, g1_tmp[0]);
	  tmp += f;
	}
	
	f.DotMEqual(g1_tmp[3], f3);
	tmp += f;

	if (mod == 0) {
	  f.DotMEqual(f3, g2_tmp[0]);
	  tmp += f;

	  f3.DotMEqual(g2_tmp[1], f2);
	  f.DotMEqual(g_tmp[2], f3);
	  tmp += f;

	} else {
	  f.DotMEqual(g2_tmp[3], f3);
	  tmp += f;

	  f3.DotMEqual(g2_tmp[4], g_tmp[0]);
	  f.DotMEqual(f1, f3);
	  tmp += f;
	}

	f3.DotMEqual(g1_tmp[1], f2);
	f.DotMEqual(g_tmp[2], f3);
	tmp += f;

	f3.DotMEqual(g1_tmp[4], g_tmp[0]);
	f.DotMEqual(f1, f3);
	tmp += f;

	// f1 = U_nu^dagger(m+mu-nu)U_mu^dagger(m-nu)
	f1.DotMEqual(g_tmp[5], g_tmp[4]);
	// f2 = U_mu^dagger(m-nu)U_nu(m-nu)
	f2.DotMEqual(g_tmp[4], g_tmp[3]); 
	// f3 = U_nu^dagger(m+mu-nu)U_mu^dagger(m-nu)U_nu(m-nu)
	f3.DotMEqual(f1, g_tmp[3]);
	
	f.DotMEqual(f3, g1_tmp[0]);
	tmp -= f;

	f.DotMEqual(g1_tmp[3], f3);
	tmp -= f;
	
	if (mod ==0) {
	  f.DotMEqual(f3, g2_tmp[0]);
	  tmp -= f;

	  f3.DotMEqual(g2_tmp[2], f2);
	  f.DotMEqual(g_tmp[5], f3);
	  tmp -= f;

	} else {
	  f.DotMEqual(g2_tmp[3], f3);
	  tmp -= f;

	  f3.DotMEqual(g2_tmp[5], g_tmp[3]);
	  f.DotMEqual(f1, f3);
	  tmp -= f;
	}

	f3.DotMEqual(g1_tmp[2], f2);
	f.DotMEqual(g_tmp[5], f3);
	tmp -= f;

	f3.DotMEqual(g1_tmp[5], g_tmp[3]);
	f.DotMEqual(f1, f3);
	tmp -= f;
      }
      
      tmp *= coeff;
      
      f.DotMEqual(*(gauge+gauge_offset), tmp) ;
 
      tmp.Dagger(f) ;
      
      f.TrLessAntiHermMatrix(tmp) ;

      *(mom+gauge_offset) += f ;

    }
    
  }
  
//------------------------------------------------------------------
// deallocate space for g_tmp, g1_tmp, g2_tmp, g1 and g2.
//------------------------------------------------------------------

  VRB.Sfree(cname, fname, str_g2t, g2_tmp) ;
  sfree(g2_tmp) ;
  
  VRB.Sfree(cname, fname, str_g1t, g1_tmp) ;
  sfree(g1_tmp) ;
  
  VRB.Sfree(cname, fname, str_gt, g_tmp) ;
  sfree(g_tmp) ;

  return;
  
}

ForceArg Fclover::EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
		      Float mass, Float dt) {
  char *fname = "EvolveMomFforce(M*,V*,V*,F,F)";
  ERR.General(cname,fname,"Not Implemented\n");
  return ForceArg(0.0,0.0,0.0);
}


CPS_END_NAMESPACE
