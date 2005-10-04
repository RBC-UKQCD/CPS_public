#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of FdwfBase class.

  $Id: f_dwf_base.C,v 1.25 2005-10-04 05:39:56 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_dwf_base/f_dwf_base.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_dwf_base.C
//
// FdwfBase is derived from FwilsonTypes and is relevant to
// domain wall fermions, moved from old Fdwf class
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <math.h>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/dwf.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/vector.h>
#include <util/random.h>
#include <util/error.h>
#include <util/time.h>
#include <util/enum_func.h>
#include <comms/scu.h> // GRF
#include <comms/glb.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// Initialize static variables.
//------------------------------------------------------------------


//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
FdwfBase::FdwfBase()
{
  cname = "FdwfBase";
  char *fname = "FdwfBase()";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Check if anisotropy is present and exit since FdwfBase has
  // not been tested for anisotropic lattices.
  //----------------------------------------------------------------
  if(GJP.XiBare() != 1 ||
     GJP.XiV()    != 1 ||
     GJP.XiVXi()  != 1   ){
    ERR.General(cname,fname,
    "XiBare=%g, XiV=%g, XiVXi=%g : FdwfBase has not been tested with anisotropy\n",
                GJP.XiBare(), GJP.XiV(), GJP.XiVXi());
  }

  //----------------------------------------------------------------
  // Do initializations before the dwf library can be used
  //----------------------------------------------------------------
  static Dwf dwf_struct;
  f_dirac_op_init_ptr = &dwf_struct;
  dwf_init((Dwf *) f_dirac_op_init_ptr);
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
FdwfBase::~FdwfBase()
{
  char *fname = "~FdwfBase()";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Un-initialize the dwf library. Memory is set free here.
  //----------------------------------------------------------------
  dwf_end((Dwf *) f_dirac_op_init_ptr);
}


//------------------------------------------------------------------
// FclassType Fclass(void):
// It returns the type of fermion class.
//------------------------------------------------------------------
FclassType FdwfBase::Fclass(void) const {
  return F_CLASS_DWF;
}


#if 0
//------------------------------------------------------------------
// int ExactFlavors() : 
// Returns the number of exact flavors of the matrix that
// is inverted during a molecular dynamics evolution.
//------------------------------------------------------------------
int FdwfBase::ExactFlavors(void) const
{
  return 2;
}


//------------------------------------------------------------------
// int SpinComponents() : 
// Returns the number of spin components.
//------------------------------------------------------------------
int FdwfBase::SpinComponents(void) const
{
  return 4;
}
#endif


//------------------------------------------------------------------
// int FsiteSize() : 
// Returns the number of fermion field components 
// (including real/imaginary) on a site of the 4-D lattice.
//------------------------------------------------------------------
int FdwfBase::FsiteSize(void) const
{
  return 2 * Colors() * SpinComponents() * GJP.SnodeSites();  
  // re/im * colors * spin_components * Ls
}

//------------------------------------------------------------------
// int FchkbEvl() :
// returns 1 => The fermion fields in the evolution
//      or the CG that inverts the evolution matrix
//      are defined on a single checkerboard (half the 
//      lattice).
//------------------------------------------------------------------
int FdwfBase::FchkbEvl(void) const
{
  return 1;
}


//------------------------------------------------------------------
// int FmatEvlInv(Vector *f_out, Vector *f_in, 
//                CgArg *cg_arg, 
//                Float *true_res,
//		  CnvFrmType cnv_frm = CNV_FRM_YES):
// It calculates f_out where A * f_out = f_in and
// A is the preconditioned fermion matrix that appears
// in the HMC evolution (even/odd  preconditioning 
// of [Dirac^dag Dirac]. The inversion is done
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
int FdwfBase::FmatEvlInv(Vector *f_out, Vector *f_in, 
		     CgArg *cg_arg, 
		     Float *true_res,
		     CnvFrmType cnv_frm)
{
  int iter;
  char *fname = "FmatEvlInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  DiracOpDwf dwf(*this, f_out, f_in, cg_arg, cnv_frm);
  
  iter = dwf.InvCg(&(cg_arg->true_rsd));
  if (true_res) *true_res = cg_arg ->true_rsd;

  // Return the number of iterations
  return iter;
}


//------------------------------------------------------------------
// Overloaded function is same as original but with true_res=0;
//------------------------------------------------------------------
int FdwfBase::FmatEvlInv(Vector *f_out, Vector *f_in, 
		     CgArg *cg_arg, 
		     CnvFrmType cnv_frm)
{ return FmatEvlInv(f_out, f_in, cg_arg, 0, cnv_frm); }


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
int FdwfBase::FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
			  int Nshift, int isz, CgArg **cg_arg, 
			  CnvFrmType cnv_frm, MultiShiftSolveType type, 
			  Float *alpha, Vector **f_out_d)
{
  int iter;
  char *fname = "FmatMInv(V**, V*, .....)";
  VRB.Func(cname,fname);

  int f_size = GJP.VolNodeSites() * FsiteSize() / (FchkbEvl()+1);
  Float dot = f_in -> NormSqGlbSum(f_size);

  Float *RsdCG = new Float[Nshift];
  for (int s=0; s<Nshift; s++) RsdCG[s] = cg_arg[s]->stop_rsd;

  //Fake the constructor
  DiracOpDwf dwf(*this, f_out[0], f_in, cg_arg[0], cnv_frm);

  int return_value= dwf.MInvCG(f_out,f_in,dot,shift,Nshift,isz,RsdCG,type,alpha);  
  for (int s=0; s<Nshift; s++) cg_arg[s]->true_rsd = RsdCG[s];  

  delete[] RsdCG;
  return return_value;
}

//------------------------------------------------------------------
// Lattice class api to the chronological inverter
//------------------------------------------------------------------
void FdwfBase::FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
			 Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm)
{

  char *fname = "FminResExt(V*, V*, V**, V**, int, CgArg *, CnvFrmType)";
  VRB.Func(cname,fname);
  
  DiracOpDwf dwf(*this, sol, source, cg_arg, cnv_frm);
  dwf.MinResExt(sol,source,sol_old,vm,degree);
  
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
int FdwfBase::FmatInv(Vector *f_out, Vector *f_in, 
		  CgArg *cg_arg, 
		  Float *true_res,
		  CnvFrmType cnv_frm,
		  PreserveType prs_f_in)
{
  int iter;
  char *fname = "FmatInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  DiracOpDwf dwf(*this, f_out, f_in, cg_arg, cnv_frm);
    
  iter = dwf.MatInv(true_res, prs_f_in);

  // Return the number of iterations
  return iter;
}


//------------------------------------------------------------------
// Overloaded function is same as original but with true_res=0;
//------------------------------------------------------------------
int FdwfBase::FmatInv(Vector *f_out, Vector *f_in, 
		  CgArg *cg_arg, 
		  CnvFrmType cnv_frm,
		  PreserveType prs_f_in)
{ return FmatInv(f_out, f_in, cg_arg, 0, cnv_frm, prs_f_in); }


//------------------------------------------------------------------
/*!
  \param five The 5-dimensional field.
  \param four The 4-dimensional field.
  \param s_u The global 5th direction (s) coordinate where the
   upper two components (right chirality) of the 5-dim. field
   take the values of those of the 4-dim. field.
  \param s_l The global 5th direction (s) coordinate where the
   lower two components (left chirality) of the 5-dim. field
   take the values of those of the 4-dim. field.
   \post The 5-dim field is zero everywhere except where the global
   5th direction coordinate (s) is \a s_l or \a s_u, where it takes the values
   explained above.
*/
//------------------------------------------------------------------
void FdwfBase::Ffour2five(Vector *five, Vector *four, int s_u, int s_l, int Ncb)
{
  int x;
  int i;
  Float *field_4D;
  Float *field_5D;
  char *fname = "Ffour2five(V*,V*,i,i)";
  VRB.Func(cname,fname);


//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------
  int f_size = GJP.VolNodeSites() * FsiteSize()*Ncb/2;
  int ls = GJP.SnodeSites();
  int vol_4d = GJP.VolNodeSites()*Ncb/2;
  int ls_stride = 24 * vol_4d;
  int s_u_local = s_u % GJP.SnodeSites();
  int s_l_local = s_l % GJP.SnodeSites();
  int s_u_node = s_u / GJP.SnodeSites();
  int s_l_node = s_l / GJP.SnodeSites();


//------------------------------------------------------------------
// Set *five using the 4D field *four. 
//------------------------------------------------------------------

  // Set all components of the 5D field to zero.
  //---------------------------------------------------------------
  field_5D  = (Float *) five;
  for(i=0; i<f_size; i++){
    field_5D[i]  = 0.0;
  }

  // Do the two upper spin components if s_u is in the node
  //---------------------------------------------------------------
  if( s_u_node == GJP.SnodeCoor() ){
    field_4D  = (Float *) four;
    field_5D  = (Float *) five;
    field_5D  = field_5D  + s_u_local * ls_stride;
    for(x=0; x<vol_4d; x++){
      for(i=0; i<12; i++){
	field_5D[i]  = field_4D[i];
      }
      field_4D  = field_4D  + 24;
      field_5D  = field_5D  + 24;
    }
  }

  // Do the two lower spin components if s_l is in the node
  //----------------------------------------------------------------
  if( s_l_node == GJP.SnodeCoor() ){
    field_4D  = (Float *) four;
    field_5D  = (Float *) five;
    field_4D  = field_4D  + 12;
    field_5D  = field_5D  + 12 + s_l_local * ls_stride;
    for(x=0; x<vol_4d; x++){
      for(i=0; i<12; i++){
	field_5D[i]  = field_4D[i];
      }
      field_4D  = field_4D  + 24;
      field_5D  = field_5D  + 24;
    }
  }

}


//------------------------------------------------------------------
/*!
  \param four The 4-dimensional field.
  \param five The 5-dimensional field.
  \param s_u The global 5th direction (s) coordinate where 
   the values of the upper two components (right chirality) of the 5-dim. field
   are taken by those of the 4-dim. field.
  \param s_l The global 5th direction coordinate (s) where the values of 
   the lower two components (left chirality) of the 5-dim. field
   are taken by  those of the 4-dim. field.
   \post The 5-dim field is zero everywhere except where the global
   5th direction coordinate is \a s_l or \a s_u, where it takes the values
   explained above.
   \post An identical 4-dim. field is reproduced on all nodes in the s
   direction.
*/
//------------------------------------------------------------------
void FdwfBase::Ffive2four(Vector *four, Vector *five, int s_u, int s_l, int Ncb)
{
  int x;
  int i;
  Float *field_4D;
  Float *field_5D;
  char *fname = "Ffive2four(V*,V*,i,i)";
  VRB.Func(cname,fname);


//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------
  int ls = GJP.SnodeSites();
  int f_size = GJP.VolNodeSites() * FsiteSize()*Ncb / (ls*2);
  int vol_4d = GJP.VolNodeSites()*Ncb/2;
  int ls_stride = 24 * vol_4d;
  int s_u_local = s_u % GJP.SnodeSites();
  int s_l_local = s_l % GJP.SnodeSites();
  int s_u_node = s_u / GJP.SnodeSites();
  int s_l_node = s_l / GJP.SnodeSites();


//------------------------------------------------------------------
// Set *four using the 5D field *five. 
//------------------------------------------------------------------

  // Set all components of the 4D field to zero.
  //---------------------------------------------------------------
  field_4D  = (Float *) four;
  for(i=0; i<f_size; i++){
    field_4D[i]  = 0.0;
  }

  // Do the two upper spin components if s_u is in the node
  //---------------------------------------------------------------
  if( s_u_node == GJP.SnodeCoor() ){
    field_4D = (Float *) four;
    field_5D = (Float *) five;
    field_5D = field_5D + s_u_local * ls_stride;
    for(x=0; x<vol_4d; x++){
      for(i=0; i<12; i++){
	field_4D[i] = field_5D[i];
      }
      field_4D = field_4D + 24;
      field_5D = field_5D + 24;
    }
  }

  // Do the two lower spin components if s_l is in the node
  //----------------------------------------------------------------
  if( s_l_node == GJP.SnodeCoor() ){
    field_4D = (Float *) four;
    field_5D = (Float *) five;
    field_4D = field_4D + 12;
    field_5D = field_5D + 12 + s_l_local * ls_stride;
    for(x=0; x<vol_4d; x++){
      for(i=0; i<12; i++){
	field_4D[i] = field_5D[i];
      }
      field_4D = field_4D + 24;
      field_5D = field_5D + 24;
    }
  }

  // Sum along s direction to get the same 4D field in all 
  // s node slices.
  //----------------------------------------------------------------
  if( GJP.Snodes() > 1) {
    Float sum;
    field_4D  = (Float *) four;
    for(i=0; i<f_size; i++){
      sum = field_4D[i];
      glb_sum_dir(&sum, 4);
      field_4D[i] = sum;    
    }
  }

}

//------------------------------------------------------------------
// int FeigSolv(Vector **f_eigenv, Float *lambda, int *valid_eig,
//              EigArg *eig_arg, 
//              CnvFrmType cnv_frm = CNV_FRM_YES):
// It solve  A * f_eigenv = lambda * f_eigenv where
// A is the fermion matrix (Dirac operator). The solution
// is done with the Ritz algorithm. eig_arg is the 
// structure that contains all the control parameters, f_eigenv
// is the fermion field eigenvectors, lambda are the
// returned eigenvalues.
// f_eigenv is defined on the whole lattice.
// The function returns the total number of Ritz iterations.
//------------------------------------------------------------------
int FdwfBase::FeigSolv(Vector **f_eigenv, Float *lambda,
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
  int N_eig = eig_arg->N_eig;
 
  if(cnv_frm == CNV_FRM_YES)
    for(int i=0; i < N_eig; ++i)
      Fconvert(f_eigenv[i], WILSON, StrOrd());
 
  // Call constructor and solve for eigenvectors.
  // Use null pointers to fake out constructor.
  Vector *v1 = (Vector *)0;
  Vector *v2 = (Vector *)0;
 
  DiracOpDwf dwf(*this, v1, v2, &cg_arg, CNV_FRM_NO);
  
  iter = dwf.RitzEig(f_eigenv, lambda, valid_eig, eig_arg);
 
  if(cnv_frm == CNV_FRM_YES)
    for(int i=0; i < N_eig; ++i)
      Fconvert(f_eigenv[i], CANONICAL, StrOrd());
 
  // rescale eigenvalues to normal convention
  Float factor = 5. - GJP.DwfHeight();
  int i;
  //for(i=0; i<N_eig; ++i)
  //lambda[i] *= factor;

  // calculate chirality
  int Ncb = NumChkb(cg_arg.RitzMatOper);
  int f_size = GJP.VolNodeSites()*2*Colors()*SpinComponents()*Ncb/2;
  Vector *four = (Vector *) smalloc (cname,fname, "four", f_size * sizeof(Float));
  Vector *fourg5 = (Vector *) smalloc (cname,fname, "fourg5", f_size * sizeof(Float));
  Float help;

  for (i=0; i<N_eig; i++) {
    Ffive2four (four, f_eigenv[i], 0, GJP.Snodes()*GJP.SnodeSites()-1,Ncb);

    // normalize four
    factor=four->NormSqNode(f_size);
    glb_sum(&factor);
    factor=1./sqrt(factor);
    four->VecTimesEquFloat(factor,f_size);
    Gamma5(fourg5,four,GJP.VolNodeSites()*Ncb/2);
    chirality[i]= four->ReDotProductNode(fourg5, f_size);
    glb_sum(&chirality[i]);
  }

  // calculate hsum
  // slice sum the eigenvector density to make a 1D vector
  if (eig_arg->print_hsum) {

    if(cnv_frm == CNV_FRM_NO)
      for(int i=0; i < N_eig; ++i)
	Fconvert(f_eigenv[i], CANONICAL, StrOrd());

    Float *f_in = (Float *) smalloc (cname, fname, "f_in", 
		GJP.VolNodeSites()*GJP.SnodeSites()*sizeof(Float));
    
    for(i=0; i < N_eig; ++i) {
      IFloat *fp= (IFloat *) f_eigenv[i];
      int j;
      for (j=0; j<GJP.VolNodeSites()*GJP.SnodeSites(); j++, fp+= 24) 
	f_in[j]= Float (dotProduct (fp,fp,24));

      if (i==0) {
	for (j=0; j<GJP.VolNodeSites()*GJP.SnodeSites(); j++)
	  printf ("%f ", f_in[j]);
      }
      printf ("\n");

      f_eigenv[i]->SliceArraySumFive (hsum[i], f_in, eig_arg->hsum_dir);
    }

    sfree(cname, fname, "f_in", f_in);
  }
  sfree(cname,fname, "four",four);
  sfree(cname,fname, "fourg5",fourg5);
  // Return the number of iterations
  return iter;
}


//------------------------------------------------------------------
// SetPhi(Vector *phi, Vector *frm1, Vector *frm2, Float mass):
// It sets the pseudofermion field phi from frm1, frm2.
// Note that frm2 is not used.
//------------------------------------------------------------------
Float FdwfBase::SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
		  Float mass){
  char *fname = "SetPhi(V*,V*,V*,F)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = mass;

  if (phi == 0)
    ERR.Pointer(cname,fname,"phi") ;

  if (frm1 == 0)
    ERR.Pointer(cname,fname,"frm1") ;

  DiracOpDwf dwf(*this, frm1, 0, &cg_arg, CNV_FRM_NO) ;
#if 0
{ IFloat *tmp = (IFloat *)frm1;
  printf("frm1[0]=%e\n",*tmp);}
#endif
  dwf.MatPcDag(phi, frm1) ;
#if 0
{ IFloat *tmp = (IFloat *)phi;
  printf("phi[0]=%e\n",*tmp);}
#endif

  return FhamiltonNode(frm1, frm1);
}

#undef PROFILE


void FdwfBase::RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
				    int isz, Float *alpha, Float mass, 
				    Float dt, Vector **sol_d) {
  char *fname = "RHMC_EvolveMomFforce";
  char *fer_force = (char*)smalloc(cname,fname,"fer_force",100*sizeof(char));

#ifdef PROFILE
  Float time = -dclock();
  ForceFlops=0;
#endif

  // Temporary fix for the moment.
  int f_size = GJP.VolNodeSites() * FsiteSize() / (FchkbEvl()+1);
  int g_size = GJP.VolNodeSites() * GsiteSize();

  Matrix *mom_old = (Matrix*)smalloc(g_size*sizeof(Float));

  for (int i=0; i<degree; i++){
    // Copy the momenta before the update
    ((Vector*)mom_old)->CopyVec((Vector*)mom,g_size);
    FdwfBase::EvolveMomFforce(mom,sol[i],mass,alpha[i]*dt);
    sprintf(fer_force, "Fermion (pole = %d):", i);
    ForceMagnitude(mom, mom_old, mass, dt, fer_force);
  }

  sfree(fer_force);
  sfree(mom_old);

#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops,time);
#endif

}

//------------------------------------------------------------------
// Float FhamiltonNode(Vector *phi, Vector *chi):
// The fermion Hamiltonian of the node sublattice
// chi must be the solution of Cg with source phi.	       
//------------------------------------------------------------------
Float FdwfBase::FhamiltonNode(Vector *phi, Vector *chi){
  char *fname = "FhamiltonNode(V*,V*)";
  VRB.Func(cname,fname);

  if (phi == 0)
    ERR.Pointer(cname,fname,"phi") ;

  if (chi == 0)
    ERR.Pointer(cname,fname,"chi") ;

  int f_size = GJP.VolNodeSites() * FsiteSize() / 2 ;

  Float ret_val;
  ret_val = phi->ReDotProductNode(chi, f_size ) ;

  // Sum accross s nodes in case Snodes() != 1
  glb_sum_dir(&ret_val, 4) ;

  return ret_val ;

}


//------------------------------------------------------------------
// Float BhamiltonNode(Vector *boson, Float mass):
// The boson Hamiltonian of the node sublattice.
//------------------------------------------------------------------
Float FdwfBase::BhamiltonNode(Vector *boson, Float mass){
  char *fname = "BhamiltonNode(V*,F)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = mass;

  if (boson == 0)
    ERR.Pointer(cname,fname,"boson");

  int f_size = GJP.VolNodeSites() * FsiteSize() / 2 ;

  char *str_tmp = "bsn_tmp" ;
  Vector *bsn_tmp = (Vector *)
    smalloc(cname,fname,str_tmp,f_size*sizeof(Float));

  DiracOpDwf dwf(*this, boson, bsn_tmp, &cg_arg, CNV_FRM_NO) ;

  dwf.MatPc(bsn_tmp,boson);

  Float ret_val = bsn_tmp->NormSqNode(f_size) ;

  sfree(cname,fname,str_tmp,bsn_tmp);

  // Sum accross s nodes in case Snodes() != 1
  glb_sum_dir(&ret_val, 4) ;

  return ret_val ;
}


//------------------------------------------------------------------
// int FsiteOffsetChkb(const int *x):
// Sets the offsets for the fermion fields on a 
// checkerboard. The fermion field storage order
// is not the canonical one but it is particular
// to the fermion type. x[i] is the 
// ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
//------------------------------------------------------------------
int FdwfBase::FsiteOffsetChkb(const int *x) const {
// ???
  ERR.NotImplemented(cname, "FsiteOffsetChkb");
  return 0; 
}


//------------------------------------------------------------------
// int FsiteOffset(const int *x):
// Sets the offsets for the fermion fields on a 
// checkerboard. The fermion field storage order
// is the canonical one. X[I] is the
// ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
//------------------------------------------------------------------
int FdwfBase::FsiteOffset(const int *x) const {
// ???
  ERR.NotImplemented(cname, "FsiteOffset");
  return 0; 
}


//--------------------------------------------------------------------
// void Freflex (Vector *out, Vector *in)
// does the reflexion in s needed for the hermitian D_dwf operator.
//
//--------------------------------------------------------------------
void FdwfBase::Freflex(Vector *out, Vector *in)
{
  char *fname = "Freflex(V*,V*)";
  VRB.Func(cname,fname);

  int i,n,node,f_size_5d,f_size_4d,half_size_5d,half_size_4d;
  Vector *send_buf, *rcv_buf;
  int numblk, blklen,s,s_reflex;

  f_size_5d= GJP.VolNodeSites()*FsiteSize();
  f_size_4d= GJP.VolNodeSites()*2*Colors()*SpinComponents();

  half_size_5d= f_size_5d/2;
  half_size_4d= f_size_4d/2;

  numblk=1;
  blklen=f_size_5d;
  while (blklen > 1023) {
    numblk*=2;
    blklen/=2;
  }

  VRB.Debug (cname, fname,"%d %d %d %d %d %d\n",
	     f_size_5d, f_size_4d, half_size_5d, half_size_4d, 
	     numblk, blklen);

  //reserve space for send and receive buffers
  send_buf= (Vector *) smalloc (f_size_5d*sizeof(IFloat));
  if (send_buf == 0) ERR.Pointer (cname, fname, "send_buf");
  VRB.Smalloc (cname, fname, "send_buf", send_buf, f_size_5d * sizeof(IFloat));
  rcv_buf= (Vector *) smalloc (f_size_5d*sizeof(IFloat));
  if (rcv_buf == 0) ERR.Pointer (cname, fname, "rcv_buf");
  VRB.Smalloc (cname, fname, "rcv_buf", rcv_buf, f_size_5d * sizeof(IFloat));

  for (n=0; n<GJP.Snodes(); n++) {

    VRB.Debug (cname,fname,"n= %d\n", n);
    if (n==0) {
      VRB.Debug (cname,fname,"copy from in to rcv_buf\n");
      // copy from in to rcv_buf
      for (i=0; i<f_size_5d; i++)
	((IFloat *)rcv_buf)[i]= ((IFloat *)in)[i];
    }
    else {
      VRB.Debug (cname,fname,
		"copy from rcv_buf to send_buf, send send_buf to rcv_buf\n");
      // copy from rcv_buf to send_buf
      // send send_buf to rcv_buf
      for (i=0; i<f_size_5d; i++)
	((IFloat *)send_buf)[i]= ((IFloat *)rcv_buf)[i];

      for (i=0; i<numblk;i++)
	getMinusData (((IFloat *) rcv_buf)+i*blklen, 
		      ((IFloat *) send_buf)+i*blklen, blklen, 4);
    }

    // node is the node where the data in rcv_buf comes from originally
    if ((node=GJP.SnodeCoor()-n)<0)
      node+= GJP.Snodes();

    VRB.Debug (cname,fname,"node= %d\n", node);

    if (node==GJP.Snodes()-1-GJP.SnodeCoor()) {
      VRB.Debug (cname,fname,"we have the right data in rcv_buf\n");
      // we have the right data in rcv_buf
      // now do the reflexion

      // in and out are in dwf checkerboard storage order

      for (s=0;s<GJP.SnodeSites();s++) {
	s_reflex= GJP.SnodeSites()-1-s;

	if ((GJP.SnodeSites()%2)==0) {
	  for (i=0; i<half_size_4d; i++) {
	    ((IFloat *)out)[s_reflex*half_size_4d+i]= 
	      ((IFloat *)rcv_buf)[s*half_size_4d+i+half_size_5d];
	    ((IFloat *)out)[s_reflex*half_size_4d+i+half_size_5d]= 
	      ((IFloat *)rcv_buf)[s*half_size_4d+i];
	  }
	}
	else {
	  for (i=0; i<half_size_4d; i++) {
	    ((IFloat *)out)[s_reflex*half_size_4d+i]= 
	      ((IFloat *)rcv_buf)[s*half_size_4d+i];
	    ((IFloat *)out)[s_reflex*half_size_4d+i+half_size_5d]= 
             ((IFloat *)rcv_buf)[s*half_size_4d+i+half_size_5d];
	  }
	}
      }
    }
  }
  
  VRB.Sfree(cname,fname, "send_buf", send_buf);
  sfree(send_buf);
  VRB.Sfree(cname,fname, "rcv_buf", rcv_buf);
  sfree(rcv_buf);

  VRB.FuncEnd (cname,fname);
}


CPS_END_NAMESPACE
