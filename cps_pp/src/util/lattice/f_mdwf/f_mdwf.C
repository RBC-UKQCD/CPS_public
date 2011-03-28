#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of Fwilson class.

  $Id: f_mdwf.C,v 1.3 2011-03-28 16:01:11 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_mdwf/f_mdwf.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_wilson.C
//
// Fwilson is derived from FwilsonTypes and is relevant to
// wilson fermions
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/wilson.h>
#include <util/verbose.h>
#include <util/gjp.h>
#include <util/error.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <util/enum_func.h>

#define BENCHMARK
#ifdef BENCHMARK
#include <util/qcdio.h>
#include <sys/time.h>
#ifndef timersub
#define timersub(a, b, result)                                                \
  do {                                                                        \
    (result)->tv_sec = (a)->tv_sec - (b)->tv_sec;                             \
    (result)->tv_usec = (a)->tv_usec - (b)->tv_usec;                          \
    if ((result)->tv_usec < 0) {                                              \
      --(result)->tv_sec;                                                     \
      (result)->tv_usec += 1000000;                                           \
    }                                                                         \
  } while (0)
#endif
#endif

CPS_START_NAMESPACE

//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
Fmdwf::Fmdwf(void)
{
  cname = "Fmdwf";
  const char *fname = "Fmdwf()";
  VRB.Func(cname,fname);

  MdwfArg *_mdwf_arg_p = GJP.GetMdwfArg();

  if(_mdwf_arg_p == NULL){
    ERR.General(cname, fname, "MdwfArg is not initialized.\n");
  }
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
Fmdwf::~Fmdwf(void)
{
  const char *fname = "~Fmdwf()";
  VRB.Func(cname,fname);
}

//! Multiplication of a lattice spin-colour vector by gamma_5.
void Fmdwf::Gamma5(Vector *v_out, Vector *v_in, int num_sites)
{
  Float * p_out = (Float *)v_out;
  Float * p_in  = (Float *)v_in;
  int half_site_size = 12 ;
  for (int site=0; site<num_sites; ++site) {
    int comp;
    for (comp=0; comp<half_site_size; ++comp) {
      *p_out++ = *p_in++ ;
    }
    for (comp=0; comp<half_site_size; ++comp) {
      *p_out++ = -*p_in++ ;
    }
  }
}


// returns the type of fermion class
FclassType Fmdwf::Fclass(void)const
{
  return F_CLASS_DWF;
}

// Sets the offsets for the fermion fields on a 
// checkerboard. The fermion field storage order
// is not the canonical one but it is particular
// to the Dwf fermion type. x[i] is the 
// ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
int Fmdwf::FsiteOffsetChkb(const int *x) const
{
  const char * fname = "FsiteOffsetChkb()";
  ERR.NotImplemented(cname, fname);
}

// Sets the offsets for the fermion fields on a 
// checkerboard. The fermion field storage order
// is the canonical one. X[I] is the
// ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
int Fmdwf::FsiteOffset(const int *x) const
{
  const char * fname = "FsiteOffset()";
  ERR.NotImplemented(cname, fname);
}

// Returns the number of fermion field 
// components (including real/imaginary) on a
// site of the 4-D lattice.
int Fmdwf::FsiteSize(void)const
{
  MdwfArg *_mdwf_arg_p = GJP.GetMdwfArg();
  return 2 * Colors() * SpinComponents() * _mdwf_arg_p->b5.b5_len;  
}

// Returns 0 => If no checkerboard is used for the evolution
//      or the CG that inverts the evolution matrix.
int Fmdwf::FchkbEvl(void)const
{
  return 1;
}

// It calculates f_out where A * f_out = f_in and
// A is the preconditioned fermion matrix that appears
// in the HMC evolution (even/odd preconditioning 
// of [Dirac^dag Dirac]). The inversion is done
// with the conjugate gradient. cg_arg is the structure
// that contains all the control parameters, f_in is the
// fermion field source vector, f_out should be set to be
// the initial guess and on return is the solution.
// f_in and f_out are defined on a checkerboard.
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// The function returns the total number of CG iterations.
int Fmdwf::FmatEvlInv(Vector *f_out, Vector *f_in, 
                      CgArg *cg_arg, 
                      Float *true_res,
                      CnvFrmType cnv_frm)
{
  const char * fname = "FmatEvlInt(V*, V*, CgArg *, ...)";
  ERR.NotImplemented(cname, fname);
  return -1;
}

int Fmdwf::FmatEvlInv(Vector *f_out, Vector *f_in, 
                      CgArg *cg_arg, 
                      CnvFrmType cnv_frm)
{
  return FmatEvlInv(f_out, f_in, cg_arg, NULL, cnv_frm);
}

int Fmdwf::FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
                       int Nshift, int isz, CgArg **cg_arg, 
                       CnvFrmType cnv_frm, MultiShiftSolveType type, Float *alpha,
                       Vector **f_out_d)
{
  const char * fname = "FmatEvlMInv(V*, V*, Float *, ...)";
  ERR.NotImplemented(cname, fname);
  return -1;
}

void Fmdwf::FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
                       Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm)
{
  const char * fname = "FminResExt(V*, V*, V**, ...)";
  ERR.NotImplemented(cname, fname);
}
    
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
int Fmdwf::FmatInv(Vector *f_out, Vector *f_in, 
                   CgArg *cg_arg, 
                   Float *true_res,
                   CnvFrmType cnv_frm, // ignored, not used by DiracOpMdwf class
                   PreserveType prs_f_in)
{
  const char * fname = "FmatInv(V*, V*, CgArg *, ...)";
  VRB.Func(cname, fname);

  MdwfArg *_mdwf_arg_p = GJP.GetMdwfArg();

  CgArg cg_arg_old = _mdwf_arg_p->cg_arg;
  _mdwf_arg_p->cg_arg = *cg_arg;  

  DiracOpMdwf mdwf(*this, _mdwf_arg_p);
  int iter = mdwf.MatInv(f_out, f_in, true_res, prs_f_in);

  _mdwf_arg_p->cg_arg = cg_arg_old;
  return iter;
}

// FmatInvMobius: same as FmatInv, except that we use mobius DWF
// formalism to speed up the CG inversion (via constructing initial guess).
// n_restart: How many restarts we perform
int Fmdwf::FmatInvMobius(Vector *f_out,
                         Vector *f_in,
                         MdwfArg *mob_l,
                         MdwfArg *mob_s,
                         Float *true_res,
                         CnvFrmType cnv_frm,
                         PreserveType prs_f_in,
                         int n_restart, Float rsd_vec[])
{
  const char *fname = "FmatInvMobius(V*, V*, ...)";
  // this implementation doesn't allow splitting in s direction(yet).
  if(GJP.Snodes() != 1){
    ERR.NotImplemented(cname, fname);
  }
  if(cnv_frm != CNV_FRM_YES){
    ERR.NotImplemented(cname, fname);
  }
  if(n_restart < 2){
    ERR.General(cname, fname, "Value %d is invalid for n_restart.\n", n_restart);
  }

  int mob_l_ls = mob_l->b5.b5_len;
  int mob_s_ls = mob_s->b5.b5_len;
  int size_4d = GJP.VolNodeSites() * SPINOR_SIZE;
  int mob_l_size_5d = size_4d * mob_l_ls;
  int mob_s_size_5d = size_4d * mob_s_ls;

  Vector *tmp_mob_l_5d = (Vector *) smalloc(cname, fname, "tmp_mob_l_5d", sizeof(Float) * mob_l_size_5d);
  Vector *tmp2_mob_l_5d = (Vector *) smalloc(cname, fname, "tmp2_mob_l_5d", sizeof(Float) * mob_l_size_5d);
  Vector *tmp_mob_s_5d = (Vector *) smalloc(cname, fname, "tmp_mob_s_5d", sizeof(Float) * mob_s_size_5d);

  CgArg *cg_arg_l = &mob_l->cg_arg;
  CgArg *cg_arg_s = &mob_s->cg_arg;

  // the first time, we solve in DWF to some degree of accuracy
  {
    Float rsd = cg_arg_l->stop_rsd;
    cg_arg_l->stop_rsd = rsd_vec[0];
    
    DiracOpMdwf mdwf(*this, mob_l);
    mdwf.MatInv(f_out, f_in, NULL, PRESERVE_YES);
    cg_arg_l->stop_rsd = rsd;
  }
  
  int restart_cnt;
  for(restart_cnt = 1; restart_cnt < n_restart; ++restart_cnt){
    //constructing the new residue
    tmp2_mob_l_5d->CopyVec(f_in, mob_l_size_5d);
    {
      DiracOpMdwf mdwf(*this, mob_l);
      mdwf.Mat(tmp_mob_l_5d, f_out);
    }
    tmp2_mob_l_5d->VecMinusEquVec(tmp_mob_l_5d, mob_l_size_5d);
    //end constructing the new residue(new residue now in tmp2_mob_l_5d)

    cg_arg_s->stop_rsd = rsd_vec[restart_cnt];

    tmp_mob_l_5d->VecZero(mob_l_size_5d);
    {
      Float mass = cg_arg_l->mass;
      Float stop_rsd = cg_arg_l->stop_rsd;
      int max_num_iter = cg_arg_l->max_num_iter;

      cg_arg_l->mass = 1.0;
      cg_arg_l->stop_rsd = rsd_vec[restart_cnt];
      cg_arg_l->max_num_iter = cg_arg_s->max_num_iter;

      DiracOpMdwf mdwf(*this, mob_l);
      mdwf.MatInv(tmp_mob_l_5d, tmp2_mob_l_5d, NULL, PRESERVE_YES);

      cg_arg_l->mass = mass;
      cg_arg_l->stop_rsd = stop_rsd;
      cg_arg_l->max_num_iter = max_num_iter;
    }
    SpinProject(tmp2_mob_l_5d, tmp_mob_l_5d, mob_l_ls, 1);
    
    // go to the small Mobius lattice
    tmp_mob_s_5d->VecZero(mob_s_size_5d);
    tmp_mob_s_5d->CopyVec(tmp2_mob_l_5d, size_4d);
    
    SpinProject(tmp_mob_l_5d, tmp_mob_s_5d, mob_s_ls, 0);
    
    {
      Float mass= cg_arg_s->mass; 
      cg_arg_s->mass = 1.0;
      
      DiracOpMdwf mdwf(*this, mob_s);
      mdwf.Mat(tmp_mob_s_5d, tmp_mob_l_5d);
      
      cg_arg_s->mass = mass;
    }
    
    tmp_mob_l_5d->VecZero(mob_s_size_5d);
    {
      DiracOpMdwf mdwf(*this, mob_s);
      mdwf.MatInv(tmp_mob_l_5d, tmp_mob_s_5d, NULL, PRESERVE_YES);
    }
    SpinProject(tmp_mob_s_5d, tmp_mob_l_5d, mob_s_ls, 1);
    
    //OV2DWF
    tmp_mob_l_5d->VecNegative(tmp_mob_s_5d, size_4d);
    {
      Vector * tmp_2nd_plane = (Vector *)((Float *)tmp_mob_l_5d + size_4d);
      Vector * src_2nd_plane = (Vector *)((Float *)tmp2_mob_l_5d + size_4d);
      tmp_2nd_plane->CopyVec(src_2nd_plane, mob_l_size_5d - size_4d);
    }
    
    SpinProject(tmp2_mob_l_5d, tmp_mob_l_5d, mob_l_ls, 0);
    
    {
      DiracOpMdwf mdwf(*this, mob_l);
      mdwf.Mat(tmp_mob_l_5d, tmp2_mob_l_5d);
    }
    
    tmp2_mob_l_5d->VecZero(mob_l_size_5d);
    {
      Float mass = cg_arg_l->mass;
      Float stop_rsd = cg_arg_l->stop_rsd;
      int max_num_iter = cg_arg_l->max_num_iter;

      cg_arg_l->mass = 1.0;
      cg_arg_l->stop_rsd = cg_arg_s->stop_rsd;
      cg_arg_l->max_num_iter = cg_arg_s->max_num_iter;

      DiracOpMdwf mdwf(*this, mob_l);
      mdwf.MatInv(tmp2_mob_l_5d, tmp_mob_l_5d, NULL, PRESERVE_NO);

      cg_arg_l->mass = mass;
      cg_arg_l->stop_rsd = stop_rsd;
      cg_arg_l->max_num_iter = max_num_iter;
    }
    SpinProject(tmp_mob_l_5d, tmp2_mob_l_5d, mob_l_ls, 1);
    
    tmp_mob_l_5d->CopyVec(tmp_mob_s_5d, size_4d);
    
    SpinProject(tmp2_mob_l_5d, tmp_mob_l_5d, mob_l_ls, 0);

    f_out->VecAddEquVec(tmp2_mob_l_5d, mob_l_size_5d);
  }

  // final solve
  int iter;
  {
    DiracOpMdwf mdwf(*this, mob_l);
    iter = mdwf.MatInv(f_out, f_in, true_res, prs_f_in);
  }

  sfree(cname, fname,  "tmp_mob_l_5d",  tmp_mob_l_5d);
  sfree(cname, fname, "tmp2_mob_l_5d", tmp2_mob_l_5d);
  sfree(cname, fname,  "tmp_mob_s_5d",  tmp_mob_s_5d);

  return iter;
}

int Fmdwf::FmatInv(Vector *f_out, Vector *f_in, 
            CgArg *cg_arg, 
            CnvFrmType cnv_frm,
            PreserveType prs_f_in)
{
  return FmatInv(f_out, f_in, cg_arg, NULL, cnv_frm, prs_f_in);
}


//!< Transforms a 4-dimensional fermion field into a 5-dimensional field.
/* The 5d field is zero */
// The 5d field is zero
// except for the upper two components (right chirality)
// at s = s_u which are equal to the ones of the 4d field
// and the lower two components (left chirality) 
// at s_l, which are equal to the ones of the 4d field
// For spread-out DWF s_u, s_l refer to the global
// s coordinate i.e. their range is from 
// 0 to [GJP.Snodes() * GJP.SnodeSites() - 1]
void Fmdwf::Ffour2five(Vector *five, Vector *four, int s_u, int s_l, int Ncb)
{
  const char * fname = "Ffour2five(V*, V*, ...)";
  VRB.Func(cname,fname);

  int x;
  int i;
  Float *field_4D;
  Float *field_5D;

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

//!< Transforms a 5-dimensional fermion field into a 4-dimensional field.
//The 4d field has
// the upper two components (right chirality) equal to the
// ones of the 5d field at s = s_u and the lower two 
// components (left chirality) equal to the
// ones of the 5d field at s = s_l, where s is the 
// coordinate in the 5th direction.
// For spread-out DWF s_u, s_l refer to the global
// s coordinate i.e. their range is from 
// 0 to [GJP.Snodes() * GJP.SnodeSites() - 1]
// The same 4D field is generarted in all s node slices.
void Fmdwf::Ffive2four(Vector *four, Vector *five, int s_u, int s_l, int Ncb)
{
  const char * fname = "Ffive2four(V*,V*,i,i)";

  int x;
  int i;
  Float *field_4D;
  Float *field_5D;
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

// It finds the eigenvectors and eigenvalues of A where
// A is the fermion matrix (Dirac operator). The solution
// uses Ritz minimization. eig_arg is the 
// structure that contains all the control parameters, f_eigenv
// are the fermion field source vectors which should be
// defined initially, lambda are the eigenvalues returned 
// on solution. f_eigenv is defined on the whole lattice.
// The function returns the total number of Ritz iterations.
int Fmdwf::FeigSolv(Vector **f_eigenv, Float *lambda,
                    Float *chirality, int *valid_eig,
                    Float **hsum,
                    EigArg *eig_arg, 
                    CnvFrmType cnv_frm)
{
  const char * fname = "FeigSolv(EigArg*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  int iter;

  CgArg cg_arg;
  cg_arg.mass = eig_arg->mass;
  cg_arg.RitzMatOper = eig_arg->RitzMatOper;
  int N_eig = eig_arg->N_eig;

  MdwfArg *_mdwf_arg_p = GJP.GetMdwfArg();
  {
    CgArg cg_arg_old = _mdwf_arg_p->cg_arg;
    _mdwf_arg_p->cg_arg = cg_arg;

    DiracOpMdwf mdwf(*this, _mdwf_arg_p);
  
    iter = mdwf.RitzEig(f_eigenv, lambda, valid_eig, eig_arg);
    _mdwf_arg_p->cg_arg = cg_arg_old;
  }
 
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
    Ffive2four(four, f_eigenv[i], 0, GJP.Snodes()*GJP.SnodeSites()-1,Ncb);

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
 
// It sets the pseudofermion field phi from frm1, frm2.
Float Fmdwf::SetPhi(Vector *phi, Vector *frm1, Vector *frm2,	       
                    Float mass, DagType dag)
{
  const char *fname = "SetPhi()";
  ERR.NotImplemented(cname, fname);
}

	
// It evolves the canonical momentum mom by step_size
// using the fermion force.
ForceArg Fmdwf::EvolveMomFforce(Matrix *mom, Vector *frm, 
                                Float mass, Float step_size)
{
  const char *fname = "EvolveMomFforce()";
  ERR.NotImplemented(cname, fname);
}

// It evolve the canonical momentum mom  by step_size
// using the bosonic quotient force.
ForceArg Fmdwf::EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
                         Float mass, Float step_size)
{
  const char *fname = "EvolveMomFforce()";
  ERR.NotImplemented(cname, fname);
}

ForceArg Fmdwf::RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
			      int isz, Float *alpha, Float mass, Float dt,
			      Vector **sol_d, ForceMeasure measure)
{
  const char *fname = "RHMC_EvolveMomFforce()";
  ERR.NotImplemented(cname, fname);
}

// The fermion Hamiltonian of the node sublattice.
// chi must be the solution of Cg with source phi.
// copied from FdwfBase
Float Fmdwf::FhamiltonNode( Vector *phi,  Vector *chi)
{
  const char *fname = "FhamiltonNode(V*, V*)";
  VRB.Func(cname,fname);

  if (phi == 0) ERR.Pointer(cname,fname,"phi") ;
  if (chi == 0) ERR.Pointer(cname,fname,"chi") ;

  int f_size = GJP.VolNodeSites() * FsiteSize() / 2 ;

  Float ret_val;

  ret_val = phi->ReDotProductNode(chi, f_size );

  // Sum accross s nodes in case Snodes() != 1
  glb_sum_dir(&ret_val, 4);

  return ret_val;
}

// Convert fermion field f_field from -> to
void Fmdwf::Fconvert(Vector *f_field,
              StrOrdType to,
              StrOrdType from)
{
  const char * fname = "Fconvert()";

  // nothing needs to be done
  //ERR.NotImplemented(cname, fname);
}


// The boson Hamiltonian of the node sublattice
Float Fmdwf::BhamiltonNode(Vector *boson, Float mass)
{
  const char * fname = "BhamiltonNode()";
  ERR.NotImplemented(cname, fname);
}

//--------------------------------------------------------------------
// void SpinProject():
//
// Does a spin projection along s direction, specifically (the
// matrices appear below apply in s direction, the example are given
// for Ls=5):
// 
// for type == 0:
//
//       [ P- P+ 0  0  0  ]
//       [ 0  P- P+ 0  0  ]
// out = [ 0  0  P- P+ 0  ] in
//       [ 0  0  0  P- P+ ]
//       [ P+ 0  0  0  P- ]
//
// for type == 1:
//
//       [ P- 0  0  0  P+ ]
//       [ P+ P- 0  0  0  ]
// out = [ 0  P+ P- 0  0  ] in
//       [ 0  0  P+ P- 0  ]
//       [ 0  0  0  P+ P- ]
//
// in and out are assumed to be in CANONICAL storage order.
//--------------------------------------------------------------------
void Fmdwf::SpinProject(Vector * out, Vector *in, int s_size, int type)
{
  const char * fname = "SpinProject(V*, V*, int)";

  // this implementation doesn't allow splitting in s direction(yet).
  if(GJP.Snodes() != 1){
    ERR.NotImplemented(cname, fname);
  }

  int lcl_size[5];

  lcl_size[0] = GJP.XnodeSites();
  lcl_size[1] = GJP.YnodeSites();
  lcl_size[2] = GJP.ZnodeSites();
  lcl_size[3] = GJP.TnodeSites();
  lcl_size[4] = s_size;

#define for_each_site_4d(x)                             \
  for(x[0] = 0; x[0] < lcl_size[0]; ++x[0])             \
    for(x[1] = 0; x[1] < lcl_size[1]; ++x[1])           \
      for(x[2] = 0; x[2] < lcl_size[2]; ++x[2])         \
        for(x[3] = 0; x[3] < lcl_size[3]; ++x[3])       

  int lcl_pos[4], color;
  int s;

  int f_size_4d= GJP.VolNodeSites()*2*Colors()*SpinComponents();

  for(s = 0; s < lcl_size[4]; ++s){
    for_each_site_4d(lcl_pos){
      int offset = s*f_size_4d;

      int offset_next;
      if(type == 0){
        offset_next = s == lcl_size[4]-1 ? 0: offset+f_size_4d;
      }else{
        offset_next = s == 0 ? (lcl_size[4]-1)*f_size_4d: offset-f_size_4d;
      }
      
      offset += SPINOR_SIZE * (lcl_pos[0] + lcl_size[0]*(lcl_pos[1] + lcl_size[1]*(lcl_pos[2] + lcl_size[2]*lcl_pos[3])));
      offset_next += SPINOR_SIZE * (lcl_pos[0] + lcl_size[0]*(lcl_pos[1] + lcl_size[1]*(lcl_pos[2] + lcl_size[2]*lcl_pos[3])));
      
      Float * out_site = (Float *)out + offset;
      Float * in_site = (Float *)in + offset;
      Float * in_site_next = (Float *)in + offset_next;
      
      memcpy((void *)out_site,
             (void *)in_site_next,
             SPINOR_SIZE/2*sizeof(Float));
      
      memcpy((void *)(out_site + SPINOR_SIZE/2),
             (void *)(in_site  + SPINOR_SIZE/2),
             SPINOR_SIZE/2*sizeof(Float));
    }
  }
  return;
}

// Reflexion in s operator, needed for the hermitian version 
// of the dirac operator in the Ritz solver.
void Fmdwf::Freflex(Vector *out, Vector *in)
{
  const char * fname = "Freflex(V*,V*)";
  VRB.Func(cname,fname);

  int i;
  MdwfArg *_mdwf_arg_p = GJP.GetMdwfArg();
  int ls = _mdwf_arg_p->b5.b5_len;
  // for MDWF, GJP.Snode() must be 1
  int slice_len = SPINOR_SIZE * GJP.VolNodeSites();
  for(i = 0; i < ls; ++i){
    Vector * out_slice = (Vector *)(((Float *)out) + slice_len * i);
    Vector * in_slice =  (Vector *)(((Float *)in) + slice_len * (ls - 1 - i));
    out_slice->CopyVec(in_slice, slice_len);
  }

  VRB.FuncEnd (cname,fname);
}

int Fmdwf::SpinComponents()const
{
  return 4;
}

int Fmdwf::ExactFlavors()const
{
  return 2;
}
    
//!< Method to ensure bosonic force works (does nothing for Wilson
//!< theories.
void Fmdwf::BforceVector(Vector *in, CgArg *cg_arg)
{
  return;
}

CPS_END_NAMESPACE
