#include <config.h>
#include <stdio.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOpMdwf class methods.

  $Id: d_op_mdwf.C,v 1.6 2013-04-05 17:51:14 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:51:14 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mdwf/d_op_mdwf.C,v 1.6 2013-04-05 17:51:14 chulwoo Exp $
//  $Id: d_op_mdwf.C,v 1.6 2013-04-05 17:51:14 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: d_op_mdwf.C,v $
//  $Revision: 1.6 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mdwf/d_op_mdwf.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// d_op_mdwf.C
//
// DiracOpMdwf is the front end for a library that contains
// all Dirac operators associated with Mobius Dwf fermions.
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
//#include <mem/p2v.h>
#include <comms/glb.h>
CPS_START_NAMESPACE
#ifdef USE_MDWF
#ifdef __cplusplus
extern "C"
{
#endif
#include <qop-mdwf3.h>
#ifdef __cplusplus
}
#endif
#endif

#ifdef USE_MDWF
// helper function used by QOP_MDWF_init only, should not be used
// otherwise.
static void sublattice(int lo[], int hi[], const int node[],
                       void *env) // env is currently ignored.
{
  const char *fname = "sublattice()";
  int xsize = GJP.XnodeSites();
  int ysize = GJP.YnodeSites();
  int zsize = GJP.ZnodeSites();
  int tsize = GJP.TnodeSites();
  lo[0] = node[0] * xsize;
  lo[1] = node[1] * ysize;
  lo[2] = node[2] * zsize;
  lo[3] = node[3] * tsize;
  hi[0] = lo[0] + xsize;
  hi[1] = lo[1] + ysize;
  hi[2] = lo[2] + zsize;
  hi[3] = lo[3] + tsize;
}

// helper function used to export gauge field to the mdwf library,
// should not be used elsewhere.
// env is used as the pointer to the instance of the lattice class.
static double read_gauge(int dir, const int pos[4], int a, int b, int re_im, void *env)
{
  const char *fname = "read_gauge()";
  Lattice *latt = (Lattice *)env;

  if(latt->StrOrd() != CANONICAL){
    ERR.General("",fname, "Storage order of the gauge field is not acceptable.\n");
  }


  int my_pos[QOP_MDWF_DIM];
  my_pos[0] = pos[0];
  my_pos[1] = pos[1];
  my_pos[2] = pos[2];
  my_pos[3] = pos[3];

  my_pos[0] %= GJP.XnodeSites();
  my_pos[1] %= GJP.YnodeSites();
  my_pos[2] %= GJP.ZnodeSites();
  my_pos[3] %= GJP.TnodeSites();

  Matrix *gauge_field = latt->GaugeField() + latt->GsiteOffset(my_pos) + dir;
  return ( re_im == 0 ) ? (*gauge_field)[a*3+b].real() : (*gauge_field)[a*3+b].imag();
}

// helper functions used to import/export fermion field from/to the mdwf library,
// should not be used elsewhere.
// env is used as the pointer to the Vector class that provides the fermion field,
// corrdinate order: pos[0] is X, pos[3] is T, pos[4] is S.
static double read_fermion_canonical(const int pos[5], int color, int dirac, int re_im, void *env)
{
  const char *fname = "read_fermion_canonical()";
  Float *fermion_field = (Float *)env;
  if(fermion_field == NULL){
    ERR.Pointer("", fname, "fermion_field");
  }
  re_im = (re_im == 0) ? 0 : 1;

  int lcl_pos[QOP_MDWF_DIM + 1];
  int lcl_size[QOP_MDWF_DIM];

  lcl_size[0] = GJP.XnodeSites();
  lcl_size[1] = GJP.YnodeSites();
  lcl_size[2] = GJP.ZnodeSites();
  lcl_size[3] = GJP.TnodeSites();

  lcl_pos[0] = pos[0] % lcl_size[0];
  lcl_pos[1] = pos[1] % lcl_size[1];
  lcl_pos[2] = pos[2] % lcl_size[2];
  lcl_pos[3] = pos[3] % lcl_size[3];
  lcl_pos[4] = pos[4]; // only one node in S dir

  int offset = re_im + 2*(color + 3*dirac)
    + SPINOR_SIZE * (lcl_pos[0] + lcl_size[0]*(lcl_pos[1] + lcl_size[1]*(lcl_pos[2] + lcl_size[2]*(lcl_pos[3] + lcl_size[3]*lcl_pos[4]))));
  return fermion_field[offset];
}

static void write_fermion_canonical(const int pos[5], int color, int dirac, int re_im, double value, void *env)
{
  const char *fname = "write_fermion_canonical()";
  Float *fermion_field = (Float *)env;
  if(fermion_field == NULL){
    ERR.Pointer("", fname, "fermion_field");
  }
  re_im = (re_im == 0) ? 0 : 1;
  if(color >= 3 || color < 0){
    ERR.General("", fname, "Unexpected color value: %d\n.", color);
  }
  if(dirac >= 4 || dirac < 0){
    ERR.General("", fname, "Unexpected spin component value: %d\n", dirac);
  }

  int lcl_pos[QOP_MDWF_DIM + 1];
  int lcl_size[QOP_MDWF_DIM];

  lcl_size[0] = GJP.XnodeSites();
  lcl_size[1] = GJP.YnodeSites();
  lcl_size[2] = GJP.ZnodeSites();
  lcl_size[3] = GJP.TnodeSites();

  lcl_pos[0] = pos[0] % lcl_size[0];
  lcl_pos[1] = pos[1] % lcl_size[1];
  lcl_pos[2] = pos[2] % lcl_size[2];
  lcl_pos[3] = pos[3] % lcl_size[3];
  lcl_pos[4] = pos[4]; // only one node in S dir

  int offset = re_im + 2*(color + 3*dirac)
    + SPINOR_SIZE * (lcl_pos[0] + lcl_size[0]*(lcl_pos[1] + lcl_size[1]*(lcl_pos[2] + lcl_size[2]*(lcl_pos[3] + lcl_size[3]*lcl_pos[4]))));
  fermion_field[offset] = value;
}

// same helper functions as above, but the source (fermion field stored
// in CPS) field is 4-D even-odd preconditioned.
static double read_fermion_precond(const int pos[5], int color, int dirac, int re_im, void *env)
{
  const char *fname = "read_fermion_precond()";
  Float *fermion_field = (Float *)env;
  if(fermion_field == NULL){
    ERR.Pointer("", fname, "fermion_field");
  }
  re_im = (re_im == 0) ? 0 : 1;

  int lcl_pos[QOP_MDWF_DIM + 1];
  int lcl_size[QOP_MDWF_DIM];

  lcl_size[0] = GJP.XnodeSites();
  lcl_size[1] = GJP.YnodeSites();
  lcl_size[2] = GJP.ZnodeSites();
  lcl_size[3] = GJP.TnodeSites();

  lcl_pos[0] = pos[0] % lcl_size[0];
  lcl_pos[1] = pos[1] % lcl_size[1];
  lcl_pos[2] = pos[2] % lcl_size[2];
  lcl_pos[3] = pos[3] % lcl_size[3];
  lcl_pos[4] = pos[4]; // only one node in S dir

  int offset = re_im + 2*(color + 3*dirac)
    + SPINOR_SIZE * (lcl_pos[0]/2 + lcl_size[0]/2*(lcl_pos[1] + lcl_size[1]*(lcl_pos[2] + lcl_size[2]*(lcl_pos[3] + lcl_size[3]*lcl_pos[4]))));
  return fermion_field[offset];
}

static void write_fermion_precond(const int pos[5], int color, int dirac, int re_im, double value, void *env)
{
  const char *fname = "write_fermion_precond()";
  Float *fermion_field = (Float *)env;
  if(fermion_field == NULL){
    ERR.Pointer("", fname, "fermion_field");
  }
  re_im = (re_im == 0) ? 0 : 1;

  int lcl_pos[QOP_MDWF_DIM + 1];
  int lcl_size[QOP_MDWF_DIM];

  lcl_size[0] = GJP.XnodeSites();
  lcl_size[1] = GJP.YnodeSites();
  lcl_size[2] = GJP.ZnodeSites();
  lcl_size[3] = GJP.TnodeSites();

  lcl_pos[0] = pos[0] % lcl_size[0];
  lcl_pos[1] = pos[1] % lcl_size[1];
  lcl_pos[2] = pos[2] % lcl_size[2];
  lcl_pos[3] = pos[3] % lcl_size[3];
  lcl_pos[4] = pos[4]; // only one node in S dir

  int offset = re_im + 2*(color + 3*dirac)
    + SPINOR_SIZE * (lcl_pos[0]/2 + lcl_size[0]/2*(lcl_pos[1] + lcl_size[1]*(lcl_pos[2] + lcl_size[2]*(lcl_pos[3] + lcl_size[3]*lcl_pos[4]))));
  fermion_field[offset] = value;
}
#endif

DiracOpMdwf::DiracOpMdwf(Lattice& latt,            // Lattice object.
                         MdwfArg *mdwf_arg_p): DiracOpWilsonTypes(latt,
                                                                   NULL,
                                                                   NULL,
                                                                   &mdwf_arg_p->cg_arg,
                                                                   CNV_FRM_NO)

{
  cname = "DiracOpMdwf";
  const char *fname = "DiracOpMdwf(L&, MdwfArg *)";
  VRB.Func(cname, fname);

#ifdef USE_MDWF
  if( GJP.Snodes() != 1 ) {
    // mdwf library doesn't allow splitting in Ls direction.
    // is this the case?
    ERR.NotImplemented(cname, fname);
  }

  latt_ptr = &latt;

  CgArg *cg_arg_p = &mdwf_arg_p->cg_arg;
  if(cg_arg_p == NULL){
    ERR.General(cname, fname, "Invalid CgArg pointer in MdwfArg.\n");
  }
  epsilon = (cg_arg_p -> stop_rsd) * (cg_arg_p -> stop_rsd);
  max_num_iter = cg_arg_p -> max_num_iter;

  if((mdwf_arg_p->c5).c5_len != (mdwf_arg_p->b5).b5_len ) {
    ERR.General(cname, fname, "in MdwfArg: b5_len doesn't match c5_len.\n");
  }

  use_single_precision = mdwf_arg_p -> use_single_precision;

  // initialize the library
  int lat_size[QOP_MDWF_DIM + 1];
  int network[QOP_MDWF_DIM];
  int node[QOP_MDWF_DIM];
  
  network[0] = GJP.Xnodes();
  network[1] = GJP.Ynodes();
  network[2] = GJP.Znodes();
  network[3] = GJP.Tnodes();
  
  lat_size[0] = network[0] * GJP.XnodeSites();
  lat_size[1] = network[1] * GJP.YnodeSites();
  lat_size[2] = network[2] * GJP.ZnodeSites();
  lat_size[3] = network[3] * GJP.TnodeSites();
  // note: We know that GJP.Snodes() must return 1.
  lat_size[4] = (mdwf_arg_p->b5).b5_len;
  
  node[0] = GJP.XnodeCoor();
  node[1] = GJP.YnodeCoor();
  node[2] = GJP.ZnodeCoor();
  node[3] = GJP.TnodeCoor();
  
  int master_p = (UniqueID() == 0) ? 1 : 0;
  
  QOP_MDWF_init(&mdwf_state, lat_size, network, node, master_p, sublattice, NULL);
  QOP_MDWF_set_generic(&mdwf_param, mdwf_state, (mdwf_arg_p->b5).b5_val, (mdwf_arg_p->c5).c5_val, -mdwf_arg_p->M5, cg_arg_p -> mass);
  VRB.Result(cname,fname,"mass(MDWF)=%g\n", cg_arg_p->mass);

  // import gauge field
  if(use_single_precision) {
    QOP_F3_MDWF_import_gauge(&mdwf_gauge_ptr.f, mdwf_state, read_gauge, latt_ptr);
  } else {
    QOP_D3_MDWF_import_gauge(&mdwf_gauge_ptr.d, mdwf_state, read_gauge, latt_ptr);
  }
#else
  ERR.NotImplemented(cname, fname);
#endif
}

DiracOpMdwf::~DiracOpMdwf()
{
  const char *fname = "~DiracOpMdwf()";
  VRB.Func(cname, fname);

#ifdef USE_MDWF
  if(use_single_precision){
    QOP_F3_MDWF_free_gauge(&mdwf_gauge_ptr.f);
  }else{
    QOP_D3_MDWF_free_gauge(&mdwf_gauge_ptr.d);
  }
  QOP_MDWF_free_parameters(&mdwf_param);
  if(mdwf_param != NULL){ // an error occured
    ERR.General(cname, fname, "An error occured while freeing MDWF parameters.\n");
  }

  QOP_MDWF_fini(&mdwf_state);
  if(mdwf_state != NULL){ //an error occured
    ERR.General(cname, fname, "An error occured while exiting MDWF library.\n");
  }
#else
  ERR.NotImplemented(cname, fname);
#endif
}

// It sets the dirac_arg pointer to arg and initializes
// the relevant parameters (kappa, m^2, ...).
void DiracOpMdwf::DiracArg(CgArg *arg)
{
  const char *fname = "DiracArg()";
  ERR.NotImplemented(cname, fname);
  
  //nothing to do here(we don't use this function in MDWF)
}

// MatPcDagMatPc is the fermion matrix that appears in the HMC 
// evolution. It is a Hermitian matrix.
// The in, out fields are defined on the checkerboard lattice.
// If dot_prd is not 0 then the dot product (on node)
// <out, in> = <MatPcDagMatPc*in, in> is returned in dot_prd.
void DiracOpMdwf::MatPcDagMatPc(Vector *out, Vector *in, Float *dot_prd)
{
  const char *fname = "MatPcDagMatPc()";
  VRB.Func(cname, fname);

#ifdef USE_MDWF
  if(use_single_precision){
    struct QOP_F3_MDWF_HalfFermion *hfermion_ptr_in = NULL;
    struct QOP_F3_MDWF_HalfFermion *hfermion_ptr_out = NULL;
    QOP_F3_MDWF_import_half_fermion(&hfermion_ptr_in, mdwf_state, read_fermion_precond, (void *)in);
    QOP_F3_MDWF_allocate_half_fermion(&hfermion_ptr_out, mdwf_state);
    if(hfermion_ptr_in == NULL || hfermion_ptr_out == NULL){
      ERR.General(cname, fname, "mdwf library initialization failed.\n");
    }

    // apply M
    QOP_F3_MDWF_M_operator(hfermion_ptr_out, mdwf_param, mdwf_gauge_ptr.f, hfermion_ptr_in);
    // calculate the inner product
    if(dot_prd != NULL){
      QOP_F3_MDWF_norm2_half_fermion(dot_prd, hfermion_ptr_out);
    }
    // apply M^\dag
    QOP_F3_MDWF_M_operator_conjugated(hfermion_ptr_in, mdwf_param, mdwf_gauge_ptr.f, hfermion_ptr_out);

    QOP_F3_MDWF_export_half_fermion(write_fermion_precond, (void *)out, hfermion_ptr_in);

    QOP_F3_MDWF_free_half_fermion(&hfermion_ptr_in);
    QOP_F3_MDWF_free_half_fermion(&hfermion_ptr_out);
  }else{
    struct QOP_D3_MDWF_HalfFermion *hfermion_ptr_in = NULL;
    struct QOP_D3_MDWF_HalfFermion *hfermion_ptr_out = NULL;
    QOP_D3_MDWF_import_half_fermion(&hfermion_ptr_in, mdwf_state, read_fermion_precond, (void *)in);
    QOP_D3_MDWF_allocate_half_fermion(&hfermion_ptr_out, mdwf_state);
    if(hfermion_ptr_in == NULL || hfermion_ptr_out == NULL){
      ERR.General(cname, fname, "mdwf library initialization failed.\n");
    }

    // apply M
    QOP_D3_MDWF_M_operator(hfermion_ptr_out, mdwf_param, mdwf_gauge_ptr.d, hfermion_ptr_in);
    // calculate the inner product
    if(dot_prd != NULL){
      QOP_D3_MDWF_norm2_half_fermion(dot_prd, hfermion_ptr_out);
    }
    // apply M^\dag
    QOP_D3_MDWF_M_operator_conjugated(hfermion_ptr_in, mdwf_param, mdwf_gauge_ptr.d, hfermion_ptr_out);

    QOP_D3_MDWF_export_half_fermion(write_fermion_precond, (void *)out, hfermion_ptr_in);

    QOP_D3_MDWF_free_half_fermion(&hfermion_ptr_in);
    QOP_D3_MDWF_free_half_fermion(&hfermion_ptr_out);
  }
#else
  ERR.NotImplemented(cname, fname);
#endif
}

// Dslash is the derivative part of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
// cb = 0/1 <--> even/odd checkerboard of in field.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
void DiracOpMdwf::Dslash(Vector *out, 
                         Vector *in,
                         ChkbType cb, 
                         DagType dag)
{
  const char *fname = "Dslash(V*, V*, ...)";
  VRB.Func(cname, fname);
  ERR.NotImplemented(cname, fname);
}

//! Multiplication by the odd-even preconditioned fermion matrix.
// MatPc is the fermion matrix.  
// The in, out fields are defined on the checkerboard lattice.
void DiracOpMdwf::MatPc(Vector *out, Vector *in)
{
  const char *fname = "MatPc()";
  VRB.Func(cname, fname);

#ifdef USE_MDWF
  if(use_single_precision){
    struct QOP_F3_MDWF_HalfFermion *hfermion_ptr_in = NULL;
    struct QOP_F3_MDWF_HalfFermion *hfermion_ptr_out = NULL;
    QOP_F3_MDWF_import_half_fermion(&hfermion_ptr_in, mdwf_state, read_fermion_precond, (void *)in);
    QOP_F3_MDWF_allocate_half_fermion(&hfermion_ptr_out, mdwf_state);
    if(hfermion_ptr_in == NULL || hfermion_ptr_out == NULL){
      ERR.General(cname, fname, "mdwf library initialization failed.\n");
    }

    // apply M
    QOP_F3_MDWF_M_operator(hfermion_ptr_out, mdwf_param, mdwf_gauge_ptr.f, hfermion_ptr_in);

    QOP_F3_MDWF_export_half_fermion(write_fermion_precond, (void *)out, hfermion_ptr_out);

    QOP_F3_MDWF_free_half_fermion(&hfermion_ptr_in);
    QOP_F3_MDWF_free_half_fermion(&hfermion_ptr_out);
  }else{
    struct QOP_D3_MDWF_HalfFermion *hfermion_ptr_in = NULL;
    struct QOP_D3_MDWF_HalfFermion *hfermion_ptr_out = NULL;
    QOP_D3_MDWF_import_half_fermion(&hfermion_ptr_in, mdwf_state, read_fermion_precond, (void *)in);
    QOP_D3_MDWF_allocate_half_fermion(&hfermion_ptr_out, mdwf_state);
    if(hfermion_ptr_in == NULL || hfermion_ptr_out == NULL){
      ERR.General(cname, fname, "mdwf library initialization failed.\n");
    }

    // apply M
    QOP_D3_MDWF_M_operator(hfermion_ptr_out, mdwf_param, mdwf_gauge_ptr.d, hfermion_ptr_in);

    QOP_D3_MDWF_export_half_fermion(write_fermion_precond, (void *)out, hfermion_ptr_out);

    QOP_D3_MDWF_free_half_fermion(&hfermion_ptr_in);
    QOP_D3_MDWF_free_half_fermion(&hfermion_ptr_out);
  }
#else
  ERR.NotImplemented(cname, fname);
#endif
}

//! Multiplication by the  hermitian conjugate odd-even preconditioned fermion matrix.
// MatPcDag is the dagger of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice.
void DiracOpMdwf::MatPcDag(Vector *out, Vector *in)
{
  const char *fname = "MatPcDag()";
  VRB.Func(cname, fname);

#ifdef USE_MDWF
  if(use_single_precision){
    struct QOP_F3_MDWF_HalfFermion *hfermion_ptr_in = NULL;
    struct QOP_F3_MDWF_HalfFermion *hfermion_ptr_out = NULL;
    QOP_F3_MDWF_import_half_fermion(&hfermion_ptr_in, mdwf_state, read_fermion_precond, (void *)in);
    QOP_F3_MDWF_allocate_half_fermion(&hfermion_ptr_out, mdwf_state);
    if(hfermion_ptr_in == NULL || hfermion_ptr_out == NULL){
      ERR.General(cname, fname, "mdwf library initialization failed.\n");
    }

    // apply M^\dag
    QOP_F3_MDWF_M_operator_conjugated(hfermion_ptr_out, mdwf_param, mdwf_gauge_ptr.f, hfermion_ptr_in);
    
    QOP_F3_MDWF_export_half_fermion(write_fermion_precond, (void *)out, hfermion_ptr_out);
    
    QOP_F3_MDWF_free_half_fermion(&hfermion_ptr_in);
    QOP_F3_MDWF_free_half_fermion(&hfermion_ptr_out);
  }else{
    struct QOP_D3_MDWF_HalfFermion *hfermion_ptr_in = NULL;
    struct QOP_D3_MDWF_HalfFermion *hfermion_ptr_out = NULL;
    QOP_D3_MDWF_import_half_fermion(&hfermion_ptr_in, mdwf_state, read_fermion_precond, (void *)in);
    QOP_D3_MDWF_allocate_half_fermion(&hfermion_ptr_out, mdwf_state);
    if(hfermion_ptr_in == NULL || hfermion_ptr_out == NULL){
      ERR.General(cname, fname, "mdwf library initialization failed.\n");
    }

    // apply M^\dag
    QOP_D3_MDWF_M_operator_conjugated(hfermion_ptr_out, mdwf_param, mdwf_gauge_ptr.d, hfermion_ptr_in);
    
    QOP_D3_MDWF_export_half_fermion(write_fermion_precond, (void *)out, hfermion_ptr_out);
    
    QOP_D3_MDWF_free_half_fermion(&hfermion_ptr_in);
    QOP_D3_MDWF_free_half_fermion(&hfermion_ptr_out);
  }
#else
  ERR.NotImplemented(cname, fname);
#endif
}

// The inverse of the unconditioned Dirac Operator 
// using Conjugate gradient.  source is *in, initial
// guess and solution is *out.
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// prs_in is used to specify if the source
// in should be preserved or not. If not the memory usage
// is less by half the size of a fermion vector.
// The function returns the total number of CG iterations.
int DiracOpMdwf::MatInv(Vector *out, 
           Vector *in, 
           Float *true_res,
           PreserveType prs_in)
{
  const char *fname = "MatInv()";

#ifdef USE_MDWF

//  printf("DiracOpMdwf::MatInv Profiled");
  struct timeval start;
  struct timeval end;
  CGflops    = 0;
  gettimeofday(&start,NULL);
  double mdwf_time = -dclock();


  int n_iter;
  Float true_res_dummy;
  Float *true_res_ptr = (true_res != NULL) ? true_res : &true_res_dummy;
    
  if(use_single_precision){
    struct QOP_F3_MDWF_Fermion *fermion_ptr_in = NULL;
    struct QOP_F3_MDWF_Fermion *fermion_ptr_out = NULL;
    struct QOP_F3_MDWF_Fermion *fermion_ptr_tmp = NULL;

    QOP_F3_MDWF_import_fermion(&fermion_ptr_in, mdwf_state, read_fermion_canonical, (void *)in);
    QOP_F3_MDWF_allocate_fermion(&fermion_ptr_out, mdwf_state);
    QOP_F3_MDWF_import_fermion(&fermion_ptr_tmp, mdwf_state, read_fermion_canonical, (void *)out);

    if(fermion_ptr_in == NULL || fermion_ptr_out == NULL){
      ERR.General(cname, fname, "mdwf library initialization failed.\n");
    }

    // is it safe to set the result and initial guess fermion fields to
    // the same one in the following function?
    QOP_F3_MDWF_DDW_CG(fermion_ptr_out, &n_iter, true_res_ptr, mdwf_param,
                       fermion_ptr_tmp, mdwf_gauge_ptr.f, fermion_ptr_in,
                       max_num_iter, epsilon,
                       QOP_MDWF_FINAL_CG_RESIDUAL);
    
    QOP_F3_MDWF_export_fermion(write_fermion_canonical, (void *)out, fermion_ptr_out);
    
    QOP_F3_MDWF_free_fermion(&fermion_ptr_in);
    QOP_F3_MDWF_free_fermion(&fermion_ptr_out);
    QOP_F3_MDWF_free_fermion(&fermion_ptr_tmp);
  }else{
    struct QOP_D3_MDWF_Fermion *fermion_ptr_in = NULL;
    struct QOP_D3_MDWF_Fermion *fermion_ptr_out = NULL;
    struct QOP_D3_MDWF_Fermion *fermion_ptr_tmp = NULL;

    QOP_D3_MDWF_import_fermion(&fermion_ptr_in, mdwf_state, read_fermion_canonical, (void *)in);
    QOP_D3_MDWF_allocate_fermion(&fermion_ptr_out, mdwf_state);
    QOP_D3_MDWF_import_fermion(&fermion_ptr_tmp, mdwf_state, read_fermion_canonical, (void *)out);

    if(fermion_ptr_in == NULL || fermion_ptr_out == NULL){
      ERR.General(cname, fname, "mdwf library initialization failed.\n");
    }

    // is it safe to set the result and initial guess fermion fields to
    // the same one in the following function?
    QOP_D3_MDWF_DDW_CG(fermion_ptr_out, &n_iter, true_res_ptr, mdwf_param,
                       fermion_ptr_tmp, mdwf_gauge_ptr.d, fermion_ptr_in,
                       max_num_iter, epsilon,
                       QOP_MDWF_FINAL_CG_RESIDUAL);
    
    QOP_D3_MDWF_export_fermion(write_fermion_canonical, (void *)out, fermion_ptr_out);
    
    QOP_D3_MDWF_free_fermion(&fermion_ptr_in);
    QOP_D3_MDWF_free_fermion(&fermion_ptr_out);
    QOP_D3_MDWF_free_fermion(&fermion_ptr_tmp);
  }

  mdwf_time +=dclock();
  print_time(cname,fname,mdwf_time);
    
  return n_iter;
#else
  ERR.NotImplemented(cname, fname);
  return -1;
#endif
}

// Same as original but true_res=0.
int DiracOpMdwf::MatInv(Vector *out, 
                        Vector *in,
                        PreserveType prs_in)
{
  const char *fname = "MatInv()";
  Float true_res;
  return MatInv(out, in, &true_res, prs_in);
}

// Same as original but in = f_in and out = f_out.
int DiracOpMdwf::MatInv(Float *true_res,
                        PreserveType prs_in)
{
  const char *fname = "MatInv()";
  ERR.NotImplemented(cname, fname);
}


// Same as original but in = f_in, out = f_out, true_res=0.
int DiracOpMdwf::MatInv(PreserveType prs_in)
{
  const char *fname = "MatInv()";
  ERR.NotImplemented(cname, fname);
}
 
// MatHerm is the hermitian version of Mat.
// MatHerm works on the full lattice.
// The in, out fields are defined on the ful.
void DiracOpMdwf::MatHerm(Vector *out, Vector *in)
{
  const char *fname = "MatHerm(V*, V*)";
  VRB.Func(cname, fname);

  //----------------------------------------------------------------
  // Implement routine
  // Warning: the following doesn't work since GJP.SnodeSites() can
  // be meaningless for DiracOpMdwf.
  //----------------------------------------------------------------
  ERR.NotImplemented(cname, fname);
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize();
  Vector *temp = (Vector *) smalloc(temp_size * sizeof(Float));

  Mat(out, in);

  lat.Freflex(temp, out);
  MultGamma(out, temp, 15, GJP.VolNodeSites()*GJP.SnodeSites() );

  VRB.Sfree(cname, fname, "temp", temp);
  sfree(temp);
}

// Mat is the unpreconditioned fermion matrix.  
// Mat works on the full lattice
// The in, out fields are defined on the full lattice.
void DiracOpMdwf::Mat(Vector *out, Vector *in)
{
  const char *fname = "Mat(V*, V*)";
  VRB.Func(cname, fname);

#ifdef USE_MDWF
  if(use_single_precision){
    struct QOP_F3_MDWF_Fermion *fermion_ptr_in = NULL;
    struct QOP_F3_MDWF_Fermion *fermion_ptr_out = NULL;
    QOP_F3_MDWF_import_fermion(&fermion_ptr_in, mdwf_state, read_fermion_canonical, (void *)in);
    QOP_F3_MDWF_allocate_fermion(&fermion_ptr_out, mdwf_state);
    if(fermion_ptr_in == NULL || fermion_ptr_out == NULL){
      ERR.General(cname, fname, "mdwf library initialization failed.\n");
    }
    
    // apply D
    QOP_F3_MDWF_DDW_operator(fermion_ptr_out, mdwf_param, mdwf_gauge_ptr.f, fermion_ptr_in);
    
    QOP_F3_MDWF_export_fermion(write_fermion_canonical, (void *)out, fermion_ptr_out);
    
    QOP_F3_MDWF_free_fermion(&fermion_ptr_in);
    QOP_F3_MDWF_free_fermion(&fermion_ptr_out);
  }else{
    struct QOP_D3_MDWF_Fermion *fermion_ptr_in = NULL;
    struct QOP_D3_MDWF_Fermion *fermion_ptr_out = NULL;
    QOP_D3_MDWF_import_fermion(&fermion_ptr_in, mdwf_state, read_fermion_canonical, (void *)in);
    QOP_D3_MDWF_allocate_fermion(&fermion_ptr_out, mdwf_state);
    if(fermion_ptr_in == NULL || fermion_ptr_out == NULL){
      ERR.General(cname, fname, "mdwf library initialization failed.\n");
    }
    
    // apply D
    QOP_D3_MDWF_DDW_operator(fermion_ptr_out, mdwf_param, mdwf_gauge_ptr.d, fermion_ptr_in);
    
    QOP_D3_MDWF_export_fermion(write_fermion_canonical, (void *)out, fermion_ptr_out);
    
    QOP_D3_MDWF_free_fermion(&fermion_ptr_in);
    QOP_D3_MDWF_free_fermion(&fermion_ptr_out);
  }
#else
  ERR.NotImplemented(cname, fname);
#endif
}

// MatDag is the dagger of the unpreconditioned fermion matrix. 
// MatDag works on the full lattice
// The in, out fields are defined on the full lattice.
void DiracOpMdwf::MatDag(Vector *out, Vector *in)
{
  const char *fname = "MatDag(V*, V*)";
  VRB.Func(cname, fname);

#ifdef USE_MDWF
  if(use_single_precision){
    struct QOP_F3_MDWF_Fermion *fermion_ptr_in = NULL;
    struct QOP_F3_MDWF_Fermion *fermion_ptr_out = NULL;
    QOP_F3_MDWF_import_fermion(&fermion_ptr_in, mdwf_state, read_fermion_canonical, (void *)in);
    QOP_F3_MDWF_allocate_fermion(&fermion_ptr_out, mdwf_state);
    if(fermion_ptr_in == NULL || fermion_ptr_out == NULL){
      ERR.General(cname, fname, "mdwf library initialization failed.\n");
    }
    
    // apply D^\dag
    QOP_F3_MDWF_DDW_operator_conjugated(fermion_ptr_out, mdwf_param, mdwf_gauge_ptr.f, fermion_ptr_in);
    
    QOP_F3_MDWF_export_fermion(write_fermion_canonical, (void *)out, fermion_ptr_out);
    
    QOP_F3_MDWF_free_fermion(&fermion_ptr_in);
    QOP_F3_MDWF_free_fermion(&fermion_ptr_out);
  }else{
    struct QOP_D3_MDWF_Fermion *fermion_ptr_in = NULL;
    struct QOP_D3_MDWF_Fermion *fermion_ptr_out = NULL;
    QOP_D3_MDWF_import_fermion(&fermion_ptr_in, mdwf_state, read_fermion_canonical, (void *)in);
    QOP_D3_MDWF_allocate_fermion(&fermion_ptr_out, mdwf_state);
    if(fermion_ptr_in == NULL || fermion_ptr_out == NULL){
      ERR.General(cname, fname, "mdwf library initialization failed.\n");
    }
    
    // apply D^\dag
    QOP_D3_MDWF_DDW_operator_conjugated(fermion_ptr_out, mdwf_param, mdwf_gauge_ptr.d, fermion_ptr_in);
    
    QOP_D3_MDWF_export_fermion(write_fermion_canonical, (void *)out, fermion_ptr_out);
    
    QOP_D3_MDWF_free_fermion(&fermion_ptr_in);
    QOP_D3_MDWF_free_fermion(&fermion_ptr_out);
  }
#else
  ERR.NotImplemented(cname, fname);
#endif
}

//! Multiplication by the square of the fermion matrix.
void DiracOpMdwf::MatDagMat(Vector *out, Vector *in)
{
  const char *fname = "MatDagMat(V*, V*)";
  VRB.Func(cname, fname);

#ifdef USE_MDWF
  if(use_single_precision){
    struct QOP_F3_MDWF_Fermion *fermion_ptr_in = NULL;
    struct QOP_F3_MDWF_Fermion *fermion_ptr_out = NULL;
    QOP_F3_MDWF_import_fermion(&fermion_ptr_in, mdwf_state, read_fermion_canonical, (void *)in);
    QOP_F3_MDWF_allocate_fermion(&fermion_ptr_out, mdwf_state);
    if(fermion_ptr_in == NULL || fermion_ptr_out == NULL){
      ERR.General(cname, fname, "mdwf library initialization failed.\n");
    }

    // apply D
    QOP_F3_MDWF_DDW_operator(fermion_ptr_out, mdwf_param, mdwf_gauge_ptr.f, fermion_ptr_in);
    // apply D^\dag
    QOP_F3_MDWF_DDW_operator_conjugated(fermion_ptr_in, mdwf_param, mdwf_gauge_ptr.f, fermion_ptr_out);
    
    QOP_F3_MDWF_export_fermion(write_fermion_canonical, (void *)out, fermion_ptr_in);
    
    QOP_F3_MDWF_free_fermion(&fermion_ptr_in);
    QOP_F3_MDWF_free_fermion(&fermion_ptr_out);
  }else{
    struct QOP_D3_MDWF_Fermion *fermion_ptr_in = NULL;
    struct QOP_D3_MDWF_Fermion *fermion_ptr_out = NULL;
    QOP_D3_MDWF_import_fermion(&fermion_ptr_in, mdwf_state, read_fermion_canonical, (void *)in);
    QOP_D3_MDWF_allocate_fermion(&fermion_ptr_out, mdwf_state);
    if(fermion_ptr_in == NULL || fermion_ptr_out == NULL){
      ERR.General(cname, fname, "mdwf library initialization failed.\n");
    }

    // apply D
    QOP_D3_MDWF_DDW_operator(fermion_ptr_out, mdwf_param, mdwf_gauge_ptr.d, fermion_ptr_in);
    // apply D^\dag
    QOP_D3_MDWF_DDW_operator_conjugated(fermion_ptr_in, mdwf_param, mdwf_gauge_ptr.d, fermion_ptr_out);
    
    QOP_D3_MDWF_export_fermion(write_fermion_canonical, (void *)out, fermion_ptr_in);
    
    QOP_D3_MDWF_free_fermion(&fermion_ptr_in);
    QOP_D3_MDWF_free_fermion(&fermion_ptr_out);
  }
#else
  ERR.NotImplemented(cname, fname);
#endif
}


//!< Computes vectors used in the HMD pseudofermionic force term.
// GRF
// chi is the solution to MatPcInv.  The user passes two full size
// CANONICAL fermion vectors with conversion enabled to the
// constructor.  Using chi, the function fills these vectors;
// the result may be used to compute the HMD fermion force.
void DiracOpMdwf::CalcHmdForceVecs(Vector *chi)
{
  const char *fname = "CalcHmdForceVecs()";
  ERR.NotImplemented(cname, fname);
}

//!< Not implemented
// Reflexion in s operator, needed for the hermitian version 
// of the dirac operator in the Ritz solver.
void DiracOpMdwf::Reflex(Vector *out, Vector *in)
{
  const char *fname = "Reflex()";
  ERR.NotImplemented(cname, fname);
}

//------------------------------------------------------------------
// DiracOpGlbSum(Float *): 
// The global sum used by InvCg. If s_nodes = 1
// it is the usual global sum. If s_nodes > 1 it
// is the 5-dimensional globals sum glb_sum_five.
//------------------------------------------------------------------
void DiracOpMdwf::DiracOpGlbSum(Float *float_p) {
  const char *fname = "DiracOpGlbSum()";
  ERR.NotImplemented(cname, fname);
//  if(GJP.Snodes() == 1) {
//    glb_sum(float_p);
//  }
//  else {
    glb_sum_five(float_p);
//  }
}

CPS_END_NAMESPACE
