#include <alg/alg_plaq.h>
#include <alg/common_arg.h>
#include <alg/no_arg.h>
#include <config.h>
#include <stdio.h>

CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOpMobiusEOFA class methods.

*/
//------------------------------------------------------------------
//
// d_op_mobius_eofa.C
//
// DiracOpMobiusEOFA is derived from the DiracOp base class. 
// DiracOpMobiusEOFA is the front end for a library that contains
// all Dirac operators associated with Dwf fermions.
//
//------------------------------------------------------------------
CPS_END_NAMESPACE

#include <comms/glb.h>
#include <util/dirac_op.h>
#include <util/dwf.h>
#include <util/error.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/mobius.h>
#include <util/time_cps.h>
#include <util/verbose.h>
#include <util/wilson.h>
  
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

DiracOpMobiusEOFA::DiracOpMobiusEOFA(Lattice& lat, Vector* f_out, Vector* f_in,
    Float _m1, Float _m2, Float _m3, Float _a, int _pm,
    CgArg* cg_arg, CnvFrmType cnv_frm_flg):
DiracOpWilsonTypes(lat, f_out, f_in, cg_arg, cnv_frm_flg),
  m1(_m1), m2(_m2), m3(_m3), a(_a), pm(_pm)
{
  cname = "DiracOpMobiusEOFA";
  const char *fname = "DiracOpMobiusEOFA(L&,V*,V*,F,F,F,F,i,CgArg*,CnvFrmType)";
  VRB.Func(cname, fname);

#ifndef USE_QUDA
  ERR.General(cname, fname, "DiracOpMobiusEOFA requires QUDA\n");
#endif

  // Check that GJP's EOFA flag is set
  if(!GJP.EOFA()){ 
    ERR.General(cname, fname, "Tried to construct DiracOpMobiusEOFA but GJP.EOFA() == false!?\n");
  }

  // Do the necessary conversions
  #undef PROFILE
  #ifdef PROFILE
  Float dtime = -dclock();
  #endif
  if(cnv_frm == CNV_FRM_YES){ lat.Convert(DWF_4D_EOPREC_EE, f_out, f_in); }
  else if(cnv_frm == CNV_FRM_NO){ lat.Convert(WILSON); }
  #ifdef PROFILE
  dtime += dclock();
  print_flops("lattice", "Convert()", 0, dtime);
  #endif

  // Initialize parameters
  DiracArg(dirac_arg);

  mobius_lib_arg = lat.FdiracOpInitPtr();

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // call to mobius_init but is done here again in case the user
  // has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the Fmobius object.
  //----------------------------------------------------------------
  Dwf* mobius_arg = static_cast<Dwf*>(mobius_lib_arg);

  const int& Ls              = GJP.SnodeSites();
  mobius_arg->ls             = Ls;
  mobius_arg->pc_type        = GJP.ZMobius_PC_Type();
  mobius_arg->mobius_kappa_b = 1.0 / ( 2.0 * (GJP.Mobius_b() * (4.0 - GJP.DwfHeight()) + GJP.DwfA5Inv()));
  mobius_arg->mobius_kappa_c = 1.0 / ( 2.0 * (GJP.Mobius_c() * (4.0 - GJP.DwfHeight()) - GJP.DwfA5Inv()));
  mobius_arg->zmobius_b.resize(Ls);
  mobius_arg->zmobius_c.resize(Ls);
  mobius_arg->zmobius_kappa_b.resize(Ls);
  mobius_arg->zmobius_kappa_c.resize(Ls);
  mobius_arg->zmobius_kappa_ratio.resize(Ls);
  for(int i=0; i<Ls; i++) 
  {
    mobius_arg->zmobius_b[i] = GJP.Mobius_b ();
    mobius_arg->zmobius_c[i] = GJP.Mobius_c ();
    mobius_arg->zmobius_kappa_b[i] = mobius_arg->mobius_kappa_b;
    mobius_arg->zmobius_kappa_c[i] = mobius_arg->mobius_kappa_c;
    mobius_arg->zmobius_kappa_ratio[i] =
    mobius_arg->zmobius_kappa_b[i] / mobius_arg->zmobius_kappa_c[i];
  }
}

void DiracOpMobiusEOFA::DiracArg(CgArg* arg){ 
 const char *fname= "DiracArg(CgArg*)";
 dirac_arg = arg; 
 VRB.Result(cname,fname,"dirac_arg=%p Inverter=%d\n",dirac_arg,dirac_arg->Inverter);
}

void DiracOpMobiusEOFA::CalcHmdForceVecs(Vector* v1, Vector* v2, Vector* phi1, Vector* phi2)
{
  const char* fname = "CalcHmdForceVecs(V*,V*,V*,V*)";
  const size_t f_size = GJP.VolNodeSites() * lat.FsiteSize() / ( lat.FchkbEvl() + 1 );

  Vector* tmp = static_cast<Vector*>( smalloc(cname, fname, "tmp", f_size*sizeof(Float)) );
  
  tmp -> CopyVec(v1, f_size);
  lat.g5R5(v1, tmp);

  sfree(tmp);
}

// Uses the chronological inversion method to forecast solutions across EOFA heatbath poles
void DiracOpMobiusEOFA::MinResExt(Vector* soln, Vector* src, Vector** soln_old, Vector** vm, int degree)
{
  const char* fname = "MinResExt(V*,V*,V**,V**,i)";

  size_t f_size = GJP.VolNodeSites() * lat.FsiteSize() / ( lat.FchkbEvl() + 1 );
  if(GJP.Gparity()){ f_size *= 2; }

  // Trivial cases
  if(degree == 0) {
    soln -> VecZero(f_size); 
    return;
  } else if(degree == 1) {
    soln -> CopyVec(soln_old[0], f_size);
    return;
  }

  Vector* r    = static_cast<Vector*>( smalloc(f_size*sizeof(Float), "r", fname, cname) );
  Vector* tmp1 = static_cast<Vector*>( smalloc(f_size*sizeof(Float), "tmp1", fname, cname) );
  Vector* tmp2 = static_cast<Vector*>( smalloc(f_size*sizeof(Float), "tmp2", fname, cname) );

  Float dot;
  Complex xp;

  // Array to hold the matrix elements
  Complex** G = static_cast<Complex**>( smalloc(degree*sizeof(Complex*), "G", fname, cname) );
  for(int i=0; i<degree; i++) {
    G[i] = static_cast<Complex*>( smalloc(degree*sizeof(Complex), "G[i]", fname, cname) );
  }

  // Solution and source vectors
  Complex* a = static_cast<Complex*>( smalloc(degree*sizeof(Complex), "a", fname, cname) );
  Complex* b = static_cast<Complex*>( smalloc(degree*sizeof(Complex), "b", fname, cname) );

  // Orthonormalize the vector basis
  for(int i=0; i<degree; i++) {
    soln_old[i] -> VecTimesEquFloat(1.0 / sqrt(soln_old[i]->NormSqGlbSum(f_size)), f_size);
    for(int j=i+1; j<degree; j++) {
      xp = soln_old[i] -> CompDotProductGlbSum(soln_old[j], f_size);
      soln_old[j] -> CTimesV1PlusV2(-xp, soln_old[i], soln_old[j], f_size);
    }
  }

  // Construct rhs
  for(int i=0; i<degree; i++)
  {
    b[i] = soln_old[i] -> CompDotProductGlbSum(src, f_size);
    tmp1 -> CopyVec(soln_old[i], f_size);
    lat.DtildeInv(tmp2, tmp1, this->m1);
    MatHerm(vm[i], tmp2);
    dot = lat.FhamiltonNode(vm[i], soln_old[i]);
    glb_sum(&dot);
    G[i][i] = Complex(dot, 0.0);
  }

  // Construct the matrix
  for(int j=0; j<degree; j++){
  for(int k=j+1; k<degree; k++){
    G[j][k] = soln_old[j] -> CompDotProductGlbSum(vm[k], f_size);
    G[k][j] = conj(G[j][k]);
  }}

  // Gauss-Jordan elimination with partial pivoting
  for(int i=0; i<degree; i++)
  {
    // Perform partial pivoting
    int k = i;
    for(int j=i+1; j<degree; j++){ if(abs(G[j][j]) > abs(G[k][k])){ k = j; } }
    if(k != i) {
      xp = b[k];
      b[k] = b[i];
      b[i] = xp;
      for(int j=0; j<degree; j++) {
        xp = G[k][j];
        G[k][j] = G[i][j];
        G[i][j] = xp;
      }
    }

    // Convert matrix to upper triangular form
    for(int j=i+1; j<degree; j++) {
      xp = G[j][i] / G[i][i];
      b[j] -= xp * b[i];
      for(int k=0; k<degree; k++){ G[j][k] -= xp * G[i][k]; }
    }
  }

  // Use Gaussian elimination to solve equations and calculate initial guesses
  tmp1 -> VecZero(f_size);
  r -> CopyVec(src, f_size);
  for(int i=degree-1; i>=0; i--) {
    a[i] = 0.0;
    for(int j=i+1; j<degree; j++){ a[i] += G[i][j] * a[j]; }
    a[i] = ( b[i] - a[i] ) / G[i][i];
    tmp1 -> CTimesV1PlusV2(a[i], soln_old[i], tmp1, f_size);
    r -> CTimesV1PlusV2(-a[i], vm[i], r, f_size);
  }
  Float error = sqrt( r->NormSqGlbSum(f_size) / src->NormSqGlbSum(f_size) );

  VRB.Result(cname, fname, "pole %d: |res|/|src| = %e\n", degree, error);

  // Now we have a forecasted guess for the unpreconditioned system.
  // Transform to the preconditioned system.
  lat.DtildeInv(soln, tmp1, this->m1);

  for(int i=degree-1; i>=0; i--){ sfree(G[i], "G[i]", fname, cname); }
  sfree(a, "a", fname, cname);
  sfree(b, "b", fname, cname);
  sfree(G, "G", fname, cname);
  sfree(r, "r", fname, cname);
  sfree(tmp1, "tmp1", fname, cname);
  sfree(tmp2, "tmp2", fname, cname);
}

DiracOpMobiusEOFA::~DiracOpMobiusEOFA()
{
  const char* fname = "~DiracOpMobiusEOFA()";
  VRB.Func(cname, fname);

  // Do the necessary conversions
  #undef PROFILE
  #ifdef PROFILE
  Float dtime = -dclock();
  #endif  
  if(cnv_frm == CNV_FRM_YES){ lat.Convert(CANONICAL, f_out, f_in); }
  else if(cnv_frm == CNV_FRM_NO){ lat.Convert(CANONICAL); }
  #ifdef PROFILE
  dtime += dclock();
  print_flops("lattice", "Convert()", 0, dtime);
  #endif
}

CPS_END_NAMESPACE
