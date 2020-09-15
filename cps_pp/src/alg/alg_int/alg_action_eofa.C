#include <config.h>
#include <math.h>
#include <string.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------------
//
// alg_action_eofa.C
//
//-------------------------------------------------------------------------
CPS_END_NAMESPACE
#include <alg/alg_hmd.h>
#include <alg/alg_int.h>
#include <alg/alg_remez.h>
#include <util/error.h>
#include <util/gjp.h>
#include <util/lattice.h>
#ifdef USE_BFM
#include <util/lattice/bfm_evo.h>
#include <util/lattice/fbfm.h>
#endif
#include <util/smalloc.h>
#include <util/time_cps.h>
#include <util/timer.h>
#include <util/vector.h>
#include <util/verbose.h>

CPS_START_NAMESPACE

// Dummy constructor - does nothing
AlgActionEOFA::AlgActionEOFA() : AlgActionBilinear(){}

AlgActionEOFA::AlgActionEOFA(AlgMomentum& mom, ActionEOFAArg& arg, bool _heatbath_forecast, bool _heatbath_test, int traj_num) : AlgActionBilinear(mom, arg.bi_arg), heatbath_forecast(_heatbath_forecast), heatbath_test(_heatbath_test)
{
  cname = "AlgActionEOFA";
  const char* fname = "AlgActionEOFA";
  VRB.Func(cname, fname);

  GJP.EnableEOFA();

  if(!UniqueID()){ printf("AlgActionEOFA constructor started\n"); fflush(stdout); }

  VRB.Result(cname, fname, "Recomputing fermion field sizes for EOFA\n");
  {
    // Recompute sizes
    Lattice& lat = LatticeFactory::Create(arg.bi_arg.fermion, G_CLASS_NONE);
    this->f_size      = GJP.VolNodeSites() * lat.FsiteSize() / ( lat.FchkbEvl() + 1 );
    this->f_vec_count = this->f_size / ( 2 * lat.Colors() );
    this->f_sites     = this->f_size / ( 2 * lat.Colors() * lat.SpinComponents() );
    this->Ncb         = (lat.FchkbEvl() == 0) ? 2 : 1;
    VRB.Result(cname, fname, "f_sites = %d, f_vec_count = %d, f_size = %d (lat.FchkbEvl() = %d)\n", 
        f_sites, f_vec_count, f_size, lat.FchkbEvl());
    LatticeFactory::Destroy();

    // Delete and reallocate phi fields
    for(int i=0; i<this->n_masses; i++) {
      sfree(phi[i], cname, fname, "phi[i]");
      phi[i] = static_cast<Vector*>( smalloc(f_size*sizeof(Float), "phi[i]", fname, cname) );
    }
  }

  int_type = INT_EOFA;
  eofa_arg = &arg;

  // First check that number of numerator masses matches number of denominator masses
  if(eofa_arg->num_mass.num_mass_len != eofa_arg->den_mass.den_mass_len){
    ERR.General(cname, fname, "Inconsistency between number of numerator and denominator masses\n");
  }

  // Also check that number of bilinear masses matches number of numerator masses
  if(eofa_arg->bi_arg.bilinears.bilinears_len != eofa_arg->num_mass.num_mass_len){
    ERR.General(cname, fname, "Inconsistency between number of bilinears and numerator masses\n");
  }

  // Also check that rational approximation parameters were supplied for each mass
  if(eofa_arg->LH_rat_approx.LH_rat_approx_len != eofa_arg->num_mass.num_mass_len){
    ERR.General(cname, fname, "Inconsistency between number of LH rational approximation parameters and numerator masses\n");
  }
  if(eofa_arg->RH_rat_approx.RH_rat_approx_len != eofa_arg->num_mass.num_mass_len){
    ERR.General(cname, fname, "Inconsistency between number of RH rational approximation parameters and numerator masses\n");
  }

  // Also check all EOFARationalDescr are specified as fermions
  for(int i=0; i<n_masses; i++){
    if(eofa_arg->LH_rat_approx.LH_rat_approx_val[i].field_type != FERMION){
      ERR.General(cname, fname, "LH_rat_approx %d not set as a fermion\n", i);
    } else if(eofa_arg->RH_rat_approx.RH_rat_approx_val[i].field_type != FERMION){
      ERR.General(cname, fname, "RH_rat_approx %d not set as a fermion\n", i);
    }
  }

  // For now only QUDA back-end has been tested with this branch
  // FIXME: BFM is probably deprecated, but should get this to work with Grid, too.
  if((fermion != F_CLASS_BFM) && (fermion != F_CLASS_MOBIUS)){
    ERR.General(cname, fname, "Fermion type not implemented for EOFA");
  }

  #ifndef USE_QUDA
  if(fermion == F_CLASS_MOBIUS){
    ERR.General(cname, fname, "EOFA with F_CLASS_MOBIUS requires QUDA\n");
  }
  #endif

  // Allocate memory for the fermion CG arguments
  if(n_masses > 0)
  {
    #ifdef USE_BFM
    // AlgActionBilinear does not set fermion field size correctly for Fbfm
    if(eofa_arg->bi_arg.fermion == F_CLASS_BFM)
    {
      int Ls = Fbfm::arg_map.at(eofa_arg->num_mass.num_mass_val[0]).Ls;
      VRB.Result(cname, fname, "Recalculating fermion field size for Fbfm based on Ls = %d\n", Ls);
      
      Lattice& lat = LatticeFactory::Create(fermion, G_CLASS_NONE);
      f_size = GJP.VolNodeSites() * Ls * (2*3*4) / ( lat.FchkbEvl() + 1 ); // (reim * color * spin) 
      f_vec_count = f_size / (2*3);
      f_sites = f_size / (2*3*4);

      VRB.Result(cname, fname, "Allocating phi fields\n");
      for(int i=0; i<n_masses; i++){
        phi[i] = static_cast<Vector*>( smalloc(this->f_size*sizeof(Float), "phi[i]", fname, cname) );
      }

      LatticeFactory::Destroy();
    }
    #endif

    num_mass = static_cast<Float*>( smalloc(n_masses * sizeof(Float), "num_mass", fname, cname) );
    den_mass = static_cast<Float*>( smalloc(n_masses * sizeof(Float), "den_mass", fname, cname) );

    for(int i=0; i<n_masses; i++)
    {
      num_mass[i] = eofa_arg -> num_mass.num_mass_val[i];
      den_mass[i] = eofa_arg -> den_mass.den_mass_val[i];

      #ifdef USE_BFM
      if(eofa_arg->bi_arg.fermion == F_CLASS_BFM)
      {
        // Make sure all quotients have the same Ls
        int Ls = Fbfm::arg_map.at(num_mass[0]).Ls;
        if(Fbfm::arg_map.at(num_mass[i]).Ls != Ls){
          ERR.General(cname, fname, "Numerator mass %d doesn't have the same Ls as numerator mass 0!\n", i);
        }
        if(Fbfm::arg_map.at(den_mass[i]).Ls != Ls){
          ERR.General(cname, fname, "Denominator mass %d doesn't have the same Ls as numerator mass 0!\n", i);
        }
      }
      #endif
    }

    // Construct approximation if necessary
    if(!loadPoles()){
      if(!UniqueID()){ printf("Generating rational approximation\n"); fflush(stdout); }
      generateApprox(num_mass, &LH_remez_arg, eofa_arg->LH_rat_approx.LH_rat_approx_val);
      generateApprox(num_mass, &RH_remez_arg, eofa_arg->RH_rat_approx.RH_rat_approx_val);
      savePoles();
      if(!UniqueID()){ printf("Finished generating rational approximation\n"); fflush(stdout); }
    }

    // Generate CG args for fg, mc, and md solves
    generateCgArg(num_mass, &LH_cg_arg_fg, &LH_cg_arg_mc, &LH_cg_arg_md, "LH_cg_arg", 
        eofa_arg->LH_stop_rsd_fg.LH_stop_rsd_fg_val,
        eofa_arg->LH_stop_rsd_mc.LH_stop_rsd_mc_val,
        eofa_arg->LH_stop_rsd_md.LH_stop_rsd_md_val);
    generateCgArg(num_mass, &RH_cg_arg_fg, &RH_cg_arg_mc, &RH_cg_arg_md, "RH_cg_arg", 
        eofa_arg->RH_stop_rsd_fg.RH_stop_rsd_fg_val,
        eofa_arg->RH_stop_rsd_mc.RH_stop_rsd_mc_val,
        eofa_arg->RH_stop_rsd_md.RH_stop_rsd_md_val);

    // Generate CG args for heatbath solves
    generateCgArg(num_mass, &LH_cg_arg_heatbath, "LH_cg_arg_heatbath", eofa_arg->LH_rat_approx.LH_rat_approx_val);
    generateCgArg(num_mass, &RH_cg_arg_heatbath, "RH_cg_arg_heatbath", eofa_arg->RH_rat_approx.RH_rat_approx_val);

    // Temporary fermion fields used in intermediate calculations for energy, heatbath, and evolution
    frmn_tmp = static_cast<Vector**>( smalloc((2*n_masses+3)*sizeof(Vector*), "frmn_tmp", fname, cname) );
    for(int i=0; i<2*n_masses+3; i++){
      frmn_tmp[i] = static_cast<Vector*>( smalloc(f_size*sizeof(Float), "frmn_tmp[i]", fname, cname) );
    }
  }

  init(traj_num);
  fg_forecast = false;

  if(eofa_arg->eigen.eigen_measure == EIGEN_MEASURE_YES){ generateEigArg(eofa_arg->eigen); }

  GJP.DisableEOFA();

  if(!UniqueID()){ printf("AlgActionEOFA constructor finished\n"); fflush(stdout); }
}

void AlgActionEOFA::init(int traj_num)
{
  const char* fname = "init(i)";
  VRB.Func(cname, fname);
  GJP.EnableEOFA();
  AlgActionBilinear::init();
  evolved = 1;
  heatbathEval = 0;
  energyEval = 0;
//  traj = traj_num - 1;
  GJP.DisableEOFA();
  VRB.FuncEnd(cname, fname);
}

// Computes largest and smallest eigenvalues of fermion matrix
std::vector<double> AlgActionEOFA::eig_range(Float m_num, Float m_den)
{
  const char* fname = "eig_range(F,F)";
  GJP.EnableEOFA();

  std::vector<double> eigs(2);

  VRB.Result(cname, fname, "Computing largest and smallest eigenves of Meofa for mass ratio = %f/%f\n", m_num, m_den);
  eigs[0] = 1.0; // analytic result
//  eigs[0] = ritz(m_num, m_den, true);//tested. measured 1!
  eigs[1] = ritz(m_num, m_den, false);

  GJP.DisableEOFA();

  return eigs;
}

#ifndef USE_BFM
static inline double quad_solve(double* ct, double* st, double a, 
    double b, double c, double d, double e, double f)
{
  double p = b * (d-f) + e * (c-a);
  double q = b * (d+f) - e * (c+a);
  double r = 2 * ( c*d - a*f );

  // solve: p + q*cos(2t) + r*sin(2t) = 0
  double den = sqrt( q*q + r*r );
  double ca = q / den;

  double ci = sqrt( 0.5 * ( 1.0 + ca) );
  double si = sqrt( 0.5 * ( 1.0 - ca) );
  if(r < 0){ si = -si; }

  double cb = -p / den;
  if(fabs(cb) > 1.0) {
    printf("Panic: cos(psi) > 1\n");
    exit(-1);
  }
  double cj = sqrt( 0.5 * ( 1.0 + cb ) );
  double sj = sqrt( 0.5 * ( 1.0 - cb ) );
  
  double ct1 = ci*cj + si*sj;
  double st1 = si*cj - ci*sj;
  double v1  = ( a*ct1*ct1 + b*st1*ct1 + c*st1*st1 ) / ( d*ct1*ct1 + e*st1*ct1 + f*st1*st1 );

  double ct2 = ci*cj - si*sj;
  double st2 = si*cj + ci*sj;
  double v2  = ( a*ct2*ct2 + b*st2*ct2 + c*st2*st2 ) / ( d*ct2*ct2 + e*st2*ct2 + f*st2*st2 );

  if(v1 < v2) {
    *ct = ct1;
    *st = st1;
    return v1;
  } else {
    *ct = ct2;
    *st = st2;
    return v2;
  }
}
#endif

// Simple ritz implementation for calculating eigenvalues of EOFA action
double AlgActionEOFA::ritz(Float m_num, Float m_den, bool compute_min)
{
  const char* fname = "ritz(F,F)";
  double eval;

  // Create an appropriate lattice
  Lattice& lat = LatticeFactory::Create(fermion, G_CLASS_NONE);

  Vector* x = static_cast<Vector*>( smalloc(f_size*sizeof(Float), "x", fname, cname) );
  Vector* y = static_cast<Vector*>( smalloc(f_size*sizeof(Float), "y", fname, cname) );
  Vector* p = static_cast<Vector*>( smalloc(f_size*sizeof(Float), "p", fname, cname) );
  Vector* z = static_cast<Vector*>( smalloc(f_size*sizeof(Float), "z", fname, cname) );
  Vector* t = static_cast<Vector*>( smalloc(f_size*sizeof(Float), "t", fname, cname) );
  Vector* u = static_cast<Vector*>( smalloc(f_size*sizeof(Float), "u", fname, cname) );

  double mu, pnorm, gnorm2, xnorm;

  // Draw normalized x
  #ifdef USE_BFM
  if(eofa_arg->bi_arg.fermion == F_CLASS_BFM){
    Fbfm::current_key_mass = num_mass;
    VRB.Result(cname, fname, "Setting Fbfm::current_key_mass = %e before calling RandGaussVector\n", Fbfm::current_key_mass);
  }
  #endif
  lat.RandGaussVector(x, 0.5, Ncb);
  double fact = x -> NormSqGlbSum(f_size);
  fact = 1.0 / sqrt(fact);
  x -> VecTimesEquFloat(fact, f_size);

  VRB.Result(cname, fname, "<x,x> = %1.10en", 1.0/(fact*fact));

  Meofa(y, x, m_num, m_den, lat); // y = Meofa*x
  mu = x -> ReDotProductGlbSum(y, f_size);
  if(!compute_min){
    mu = -mu;
    y -> VecTimesEquFloat(-1.0, f_size);
  }

  p -> FTimesV1PlusV2(-mu, x, y, f_size); // p = y - mu*x
  gnorm2 = p -> NormSqGlbSum(f_size);
  pnorm = sqrt(gnorm2);

  int i;
  double stop_tol = eofa_arg->eigen.stop_rsd * eofa_arg->eigen.stop_rsd;
  VRB.Result(cname, fname, "stopping tolerance: gnorm2 = %1.8e\n", stop_tol);
  for(i=0; i<eofa_arg->eigen.max_num_iter; i++)
  {
    // current state
    xnorm = x -> NormSqGlbSum(f_size);
    eval = (compute_min) ? (mu/xnorm) : (-mu/xnorm);
    VRB.Result(cname, fname, "iter %d: gnorm2 = %1.8e, mu = %1.8e, eval = %1.8e\n", i, gnorm2, mu, eval);

    // check for convergence
    if(gnorm2 < stop_tol){ break; }

    Meofa(z, p, m_num, m_den, lat);
    double pap = p -> ReDotProductGlbSum(z, f_size);
    if(!compute_min){
      pap = -pap;
      z -> VecTimesEquFloat(-1.0, f_size);
    }

    // minimize x cos(theta) + p / pnorm * sin(theta) in theta
    double d = x -> NormSqGlbSum(f_size);
    double e = (2.0/pnorm) * ( x -> ReDotProductGlbSum(p, f_size) );
    double f = 1.0;
    double a = mu * d;
    double b = (2.0/pnorm) * ( x -> ReDotProductGlbSum(z, f_size) );
    double c = pap / pnorm / pnorm;

    double ct, st;
    mu = quad_solve(&ct, &st, a, b, c, d, e, f);

    x -> VecTimesEquFloat(ct, f_size);
    x -> FTimesV1PlusV2(st/pnorm, p, x, f_size); // x = ct*x + st/pnorm*p
    y -> VecTimesEquFloat(ct, f_size);
    y -> FTimesV1PlusV2(st/pnorm, z, y, f_size); // y = ct*y + st/pnorm*z

    t -> FTimesV1PlusV2(-mu, x, y, f_size);
    double gnew = t -> NormSqGlbSum(f_size);
    double beta = ct * gnew / gnorm2;
    gnorm2 = gnew;

    double xpp = x -> ReDotProductGlbSum(p, f_size);
    u -> FTimesV1PlusV2(-xpp, x, p, f_size);

    p -> FTimesV1PlusV2(beta, u, t, f_size);
    pnorm = sqrt( p -> NormSqGlbSum(f_size) );
  }

  if(i < eofa_arg->eigen.max_num_iter) {
    VRB.Result(cname, fname, "converged at iteration %d.\n", i);
  } else {
    VRB.Result(cname, fname, "maximum iteration number reached!\n");
  }

  if(!compute_min){ mu = -mu; }
  xnorm = x -> NormSqGlbSum(f_size);
  eval = mu / xnorm;
  VRB.Result(cname, fname, "result: eval = %1.8e\n", eval);

  sfree(u, "u", fname, cname);
  sfree(t, "t", fname, cname);
  sfree(z, "z", fname, cname);
  sfree(p, "p", fname, cname);
  sfree(y, "y", fname, cname);
  sfree(x, "x", fname, cname);
  LatticeFactory::Destroy();

  return eval;
}

// Computes out = Meofa * in
void AlgActionEOFA::Meofa(Vector* out, const Vector* in, Float m_num, Float m_den, Lattice& lat)
{
  const char* fname = "Meofa(V*,V*,F,F,Lattice*)";
  
  GJP.EnableEOFA();

  // Zero tmp vectors
  for(int j=0; j<3; j++){ frmn_tmp[j] -> VecZero(f_size); }

  // out = in
  out -> CopyVec(in, f_size);

  lat.SetEOFAParams(m_num, 0.0, 0.0, 0.0, 1);
  Float kt = lat.kt(m_num, m_den);

  // Contribution from LH term
  lat.ChiralProj(frmn_tmp[0], in, -1);
  lat.Omega(frmn_tmp[1], frmn_tmp[0], -1, 0);
  lat.FmatEvlMeofa(frmn_tmp[2], frmn_tmp[1], LH_cg_arg_mc[0], m_num, 0.0, 0.0, 0.0, -1);
  lat.Omega(frmn_tmp[1], frmn_tmp[2], -1, 1);
  lat.ChiralProj(frmn_tmp[0], frmn_tmp[1], -1);
  frmn_tmp[1] -> VecEqualsVecTimesEquFloat(frmn_tmp[0], -kt, f_size);
  out -> VecAddEquVec(frmn_tmp[1], f_size); // out = out - k*Pm*Om^\dag*H(mf)^{-1}*Om*Pm*in

  // Contribution from RH term
  lat.ChiralProj(frmn_tmp[0], in, 1);
  lat.Omega(frmn_tmp[1], frmn_tmp[0], 1, 0);
  lat.FmatEvlMeofa(frmn_tmp[2], frmn_tmp[1], RH_cg_arg_mc[0], m_den, m_num, m_den, -1.0, 1);
  lat.Omega(frmn_tmp[1], frmn_tmp[2], 1, 1);
  lat.ChiralProj(frmn_tmp[0], frmn_tmp[1], 1);
  frmn_tmp[1] -> VecEqualsVecTimesEquFloat(frmn_tmp[0], kt, f_size);
  out -> VecAddEquVec(frmn_tmp[1], f_size); // out = out + k*Pp*Op^\dag*(H(mb)-Dp(mf,mb)*Pp)^{-1}*Op*Pp*in

  GJP.DisableEOFA();
}

// Calculate EOFA pseudofermion contribution to the Hamiltonian
Float AlgActionEOFA::energy()
{
  char fname[20 + strlen(force_label)];
  sprintf(fname, "energy()[%s]", force_label);

  if(energyEval) {
    return 0.0;
  } 

  else if(!evolved) 
  {
    energyEval = 1;
    
    {
      Float gsum_h(h_init);
      glb_sum(&gsum_h);
      if(UniqueID() == 0){ printf("AlgActionEOFA::energy() [%s] %1.8e\n", force_label, gsum_h); }
    }
    
    return h_init;
  } 

  else 
  {
    static Timer timer(cname, fname);
    timer.start(true);
    Float dtime = -dclock();
    int shift(0);
    Float h(0.0);

    GJP.EnableEOFA();

    // Create an appropriate lattice
    Lattice& lat = LatticeFactory::Create(fermion, G_CLASS_NONE);

    Float total_h_i;
    for(int i=0; i<n_masses; i++)
    {
      // Zero tmp vectors
      for(int j=0; j<3; j++){ frmn_tmp[j] -> VecZero(f_size); }

      lat.SetEOFAParams(num_mass[i], 0.0, 0.0, 0.0, 1);
      Float kt = lat.kt(num_mass[i], den_mass[i]);

      Float h_i = lat.FhamiltonNode(phi[i], phi[i]);

      // LH term
      lat.ChiralProj(frmn_tmp[0], phi[i], -1);
      lat.Omega(frmn_tmp[1], frmn_tmp[0], -1, 0);
      frmn_tmp[2] -> VecZero(f_size);
      lat.FmatEvlMeofa(frmn_tmp[2], frmn_tmp[1], LH_cg_arg_mc[i], num_mass[i], 0.0, 0.0, 0.0, -1);
      lat.Omega(frmn_tmp[1], frmn_tmp[2], -1, 1);
      frmn_tmp[2] -> VecEqualsVecTimesEquFloat(frmn_tmp[1], -kt, f_size);
      h_i += lat.FhamiltonNode(frmn_tmp[0], frmn_tmp[2]);

      // RH term
      lat.ChiralProj(frmn_tmp[0], phi[i], 1);
      lat.Omega(frmn_tmp[1], frmn_tmp[0], 1, 0);
      frmn_tmp[2] -> VecZero(f_size);
      lat.FmatEvlMeofa(frmn_tmp[2], frmn_tmp[1], RH_cg_arg_mc[i], den_mass[i], num_mass[i], den_mass[i], -1.0, 1);
      lat.Omega(frmn_tmp[1], frmn_tmp[2], 1, 1);
      frmn_tmp[2] -> VecEqualsVecTimesEquFloat(frmn_tmp[1], kt, f_size);
      h_i += lat.FhamiltonNode(frmn_tmp[0], frmn_tmp[2]);

      updateCgStats(LH_cg_arg_mc[i]);
      updateCgStats(RH_cg_arg_mc[i]);

      h += h_i;
      total_h_i = h_i;
      glb_sum(&total_h_i);
      VRB.Result(cname, fname, "mass ratio %f/%f: final H = %1.15e\n", num_mass[i], den_mass[i], total_h_i);
    }

    LatticeFactory::Destroy();

    energyEval = 1;

    dtime += dclock();
    print_flops(cname, fname, 0, dtime);
    timer.stop(true);

    total_h_i = h - h_init;
    glb_sum(&total_h_i);
    VRB.Result(cname, fname, "dH = %1.15e\n", total_h_i);

    GJP.DisableEOFA();

    return h;
  }
}

void AlgActionEOFA::reweight(Float* rw_fac, Float* norm)
{
  const char* fname = "reweight(F*,F*)";

  GJP.EnableEOFA();

  Lattice& lat = LatticeFactory::Create(fermion, G_CLASS_NONE);

  for(int i=0; i<n_masses; i++)
  {
    #ifdef USE_BFM
    if(eofa_arg->bi_arg.fermion == F_CLASS_BFM){
      Fbfm::current_key_mass = num_mass[i];
      VRB.Result(cname, fname, "Setting Fbfm::current_key_mass = %e before calling RandGaussVector\n", Fbfm::current_key_mass);
    }
    #endif   
    lat.RandGaussVector(phi[i], 0.5, Ncb);
    // IFloat* pi = reinterpret_cast<IFloat*>(phi[i]);
    // for(int idx=0; idx<f_size; idx++){ pi[idx] = 1.0; }

    // Zero tmp vectors
    for(int j=0; j<3; j++){ frmn_tmp[j] -> VecZero(f_size); }

    if(GJP.Gparity1fX() && GJP.Gparity1fY())
    {
      if(!UniqueID()){ printf("Putting minus sign on fermion source in UR quadrant\n"); fflush(stdout); }

      // Make source on upper-right quadrant negative (RNGs should be correct)
      for(int s=0; s<GJP.SnodeSites(); s++){
      for(int t=0; t<GJP.TnodeSites(); t++){
      for(int z=0; z<GJP.ZnodeSites(); z++){
      for(int y=0; y<GJP.YnodeSites(); y++){
      for(int x=0; x<GJP.XnodeSites(); x++){

        int gx = x + GJP.XnodeCoor() * GJP.XnodeSites();
        int gy = y + GJP.YnodeCoor() * GJP.YnodeSites();

        if((gx >= GJP.Xnodes()*GJP.XnodeSites()/2) && (gy >= GJP.Ynodes()*GJP.YnodeSites()/2))
        {
          int pos[5] = {x, y, z, t, s};
          int f_off  = lat.FsiteOffsetChkb(pos) * lat.SpinComponents();
          for(int spn=0; spn<lat.SpinComponents(); spn++){ *(frmn_tmp[0]+f_off+spn) *= -1.0; }
        }
      }}}}}
    }

    norm[i] = lat.FhamiltonNode(phi[i], phi[i]);

    lat.SetEOFAParams(num_mass[i], 0.0, 0.0, 0.0, 1);
    Float kt = lat.kt(num_mass[i], den_mass[i]);

    // Contribution from LH term
    lat.ChiralProj(frmn_tmp[0], phi[i], -1);
    lat.Omega(frmn_tmp[1], frmn_tmp[0], -1, 0);
    frmn_tmp[2] -> VecZero(f_size);
    lat.FmatEvlMeofa(frmn_tmp[2], frmn_tmp[1], LH_cg_arg_mc[i], num_mass[i], 0.0, 0.0, 0.0, -1);
    lat.Omega(frmn_tmp[1], frmn_tmp[2], -1, 1);
    frmn_tmp[2] -> VecEqualsVecTimesEquFloat(frmn_tmp[1], -kt, f_size);
    rw_fac[i] = lat.FhamiltonNode(frmn_tmp[0], frmn_tmp[2]);

    // Contribution from RH term
    lat.ChiralProj(frmn_tmp[0], phi[i], 1);
    lat.Omega(frmn_tmp[1], frmn_tmp[0], 1, 0);
    frmn_tmp[2] -> VecZero(f_size);
    lat.FmatEvlMeofa(frmn_tmp[2], frmn_tmp[1], RH_cg_arg_mc[i], den_mass[i], num_mass[i], den_mass[i], -1.0, 1);
    lat.Omega(frmn_tmp[1], frmn_tmp[2], 1, 1);
    frmn_tmp[2] -> VecEqualsVecTimesEquFloat(frmn_tmp[1], kt, f_size);
    rw_fac[i] += lat.FhamiltonNode(frmn_tmp[0], frmn_tmp[2]);

    updateCgStats(LH_cg_arg_mc[i]);
    updateCgStats(RH_cg_arg_mc[i]);

    glb_sum(rw_fac+i);
    glb_sum(norm+i);
  }

  LatticeFactory::Destroy();

  GJP.DisableEOFA();
}

void AlgActionEOFA::heatbath()
{
  char fname[20 + strlen(force_label)];
  sprintf(fname, "heatbath()[%s]", force_label);

  static Timer timer(cname, fname);
  timer.start(true);
  Float dtime2 = -dclock(true);

  std::string veloc_label;
  veloc_label ="Phi_traj"+std::to_string(traj_num);
  int veloc_v = getVer(veloc_label.c_str());
  VRB.Result(cname,fname,"label version=%s %d\n",veloc_label.c_str(),veloc_v);

  Float dtime = -dclock(true);
  GJP.EnableEOFA();
    dtime +=dclock(true);
    print_flops(fname,"EnableEOFA()",0,dtime);
    dtime =-dclock(true);

  if(!heatbathEval)
  {
    // If enabled, check that the spectral range of Meofa lies inside the bounds of our rational approximation
    if(eofa_arg->eigen.eigen_measure == EIGEN_MEASURE_YES)
    {
      std::vector<double> eigs(2);
      for(int i=0; i<n_masses; i++)
      {
        eigs = eig_range(num_mass[i], den_mass[i]);

        if((eigs[0] < LH_remez_arg[i].lambda_low) || (eigs[0] < RH_remez_arg[i].lambda_low)){
          ERR.General(cname, fname, "Lower bound exceeded: mass ratio = %f/%f, %e < %e\n",
              num_mass[i], den_mass[i], eigs[0], LH_remez_arg[i].lambda_low);
        } else {
          VRB.Result(cname, fname, "Lower bound valid: mass ratio = %f/%f, %e > %e\n",
              num_mass[i], den_mass[i], eigs[0], LH_remez_arg[i].lambda_low);
        }

        if((eigs[1] > LH_remez_arg[i].lambda_high) || (eigs[1] > RH_remez_arg[i].lambda_high)){
          ERR.General(cname, fname, "Upper bound exceeded: mass ratio = %f/%f, %e > %e\n",
              num_mass[i], den_mass[i], eigs[0], LH_remez_arg[i].lambda_high);
        } else {
          VRB.Result(cname, fname, "Upper bound valid: mass ratio = %f/%f, %e < %e\n",
              num_mass[i], den_mass[i], eigs[0], LH_remez_arg[i].lambda_high);
        }
      }
    }
    dtime +=dclock(true);
    print_flops(fname,"eigen_measure()",0,dtime);
    dtime =-dclock(true);


    Lattice& lat = LatticeFactory::Create(fermion, G_CLASS_NONE);
    h_init = 0.0;

    for(int i=0; i<n_masses; i++)
    {

      lat.RandGaussVector(phi[i], 0.5, Ncb);

      // Zero tmp vectors
      for(int j=0; j<3; j++){ frmn_tmp[j] -> VecZero(f_size); }

      if(GJP.Gparity1fX() && GJP.Gparity1fY())
      {
        if(!UniqueID()){ printf("Putting minus sign on fermion source in UR quadrant\n"); fflush(stdout); }

        // Make source on upper-right quadrant negative (RNGs should be correct)
        for(int s=0; s<GJP.SnodeSites(); s++){
        for(int t=0; t<GJP.TnodeSites(); t++){
        for(int z=0; z<GJP.ZnodeSites(); z++){
        for(int y=0; y<GJP.YnodeSites(); y++){
        for(int x=0; x<GJP.XnodeSites(); x++){

          int gx = x + GJP.XnodeCoor() * GJP.XnodeSites();
          int gy = y + GJP.YnodeCoor() * GJP.YnodeSites();

          if((gx >= GJP.Xnodes()*GJP.XnodeSites()/2) && (gy >= GJP.Ynodes()*GJP.YnodeSites()/2))
          {
            int pos[5] = {x, y, z, t, s};
            int f_off  = lat.FsiteOffsetChkb(pos) * lat.SpinComponents();
            for(int spn=0; spn<lat.SpinComponents(); spn++){ *(frmn_tmp[0]+f_off+spn) *= -1.0; }
          }
        }}}}}
      }
    dtime +=dclock(true);
    print_flops(fname,"VecZero()",0,dtime);
    dtime =-dclock(true);
      
      Float delta_h = lat.FhamiltonNode(phi[i], phi[i]);
      Float gsum_h(delta_h);
      {
        glb_sum(&gsum_h);
        VRB.Result(cname, fname, "mass ratio %f/%f: initial H = %1.15e\n", num_mass[i], den_mass[i], gsum_h);
      }
    dtime +=dclock(true);
    print_flops(fname,"FhamiltonNode()",0,dtime);
    dtime =-dclock(true);

      h_init += delta_h;

      // Stored solutions (CG_solns) and M*solutions (vm) for heatbath forecasting
      Vector** CG_solns_LH = static_cast<Vector**>( smalloc(LH_remez_arg[i].degree*sizeof(Vector*), "CG_solns_LH", fname, cname) );
      Vector** CG_solns_RH = static_cast<Vector**>( smalloc(RH_remez_arg[i].degree*sizeof(Vector*), "CG_solns_RH", fname, cname) );
      Vector** vm_LH       = static_cast<Vector**>( smalloc(LH_remez_arg[i].degree*sizeof(Vector*), "vm_LH", fname, cname) );
      Vector** vm_RH       = static_cast<Vector**>( smalloc(RH_remez_arg[i].degree*sizeof(Vector*), "vm_RH", fname, cname) );

      for(int j=0; j<LH_remez_arg[i].degree; j++){
        CG_solns_LH[j] = static_cast<Vector*>( smalloc(f_size*sizeof(Float), "CG_solns_LH[j]", fname, cname) );
        vm_LH[j]       = static_cast<Vector*>( smalloc(f_size*sizeof(Float), "vm_LH[j]", fname, cname) );
      }
      for(int j=0; j<RH_remez_arg[i].degree; j++){
        CG_solns_RH[j] = static_cast<Vector*>( smalloc(f_size*sizeof(Float), "CG_solns_RH[j]", fname, cname) );
        vm_RH[j]       = static_cast<Vector*>( smalloc(f_size*sizeof(Float), "vm_RH[j]", fname, cname) );
      }

      lat.SetEOFAParams(num_mass[i], 0.0, 0.0, 0.0, -1);
      Float kt = lat.kt(num_mass[i], den_mass[i]);

      // LH term
      lat.ChiralProj(frmn_tmp[0], phi[i], -1);
      lat.Omega(frmn_tmp[1], frmn_tmp[0], -1, 0);
      Float norm = LH_remez_arg[i].norm;
      for(int j=0; j<LH_remez_arg[i].degree; j++)
      {
        int pole_idx = LH_remez_arg[i].degree - j - 1; // reverse order so we are forecasting on largest shifts
        Float gamma = 1.0 / ( 1.0 + LH_remez_arg[i].pole[pole_idx] );
        norm += LH_remez_arg[i].residue[pole_idx] * gamma;
        lat.SetEOFAParams(num_mass[i], num_mass[i], den_mass[i], -gamma, -1);
        if(heatbath_forecast) {
          lat.FminResExt(CG_solns_LH[j], frmn_tmp[1], CG_solns_LH, vm_LH, j, num_mass[i], num_mass[i], den_mass[i], -gamma, -1, LH_cg_arg_heatbath[i][pole_idx], CNV_FRM_YES);
        } else {
          CG_solns_LH[j] -> VecZero(f_size);
        }
//if(veloc_v<1)
        lat.FmatEvlMeofa(CG_solns_LH[j], frmn_tmp[1], LH_cg_arg_heatbath[i][pole_idx], num_mass[i], num_mass[i], den_mass[i], -gamma, -1);
        frmn_tmp[2] -> FTimesV1PlusV2(LH_remez_arg[i].residue[pole_idx]*gamma*gamma*kt, CG_solns_LH[j], frmn_tmp[2], f_size);
      }
      lat.Omega(frmn_tmp[1], frmn_tmp[2], -1, 1);
      lat.ChiralProj(frmn_tmp[2], frmn_tmp[1], -1);
    dtime +=dclock(true);
    print_flops(fname,"LHterm",0,dtime);
    dtime =-dclock(true);

      // RH term
      lat.ChiralProj(frmn_tmp[0], phi[i], 1);
      frmn_tmp[2] -> FTimesV1PlusV2(norm, phi[i], frmn_tmp[2], f_size);
      phi[i] -> VecZero(f_size);
      lat.Omega(frmn_tmp[1], frmn_tmp[0], 1, 0);
      for(int j=0; j<RH_remez_arg[i].degree; j++)
      {
        Float gamma = 1.0 / ( 1.0 + RH_remez_arg[i].pole[j] );
        lat.SetEOFAParams(den_mass[i], num_mass[i], den_mass[i], -gamma*RH_remez_arg[i].pole[j], 1);
        if(heatbath_forecast) {
          lat.FminResExt(CG_solns_RH[j], frmn_tmp[1], CG_solns_RH, vm_RH, j, den_mass[i], num_mass[i], den_mass[i], -gamma*RH_remez_arg[i].pole[j], 1, RH_cg_arg_heatbath[i][j], CNV_FRM_YES);
        } else {
          CG_solns_RH[j] -> VecZero(f_size);
        }
//if(veloc_v<1)
        lat.FmatEvlMeofa(CG_solns_RH[j], frmn_tmp[1], RH_cg_arg_heatbath[i][j], den_mass[i], num_mass[i], den_mass[i], -gamma*RH_remez_arg[i].pole[j], 1);
        phi[i] -> FTimesV1PlusV2(-RH_remez_arg[i].residue[j]*gamma*gamma*kt, CG_solns_RH[j], phi[i], f_size);
      }
      lat.Omega(frmn_tmp[1], phi[i], 1, 1);
      lat.ChiralProj(phi[i], frmn_tmp[1], 1);
    dtime +=dclock(true);
    print_flops(fname,"RHterm",0,dtime);
    dtime =-dclock(true);

      phi[i] -> FTimesV1PlusV2(1.0, frmn_tmp[2], phi[i], f_size);
    dtime +=dclock(true);
    print_flops(fname,"FTimesV1PlusV2()",0,dtime);
    dtime =-dclock(true);

      for(int j=0; j<LH_remez_arg[i].degree; j++){ updateCgStats(LH_cg_arg_heatbath[i][j]); }
      for(int j=0; j<RH_remez_arg[i].degree; j++){ updateCgStats(RH_cg_arg_heatbath[i][j]); }

      for(int j=LH_remez_arg[i].degree-1; j>=0; j--){
        sfree(CG_solns_LH[j], "CG_solns_LH[j]", fname, cname);
        sfree(vm_LH[j], "vm_LH[j]", fname, cname);
      }
      for(int j=RH_remez_arg[i].degree-1; j>=0; j--){
        sfree(CG_solns_RH[j], "CG_solns_RH[j]", fname, cname);
        sfree(vm_RH[j], "vm_RH[j]", fname, cname);
      }
      sfree(vm_RH, "vm_RH", fname, cname);
      sfree(vm_LH, "vm_LH", fname, cname);
      sfree(CG_solns_RH, "CG_solns_RH", fname, cname);
      sfree(CG_solns_LH, "CG_solns_LH", fname, cname);

      // Sanity check: does the action of the new pseudofermion field 
      // agree with the norm of the noise vector used to seed it?
      // This should fail if the rational approximation was not tuned correctly.
//if(veloc_v<1)
      if(heatbath_test)
      {
        // Zero tmp vectors
        for(int j=0; j<3; j++){ frmn_tmp[j] -> VecZero(f_size); }

        lat.SetEOFAParams(num_mass[i], 0.0, 0.0, 0.0, 1);

        Float h_i = lat.FhamiltonNode(phi[i], phi[i]);

        // LH term
        lat.ChiralProj(frmn_tmp[0], phi[i], -1);
        lat.Omega(frmn_tmp[1], frmn_tmp[0], -1, 0);
        frmn_tmp[2] -> VecZero(f_size);
        lat.FmatEvlMeofa(frmn_tmp[2], frmn_tmp[1], LH_cg_arg_mc[i], num_mass[i], 0.0, 0.0, 0.0, -1);
        lat.Omega(frmn_tmp[1], frmn_tmp[2], -1, 1);
        frmn_tmp[2] -> VecEqualsVecTimesEquFloat(frmn_tmp[1], -kt, f_size);
        h_i += lat.FhamiltonNode(frmn_tmp[0], frmn_tmp[2]);

        // RH term
        lat.ChiralProj(frmn_tmp[0], phi[i], 1);
        lat.Omega(frmn_tmp[1], frmn_tmp[0], 1, 0);
        frmn_tmp[2] -> VecZero(f_size);
        lat.FmatEvlMeofa(frmn_tmp[2], frmn_tmp[1], RH_cg_arg_mc[i], den_mass[i], num_mass[i], den_mass[i], -1.0, 1);
        lat.Omega(frmn_tmp[1], frmn_tmp[2], 1, 1);
        frmn_tmp[2] -> VecEqualsVecTimesEquFloat(frmn_tmp[1], kt, f_size);
        h_i += lat.FhamiltonNode(frmn_tmp[0], frmn_tmp[2]);

        // updateCgStats(LH_cg_arg_mc[i]);
        // updateCgStats(RH_cg_arg_mc[i]);

        glb_sum(&h_i);

        Float rel_err = std::fabs(gsum_h - h_i) / std::fabs(gsum_h);

        VRB.Result(cname, fname, "mass ratio %f/%f: heatbath test rel. error = %1.4e\n", 
            num_mass[i], den_mass[i], rel_err);
        
        Float stop_rsd = std::max(LH_cg_arg_mc[i]->stop_rsd, RH_cg_arg_mc[i]->stop_rsd); // use sloppier
        if(rel_err > 100.0 * stop_rsd) {
          ERR.General(cname, fname, "heatbath test failed. rel_error(%e) >100*stop_rsd(%e)\n",rel_err,stop_rsd);
        }
      }
    }

    LatticeFactory::Destroy();

    evolved = 0;
    heatbathEval = 1;
    energyEval = 0;
    dtime +=dclock(true);
    print_flops(fname,"eigen_measure()",0,dtime);
    dtime =-dclock(true);
//    traj++;
  }

  dtime2 += dclock();
  print_flops(cname, fname, 0, dtime2);
  timer.stop(true);

  GJP.DisableEOFA();
}

void AlgActionEOFA::prepare_fg(Matrix* force, Float dt_ratio)
{
  char fname[30 + strlen(force_label)];
  sprintf(fname, "prepare_fg(M*,F)[%s]", force_label);

  static Timer timer(cname, fname);
  timer.start(true);
  Float dtime = -dclock();
  Float dtime_cg(0.0);
  Float dtime_force(0.0);

  if(skip_force)
  {
    VRB.Result(cname, fname, "WARNING: skipping prepare_fg() because AlgActionEOFA::skip_force is true!\n");
    evolved = 1;
    timer.stop(true);
    return;
  }

  GJP.EnableEOFA();

  Lattice& lat = LatticeFactory::Create(fermion, G_CLASS_NONE);

  for(int i=0; i<n_masses; i++)
  {
    lat.SetEOFAParams(num_mass[i], 0.0, 0.0, 0.0, -1);
    Float kt = lat.kt(num_mass[i], den_mass[i]);

    dtime_cg -= dclock();

    for(int j=0; j<3; j++){ frmn_tmp[j] -> VecZero(f_size); }

    // Zero FG forecasting vectors
    frmn_tmp[2*i+3] -> VecZero(f_size);
    frmn_tmp[2*i+4] -> VecZero(f_size);

    // LH term
    // \chi_{L} = Heofa(m_f)^{-1} * \Omega_{-} * P_{-} * \Phi
    lat.ChiralProj(frmn_tmp[0], phi[i], -1);
    lat.Omega(frmn_tmp[2], frmn_tmp[0], -1, 0);
    cg_iter = lat.FmatEvlMeofa(frmn_tmp[0], frmn_tmp[2], LH_cg_arg_fg[i], num_mass[i], 0.0, 0.0, 0.0, -1, frmn_tmp[2*i+3]);
    VRB.Result(cname, fname, "mass ratio %f/%f: cg_iter = %d (LH)\n", num_mass[i], den_mass[i], cg_iter);

    // RH term
    // \chi_{R} = ( Heofa(m_b) - \Delta_{+}(m_f,m_b)*P_{+} )^{-1} * \Omega_{+} * P_{+} * \Phi
    lat.ChiralProj(frmn_tmp[1], phi[i], 1);
    lat.Omega(frmn_tmp[0], frmn_tmp[1], 1, 0);
    cg_iter = lat.FmatEvlMeofa(frmn_tmp[1], frmn_tmp[0], RH_cg_arg_fg[i], den_mass[i], num_mass[i], den_mass[i], -1.0, 1, frmn_tmp[2*i+4]);
    VRB.Result(cname, fname, "mass ratio %f/%f: cg_iter = %d (RH)\n", num_mass[i], den_mass[i], cg_iter);

    dtime_cg += dclock();
    updateCgStats(LH_cg_arg_fg[i]);
    updateCgStats(RH_cg_arg_fg[i]);

    Matrix* mom_tmp = force;
    if(force_measure == FORCE_MEASURE_YES){
      mom_tmp = static_cast<Matrix*>( smalloc(g_size*sizeof(Float), cname, fname, "mom_tmp") );
      reinterpret_cast<Vector*>(mom_tmp) -> VecZero(g_size);
    }

    dtime_force -= dclock();

    // Evolve momentum using force from LH term
    Fdt = lat.EvolveMomFforce(mom_tmp, frmn_tmp[0], num_mass[i], 0.5*kt*dt_ratio);
    if(force_measure == FORCE_MEASURE_YES) {
      char label[200];
      sprintf(label, "%s (LH), mass ratio %f/%f:", force_label, num_mass[i], den_mass[i]);
      Fdt.print(dt_ratio, label);
    }

    // Evolve momentum using force from RH term
    Fdt = lat.EvolveMomFforce(mom_tmp, frmn_tmp[1], den_mass[i], -0.5*kt*dt_ratio);
    dtime_force += dclock();
    if(force_measure == FORCE_MEASURE_YES) {
      char label[200];
      sprintf(label, "%s (RH), mass ratio %f/%f:", force_label, num_mass[i], den_mass[i]);
      Fdt.print(dt_ratio, label);
    }

    // Monitor total force contribution
    if(force_measure == FORCE_MEASURE_YES)
    {
      Fdt.measure(mom_tmp);
      Fdt.glb_reduce();

      reinterpret_cast<Vector*>(force) -> VecAddEquVec(reinterpret_cast<Vector*>(mom_tmp), g_size);

      char label[200];
      sprintf(label, "%s (total), mass ratio %f/%f:", force_label, num_mass[i], den_mass[i]);
      Fdt.print(dt_ratio, label);

      sfree(mom_tmp, "mom_tmp", fname, cname);
    }
  }

  fg_forecast = true; // We now have solutions to forecast the next CG solve
#ifdef USE_QUDA
  fg_forecast = false; // CJ:DtildeInv() is slow at the moment.
#endif

  LatticeFactory::Destroy();

  dtime += dclock();

  char fname_cg[30 + strlen(force_label)];
  sprintf(fname_cg, "prepare_fg::cg()[%s]", force_label);
  char fname_force[30 + strlen(force_label)];
  sprintf(fname_force, "prepare_fg::force()[%s]", force_label);

  print_flops(cname, fname, 0, dtime);
  print_flops(cname, fname_cg, 0, dtime_cg);
  print_flops(cname, fname_force, 0, dtime_force);
  timer.stop(true);

  GJP.DisableEOFA();
}

void AlgActionEOFA::evolve(Float dt, int nsteps)
{
  char fname[30 + strlen(force_label)];
  sprintf(fname, "evolve(F,i)[%s]", force_label);
  
  static Timer timer(cname, fname);
  timer.start(true);
  Float dtime = -dclock();
  Float dtime_cg(0.0);
  Float dtime_force(0.0);

  if(skip_force){
    VRB.Result(cname, fname, "WARNING: skipping evolve() because AlgActionEOFA::skip_force is true\n");
    evolved = 1;
    heatbathEval = 0;
    energyEval = 0;
    timer.stop(true);
    return;
  }

  GJP.EnableEOFA();

  Lattice& lat = LatticeFactory::Create(fermion, G_CLASS_NONE);

  for(int steps=0; steps<nsteps; steps++)
  {
    for(int i=0; i<n_masses; i++)
    {
      if(UniqueID() == 0){ printf("AlgActionEOFA::evolve()[%s] step %d mass %d\n", force_label, steps, i); }

      lat.SetEOFAParams(num_mass[i], 0.0, 0.0, 0.0, -1);
      Float kt = lat.kt(num_mass[i], den_mass[i]);

      dtime_cg -= dclock();

      for(int j=0; j<3; j++){ frmn_tmp[j] -> VecZero(f_size); }

      if(!fg_forecast){ // Zero solution vectors if we don't have an initial guess
        frmn_tmp[2*i+3] -> VecZero(f_size);
        frmn_tmp[2*i+4] -> VecZero(f_size);
      } else {          // Otherwise we use forecasted guesses
        VRB.Result(cname, fname, "Using force gradient forecasting.\n");
        // lat.DtildeInv(frmn_tmp[2*i+3], num_mass[i]);
        // lat.DtildeInv(frmn_tmp[2*i+4], den_mass[i]);
      }

      // LH term
      // \chi_{L} = Heofa(m_f)^{-1} * \Omega_{-} * P_{-} * \Phi
      lat.ChiralProj(frmn_tmp[0], phi[i], -1);
      lat.Omega(frmn_tmp[2], frmn_tmp[0], -1, 0);
      cg_iter = lat.FmatEvlMeofa(frmn_tmp[2*i+3], frmn_tmp[2], LH_cg_arg_md[i], num_mass[i], 0.0, 0.0, 0.0, -1);
      VRB.Result(cname, fname, "mass ratio %f/%f: cg_iter = %d (LH)\n");

      // RH term
      // \chi_{R} = ( Heofa(m_b) - \Delta_{+}(m_f,m_b)*P_{+} )^{-1} * \Omega_{+} * P_{+} * \Phi
      lat.ChiralProj(frmn_tmp[2], phi[i], 1);
      lat.Omega(frmn_tmp[0], frmn_tmp[2], 1, 0);
      cg_iter = lat.FmatEvlMeofa(frmn_tmp[2*i+4], frmn_tmp[0], RH_cg_arg_md[i], den_mass[i], num_mass[i], den_mass[i], -1.0, 1);
      VRB.Result(cname, fname, "mass ratio %f/%f: cg_iter = %d (RH)\n");

      dtime_cg += dclock();
      updateCgStats(LH_cg_arg_md[i]);
      updateCgStats(RH_cg_arg_md[i]);

      Matrix* mom_tmp;
      if(force_measure == FORCE_MEASURE_YES){
        mom_tmp = static_cast<Matrix*>( smalloc(g_size*sizeof(Float), cname, fname, "mom_tmp") );
        reinterpret_cast<Vector*>(mom_tmp) -> VecZero(g_size);
      } else {
        mom_tmp = mom;
      }

      dtime_force -= dclock();

      // Evolve momentum using force from LH spinor
      Fdt = lat.EvolveMomFforce(mom_tmp, frmn_tmp[2*i+3], num_mass[i], 0.5*kt*dt);
      if(force_measure == FORCE_MEASURE_YES) {
        char label[200];
        sprintf(label, "%s (LH), mass ratio %f/%f:", force_label, num_mass[i], den_mass[i]);
        Fdt.print(dt, label);
      }

      // Evolve momentum using force from RH spinor
      Fdt = lat.EvolveMomFforce(mom_tmp, frmn_tmp[2*i+4], den_mass[i], -0.5*kt*dt);
      dtime_force += dclock();
      if(force_measure == FORCE_MEASURE_YES) {
        char label[200];
        sprintf(label, "%s (RH), mass ratio %f/%f:", force_label, num_mass[i], den_mass[i]);
        Fdt.print(dt, label);
      }

      // Monitor total force contribution
      if(force_measure == FORCE_MEASURE_YES) 
      {
        Fdt.measure(mom_tmp);
        Fdt.glb_reduce();

        fTimesV1PlusV2(reinterpret_cast<IFloat*>(mom), 1.0, reinterpret_cast<IFloat*>(mom_tmp), reinterpret_cast<Float*>(mom), g_size);
        sfree(mom_tmp, "mom_tmp", fname, cname);

        char label[200];
        sprintf(label, "%s (total), mass ratio %f/%f:", force_label, num_mass[i], den_mass[i]);
        Fdt.print(dt, label);
      }

      if(UniqueID() == 0){
        printf("AlgActionEOFA::evolve() [%s] end of step %d, mass ratio %f/%f\n", force_label, steps, num_mass[i], den_mass[i]);
      }
    }

    fg_forecast = false;
    evolved = 1;
    heatbathEval = 0;
    energyEval = 0;
    md_steps++;
  }

  if(UniqueID() == 0){
    printf("AlgActionEOFA::evolve() [%s] end\n", force_label);
  }

  LatticeFactory::Destroy();

  char fname_cg[30 + strlen(force_label)];
  sprintf(fname_cg, "evolve::cg()[%s]", force_label);
  char fname_force[30 + strlen(force_label)];
  sprintf(fname_force, "evolve::force()[%s]", force_label);

  dtime += dclock();
  print_flops(cname, fname, 0, dtime);
  print_flops(cname, fname_cg, 0, dtime_cg);
  print_flops(cname, fname_force, 0, dtime_force);
  timer.stop(true);

  GJP.DisableEOFA();
}

void AlgActionEOFA::generateApprox(Float* mass, RemezArg** remez_arg, EOFARationalDescr* rat)
{
  const char* fname = "generateApprox(F*,RemezArg**,EOFARationalDescr*)";

  *remez_arg = new RemezArg[n_masses];

  for(int i=0; i<n_masses; i++) 
  {
    (*remez_arg)[i].approx_type = rat[i].rat_approx.approx_type;
    
    switch(rat[i].rat_approx.bounds_type) 
    {
      case RATIONAL_BOUNDS_MANUAL:
        (*remez_arg)[i].lambda_low  = rat[i].rat_approx.lambda_low;
        (*remez_arg)[i].lambda_high = rat[i].rat_approx.lambda_high;
        break;
      case RATIONAL_BOUNDS_AUTOMATIC:
        ERR.General(cname, fname, "RationalApproxType %d not implemented for FclassType %d\n",
            RATIONAL_BOUNDS_AUTOMATIC, fermion);
        break;
      default:
        ERR.General(cname, fname, "RationalApproxType %d not implemented\n", rat[i].rat_approx.bounds_type);
    }

    (*remez_arg)[i].degree       = rat[i].rat_approx.stop_rsd.stop_rsd_len;
    (*remez_arg)[i].field_type   = rat[i].field_type;
    (*remez_arg)[i].power_num    = rat[i].power_num;
    (*remez_arg)[i].power_den    = rat[i].power_den;
    (*remez_arg)[i].precision    = rat[i].precision;
    (*remez_arg)[i].valid_approx = 0;
    (*remez_arg)[i].delta_m      = -(*remez_arg)[i].lambda_low;
  }

  // Construct approximations
  for(int i=0; i<n_masses; i++) {
    for(int j=0; j<i; j++){ (*remez_arg)[i].valid_approx = compareApprox((*remez_arg)[i], (*remez_arg)[j]); }
    if(!(*remez_arg)[i].valid_approx) {
      AlgRemez remez((*remez_arg)[i]);
      remez.generateApprox();
      (*remez_arg)[i].valid_approx = 1;
    }
  }
}

void AlgActionEOFA::destroyApprox(RemezArg* remez_arg)
{
  const char* fname = "destroyApprox(RemezArg*)";
  delete[] remez_arg;
}

int AlgActionEOFA::compareApprox(RemezArg& arg1, RemezArg& arg2)
{
  int result(0);

  if( ( arg1.field_type  == arg2.field_type  ) &&
      ( arg1.power_num   == arg2.power_num   ) &&
      ( arg1.power_den   == arg2.power_den   ) &&
      ( arg1.lambda_low  == arg2.lambda_low  ) &&
      ( arg1.lambda_high == arg2.lambda_high ) &&
      ( arg1.precision   == arg2.precision   ) &&
      ( arg1.degree      == arg2.degree      ) )
  {
    arg1.norm     = arg2.norm;
    arg1.norm_inv = arg2.norm_inv;

    for(int k=0; k<arg1.degree; k++) {
      arg1.residue[k]     = arg2.residue[k];
      arg1.pole[k]        = arg2.pole[k];
      arg1.residue_inv[k] = arg2.residue_inv[k];
      arg1.pole_inv[k]    = arg2.pole_inv[k];
    }

    result = 1;
  }

  return result;
}

void AlgActionEOFA::generateCgArg(Float* mass, CgArg**** cg_arg, const char* l, EOFARationalDescr* rat)
{
  const char* fname = "generateCgArg(F*,CgArg****,char*,EOFARationalDescr*)";

  char label[100], label_i[100], label_ij[100];
  sprintf(label, "%s", l);
  sprintf(label_i, "%s[i]", l);
  sprintf(label_ij, "s[i][j]", l);

  (*cg_arg) = static_cast<CgArg***>( smalloc(n_masses*sizeof(CgArg**), label, fname, cname) );

  for(int i=0; i<n_masses; i++)
  {
    (*cg_arg)[i] = static_cast<CgArg**>( smalloc(rat[i].rat_approx.stop_rsd.stop_rsd_len*sizeof(CgArg*), label_i, fname, cname) );
    
    for(int j=0; j<rat[i].rat_approx.stop_rsd.stop_rsd_len; j++) 
    {
//      (*cg_arg)[i][j]               = static_cast<CgArg*>( smalloc(sizeof(CgArg), label_ij, fname, cname) );
      (*cg_arg)[i][j]               = new CgArg;
      (*cg_arg)[i][j]->mass         = mass[i];
      (*cg_arg)[i][j]->max_num_iter = max_num_iter[i];
      (*cg_arg)[i][j]->stop_rsd     = rat[i].rat_approx.stop_rsd.stop_rsd_val[j];
    }
  }
}

void AlgActionEOFA::generateCgArg(Float* mass, CgArg*** cg_arg_fg, CgArg*** cg_arg_mc, 
    CgArg*** cg_arg_md, const char* l, Float* stop_rsd_fg, Float* stop_rsd_mc, Float* stop_rsd_md)
{
  const char* fname = "generateCgArg(F*,CA***,CA***,CA***,char*,F*,F*,F*)";

  char fg_label[100], fg_label_i[100];
  char mc_label[100], mc_label_i[100];
  char md_label[100], md_label_i[100];
  sprintf(fg_label, "%s_fg", l);
  sprintf(fg_label_i, "%s_fg[i]", l);
  sprintf(mc_label, "%s_mc", l);
  sprintf(mc_label_i, "%s_mc[i]", l);
  sprintf(md_label, "%s_md", l);
  sprintf(md_label_i, "%s_md[i]", l);

  (*cg_arg_fg) = static_cast<CgArg**>( smalloc(n_masses*sizeof(CgArg*), fg_label, fname, cname) );
  (*cg_arg_mc) = static_cast<CgArg**>( smalloc(n_masses*sizeof(CgArg*), mc_label, fname, cname) );
  (*cg_arg_md) = static_cast<CgArg**>( smalloc(n_masses*sizeof(CgArg*), md_label, fname, cname) );

  for(int i=0; i<n_masses; i++)
  {
//    (*cg_arg_fg)[i]               = static_cast<CgArg*>( smalloc(sizeof(CgArg), fg_label_i, fname, cname) );
    (*cg_arg_fg)[i]              = new CgArg;
    (*cg_arg_fg)[i]->mass         = mass[i];
    (*cg_arg_fg)[i]->max_num_iter = max_num_iter[i];
    (*cg_arg_fg)[i]->stop_rsd     = stop_rsd_fg[i];

//    (*cg_arg_mc)[i]               = static_cast<CgArg*>( smalloc(sizeof(CgArg), mc_label_i, fname, cname) );
    (*cg_arg_mc)[i]              = new CgArg;
    (*cg_arg_mc)[i]->mass         = mass[i];
    (*cg_arg_mc)[i]->max_num_iter = max_num_iter[i];
    (*cg_arg_mc)[i]->stop_rsd     = stop_rsd_mc[i];

//    (*cg_arg_md)[i]               = static_cast<CgArg*>( smalloc(sizeof(CgArg), md_label_i, fname, cname) );
    (*cg_arg_md)[i]              = new CgArg;
    (*cg_arg_md)[i]->mass         = mass[i];
    (*cg_arg_md)[i]->max_num_iter = max_num_iter[i];
    (*cg_arg_md)[i]->stop_rsd     = stop_rsd_md[i];
  }
}

void AlgActionEOFA::destroyCgArg(CgArg*** cg_arg, const char* l, RemezArg* remez_arg)
{
  const char* fname = "destroyCgArg(CgArg***,char*,RemezArg*)";
  char label[100], label_i[100], label_ij[100];
  sprintf(label, "%s", l);
  sprintf(label_i, "%s[i]", l);
  sprintf(label_ij, "%s[i][j]", l);

  for(int i=0; i<n_masses; i++) {
//    for(int j=0; j<remez_arg[i].degree; j++){ sfree(cg_arg[i][j], label_ij, fname, cname); }
    for(int j=0; j<remez_arg[i].degree; j++){  delete cg_arg[i][j]; }
    sfree(cg_arg[i], label_i, fname, cname);
  }

  sfree(cg_arg, label, fname, cname);
}

void AlgActionEOFA::destroyCgArg(CgArg** cg_arg_fg, CgArg** cg_arg_mc, CgArg** cg_arg_md, const char* l)
{
  const char* fname = "destroyCgArg(CgArg**,CgArg**,CgArg**,char*)";
  char fg_label[100], fg_label_i[100];
  char mc_label[100], mc_label_i[100];
  char md_label[100], md_label_i[100];
  sprintf(fg_label, "%s_fg", l);
  sprintf(fg_label_i, "%s_fg[i]", l);
  sprintf(mc_label, "%s_mc", l);
  sprintf(mc_label_i, "%s_mc[i]", l);
  sprintf(md_label, "%s_md", l);
  sprintf(md_label_i, "%s_md[i]", l);

  for(int i=0; i<n_masses; i++) {
//    sfree(cg_arg_fg[i], fg_label_i, fname, cname);
//    sfree(cg_arg_mc[i], mc_label_i, fname, cname);
//    sfree(cg_arg_md[i], md_label_i, fname, cname);
     delete cg_arg_fg[i];
     delete cg_arg_mc[i];
     delete cg_arg_md[i];
  }
  
  sfree(cg_arg_fg, fg_label, fname, cname);
  sfree(cg_arg_mc, mc_label, fname, cname);
  sfree(cg_arg_md, md_label, fname, cname);
}

void AlgActionEOFA::generateEigArg(const EigenDescr& eigen)
{
  const char* fname = "generateEigArg(EigenDescr)";

  eig_arg.pattern_kind  = ARRAY;
  eig_arg.Mass.Mass_len = n_masses;
  eig_arg.Mass.Mass_val = static_cast<Float*>( smalloc(n_masses*sizeof(Float), "Mass_val", fname, cname) );

  eig_arg.N_eig      = 1;
  eig_arg.Kalk_Sim   = 0;
  eig_arg.MaxCG      = eigen.max_num_iter;
  eig_arg.RsdR_a     = eigen.stop_rsd;
  eig_arg.RsdR_r     = eigen.stop_rsd;
  eig_arg.Rsdlam     = eigen.stop_rsd;
  eig_arg.Cv_fact    = 0.0;
  eig_arg.N_min      = 0;
  eig_arg.N_max      = 0;
  eig_arg.N_KS_max   = 0;
  eig_arg.n_renorm   = 100;
  eig_arg.ProjApsiP  = 0;
  eig_arg.print_hsum = 0;
  eig_arg.hsum_dir   = 0;
  eig_arg.ncorr      = 0;
  eig_arg.fname      = 0;

  lambda_low  = static_cast<Float**>( smalloc(eig_arg.N_eig*sizeof(Float*), "lambda_low", fname, cname) );
  lambda_high = static_cast<Float**>( smalloc(eig_arg.N_eig*sizeof(Float*), "lambda_high", fname, cname) );
  for(int i=0; i<eig_arg.N_eig; i++){
    lambda_low[i]  = static_cast<Float*>( smalloc(n_masses*sizeof(Float), "lambda_low[i]", fname, cname) );
    lambda_high[i] = static_cast<Float*>( smalloc(n_masses*sizeof(Float), "lambda_high[i]", fname, cname) );
  }
}

void AlgActionEOFA::destroyEigArg()
{
  const char* fname = "destroyEigArg()";

  for(int i=0; i<eig_arg.N_eig; i++) {
    sfree(lambda_low[i], "lambda_low[i]", fname, cname);
    sfree(lambda_high[i], "lambda_high[i]", fname, cname);
  }
  sfree(lambda_low, "lambda_low", fname, cname);
  sfree(lambda_high, "lambda_high", fname, cname);

  sfree(eig_arg.Mass.Mass_val, "Mass_val", fname, cname);
}

bool AlgActionEOFA::checkPolesFile(const RemezArg& ra, const EOFARationalDescr& r)
{
  if( ra.field_type != r.field_type                       ) { return false; }
  if( ra.power_num  != r.power_num                        ) { return false; }
  if( ra.degree     != r.rat_approx.stop_rsd.stop_rsd_len ) { return false; }
  
  double r_lo = fabs( ra.lambda_low - r.rat_approx.lambda_low ) / fabs( r.rat_approx.lambda_low );
  double r_hi = fabs( ra.lambda_high - r.rat_approx.lambda_high ) / fabs( r.rat_approx.lambda_high );
  if((r_lo > 1.0e-03) || (r_hi > 1.0e-03)){ return false; }

  return true;
}

bool AlgActionEOFA::loadPoles()
{
  const char* fname = "loadPoles()";
  
  if(eofa_arg->remez_generate){ return false; }
  
  if(strlen(eofa_arg->rat_poles_file) == 0){ return false; }

  FILE* fp = fopen(eofa_arg->rat_poles_file, "r");
  if(fp == nullptr){ return false; }
  fclose(fp);

  EOFARemezArg ea;
  if(!ea.Decode(eofa_arg->rat_poles_file, "ea")){ return false; }

  if(ea.LH.LH_len != n_masses){ return false; }
  if(ea.RH.RH_len != n_masses){ return false; }

  LH_remez_arg = new RemezArg[n_masses];
  RH_remez_arg = new RemezArg[n_masses];

  for(int i=0; i<n_masses; i++)
  {
    LH_remez_arg[i] = ea.LH.LH_val[i];
    RH_remez_arg[i] = ea.RH.RH_val[i];

    if(!checkPolesFile(LH_remez_arg[i], eofa_arg->LH_rat_approx.LH_rat_approx_val[i])){ return false; }
    if(!checkPolesFile(RH_remez_arg[i], eofa_arg->RH_rat_approx.RH_rat_approx_val[i])){ return false; }
  }

  VRB.Result(cname, fname, "Successfully loaded poles file %s.\n", eofa_arg->rat_poles_file);
  return true;
} 

bool AlgActionEOFA::savePoles()
{
  const char* fname = "savePoles()";

  if(strlen(eofa_arg->rat_poles_file) == 0){ return false; }

  EOFARemezArg ea;

  ea.LH.LH_len = n_masses;
  ea.RH.RH_len = n_masses;

  ea.LH.LH_val = LH_remez_arg;
  ea.RH.RH_val = RH_remez_arg;

  return ea.Encode(eofa_arg->rat_poles_file, "ea");
}

AlgActionEOFA::~AlgActionEOFA()
{
  const char* fname = "~AlgActionEOFA()";
  VRB.Func(cname, fname);

  if(n_masses > 0)
  {
    // Free temporary pseudofermion fields
    for(int i=2*n_masses+2; i>=0; i--){
      sfree(frmn_tmp[i], "frmn_tmp[i]", fname, cname);
    }
    sfree(frmn_tmp, "frmn_tmp", fname, cname);

    // Free memory for the fermion CG arguments
    destroyCgArg(LH_cg_arg_heatbath, "LH_cg_arg_heatbath", LH_remez_arg);
    destroyCgArg(LH_cg_arg_fg, LH_cg_arg_mc, LH_cg_arg_md, "LH_cg_arg");
    destroyCgArg(RH_cg_arg_heatbath, "RH_cg_arg_heatbath", RH_remez_arg);
    destroyCgArg(RH_cg_arg_fg, RH_cg_arg_mc, RH_cg_arg_md, "RH_cg_arg");

    // Must not free these until CG args are freed
    destroyApprox(LH_remez_arg);
    destroyApprox(RH_remez_arg);

    sfree(num_mass, "num_mass", fname, cname);
    sfree(den_mass, "den_mass", fname, cname);

    if(eofa_arg->eigen.eigen_measure == EIGEN_MEASURE_YES){ destroyEigArg(); }

    // Delete and reallocate phi fields
    for(int i=0; i<this->n_masses; i++) {
      sfree(phi[i], cname, fname, "phi[i]");
      phi[i] = static_cast<Vector*>( smalloc(this->f_size*sizeof(Float)/this->Ncb, "phi[i]", fname, cname) );
    }
  }
}

CPS_END_NAMESPACE
