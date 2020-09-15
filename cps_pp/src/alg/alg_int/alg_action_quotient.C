#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_action_quotient.C
//
// AlgActionQuotient is a class describing a bilinear action, with a
// kernel given as a matrix quotient, i.e.,
//
//   S = phi^dagger M_1 (M_2^dagger M_2)^{-1} M_1^dagger phi.
//
// This action covers 2 flavour domain wall theories with the
// Pauli-Villars cancellation and also Hasenbusch type actions.  The
// mass parameters for M_1 and M_2 are specified separately.
//
// This action also uses a chronological predictor to reduce the
// number of cg iterations.
//------------------------------------------------------------------
//CK: In order to directly reproduce the Quotient action results for WilsonTm fermions from CPS using Fbfm, we need to apply a unitary
//    transformation to the random vectors due to the differing normalization choices. If disabled the Fbfm results will differ slightly from the CPS results
#define CPS_FBFM_WILSONTM_COMPAT_MODE
  CPS_END_NAMESPACE
#include<alg/alg_hmd.h>
#include<util/lattice.h>
#include<util/vector.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/time_cps.h>
#include<alg/alg_int.h>
#include<util/dirac_op.h>
#include <util/timer.h>
#include <vector>
 
CPS_START_NAMESPACE
AlgActionQuotient::AlgActionQuotient (AlgMomentum & mom,
					ActionQuotientArg & q_arg)
:AlgActionBilinear (mom, q_arg.bi_arg), cname ("AlgActionQuotient")
//      skip_force(false)
{
  const char *fname = "AlgActionQuotient()";

  int_type = INT_QUOTIENT;
  quo_arg = &q_arg;

  //!< First check n_masses bilinear = n_masses quotient
  if (quo_arg->quotients.quotients_len !=
      quo_arg->bi_arg.bilinears.bilinears_len)
    ERR.General (cname, fname,
		 "Inconsistency between QuotientArg and BilinearArg n_masses\n");

  if (n_masses > 0) {

    bsn_cg_arg.resize (n_masses);
    frm_cg_arg_fg.resize (n_masses);
    frm_cg_arg_md.resize (n_masses);
    frm_cg_arg_mc.resize (n_masses);

    //!< Initialize the CG arguments
    for (int i = 0; i < n_masses; i++) {
      const QuotientDescr & qi = quo_arg->quotients.quotients_val[i];

      bsn_mass.push_back (qi.bsn_mass);
      frm_mass.push_back (qi.frm_mass);


      //~~ added for twisted mass Wilson fermions
      bsn_mass_epsilon.push_back (qi.bsn_mass_epsilon);
      frm_mass_epsilon.push_back (qi.frm_mass_epsilon);

      bsn_cg_arg[i].mass = qi.bsn_mass;
      //~~ added for twisted mass Wilson fermions
      bsn_cg_arg[i].epsilon = qi.bsn_mass_epsilon;
      bsn_cg_arg[i].max_num_iter = max_num_iter[i];
      bsn_cg_arg[i].stop_rsd = qi.stop_rsd_hb;

      frm_cg_arg_md[i].mass = qi.frm_mass;
      //~~ added for twisted mass Wilson fermions
      frm_cg_arg_md[i].epsilon = qi.frm_mass_epsilon;
      frm_cg_arg_md[i].max_num_iter = max_num_iter[i];
      frm_cg_arg_md[i].stop_rsd = qi.stop_rsd_md;

      frm_cg_arg_fg[i] = frm_cg_arg_md[i];
      frm_cg_arg_fg[i].stop_rsd = qi.stop_rsd_md * qi.stop_rsd_fg_mult;

      frm_cg_arg_mc[i] = frm_cg_arg_md[i];
      frm_cg_arg_mc[i].stop_rsd = qi.stop_rsd_mc;

      chrono.push_back (qi.chrono);
    }

    evolved = 1;

    //!< Vectors used to store solution history
    v = (Vector ***) smalloc (n_masses * sizeof (Vector **), "v", fname, cname);
    cg_sol_old = (Vector ***) smalloc (n_masses * sizeof (Vector **),
				       "cg_sol_old", fname, cname);
    vm = (Vector ***) smalloc (n_masses * sizeof (Vector **),
			       "vm", fname, cname);

    tmp1 = (Vector *) smalloc (2*f_size * sizeof (Float), "tmp1", fname, cname);
    tmp2 = (Vector *) smalloc (2*f_size * sizeof (Float), "tmp2", fname, cname);
    tmp3 = (Vector *) smalloc (2*f_size * sizeof (Float), "tmp3", fname, cname);

    VRB.Result (cname, fname, "allocating fermion fields of size %d Floats\n",
		f_size);

    for (int i = 0; i < n_masses; i++) {
      int deg = 0;
      if (chrono[i] > 0)
	deg = chrono[i];
      else if (chrono[i] == 0)
	deg = 1;
      else
	ERR.General (cname, fname, "Cannot have negative chronology\n");

      v[i] =
	(Vector **) smalloc (deg * sizeof (Vector *), "v[i]", fname, cname);
      cg_sol_old[i] =
	(Vector **) smalloc (deg * sizeof (Vector *), "cg_sol_old[i]", fname,
			     cname);
      vm[i] =
	(Vector **) smalloc (deg * sizeof (Vector *), "vm[i]", fname, cname);
      for (int j = 0; j < deg; j++) {
	v[i][j] = (Vector *) smalloc (f_size * sizeof (Float),
				      "v[i][j]", fname, cname);
	vm[i][j] = (Vector *) smalloc (f_size * sizeof (Float),
				       "vm[i][j]", fname, cname);
      }
    }

  }

  init ();

  fg_forecast = false;
}

void AlgActionQuotient::init ()
{
  AlgActionBilinear::init ();
  evolved = 1;
  for (int i = 0; i < n_masses; i++)
    v[i][0]->VecZero (f_size);
}

AlgActionQuotient::~AlgActionQuotient ()
{
  const char *fname = "~AlgActionQuotient()";

  if (n_masses > 0) {
    //!< Free chronology
    for (int i = 0; i < n_masses; i++) {
      int deg = 0;
      if (chrono[i] > 0)
	deg = chrono[i];
      else if (chrono[i] == 0)
	deg = 1;

      for (int j = 0; j < deg; j++) {
	sfree (vm[i][j], "vm[i][j]", fname, cname);
	sfree (v[i][j], "v[i][j]", fname, cname);
      }
      sfree (cg_sol_old[i], "cg_sol_old[i]", fname, cname);
      sfree (vm[i], "vm[i]", fname, cname);
      sfree (v[i], "v[i]", fname, cname);
    }
    sfree (tmp3, "tmp3", fname, cname);
    sfree (tmp2, "tmp2", fname, cname);
    sfree (tmp1, "tmp1", fname, cname);
    sfree (cg_sol_old, "cg_sol_old", fname, cname);
    sfree (vm, "vm", fname, cname);
    sfree (v, "v", fname, cname);
  }

}

//!< Heat Bath for quotients
void AlgActionQuotient::reweight (Float * rw_fac, Float * norm)
{

  char *fname = "reweight()";
  VRB.Func (cname, fname);

  if (n_masses > 0) {
    Lattice & lat = LatticeFactory::Create (fermion, G_CLASS_NONE);

    // tmp1, tmp2 < - random Gaussian vector (RGV)
    for (int i = 0; i < n_masses; i++) {
      lat.RandGaussVector (tmp1, 0.5, Ncb);
      lat.RandGaussVector (tmp2, 0.5, Ncb);


      {
	// everything besides Fdwf4d:

	//~~ changed for twisted mass Wilson fermions
	// phi <- M_b^\dag (RGV)
	norm[i] = (lat.Fclass () == F_CLASS_WILSON_TM) ?
	  lat.SetPhi (phi[i], tmp1, tmp2, bsn_mass[i], bsn_mass_epsilon[i], DAG_YES) : 
	  lat.SetPhi (phi[i], tmp1, tmp2, bsn_mass[i], DAG_YES);

	// tmp2 <- (M_f^\dag M_f)^{-1} M_b^\dag (RGV)
	tmp2->VecZero (f_size);
	cg_iter = lat.FmatEvlInv (tmp2, phi[i], &frm_cg_arg_mc[i], CNV_FRM_NO);

	rw_fac[i] = lat.FhamiltonNode (phi[i], tmp2);
	rw_fac[i] -= norm[i];
      }

      VRB.Result (cname, fname, "rw_fac=%e norm=%e\n", rw_fac[i], norm[i]);
      glb_sum (rw_fac + i);
      glb_sum (norm + i);

      updateCgStats (&bsn_cg_arg[i]);
    }

    LatticeFactory::Destroy ();
    //    evolved = 0;
  }

}

//!< Heat Bath for quotients
void AlgActionQuotient::heatbath ()
{

  char *fname = "heatbath()";
  VRB.Func (cname, fname);
  static Timer timer (cname, fname);
  timer.start ();
//  Float dtime = -dclock ();
#if 0
std::string veloc_label;
//veloc_label <<"Phi_traj"<<traj_num;
veloc_label ="Phi_traj"+std::to_string(traj_num);
int veloc_v = getVer(veloc_label.c_str());
VRB.Result(cname,fname,"label version=%s %d\n",veloc_label.c_str(),veloc_v);
#endif


  if (n_masses > 0) {
    Lattice & lat = LatticeFactory::Create (fermion, G_CLASS_NONE);

    h_init = 0.;
    VRB.Debug (cname, fname, "h_init=%0.14e\n", h_init);

    // tmp1, tmp2 < - random Gaussian vector (RGV)
    for (int i = 0; i < n_masses; i++) {
      lat.RandGaussVector (tmp1, 0.5, Ncb);
      lat.RandGaussVector (tmp2, 0.5, Ncb);


      Float h_i;
      {
	// everything besides Fdwf4d:

	//~~ changed for twisted mass Wilson fermions
	// phi <- M_f^\dag (RGV)
	VRB.Result (cname,fname,"|R>=%e\n", lat.FhamiltonNode (tmp1, tmp1));
	if (lat.Fclass () == F_CLASS_WILSON_TM ) {
	  h_i = lat.SetPhi (tmp3, tmp1, tmp2, frm_mass[i], frm_mass_epsilon[i], DAG_YES);
	} else {
	  h_i = lat.SetPhi (tmp3, tmp1, tmp2, frm_mass[i], DAG_YES);
	}
	h_init += h_i;
	VRB.Result (cname, fname, " %d:h_init=%0.14e\n", i, h_init);
	VRB.Result (cname,fname,"Mdag(m_f)|R>=%e\n", lat.FhamiltonNode (tmp3,tmp3));

	// tmp2 <- (M_b^\dag M_b)^{-1} M_f^\dag (RGV)
	tmp2->VecZero (f_size);
{
#ifdef PHI_SAVE
std::stringstream phi_file;
phi_file <<"Phi_traj"<<traj_num<<"."<<phi_veloc[i]<<"."<<UniqueID();
std::string phi_str=phi_file.str();
FILE *fp = fopen(phi_str.c_str(),"r");
printf ("%s::%s %s fp=%p\n", cname,fname,phi_file.str().c_str(),fp);
//exit(4);
if (!fp){
//fclose(fp);
#endif
//        if(veloc_v <1)
	cg_iter = lat.FmatEvlInv (tmp2, tmp3, &bsn_cg_arg[i], CNV_FRM_NO);
#ifdef PHI_SAVE
fp = fopen(phi_str.c_str(),"w");
fwrite (tmp2,sizeof(Float),f_size,fp); 
printf ("%s::%s %s fwrite done\n", cname,fname,phi_str.c_str());
fclose(fp);
} else {
fread (tmp2,sizeof(Float),f_size,fp); 
printf ("%s::%s %s fread done\n", cname,fname,phi_str.c_str());
fclose(fp);
}
#endif
}
	VRB.Result (cname,fname,"(1/MdagM(m_b)) Mdag(m_f)|R>=%e\n", lat.FhamiltonNode (tmp2,tmp2) );
//exit(-4);

	//~~ changed for twisted mass Wilson fermions
	// phi <- M_b (M_b^\dag M_b)^{-1} M_f^\dag (RGV)
	if (lat.Fclass () == F_CLASS_WILSON_TM)
	  lat.SetPhi (phi[i], tmp2, tmp1, bsn_mass[i], bsn_mass_epsilon[i], DAG_NO);
	else
	  lat.SetPhi (phi[i], tmp2, tmp1, bsn_mass[i], DAG_NO);
      }
#ifndef USE_VELOC
	VRB.Result (cname,fname,"%p %p: phi[%d] = M(m_b) (1/MdagM(m_b)) Mdag(m_f)|R>=%e\n", phi[i],phi[i], i, lat.FhamiltonNode (phi[i],phi[i]) );
	VRB.Result (cname,fname,"%p %p: <R|M(m_b) (1/MdagM(m_b)) Mdag(m_f)|R>=%e\n",tmp1,phi[i],  lat.FhamiltonNode (tmp1,phi[i]) );
#endif
      Float total_h_i = h_i;
      glb_sum (&total_h_i);
      VRB.Result (cname, fname,
		  "heatbath: mass ratio %0.4f/(%0.4f) initial ham = %0.16e\n",
		  frm_mass[i], bsn_mass[i], total_h_i);
      VRB.Result (cname, fname,
		  "heatbath: mass ratio %0.4f/(%0.4f) cg_iter = %d\n",
		  frm_mass[i], bsn_mass[i], cg_iter);
      updateCgStats (&bsn_cg_arg[i]);
    }
Float f_temp=0;
glb_sum(&f_temp);
//exit(-42);

    LatticeFactory::Destroy ();
    evolved = 0;
  }

//  dtime += dclock ();
//  print_flops (cname, fname, 0, dtime);
  timer.stop (true);
  VRB.Func (cname, fname);
}

//!< Calculate fermion contribution to the Hamiltonian
Float AlgActionQuotient::energy ()
{

  char *fname = "energy()";
  VRB.Func (cname, fname);
  static Timer timer (cname, fname);
//  Float dtime = -dclock ();
  Float h = 0.0;
  VRB.Result (cname, fname, "evolved h_init= %d %0.14e\n", evolved, h_init);

  if (n_masses > 0) {
//    if (!evolved && h_init != 0.) {
    if (!evolved) {
      return h_init;
    } else {
      timer.start ();
      Lattice & lat = LatticeFactory::Create (fermion, G_CLASS_NONE);

      for (int i = 0; i < n_masses; i++) {
	Float h_i;

	{
	  // everything besides Fdwf4d:

	  //~~ changed for twisted mass Wilson fermions
	VRB.Result (cname,fname,"%p %p : phi[%d]: M(m_b) (1/MdagM(m_b)) Mdag(m_f)|R>=%e\n", phi[i],phi[i],i,lat.FhamiltonNode (phi[i],phi[i]) );
	VRB.Result (cname,fname,"%p %p : <R|M(m_b) (1/MdagM(m_b)) Mdag(m_f)|R>=%e\n", tmp1, phi[i],lat.FhamiltonNode (tmp1,phi[i]) );
	  tmp2->VecZero (f_size);
	  (lat.Fclass () == F_CLASS_WILSON_TM) ?
	    lat.SetPhi (tmp1, phi[i], tmp2, bsn_mass[i], bsn_mass_epsilon[i], DAG_YES) :
	    lat.SetPhi (tmp1, phi[i], tmp2, bsn_mass[i], DAG_YES);
	VRB.Result (cname,fname,"%p %p: : Mdag(m_b)M(m_b) (1/MdagM(m_b)) Mdag(m_f)|R>=%e\n", tmp1,tmp1,lat.FhamiltonNode (tmp1,tmp1));
	VRB.Result (cname,fname,"%p %p: Mdag(m_b)M(m_b) (1/MdagM(m_b)) Mdag(m_f)|R>=%e\n", tmp3,tmp1,lat.FhamiltonNode (tmp3,tmp1));

	  cg_iter = lat.FmatEvlInv (tmp2, tmp1, &frm_cg_arg_mc[i], CNV_FRM_NO);
	VRB.Result (cname,fname,"(1/MdagM(m_f))Mdag(m_b)M(m_b) (1/MdagM(m_b)) Mdag(m_f)|R>=%e\n", lat.FhamiltonNode (tmp2,tmp2));

	  h_i = lat.FhamiltonNode (tmp1, tmp2);
	VRB.Result (cname,fname,"<R|M(m_f)(1/MdagM(m_f))Mdag(m_b)M(m_b) (1/MdagM(m_b)) Mdag(m_f)|R>=%e\n", lat.FhamiltonNode (tmp1,tmp2));
	  h += h_i;
	}
	Float total_h_i = h_i;
	glb_sum (&total_h_i);
//        h_init = h_i - h_init;
//	glb_sum (&h_init);
	VRB.Result (cname, fname,
		    "energy: mass ratio (%0.4f)/%0.4f final ham = %0.16e\n",
		    frm_mass[i], bsn_mass[i], total_h_i);
	VRB.Result (cname, fname,
		    "energy: mass ratio (%0.4f)/%0.4f cg_iter = %d\n",
		    frm_mass[i], bsn_mass[i], cg_iter);
	updateCgStats (&frm_cg_arg_mc[i]);

//      h += lat.FhamiltonNode (tmp1, tmp2);
      }

      LatticeFactory::Destroy ();
    }
  }
//  dtime += dclock ();
//  print_flops (cname, fname, 0, dtime);

  Float total_dh = h - h_init;
  glb_sum (&total_dh);
  VRB.Result (cname, fname,
	      "total delta_ham from this set of quotients = %0.16e\n", total_dh);

  timer.stop (true);
  return h;
}

void AlgActionQuotient::prepare_fg (Matrix * force, Float dt_ratio)
{
  char *fname = "prepare_fg(M*,F)";
  VRB.Func (cname, fname);
  static Timer timer (cname, fname);
  static Timer timer_cg (fname,"prepare_fg::cg()");
  static Timer timer_force (fname,"prepare_fg::force()");
  timer.start ();
//  Float dtime = -dclock ();
//  Float dtime_cg = 0.;
//  Float dtime_force = 0.;

  if (skip_force) {
    VRB.Result (cname, fname,
		"WARNING! skipping prepare_fg() because AlgActionQuotient::skip_force is true!\n");
    evolved = 1;
    timer.stop ();
    return;
  }

  Lattice & lat = LatticeFactory::Create (fermion, G_CLASS_NONE);

  int chronoDeg;
  ForceArg Fdt;
  for (int i = 0; i < n_masses; i++) {

#ifdef USE_BFM
    Fdwf4d::pauli_villars_resid = Fdwf4d::pauli_villars_resid_md;
#endif

    //~~ changed for twisted mass Wilson fermions
    // tmp1 <- (M_b^\dag M_b) (M_b^\dag M_b)^{-1} M_f^\dag (RGV) = M_f^\dag (RGV)
    (lat.Fclass () == F_CLASS_WILSON_TM || lat.Fclass () == F_CLASS_BFM) ?
      lat.SetPhi (tmp1, phi[i], tmp1, bsn_mass[i], bsn_mass_epsilon[i], DAG_YES) : 
      lat.SetPhi (tmp1, phi[i], tmp1, bsn_mass[i], DAG_YES);
     VRB.Result (cname,fname,"%p %p: phi[%d] = M(m_b) (1/MdagM(m_b)) Mdag(m_f)|R>=%e\n", phi[i],phi[i], i, lat.FhamiltonNode (phi[i],phi[i]) );
     moveFloat((Float*)tmp3,(Float*)phi[i],lat.half_size);

    chronoDeg = (md_steps > chrono[i]) ? chrono[i] : md_steps;

    //!< Perform pointer arithmetic to avoid unnecessary copying
    int isz = (chrono[i] > 0) ? (md_steps % chrono[i]) : 0;
    cg_sol = v[i][isz];

    for (int j = 0; j < chrono[i]; j++) {
      int shift = isz - (j + 1);
      if (shift < 0)
	shift += chrono[i];
      cg_sol_old[i][j] = v[i][shift];
    }

    //!< Construct the initial guess
    lat.FminResExt (cg_sol, tmp1, cg_sol_old[i], vm[i],
		    chronoDeg, &frm_cg_arg_fg[i], CNV_FRM_NO);

//    dtime_cg -= dclock ();
    // cg_sol = (M_f^\dag M_f)^{-1} M_f^\dag (RGV)
    timer_cg.start ();
    cg_iter = lat.FmatEvlInv (cg_sol, tmp1, &frm_cg_arg_fg[i], CNV_FRM_NO);
    VRB.Result (cname, fname,
		"prepare_fg: mass ratio (%0.4f)/%0.4f cg_iter = %d\n",
		frm_mass[i], bsn_mass[i], cg_iter);
//    dtime_cg += dclock ();
    timer_cg.stop ();

    updateCgStats (&frm_cg_arg_fg[i]);

    //int g_size = GJP.VolNodeSites() * lat.GsiteSize();

    Matrix *mom_tmp = force;

    if (force_measure == FORCE_MEASURE_YES) {
      mom_tmp =
	(Matrix *) smalloc (g_size * sizeof (Float), "mom_tmp", fname, cname);
      ((Vector *) mom_tmp)->VecZero (g_size);
    }

//    dtime_force -= dclock ();
    timer_force.start ();

    //!< Evolve mom using fermion force
    //~~ changed for twisted mass Wilson fermions
    // cg_sol is aka \chi
    Fdt = (lat.Fclass () == F_CLASS_WILSON_TM || lat.Fclass () == F_CLASS_BFM) ?
      lat.EvolveMomFforce (mom_tmp, cg_sol, frm_mass[i], frm_mass_epsilon[i], dt_ratio) : 
      lat.EvolveMomFforce (mom_tmp, cg_sol, frm_mass[i], dt_ratio);

    if (force_measure == FORCE_MEASURE_YES) {
      char label[200];
      sprintf (label, "%s (fermion), mass = %e:", force_label, frm_mass[i]);
      Fdt.print (dt_ratio, label);
    }
    //!< Evolve mom using boson force
    //~~ changed for twisted mass Wilson fermions
    // cg_sol is aka \chi
    // phi <- M_b (M_b^\dag M_b)^{-1} M_f^\dag (RGV)
    //Twisted mass guy has argument order backwards (nevertheless it is correct)
    if (lat.Fclass () == F_CLASS_WILSON_TM || lat.Fclass () == F_CLASS_BFM)
      Fdt =
	lat.EvolveMomFforce (mom_tmp, cg_sol, phi[i], bsn_mass[i], bsn_mass_epsilon[i], dt_ratio);
    else
      Fdt =
	lat.EvolveMomFforce (mom_tmp, phi[i], cg_sol, bsn_mass[i], dt_ratio);
    moveFloat((Float*)phi[i],(Float*)tmp3,lat.half_size);
	VRB.Result (cname,fname,"%p %p: phi[%d] = M(m_b) (1/MdagM(m_b)) Mdag(m_f)|R>=%e\n", phi[i],phi[i], i, lat.FhamiltonNode (phi[i],phi[i]) );

    //CK: for DWF above does:
    //v1 = (cg_sol, Deo cg_sol)
    //v2 = (-kappa^2 phi[i], -kappa^2 Deo^dag phi[i])
    //where  cg_sol = (M_f^\dag M_f)^{-1} M_f^\dag (RGV)

    //for WilsonTM
    //v1 = f_field_out = (phi[i], g5theta(ctheta,-stheta)D_eo phi[i])
    //v2 = f_field_in  = ( -kappa^2 g5theta(ctheta,stheta)cg_sol, -kappa^2 g5theta(ctheta,stheta) Deo^dag g5theta(ctheta,stheta)cg_sol )


//    dtime_force += dclock ();

    if (force_measure == FORCE_MEASURE_YES) {
      char label[200];
      sprintf (label, "%s (boson), mass = %e:", force_label, bsn_mass[i]);
      Fdt.print (dt_ratio, label);

      // If measuring the force, need to measure and then sum to mom
      Fdt.measure (mom_tmp);
      Fdt.glb_reduce ();

      ((Vector *) force)->VecAddEquVec ((Vector *) mom_tmp, g_size);

      sprintf (label, "%s (total), mass = (%e,%e):", force_label, frm_mass[i],
	       bsn_mass[i]);
      Fdt.print (dt_ratio, label);

      sfree (mom_tmp, "mom_tmp", fname, cname);
    }
    timer_force.stop ();
  }
  // We now have a solution to forecast the next normal solve.
  fg_forecast = true;

  md_steps++;
  LatticeFactory::Destroy ();

//  dtime += dclock ();
//  print_flops (cname, fname, 0, dtime);
//  print_flops (cname, "prepare_fg::cg()", 0, dtime_cg);
//  print_flops (cname, "prepare_fg::force()", 0, dtime_force);
  timer.stop (true);
}

//!< run method evolves the momentum due to the fermion force
void AlgActionQuotient::evolve (Float dt, int nsteps)
{
  char *fname = "evolve(Float,int)";
  VRB.Func (cname, fname);
  static Timer timer (cname, fname);
  static Timer timer_cg (fname,"evolve::cg()");
  static Timer timer_force (fname,"evolve::force()");
  timer.start ();
//  Float dtime = -dclock ();
//  Float dtime_cg = 0.;
//  Float dtime_force = 0.;

  if (skip_force) {
    VRB.Result (cname, fname,
		"WARNING! skipping evolve() because AlgActionQuotient::skip_force is true!\n");
    evolved = 1;
    timer.stop (true);
    return;
  }

  int chronoDeg;
  ForceArg Fdt;
  if (n_masses <= 0)
    return;
  Lattice & lat = LatticeFactory::Create (fermion, G_CLASS_NONE);

  for (int steps = 0; steps < nsteps; steps++) {
    for (int i = 0; i < n_masses; i++) {
      //phi[i] = M_b (M_b^\dag M_b)^{-1} M_f^\dag (RGV)

      // tmp1 <- M_b^dag phi[i] = (M_b^\dag M_b) (M_b^\dag M_b)^{-1} M_f^\dag (RGV) = M_f^\dag (RGV)
      if (lat.Fclass () == F_CLASS_WILSON_TM || lat.Fclass () == F_CLASS_BFM)
	lat.SetPhi (tmp1, phi[i], tmp3, bsn_mass[i], bsn_mass_epsilon[i],
		    DAG_YES);
      else
	lat.SetPhi (tmp1, phi[i], tmp3, bsn_mass[i], DAG_YES);
	VRB.Result (cname,fname,"%p %p: phi[%d] = M(m_b) (1/MdagM(m_b)) Mdag(m_f)|R>=%e\n", phi[i],phi[i], i, lat.FhamiltonNode (phi[i],phi[i]) );
// copying phi[i]

      chronoDeg = (md_steps > chrono[i]) ? chrono[i] : md_steps;

      //!< Perform pointer arithmetic to avoid unnecessary copying
      int isz = (chrono[i] > 0) ? (md_steps % chrono[i]) : 0;
      cg_sol = v[i][isz];

      for (int j = 0; j < chrono[i]; j++) {
	int shift = isz - (j + 1);
	if (shift < 0)
	  shift += chrono[i];
	cg_sol_old[i][j] = v[i][shift];
      }

      //!< Construct the initial guess
      //
      // Branch condition added by Hantao: only if we are NOT using
      // chronological inverter as well as we have a previous
      // fg-solution in hand, we skip the function FminResExt().
      //
      // If chronoDeg == 0 and we don't have a fg forecast, we still
      // want this function to zero the initial guess for us.
      if (fg_forecast == false || chronoDeg != 0) {
	lat.FminResExt (cg_sol, tmp1, cg_sol_old[i], vm[i], chronoDeg, &frm_cg_arg_md[i], CNV_FRM_NO);	//find guess for cg_sol
      } else {
	VRB.Result (cname, fname, "Using force gradient forecasting.\n");
      }

//      dtime_cg -= dclock ();
      // cg_sol = (M_f^\dag M_f)^{-1} tmp1 = (M_f^\dag M_f)^{-1} M_f^\dag (RGV)
      timer_cg.start ();
      cg_iter = lat.FmatEvlInv (cg_sol, tmp1, &frm_cg_arg_md[i], CNV_FRM_NO);
      VRB.Result (cname, fname,
		  "evolve: mass ratio (%0.4f)/%0.4f cg_iter = %d\n",
		  frm_mass[i], bsn_mass[i], cg_iter);
//      dtime_cg += dclock ();
      timer_cg.stop ();

      updateCgStats (&frm_cg_arg_md[i]);

      int g_size = GJP.VolNodeSites () * lat.GsiteSize ();
      if (GJP.Gparity ())
	g_size *= 2;
      Matrix *mom_tmp;

      if (force_measure == FORCE_MEASURE_YES) {
	mom_tmp =
	  (Matrix *) smalloc (g_size * sizeof (Float), "mom_tmp", fname, cname);
	((Vector *) mom_tmp)->VecZero (g_size);
      } else {
	mom_tmp = mom;
      }

//      dtime_force -= dclock ();
      timer_force.start ();
      //!< Evolve mom using fermion force
      if (lat.Fclass () == F_CLASS_WILSON_TM || lat.Fclass () == F_CLASS_BFM)
	Fdt = lat.EvolveMomFforce (mom_tmp, cg_sol, frm_mass[i], frm_mass_epsilon[i], dt);
      else
	Fdt = lat.EvolveMomFforce (mom_tmp, cg_sol, frm_mass[i], dt);

      if (force_measure == FORCE_MEASURE_YES) {
	char label[200];
	sprintf (label, "%s (fermion), mass = %e:", force_label, frm_mass[i]);
	Fdt.print (dt, label);
      }
      //!< Evolve mom using boson force
      //Arguments for general case:   (M_f^\dag M_f)^{-1} M_f^\dag (RGV)   and   M_b (M_b^\dag M_b)^{-1} M_f^\dag (RGV)
      //Arguments for WilsonTM case:  M_b (M_b^\dag M_b)^{-1} M_f^\dag (RGV)   and    (M_f^\dag M_f)^{-1} M_f^\dag (RGV)

      moveFloat((Float*)tmp3,(Float*)phi[i],lat.half_size);
      if (lat.Fclass () == F_CLASS_WILSON_TM || (lat.Fclass () == F_CLASS_BFM))
	Fdt = lat.EvolveMomFforce (mom_tmp, cg_sol, phi[i], bsn_mass[i], bsn_mass_epsilon[i], dt);
      else
	Fdt = lat.EvolveMomFforce (mom_tmp, phi[i], cg_sol, bsn_mass[i], dt);
      moveFloat((Float*)phi[i],(Float*)tmp3,lat.half_size);
	VRB.Result (cname,fname,"%p %p: phi[%d] = M(m_b) (1/MdagM(m_b)) Mdag(m_f)|R>=%e\n", phi[i],phi[i], i, lat.FhamiltonNode (phi[i],phi[i]) );

//      dtime_force += dclock ();

      if (force_measure == FORCE_MEASURE_YES) {
	char label[200];
	sprintf (label, "%s (boson), mass = %e:", force_label, bsn_mass[i]);
	Fdt.print (dt, label);

	// If measuring the force, need to measure and then sum to mom
	Fdt.measure (mom_tmp);
	Fdt.glb_reduce ();

	fTimesV1PlusV2 ((IFloat *) mom, 1.0, (IFloat *) mom_tmp, (IFloat *) mom,
			g_size);
	sprintf (label, "%s (total), mass = (%e,%e):", force_label, frm_mass[i],
		 bsn_mass[i]);
	Fdt.print (dt, label);

	sfree (mom_tmp, "mom_tmp", fname, cname);
      }
      timer_force.stop ();
    }
    // Note that as long as the last solve in a trajectory is NOT a
    // force gradient solve (which should always be the case), we
    // won't bring our solution across trajectories by using the
    // statement here.
    fg_forecast = false;

    md_steps++;
  }

  LatticeFactory::Destroy ();
  evolved = 1;

//  dtime += dclock ();
//  print_flops (cname, fname, 0, dtime);
//  print_flops (cname, "evolve::cg()", 0, dtime_cg);
//  print_flops (cname, "evolve::force()", 0, dtime_force);
  timer.stop (true);
}

CPS_END_NAMESPACE
