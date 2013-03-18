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
CPS_START_NAMESPACE

AlgActionQuotient::AlgActionQuotient(AlgMomentum &mom,
                                     ActionQuotientArg &q_arg)
    : AlgActionBilinear(mom, q_arg.bi_arg),
      cname("AlgActionQuotient")
{
    const char *fname = "AlgActionQuotient()";

    int_type = INT_QUOTIENT;
    quo_arg = &q_arg;

    //!< First check n_masses bilinear = n_masses quotient
    if (quo_arg->quotients.quotients_len != quo_arg->bi_arg.bilinears.bilinears_len)
        ERR.General(cname, fname,
                    "Inconsistency between QuotientArg and BilinearArg n_masses\n");

    if(n_masses > 0) {
        bsn_cg_arg.resize(n_masses);
        frm_cg_arg_fg.resize(n_masses);
        frm_cg_arg_md.resize(n_masses);
        frm_cg_arg_mc.resize(n_masses);

        //!< Initialize the CG arguments
        for(int i=0; i<n_masses; i++) {
            const QuotientDescr &qi = quo_arg->quotients.quotients_val[i];

            bsn_mass.push_back(qi.bsn_mass);
            frm_mass.push_back(qi.frm_mass);

            //~~ added for twisted mass Wilson fermions
            bsn_mass_epsilon.push_back(qi.bsn_mass_epsilon);
            frm_mass_epsilon.push_back(qi.frm_mass_epsilon);

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

            chrono.push_back(qi.chrono);
        }

        evolved = 1;

        //!< Vectors used to store solution history
        v = (Vector***) smalloc(n_masses*sizeof(Vector**),
                                "v", fname, cname);
        cg_sol_old = (Vector***) smalloc(n_masses*sizeof(Vector**),
                                         "cg_sol_old", fname, cname);
        vm = (Vector***) smalloc(n_masses*sizeof(Vector**),
                                 "vm", fname, cname);

        tmp1 = (Vector*)smalloc(f_size*sizeof(Float),"tmp1",fname,cname);
        tmp2 = (Vector*)smalloc(f_size*sizeof(Float),"tmp2",fname,cname);
    
        for (int i=0; i<n_masses; i++) {
            int deg=0;
            if (chrono[i] > 0) deg = chrono[i];
            else if (chrono[i] == 0) deg = 1;
            else ERR.General(cname,fname,"Cannot have negative chronology\n");

            v[i] = (Vector**) smalloc(deg*sizeof(Vector*),"v[i]", fname, cname);
            cg_sol_old[i] = (Vector**) smalloc(deg*sizeof(Vector*),
                                               "cg_sol_old[i]", fname, cname);
            vm[i] = (Vector**) smalloc(deg*sizeof(Vector*), "vm[i]", fname, cname);
            for (int j=0; j<deg; j++) {
                v[i][j] = (Vector*) smalloc(f_size*sizeof(Float),
                                            "v[i][j]", fname, cname);
                vm[i][j] = (Vector*) smalloc(f_size*sizeof(Float),
                                             "vm[i][j]", fname, cname);
            }
        }

    }

    init();

    fg_forecast = false;
}

void AlgActionQuotient::init()
{
    AlgActionBilinear::init();
    evolved = 1;
    for (int i=0; i<n_masses; i++) v[i][0]->VecZero(f_size);
}

AlgActionQuotient::~AlgActionQuotient()
{
    const char *fname = "~AlgActionQuotient()";

    if(n_masses > 0) {
        //!< Free chronology
        for (int i=0; i<n_masses; i++) {
            int deg=0;
            if (chrono[i] > 0) deg = chrono[i];
            else if (chrono[i] == 0) deg = 1;

            for (int j=0; j<deg; j++) {
                sfree(vm[i][j],"vm[i][j]", fname, cname);
                sfree(v[i][j],"v[i][j]", fname, cname);
            }
            sfree(cg_sol_old[i], "cg_sol_old[i]", fname, cname);
            sfree(vm[i],"vm[i]", fname, cname);
            sfree(v[i],"v[i]", fname, cname);
        }
        sfree(tmp2, "tmp2", fname, cname);
        sfree(tmp1, "tmp1", fname, cname);
        sfree(cg_sol_old, "cg_sol_old", fname, cname);
        sfree(vm, "vm", fname, cname);
        sfree(v, "v", fname, cname);
    }

}

//!< Heat Bath for quotients
void AlgActionQuotient::reweight(Float *rw_fac, Float *norm) {

    char *fname = "reweight()";
  
    if (n_masses > 0) {
        Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);
    
        //    h_init = 0.0;

        // tmp1, tmp2 < - random Gaussian vector (RGV)
        for(int i=0; i<n_masses; i++){
            lat.RandGaussVector(tmp1, 0.5, Ncb);
            lat.RandGaussVector(tmp2, 0.5, Ncb);

            //~~ changed for twisted mass Wilson fermions
            // phi <- M_f^\dag (RGV)
            norm[i] = (lat.Fclass() == F_CLASS_WILSON_TM) ?
                lat.SetPhi(phi[i], tmp1, tmp2, bsn_mass[i], bsn_mass_epsilon[i], DAG_YES) :
                lat.SetPhi(phi[i], tmp1, tmp2, bsn_mass[i], DAG_YES);
	 
            // tmp2 <- (M_b^\dag M_b)^{-1} M_f^\dag (RGV)
            tmp2 -> VecZero(f_size);
            cg_iter = lat.FmatEvlInv(tmp2, phi[i], &frm_cg_arg_mc[i], CNV_FRM_NO);

            rw_fac[i] = lat.FhamiltonNode(phi[i], tmp2);
            rw_fac[i] -= norm[i];
            VRB.Result(cname,fname,"rw_fac=%e norm=%e\n",rw_fac[i],norm[i]);
            glb_sum(rw_fac+i);
            glb_sum(norm+i);
      
            updateCgStats(&bsn_cg_arg[i]);          
        }

        LatticeFactory::Destroy();
        //    evolved = 0;
    }

}

//!< Heat Bath for quotients
void AlgActionQuotient::heatbath() {

    char *fname = "heatbath()";
    Float dtime = -dclock();
  
    if (n_masses > 0) {
        Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);
    
        h_init = 0.0;

        // tmp1, tmp2 < - random Gaussian vector (RGV)
        for(int i=0; i<n_masses; i++){
            lat.RandGaussVector(tmp1, 0.5, Ncb);
            lat.RandGaussVector(tmp2, 0.5, Ncb);

            //~~ changed for twisted mass Wilson fermions
            // phi <- M_f^\dag (RGV)
            h_init += (lat.Fclass() == F_CLASS_WILSON_TM) ?
                lat.SetPhi(phi[i], tmp1, tmp2, frm_mass[i], frm_mass_epsilon[i], DAG_YES) :
                lat.SetPhi(phi[i], tmp1, tmp2, frm_mass[i], DAG_YES);
	 
            // tmp2 <- (M_b^\dag M_b)^{-1} M_f^\dag (RGV)
            tmp2 -> VecZero(f_size);
            cg_iter = lat.FmatEvlInv(tmp2, phi[i], &bsn_cg_arg[i], CNV_FRM_NO);
      
            //~~ changed for twisted mass Wilson fermions
            // phi <- M_b (M_b^\dag M_b)^{-1} M_f^\dag (RGV)
            (lat.Fclass() == F_CLASS_WILSON_TM) ?
                lat.SetPhi(phi[i], tmp2, tmp1, bsn_mass[i], bsn_mass_epsilon[i], DAG_NO) :
                lat.SetPhi(phi[i], tmp2, tmp1, bsn_mass[i], DAG_NO);
            updateCgStats(&bsn_cg_arg[i]);

            //      h_init2 += lat.FhamiltonNode(phi[i], phi[i]);
        }
        //    h_init2 = h_init - h_init2;
        //    glb_sum(&h_init2);

        //  VRB.Result(cname, fname, "Hamiltonian = %17.12e\n", h_init2);

        LatticeFactory::Destroy();
        evolved = 0;
    }

    dtime += dclock();
    print_flops(cname, fname, 0, dtime);
}

//!< Calculate fermion contribution to the Hamiltonian
Float AlgActionQuotient::energy() {

    char *fname = "energy()";
    Float dtime = -dclock();
    Float h = 0.0;

    if (n_masses > 0) {
        if (!evolved && h_init != 0.0) {
            return h_init;
        } else {
            Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);

            for(int i=0; i<n_masses; i++) {
   
                //~~ changed for twisted mass Wilson fermions
		(lat.Fclass() == F_CLASS_WILSON_TM) ?
                    lat.SetPhi(tmp1, phi[i], tmp2, bsn_mass[i], bsn_mass_epsilon[i], DAG_YES) :
                    lat.SetPhi(tmp1, phi[i], tmp2, bsn_mass[i], DAG_YES);

                tmp2 -> VecZero(f_size);
                cg_iter = 
                    lat.FmatEvlInv(tmp2, tmp1, &frm_cg_arg_mc[i], CNV_FRM_NO);
	
                updateCgStats(&frm_cg_arg_mc[i]);
	  
                h += lat.FhamiltonNode(tmp1, tmp2);
            }
      
            LatticeFactory::Destroy();
        }
    } 
    dtime += dclock();
    print_flops(cname, fname, 0, dtime);

    return h;
}

void AlgActionQuotient::prepare_fg(Matrix * force, Float dt_ratio)
{
    char * fname = "prepare_fg(M*,F)";
    Float dtime = -dclock();
    Float dtime_cg = 0.;
    Float dtime_force = 0.;

    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);

    int chronoDeg;
    ForceArg Fdt;
    for(int i=0; i<n_masses; i++) {
        //~~ changed for twisted mass Wilson fermions
        // tmp1 <- (M_b^\dag M_b) (M_b^\dag M_b)^{-1} M_f^\dag (RGV) = M_f^\dag (RGV)
        (lat.Fclass() == F_CLASS_WILSON_TM) ?
            lat.SetPhi(tmp1, phi[i], tmp1, bsn_mass[i], bsn_mass_epsilon[i], DAG_YES) :
            lat.SetPhi(tmp1, phi[i], tmp1, bsn_mass[i], DAG_YES);

        chronoDeg = ( md_steps > chrono[i] ) ? chrono[i] : md_steps ;

        //!< Perform pointer arithmetic to avoid unnecessary copying
        int isz = ( chrono[i] > 0 ) ? ( md_steps % chrono[i] ) : 0 ;
        cg_sol = v[i][isz];
	
        for (int j=0; j<chrono[i]; j++) {
            int shift = isz - (j+1);
            if (shift<0) shift += chrono[i];
            cg_sol_old[i][j] = v[i][shift];
        }

        //!< Construct the initial guess
        lat.FminResExt(cg_sol, tmp1, cg_sol_old[i], vm[i], 
                       chronoDeg, &frm_cg_arg_fg[i], CNV_FRM_NO);

        dtime_cg -= dclock();
        // cg_sol = (M_f^\dag M_f)^{-1} M_f^\dag (RGV)
        cg_iter = 
            lat.FmatEvlInv(cg_sol, tmp1, &frm_cg_arg_fg[i], CNV_FRM_NO);
        dtime_cg += dclock();

        updateCgStats(&frm_cg_arg_fg[i]);

        //int g_size = GJP.VolNodeSites() * lat.GsiteSize();

        Matrix * mom_tmp = force;

        if (force_measure == FORCE_MEASURE_YES) {
            mom_tmp = (Matrix*)smalloc(g_size*sizeof(Float),"mom_tmp", fname, cname);
            ((Vector*)mom_tmp) -> VecZero(g_size);
        }
    
        dtime_force -= dclock();

        //!< Evolve mom using fermion force
        //~~ changed for twisted mass Wilson fermions
        // cg_sol is aka \chi
        Fdt = (lat.Fclass() == F_CLASS_WILSON_TM) ?
            lat.EvolveMomFforce(mom_tmp, cg_sol, frm_mass[i], frm_mass_epsilon[i], dt_ratio) :
            lat.EvolveMomFforce(mom_tmp, cg_sol, frm_mass[i], dt_ratio);

        if (force_measure == FORCE_MEASURE_YES) {
            char label[200];
            sprintf(label, "%s (fermion), mass = %e:", force_label, frm_mass[i]);
            Fdt.print(dt_ratio, label);
        }
	
        //!< Evolve mom using boson force
        //~~ changed for twisted mass Wilson fermions
        // cg_sol is aka \chi
        // phi <- M_b (M_b^\dag M_b)^{-1} M_f^\dag (RGV)
        Fdt = (lat.Fclass() == F_CLASS_WILSON_TM) ?
            lat.EvolveMomFforce(mom_tmp, cg_sol, phi[i], bsn_mass[i], bsn_mass_epsilon[i], dt_ratio) :
            // TB flipped cg_sol <---> phi since this is what was in original
            lat.EvolveMomFforce(mom_tmp, phi[i], cg_sol, bsn_mass[i], dt_ratio) ;
	
        dtime_force += dclock();

        if (force_measure == FORCE_MEASURE_YES) {
            char label[200];
            sprintf(label, "%s (boson), mass = %e:", force_label, bsn_mass[i]);
            Fdt.print(dt_ratio, label);
	  
            // If measuring the force, need to measure and then sum to mom
            Fdt.measure(mom_tmp);
            Fdt.glb_reduce();

            ((Vector *)force)->VecAddEquVec((Vector *)mom_tmp, g_size);

            sprintf(label, "%s (total), mass = (%e,%e):", force_label, frm_mass[i], bsn_mass[i]);
            Fdt.print(dt_ratio, label);

            sfree(mom_tmp, "mom_tmp", fname, cname);
        }
    }
    // We now have a solution to forecast the next normal solve.
    fg_forecast = true;
  
    md_steps++;
    LatticeFactory::Destroy();

    dtime += dclock();
    print_flops(cname, fname, 0, dtime);
    print_flops(cname, "prepare_fg::cg()", 0, dtime_cg);
    print_flops(cname, "prepare_fg::force()", 0, dtime_force);
}

//!< run method evolves the momentum due to the fermion force
void AlgActionQuotient::evolve(Float dt, int nsteps) 
{
    char *fname = "evolve(Float, int)";
    Float dtime = -dclock();
    Float dtime_cg = 0.;
    Float dtime_force = 0.;

    int chronoDeg;
    ForceArg Fdt;
    if (n_masses <= 0) return;
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);

    for (int steps=0; steps<nsteps; steps++) {
        for(int i=0; i<n_masses; i++) {

            //~~ changed for twisted mass Wilson fermions
            // tmp1 <- (M_b^\dag M_b) (M_b^\dag M_b)^{-1} M_f^\dag (RGV) = M_f^\dag (RGV)
            (lat.Fclass() == F_CLASS_WILSON_TM) ?
                lat.SetPhi(tmp1, phi[i], tmp1, bsn_mass[i], bsn_mass_epsilon[i], DAG_YES) :
                lat.SetPhi(tmp1, phi[i], tmp1, bsn_mass[i], DAG_YES);

            chronoDeg = ( md_steps > chrono[i] ) ? chrono[i] : md_steps ;

            //!< Perform pointer arithmetic to avoid unnecessary copying
            int isz = ( chrono[i] > 0 ) ? ( md_steps % chrono[i] ) : 0 ;
            cg_sol = v[i][isz];
	
            for (int j=0; j<chrono[i]; j++) {
                int shift = isz - (j+1);
                if (shift<0) shift += chrono[i];
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
                lat.FminResExt(cg_sol, tmp1, cg_sol_old[i], vm[i], 
                               chronoDeg, &frm_cg_arg_md[i], CNV_FRM_NO);
            }else{
                VRB.Result(cname, fname, "Using force gradient forecasting.\n");
            }

            dtime_cg -= dclock();
            // cg_sol = (M_f^\dag M_f)^{-1} M_f^\dag (RGV)
            cg_iter = lat.FmatEvlInv(cg_sol, tmp1, &frm_cg_arg_md[i], CNV_FRM_NO);
            dtime_cg += dclock();

            updateCgStats(&frm_cg_arg_md[i]);

            int g_size = GJP.VolNodeSites() * lat.GsiteSize();
            Matrix *mom_tmp;

            if (force_measure == FORCE_MEASURE_YES) {
                mom_tmp = (Matrix*)smalloc(g_size*sizeof(Float),"mom_tmp", fname, cname);
                ((Vector*)mom_tmp) -> VecZero(g_size);
            } else {
                mom_tmp = mom;
            }

            dtime_force -= dclock();
            //!< Evolve mom using fermion force
            //~~ changed for twisted mass Wilson fermions
            // cg_sol is aka \chi
            Fdt = (lat.Fclass() == F_CLASS_WILSON_TM) ?
                lat.EvolveMomFforce(mom_tmp, cg_sol, frm_mass[i], frm_mass_epsilon[i], dt) :
                lat.EvolveMomFforce(mom_tmp, cg_sol, frm_mass[i], dt);
            if (force_measure == FORCE_MEASURE_YES) {
                char label[200];
                sprintf(label, "%s (fermion), mass = %e:", force_label, frm_mass[i]);
                Fdt.print(dt, label);
            }
	
            //!< Evolve mom using boson force
            //~~ changed for twisted mass Wilson fermions
            // cg_sol is aka \chi
            // phi <- M_b (M_b^\dag M_b)^{-1} M_f^\dag (RGV)
            Fdt = (lat.Fclass() == F_CLASS_WILSON_TM) ?
                lat.EvolveMomFforce(mom_tmp, cg_sol, phi[i], bsn_mass[i], bsn_mass_epsilon[i], dt) :
                // TB flipped cg_sol <---> phi since this is what was in original
                lat.EvolveMomFforce(mom_tmp, phi[i], cg_sol, bsn_mass[i], dt) ;
            dtime_force += dclock();
	
            if (force_measure == FORCE_MEASURE_YES) {
                char label[200];
                sprintf(label, "%s (boson), mass = %e:", force_label, bsn_mass[i]);
                Fdt.print(dt, label);
	  
                // If measuring the force, need to measure and then sum to mom
                Fdt.measure(mom_tmp);
                Fdt.glb_reduce();

                fTimesV1PlusV2((IFloat*)mom, 1.0, (IFloat*)mom_tmp, (IFloat*)mom, g_size);
                sprintf(label, "%s (total), mass = (%e,%e):", force_label, frm_mass[i], bsn_mass[i]); 
                Fdt.print(dt, label);
	  
                sfree(mom_tmp, "mom_tmp", fname, cname);
            }

        }
        // Note that as long as the last solve in a trajectory is NOT a
        // force gradient solve (which should always be the case), we
        // won't bring our solution across trajectories by using the
        // statement here.
        fg_forecast = false;

        md_steps++;
    }

    LatticeFactory::Destroy();
    evolved = 1;

    dtime += dclock();
    print_flops(cname, fname, 0, dtime);
    print_flops(cname, "evolve::cg()", 0, dtime_cg);
    print_flops(cname, "evolve::force()", 0, dtime_force);
}

CPS_END_NAMESPACE
