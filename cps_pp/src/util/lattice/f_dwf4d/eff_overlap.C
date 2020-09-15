#ifdef USE_BFM

#include <util/lattice/eff_overlap.h>
#include <util/lattice/bfm_evo.h>
#include <util/lattice/bfm_eigcg.h>
#include <util/lattice/bfm_mixed_solver.h>
#ifndef BFM_GPARITY
#include <util/lattice/hdcg_controller.h>
#endif
#include <util/timer.h>
#include <util/gjp.h>
#include <util/smalloc.h>

CPS_START_NAMESPACE


static const char* cname = "";

// Applies Dov
//
// How to apply the 4D overlap operator:
//
// 1. Promote the 4D input vector to a 5D vector, putting the left-handed part
//    at s=0 and the right-handed part at s=Ls-1
//
// 3. Apply D_DW(m)    
//
// 4. Apply D_DW(1)^-1          
//
// 5. Reduce the 5D vector to a 4D vector, taking the left-handed part from s=0
//    and the right-handed part from s=Ls-1
//
// We don't have to worry about issues with D_- here because 
//
// D_DW(1)^-1 D_DW(m) = (D_-^-1 D_DW(1))^-1 (D_-^-1 D_DW(m))
//
// (the D_-'s cancel in the product we compute).
int ApplyOverlap(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Vector *out, Vector *in, Float mass, Float pv_stop_rsd)
{
    const char* fname = "ApplyOverlap()";

    static Timer timer(fname);
    static std::map<Float, Timer*> timers;
    if (timers.count(mass) == 0) {
	char timer_mass_name[512];
	sprintf(timer_mass_name, "ApplyOverlap(mass=%0.4f)", mass);
	timers[mass] = new Timer(timer_mass_name);
    }
    timer.start(true);
    timers[mass]->start(true);

    Fermion_t vec_5d[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };
    Fermion_t tmp_5d[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };

    // vec_5d = 5D version of input vector
    bfm_d.cps_impexFermion_4d((Float *)in, vec_5d, Import);

    bfm_d.set_mass(mass);
#pragma omp parallel
    {
	// tmp_5d = D_DW(m) vec_5d
	bfm_d.G5D_Munprec(vec_5d, tmp_5d, DaggerNo);
    }

    int iters;
    bfm_d.set_mass(1.0);
    if (use_mixed_solver) bfm_f.set_mass(1.0);
    bfm_d.residual = pv_stop_rsd;
    if (use_mixed_solver) bfm_f.residual = 1e-5;
#pragma omp parallel
    {
	// vec_5d = D_DW(1)^-1 tmp_5d = D_DW(1)^-1 D_DW(m) in
	bfm_d.set_zero(vec_5d[Even]);
	bfm_d.set_zero(vec_5d[Odd]);
	iters = use_mixed_solver ?
	    mixed_cg::threaded_cg_mixed_M(vec_5d, tmp_5d, bfm_d, bfm_f, 5) :
	    bfm_d.CGNE_M(vec_5d, tmp_5d);
    }

    // out = 4D version of vec_5d
    bfm_d.cps_impexFermion_4d((Float *)out, vec_5d, Export);

    bfm_d.freeFermion(vec_5d[Even]);
    bfm_d.freeFermion(vec_5d[Odd]);
    bfm_d.freeFermion(tmp_5d[Even]);
    bfm_d.freeFermion(tmp_5d[Odd]);

    VRB.Result("", fname, "Applied overlap operator for mass = %e, Ls = %d: PV iterations = %d\n", mass, bfm_d.Ls, iters);

    timers[mass]->stop(true);
    timer.stop(true);
    
    return iters;
}

// Invert the 4D effective overlap operator. 
// Can use EigCG for the inner 5D solve if itype is set to EIGCG
int ApplyOverlapInverse(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Vector *out, Vector *in, Float mass, Float stop_rsd, InverterType itype)
{
    const char* fname = "ApplyOverlapInverse()";

    static Timer timer(fname);
    static std::map<Float, Timer*> timers;
    if (timers.count(mass) == 0) {
	char timer_mass_name[512];
	sprintf(timer_mass_name, "ApplyOverlapInverse(mass=%0.4f)", mass);
	timers[mass] = new Timer(timer_mass_name);
    }
    timer.start(true);
    timers[mass]->start(true);

    Fermion_t vec_5d[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };
    Fermion_t tmp_5d[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };

    bfm_d.cps_impexFermion_4d((Float *)in, vec_5d, Import);

    // Do
    //   tmp_5d = D_DW(1) vec_5d = D_DW(1) P in
    bfm_d.set_mass(1.0);
#pragma omp parallel
    {
	bfm_d.G5D_Munprec(vec_5d, tmp_5d, DaggerNo);
    }

    // Now do
    //   vec_5d = D_DW(m)^{-1} tmp_5d = D_DW(m)^{-1} D_DW(1) P in
    int iters;
    bfm_d.set_mass(mass);
    if (use_mixed_solver) bfm_f.set_mass(mass);
    bfm_d.residual = stop_rsd;
    if (use_mixed_solver) bfm_f.residual = 1e-5;
#pragma omp parallel
    {
	bfm_d.set_zero(vec_5d[Even]);
	bfm_d.set_zero(vec_5d[Odd]);
    }

#ifndef BFM_GPARITY
    if (itype == HDCG) {
	HDCG_wrapper *control = HDCGInstance::getInstance();
	assert(control != NULL);
	control->HDCG_set_mass(mass);
	control->HDCG_invert(vec_5d, tmp_5d, stop_rsd, bfm_d.max_iter);
    } else 
#endif
{
#pragma omp parallel
	{
	    if (use_mixed_solver) {
		iters = mixed_cg::threaded_cg_mixed_M(vec_5d, tmp_5d, bfm_d, bfm_f, 5, itype);
	    } else {
		switch (itype) {
		    case CG:
			iters = bfm_d.CGNE_M(vec_5d, tmp_5d);
			break;
		    case EIGCG:
			iters = bfm_d.EIG_CGNE_M(vec_5d, tmp_5d);
			break;
		    default:
			if (bfm_d.isBoss()) {
			    printf("%s::%s: Inverter type %d not implemented\n", cname, fname, itype);
			}
			exit(-1);
			break;
		}
	    }
	}
    }

    // out = [P^{-1} vec_5d]_0 = [P^{-1} D_DW(m)^{-1} D_DW(1) P]_00 in
    bfm_d.cps_impexFermion_4d((Float *)out, vec_5d, Export);

    bfm_d.freeFermion(vec_5d[Even]);
    bfm_d.freeFermion(vec_5d[Odd]);
    bfm_d.freeFermion(tmp_5d[Even]);
    bfm_d.freeFermion(tmp_5d[Odd]);

    timers[mass]->stop(true);
    timer.stop(true);

    return iters;
}

int ApplyOverlapDag(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Vector *out, Vector *in, Float mass, Float pv_stop_rsd)
{
    const char* fname = "ApplyOverlapDag()";

    static Timer timer(fname);
    static std::map<Float, Timer*> timers;
    if (timers.count(mass) == 0) {
	char timer_mass_name[512];
	sprintf(timer_mass_name, "ApplyOverlapDag(mass=%0.4f)", mass);
	timers[mass] = new Timer(timer_mass_name);
    }
    timer.start(true);
    timers[mass]->start(true);


    Fermion_t vec_5d[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };
    Fermion_t tmp_5d[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };

    // vec_5d = 5D version of input vector
    bfm_d.cps_impexFermion_4d((Float *)in, vec_5d, Import);

    bfm_d.set_mass(1.0);
    if (use_mixed_solver) bfm_f.set_mass(1.0);
    bfm_d.residual = pv_stop_rsd;
    if (use_mixed_solver) bfm_f.residual = 1e-5;
    int iters;
#ifdef BFM_GPARITY
   ERR.General(cname,fname,"mixed_cg::threaded_cg_mixed_Mdag_guess not implemented for BFM with Gparity\n");
#else
#pragma omp parallel
    {
	// tmp_5d = D_DW(1)^dag^-1 vec_5d
	bfm_d.set_zero(tmp_5d[Even]);
	bfm_d.set_zero(tmp_5d[Odd]);
	iters = use_mixed_solver ?
	    mixed_cg::threaded_cg_mixed_Mdag_guess(tmp_5d, vec_5d, bfm_d, bfm_f, 5) :
	    bfm_d.CGNE_Mdag(tmp_5d, vec_5d);
    }
#endif

    bfm_d.set_mass(mass);
#pragma omp parallel
    {
	// vec_5d = D_DW(m)^-1^dag tmp_5d = D_DW(m)^-1^dag D_DW(1)^dag in
	bfm_d.G5D_Munprec(tmp_5d, vec_5d, DaggerYes);
    }

    // out = 4D version of vec_5d
    bfm_d.cps_impexFermion_4d((Float *)out, vec_5d, Export);

    bfm_d.freeFermion(vec_5d[Even]);
    bfm_d.freeFermion(vec_5d[Odd]);
    bfm_d.freeFermion(tmp_5d[Even]);
    bfm_d.freeFermion(tmp_5d[Odd]);

    VRB.Result("", fname, "Applied overlap-dagger operator for mass = %e, Ls = %d: PV iterations = %d\n", mass, bfm_d.Ls, iters);

    timers[mass]->stop(true);
    timer.stop(true);

    return iters;
}

int ApplyOverlapDagInverse(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Vector *out, Vector *in, Float mass, Float stop_rsd)
{
    const char* fname = "ApplyOverlapDagInverse()";

    static Timer timer(fname);
    static std::map<Float, Timer*> timers;
    if (timers.count(mass) == 0) {
	char timer_mass_name[512];
	sprintf(timer_mass_name, "ApplyOverlapDagInverse(mass=%0.4f)", mass);
	timers[mass] = new Timer(timer_mass_name);
    }
    timer.start(true);
    timers[mass]->start(true);

    Fermion_t vec_5d[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };
    Fermion_t tmp_5d[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };

    bfm_d.cps_impexFermion_4d((Float *)in, vec_5d, Import);

    bfm_d.set_mass(mass);
    if (use_mixed_solver) bfm_f.set_mass(mass);
    bfm_d.residual = stop_rsd;
    if (use_mixed_solver) bfm_f.residual = 1e-5;
    int iters;
#ifdef BFM_GPARITY
   ERR.General(cname,fname,"mixed_cg::threaded_cg_mixed_Mdag_guess not implemented for BFM with Gparity\n");
#else
#pragma omp parallel
    {
	// TODO: optionally make use of initial guess
	bfm_d.set_zero(tmp_5d[Even]);
	bfm_d.set_zero(tmp_5d[Odd]);
	iters = use_mixed_solver ?
	    mixed_cg::threaded_cg_mixed_Mdag_guess(tmp_5d, vec_5d, bfm_d, bfm_f, 5) :
	    bfm_d.CGNE_Mdag(tmp_5d, vec_5d);
    }
#endif

    bfm_d.set_mass(1.0);
#pragma omp parallel
    {
	bfm_d.G5D_Munprec(tmp_5d, vec_5d, DaggerYes);
    }

    bfm_d.cps_impexFermion_4d((Float *)out, vec_5d, Export);

    bfm_d.freeFermion(vec_5d[Even]);
    bfm_d.freeFermion(vec_5d[Odd]);
    bfm_d.freeFermion(tmp_5d[Even]);
    bfm_d.freeFermion(tmp_5d[Odd]);

    timers[mass]->stop(true);
    timer.stop(true);

    return iters;
}

// Like ApplyOverlapInverse, but uses the initial value of out
// as a guess. This requires an extra inversion of D_DW(1)
int ApplyOverlapInverseGuess(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Vector *out, Vector *in, Float mass, Float stop_rsd)
{
    const char* fname = "ApplyOverlapInverseGuess";

    Fermion_t out_5d[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };
    Fermion_t Dm_Pout[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };
    Fermion_t in_5d[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };
    Fermion_t rhs_5d[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };

    // out_5d = P out
    bfm_d.cps_impexFermion_4d((Float *)out, out_5d, Import);

    // in_5d = P in
    bfm_d.cps_impexFermion_4d((Float *)in, in_5d, Import);

    // Contruct initial guess for 5D solve

    // Dm_Pout = -D_DW(m) P out
    bfm_d.set_mass(mass);
#pragma omp parallel
    {
	bfm_d.scale(out_5d, -1.0);
	bfm_d.G5D_Munprec(out_5d, Dm_Pout, DaggerNo);
    }
    
    VRB.Result(cname, fname, "going to do Pauli-Villars solve to get 5D guess\n");
    bfm_d.set_mass(1.0);
    if (use_mixed_solver) bfm_f.set_mass(1.0);
    bfm_d.residual = stop_rsd;
    if (use_mixed_solver) bfm_f.residual = 1e-5;
    int iters_PV;
#pragma omp parallel 
    {
	// out_5d = -D_DW(1)^{-1} D_DW(m) P out
	bfm_d.set_zero(out_5d[Even]);
	bfm_d.set_zero(out_5d[Odd]);
	iters_PV = use_mixed_solver ?
	    mixed_cg::threaded_cg_mixed_M(out_5d, Dm_Pout, bfm_d, bfm_f, 5) :
	    bfm_d.CGNE_M(out_5d, Dm_Pout);

	// Construct RHS for 5D solve
	// rhs_5d = D_DW(1) P in
	bfm_d.G5D_Munprec(in_5d, rhs_5d, DaggerNo);
    }

    // get the zero component of the initial guess.
    // We do this by again importing out into out_5d. out_5d now has all
    // the other components of the initial guess, so we set the "prezero"
    // argument to false and just import to the zero component.
    bfm_d.cps_impexFermion_4d((Float *)out, out_5d, Import, false); 

    // Do the main 5D solve
    VRB.Result(cname, fname, "Now doing main solve\n");
    bfm_d.set_mass(mass);
    if (use_mixed_solver) bfm_f.set_mass(mass);
    bfm_d.residual = stop_rsd;
    if (use_mixed_solver) bfm_f.residual = 1e-5;
    int iters_main;
#pragma omp parallel
    {
	iters_main = use_mixed_solver ?
	    mixed_cg::threaded_cg_mixed_M(out_5d, rhs_5d, bfm_d, bfm_f, 5) :
	    bfm_d.CGNE_M(out_5d, rhs_5d);
    }

    bfm_d.cps_impexFermion_4d((Float *)out, out_5d, Export);

    bfm_d.freeFermion(out_5d[Even]);
    bfm_d.freeFermion(out_5d[Odd]);
    bfm_d.freeFermion(Dm_Pout[Even]);
    bfm_d.freeFermion(Dm_Pout[Odd]);
    bfm_d.freeFermion(in_5d[Even]);
    bfm_d.freeFermion(in_5d[Odd]);
    bfm_d.freeFermion(rhs_5d[Even]);
    bfm_d.freeFermion(rhs_5d[Odd]);

    int iters = iters_PV + iters_main;
    VRB.Result(cname, fname, "iters_PV   = %d\n", iters_PV);
    VRB.Result(cname, fname, "iters_main = %d\n", iters_main);
    VRB.Result(cname, fname, "final total iters = %d at Ls = %d\n", iters, bfm_d.Ls);
    VRB.Result(cname, fname, "Done.\n");
    return iters;
}


// Like ApplyOverlapDagInverse, but uses the initial value of out
// as a guess. This requires an extra inversion of D_DW(1)^\dag
int ApplyOverlapDagInverseGuess(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Vector *out, Vector *in, Float mass, Float stop_rsd)
{
    const char* fname = "ApplyOverlapDagInverseGuess";

    Fermion_t out_5d[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };
    Fermion_t in_5d[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };
    Fermion_t tmp_5d[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };

    // out_5d = P out
    bfm_d.cps_impexFermion_4d((Float *)out, out_5d, Import);

    // in_5d = P in
    bfm_d.cps_impexFermion_4d((Float *)in, in_5d, Import);

    // Contruct initial guess for 5D solve
    // tmp_5d = D_DW(1)^{\dag-1} P out
    bfm_d.set_mass(1.0);
    if (use_mixed_solver) bfm_f.set_mass(1.0);
    bfm_d.residual = stop_rsd;
    if (use_mixed_solver) bfm_f.residual = 1e-5;
    int iters_PV;
    VRB.Result(cname, fname, "Doing Pauli-Villars solve to get 5D guess\n");
#ifdef BFM_GPARITY
   ERR.General(cname,fname,"mixed_cg::threaded_cg_mixed_Mdag_guess not implemented for BFM with Gparity\n");
#else
#pragma omp parallel
    {
        bfm_d.set_zero(tmp_5d[Even]);
        bfm_d.set_zero(tmp_5d[Odd]);
	iters_PV = use_mixed_solver ?
	    mixed_cg::threaded_cg_mixed_Mdag_guess(tmp_5d, out_5d, bfm_d, bfm_f, 5) :
	    bfm_d.CGNE_Mdag(tmp_5d, out_5d);
    }
#endif

    // tmp_5d = D_DW(m)^{\dag-1} P in
    // Above we constructed the initial guess for the inversion.
    bfm_d.set_mass(mass);
    if (use_mixed_solver) bfm_f.set_mass(mass);
    bfm_d.residual = stop_rsd;
    if (use_mixed_solver) bfm_f.residual = 1e-5;
    int iters_main;
    VRB.Result(cname, fname, "Doing main solve.\n");
#ifdef BFM_GPARITY
   ERR.General(cname,fname,"mixed_cg::threaded_cg_mixed_Mdag_guess not implemented for BFM with Gparity\n");
#else
#pragma omp parallel 
    {
	iters_main = use_mixed_solver ?
	    mixed_cg::threaded_cg_mixed_Mdag_guess(tmp_5d, in_5d, bfm_d, bfm_f, 5) :
	    bfm_d.CGNE_Mdag(tmp_5d, in_5d);
    }
#endif

    // Finally, out_5d = D_DW(1)^\dag D_DW(m)^{\dag-1} P in
    bfm_d.set_mass(1.0);
#pragma omp parallel
    {
        bfm_d.G5D_Munprec(tmp_5d, out_5d, DaggerYes);
    }

    // Extract the 4D part of the 5D result
    bfm_d.cps_impexFermion_4d((Float *)out, out_5d, Export);

    bfm_d.freeFermion(out_5d[Even]);
    bfm_d.freeFermion(out_5d[Odd]);
    bfm_d.freeFermion(in_5d[Even]);
    bfm_d.freeFermion(in_5d[Odd]);
    bfm_d.freeFermion(tmp_5d[Even]);
    bfm_d.freeFermion(tmp_5d[Odd]);

    int iters = iters_PV + iters_main;
    VRB.Result(cname, fname, "iters_PV   = %d\n", iters_PV);
    VRB.Result(cname, fname, "iters_main = %d\n", iters_main);
    VRB.Result(cname, fname, "final total iters = %d at Ls = %d\n", iters, bfm_d.Ls);
    VRB.Result(cname, fname, "Done.\n");
    return iters;
}

void AutofillBfmarg(bfmarg &arg); // defined in f_bfm.C

int InvertOverlapDefectCorrection(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Vector *out, Vector *in, Float mass, DagType dag, bfmarg cheap_approx, Matrix *gauge_field, int num_dc_steps,
    Float cheap_solve_stop_rsd, Float exact_solve_stop_rsd)
{
    const char* fname = "InvertOverlapDefectCorrection()";

    VRB.Result("", fname, "Start!\n");

    bfm_d.comm_end();

    bfm_evo<double> bfm_cheap_d;
    bfm_evo<float> bfm_cheap_f;
    AutofillBfmarg(cheap_approx); // make sure various fields are filled in
    bfm_cheap_d.init(cheap_approx);
    bfm_cheap_d.cps_importGauge((Float *)gauge_field);
    if (use_mixed_solver) {
	bfm_cheap_d.comm_end();
	bfm_cheap_f.init(cheap_approx);
	bfm_cheap_f.cps_importGauge((Float *)gauge_field);
	bfm_cheap_f.comm_end();
	bfm_cheap_d.comm_init();
    }

    int iters;
    int total_iters = 0;
    int total_iters_times_Ls = 0;

    size_t f_size = 24 * GJP.VolNodeSites();
    Vector *approx_sol = (Vector*)smalloc(f_size * sizeof(Float), "approx_sol", fname, "");
    Vector *residual = (Vector*)smalloc(f_size * sizeof(Float), "residual", fname, "");
    Vector *tmp = (Vector*)smalloc(f_size * sizeof(Float), "tmp", fname, "");

    // TODO: try to use initial guess?
    out->VecZero(f_size);
    residual->CopyVec(in, f_size);
    Float norm2_residual = residual->NormSqGlbSum(f_size);

    VRB.Result("", fname, "Starting norm2_residual = %e\n", norm2_residual);

    for (int dc_step = 0; dc_step < num_dc_steps; dc_step++) {
	// Approximately invert Dov on residual
	// approx_sol = Dov'^{-1} res
	VRB.Result("", fname, "Doing cheap solve (dc_step = %d of %d)\n", dc_step, num_dc_steps);
        if (dag == DAG_NO) {
	    iters = ApplyOverlapInverse(bfm_cheap_d, bfm_cheap_f, use_mixed_solver,
	        approx_sol, residual, cheap_approx.mass, cheap_solve_stop_rsd);
        } else {
	    iters = ApplyOverlapDagInverse(bfm_cheap_d, bfm_cheap_f, use_mixed_solver,
	        approx_sol, residual, cheap_approx.mass, cheap_solve_stop_rsd);
        }
        VRB.Result(cname, fname, "cheap inverse took %d iters at Ls = %d\n", iters, bfm_cheap_d.Ls);
        total_iters += iters;
        total_iters_times_Ls += iters * bfm_cheap_d.Ls;

	out->VecAddEquVec(approx_sol, f_size);

        if (dc_step < num_dc_steps - 1) {
	    // compute new true residual for use in next defect correction step

	    bfm_cheap_d.comm_end();
	    bfm_d.comm_init();

	    // tmp = Dov out
            if (dag == DAG_NO) {
	        iters = ApplyOverlap(bfm_d, bfm_f, use_mixed_solver,
	            tmp, out, mass, exact_solve_stop_rsd);
            } else {
	        iters = ApplyOverlapDag(bfm_d, bfm_f, use_mixed_solver,
	            tmp, out, mass, exact_solve_stop_rsd);
            }
            VRB.Result(cname, fname, "ApplyOverlap took %d iters at Ls = %d\n", iters, bfm_d.Ls);
            total_iters += iters;
            total_iters_times_Ls += iters * bfm_d.Ls;

	    bfm_d.comm_end();
	    bfm_cheap_d.comm_init();

	    // residual = in - Dov out
	    residual->FTimesV1MinusV2(1.0, in, tmp, f_size);

	    norm2_residual = residual->NormSqGlbSum(f_size);
	    VRB.Result("", fname, "After defect correction step #%d, norm2_residual = %e\n", dc_step, norm2_residual);
        }
    }

    bfm_cheap_d.end();
    if (use_mixed_solver) {
	bfm_cheap_f.end();
    }
    bfm_d.comm_init();

    // Final cleanup solve. We use the approximate solution
    // we have computed as the initial guess.
    VRB.Result("", fname, "Doing final cleanup solve\n");
    if (dag == DAG_NO) {
        iters = ApplyOverlapInverseGuess(bfm_d, bfm_f, use_mixed_solver,
	    out, in, mass, exact_solve_stop_rsd);
    } else {
        iters = ApplyOverlapDagInverseGuess(bfm_d, bfm_f, use_mixed_solver,
	    out, in, mass, exact_solve_stop_rsd);
    }
    total_iters += iters;
    total_iters_times_Ls += iters * bfm_d.Ls;
    VRB.Result("", fname, "Final cleanup solve took %d iters at Ls = %d\n", iters, bfm_d.Ls);

    int total_normalized_iters = total_iters_times_Ls / bfm_d.Ls;
    VRB.Result(cname, fname, "total iters = %d\n", total_iters);
    VRB.Result(cname, fname, "total iters*Ls = %d\n", total_iters_times_Ls);
    VRB.Result(cname, fname, "total iters normalized to Ls of %d = %d\n", bfm_d.Ls, total_normalized_iters);

    sfree(approx_sol, "approx_sol", fname, "");
    sfree(residual, "residual", fname, "");
    sfree(tmp, "tmp", fname, "");

    VRB.Result("", fname, "Done!\n");
    return total_iters;
}

// For use in MADWF
//
// Does 
//   Pc = D_DW(1)^{-1} b
//   c0 = [P^{-1} Pc]_0 = [P^{-1} D_DW(1)^{-1} b]_0
static void Convert5dRhsTo4dRhs(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Vector *c0, Fermion_t Pc[2], Fermion_t b[2], Float pv_stop_rsd)
{
    bfm_d.set_mass(1.0);
    if (use_mixed_solver) bfm_f.set_mass(1.0);
    bfm_d.residual = pv_stop_rsd;
    if (use_mixed_solver) bfm_f.residual = 1e-5;
#pragma omp parallel
    {
        // D1inv_in = D_DW(1)^{-1} rhs_5d
	bfm_d.set_zero(Pc[Even]);
	bfm_d.set_zero(Pc[Odd]);
	use_mixed_solver ?
	    mixed_cg::threaded_cg_mixed_M(Pc, b, bfm_d, bfm_f, 5) :
	    bfm_d.CGNE_M(Pc, b);

    }

    // rhs_4d = [P^{-1} D_DW(1)^{-1} in]_0
    bfm_d.cps_impexFermion_4d((Float *)c0, Pc, Export);
}


// For use in MADWF.
// From the approximate 4D solution, derive the full approximate 5D solution.
// Computes 
//     ( -Dov(m) y0 )                          ( -y0 )
//   P (     y1     ) = D_DW(1)^{-1} D_DW(m) P (  c1 )
//     (     y2     )                          (  c2 )
//     (     ...    )                          ( ... )
// and then sets
//         ( y0 )
//   x = P ( y1 )
//         ( y2 )
//         ( y3 )
static void Reconstruct5dSol(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Fermion_t x[2], Vector *y0, Fermion_t Pc[2], 
    Float mass, Float pv_stop_rsd)
{
    Fermion_t tmp[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };
    Fermion_t Dm_tmp[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };
    
    // Construct tmp = P transpose(-y0, c1, c2, c3, ...)
#pragma omp parallel 
    {
	bfm_d.copy(tmp[Even], Pc[Even]);
	bfm_d.copy(tmp[Odd], Pc[Odd]); // now tmp = c
	bfm_d.scale(tmp[Even], -1.0); 
	bfm_d.scale(tmp[Odd], -1.0); // now tmp = -c
    }
    bfm_d.cps_impexFermion_4d((Float *)y0, tmp, Import, false); // now tmp = P transpose(y0, -c1, -c2, ...)
#pragma omp parallel 
    {
	bfm_d.scale(tmp[Even], -1.0); 
	bfm_d.scale(tmp[Odd], -1.0); // now tmp = P transpose(-y0, c1, c2, ... )
    }

    // Compute Dm_tmp = D_DW(m) tmp
    bfm_d.set_mass(mass);
    if (use_mixed_solver) bfm_f.set_mass(mass);
#pragma omp parallel
    {
        bfm_d.G5D_Munprec(tmp, Dm_tmp, DaggerNo);
    }

    // Compute D1inv_Dm_tmp = D_DW(1)^{-1} Dm_tmp = D_DW(1)^{-1} D_DW(m) tmp
    bfm_d.set_mass(1.0);
    if (use_mixed_solver) bfm_f.set_mass(1.0);
    bfm_d.residual = pv_stop_rsd;
    if (use_mixed_solver) bfm_f.residual = 1e-5;
#pragma omp parallel
    {
        bfm_d.set_zero(x[Even]);
        bfm_d.set_zero(x[Odd]);
	use_mixed_solver ?
	    mixed_cg::threaded_cg_mixed_M(x, Dm_tmp, bfm_d, bfm_f, 5) :
	    bfm_d.CGNE_M(x, Dm_tmp);
    }

    // Now x is equal to P y except that the part corresponding to 
    // y0 is wrong. Fix that by overwriting the y0 part of x:
    bfm_d.cps_impexFermion_4d((Float *)y0, x, Import, false);
    // Now x = P y

    bfm_d.freeFermion(tmp[Even]);
    bfm_d.freeFermion(tmp[Odd]);
    bfm_d.freeFermion(Dm_tmp[Even]);
    bfm_d.freeFermion(Dm_tmp[Odd]);
}
    

// See Hantao's thesis
// 
// If we need to solve
//     D_DW(m) x = b
// We construct the equivalent problem
//     P^{-1} D_DW(1)^{-1} D_DW(m) P [P^{-1} x] x = P^{-1} D_DW(1)^{-1} b
// Then we solve the first row of this 5d equation, which is
//     Dov(m) y_0 = c_0
// where
//     y = P^{-1} x
//     c = P^{-1} D_DW(1)^{-1} b
// This 4D problem can be solved with InvertOverlap. To speed things up we can
// use a cheap approximate version of Dov(m) in InvertOverlap. Then we need to recover
// y from y_0 and x from y. 
int MADWF_CG_M(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Fermion_t x[2],
    Fermion_t b[2],
    Float mass,
    Matrix *gauge_field,
    Float exact_solve_stop_rsd,
    MADWFParams madwf_params,
    InverterType itype)
{
    const char* fname = "MADWF_CG_M()";

    int iters = 0;

    VRB.Result("", fname, "Start!\n");

    // Make bfm objects for the cheap operator
    bfm_d.comm_end();
    bfm_evo<double> bfm_cheap_d;
    bfm_evo<float> bfm_cheap_f;
    AutofillBfmarg(madwf_params.cheap_approx); // make sure various fields are filled in
    bfm_cheap_d.init(madwf_params.cheap_approx);
    bfm_cheap_d.cps_importGauge((Float *)gauge_field);
    if (use_mixed_solver) {
	bfm_cheap_d.comm_end();
	bfm_cheap_f.init(madwf_params.cheap_approx);
	bfm_cheap_f.cps_importGauge((Float *)gauge_field);
	bfm_cheap_f.comm_end();
	bfm_cheap_d.comm_init();
    }
    bfm_cheap_d.comm_end();
    bfm_d.comm_init();

    Fermion_t Pc[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };
    Fermion_t residual[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };
    Fermion_t z[] = { bfm_d.allocFermion(), bfm_d.allocFermion() };

    size_t f_size_4d = 24 * GJP.VolNodeSites();
    Vector *c0 = (Vector*)smalloc(f_size_4d * sizeof(Float), "c0", fname, ""); 
    Vector *guess_y0 = (Vector*)smalloc(f_size_4d * sizeof(Float), "guess_y0", fname, ""); 

    // Compute the initial 5D residual
    bfm_d.set_mass(mass);
#pragma omp parallel
    {
	bfm_d.G5D_Munprec(x, z, DaggerNo); // z = D_DW(m) * x
	bfm_d.axpy(residual, z, b, -1.0); // residual = b - D_DW(m) * x
    }
    
    // make sure certain residuals aren't pointlessly strict
    // if (madwf_params.cheap_solve_stop_rsd < exact_solve_stop_rsd) madwf_params.cheap_solve_stop_rsd = exact_solve_stop_rsd;
    // if (madwf_params.exact_pv_stop_rsd < exact_solve_stop_rsd) madwf_params.exact_pv_stop_rsd = exact_solve_stop_rsd;

    // Now we will use the cheap operator to iteratively compute an approximation to 
    //   D_DW(m)^{-1} residual
    for (int dc_step = 0; dc_step < madwf_params.num_dc_steps; dc_step++) {
	VRB.Result(cname, fname, "Starting defect correction step #%d of %d\n", dc_step + 1, madwf_params.num_dc_steps);

	Float norm2_residual;
#pragma omp parallel
	{
	    norm2_residual = bfm_d.norm(residual);
	}
	VRB.Result(cname, fname, "Current norm2 of residual is %0.16e\n", norm2_residual);

	// Translate the 5D problem
	//    D_DW(m) z = residual
	// into a 4D problem.
	// First construct the RHS of the 4D problem.
	//   Pc = D_DW(1)^{-1} residual
	//   c0 = [P^{-1} Pc]_0 = [P^{-1} D_DW(1)^{-1} residual]_0
	VRB.Result(cname, fname, "Converting 5D RHS to 4D RHS.\n");
	Convert5dRhsTo4dRhs(bfm_d, bfm_f, use_mixed_solver, c0, Pc, residual, madwf_params.exact_pv_stop_rsd);

	// now invert the cheap overlap approximation
	VRB.Result(cname, fname, "Doing cheap overlap inversion.\n");
	bfm_d.comm_end();
	bfm_cheap_d.comm_init();
	// guess_y0 = Dov'(m)^{-1} c0
	iters += ApplyOverlapInverse(bfm_cheap_d, bfm_cheap_f, use_mixed_solver, 
	    guess_y0, c0, madwf_params.cheap_approx.mass, madwf_params.cheap_solve_stop_rsd, itype);
	bfm_cheap_d.comm_end();
	bfm_d.comm_init();
	
	// Now reconstruct the guess for z, the solution to the 5D problem
	VRB.Result(cname, fname, "Reconstructing 5D guess.\n");
	Reconstruct5dSol(bfm_d, bfm_f, use_mixed_solver, z, guess_y0, Pc, mass, madwf_params.exact_pv_stop_rsd);

	// z is approximately the vector we need to add to x to get the full solution
	// to the overall 5D problem. So update x with
	//   x = x + z
#pragma omp parallel
	{
	    bfm_d.axpy(x, x, z, 1.0);
	}	

	// If we are going to do more steps of defect correction, compute the new residual
	if (dc_step < madwf_params.num_dc_steps - 1) {
	    bfm_d.set_mass(mass);
#pragma omp parallel
	    {
		bfm_d.G5D_Munprec(x, z, DaggerNo); // z = D_DW(m) * x
		bfm_d.axpy(residual, z, b, -1.0); // residual = b - D_DW(m) * x
	    }
	} 
    }

    sfree(c0, "c0", fname, cname);
    sfree(guess_y0, "guess_y0", fname, cname);

    bfm_d.freeFermion(Pc[Even]);
    bfm_d.freeFermion(Pc[Odd]);
    bfm_d.freeFermion(residual[Even]);
    bfm_d.freeFermion(residual[Odd]);
    bfm_d.freeFermion(z[Even]);
    bfm_d.freeFermion(z[Odd]);

    bfm_cheap_d.end();
    if (use_mixed_solver) {
	bfm_cheap_f.end();
    }

    // Now we clean up our approximate solution x. The approximate solution serves
    // as a guess for a regular 5D inversion. 
    VRB.Result(cname, fname, "Doing final 5D solve.\n");
    bfm_d.set_mass(mass);
    if (use_mixed_solver) bfm_f.set_mass(mass);
    bfm_d.residual = exact_solve_stop_rsd;
    if (use_mixed_solver) bfm_f.residual = 1e-5;
    Float final_solve_iters;
#pragma omp parallel
    {
	final_solve_iters = use_mixed_solver ?
	    mixed_cg::threaded_cg_mixed_M(x, b, bfm_d, bfm_f, 5) :
	    bfm_d.CGNE_M(x, b);
    }
    iters += final_solve_iters;

    return iters;
}



    



    


CPS_END_NAMESPACE

#endif
