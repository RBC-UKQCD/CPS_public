#include<config.h>

#ifdef USE_BFM
#include <util/lattice/bfm_evo.h>
#include <util/lattice/bfm_eigcg.h>
#include <util/lattice/fbfm.h>
#endif


CPS_START_NAMESPACE
//Controls the version of the multi-shift algorithm used. The user can define different versions to use in different environments, for example if you
//want to use an approximate method within the molecular dynamics evolution.

//The current environment must be set manually. It defaults to Generic
class MultiShiftCGcontroller{
public:
    enum Mode { SINGLE_PREC,
	       DOUBLE_PREC,
	       SINGLE_PREC_PLUS_OUTER_DEFECT_CORRECTION_LOOP, //Single precision multi-shift followed by single precision restarted defect correction loop over poles
	       SINGLE_PREC_AS_DOUBLE_PREC_GUESS, //Single precision multi-shift followed by double precision multi-shift using the single prec results as a guess (using Osborn's method) 
	       SINGLE_PREC_RESTARTED_AS_DOUBLE_PREC_GUESS, //Restarted single precision multi-shift with defect correction followed by double precision multi-shift using the single prec results as a guess (also using Osborn's method)
	       SINGLE_PREC_RELIABLE_UPDATE_PLUS_OUTER_DEFECT_CORRECTION_LOOP, //Single precision multi-shift with reliable update followed by single precision restarted defect correction loop over poles
	       NMultiShiftCGMode };
    enum Environment { MolecularDynamics, EnergyCalculation, Heatbath, Generic, NMultiShiftEnvironment };
private:
    Mode environ_mode[(int)NMultiShiftEnvironment];
    Environment current_environment;

    double minimum_single_prec_residual;  //For variants with an initial single precision solve, the stopping conditions are set equal to the larger of the double precision residual and this bound. Does not apply
                                          //to the reliable update version.
    int reliable_update_freq; //Used in versions with reliable update
    int max_defect_correction_cycles;
public:
    MultiShiftCGcontroller(): current_environment(Generic), 
			 minimum_single_prec_residual(1e-08), 
			 reliable_update_freq(100),
			 max_defect_correction_cycles(500){
	for(int i=0;i<(int)NMultiShiftEnvironment;i++) environ_mode[i] = DOUBLE_PREC;
    }

    void setEnvironmentMode(const Environment & environ,  const Mode & mode){ environ_mode[(int)environ] = mode; }
    void setEnvironment(const Environment & environ){ current_environment = environ; }
    const Mode& getMode() const{ return environ_mode[(int)current_environment]; }

    void setMinimumSinglePrecResidual(const double &r){ minimum_single_prec_residual = r; }
    void setReliableUpdateFreq(const int &f){ reliable_update_freq = f; }
    void setMaximumDefectCorrectionCycles(const int &c){ max_defect_correction_cycles = c; }

#ifdef USE_BFM
    int MInv(Fermion_t *sol_multi, Fermion_t src, 
	     Float *shift, int Nshift, 
	     Float *mresidual, Float *alpha, int single,
	     bfm_evo<double> &bd,
	     bfm_evo<float> &bf){
	
	const Mode& mode = getMode();
	int iter;

	if(mode == SINGLE_PREC){ //Note, this uses the residuals specified in the cg_arg without modification
#pragma omp parallel
	    {
		Fermion_t src_f = bf.threadedAllocFermion();
		Fermion_t sol_f[Nshift];
		for(int i=0;i<Nshift;i++) sol_f[i] = bf.threadedAllocFermion();

		mixed_cg::threaded_convFermion(src_f,src,bf,bd);
		mixed_cg::switch_comm(bf,bd);
		iter = bf.CGNE_prec_MdagM_multi_shift(sol_f, src_f, shift, alpha, Nshift, mresidual, single);
		mixed_cg::switch_comm(bd,bf);
		for(int i=0;i<Nshift;i++){
		    mixed_cg::threaded_convFermion(sol_multi[i],sol_f[i],bd,bf);
		    bf.threadedFreeFermion(sol_f[i]);
		}
		bf.threadedFreeFermion(src_f);
	    }       
	}else if(mode == DOUBLE_PREC){
#pragma omp parallel
	    {
		iter = bd.CGNE_prec_MdagM_multi_shift(sol_multi, src, shift, alpha, Nshift, mresidual, single);
	    }
	}else if(mode == SINGLE_PREC_PLUS_OUTER_DEFECT_CORRECTION_LOOP){
	    double fresidual[Nshift]; //residuals for initial single prec solve
	    for(int s=0;s<Nshift;s++) fresidual[s] = (mresidual[s] >= minimum_single_prec_residual ? mresidual[s] : minimum_single_prec_residual);
#pragma omp parallel
	    {
		iter = mixed_cg::threaded_cg_mixed_defect_correction_multi_shift_MdagM(sol_multi,src, shift,alpha, bd,bf, Nshift, mresidual, fresidual, single, max_defect_correction_cycles);
	    }
	}else if(mode == SINGLE_PREC_AS_DOUBLE_PREC_GUESS){
	    double fresidual[Nshift]; //residuals for initial single prec solve
	    for(int s=0;s<Nshift;s++) fresidual[s] = (mresidual[s] >= minimum_single_prec_residual ? mresidual[s] : minimum_single_prec_residual);
#pragma omp parallel
	    {
		iter = mixed_cg::threaded_cg_mixed_single_prec_as_guess_multi_shift_MdagM(sol_multi,src, shift,alpha, Nshift, mresidual, fresidual, single, bd, bf);
	    }
	}else if(mode == SINGLE_PREC_RESTARTED_AS_DOUBLE_PREC_GUESS){
	    double fresidual[Nshift]; //residuals for initial single prec solve
	    for(int s=0;s<Nshift;s++) fresidual[s] = (mresidual[s] >= minimum_single_prec_residual ? mresidual[s] : minimum_single_prec_residual);
#pragma omp parallel
	    {
		iter = mixed_cg::threaded_cg_mixed_restarted_multi_shift_MdagM(sol_multi,src, shift,alpha, Nshift, mresidual, fresidual, single, bd, bf, max_defect_correction_cycles);
	    }
	}else if(mode == SINGLE_PREC_RELIABLE_UPDATE_PLUS_OUTER_DEFECT_CORRECTION_LOOP){
	    #pragma omp parallel
	    {
		iter = mixed_cg::threaded_cg_mixed_multi_shift_MdagM_sp_relup_dp_defect_correction(sol_multi,src, shift,alpha, Nshift, mresidual, single, bf, bd, reliable_update_freq, -1, max_defect_correction_cycles);
	    }
	}else ERR.General("_MultiShiftCGargs","MInv(..)","Unknown multi-shift mode\n");
	return iter;
    }
#endif
};

extern MultiShiftCGcontroller MultiShiftController; //global instance (created in fbfm.C)
CPS_END_NAMESPACE

