#include<config.h>

#ifdef USE_BFM
#include <util/lattice/bfm_evo.h>
#include <util/lattice/bfm_eigcg.h>
#include <util/lattice/fbfm.h>
#endif

// Promoting the class to work without BFM. BFM implementation in bfm_mixed_solver.
// Has to be manually reconciled after merging Gparity branch is done!!
#ifndef USE_BFM

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

};

extern MultiShiftCGcontroller MultiShiftController; //global instance (created in fbfm.C)
CPS_END_NAMESPACE
#endif
