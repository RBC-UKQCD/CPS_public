#ifndef INCLUDED_EFF_OVERLAP_H__
#define INCLUDED_EFF_OVERLAP_H__

#ifdef USE_BFM

#include <util/lattice/bfm_evo.h>
#include <util/vector.h>

CPS_START_NAMESPACE

int ApplyOverlap(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Vector *out, Vector *in, Float mass, Float pv_stop_rsd);

int ApplyOverlapInverse(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Vector *out, Vector *in, Float mass, Float stop_rsd,
    InverterType itype = CG);

int ApplyOverlapDag(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Vector *out, Vector *in, Float mass, Float pv_stop_rsd);

int ApplyOverlapDagInverse(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Vector *out, Vector *in, Float mass, Float stop_rsd);


int ApplyOverlapInverseGuess(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Vector *out, Vector *in, Float mass, Float stop_rsd);

int ApplyOverlapDagInverseGuess(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Vector *out, Vector *in, Float mass, Float stop_rsd);

int InvertOverlapDefectCorrection(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Vector *out, Vector *in, Float mass, DagType dag, bfmarg cheap_approx, Matrix *gauge_field, int num_dc_steps,
    Float cheap_solve_stop_rsd, Float exact_solve_stop_rsd);



/*void Convert5dRhsTo4dRhs(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Vector *c0, Fermion_t Pc[2], Fermion_t b[2], Float pv_stop_rsd);

void Reconstruct5dSol(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Fermion_t x[2], Vector *y0, Fermion_t Pc[2],
    Float mass, Float pv_stop_rsd);*/

struct MADWFParams
{
    bfmarg cheap_approx;
    int num_dc_steps;
    Float cheap_solve_stop_rsd;
    Float exact_pv_stop_rsd;
};
typedef struct MADWFParams MADWFParams;

int MADWF_CG_M(bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, bool use_mixed_solver,
    Fermion_t x[2], 
    Fermion_t b[2], 
    Float mass, 
    Matrix *gauge_field, 
    Float exact_solve_stop_rsd, 
    MADWFParams madwf_params,
    InverterType itype = CG);




CPS_END_NAMESPACE

#endif 

#endif
