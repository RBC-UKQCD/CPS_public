#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_hmd.h
//
// Header file for the AlgHmd class and its derived classes
// AlgHmdR and AlgHmcPhi. 
// The type of glue or fermion is given as
// an argument of type Lattice& to the constructor.
//
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_HMD_H
#define INCLUDED_ALG_HMD_H

CPS_END_NAMESPACE
#include<util/lattice.h>
#include<util/smalloc.h>
#include<util/pmalloc.h>
#include<alg/alg_base.h>
#include<alg/common_arg.h>
#include<alg/hmd_arg.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
//
// AlgHmd is derived from Alg and is relevant to the Hybrid  
// Molecular Dynamics algorithms.
//
//------------------------------------------------------------------
class AlgHmd : public Alg
{
 private:
    char *cname;

 protected:

    HmdArg *hmd_arg;
        // The argument structure for the hmc phi algorithm
  
    int g_size;       
        // Node size of the gauge field

    int Ncb;
        // Number of checkerboards on which the fermion field is defined

    Matrix* mom;
        // Conjugate momentum (traceless antihermitian) 

 public:

  AlgHmd(Lattice& latt, CommonArg *c_arg, HmdArg *arg);

  virtual ~AlgHmd();

};


//------------------------------------------------------------------
//
// AlgHmcPhi is derived from AlgHmd and is relevant to the Hybrid  
// Molecular Dynamics HMC Phi algorithm.
//
//------------------------------------------------------------------
class AlgHmcPhi : public AlgHmd
{
 private:
    char *cname;


 protected:

    int n_frm_masses;     
        // The number of dynamical fermion masses.

    int n_bsn_masses;     
        // The number of dynamical boson masses.

    int f_size;       
        // Node checkerboard size of the fermion field

    CgArg **frm_cg_arg;
        // pointer to an array of CG argument structures
	// relevant to fermions.

    CgArg **bsn_cg_arg;
        // pointer to an array of CG argument structures
	// relevant to bosons.

    Vector** phi;
        // Pseudo fermion field phi (checkerboarded).

    Vector** bsn;
        // Boson field bsn (checkerboarded).

    Matrix* gauge_field_init;
        // Initial gauge field needed if the evolved one is rejected

    Vector** frm1;
    Vector** frm2;
        // 2 general purpose fermion/boson field arrays (checkerboarded).

    Vector** cg_sol_cur; 
    Vector** cg_sol_prev;
        // Temporary pointers to the current and previous CG solution.

    Float *h_f_init;    
        // Initial fermion Hamiltonian (one for each mass)

    Float *h_f_final;   
        // Final fermion Hamiltonian (one for each mass)

    Float *delta_h_f;   
        // Final-Init fermion Hamiltonian (one for each mass)

    Float *h_b_init;    
        // Initial boson Hamiltonian (one for each mass)

    Float *h_b_final;   
        // Final boson Hamiltonian (one for each mass)

    Float *delta_h_b;   
        // Final-Init boson Hamiltonian (one for each mass)

 public:

  AlgHmcPhi(Lattice& latt, CommonArg *c_arg, HmdArg *arg);

  virtual ~AlgHmcPhi();

  void run(void);
};


//------------------------------------------------------------------
//
// AlgHmdR is derived from AlgHmd and is relevant to the Hybrid  
// Molecular Dynamics R algorithm. Boson fields are simulated as
// fermion fields with negative flavor number.
//
//------------------------------------------------------------------
class AlgHmdR : public AlgHmd
{
 private:
    char *cname;

 protected:
    int n_frm_masses;     
        // The number of dynamical fermion masses.

    Float *flavor_time_step;
        // Pointer to an array of size n_frm_masses + 1.
        // Each entry contains the time_step for each 
        // intermediate (one for each flavor plus the 
        // midpoint one) evolution of the gauge field.

    int f_size;       
        // Node checkerboard size of the fermion field

    CgArg **frm_cg_arg;
        // pointer to an array of CG argument structures
	// relevant to fermions.

    Vector** phi;
        // Pseudo fermion field phi (checkerboarded).

    Vector* frm1;
    Vector* frm2;
        // 2 general purpose fermion fields (checkerboarded).

 public:

  AlgHmdR(Lattice& latt, CommonArg *c_arg, HmdArg *arg);

  virtual ~AlgHmdR();

  void run(void);
};


#endif




CPS_END_NAMESPACE
