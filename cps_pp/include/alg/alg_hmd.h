#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief  Definitions of the AlgHmd class and derived classes.

  $Id: alg_hmd.h,v 1.3 2003-08-11 22:15:45 mike Exp $
*/
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_HMD_H
#define INCLUDED_ALG_HMD_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/hmd_arg.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
//! A virtual base class for implementing Molecular Dynamics algorithms.
/*!
  \ingroup alg
  \todo The derived classes should inherit more.
*/
//------------------------------------------------------------------
class AlgHmd : public Alg
{
 private:
    char *cname;
   
 protected:

    HmdArg *hmd_arg;
        //!< A structure containing the algorithm parameters
  
    int g_size;       
    //!< The size of the gauge field.
    /*!<
      The size of the gauge field on the local lattice in terms of the
      total number of floating point numbers.
    */

    int Ncb;
    //!< The number of parities on which the pseudofermion field is defined.

    Matrix* mom;
    //!< The (traceless antihermitian) conjugate momentum field.
    /*!<
      This is \e i times the actual conjugate momentum.
     */	

 public:

  AlgHmd(Lattice& latt, CommonArg *c_arg, HmdArg *arg);

  virtual ~AlgHmd();

};


//------------------------------------------------------------------
//! A class implementing the Hybrid Monte Carlo algorithm.
/*!
  This evolves a gauge field by a single iteration of
  the standard HMC algorithm (phi algorithm), \e i.e. a molecular dynamics
  trajectory followed by a metropolis accept/reject step.

  This implementation uses the leapfrog integration scheme and a
  two-step chronological method to choose the starting guess for the solver.  

  The algorithm is configurable to include dynamical fermions
  of several masses, each mass having its own set of pseudofermion fields
  and solver parameters. The number of degenerate flavours at each mass
  depends on the type of fermion action used. Similarly
  one can have in addition pseudobosonic 
  \f$ \chi^\dagger M^\dagger M \chi \f$ terms in the hamiltonian contributing
  to the total force on the momentum field. 

  One also has the option of forcing acceptance at the metropolis step.

  A bunch of statistics relating to the performance of the solver, the
  conservation of the hamiltonian and the metropolis acceptance are collected.
  
  \ingroup alg
*/
//------------------------------------------------------------------
class AlgHmcPhi : public AlgHmd
{
 private:
    char *cname;


 protected:

    int n_frm_masses;     
    //!< The number of dynamical fermion masses.
    /*!
      This is not necessarily the number of dynamical flavours.
     */

    int n_bsn_masses;     
    //!< The number of dynamical boson masses.

    int f_size;       
    //!< The size of the pseudofermion (and similar) fields.
    /*!< The size is given in terms of the total number of floating point
      numbers in the field on the local lattice, taking into account whether
      or not the field is defined on just a single parity.
    */

    CgArg **frm_cg_arg;
    //!< Pointer to an array of structures containing solver parameters.
    /*!<
      These are the parameters corresponding to each of the dynamical fermion
      masses.      
     */

    CgArg **bsn_cg_arg;
    //!< Pointer to an array of structures containing solver parameters.
    /*!<
      These are the parameters corresponding to each of the dynamical boson
      masses.      
     */

    Vector** phi;
    //!< Pseudofermion fields
    /*!< One for each mass */

    Vector** bsn;
    //!< Boson (pseudoboson?) fields 
    /*!< One for each mass */
    
    Matrix* gauge_field_init;
    //!< The initial gauge field configuration

    Vector** frm1;
    //!< Array of general purpose fields
    Vector** frm2;
    //!< Another array of general purpose fields    

    Vector** cg_sol_cur;
    //!< Pointer to the most recent solution produced by the solver.
    /*!< One for each mass */    
    Vector** cg_sol_prev;
    //!< Pointer to the most recent solution produced by the solver.
    /*!< One for each mass */    


    Float *h_f_init;    
    //!< The initial value of the pseudofermion action.
    /*!< The value at the start of the trajectory. One for each mass */

    Float *h_f_final;   
    //!< The final value of the pseudofermion action.
    /*!< The value at the end of the trajectory. One for each mass */

    Float *delta_h_f;   
    //!< The change in the value of the pseudofermion action.
    /*!< The final value - the initial value. One for each mass */

    Float *h_b_init;    
    //!< The initial value of the boson action.
    /*!< The value at the start of the trajectory. One for each mass */

    Float *h_b_final;   
    //!< The final value of the boson action.
    /*!< The value at the end of the trajectory. One for each mass */

    Float *delta_h_b;   
    //!< The change in the value of the boson action.
    /*!< The final value - the initial value. One for each mass */

    int *FRatDeg;
    Float *FRatNorm;
    Float **FRatConst;
    Float **FRatPole;

    int *HBRatDeg;
    Float *HBRatNorm;
    Float **HBRatConst;
    Float **HBRatPole;

 public:

  AlgHmcPhi(Lattice& latt, CommonArg *c_arg, HmdArg *arg);

  virtual ~AlgHmcPhi();

  void run(void);
};

//------------------------------------------------------------------
//
// AlgHmcRHMC is derived from AlgHmd and is relevant to the Rational Hybrid  
// Monte Carlo algorithm.
//
//------------------------------------------------------------------
class AlgHmcRHMC : public AlgHmd
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

    int *FRatDeg;
    Float *FRatNorm;
    Float **FRatConst;
    Float **FRatPole;

    int *HBRatDeg;
    Float *HBRatNorm;
    Float **HBRatConst;
    Float **HBRatPole;

 public:

  AlgHmcRHMC(Lattice& latt, CommonArg *c_arg, HmdArg *arg);

  virtual ~AlgHmcRHMC();

  void run(void);
};

//------------------------------------------------------------------
//! A class implementing the Hybrid Molecular Dynamics (R) algorithm.
/*!
  This evolves a gauge field by a single iteration of
  the standard HMD algorithm, \e i.e. a molecular dynamics trajectory
  with intermediate gauge field updates to take care of the arbitrary numbers
  of dynamical flavours. The initial guess in the solver is a random gaussian
  vector.
  
  The algorithm is configurable to include dynamical fermions
  of several masses, each with a different number of flavours.
  For each mass there is a set of solver parameters. 

  A bunch of statistics relating to the performance of the solver, the
  conservation of the hamiltonian and the metropolis acceptance are collected.
  
  \ingroup alg
*/
//------------------------------------------------------------------
class AlgHmdR : public AlgHmd
{
 private:
    char *cname;

 protected:
    int n_frm_masses;     
    //!< The number of dynamical fermion masses.
    /*!<
      This is not necessarily the number of dynamical flavours.
     */

    Float *flavor_time_step;
    //!< A tricky thing to describe succinctly.
    /*!<
      This is an array of the time steps used for the intermediate gauge
      field updates for each dynamical mass.
      Actually it is the differences between them. At least most of it is.
    */

    int f_size;       
    //!< The size of the pseudofermion (and similar) fields.
    /*!< The size is given in terms of the total number of floating point
      numbers in the field on the local lattice, taking into account whether
      or not the field is defined on just a single parity.
    */

    CgArg **frm_cg_arg;
    //!< Pointer to an array of structures containing solver parameters.
    /*!<
      These are the parameters corresponding to each of the dynamical fermion
      masses.      
     */

    Vector** phi;
    //!< Pseudofermion fields
    /*!< One for each mass */

    Vector* frm1;
    //!< Array of general purpose fields

    Vector* frm2;
    //!< Another array of general purpose fields
    

 public:

  AlgHmdR(Lattice& latt, CommonArg *c_arg, HmdArg *arg);

  virtual ~AlgHmdR();

  //! Performs a single HMD trajectory.
  void run(void);
};


#endif





CPS_END_NAMESPACE
