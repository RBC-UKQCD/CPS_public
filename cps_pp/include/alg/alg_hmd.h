#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief  Definitions of the AlgHmd class and derived classes.

  $Id: alg_hmd.h,v 1.18 2005-11-15 06:25:20 chulwoo Exp $
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

    int alloc_sum;

 public:

  AlgHmd(Lattice& latt, CommonArg *c_arg, HmdArg *arg);


  virtual ~AlgHmd();

  //! Performs a single trajectory of the HMD algorithm. 
  virtual Float run() = 0;

};


//------------------------------------------------------------------
//! A class implementing the Hybrid Monte Carlo algorithm.
/*!
  This evolves a gauge field by a single iteration of
  the standard HMC algorithm (phi algorithm), \e i.e. a molecular dynamics
  trajectory followed by a metropolis accept/reject step.

  This implementation uses a Upqp leapfrog integration scheme and a
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

 public:

  AlgHmcPhi(Lattice& latt, CommonArg *c_arg, HmdArg *arg);

  virtual ~AlgHmcPhi();

  //! Performs a single  HMC trajectory.  
  Float run(void);
};

//------------------------------------------------------------------
//! A class implementing the Hybrid Monte Carlo algorithm.
/*!
  This evolves a gauge field by a single iteration of
  the HMC algorithm, \e i.e. a molecular dynamics
  trajectory followed by a metropolis accept/reject step.

  This implementation uses a Uqpq leapfrog integration scheme and a minimum 
  residual chronological method to choose the starting guess for the solver.  

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
class AlgHmcQPQ : public AlgHmd
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

    Vector** cg_sol_prev;
    Vector*** cg_sol;
    //!< Array holding the conjugate gradient solutions
    Vector*** vm;
    //!< Array containing the orthogonal basis used by the chronological preconditioner

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

 public:

  AlgHmcQPQ(Lattice& latt, CommonArg *c_arg, HmdArg *arg);

  virtual ~AlgHmcQPQ();

  //! Performs a HMC trajectory.
  Float run(void);
};

//------------------------------------------------------------------
//! A class implementing the Rational Hybrid Monte Carlo algorithm.
/*!
  The main method evolves a gauge field along a molecular dynamics
  trajectory followed by a metropolis accept/reject step using the RHMC
  algorithm.

  The algorithm is configurable to include dynamical fermions
  of several masses, each mass having its own set of pseudofermion fields
  and solver parameters. The number of degenerate flavours at each mass
  depends on the type of fermion action used.

  As with AlgHmcPhi, one can also have bosons in the action.

  The rational approximation can be generated dynamically, but the
  correct constructor must be used and the code should be configured with
  the --enable-gmp flag to use the GNU multiprecision library.

  \note When used with Asqtad fermions, the local lattice size has
  to be greater than 2 in all directions.
*/
//------------------------------------------------------------------
class AlgHmcRHMC : public AlgHmd
{
 private:
    char *cname;

    //! Does the main work of the constructors.
    void init();

    //! Renormalises the shifts into the mass for staggered fermions
    void massRenormalise(Float *mass, Float *trueMass, int degree, 
			 Float *shift, MassRenormaliseDir direction);

 protected:

    int n_frm_masses;     
        //!< The number of dynamical fermion masses.

    int n_bsn_masses;     
        //!< The number of dynamical boson masses.

    int f_sites;       
    int f_vec_count;       
    int f_count;       
    int f_size;       
    //!< The size of a fermion field.
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
    //!< The final gauge field configuration

    unsigned int **rng4d_init;
    unsigned int **rng5d_init;
    //!< The initial random numbers

    Matrix* gauge_field_final;
    //!< The initial gauge field configuration

    Vector* frm;

    Vector** frmn;
    //!< Array of vectors
    /*!< These will hold the solutions from the solves. */

    Vector** frmn_d;
    //!< Array of vectors
    /*!< These will hold the solutions from the solves multiplied by
      the D-slash operator. */ 

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

    int total_size;
    //!< The sum of the rational approximation degrees used for the force

    Float *all_res;
    //!< An array holding the residues - used for asqtad force

    EigArg *eig_arg;
    //!< ?

    Float **alpha;

 public:

  //!  Constructor for when the rational approximation is fixed.
  AlgHmcRHMC(Lattice& latt, CommonArg *c_arg, HmdArg *arg);

  //!   Constructor for when the rational approximations are to be generated dynamically.
  AlgHmcRHMC(Lattice& latt, CommonArg *c_arg, HmdArg *arg, EigArg *e_arg);

  virtual ~AlgHmcRHMC();

  //! Performs a single RHMC trajectory
  Float run(void);

  //! Automatic generation of the rational approximation.
  void generateApprox(HmdArg*);

  //! Dynamical generation of the rational approximation.
  void dynamicalApprox();

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
      field updates in the R algorithm for each dynamical mass.
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
    Vector* frm2;
    //!< Fermion work vectors
    
 public:

  AlgHmdR(Lattice& latt, CommonArg *c_arg, HmdArg *arg);

  virtual ~AlgHmdR();

  //! Performs a single HMD trajectory.
  Float run(void);
};

//------------------------------------------------------------------
//! A class implementing the R2 algorithm.
/*!
  This evolves a gauge field by a single iteration of the standard HMD
  algorithm, \e i.e. a molecular dynamics trajectory with intermediate
  gauge field updates to take care of the arbitrary numbers of
  dynamical flavours. The initial guess in the solver is a random
  gaussian vector.
  
  The R2 algorithm is specifically for simulating 2+1 flavours of
  dynamical staggered fermions of different masses.  It uses the
  generalised multimass solver to minimise the number of cg iterations
  performed.  The algorithm can only be used when the mass matrix is a
  scalar multiple of the identity.

  \ingroup alg
*/
//------------------------------------------------------------------
class AlgHmdR2 : public AlgHmd
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

    Float *force_coeff;
    //!< The coefficient of each force passed to RHMC_EvolveFforce
    /*!<
      Array of coefficients = flavours * flavour_coeff
    */

    int f_sites;       
    int f_vec_count;       
    int f_count;       
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

    Vector** frmn;
    //!< Array of vectors
    /*!< These will hold the solutions from the solves. */

    Vector** frmn_d;
    //!< Array of vectors
    /*!< These will hold the solutions from the solves multiplied by
      the D-slash operator. */ 

    Float *shift;
    //!< The renormalised shift

    int light;
    //!< The index of the lightest mass

    int heavy;
    //!< The index of the heaviest mass

 public:

  AlgHmdR2(Lattice& latt, CommonArg *c_arg, HmdArg *arg);

  virtual ~AlgHmdR2();

  //! Performs a single HMD trajectory.
  Float run(void);
};

#endif

CPS_END_NAMESPACE
