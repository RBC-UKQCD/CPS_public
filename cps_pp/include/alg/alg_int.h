#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief  Definitions of the AlgInt class and derived classes.

*/
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_INT_H
#define INCLUDED_ALG_INT_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/alg_meas.h>
#include <alg/cg_stats.h>
#include <util/checksum.h>
#include <alg/int_arg.h>
#include <alg/remez_arg.h>
CPS_START_NAMESPACE

class AlgInt {

 private:
  char *cname;

 public:
  
  AlgInt();
  virtual ~AlgInt();

  //!< method to do heatbath (if necessary)
  virtual void heatbath() = 0;

  //!< run method evolves the integrator
  virtual void evolve (Float dt, int steps) = 0;

  //!< method to calculate energy (if necessary)
  virtual Float energy() = 0;

  //!< method to return the cost of the integrator (wrt cg iterations)
  virtual void cost(CgStats*) = 0;

  //!< method to reverse the direction of evolution (i.e. flip momenta)
  virtual void reverse() = 0;

  //!< method used to reinitialise the integrator
  virtual void init() = 0;

};

//!< 
class AlgIntAB : public AlgInt {

 private:
  char *cname;

 protected:
  AlgInt *A;
  AlgInt *B;
  int A_steps, B_steps;
  int A_calls, B_calls;
  
  IntegratorLevel level; //!< Is this the top level integrator?
  unsigned long step_cnt;
  IntABArg *ab_arg;

 public:
  AlgIntAB(AlgInt &A, AlgInt &B, IntABArg &);
  virtual ~AlgIntAB();

  void heatbath();

  Float energy();

  //!< evolve method evolves the integrator
  virtual void evolve(Float dt, int steps) = 0;

  void cost(CgStats*);
  void reverse();

  //!< Dummy method
  void init();

  //!< AlgIntAB factory
  static AlgIntAB& Create(AlgInt &A, AlgInt &B, IntABArg &ab_arg);
  static void Destroy(AlgIntAB&);

};

//!< An implementation of a leapfrog integrater
class AlgIntLeap : public AlgIntAB {

private:
  char *cname;

public:
  AlgIntLeap(AlgInt &A, AlgInt &B, IntABArg &);
  virtual ~AlgIntLeap();

  //!< evolve method evolves the integrator
  void evolve(Float dt, int steps);
  
};

//!< An implementation of the Omelyan integrater
class AlgIntOmelyan : public AlgIntAB {

private:
  char *cname;
  Float lambda;

public:
  AlgIntOmelyan(AlgInt &A, AlgInt &B, IntABArg &);
  virtual ~AlgIntOmelyan();

  //!< evolve method evolves the integrator
  void evolve(Float dt, int steps);
  
};

//!< Produces the composite integrator of A and B
//!< If [A,B] != [B,A] then obviously order matters!
class AlgIntSum : public AlgIntAB {

private:
  char *cname;

public:
  AlgIntSum(AlgInt &A, AlgInt &B, IntABArg &);
  virtual ~AlgIntSum();

  //!< evolve method evolves the integrator
  void evolve(Float dt, int steps);

};

//!< Super class of all Hamiltonian constituents
class AlgHamiltonian : public AlgInt {

 private:
  char *cname;

 public:
  AlgHamiltonian();
  virtual ~AlgHamiltonian();

  virtual void heatbath() = 0;
  virtual Float energy() = 0;
  virtual void evolve(Float dt, int steps) = 0;
  virtual void cost(CgStats*) = 0;

};

//!< Class describing momentum
class AlgMomentum : public AlgHamiltonian {

 private:
  char *cname;

  //!< This class tracks the MD time since it updates the gauge field
  char *md_time_str;
  Matrix *mom;
  int g_size;

 public:
  AlgMomentum();
  virtual ~AlgMomentum();

  void heatbath();

  Float energy();

  //!< evolve method evolves the integrator
  void evolve(Float dt, int steps);

  void cost(CgStats*);

  Matrix* getMom();

  void reverse();

  //!< Dummy method
  void init();
};

//!< Super class of all possible actions
class AlgAction : public AlgHamiltonian {

 private:
  char *cname;

 protected:
  Matrix *mom;

 public:
  AlgAction();
  AlgAction(AlgMomentum &mom);
  ~AlgAction();

  virtual void heatbath() = 0;
  virtual Float energy() = 0;
  virtual void evolve(Float dt, int steps) = 0;
  virtual void cost(CgStats*) = 0;

  void reverse();

};

//!< Super class of all possible bilinear actions
class AlgActionBilinear : public AlgAction {

 private:
  char *cname;

 protected:
  ActionBilinearArg *bi_arg;

  int n_masses;
  FclassType fermion;

  //!< An array which stores the values of the masses
  Float *mass;

  //!< Maximum number of cg iterations
  int *max_num_iter;

  //!< Number of lattice sites
  int f_sites;
  //!< Number of Vectors in a Vector array
  int f_vec_count;
  //!< Number of Floats in a Vector array
  int f_size;
  //!< Number of checkerboards
  int Ncb;

  //!< The conjugate gradient statistics
  CgStats cg_stats;
  int cg_iter;

  //!< Pseudofermion fields, one for each mass
  Vector **phi;

  int md_steps;

 public:
  AlgActionBilinear();
  AlgActionBilinear(AlgMomentum &, ActionBilinearArg &);
  virtual ~AlgActionBilinear();
  void cost(CgStats*);
  void updateCgStats(CgArg*);

  int getNmass();
  Float getMass(int);
  FclassType getFermion();

  virtual void heatbath() = 0;
  virtual Float energy() = 0;
  virtual void evolve(Float dt, int steps) = 0;

  void init();

};

//!< Class describing bilinear action with rational matrix function
class AlgActionRational : public AlgActionBilinear {

 private:
  char *cname;
  int **fractionSplit;
  ActionRationalArg *rat_arg;

 protected:

  //!< Has any evolution taken place?
  int evolved;

  //!< Has the heatbath been evaluated?
  int heatbathEval;

  //!< Has the energy been evaluated?
  int energyEval;

  //!< This is where the rational parameters are stored
  RemezArg *remez_arg_md;
  RemezArg *remez_arg_mc;

  Float h_init;

  CgArg ***frm_cg_arg_md;
  CgArg ***frm_cg_arg_mc;
  //!< Pointer to an array of structures containing solver parameters.
  /*!<
    These are the parameters corresponding to each of the dynamical fermion
    masses.      
  */

  Vector** frmn;
  //!< Array of vectors
  /*!< These will hold the solutions from the solves. */
  
  Vector** frmn_d;
  //!< Array of vectors
  /*!< These will hold the solutions from the solves multiplied by
    the D-slash operator. */ 
  
  Vector** frmn_tmp;
  //!< Used for Asqtad fraction splitting

  int total_size;
  //!< The sum of the rational approximation degrees used for the force

  int max_size;
  //!< The maximum degree of rational approximation used for the force
  
  Float *all_res;
  //!< An array holding the residues - used for asqtad force


  //!< Renormalise the mass (optimisation)
  void massRenormalise(Float *mass, Float *trueMass, int degree, 
		       Float *shift, MassRenormaliseDir direction);

  //!< Automatic generation of the rational approximation.
  void generateApprox();

 public:

  AlgActionRational();
  AlgActionRational(AlgMomentum &mom, ActionRationalArg &rat_arg);
  virtual ~AlgActionRational();

  void heatbath();
  Float energy();

  //!< evolve method evolves the integrator
  void evolve(Float dt, int steps);
  //!< evolve method relevant to split timescales
  void evolve(Float dt, int steps, int **fractionSplit);

  //!< Compare two approximations to avoid recalculation if possible
  static int compareApprox(RemezArg &, RemezArg &);

  void init();

};

//!< Derived from AlgActionRational - allows independent evolution of
//!< partial fractions.
class AlgActionRationalSplit : public AlgActionRational {

 private:
  char *cname;
  AlgActionRational *rat;
  ActionRationalSplitArg *rat_split_arg;
  int **fractionSplit;

 public:
  AlgActionRationalSplit(AlgActionRational &rat, 
			 ActionRationalSplitArg &rat_split_arg); 
			 
  virtual ~AlgActionRationalSplit();

  void heatbath();
  Float energy();

  //!< evolve method evolves the integrator
  void evolve(Float dt, int steps);
  
  int getNmass();
  Float getMass(int);

  void cost(CgStats*);
};

//!< Class describing bosonic bilinear action
class AlgActionBoson : public AlgActionBilinear {

 private:
  char *cname;

  ActionBosonArg *bsn_arg;
  CgArg **bsn_cg_arg;   //!< Pointer to an array of solver parameters.

 public:

  AlgActionBoson(AlgMomentum &mom, ActionBosonArg &boson_arg);
  virtual ~AlgActionBoson();
  
  void heatbath();

  Float energy();

  //!< evolve method evolves the integrator
  void evolve(Float dt, int steps);

};

//!< Class describing fermionic bilinear action
class AlgActionFermion : public AlgActionBilinear {

 private:
  char *cname;

  ActionFermionArg *frm_arg;

  CgArg **frm_cg_arg_md;   //!< Pointer to an array of solver parameters.
  CgArg **frm_cg_arg_mc;   //!< Pointer to an array of solver parameters.

  int evolved;
  Float h_init;

  //!< Stores the history of cg solutions - used by chronological inversion
  Vector ***v;
  Vector ***cg_sol_old;
  Vector *cg_sol;

  //!< Stores the orthogonalised vectors multiplied by MatPcDagMatPc
  //!< These currently live in AlgActionFermion for a future
  //!< optimisation (chronological preconditioner)
  Vector ***vm;

  int *chrono;

 public:

  AlgActionFermion(AlgMomentum &mom, ActionFermionArg &frm_arg);
  virtual ~AlgActionFermion();
  
  void heatbath();

  Float energy();

  //!< evolve method evolves the integrator
  void evolve(Float dt, int steps);

  void init();
};

//!< Class describing pure gauge action
class AlgActionGauge : public AlgAction {

 private:
  char *cname;
  ActionGaugeArg *gauge_arg;
  GclassType gluon;

 public:

  AlgActionGauge(AlgMomentum &mom, ActionGaugeArg &gauge_arg);
  virtual ~AlgActionGauge();

  void heatbath();

  Float energy();

  //!< evolve method evolves the integrator
  void evolve(Float dt, int steps);

  void cost(CgStats*);

  void init();

};

#endif

CPS_END_NAMESPACE
