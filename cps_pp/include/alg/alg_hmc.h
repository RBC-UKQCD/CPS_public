#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief  Definition of the AlgHmc class.

*/
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_HMC_H
#define INCLUDED_ALG_HMC_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/hmc_arg.h>
#include <alg/alg_meas.h>
#include <alg/cg_stats.h>
#include <util/checksum.h>
#include <alg/alg_int.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
//! A class implementing the Hybrid Monte Carlo algorithm.
/*!

  This is an abstract implementation of HMC, the exact type of hmc is
  determined by the integrator that is passed to it.

  The main method evolves a gauge field along a molecular dynamics
  trajectory followed by a metropolis accept/reject step using the
  integrator which is passed to it.

*/
//------------------------------------------------------------------

class AlgHmc
{
 private:
    char *cname;

    //!< Save the gauge field and rngs (for accept/reject and repro)
    int saveInitialState();
    //!< Save the final gauge field (for repro test)
    void saveFinalState();
    //!< Restore the initial gauge field and rngs
    void restoreInitialState();
    //!< Restore the final gauge field 
    void restoreFinalState();
    //!< Test reproducibilty
    int reproduceTest(int attempt);
    //!< Used to test reversiblity
    void reverseTest();
    //!< used to shift lattice/RNG for a stronger reproducibility testing
    void shiftStates(const int x,const int y,const int z, const int t);

 protected:

    //!< The size of the gauge field.
    int g_size;
    /*!<
      The size of the gauge field on the local lattice in terms of the
      total number of floating point numbers.
    */

    Matrix* gauge_field_init;
    //!< The final gauge field configuration

    Matrix* gauge_field_final;
    //!< The initial gauge field configuration

//    unsigned int **rng4d_init;
//    unsigned int **rng5d_init;
      LRGState lrg_state;
    //!< The initial random numbers

    int config_no;

    Float h_init;   //!< Initial Hamiltonian
    Float h_final;  //!< Final Hamiltonian
    Float delta_h;  //!< Final-Initial Hamiltonian
    Float h_delta;  //!< Final-Initial Hamiltonian (reversibility)

    Matrix* mom;
    //!< The (traceless antihermitian) conjugate momentum field.
    /*!<
      This is \e i times the actual conjugate momentum.
     */	

    HmcArg *hmc_arg;
    CommonArg *common_arg;

    AlgIntAB *integrator;
    //!< The integrator which defines the algorithm

    unsigned int checksum[2];
    //!< Store the checksums of the final lattice (used for repro test)

 public:

  AlgHmc(AlgIntAB &integrator, CommonArg &c_arg, HmcArg &arg);

  virtual ~AlgHmc();

  //!< Performs a single HMC trajectory
  Float run(void);
};

#endif

CPS_END_NAMESPACE
