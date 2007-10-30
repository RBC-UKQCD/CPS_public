/*! A structure holding the parameters relevant to the HMC algorithm.*/

class HmcArg {

    int steps_per_traj;		    /*!<
				      For the R algorithm, this is
				      the number of steps per trajectory. For
				      the HMC/phi algorithm this is
				      <em> one less than </em> the number of
				      steps per trajectory.
				    */
  
    Float step_size;		    /*!< The molecular dynamics step size. */
  
    MetropolisType metropolis;        /*!< Whether to do the metropolis
					accept/reject step.
					This is relevant to the Phi
					algorithm only. */ 
  
    ReunitarizeType reunitarize;      /*!< Whether to reunitarize at the end of
					the trajectory before the metropolis
					step. */

    //! Reversibility test? 1 = TRUE, 0 = FALSE
    ReverseTest reverse;

    //! Reproduce test? 1 = TRUE, 0 = FALSE
    ReproduceTest reproduce;

    //! How many attempts do we try to reproduce?
    int reproduce_attempt_limit;

    //! whether sloppy precision is used for MD
    bool wfm_md_sloppy;
};
