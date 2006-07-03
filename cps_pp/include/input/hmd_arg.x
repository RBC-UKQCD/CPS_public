/*! A structure holding the parameters relevant to the HMD algorithms.*/
/*!
  Not all parameters are relevant to all algorithms. An algorithm should
  safely ignore  parameters irrelevant to it.
 \ingroup algargs 
*/
typedef Float FRatVec[MAX_RAT_DEGREE];
typedef int IMassVec[MAX_HMD_MASSES];
typedef Float FMassVec[MAX_HMD_MASSES];
typedef FRatVec FRatMassVec[MAX_HMD_MASSES];

class HmdArg {


    int n_frm_masses;                /*!< The number of dynamical fermion
				       masses. This is not necessarily the
				       number of flavours: At each mass the
				       number of flavours is given by
				       Lattice::ExactFlavors. */
				
    int n_bsn_masses;                 /*!< The number of dynamical boson 
					masses. This is not necessarily the
					number of flavours: At each mass the
					number of flavours is given by
					Lattice::ExactFlavors. This is relevant
					to the Phi algorithm only.*/

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

    //! Whether the RHMC approximation is static or dynamic.
    RatApproxType approx_type;

    //! the allowed sprectral spread in the rational approximation
    Float spread;

    //! The precision used in the Remez algorithm of the RHMC approximation.
    long precision;

    //! The location of the smallest polar shift in the RHMC rational approximations.
    int isz;

    //! The Sexton-Weingarten factor
    /*! RHMC only.*/
    int sw;

    //! Chonological inversion parameter
    /*!
      The number of previous solutions used to form the initial solver
      guess (HMC only).
    */
    int chrono;

    //! Reproduce test? 1 = TRUE, 0 = FALSE
    int reproduce;

    //! How many attempts do we try to reproduce?
    int reproduce_attempt_limit;

    FMassVec frm_mass;   /*!< The array of dynamical fermions
					masses . */

    IMassVec frm_flavors;  /*!<
					The flavour number of each dynamical
					mass:
					This is relevant to the R
					algorithm only. Bosonic fields
					are simulated by setting 
					a negative flavor number. */

				
    FMassVec bsn_mass;   /* The array of dynamical bosons masses.
					 This is relevant to the Phi algorithm
					 only.*/

    IMassVec max_num_iter; /*!< The maximum number of iterations
					of the solver for each mass 
					of dynamical fermions/bosons */
                                   

    FMassVec stop_rsd;   /*!< The target residual for the solver
					stopping condition for each mass 
					of dynamical fermions/bosons */

    FMassVec stop_rsd_md;  /*!< The target residual for the solver
					   stopping condition for each mass 
					   of dynamical fermions/bosons */

    FMassVec stop_rsd_mc;  /*!< The target residual for the solver
					   stopping condition for each mass 
					   of dynamical fermions/bosons */

  
    //! Type of fields which are being simulated (FERMION OR BOSON)
    FieldType field_type[MAX_HMD_MASSES];

    //! Has a valid approximation been constructed?
    IMassVec valid_approx;

    FMassVec lambda_low;
    FMassVec lambda_high;
    FMassVec lambda_min;
    FMassVec lambda_max;

    /*Placeholder for splitting into two types*/
    RhmcPolesAction rhmc_poles_action;
    string rhmc_poles_file<>;

    //! The approximation represented by the RHMC rational functions.
    IMassVec frm_power_num;
    //! The approximation represented by the RHMC rational functions.
    IMassVec frm_power_den;

    //! Parameters for the RHMC force rational approximations.
    IMassVec FRatDeg;

    //! The new degrees to use for dynamic RHMC.
    IMassVec FRatDegNew;

    //! The new degrees to use for dynamic RHMC.
    IMassVec SRatDegNew;

    //! Parameters for the RHMC rational approximations.    
    FMassVec FRatError;
    //! Parameters for the RHMC rational approximations.        
    FMassVec FRatNorm;
    //! Parameters for the RHMC rational approximations.
    FRatMassVec FRatRes;
    //! Parameters for the RHMC rational approximations..
    FRatMassVec FRatPole;

    //! Parameters for the RHMC action rational approximations. 
    IMassVec SRatDeg;
    //! Parameters for the RHMC action rational approximations.     
    FMassVec SRatError;
    //! Parameters for the RHMC action rational approximations.     
    FMassVec SRatNorm;
    //! Parameters for the RHMC action rational approximations.         
    FRatMassVec SRatRes;
    //! Parameters for the RHMC action rational approximations.             
    FRatMassVec SRatPole;

    //! Parameters for the RHMC inverse action rational approximations.
    FMassVec SIRatNorm;
    //! Parameters for the RHMC inverse action rational approximations.
    FRatMassVec SIRatRes;
    //! Parameters for the RHMC inverse action rational approximations.
    FRatMassVec SIRatPole;

    

};
/*
 * PAB. While not an alg, it is convenient to place controls for evolution in
 * the VML portion of the HMD algorithm.
 */
class EvoArg {
  int traj_start;
  int gauge_unload_period;
  int gauge_configurations;
  int io_concurrency;
  int hdw_xcsum;
  int hdw_rcsum;
  int reproduce_interval;
  string ensemble_id<>;
  string ensemble_label<>;
  string creator<>;
  string gauge_file_stem<>;
  string rng_file_stem<>;
  string plaquette_stem<>;
  string pbp_stem<>;
  string evo_stem<>;
  string w_spect_directory<>;
  string work_directory<>;
  int measure_pbp;
//  int measure_eig;
  int measure_w_spect_interval;
//  string eig_lo_stem<>;
//  string eig_hi_stem<>;
};

class RhmcPolesState {

    //! The approximation represented by the RHMC rational functions.
    IMassVec frm_power_num;
    //! The approximation represented by the RHMC rational functions.
    IMassVec frm_power_den;

    //! Parameters for the RHMC force rational approximations.
    IMassVec FRatDeg;

    //! The new degrees to use for dynamic RHMC.
    IMassVec FRatDegNew;

    //! The new degrees to use for dynamic RHMC.
    IMassVec SRatDegNew;

    //! Parameters for the RHMC rational approximations.    
    FMassVec FRatError;
    //! Parameters for the RHMC rational approximations.        
    FMassVec FRatNorm;
    //! Parameters for the RHMC rational approximations.
    FRatMassVec FRatRes;
    //! Parameters for the RHMC rational approximations..
    FRatMassVec FRatPole;

    //! Parameters for the RHMC action rational approximations. 
    IMassVec SRatDeg;
    //! Parameters for the RHMC action rational approximations.     
    FMassVec SRatError;
    //! Parameters for the RHMC action rational approximations.     
    FMassVec SRatNorm;
    //! Parameters for the RHMC action rational approximations.         
    FRatMassVec SRatRes;
    //! Parameters for the RHMC action rational approximations.             
    FRatMassVec SRatPole;

    //! Parameters for the RHMC inverse action rational approximations.
    FMassVec SIRatNorm;
    //! Parameters for the RHMC inverse action rational approximations.
    FRatMassVec SIRatRes;
    //! Parameters for the RHMC inverse action rational approximations.
    FRatMassVec SIRatPole;

};    
