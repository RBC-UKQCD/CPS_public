#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief  Definitions of the HmdArg structure.

  $Id: hmd_arg.h,v 1.11 2005-03-07 00:03:21 chulwoo Exp $
*/
//------------------------------------------------------------------


#ifndef INCLUDED_HMD_ARG_H
#define INCLUDED_HMD_ARG_H

CPS_END_NAMESPACE
#include <util/vector.h>
CPS_START_NAMESPACE


enum ReunitarizeType { REUNITARIZE_NO    = 0,
		       REUNITARIZE_YES   = 1 };

enum MetropolisType { METROPOLIS_NO    = 0,
		      METROPOLIS_YES   = 1 };



#define MAX_HMD_MASSES  8       //!< The maximum number of dynamical masses.

#define MAX_RAT_DEGREE  20          /* The maximum degree of the rational
				       approximation. */

//! A structure holding the parameters relevant to the HMD algorithms.
/*!
  Not all parameters are relevant to all algorithms. An algorithm should
  safely ignore  parameters irrelevant to it.
 \ingroup algargs 
*/

struct HmdArg {


    int n_frm_masses;                /*!< The number of dynamical fermion
				       masses. This is not necessarily the
				       number of flavours: At each mass the
				       number of flavours is given by
				       Lattice::ExactFlavors. */
				
    Float frm_mass[MAX_HMD_MASSES];   /*!< The array of dynamical fermions
					masses . */

    int frm_flavors[MAX_HMD_MASSES];  /*!<
					The flavour number of each dynamical
					mass:
					This is relevant to the R
					algorithm only. Bosonic fields
					are simulated by setting 
					a negative flavor number. */

    int n_bsn_masses;                 /*!< The number of dynamical boson 
					masses. This is not necessarily the
					number of flavours: At each mass the
					number of flavours is given by
					Lattice::ExactFlavors. This is relevant
					to the Phi algorithm only.*/
				
    Float bsn_mass[MAX_HMD_MASSES];   /* The array of dynamical bosons masses.
					 This is relevant to the Phi algorithm
					 only.*/

    int max_num_iter[MAX_HMD_MASSES]; /*!< The maximum number of iterations
					of the solver for each mass 
					of dynamical fermions/bosons */
                                   

    Float stop_rsd[MAX_HMD_MASSES];   /*!< The target residual for the solver
					stopping condition for each mass 
					of dynamical fermions/bosons */

    Float stop_rsd_md[MAX_HMD_MASSES];  /*!< The target residual for the solver
					   stopping condition for each mass 
					   of dynamical fermions/bosons */

    Float stop_rsd_mc[MAX_HMD_MASSES];  /*!< The target residual for the solver
					   stopping condition for each mass 
					   of dynamical fermions/bosons */

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
  
    //! Type of fields which are being simulated (FERMION OR BOSON)
    FieldType field_type[MAX_HMD_MASSES];

    //! Has a valid approximation been constructed?
    int valid_approx[MAX_HMD_MASSES];

    //! Parameters for the RHMC force rational approximations.
    int FRatDeg[MAX_HMD_MASSES];
    //! Parameters for the RHMC rational approximations.    
    Float FRatError[MAX_HMD_MASSES];
    //! Parameters for the RHMC rational approximations.        
    Float FRatNorm[MAX_HMD_MASSES];
    //! Parameters for the RHMC rational approximations.
    Float FRatRes[MAX_HMD_MASSES][MAX_RAT_DEGREE];
    //! Parameters for the RHMC rational approximations..
    Float FRatPole[MAX_HMD_MASSES][MAX_RAT_DEGREE];
    

    //! Parameters for the RHMC action rational approximations. 
    int SRatDeg[MAX_HMD_MASSES];
    //! Parameters for the RHMC action rational approximations.     
    Float SRatError[MAX_HMD_MASSES];
    //! Parameters for the RHMC action rational approximations.     
    Float SRatNorm[MAX_HMD_MASSES];
    //! Parameters for the RHMC action rational approximations.         
    Float SRatRes[MAX_HMD_MASSES][MAX_RAT_DEGREE];
    //! Parameters for the RHMC action rational approximations.             
    Float SRatPole[MAX_HMD_MASSES][MAX_RAT_DEGREE];

    //! Parameters for the RHMC inverse action rational approximations.
    Float SIRatNorm[MAX_HMD_MASSES];
    //! Parameters for the RHMC inverse action rational approximations.
    Float SIRatRes[MAX_HMD_MASSES][MAX_RAT_DEGREE];
    //! Parameters for the RHMC inverse action rational approximations.
    Float SIRatPole[MAX_HMD_MASSES][MAX_RAT_DEGREE];

    //! The approximation represented by the RHMC rational functions.
    int frm_power_num[MAX_HMD_MASSES];
    //! The approximation represented by the RHMC rational functions.
    int frm_power_den[MAX_HMD_MASSES];

    //! Whether the RHMC approximation is static or dynamic.
    int approx_type;

    //! The new degrees to use for dynamic RHMC.
    int FRatDegNew[MAX_HMD_MASSES];
    //! The new degrees to use for dynamic RHMC.
    int SRatDegNew[MAX_HMD_MASSES];
    
    Float lambda_low[MAX_HMD_MASSES];
    Float lambda_high[MAX_HMD_MASSES];
    Float lambda_min[MAX_HMD_MASSES];
    Float lambda_max[MAX_HMD_MASSES];

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

};


#endif /* INCLUDED_HMD_ARG_H */



CPS_END_NAMESPACE
