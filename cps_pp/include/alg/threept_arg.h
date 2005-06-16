#include <config.h>
/*!\file
  \brief Definition of the ThreePtArg structure.

  $Id: threept_arg.h,v 1.5 2005-06-16 07:15:06 chulwoo Exp $ 
*/
//---------------------------------------------------------------------------

CPS_START_NAMESPACE
#ifndef INCLUDED_3PT_ARG_H
#define INCLUDED_3PT_ARG_H          //!< Prevent multiple inclusion.
CPS_END_NAMESPACE
#include <alg/cg_arg.h>
#include <util/vector.h>
#include <alg/qpropw_arg.h>
CPS_START_NAMESPACE

//! The structure holds parameters relevant to the three-point functions measurement.
/*!  \ingroup algargs */
struct ThreePtArg {

  CgArg cg;		/*!< Parameters for fermion matrix inversion to
				  compute the quark propagator. */

  RandomType rng;      /*!< Type of RNG for random source */
  int seed;		/*!< Seed for random source */

  int gauge_fix; /*!< Type of gauge fixing to apply */

  int t_src;	/*!< Source timeslice */
  int t_snk;	/*!< Sink timeslice */  

  int t_op;		/*!< The location of the random propagator */
  int width;    /*!< Width of random propagator slab (starting at t_op) */
  int num_hits; /*!< Number of random propagators to do */
  
  int num_light;	/*!< The number of light masses */
  Float l_mass[20];	/*!< The list of light masses */
  int num_heavy;    /*!< The number of heavy masses */
  Float h_mass[20]; /*!< The list of heavy masses */
 
};

#endif /* !INCLUDED_3PT_ARG_H */
CPS_END_NAMESPACE
