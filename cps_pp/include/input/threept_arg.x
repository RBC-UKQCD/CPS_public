#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Definition of the ThreePtArg structure.

  $Id: threept_arg.x,v 1.2 2004/12/11 20:57:38 chulwoo Exp $ 
*/
/*---------------------------------------------------------------------------*/

#ifndef INCLUDED_3PT_ARG_H
#define INCLUDED_3PT_ARG_H          /*!< Prevent multiple inclusion.*/

CPS_END_NAMESPACE
#include <alg/cg_arg.h>
#include <util/vector.h>
CPS_START_NAMESPACE

/*! The structure holds parameters relevant to the three-point functions measurement.*/
/*!  \ingroup algargs */
struct ThreePtArg {
  /* ??? */

  CgArg cg;		/*!< Parameters for fermion matrix inversion to
			   compute the quark propagator. */
  int seed;		/*!< Seed for random source */

    int t_src;		/*!< Source timeslice for quarks in I graphs. */

    int t_Op;		/*!< The operator timeslice */

    int t_Op_2;		/*!< The 2nd operator timeslice. */

    int t_sink;		/*!< Sink timeslice for spectator quark in I graphs. */

    int num_masses;	/*!< The number of masses to do.*/

    Float mass[20];	/*!< The list of masses. */
    

};

#endif /* !INCLUDED_3PT_ARG_H */

CPS_END_NAMESPACE
