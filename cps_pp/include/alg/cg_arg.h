#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of the CgArg structure..

  $Id: cg_arg.h,v 1.4 2004-08-18 11:57:36 zs Exp $
*/
//---------------------------------------------------------------------------

#ifndef INCLUDED_CG_ARG_H
#define INCLUDED_CG_ARG_H		//!< Prevent multiple inclusion.

CPS_END_NAMESPACE
#include <util/enum.h>
#include <util/vector.h>
CPS_START_NAMESPACE

//! A structure to hold the solver parameters.
struct CgArg {

  Float mass;			/*!<  The mass parameter. */

  int max_num_iter;		/*!<  The maximum number of solver
				 iterations to do. */

  Float stop_rsd;		/*!<  The target residual. */
  Float true_rsd;               /*!<  The true residual. */

  enum RitzMatType RitzMatOper; /*!< Which operator to determine eigenvalues
				  of, if any. */
};

#endif /* !INCLUDED_CG_ARG_H */

CPS_END_NAMESPACE
