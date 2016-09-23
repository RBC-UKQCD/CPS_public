/*!\file
  \brief Definition of EqStateArg structure.

  $Id: eq_state_arg.h,v 1.3 2004/09/02 17:00:08 zs Exp $
*/
	
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/09/02 17:00:08 $
//  $Header: /space/cvs/cps/cps++/include/alg/eq_state_arg.h,v 1.3 2004/09/02 17:00:08 zs Exp $
//  $Id: eq_state_arg.h,v 1.3 2004/09/02 17:00:08 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: eq_state_arg.h,v $
//  $Revision: 1.3 $
//  $Source: /space/cvs/cps/cps++/include/alg/eq_state_arg.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

#include<config.h>
CPS_START_NAMESPACE

#ifndef INCLUDED_EQ_STATE_ARG_H
#define INCLUDED_EQ_STATE_ARG_H        //!< Prevent multiple inclusion.


//! Container of parameters for AlgEqState. 

struct EqStateArg {
  int dir;      /*!< The direction  determining what the hyperplanes
		  in which the plaquette is computed: The plaquette is
		  computed in the hyperplanes perpendicular and parallel
		  to this direction.	
		 */
};

#endif /* !INCLUDED_EQ_STATE_ARG_H */

CPS_END_NAMESPACE
