#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of the OverRelaxArg structure.

  $Id: overrelax_arg.h,v 1.2 2004-09-21 18:07:14 chulwoo Exp $
*/
//-----------------------------------------------------------------------------

/*  The structure type OverrelaxArg holds the parameters specific to
    the gauge field heat bath routines. */

#ifndef INCLUDED_OVERRELAX_ARG_H
#define INCLUDED_OVERRELAX_ARG_H          //!< Prevent multiple inclusion.


//!  The structure holds parameters specific to the gauge field heat bath routines. 
struct OverRelaxArg {
  int  num_iter;		/*!< The number of overrelaxation iterations to do.*/
};

#endif /* !INCLUDED_OVERRELAX_ARG_H */

CPS_END_NAMESPACE
