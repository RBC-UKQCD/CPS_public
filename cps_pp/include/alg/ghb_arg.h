#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of the GhbArg structure.

  $Id: ghb_arg.h,v 1.2 2003-07-24 16:53:53 zs Exp $
*/
//-----------------------------------------------------------------------------

/*  The structure type GhbArg holds the parameters specific to
    the gauge field heat bath routines. */

#ifndef INCLUDED_GHB_ARG_H
#define INCLUDED_GHB_ARG_H          //!< Prevent multiple inclusion.


//!  The structure holds parameters specific to the gauge field heat bath routines. 
struct GhbArg {
  int  num_iter;		/*!< The number of heat bath iterations to do.*/
};

#endif /* !INCLUDED_GHB_ARG_H */

CPS_END_NAMESPACE
