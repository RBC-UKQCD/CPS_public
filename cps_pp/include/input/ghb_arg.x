#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of the GhbArg structure.

  $Id: ghb_arg.x,v 1.2 2004-12-11 20:57:34 chulwoo Exp $
*/
/*-----------------------------------------------------------------------------*/

/*  The structure type GhbArg holds the parameters specific to
    the gauge field heat bath routines. */

#ifndef INCLUDED_GHB_ARG_H
#define INCLUDED_GHB_ARG_H          /*!< Prevent multiple inclusion.*/


/*!  The structure holds parameters specific to the gauge field heat bath routines. */
struct GhbArg {
  int  num_iter;		/*!< The number of heat bath iterations to do.*/
};

#endif /* !INCLUDED_GHB_ARG_H */

CPS_END_NAMESPACE
