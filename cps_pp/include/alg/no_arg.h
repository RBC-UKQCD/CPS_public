#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of a dummy (empty) argument structure.

  $Id: no_arg.h,v 1.4 2004-09-02 16:56:40 zs Exp $
*/
//---------------------------------------------------------------------------
/*  no_arg.h */


#ifndef INCLUDED_NO_ARG_H
#define INCLUDED_NO_ARG_H          //!< Prevent multiple inclusion.

//! An absence of parameters
/*!   Is an empty structure,
    needed for the cases where no arguments are necessary */
struct NoArg {

};

#endif /* !INCLUDED_NO_ARG_H */

CPS_END_NAMESPACE
