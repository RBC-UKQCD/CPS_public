#include<config.h>
CPS_START_NAMESPACE
/*  ghb_arg.h */

/*  The structure type GhbArg holds the parameters specific to
    the gauge field heat bath routines. */

#ifndef INCLUDED_GHB_ARG_H
#define INCLUDED_GHB_ARG_H



struct GhbArg {
  int  num_iter;		/*  The number of heat bath
				    iterations to do. */
};

#endif /* !INCLUDED_GHB_ARG_H */
CPS_END_NAMESPACE
