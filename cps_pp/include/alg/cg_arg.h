#include<config.h>
CPS_START_NAMESPACE
/*  cg_arg.h */

/*  The structure type CgArg holds the parameters specific to
    the subroutine cg(). */

#ifndef INCLUDED_CG_ARG_H
#define INCLUDED_CG_ARG_H

CPS_END_NAMESPACE
#include<util/enum.h>
#include<util/vector.h>
CPS_START_NAMESPACE

struct CgArg {
  /* ??? */

  Float mass;			/*  The mass to use in the conjugate
				gradient solution. */

  int max_num_iter;		/*  The maximum number of conjugate
				gradient iterations to do. */


  Float stop_rsd;		/*  The residual for the stopping 
                                condition. */

  enum RitzMatType RitzMatOper; // Which operator to determine eigenvalues of
                                // if any
};

#endif /* !INCLUDED_CG_ARG_H */
CPS_END_NAMESPACE
