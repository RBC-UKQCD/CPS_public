#include<config.h>
CPS_START_NAMESPACE
/*  eq_state_arg.h */

/*  The structure type EqStateArg is the arguments structure
    for AlgEqState. */

#ifndef INCLUDED_EQ_STATE_ARG_H
#define INCLUDED_EQ_STATE_ARG_H

struct EqStateArg {
  int dir;      /* the special direction which determines what hyperplane(s)
		 * the sum of plaq is applied on. There are two sum
		 * operations, one on the hyperplane perpendicular to this
		 * direction, the other on all the hyperplanes parallel to
		 * this direction.
		 */
};

#endif /* !INCLUDED_EQ_STATE_ARG_H */

CPS_END_NAMESPACE
