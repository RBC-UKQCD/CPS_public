#include<config.h>
CPS_START_NAMESPACE
/*  noise_arg.h */

/*  The structure type NoiseArg holds the parameters relevant
    to the generation of a "noise" configuration. */

#ifndef INCLUDED_NOISE_ARG_H
#define INCLUDED_NOISE_ARG_H

CPS_END_NAMESPACE
#include<util/data_types.h>
CPS_START_NAMESPACE



enum NoiseType { GAUSSIAN = 0,
                 FLAT	  = 1};


struct NoiseArg {

  NoiseType noise_kind;  /* The kind of noise. */

  Float size;            /* The noise magnitude */

};

#endif /* !INCLUDED_NOISE_ARG_H */
CPS_END_NAMESPACE
