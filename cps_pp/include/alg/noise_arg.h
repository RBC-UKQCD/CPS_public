#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of the NoiseArg structure.

  $Id: noise_arg.h,v 1.3 2004-08-18 11:57:36 zs Exp $
*/
//---------------------------------------------------------------------------
/*  noise_arg.h */

/*  The structure type NoiseArg holds the parameters relevant
    to the generation of a "noise" configuration. */

#ifndef INCLUDED_NOISE_ARG_H
#define INCLUDED_NOISE_ARG_H          //!< Prevent multiple inclusion.

CPS_END_NAMESPACE
#include <util/data_types.h>
CPS_START_NAMESPACE


//! Types of noise distribution.

enum NoiseType { GAUSSIAN = 0, /*!< Gaussian distribution with zero mean.  */
                 FLAT	  = 1  /*!< Uniform distribution with zero mean.  */
};

//! A structure holding the parameters relevant to the noisy gauge field calculation.
/*!  \ingroup algargs */
struct NoiseArg {

    NoiseType noise_kind;  /*!< The type of distribution from which
			     to draw the random group elements. */

    Float size;            /* The noise magnitude:
			      For the (zero mean) gaussian distribution this is
			      the standard deviation.
			      For the (zero mean) uniform distribution this is
			      half of the width.
			   */

};

#endif /* !INCLUDED_NOISE_ARG_H */

CPS_END_NAMESPACE
