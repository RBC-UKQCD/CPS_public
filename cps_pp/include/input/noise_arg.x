/*!\file
  \brief  Definition of the NoiseArg structure.

  $Id: noise_arg.x,v 1.3 2006/05/22 21:12:04 chulwoo Exp $
*/
/*---------------------------------------------------------------------------*/
/*  noise_arg.h */

/*  The structure type NoiseArg holds the parameters relevant
    to the generation of a "noise" configuration. */


/*! Types of noise distribution.*/

enum NoiseType { GAUSSIAN = 0, /*!< Gaussian distribution with zero mean.  */
                 FLAT	  = 1  /*!< Uniform distribution with zero mean.  */
};

/*! A structure holding the parameters relevant to the noisy gauge field calculation.*/
/*!  \ingroup algargs */
class NoiseArg {

    NoiseType noise_kind;  /*!< The type of distribution from which
			     to draw the random group elements. */

    Float size;            /* The noise magnitude:
			      For the (zero mean) gaussian distribution this is
			      the standard deviation.
			      For the (zero mean) uniform distribution this is
			      half of the width.
			   */

};

