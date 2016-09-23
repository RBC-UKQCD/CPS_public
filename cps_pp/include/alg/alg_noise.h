#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief  Definitions of the AlgNoise class.

  $Id: alg_noise.h,v 1.4 2004/09/02 16:56:35 zs Exp $
*/
//------------------------------------------------------------------

#ifndef INCLUDED_ALG_NOISE_H
#define INCLUDED_ALG_NOISE_H           //!< Prevent multiple inclusion.

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/noise_arg.h>
CPS_START_NAMESPACE


//! Class for noising up a gauge configuration.
/*!
  Each gauge field link is left multiplied by an SU(3) group element
  taken at random from some distribution with unit mean. The type and width
  of the distribution is specified by the user.

  \ingroup alg
 */
class AlgNoise : public Alg
{
 private:
    char *cname;

    NoiseArg *alg_noise_arg;
        // The argument structure for the noise algorithm

    Matrix Exponentiate_Matrix( Matrix, int );
    
 public:
    AlgNoise(Lattice & latt, CommonArg *c_arg, NoiseArg *arg);

    virtual ~AlgNoise();

    //! Run the algorithm.
    void run(void);
};

#endif







CPS_END_NAMESPACE
