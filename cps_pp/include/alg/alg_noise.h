#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_noise.h
//
// Header file for the AlgNoise class.
//
// AlgNoise is derived from Alg and is relevant to the
// multiplication of the current configuration with group elements 
// that take values that fluctuate according to some random 
// distribution (specified in NoiseArg) away from the identity.
// The magnitude of the fluctuations is specified in NoiseArg.
// 
//------------------------------------------------------------------

#ifndef INCLUDED_ALG_NOISE_H
#define INCLUDED_ALG_NOISE_H

CPS_END_NAMESPACE
#include<util/lattice.h>
#include<util/smalloc.h>
#include<util/pmalloc.h>
#include<alg/alg_base.h>
#include<alg/common_arg.h>
#include<alg/noise_arg.h>
CPS_START_NAMESPACE


class AlgNoise : public Alg
{
 private:
    char *cname;

    NoiseArg *alg_noise_arg;
        // The argument structure for the noise algorithm

 public:
    AlgNoise(Lattice & latt, CommonArg *c_arg, NoiseArg *arg);

    virtual ~AlgNoise();

    void run(void);
};

#endif






CPS_END_NAMESPACE
