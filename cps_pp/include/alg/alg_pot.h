#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_pot.h
//
// Header file for the AlgPot class.
//
// AlgPot is derived from Alg and it measures the potential
// into different propagation directions
//
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_POT_H
#define INCLUDED_ALG_POT_H

CPS_END_NAMESPACE
#include <math.h>    // for exponential function
#include <alg/pot_arg.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <util/lattice.h> // for PathOrderedProduct
CPS_START_NAMESPACE


const Float le=log(exp(1.));

class AlgPot : public Alg
{
 private:
    char *cname;
    
    PotArg *alg_pot_arg; // The argument structure for the potential
    Float norm_fac;     // normalization factor for colour and local volume
    Float xiB2;         // bare anisotropy squared
    // for the purpose of the (anisotropic) heat bath all 
    // temporal links were multiplied by the bare anisotropy
    // in this manner the update could be done using
    // the isotropic code --> we have to convert the temporal links back
    // aniso_factor = xi0^{-2* (number of temporal links) }

    static Float power(Float x, Float y) {return exp(y * log(x)/le );}
    // define x^y through  x^y = exp [y* ln(x)] and ln(x) = log(x)/log(e)

 public:
    AlgPot(Lattice & latt, CommonArg *c_arg, PotArg *arg);

    virtual ~AlgPot();

    void run(void);
    
    // historical junk -- uncommented by manke for v4.1.0 
    // Float Wzx[100];  
    // Float Wzt[100];
    
};



#endif





CPS_END_NAMESPACE
