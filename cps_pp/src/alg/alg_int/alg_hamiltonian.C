#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_hamiltonian.C
//
// AlgHamiltonian is a class defining methods common to all elements
// of the Hamiltonian
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<alg/alg_hmd.h>
#include<util/lattice.h>
#include<util/vector.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/alg_int.h>
CPS_START_NAMESPACE

AlgHamiltonian::AlgHamiltonian() : AlgInt()
{

}

AlgHamiltonian::~AlgHamiltonian()
{

}

CPS_END_NAMESPACE
