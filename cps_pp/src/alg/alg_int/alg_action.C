#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_action.C
//
// AlgAction is a class defining methods common to all actions
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

AlgAction::AlgAction(AlgMomentum &momentum) : AlgHamiltonian()
{

  cname = "AlgAction(Matrix*)";
  mom = momentum.getMom();

}

AlgAction::AlgAction() : AlgHamiltonian()
{

}

AlgAction::~AlgAction()
{

}

//!< Dummy method
void AlgAction::reverse()
{
  
}

CPS_END_NAMESPACE
