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

AlgAction::AlgAction(AlgMomentum &momentum, ActionArg &action_arg) 
  : AlgHamiltonian()
{

  cname = "AlgAction(Matrix*)";
  char *fname = "AlgAction(AlgMomentum&,ActionArg&)";
  VRB.Func(cname,fname);
  mom = momentum.getMom();
  force_measure = action_arg.force_measure;
  force_label = action_arg.force_label;
}

AlgAction::AlgAction() : AlgHamiltonian()
{

}

AlgAction::~AlgAction()
{
  char *fname = "~AlgAction()";
}

//!< Dummy method
void AlgAction::reverse()
{
  
}

void AlgAction::printForce(Float Fdt, Float dt, char *label) {

  char *fname = "printForce(Float, Float, char*)";
  FILE *fp;

#if TARGET==cpsMPI
  using MPISCU::fprintf;
#endif

  // Print out monitor info
  //---------------------------------------------------------------
  if( (fp = Fopen("force.dat", "a")) == NULL ) {
    ERR.FileA(cname,fname, "force.dat");
  }
  Fprintf(fp,"%s %e (L2) dt = %f\n", label, (IFloat)Fdt, (IFloat)dt);
  Fclose(fp);

}

CPS_END_NAMESPACE
