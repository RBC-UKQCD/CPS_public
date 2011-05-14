#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_int_ab.C
//
// AlgIntAB is derived from AlgInt, it is the class from which all
// intgerators with two operators are derived from, e.g., leapfrog.
// 
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<math.h>
#include<alg/alg_hmd.h>
#include<util/lattice.h>
#include<util/vector.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/alg_int.h>
CPS_START_NAMESPACE


AlgIntAB::AlgIntAB(AlgInt &a, AlgInt &b,  IntABArg &a_arg) 
  : AlgInt()
{

  cname = "AlgIntAB";
  VRB.Func(cname,"AlgIngAB()");
  A = &a;
  B = &b;
  ab_arg = &a_arg;
  A_steps = ab_arg->A_steps;
  B_steps = ab_arg->B_steps;
  level = ab_arg->level;
}

AlgIntAB::~AlgIntAB() {

}

// Maybe use tmp in future but can ignore for now
void AlgIntAB::heatbath() {
  traj++;
  A->heatbath();
  B->heatbath();
}

Float AlgIntAB::energy() {
  return A->energy() + B->energy();
}

void AlgIntAB::cost(CgStats *cg_stats) {
  A->cost(cg_stats);
  B->cost(cg_stats);
}

void AlgIntAB::reverse() {
  A->reverse();
  B->reverse();
}

void AlgIntAB::init() {
  A->init();
  B->init();
}

AlgIntAB& AlgIntAB::Create(AlgInt &A, AlgInt &B, IntABArg &ab_arg) {
  VRB.Func("AlgIntAB","Create()");
  if (ab_arg.type == INT_LEAP) {
    AlgIntLeap *leap = new AlgIntLeap(A, B, ab_arg);
    return *leap;
  } else if (ab_arg.type == INT_OMELYAN) {
    AlgIntOmelyan *ome = new AlgIntOmelyan(A, B, ab_arg);
    return *ome;
  } else if (ab_arg.type == INT_CAMPOSTRINI) {
    AlgIntCampostrini *ome = new AlgIntCampostrini(A, B, ab_arg);
    return *ome;
  } else if (ab_arg.type == INT_OMELYAN_44) {
    AlgIntOmelyan44 *ome = new AlgIntOmelyan44(A, B, ab_arg);
    return *ome;
  } else if (ab_arg.type == INT_OMELYAN_45) {
    AlgIntOmelyan45 *ome = new AlgIntOmelyan45(A, B, ab_arg);
    return *ome;
  } else if (ab_arg.type == INT_SUM) {
    AlgIntSum *sum = new AlgIntSum(A, B, ab_arg);
    return *sum;
  } else if (ab_arg.type == INT_FORCE_GRAD_QPQPQ || 
	     ab_arg.type == INT_FORCE_GRAD_PQPQP) {
    AlgIntForceGrad *fg = new AlgIntForceGrad(A, B, ab_arg);
    return *fg;
  } else {
    ERR.General("AlgIntFactory","CreateAB()","Integrator type not defined.\n");
  }
}

void AlgIntAB::Destroy(AlgIntAB &ab) {

  delete &ab;

}

CPS_END_NAMESPACE
