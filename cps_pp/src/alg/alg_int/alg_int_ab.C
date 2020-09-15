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
#include<util/time_cps.h>
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
//  traj++;
  const char *fname="heatbath()";
  std::string  veloc_label = "Phi_traj" +std::to_string(traj_num);
//  veloc_label<<"Phi_traj"<<traj_num;
//  char veloc_tmp[256];
//  snprintf(veloc_tmp,256,"Phi_traj%d",traj_num);

//  char *veloc_p = veloc_label.c_str();
  A->heatbath();
  B->heatbath();
#ifdef HAVE_VELOC
  if (level == TOP_LEVEL_INTEGRATOR){
     int veloc_v = getVer(veloc_label.c_str());
     VRB.Result(cname,fname,"veloc_v=%d\n",veloc_v);
     Float dtime = -dclock();
     if (veloc_v<1){
	VELOC_Checkpoint_begin(veloc_label.c_str(),1);
//	VELOC_Checkpoint_selective(VELOC_CKPT_SOME,phi_veloc_all.data(),phi_veloc_all.size());
	VELOC_Checkpoint_mem();
	VELOC_Checkpoint_end(1);
     } else {
        VELOC_Restart_begin (veloc_label.c_str(),1);
//	VELOC_Recover_selective(VELOC_CKPT_SOME,phi_veloc_all.data(),phi_veloc_all.size());
	VELOC_Recover_mem();
	VELOC_Restart_end(1);
        exit(-42);
     }
//     VELOC_Finalize(0);
     dtime +=dclock();
     print_time(fname,"VeloC()",dtime);
  }
#endif
}

Float AlgIntAB::energy() {
  Float out = A->energy() + B->energy();
  return out;
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
