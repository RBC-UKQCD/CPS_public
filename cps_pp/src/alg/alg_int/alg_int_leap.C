#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_int_leap.C
//
// AlgIntLeap is derived from AlgInt, it is an implementation
// of the leapfrog integrator using abstract operators, i.e.
//
// To construct a QPQ integrator, the update to the coordinate
// must be the first argument, the momentum the second.
//
// To construct a PQP integrator, the update to the momentum
// must be the second argument, the coordinate the second.
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
#include<util/checksum.h>
#include<util/time_cps.h>
CPS_START_NAMESPACE


AlgIntLeap::AlgIntLeap(AlgInt &a, AlgInt &b, IntABArg &ab_arg) 
  : AlgIntAB(a,b,ab_arg)
{
  int_type = INT_LEAP;
  A_calls = 2;
  B_calls = 1;
}

AlgIntLeap::~AlgIntLeap() {

}

void AlgIntLeap::evolve(Float dt, int steps) 
{
  char * fname = "evolve(Float,int)";

  step_cnt = 0;
  if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(step_cnt);
//  int veloc_interval = this->s_ckpt;
  int veloc_interval = 1;
  if(veloc_interval<1) veloc_interval=steps;

// if( getVer(cname,fname)  <2) 
  A->evolve(0.5*dt/(Float)A_steps, A_steps);

  for (int i=0; i<steps; i++) {
    if (level == TOP_LEVEL_INTEGRATOR && (i%veloc_interval==0))
#if 1
    checkpoint(cname,fname,i,steps); // i will be modified if checkpoint exist
#else
     {
      std::stringstream veloc_label, veloc_label2;
#ifdef HAVE_VELOC
       Float dtime = -dclock();
// for debugging. Needs to be fixed
//       int traj_num=1;
       veloc_label <<"Veloc_traj"<<traj_num;
       std::cout << "Veloc label : "<< (veloc_label.str()).c_str() << std::endl;  
//       VRB.Result(cname,fname,"Veloc label : %s\n",(veloc_label.str()).c_str());
       int veloc_v = VELOC_Restart_test((veloc_label.str()).c_str(),steps);
       VRB.Result(cname,fname,"Veloc label version: %s %d\n",(veloc_label.str().c_str()),veloc_v);
       if(veloc_v == VELOC_FAILURE  || (veloc_v < i ) ){ // no restart, checkpoint
        VELOC_Checkpoint_begin((veloc_label.str()).c_str(), i );
        VELOC_Checkpoint_mem();
        VELOC_Checkpoint_end(1 );
       veloc_label2 <<"Veloc_traj"<<traj_num<<"_"<<i;
       VRB.Result(cname,fname, "Veloc Checkpoint label2: %s \n",(veloc_label2.str()).c_str());
        LRG.Write((veloc_label2.str()).c_str());
        } else {
       assert(veloc_v <steps);
        i = veloc_v;
        VELOC_Restart_begin (veloc_label.str().c_str(), i );
        VELOC_Recover_mem();
        VELOC_Restart_end(1 );
       veloc_label2 <<"Veloc_traj"<<traj_num<<"_"<<i;
       VRB.Result(cname,fname, "Veloc Recover label2: %s .\n",(veloc_label2.str()).c_str());
        LRG.Read(veloc_label2.str().c_str());
       }
       dtime += dclock();
       print_flops(fname,"VeloC()",0,dtime);
//       exit(-43);
#endif
//	 CSM.SaveComment(++step_cnt);
    }
#endif

    B->evolve(dt/(Float)B_steps, B_steps);    
    if (i < steps-1) A->evolve(dt/(Float)A_steps, A_steps);
    else A->evolve(0.5*dt/(Float)A_steps, A_steps);
  }

  if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(++step_cnt);

}

CPS_END_NAMESPACE
