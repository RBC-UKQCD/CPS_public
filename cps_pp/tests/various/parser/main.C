#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /space/cvs/cps/cps++/tests/various/parser/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#define MAX_QUEFILE_LEN 10000 

CPS_END_NAMESPACE
#include<parser.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<task_func.h>
#include<alg/hmd_arg.h>
#include <strings.h>
CPS_START_NAMESPACE

GlobalJobParameter GJP;
LatRanGen LRG;
Verbose VRB;
Error ERR;



int main()
{
  char quefile[MAX_QUEFILE_LEN];
  strcpy(quefile, " ");

// PBP Gwilson Fstag
//------------------------------------------------------------------
//  strcat(quefile, " >PBPGWFS -mass 1.0E-05 -max_num_iter 500 -stop_rsd 1e-8 ");
//  strcat(quefile, " -src_u_s 0 -src_l_s 3 -snk_u_s 3 -snk_l_s 0 ");
//
// HMC Phi Gwilson Fstag
//------------------------------------------------------------------
  strcat(quefile, " >HMCPHIGWFS -n_frm_masses 2 -frm_mass0 0.1 -frm_mass1 0.2 ");
  strcat(quefile, " -frm_flavors0 2 -frm_flavors1 0 ");
  strcat(quefile, " -max_num_iter0 9 -max_num_iter1 15 ");
  strcat(quefile, " -n_bsn_masses 2 -bsn_mass0 0.8 -bsn_mass1 0.9 ");
  strcat(quefile, " -bsn_flavors0 2 -bsn_flavors1 0 ");
  strcat(quefile, " -max_num_iter0 9 -max_num_iter1 9 ");
  strcat(quefile, " -stop_rsd0 1e-10 -stop_rsd1 1e-10 ");
  strcat(quefile, " -steps_per_traj 4 -step_size 1e-02 ");
  strcat(quefile, " -metropolis 1 ");

//
// GHB Gwilson Fstag
//------------------------------------------------------------------
//  strcat(quefile, " >GHBGWFS -num_iter 500 "); 
//
// DO LOOP
//------------------------------------------------------------------
  strcat(quefile, " #DO -x_node_sites 2 -y_node_sites 2 ");
  strcat(quefile, " -z_node_sites 2 -t_node_sites 2 -s_node_sites 1 ");
  strcat(quefile, " -x_nodes 2 -y_nodes 2 -z_nodes 4 -t_nodes 4 ");
  strcat(quefile, " -x_bc 1 -y_bc 1 -z_bc 1 -t_bc 1 "); 
  strcat(quefile, " -start_conf_kind 1 -start_conf_load_addr 0 ");
  strcat(quefile, " -start_seed_kind 0 -start_seed_value 0 ");
  strcat(quefile, " -colors 3 -beta 2.5 -dwf_height 0.9 "); 
  strcat(quefile, " -verbose_level 100 -exec_task_list 2 ");

   Quefile q(quefile);
   
   Job job(q);

   VRB.Level(GJP.VerboseLevel());

   job.run();
   
   return 1;
}

CPS_END_NAMESPACE
