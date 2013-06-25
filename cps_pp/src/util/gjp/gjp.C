#include<config.h>
#include<math.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of GlobalJobParameter class methods.

  $Id: gjp.C,v 1.47 2013-06-25 12:51:12 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-06-25 12:51:12 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/gjp/gjp.C,v 1.47 2013-06-25 12:51:12 chulwoo Exp $
//  $Id: gjp.C,v 1.47 2013-06-25 12:51:12 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: gjp.C,v $
//  $Revision: 1.47 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/gjp/gjp.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// gjp.C
//
// GlobalJobParameter is the base class. The constructor of this
// class sets the values of the global parameters. These values
// are accessible through function calls. An object of this class
// called GJP should be created at the highest scope and it should
// be made global. The header file declares GJP as external.
//
//
// NOTE that the GJP.Xnodes, ... functions do not necessarily return 
// the same value as their qos sister functions. The GJP.Xnodes, ...
// return the values set by the do_arg structure. Because
// of this one can "divide" the machine into a number of identical
// hypercubes. Similarly the GJP.XnodeCoor, ... functions return
// the coordinate of the node in the divided section i.e.
// GJP.XnodeCoor = CoorX % GJP.Xnodes, where CoorX is the
// qos sytem function call.
// 
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdlib.h>
#include <util/qcdio.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/vector.h>
#include <util/error.h>
#include <util/checksum.h>
#include <alg/do_arg.h>
//#include <mem/p2v.h>

#ifdef USE_QUDA
#include <alg/quda_arg.h>
#endif

CPS_START_NAMESPACE

static const double SMALL = 1e-10;

#ifdef PARALLEL
int gjp_local_axis[6] = {0, 0, 0, 0, 1, 1}; 
     // For gjp_local_axis[n], n = {0,1,2,3,4}
     // corresponds to {x,y,z,t,s}. It is 1 if *_nodes = 1
     // and it is 0 otherwise.
     // gjp_local_axis[5] = 0 indicates that none of the
     // x,y,z,t directions are local. If = 1 it indicates that
     // at least one of the x,y,z,t directions is local.
     // Needed for fast access by communication routines
     // (some written in assembly).
     // It is set by GJP.Initialize.

#if 0
SCUDir gjp_scu_dir[10] = { SCU_XP, SCU_XM, SCU_YP, SCU_YM,	
                           SCU_ZP, SCU_ZM, SCU_TP, SCU_TM,
                           SCU_TP, SCU_TM };
     // set to:  SCU_XP, SCU_XM, SCU_YP, SCU_YM,
     // SCU_ZP, SCU_ZM, SCU_TP, SCU_TM, s_p, s_m
     // where s_p, s_m is one of the SCU_*P, SCU_*M.
     // Needed by get_plus_data, get_minus_data and glb_sum.
     // This in combination with gjp_local_axis determines
     // the direction for communication.
     // It is set by GJP.Initialize.
     // {0,1,2,3,4} corresponds to {x,y,z,t,s}
#endif

int gjp_scu_wire_map[10] = {0, 1, 2, 3, 4, 5, 6, 7, 0, 0};
     // it gives the wire number for directions
     // 0-9 corresponding to
     // x+, x-, y+, y-, z+, z-, t+, t-, s+, s-
     // The local wires are set to 0 but it is
     // assumed that gjp_local_axis is used in conjunction
     // so that the local direction wire number is not
     // used. 
#if TARGET == BGL
int bgl_machine_dir[8];
     // This array is set by GJP.Initialize to:
     // bgl_machine_dir[0] = 2*bgl_machine_dir_x;
     // bgl_machine_dir[1] = 2*bgl_machine_dir_x+1;
     // bgl_machine_dir[2] = 2*bgl_machine_dir_y;
     // bgl_machine_dir[3] = 2*bgl_machine_dir_y+1;
     // bgl_machine_dir[4] = 2*bgl_machine_dir_z;
     // bgl_machine_dir[5] = 2*bgl_machine_dir_z+1;
     // bgl_machine_dir[6] = 2*bgl_machine_dir_t;
     // bgl_machine_dir[7] = 2*bgl_machine_dir_t+1;
     // This array is for convenience when translating
     // from the physics system directions to the processor
     // grid directions.

int bgl_cps_dir[8];
     // This array is set by GJP.Initialize to be the
     // "reverse" array of bgl_machine_dir. 
     // This array is for convenience when translating
     // from the the processor grid directions to the
     // physics system directions
#endif //TARGET == BGL
#endif

GlobalJobParameter GJP;

#ifdef USE_QUDA
QudaArg QudaParam;
#endif

//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
GlobalJobParameter::GlobalJobParameter() 
{
  cname = "GlobalJobParameter";
  char *fname = "GlobalJobParameter()";
//  printf("%s::%s Entered\n",cname,fname);
//  VRB.Func(cname,fname);
  doext_p = NULL;
  arg_set=0;
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
GlobalJobParameter::~GlobalJobParameter() {
  char *fname = "~GlobalJobParameter()";
//  printf("%s::%s Entered\n",cname,fname);
  VRB.Func(cname,fname);
}

//------------------------------------------------------------------
/*!
  Initializes all global variables using the DoArg structure, and performs
  checks to make sure the values are suitable.
  \param rda Structure containing the initial values of the global variables
*/

void GlobalJobParameter::Initialize(char *filename, char *instname) {
  doarg_int.Decode(filename,instname);
  Initialize();
}

void GlobalJobParameter::Initialize(const DoArg& rda) {
#if 0
  ERR.General(cname,"Initizlize(DoArg)", 
    "No longer allowed. Please use GlobalJobParameter::Initialize(char *filename, char *instname) instead");
#else
  doarg_int = rda;
  Initialize();
#endif
}

void GlobalJobParameter::InitializeExt(const DoArgExt& rda) {
  doext_int = rda;
  doext_p = &doext_int;
//  Initialize();
}

void GlobalJobParameter::Initialize() {
  char *fname = "Initialize()";
  VRB.Func(cname,fname);
  int i, j;
  char *dim_name[5] = {"X","Y","Z","T","S"};

#if 0
  bgl_machine_dir[0] = 2*doarg_int.bgl_machine_dir_x;
  bgl_machine_dir[1] = 2*doarg_int.bgl_machine_dir_x+1;
  bgl_machine_dir[2] = 2*doarg_int.bgl_machine_dir_y;
  bgl_machine_dir[3] = 2*doarg_int.bgl_machine_dir_y+1;
  bgl_machine_dir[4] = 2*doarg_int.bgl_machine_dir_z;
  bgl_machine_dir[5] = 2*doarg_int.bgl_machine_dir_z+1;
  bgl_machine_dir[6] = 2*doarg_int.bgl_machine_dir_t;
  bgl_machine_dir[7] = 2*doarg_int.bgl_machine_dir_t+1;
  
  for(i = 0; i<8 ; i++)
    bgl_cps_dir[bgl_machine_dir[i]] = i;
  if (UniqueID()==0)
  for(i = 0; i<8 ; i++){
    printf("bgl_machine_dir[%d]=%d bgl_cps_dir[%d]=%d\n",
    i, bgl_machine_dir[i], i, bgl_cps_dir[i]);
  }
#endif



  // Set the number of nodes
  //----------------------------------------------------------------
  nodes[0] = doarg_int.x_nodes;
  nodes[1] = doarg_int.y_nodes;
  nodes[2] = doarg_int.z_nodes;
  nodes[3] = doarg_int.t_nodes;
  nodes[4] = doarg_int.s_nodes;

  if( nodes[0]==0 && nodes[1]==0 && nodes[2]==0 && nodes[3]==0 && nodes[4] == 0 )
  {
    nodes[0] = SizeX();
    nodes[1] = SizeY();
    nodes[2] = SizeZ();
    nodes[3] = SizeT();
    nodes[4] = SizeS();
  }
  VRB.Result(cname,fname, "nodes= %d %d %d %d %d\n",
nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]);

  // Check that the number of nodes divides the machine into
  // same size partitions.
  //----------------------------------------------------------------
  int size[5];
  size[0] = SizeX(); 
  size[1] = SizeY(); 
  size[2] = SizeZ(); 
  size[3] = SizeT(); 
  size[4] = SizeS(); 
//  VRB.Result(cname,fname, "size= %d %d %d %d %d\n",
//size[0], size[1], size[2], size[3], size[4]);
  for(i = 0; i<5 ; i++)
  if( nodes[i] == 0 || size[i]%nodes[i] != 0) 
      ERR.General(cname,fname,	
		  "Illegal machine partition in %s direction; physical grid size = %d must be a multiple of DoArg::%s_nodes = %d\n",
		  dim_name[i], size[i], dim_name[i], nodes[i] );
  

  // Set the number of sites of a single node
  //----------------------------------------------------------------
  node_sites[0] = doarg_int.x_node_sites;
  node_sites[1] = doarg_int.y_node_sites;
  node_sites[2] = doarg_int.z_node_sites;
  node_sites[3] = doarg_int.t_node_sites;
  node_sites[4] = doarg_int.s_node_sites;

  if( node_sites[0]==0 && node_sites[1]==0 && node_sites[2]==0 && node_sites[3]==0 && node_sites[4] == 0 )
    {
      node_sites[0] = doarg_int.x_sites/nodes[0];
      node_sites[1] = doarg_int.y_sites/nodes[1];
      node_sites[2] = doarg_int.z_sites/nodes[2];
      node_sites[3] = doarg_int.t_sites/nodes[3];
      node_sites[4] = doarg_int.s_sites/nodes[4];
    }

  for(i = 0; i<5 ; i++)
  if (node_sites[i] ==0 ) node_sites[i] = 1;
  for(i = 0; i<4 ; i++)
  if (node_sites[i]<=0 ||node_sites[i]%2!=0)
      ERR.General(cname,fname,
	"Bad value %d for %s_node_sites; must be divisible by 2\n", node_sites[i], dim_name[i]);

  // Set the volume values
  //----------------------------------------------------------------
 
  vol_node_sites = 1;
  for(i = 0; i<4 ; i++) vol_node_sites *= node_sites[i];
  vol_sites = vol_node_sites;
  for(i = 0; i<4 ; i++) vol_sites *= nodes[i];


  // Set the coordinates of the node
  //----------------------------------------------------------------

  for(i = 0; i<5 ; i++) node_coor[i] = 0;
#ifdef PARALLEL
  int coor[5];
  coor[0] = CoorX();
  coor[1] = CoorY();
  coor[2] = CoorZ();
  coor[3] = CoorT();
  coor[4] = CoorS();
	
  for(i = 0; i<5 ; i++){
    if(nodes[i] != 1) node_coor[i] = coor[i] % nodes[i];
    else node_coor[i]=0;
  }
  VRB.Result(cname,fname, "node_sites= %d %d %d %d %d\n",
node_sites[0], node_sites[1], node_sites[2], node_sites[3], node_sites[4]);
  VRB.Result(cname,fname, "coor= %d %d %d %d %d\n",
node_coor[0], node_coor[1], node_coor[2], node_coor[3], node_coor[4]);
//  printf("Node %d: coor= %d %d %d %d %d\n",UniqueID(),
//node_coor[0], node_coor[1], node_coor[2], node_coor[3], node_coor[4]);

  // Set the static arrays gjp_local_axis[5], gjp_scu_dir[10],
  // and gjp_scu_wire_map[10].
  //----------------------------------------------------------------

  for(int la=0; la<6; la++) gjp_local_axis[la] = 0;
  for(int la=0; la<4; la++){
  if(nodes[la] == 1) gjp_local_axis[la] = gjp_local_axis[5] = 1;
  }
  if(nodes[4] == 1) gjp_local_axis[4] = 1;
  if (!UniqueID())
  for(int la=0; la<4; la++){
    printf("dim %d: nodes=%d gjp_local_axis=%d\n",la,nodes[la],gjp_local_axis[la]);
  }

#if 0
  gjp_scu_dir[0] = SCU_XP;
  gjp_scu_dir[1] = SCU_XM;
  gjp_scu_dir[2] = SCU_YP;
  gjp_scu_dir[3] = SCU_YM;
  gjp_scu_dir[4] = SCU_ZP;
  gjp_scu_dir[5] = SCU_ZM;
  gjp_scu_dir[6] = SCU_TP;
  gjp_scu_dir[7] = SCU_TM;
  gjp_scu_dir[8] = SCU_SP;
  gjp_scu_dir[9] = SCU_SM;
#endif

#if TARGET == QCDOC
  for(int i = 0;i<5;i++)
  if(nodes[i] > 1){
  gjp_scu_wire_map[2*i]   = SCURemap(gjp_scu_dir[2*i]);
  gjp_scu_wire_map[2*i+1] = SCURemap(gjp_scu_dir[2*i+1]);
  }

#endif // TARGET == QCDOC
  
#endif //PARALLEL


  
  // Set the boundary conditions for the whole lattice
  //----------------------------------------------------------------
#ifndef UNIFORM_SEED_TESTING
  bc[0] = doarg_int.x_bc;
  bc[1] = doarg_int.y_bc;
  bc[2] = doarg_int.z_bc;
  bc[3] = doarg_int.t_bc;
#else
  for(i = 0; i<4 ; i++)
    bc[i] = BND_CND_PRD;
#endif

  // Set the boundary conditions for the sub-lattice on this node
  //----------------------------------------------------------------
  for(i = 0; i<4 ; i++){
  node_bc[i] = BND_CND_PRD;
  if(bc[i] == BND_CND_APRD) 
    node_bc[i] = ( node_coor[i] == (nodes[i]-1) ) ? BND_CND_APRD : BND_CND_PRD;
  }

  // Set the initial configuration load address
  //----------------------------------------------------------------
  StartConfType conf_kind = doarg_int.start_conf_kind;
   if(conf_kind != START_CONF_MEM && 
      conf_kind != START_CONF_LOAD )
//      conf_kind != START_CONF_FILE)
   doarg_int.start_conf_load_addr = 0;

    VRB.Result(cname,fname,"start_conf_alloc_flag=%d\n",doarg_int.start_conf_alloc_flag);
    
  // Set parameters for anisotropic lattices and clover improvement.
  // MUST BE AFTER THE SETTING OF BETA [which is re-adjusted for
  // anisotropic implementations]. 
  //----------------------------------------------------------------
  if (fabs(doarg_int.xi_bare - 1.0)> SMALL ) {
//    doarg_int.clover_coeff_xi = rda.clover_coeff_xi;
    doarg_int.beta /= doarg_int.xi_bare;    
    doarg_int.xi_gfix *= doarg_int.xi_bare;
  } else {
    doarg_int.clover_coeff_xi = doarg_int.clover_coeff;
  }

  //================================================================
  // Other initializations
  //================================================================

  VRB.DeactivateAll();
//  printf("verbose_level =%d\n",doarg_int.verbose_level);
  VRB.Level(doarg_int.verbose_level);
if (!UniqueID())
  printf("verbose_level =%d\n",VRB.Level());
  CSM.Initialize(1000);
  CSM.Activate(doarg_int.checksum_level);
if (!UniqueID())
  printf("checksum_level =%d\n",doarg_int.checksum_level);

 mdwf_arg = NULL;
 mdwf_tuning = NULL;
  VRB.FuncEnd(cname,fname);
}


//------------------------------------------------------------------
/*!
  Sets the type of global lattice boundary condition in the X direction
  and also adjusts the local X direction boundary conditions to match.
  \param bc The type of boundary condition.
*/


void GlobalJobParameter::Bc(int dir, BndCndType cond){

  // Set the x boundary condition for the whole lattice
  //----------------------------------------------------------------
  bc[dir] = cond;

  // Set the x boundary condition for the sub-lattice on this node
  //----------------------------------------------------------------
  node_bc[dir] = BND_CND_PRD;
  if(bc[dir] == BND_CND_APRD) 
    node_bc[dir] = ( node_coor[dir] == (nodes[dir]-1) ) ? BND_CND_APRD : BND_CND_PRD;
}


int GlobalJobParameter::argc(void){

  char *fname="argc()";
  if (!arg_set){
    ERR.General(cname,fname,"ERROR: GJP.argc not initialized\n");
    exit(-1);
  }
  return *argc_int;
}

int *GlobalJobParameter::argc_p(void){

  char *fname="argc_p()";
  if (!arg_set){
    ERR.General(cname,fname,"ERROR: GJP.argc not initialized\n");
    exit(-1);
  }
  return argc_int;
}


char** GlobalJobParameter::argv(void){

  char *fname="argv()";
  if (!arg_set){
    ERR.General(cname,fname,"ERROR: GJP.argv not initialized\n");
    exit(-1);
  }
  return *argv_int;
}


char*** GlobalJobParameter::argv_p(void){

  char *fname="argv_p()";
  if (!arg_set){
    ERR.General(cname,fname,"ERROR: GJP.argv not initialized\n");
    exit(-1);
  }
  return argv_int;
}

void GlobalJobParameter::setArg(int* argc, char*** argv){

  argc_int = argc;
  argv_int = argv;

  arg_set=1;
}

void GlobalJobParameter::SetMdwfArg(const MdwfArg *_mdwf_arg)
{
  const char *fname = "SetMdwfArg()";

  // free allocated memory.
  if(mdwf_arg != NULL) FreeMdwfArg();
  if(_mdwf_arg == NULL) return;

  mdwf_arg = (MdwfArg *)smalloc(cname, fname, "mdwf_arg", sizeof(MdwfArg));
  *mdwf_arg = *_mdwf_arg;

  int ls = mdwf_arg->b5.b5_len;
  int rsd_len = mdwf_arg->rsd_vec.rsd_vec_len;
  mdwf_arg->b5.b5_val = (Float *)smalloc(cname, fname, "b5_val", sizeof(Float)*ls);
  mdwf_arg->c5.c5_val = (Float *)smalloc(cname, fname, "c5_val", sizeof(Float)*ls);
  mdwf_arg->rsd_vec.rsd_vec_val = (Float *)smalloc(cname, fname, "rsd_vec_val", sizeof(Float)*rsd_len);

  int i;
  for(i = 0; i < ls; ++i){
    mdwf_arg->b5.b5_val[i] = _mdwf_arg->b5.b5_val[i];
    mdwf_arg->c5.c5_val[i] = _mdwf_arg->c5.c5_val[i];
  }
  for(i = 0; i < rsd_len; ++i){
    mdwf_arg->rsd_vec.rsd_vec_val[i] = _mdwf_arg->rsd_vec.rsd_vec_val[i];
  }
}

void GlobalJobParameter::FreeMdwfArg(void)
{
  const char *fname = "FreeMdwfArg()";
  if(mdwf_arg == NULL) return;

  sfree(cname, fname, "b5_val", mdwf_arg->b5.b5_val);
  sfree(cname, fname, "c5_val", mdwf_arg->c5.c5_val);
  sfree(cname, fname, "rsd_vec_val", mdwf_arg->rsd_vec.rsd_vec_val);
  sfree(cname, fname, "mdwf_arg", mdwf_arg);

  mdwf_arg = NULL;
}

static void InitMdwfArg(MdwfArg *mdwf_arg, int ls, int rc_max, int use_mdwf_for_dwf, int use_single_precision);

bool GlobalJobParameter::InitMdwfTuning(const MdwfTuningInitArg &mti_arg)
{
  const char *fname = "InitMdwfTuning(...)";

  static int initialized = 0;

  if(initialized != 0) return true;
  initialized = 1;

  int len1 = strlen(mti_arg.tuning_fn);
  int len2 = strlen(mti_arg.tuning_record_fn);

  mdwf_tuning_fn = (char *)smalloc(cname, fname, "mdwf_tuning_fn", (len1+1)*sizeof(char));
  mdwf_tuning_record_fn = (char *)smalloc(cname, fname, "mdwf_tuning_record_fn", (len2+1)*sizeof(char));

  strncpy(mdwf_tuning_fn, mti_arg.tuning_fn, len1);
  strncpy(mdwf_tuning_record_fn, mti_arg.tuning_record_fn, len2);

  // Check if the record file exists.
  // If the record file doesn't exist we create it and write some header information to the record file.
  if(UniqueID() == 0){
    FILE *fp = fopen(mdwf_tuning_record_fn, "r");

    if(fp == NULL){
      fp = fopen(mdwf_tuning_record_fn, "w");
      fprintf(fp, "#Mobius tuning test record. Automatically generated by tuning program, please do not modify this file.\n");
      fprintf(fp, "#%4s %17s %17s %17s %5s %17s %10s\n", "#ls", "b5", "c5", "residual", "rc", "time", "DWF CG");
    }
    fclose(fp);
  }

  mdwf_tuning = (MdwfTuning *)smalloc(cname, fname, "mwdf_tuning", sizeof(MdwfTuning));

  FILE *fp = fopen(mdwf_tuning_fn, "r");
  if(fp != NULL){
    fclose(fp);
    return mdwf_tuning->Decode(mdwf_tuning_fn, "mdwf_tuning");
  }else{
    int ls_min = mti_arg.ls_min;
    int ls_max = mti_arg.ls_max;

    // make sure we get even numbers,
    // Mdwf library doesn't like odd Ls.
    ls_min -= ls_min & 1;
    ls_max += ls_max & 1;

    mdwf_tuning->ls_min = ls_min;
    mdwf_tuning->ls_max = ls_max;
    mdwf_tuning->index = 0;
    mdwf_tuning->stage = 1;
    mdwf_tuning->rc_max = mti_arg.rc_max;
    mdwf_tuning->c5_range = mti_arg.c5_range;
    mdwf_tuning->rsd_granularity = mti_arg.rsd_granularity;

    // fill in 'poisoned' values
    mdwf_tuning->opti_time = -1.0;
    mdwf_tuning->opti_index = -2;
    mdwf_tuning->rsd_val = -1.0;
    mdwf_tuning->rsd_time = -1.0;
    mdwf_tuning->rc_val = -1;
    mdwf_tuning->rc_time = -1.0;

    int len = (ls_max - ls_min + 2)>>1;
    
    mdwf_tuning->mdwf_arg.mdwf_arg_len = len;
    mdwf_tuning->mdwf_arg.mdwf_arg_val = (MdwfArg *)smalloc(cname, fname, "mdwf_arg_val", sizeof(MdwfArg)*len);
    
    for(int i = 0; i < len; ++i){
      InitMdwfArg(&(mdwf_tuning->mdwf_arg.mdwf_arg_val[i]),
                  ls_min+2*i,
                  mdwf_tuning->rc_max,
                  mti_arg.use_mdwf_for_dwf,
                  mti_arg.use_single_precision);
    }
    
    mdwf_tuning->c5.c5_len = 4;
    mdwf_tuning->c5.c5_val = (C5State *)smalloc(cname, fname, "c5_val", sizeof(C5State)*4);
    memset(mdwf_tuning->c5.c5_val, 0, sizeof(C5State)*4);

    return true;
  }
}

static void InitMdwfArg(MdwfArg *mdwf_arg, int ls, int rc_max, int use_mdwf_for_dwf, int use_single_precision)
{
  const char *cname = "";
  const char *fname = "InitMdwfArg()";

  CgArg cg_arg_template;
  cg_arg_template.mass = -1.0;
  cg_arg_template.epsilon = -1.0;
  cg_arg_template.max_num_iter = -1;
  cg_arg_template.stop_rsd = -1.0;
  cg_arg_template.true_rsd = -1.0;
  cg_arg_template.RitzMatOper = NONE;
  cg_arg_template.Inverter = CG;
  cg_arg_template.bicgstab_n = -1;

  mdwf_arg->cg_arg = cg_arg_template;

  mdwf_arg->b5.b5_len = ls;
  mdwf_arg->c5.c5_len = ls;
  mdwf_arg->rsd_vec.rsd_vec_len = rc_max;

  mdwf_arg->b5.b5_val = (Float *)smalloc(cname, fname, "b5_val", sizeof(Float)*ls);
  mdwf_arg->c5.c5_val = (Float *)smalloc(cname, fname, "c5_val", sizeof(Float)*ls);
  mdwf_arg->rsd_vec.rsd_vec_val = (Float *)smalloc(cname, fname, "rsd_vec_val", sizeof(Float)*rc_max);

  mdwf_arg->use_single_precision = use_single_precision;
  mdwf_arg->use_mdwf_for_dwf = use_mdwf_for_dwf;
  mdwf_arg->M5 = GJP.DwfHeight();

  memset(mdwf_arg->b5.b5_val, 0, sizeof(Float)*ls);
  memset(mdwf_arg->c5.c5_val, 0, sizeof(Float)*ls);
  memset(mdwf_arg->rsd_vec.rsd_vec_val, 0, sizeof(Float)*rc_max);
}

CPS_END_NAMESPACE
