#include<config.h>
#include<math.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of GlobalJobParameter class methods.

  $Id: gjp.C,v 1.24 2004-12-21 19:45:15 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-12-21 19:45:15 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/gjp/gjp.C,v 1.24 2004-12-21 19:45:15 chulwoo Exp $
//  $Id: gjp.C,v 1.24 2004-12-21 19:45:15 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: gjp.C,v $
//  $Revision: 1.24 $
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
#include <alg/do_arg.h>
#include <mem/p2v.h>
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

int gjp_scu_wire_map[10] = {0, 1, 2, 3, 4, 5, 6, 7, 0, 0};
     // it gives the wire number for directions
     // 0-9 corresponding to
     // x+, x-, y+, y-, z+, z-, t+, t-, s+, s-
     // The local wires are set to 0 but it is
     // assumed that gjp_local_axis is used in conjunction
     // so that the local direction wire number is not
     // used. 
#endif

GlobalJobParameter GJP;
//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
GlobalJobParameter::GlobalJobParameter() 
{
  cname = "GlobalJobParameter";
  char *fname = "GlobalJobParameter()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
GlobalJobParameter::~GlobalJobParameter() {
  char *fname = "~GlobalJobParameter()";
  VRB.Func(cname,fname);
}

//------------------------------------------------------------------
/*!
  Initializes all global variables using the DoArg structure, and performs
  checks to make sure the values are suitable.
  \param rda Structure containing the initial values of the global variables
*/
void GlobalJobParameter::Initialize(const DoArg& rda) {
  char *fname = "Initialize()";
  VRB.Func(cname,fname);
  int i, j;
  char *dim_name[5] = {"X","Y","Z","T","S"};

  doarg_int = rda;

  // Set the number of sites of a single node
  //----------------------------------------------------------------
  node_sites[0] = doarg_int.x_node_sites;
  node_sites[1] = doarg_int.y_node_sites;
  node_sites[2] = doarg_int.z_node_sites;
  node_sites[3] = doarg_int.t_node_sites;
  node_sites[4] = doarg_int.s_node_sites;

  // Set the number of nodes
  //----------------------------------------------------------------
  nodes[0] = doarg_int.x_nodes;
  nodes[1] = doarg_int.y_nodes;
  nodes[2] = doarg_int.z_nodes;
  nodes[3] = doarg_int.t_nodes;
  nodes[4] = doarg_int.s_nodes;
  if (nodes[4] == 0) nodes[4] = 1;

  for(i = 0; i<4 ; i++)
  if (node_sites[i]<=0 ||node_sites[i]%2!=0)
      ERR.General(cname,fname,
	"Bad value %d for %s_node_sites; must be divisible by 2\n", node_sites[i], dim_name[i]);

  // Check that the number of nodes divides the machine into
  // same size partitions.
  //----------------------------------------------------------------
  int size[5];
#if 0
  for(i = 0; i<5 ; i++)
  size[i] = 1;
#endif
  size[0] = SizeX(); 
  size[1] = SizeY(); 
  size[2] = SizeZ(); 
  size[3] = SizeT(); 
  size[4] = SizeS(); 
  for(i = 0; i<5 ; i++)
  if( nodes[i] == 0 || size[i]%nodes[i] != 0) 
      ERR.General(cname,fname,	
		  "Illegal machine partition in %s direction; physical grid size = %d must be a multiple of DoArg::x_nodes = %d\n",
		  dim_name[i], size[i], nodes[i] );
  

  if( node_sites[0]==0 && node_sites[1]==0 && node_sites[2]==0 && node_sites[3]==0 && node_sites[4] == 0 )
    {
      node_sites[0] = rda.x_sites/nodes[0];
      node_sites[1] = rda.y_sites/nodes[1];
      node_sites[2] = rda.z_sites/nodes[2];
      node_sites[3] = rda.t_sites/nodes[3];
      node_sites[4] = rda.s_sites/nodes[4];
    }

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
  }
  VRB.Result(cname,fname, "node_sites= %d %d %d %d %d\n",
node_sites[0], node_sites[1], node_sites[2], node_sites[3], node_sites[4]);
  VRB.Result(cname,fname, "nodes= %d %d %d %d %d\n",
nodes[0], nodes[1], nodes[2], nodes[3], nodes[4]);

  // Set the static arrays gjp_local_axis[5], gjp_scu_dir[10],
  // and gjp_scu_wire_map[10].
  //----------------------------------------------------------------

  for(int la=0; la<6; la++) gjp_local_axis[la] = 0;
  for(int la=0; la<4; la++){
  if(nodes[la] == 1) gjp_local_axis[la] = gjp_local_axis[5] = 1;
  }
  if(nodes[4] == 1) gjp_local_axis[4] = 1;

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

  for(int i = 0;i<5;i++)
  if(nodes[i] > 1){
  gjp_scu_wire_map[2*i]   = SCURemap(gjp_scu_dir[2*i]);
  gjp_scu_wire_map[2*i+1] = SCURemap(gjp_scu_dir[2*i+1]);
  }
  
#endif //PARALLEL


  
  // Set the boundary conditions for the whole lattice
  //----------------------------------------------------------------
  bc[0] = doarg_int.x_bc;
  bc[1] = doarg_int.y_bc;
  bc[2] = doarg_int.z_bc;
  bc[3] = doarg_int.t_bc;

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
      conf_kind != START_CONF_FILE)
   doarg_int.start_conf_load_addr = 0;

#if 0
  if(start_conf_kind == START_CONF_FILE){
    if(strlen(rda.start_conf_filename)<1){
      ERR.General(cname,fname,
      "start_conf_filename is not set correctly");
    } else
    strcpy(start_conf_filename,rda.start_conf_filename);
  }
#endif

    VRB.Flow(cname,fname,"start_conf_alloc_flag=%d\n",doarg_int.start_conf_alloc_flag);
    

#if 0
  xi_bare = rda.xi_bare;
  xi_dir  = rda.xi_dir;
  xi_v    = rda.xi_v;
  xi_v_xi = rda.xi_v_xi;
  xi_gfix = rda.xi_gfix;  // for Landau gauge
  clover_coeff_xi = clover_coeff = rda.clover_coeff;
#endif

  // Set parameters for anisotropic lattices and clover improvement.
  // MUST BE AFTER THE SETTING OF BETA [which is re-adjusted for
  // anisotropic implementations]. 
  //----------------------------------------------------------------
  if (fabs(doarg_int.xi_bare - 1.0)> SMALL ) {
    doarg_int.clover_coeff_xi = rda.clover_coeff_xi;
    doarg_int.beta /= doarg_int.xi_bare;    
    doarg_int.xi_gfix *= doarg_int.xi_bare;
  } else {
    doarg_int.clover_coeff_xi = doarg_int.clover_coeff;
  }

  //================================================================
  // Other initializations
  //================================================================
    

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

#if 0
void GlobalJobParameter::Xbc(BndCndType bc){ Bc(0,bc);}
void GlobalJobParameter::Ybc(BndCndType bc){ Bc(1,bc);}
void GlobalJobParameter::Zbc(BndCndType bc){ Bc(2,bc);}
void GlobalJobParameter::Tbc(BndCndType bc){ Bc(3,bc);}

void GlobalJobParameter::Xbc(BndCndType bc){

  // Set the x boundary condition for the whole lattice
  //----------------------------------------------------------------
  x_bc = bc;

  // Set the x boundary condition for the sub-lattice on this node
  //----------------------------------------------------------------
  x_node_bc = BND_CND_PRD;
  if(x_bc == BND_CND_APRD) 
    x_node_bc = ( x_node_coor == (x_nodes-1) ) ? BND_CND_APRD : BND_CND_PRD;
}


//------------------------------------------------------------------
/*!
  Sets the type of global lattice boundary condition in the Y direction
  and also adjusts the local Y direction boundary conditions to match.
  \param bc The type of boundary condition.
*/

void GlobalJobParameter::Ybc(BndCndType bc){

  // Set the y boundary condition for the whole lattice
  //----------------------------------------------------------------
  y_bc = bc;

  // Set the y boundary condition for the sub-lattice on this node
  //----------------------------------------------------------------
  y_node_bc = BND_CND_PRD;
  if(y_bc == BND_CND_APRD) 
    y_node_bc = ( y_node_coor == (y_nodes-1) ) ? BND_CND_APRD : BND_CND_PRD;
}


//------------------------------------------------------------------
/*!
  Sets the type of global lattice boundary condition in the Z direction
  and also adjusts the local Z direction boundary conditions to match.
  \param bc The type of boundary condition.
*/

void GlobalJobParameter::Zbc(BndCndType bc){

  // Set the z boundary condition for the whole lattice
  //----------------------------------------------------------------
  z_bc = bc;

  // Set the z boundary condition for the sub-lattice on this node
  //----------------------------------------------------------------
  z_node_bc = BND_CND_PRD;
  if(z_bc == BND_CND_APRD) 
    z_node_bc = ( z_node_coor == (z_nodes-1) ) ? BND_CND_APRD : BND_CND_PRD;
}


//------------------------------------------------------------------
/*!
  Sets the type of global lattice boundary condition in the T direction
  and also adjusts the local T direction boundary conditions to match.
  \param bc The type of boundary condition.
*/

void GlobalJobParameter::Tbc(BndCndType bc){

  // Set the t boundary condition for the whole lattice
  //----------------------------------------------------------------
  t_bc = bc;

  // Set the t boundary condition for the sub-lattice on this node
  //----------------------------------------------------------------
  t_node_bc = BND_CND_PRD;
  if(t_bc == BND_CND_APRD) 
    t_node_bc = ( t_node_coor == (t_nodes-1) ) ? BND_CND_APRD : BND_CND_PRD;
}
#endif


CPS_END_NAMESPACE
