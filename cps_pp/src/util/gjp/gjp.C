#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of GlobalJobParameter class methods.

  $Id: gjp.C,v 1.15 2004-08-18 13:10:01 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 13:10:01 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/gjp/gjp.C,v 1.15 2004-08-18 13:10:01 zs Exp $
//  $Id: gjp.C,v 1.15 2004-08-18 13:10:01 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: gjp.C,v $
//  $Revision: 1.15 $
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

// CPS_END_NAMESPACE
// #ifdef PARALLEL
// #include <comms/sysfunc.h>
// #else
// #include <time.h>
// #endif
// CPS_START_NAMESPACE

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

  // Set the number of sites of a single node
  //----------------------------------------------------------------
  x_node_sites = rda.x_node_sites;
  y_node_sites = rda.y_node_sites;
  z_node_sites = rda.z_node_sites;
  t_node_sites = rda.t_node_sites;
  s_node_sites = rda.s_node_sites;   
  if(s_node_sites == 0) s_node_sites=1;
  
  if (x_node_sites<=0 ||x_node_sites%2!=0)
      ERR.General(cname,fname,
	"Bad value %d for x_node_sites; must be divisible by 2\n", x_node_sites);
  if (y_node_sites<=0 ||y_node_sites%2!=0)
      ERR.General(cname,fname,
	"Bad value %d for y_node_sites; must be divisible by 2\n", y_node_sites);
  if (z_node_sites<=0 ||z_node_sites%2!=0)
      ERR.General(cname,fname,
	"Bad value %d for z_node_sites; must be divisible by 2\n", z_node_sites);
  if (t_node_sites<=0 ||t_node_sites%2!=0)
      ERR.General(cname,fname,
	"Bad value %d for t_node_sites; must be divisible by 2\n", t_node_sites);
  if (s_node_sites<=0 || (s_node_sites>1 && s_node_sites%2!=0))
      ERR.General(cname,fname,
	"Bad value %d for s_node_sites; must be 1 or divisible by 2\n", t_node_sites);

  // Set the number of nodes
  //----------------------------------------------------------------
  x_nodes = rda.x_nodes;
  y_nodes = rda.y_nodes;
  z_nodes = rda.z_nodes;
  t_nodes = rda.t_nodes;
  s_nodes = rda.s_nodes;

  // On the QCDSP, check that at least one number_of_nodes is equal to 1
  //----------------------------------------------------------------
#if TARGET == QCDSP
  if( x_nodes != 1 &&
      y_nodes != 1 &&
      z_nodes != 1 &&
      t_nodes != 1 &&
      s_nodes != 1 ) {
      ERR.General(cname,fname,
		  "The processor grid must be of dimension 1 in at least one direction, but your grid is (X, Y, Z, T, S) = (%d, %d, %d, %d, %d).\n",
		  x_nodes, y_nodes, z_nodes, t_nodes, s_nodes);
  }
#endif
  

  // Check and set s_axis
  //-----------------------------------------------------------------
#if TARGET == QCDSP
  if(s_nodes != 1) {
      s_axis = rda.s_axis;
      if(s_axis != SCU_X &&
	 s_axis != SCU_Y &&
	 s_axis != SCU_Z &&
	 s_axis != SCU_T )
	  ERR.General(cname, fname,
		      "Illegal direction %d for the parallel S axis. Set DoArg::s_axis to one of 0, 1, 2 or 3.\n",
		      s_axis);

      if(s_axis == SCU_X && x_nodes != 1)
	  ERR. General(cname, fname, "You cannot choose %d as the direction for the parallel S axis (DoArg::s_axis = %d) AND parallelise along the X axis too (DoArg::x_nodes = %d).\n", s_axis, s_axis, x_nodes);
      
      if(s_axis == SCU_Y && y_nodes != 1) 
	  ERR. General(cname, fname, "You cannot choose %d as the direction for the parallel S axis (DoArg::s_axis = %d) AND parallelise along the Y axis too (DoArg::y_nodes = %d).\n", s_axis,s_axis,  y_nodes);
      
      if(s_axis == SCU_Z && z_nodes != 1) 
	  ERR. General(cname, fname, "You cannot choose %d as the direction for the parallel S axis (DoArg::s_axis = %d) AND parallelise along the Z axis too (DoArg::z_nodes = %d).\n", s_axis, s_axis, t_nodes);
      
      if(s_axis == SCU_T && t_nodes != 1) 
	  ERR. General(cname, fname, "You cannot chosen %d as the direction for the parallel S axis (DoArg::s_axis = %d) AND parallelise along the T axis too (DoArg::t_nodes = %d).\n", s_axis, s_axis, t_nodes);

      
  }
#endif

  // Check that the number of nodes divides the machine into
  // same size partitions.
  //----------------------------------------------------------------
  int size_x = 1;
  int size_y = 1;
  int size_z = 1;
  int size_t = 1;
  int size_s = 1;
#ifdef PARALLEL
  if( x_nodes != 1) size_x = SizeX(); 
  if( y_nodes != 1) size_y = SizeY(); 
  if( z_nodes != 1) size_z = SizeZ(); 
  if( t_nodes != 1) size_t = SizeT(); 
  if( s_nodes != 1) size_s = SizeS(); 
#endif
#if TARGET==QCDSP
  if( s_nodes != 1) {
    if (s_axis == SCU_X) size_s = SizeX();
    if (s_axis == SCU_Y) size_s = SizeY();
    if (s_axis == SCU_Z) size_s = SizeZ();
    if (s_axis == SCU_T) size_s = SizeT();
  } 
#endif

  if( x_nodes == 0 || size_x%x_nodes != 0) 
      ERR.General(cname,fname,	
		  "Illegal machine partition in X direction; physical grid size = %d must be a multiple of DoArg::x_nodes = %d\n",
		  size_x, x_nodes );
  
  if( y_nodes == 0 || size_y%y_nodes != 0) 
      ERR.General(cname,fname,
		  "Illegal machine partition in Y direction; physical grid size = %d must be a multiple of DoArg::y_nodes = %d\n",
		  size_y, y_nodes);	
  
  if( z_nodes == 0 || size_z%z_nodes != 0) 
       ERR.General(cname,fname,
		  "Illegal machine partition in Z direction; physical grid size = %d must be a multiple of DoArg::z_nodes = %d\n",
		   size_z, z_nodes);

  if( t_nodes == 0 || size_t%t_nodes != 0) 
      ERR.General(cname,fname,
		  "Illegal machine partition in T direction; physical grid size = %d must be a multiple of DoArg::t_nodes = %d\n",
		  size_t, t_nodes);
  
  if( s_nodes == 0 || size_s%s_nodes != 0) 
      ERR.General(cname,fname,
		  "Illegal machine partition in S direction; physical grid size = %d must be a multiple of DoArg::s_nodes = %d\n",
		  size_s, s_nodes);
  

  // Set the volume values
  //----------------------------------------------------------------
  vol_node_sites = x_node_sites * 
                   y_node_sites * 
                   z_node_sites * 
                   t_node_sites; 
  vol_sites = x_nodes * 
              y_nodes * 
              z_nodes * 
              t_nodes *
              vol_node_sites;

  // Set the coordinates of the node
  //----------------------------------------------------------------
  x_node_coor = 0;
  y_node_coor = 0;
  z_node_coor = 0;
  t_node_coor = 0;
  s_node_coor = 0;
#ifdef PARALLEL
  if(x_nodes != 1) x_node_coor = CoorX() % x_nodes;
  if(y_nodes != 1) y_node_coor = CoorY() % y_nodes;
  if(z_nodes != 1) z_node_coor = CoorZ() % z_nodes;
  if(t_nodes != 1) t_node_coor = CoorT() % t_nodes;
#if TARGET == QCDSP
  if(s_nodes != 1) {
    if (s_axis == SCU_X) s_node_coor = CoorX() % s_nodes;
    if (s_axis == SCU_Y) s_node_coor = CoorY() % s_nodes;
    if (s_axis == SCU_Z) s_node_coor = CoorZ() % s_nodes;
    if (s_axis == SCU_T) s_node_coor = CoorT() % s_nodes;
  }
#else
  if(s_nodes != 1) s_node_coor = CoorS() % s_nodes;
#endif
#endif

  // Set the static arrays gjp_local_axis[5], gjp_scu_dir[10],
  // and gjp_scu_wire_map[10].
  //----------------------------------------------------------------
#ifdef PARALLEL

  for(int la=0; la<6; la++) gjp_local_axis[la] = 0;
  if(x_nodes == 1) gjp_local_axis[0] = gjp_local_axis[5] = 1;
  if(y_nodes == 1) gjp_local_axis[1] = gjp_local_axis[5] = 1;
  if(z_nodes == 1) gjp_local_axis[2] = gjp_local_axis[5] = 1;
  if(t_nodes == 1) gjp_local_axis[3] = gjp_local_axis[5] = 1;
  if(s_nodes == 1) gjp_local_axis[4] = 1;

  gjp_scu_dir[0] = SCU_XP;
  gjp_scu_dir[1] = SCU_XM;
  gjp_scu_dir[2] = SCU_YP;
  gjp_scu_dir[3] = SCU_YM;
  gjp_scu_dir[4] = SCU_ZP;
  gjp_scu_dir[5] = SCU_ZM;
  gjp_scu_dir[6] = SCU_TP;
  gjp_scu_dir[7] = SCU_TM;

  gjp_scu_wire_map[0] = SCURemap( SCU_XP );
  gjp_scu_wire_map[1] = SCURemap( SCU_XM );
  gjp_scu_wire_map[2] = SCURemap( SCU_YP );
  gjp_scu_wire_map[3] = SCURemap( SCU_YM );
  gjp_scu_wire_map[4] = SCURemap( SCU_ZP );
  gjp_scu_wire_map[5] = SCURemap( SCU_ZM );
  gjp_scu_wire_map[6] = SCURemap( SCU_TP );
  gjp_scu_wire_map[7] = SCURemap( SCU_TM );

#if TARGET==QCDSP 
  if(s_nodes != 1) {
    if (s_axis == SCU_X) {
      gjp_scu_dir[8] = SCU_XP;
      gjp_scu_dir[9] = SCU_XM;
      gjp_scu_wire_map[8] = SCURemap( SCU_XP );
      gjp_scu_wire_map[9] = SCURemap( SCU_XM );
    }
    if (s_axis == SCU_Y) {
      gjp_scu_dir[8] = SCU_YP;
      gjp_scu_dir[9] = SCU_YM;
      gjp_scu_wire_map[8] = SCURemap( SCU_YP );
      gjp_scu_wire_map[9] = SCURemap( SCU_YM );
    }
    if (s_axis == SCU_Z) {
      gjp_scu_dir[8] = SCU_ZP;
      gjp_scu_dir[9] = SCU_ZM;
      gjp_scu_wire_map[8] = SCURemap( SCU_ZP );
      gjp_scu_wire_map[9] = SCURemap( SCU_ZM );
    }
    if (s_axis == SCU_T) {
      gjp_scu_dir[8] = SCU_TP;
      gjp_scu_dir[9] = SCU_TM;
      gjp_scu_wire_map[8] = SCURemap( SCU_TP );
      gjp_scu_wire_map[9] = SCURemap( SCU_TM );
    }
  }
#else
      gjp_scu_dir[8] = SCU_SP;
      gjp_scu_dir[9] = SCU_SM;
      gjp_scu_wire_map[8] = SCURemap( SCU_SP );
      gjp_scu_wire_map[9] = SCURemap( SCU_SM );
#endif // QCDSP


#endif //PARALLEL
  
  // Set the boundary conditions for the whole lattice
  //----------------------------------------------------------------
  x_bc = rda.x_bc;
  y_bc = rda.y_bc;
  z_bc = rda.z_bc;
  t_bc = rda.t_bc;

  // Set the boundary conditions for the sub-lattice on this node
  //----------------------------------------------------------------
  x_node_bc = BND_CND_PRD;
  if(x_bc == BND_CND_APRD) 
    x_node_bc = ( x_node_coor == (x_nodes-1) ) ? BND_CND_APRD : BND_CND_PRD;
  y_node_bc = BND_CND_PRD;
  if(y_bc == BND_CND_APRD) 
    y_node_bc = ( y_node_coor == (y_nodes-1) ) ? BND_CND_APRD : BND_CND_PRD;
  z_node_bc = BND_CND_PRD;
  if(z_bc == BND_CND_APRD) 
    z_node_bc = ( z_node_coor == (z_nodes-1) ) ? BND_CND_APRD : BND_CND_PRD;
  t_node_bc = BND_CND_PRD;
  if(t_bc == BND_CND_APRD) 
    t_node_bc = ( t_node_coor == (t_nodes-1) ) ? BND_CND_APRD : BND_CND_PRD;

  // Set the initial configuration kind.
  //----------------------------------------------------------------
  start_conf_kind = rda.start_conf_kind;

  // Set the initial configuration load address
  //----------------------------------------------------------------
  if(start_conf_kind == START_CONF_MEM || 
     start_conf_kind == START_CONF_LOAD ||
     start_conf_kind == START_CONF_FILE)
    start_conf_load_addr = rda.start_conf_load_addr;
  else
    start_conf_load_addr = 0;

  if(start_conf_kind == START_CONF_FILE){
    if(strlen(rda.start_conf_filename)<1){
      ERR.General(cname,fname,
      "start_conf_filename is not set correctly");
    } else
    strcpy(start_conf_filename,rda.start_conf_filename);
  }
    

  // Set the initial seed type.
  //----------------------------------------------------------------
  start_seed_kind = rda.start_seed_kind;

  // Set the initial seed value
  //----------------------------------------------------------------

  if(start_seed_kind == START_SEED_INPUT ||
     start_seed_kind == START_SEED_INPUT_UNIFORM ||
     start_seed_kind == START_SEED_INPUT_NODE) 
      start_seed_value = rda.start_seed_value;
  

  // Set beta.
  //----------------------------------------------------------------
  beta = rda.beta;

  // Set c_1
  //----------------------------------------------------------------
  c_1 = rda.c_1;

  //Set u_0
  //----------------------------------------------------------------
  u_0 = rda.u0;

  // Set dwf_height.
  //----------------------------------------------------------------
  dwf_height = rda.dwf_height;

  // Set the inverse of the dwf 5th dir. lattice spacing 
  //----------------------------------------------------------------
  dwf_a5_inv = rda.dwf_a5_inv;


  //------------------------------------------------------------------
  // Added in by Ping for anisotropic lattices and clover improvement.
  //------------------------------------------------------------------

  // Set parameters for anisotropic lattices and clover improvement.
  // MUST BE AFTER THE SETTING OF BETA [which is re-adjusted for
  // anisotropic implementations]. 
  //----------------------------------------------------------------
  xi_bare = rda.xi_bare;
  xi_dir  = rda.xi_dir;
  xi_v    = rda.xi_v;
  xi_v_xi = rda.xi_v_xi;
  xi_gfix = rda.xi_gfix;  // for Landau gauge
  clover_coeff_xi = clover_coeff = rda.clover_coeff;
  if (xi_bare != 1.0) {
    clover_coeff_xi = rda.clover_coeff_xi;
    beta /= xi_bare;    
    xi_gfix *= xi_bare;
  }

  //------------------------------------------------------------------
  // Added in by Ping for global sum
  //------------------------------------------------------------------
  // The following two parameters are relevant to scu transfer frequency
  // and hardware global sum.

  gsum_fast_mode = rda.gsum_fast_mode;
  gsum_max_try = rda.gsum_max_try;


  // Set the power plaquette cutoff.
  //----------------------------------------------------------------
  power_plaq_cutoff = rda.power_plaq_cutoff;

  // Set the power plaquette exponent.
  //----------------------------------------------------------------
  power_plaq_exponent = rda.power_plaq_exponent;

  // Set the power rectangle cutoff.
  //----------------------------------------------------------------
  power_rect_cutoff = rda.power_rect_cutoff;

  // Set the power rectangle exponent.
  //----------------------------------------------------------------
  power_rect_exponent = rda.power_rect_exponent;


  // Set the asqtad improved  staggered fermion action parameters.

  asqtad_KS = rda.asqtad_KS;	  
  asqtad_naik = rda.asqtad_naik;	  
  asqtad_3staple = rda.asqtad_3staple; 
  asqtad_5staple = rda.asqtad_5staple; 
  asqtad_7staple = rda.asqtad_7staple; 
  asqtad_lepage = rda.asqtad_lepage;  
 
  //================================================================
  // Other initializations
  //================================================================
    
  //----------------------------------------------------------------
  // Copy into cram the relevant vector library segment 
  //----------------------------------------------------------------
  p2vVector();

}


//------------------------------------------------------------------
/*!
  Sets the type of global lattice boundary condition in the X direction
  and also adjusts the local X direction boundary conditions to match.
  \param bc The type of boundary condition.
*/

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









CPS_END_NAMESPACE
