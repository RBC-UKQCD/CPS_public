#ifndef INCLUDED_GLOBAL_JOB_PARAMETER_H
#define INCLUDED_GLOBAL_JOB_PARAMETER_H

#include<config.h>
/*!\file
  \brief  Definitions of global job parameters.

  $Id: gjp.h,v 1.47 2013-06-25 12:51:12 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-06-25 12:51:12 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/gjp.h,v 1.47 2013-06-25 12:51:12 chulwoo Exp $
//  $Id: gjp.h,v 1.47 2013-06-25 12:51:12 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: gjp.h,v $
//  $Revision: 1.47 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/gjp.h,v $
//  $State: Exp $
//--------------------------------------------------------------------
//
// gjp.h
//
// GlobalJobParameter is the base class. The constructor of this
// class sets the values of the global parameters. These values
// are accessible through function calls. An object of this class
// called GJP should be created at the highest scope (outside main).
// The header file declares GJP as external.
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


#include <util/lattice.h>
#include <util/vector.h>
#include <comms/sysfunc_cps.h>
#include <alg/do_arg.h>
#include <alg/cg_arg.h>

#ifdef USE_QUDA
#include <alg/quda_arg.h>
#endif

CPS_START_NAMESPACE

#ifdef PARALLEL
extern int gjp_local_axis[];
     // 1 if *_nodes = 1
     // 0 otherwise
     // Needed for fast access by communication routines
     // (some written in assembly).
     // It is set by GJP.Initialize.
     // {0,1,2,3,4} corresponds to {x,y,z,t,s}

//#if TARGET==QCDOC
#if 0
extern SCUDir gjp_scu_dir[];
     // set to:  SCU_XP, SCU_XM, SCU_YP, SCU_YM,
     // SCU_ZP, SCU_ZM, SCU_TP, SCU_TM, s_p, s_m
     // where s_p, s_m is one of the SCU_*P, SCU_*M.
     // Needed by get_plus_data, get_minus_data and glb_sum.
     // This in combination with gjp_local_axis determines
     // the direction for communication.
     // It is set by GJP.Initialize.
     // {0,1,2,3,4} corresponds to {x,y,z,t,s}
#endif

extern int gjp_scu_wire_map[];
     // it gives the wire number for directions
     // 0-9 corresponding to
     // x+, x-, y+, y-, z+, z-, t+, t-, s+, s-
     // The local wires are set to 0 but it is
     // assumed that gjp_local_axis is used in conjunction
     // so that the local direction wire number is not
     // used. 

#if TARGET == BGL
extern int bgl_machine_dir[8];
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

extern int bgl_cps_dir[8];
     // This array is set by GJP.Initialize to be the
     // "reverse" array of bgl_machine_dir. 
     // This array is for convenience when translating
     // from the the processor grid directions to the
     // physics system directions

#endif //TARGET == BGL

#endif //PARALLEL

//! A container class for global parameters.
/*! An object of this class, called GJP, should be created at the highest
  scope (outside main). The header file declares GJP as external. All global
  values are then accessible via the methods of this object.  
*/

//! Maximum filname for START_CONF_FILE.
const int MAX_FILENAME_LEN = 100;

//! Global parameters.
class GlobalJobParameter
{
 private:
  const char *cname;    // Class name.
  DoArg doarg_int;
  DoArgExt doext_int;
  DoArgExt *doext_p;

  int* argc_int;
  char*** argv_int;
  int arg_set;


  int node_sites[5]; // sites of a single node along {X,Y,Z,T,S} direction
  int nodes[5];      // number of nodes along {X,Y,Z,T,S} direction
  BndCndType bc[5];       // sites of a single node along {X,Y,Z,T,S} direction
  BndCndType node_bc[5];  // sites of a single node along {X,Y,Z,T,S} direction
  int node_coor[5];  // sites of a single node along {X,Y,Z,T,S} direction

  int vol_node_sites;  // The number of sites (4-D) of a single node.
  int vol_sites;       // The number of sites (4-D) of the whole lattice

  void Initialize();


  MdwfArg *mdwf_arg;
  MdwfTuning *mdwf_tuning;
  char *mdwf_tuning_fn;
  char *mdwf_tuning_record_fn;
public:
  GlobalJobParameter();

  ~GlobalJobParameter();

#if TARGET == BGL

  int BglMachineDirX(void)
    {return bgl_machine_dir[0]/2;}
     // bgl_machine_dir[0] = 2*bgl_machine_dir_x;
  //!< Gets the direction of the processor grid that X is "mapped" on. 
  /*!<
    \return The direction of the processor grid that X is "mapped" on. 
  */

  int BglMachineDirY(void)
    {return bgl_machine_dir[2]/2;}
  //!< Gets the direction of the processor grid that X is "mapped" on. 
  /*!<
    \return The direction of the processor grid that X is "mapped" on. 
  */

  int BglMachineDirZ(void)
    {return bgl_machine_dir[4]/2;}
  //!< Gets the direction of the processor grid that X is "mapped" on. 
  /*!<
    \return The direction of the processor grid that X is "mapped" on. 
  */

  int BglMachineDirT(void)
    {return bgl_machine_dir[6]/2;}
  //!< Gets the direction of the processor grid that X is "mapped" on. 
  /*!<
    \return The direction of the processor grid that X is "mapped" on. 
  */

#endif



  /*!\defgroup gjp_get_methods Methods that return the value of a global variable
    @{ */

  int Sites(int dir) const { return node_sites[dir]*nodes[dir]; }
   
  int NodeSites(int dir) const { return node_sites[dir]; }
  //!< Gets the dimension of the local lattice in a given direction.
  /*!<
    \param dir The direction in which to obtain the local lattice
    size; 0, 1, 2, 3 or 4 corresponding to X, Y, Z, T or S (the latter
    is only relevant for Domain Wall Fermions).
    \return The size of the local lattice in direction \a dir.
  */
  
    //making use of the storage order of the elements in the class, ie
    //the fact that x_node_sites, y_node_sites, z_node_sites, t_node_sites, 
    //s_node_sites, are in the order, x, y, z, t, s, with nothing in between. 
    //should be very careful when changing the storage order of these 
    //variables.
    //
    //This is not a good practice in general, but useful in many cases, where 
    //we need to know properties in a direction, but the direction 
    //parameter is a variable. This tweak requires the least modification
    //of the original code. 

  // return positions from index number, assuming xyzt ordering
  int LocalIndex(int index, int pos[]){
    int rest=index;
    for(int i=0;i<4;i++){
      pos[i] = rest%node_sites[i]; rest = rest/node_sites[i];
    }
    return rest;
  }

  int XnodeSites() const
      {return node_sites[0];}
  //!< Gets the dimension of the local lattice in the X direction.
  /*!<
    \return The size of the local lattice in the X direction.
  */

  int YnodeSites() const
      {return node_sites[1];}
  //!< Gets the dimension of the local lattice in the Y direction.  
  /*!<
    \return The size of the local lattice in the Y direction.
  */

  int ZnodeSites() const
      {return node_sites[2];}
  //!< Gets the dimension of the local lattice in the Z direction.
  /*!<
    \return The size of the local lattice in the Z direction.
  */

  int TnodeSites() const
      {return node_sites[3];}
  //!< Gets the dimension of the local lattice in the T direction.
  /*!<
    \return The size of the local lattice in the T direction.
  */

  int SnodeSites() const
      {return node_sites[4];}
  //!< Gets the dimension of the local lattice in the 5th direction.
  /*!<
    This is only relevant for Domain Wall Fermions.
    \return The size of the local lattice in the 5th direction.
  */

  int Nodes(int dir) const { return nodes[dir];}
  //!<Gets the dimension of the processor grid in a given direction.
  /*!<
    \param dir The direction in which to obtain the node grid
    size; 0, 1, 2, 3 or 4 corresponding to X, Y, Z, T or S (the latter
    is only relevant for Domain Wall Fermions).
    \return The size of the grid in direction \a dir.
  */

     // refer to the comments of NodeSites(int dir) above

  int Xnodes() const
      {return nodes[0];}
  //!< Gets the dimension of the node grid in the X direction.  
  /*!<
    \return The size of the grid in the X direction.
  */

  int Ynodes() const
      {return nodes[1];}
  //!< Gets the dimension of the node grid in the Y direction.  
  /*!<
    \return The size of the grid in the Y direction.
  */

  int Znodes() const
      {return nodes[2];}
  //!< Gets the dimension of the node grid in the Z direction.  
  /*!<
    \return The size of the grid in the Z direction.
  */

  int Tnodes() const
      {return nodes[3];}
  //!< Gets the dimension of the node grid in the T direction.  
  /*!<
    \return The size of the grid in the T direction.
  */

  int Snodes() const
      {return nodes[4];}
  //!< Gets the dimension of the node grid in the 5th direction.  
  /*!<
    This is only relevant for Domain Wall Fermions.
    \return The size of the grid in the 5th direction.
  */


  int VolNodeSites() const
      {return vol_node_sites;}
  //!< Gets the local lattice volume.
  /*!<
    In a domain wall fermion context, where the lattice is 5-dimensional,
    this is the volume of the 4-dimensional slice perpendicular to the 5th
    direction.
    \return The number of lattice sites on the node.
   */

  int VolSites() const
      {return vol_sites;}
  //!< Gets the global lattice volume.
  /*!<
    In a domain wall fermion context, where the lattice is 5-dimensional,
    this is the volume of the 4-dimensional slice perpendicular to the 5th
    direction.
    \return The number of lattice sites in the entire lattice..
   */

  int NodeCoor(int dir) {return node_coor[dir];}
  //!< Gets the grid coordinate of this node in a given direction.
  /*!<
    \param dir The direction in which to obtain the node grid
    coordinate; 0, 1, 2, 3 or 4 corresponding to X, Y, Z, T or S (the latter
    is only relevant for Domain Wall Fermions).
    \return The grid coordinate of this node in direction \a dir.
  */    
  // refer to the comments of NodeSites(int dir) above

  int XnodeCoor() const
      {return node_coor[0];}
  //!< Gets this nodes X direction grid coordinate.
  /*!<
    \return The grid coordinate of this node in the X direction.
  */

  int YnodeCoor() const
      {return node_coor[1];}
  //!< Gets this nodes Y direction grid coordinate.
  /*!<
    \return The grid coordinate of this node in the Y direction.
  */

  int ZnodeCoor() const
      {return node_coor[2];}
  //!< Gets this nodes Z direction grid coordinate.
  /*!<
    \return The grid coordinate of this node in the Z direction.
  */

  int TnodeCoor() const
      {return node_coor[3];}
  //!< Gets this nodes T direction grid coordinate.
  /*!<
    \return The grid coordinate of this node in the T direction.
  */

  int SnodeCoor()
      {return node_coor[4];}
  //!< Gets this nodes 5th direction grid coordinate.
  /*!<
    This is only relevant for Domain Wall Fermions.    
    \return The grid coordinate of this node in the 5th direction.
  */

  Float TwistBc(int dir) const
  { 
    switch(dir){
    case 0: return doext_p->twist_bc_x;
    case 1: return doext_p->twist_bc_y;
    case 2: return doext_p->twist_bc_z;
    case 3: return doext_p->twist_bc_t;
    default: printf("GJP::TwistBc(): Incorrect dir for twist\n"); 
      exit(0);
    }
    return 0.;
  }

  int Traj(){ return doext_p->trajectory; }

  BndCndType Bc(int dir) const
      { return bc[dir];}
  //!< Gets the global lattice boundary condition in a given direction.
  /*!< 
    \param dir The direction in which to obtain the boundary 
    condition; 0, 1, 2 or 3 corresponding to X, Y, Z or T.
    \return The type of boundary condition in direction \a dir.
  */
  // refer to the comments of NodeSites(int dir) above
  
  BndCndType Xbc() const
      {return bc[0];}
  //!< Gets the global lattice boundary condition in the X direction.
  /*!<
    \return The type of global boundary condition along the X axis.
  */

  BndCndType Ybc() const
      {return bc[1];}
  //!< Gets the global lattice boundary condition in the Y direction.
  /*!<
    \return The type of global boundary condition along the Y axis.
  */

  
  BndCndType Zbc() const
      {return bc[2];}
  //!< Gets the global lattice boundary condition in the Z direction.
  /*!<
    \return The type of global boundary condition along the Z axis.
  */

  BndCndType Tbc() const
      {return bc[3];}
  //!< Gets the global lattice boundary condition in the T direction.
  /*!<
    \return The type of global boundary condition along the T axis.
  */

  BndCndType NodeBc(int dir) const
      { return node_bc[dir];}
  //!< Gets the local lattice boundary condition in a given direction.
  /*!< 
    \param dir The direction in which to obtain the local boundary 
    condition; 0, 1, 2 or 3 corresponding to X, Y, Z or T.
    \return The type of boundary condition on this node in direction \a dir.
  */
  // refer to the comments of NodeSites(int dir) above

  BndCndType XnodeBc() const
      { return node_bc[0];}
  //!< Gets the local lattice boundary condition in the X direction.
  /*!<
    \return The type of local boundary condition along the X axis.
  */
    
  BndCndType YnodeBc() const
      { return node_bc[1];}
  //!< Gets the local lattice boundary condition in the Y direction.
  /*!<
    \return The type of local boundary condition along the Y axis.
  */
    
  BndCndType ZnodeBc() const
      { return node_bc[2];}
  //!< Gets the local lattice boundary condition in the Z direction.
  /*!<
    \return The type of local boundary condition along the Z axis.
  */
    
  BndCndType TnodeBc() const
      { return node_bc[3];}
  //!< Gets the local lattice boundary condition in the T direction.
  /*!<
    \return The type of local boundary condition along the T axis.
  */

  int CGreprodFreq() const
      {return doarg_int.cg_reprod_freq;}
  //!< Gets the frequency for CG reproducibility test

  StartConfType StartConfKind() const
      {return doarg_int.start_conf_kind;}
  //!< Gets the type of initial  gauge configuration.
  /*!<
    \return The type of initial gauge configuration.
  */    
  StartConfType StartU1ConfKind() const
      {return doext_p->start_u1_conf_kind;}
  //!< Gets the type of initial u1 gauge configuration.
  /*!<
    \return The type of initial u1 gauge configuration.
  */

  Matrix *StartConfLoadAddr() const
      {return (Matrix *)doarg_int.start_conf_load_addr;}
  void StartConfLoadAddr( Matrix * addr) 
      { doarg_int.start_conf_load_addr = (unsigned long) addr;}
  Float *StartU1ConfLoadAddr() const
      { if (!doext_p) return NULL;
        else return (Float *)doext_p->start_u1_conf_load_addr;}
  void StartU1ConfLoadAddr( Float * addr)
      { doext_p->start_u1_conf_load_addr = (unsigned long) addr;}
  //!< Gets the initial configuration.
  /*!<
    \return The address of the starting configuration
    if the gauge field starting type is \c START_CONF_MEM, 0 otherwise.
  */

  const char * StartConfFilename() const
      {return doarg_int.start_conf_filename;}
  const char * StartU1ConfFilename() const
      {return doext_p->start_u1_conf_filename;}

  const char * StartSeedFilename() const
      {return doarg_int.start_seed_filename;}

  const int StartConfAllocFlag() 
      {return doarg_int.start_conf_alloc_flag;}
  const int StartU1ConfAllocFlag()
      {return doext_p->start_u1_conf_alloc_flag;}
  const int mult_u1()
      {return doext_p->mult_u1_conf_flag;}

  const int WfmSendAllocFlag() 
      {return doarg_int.wfm_send_alloc_flag;}
  const int WfmAllocFlag() 
      {return doarg_int.wfm_alloc_flag;}

  StartSeedType StartSeedKind() const
      {return doarg_int.start_seed_kind;}
  //!< Gets the type of the initial RNG seed.
  /*!<
    \return The type of the initial RNG seed.
  */

  int StartSeedValue() const
      {return doarg_int.start_seed_value;}
  //!< Gets the value of the starting seed.
  /*!<
    \return The value of the starting seed.
   */
  
  int Colors() const {return 3;} 
  //!< Gets the number of colours.
  /*!<
    \return The number of colours.
  */

  int VerboseLevel() const   {return doarg_int.verbose_level;}


  Float Beta() const
      {return doarg_int.beta;}
  //!< Gets the "beta" parameter in the pure gauge action.
 
  int SaveStride() const
      {return doext_p->save_stride;}
  //!< Gets the stride (number) of eigenvectors to save at-a-time in Eigencontainer
 
  /*!
    \return The coefficient of the plaquette term in the pure gauge action. .
  */
  Float C1() const
      {return doarg_int.c_1;}
  //!< Gets c_1, the coefficient of the rectangle term in the pure gauge action.
  /*!<
- c_1 = 0 is the Wilson gauge action.
- c_1 = -0.05 is the tree level value.
- c_1 = -0.331 is the Iwasaki action.
    
    \return The coefficient of the rectangle term in the pure gauge action.
   */  

  Float u0() const
      {return doarg_int.u0;}
  //!< Gets the tadpole coefficient (the mean link).
  /*
    \return The tadpole coefficient.
  */


  Float DwfHeight() const
      {return doarg_int.dwf_height;}
  //!< Gets the height of the domain wall.
  /*!<
    Obviously, only relevant for Domain Wall Fermions.
    \return The height of the domain wall.
  */

  Float DwfA5Inv() const
      {return doarg_int.dwf_a5_inv;}
  //!< Gets the inverse of the 5th direction lattice spacing.
  /*!<
    Obviously, only relevant for Domain Wall Fermions.
    \return The inverse of the 5th direction lattice spacing.
  */
  Float Mobius_b() const
      {return doext_p->mobius_b_coeff;}
  Float Mobius_c() const
      {return doext_p->mobius_c_coeff;}
  

  //------------------------------------------------------------------
  // Added in by Ping for anisotropic lattices and clover improvement
  //------------------------------------------------------------------
  // The following parameters are relevant to anisotropic lattices and
  // clover improvement.
  // xi as a prefix  indicates relevancy to anisotropic lattices.
  // xi as a postfix indicates relevancy to the special anisotropic direction.

  Float XiBare() const            {return doarg_int.xi_bare;}
  //!< Gets the bare lattice anisotropy.
  /*!<
    The anisotropy is 1 for an isotropic lattice.
    \return The bare anisotropy,
   */

  int   XiDir()       const       {return doarg_int.xi_dir;}
  //!< Gets the anisotropic direction.
  /*!<
    This will be one of 0, 1, 2 or 3 corresponding to the X, Y, Z and
    T directions. The default is 3.
    \return The anisotropic direction.
  */
 
  Float XiV() const               {return doarg_int.xi_v;}
  //!< Gets the  bare speed of light.
  /*!<
    This is 1 for an isotropic lattice.
    \return The bare velocity of light.
  */
 
  Float XiVXi()    const          {return doarg_int.xi_v_xi;}
  //!<  Gets the bare speed of light in the anisotropic direction.
  /*
    This is 1 for an isotropic lattice.
    \return the bare speed of light  in the anisotropic direction.
  */
  
  Float CloverCoeff()   const     {return doarg_int.clover_coeff;}
  //!< Gets the clover coefficient.
  /*!<
    The coefficient of the clover term in the Sheikoleslami-Wohlert improved
    fermion action.
    This 1 for tree level improvement.
    \return The clover coefficient.
  */
 
  Float CloverCoeffXi() const     {return doarg_int.clover_coeff_xi;}
  //!< Gets the anisotropic clover coefficient.
  /*!<
    The coefficient of the clover term with links in the anisotropic
    direction Sheikoleslami-Wohlert improved anisotropic fermion action.
    \return The anisotropic clover coefficient.
  */
  void  XiV( Float xi_v ) {doarg_int.xi_v = xi_v;}
  void CloverCoeff(Float clover_coeff)
         {doarg_int.clover_coeff = clover_coeff;}
  void CloverCoeffXi(Float clover_coeff_xi)
         {doarg_int.clover_coeff_xi = clover_coeff_xi;}

  Float XiGfix() const            {return doarg_int.xi_gfix;}
  int GfixChkb() const 		  {return doarg_int.gfix_chkb;}
  //!< Gets the Landau gauge coefficient
  /*!<
    The coefficient for fixing to the Landau gauge on anisotropic
    lattices.  
    This is 1 for an isotropic lattice.
    \return The Landau gauge coefficient,
  */


  Float PowerPlaqCutoff() const
      {return doarg_int.power_plaq_cutoff;}
  //!< Gets the cut-off parameter of the power plaquette term in the pure gauge action.
  /*!< 
    \return The cut-off parameter.
  */

  int PowerPlaqExponent() const
      {return doarg_int.power_plaq_exponent;} 
  //!< Gets the exponent of the power plaquette term in the pure gauge action.
  /*!<
    \return The exponent.
  */


  Float PowerRectCutoff() const
      {return doarg_int.power_rect_cutoff;}
  //!< Gets the cut-off parameter of the power rectangle term in the pure gauge action.
  /*!<
    \return The cutoff parameter.
   */

  int PowerRectExponent() const
      {return doarg_int.power_rect_exponent;}
  //!< Gets the exponent of the power rectangle term in the pure gauge action.
  /*!<
    \return The exponent.
  */

  // Asqtad improved staggered action parameters

  //! Gets the coefficient of the Kogut-Susskind term in the Asqtad improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float KS_coeff() const { return doarg_int.asqtad_KS; }

  //! Gets the coefficient of the Naik term in the Asqtad improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float Naik_coeff() const { return doarg_int.asqtad_naik; }

  //! Gets the coefficient of the 3-staple term in the Asqtad improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float staple3_coeff() const { return doarg_int.asqtad_3staple; }

  //! Gets the coefficient of the 5-staple term in the Asqtad improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float staple5_coeff() const { return doarg_int.asqtad_5staple; }

  //! Gets the coefficient of the 7-staple term in the Asqtad improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float staple7_coeff() const { return doarg_int.asqtad_7staple; }
  
  //! Gets the coefficient of the Lepage term in the Asqtad improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float Lepage_coeff() const { return doarg_int.asqtad_lepage; }

  

  
   // P4 improved staggered action parameters

  //! Gets the coefficient of the Kogut-Susskind term in the P4 improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float p4_KS_coeff() const { return doarg_int.p4_KS; }

  //! Gets the coefficient of the Naik term in the P4 improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float p4_knight_coeff() const { return doarg_int.p4_knight; }

  //! Gets the coefficient of the 3-staple term in the P4 improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float p4_staple3_coeff() const { return doarg_int.p4_3staple; }

  //! Gets the coefficient of the 5-staple term in the P4 improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float p4_staple5_coeff() const { return doarg_int.p4_5staple; }

  //! Gets the coefficient of the 7-staple term in the P4 improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float p4_staple7_coeff() const { return doarg_int.p4_7staple; }
  
  //! Gets the coefficient of the Lepage term in the P4 improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float p4_Lepage_coeff() const { return doarg_int.p4_lepage; }

  
  /*! @} */

  /*!\defgroup gjp_set_methods Methods that set the value of a global variable
    @{ */

  //!  Initializes all global variables
  void Initialize(char *filename, char *instname);
  //! to be deprecated. Will point to the routine above
  void Initialize(const DoArg& do_arg);
  void InitializeExt(const DoArgExt& do_ext);
  int  ExtInitialized(){
    if (!doext_p) return 0;
    else return 1;
  }

  //PAB... Need to serialise the do arg as a means of meta-data preservation
  DoArg *GetDoArg(void) { return &doarg_int;};

  void SnodeSites(int sites)
      {node_sites[4] = sites;}
  //!< Sets the value of the dimension of the local lattice in the 5th direction.
  /*!<
    This is only relevant for Domain Wall Fermions.
    \param sites The dimension of the local lattice in the 5th direction.
  */


  //! Sets the global lattice boundary condition in the (dir) direction.
  void Bc(int dir, BndCndType cond);

  //! Sets the global lattice boundary condition in the X direction.
  void Xbc(BndCndType bc) { Bc(0,bc);}

  //! Sets the global lattice boundary condition in the Y direction.
  void Ybc(BndCndType bc) { Bc(1,bc);}

  //! Sets the global lattice boundary condition in the Z direction.
  void Zbc(BndCndType bc) { Bc(2,bc);}

  //! Sets the global lattice boundary condition in the T direction.
  void Tbc(BndCndType bc) { Bc(3,bc);}

  void StartConfKind(StartConfType sc)
      {doarg_int.start_conf_kind = sc;}
  //!< Sets the type of initial  gauge configuration.
  /*!<
    \param sc The type of initial gauge configuration.
  */
  void StartU1ConfKind(StartConfType sc)
      {doext_p->start_u1_conf_kind = sc;}
  //!< Sets the type of initial  gauge configuration.
  /*!<
    \param sc The type of initial gauge configuration.
  */

  void StartSeedKind(StartSeedType ss)
      {doarg_int.start_seed_kind = ss;}
  //!< Sets the type of the initial RNG seed.
  /*!<
    \param sc The type of the initial RNG seed.
  */

  void DwfHeight(Float height)
      {doarg_int.dwf_height = height;}
  //!< Sets the height of the domain wall.
  /*!<
    Obviously, only relevant for Domain Wall Fermions.
    \param height The height of the domain wall.
  */

  void DwfA5Inv(Float a5_inv)
      {doarg_int.dwf_a5_inv = a5_inv;}
  //!< Sets the inverse of the 5th direction lattice spacing.
  /*!<
    Obviously, only relevant for Domain Wall Fermions.
    \param a5_inv The inverse of the 5th direction lattice spacing.
  */

     // Sets the inverse of the dwf 5th dir. lattice spacing.


  // accessor for the main-arguments argc, argv

  int argc(void);
  char** argv(void);
  int *argc_p(void);
  char*** argv_p(void);
  
  void setArg(int* argc, char*** argv);


  void SetMdwfArg(const MdwfArg *);

  MdwfArg *GetMdwfArg(void){return mdwf_arg;}

  void FreeMdwfArg(void);

  bool InitMdwfTuning(const MdwfTuningInitArg &);

  MdwfTuning *GetMdwfTuning(void){return mdwf_tuning;}

  char *GetMdwfTuningFN(void){return mdwf_tuning_fn;}

  char *GetMdwfTuningRecordFN(void){return mdwf_tuning_record_fn;}

  /*! @} */
};

/*! An instance of the GlobalJobParameter class, named GJP, should be
  created at the highest scope (outside main). This external declaration
  allows access to all global variables.
*/
extern GlobalJobParameter GJP;


#ifdef USE_QUDA
  extern QudaArg QudaParam;
#endif

/*! declaration for Start() and End() which should be called at the
start and end of main()
*/

#if TARGET == NOARCH && (!defined USE_QMP)
inline void Start(){}
#endif
//inline void End(){}
//inline void Start(int * argc, char ***argv){GJP.setArg(argc, argv);}
void End();
void Start(int * argc, char ***argv);

#if TARGET == QCDOC 
extern "C" {
  void _mcleanup(void);
}
inline void Start(int * argc, char ***argv){Start(); GJP.setArg(argc, argv);}
#elif USE_QMP
namespace QMPSCU {
  void init_qmp();
  void init_qmp(int * argc, char *** argv);
  void destroy_qmp();
}
#elif TARGET == BGL
void Start(const BGLAxisMap *);
#endif

#ifdef USE_BFM
int cps_qdp_init(int *argc, char ***argv);
int cps_qdp_finalize();
#endif



CPS_END_NAMESPACE


#endif
