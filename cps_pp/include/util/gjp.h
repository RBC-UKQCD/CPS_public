#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definitions of global job parameters.

  $Id: gjp.h,v 1.4 2003-10-21 15:34:42 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-10-21 15:34:42 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/gjp.h,v 1.4 2003-10-21 15:34:42 zs Exp $
//  $Id: gjp.h,v 1.4 2003-10-21 15:34:42 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.3.4.2  2003/10/21 15:26:55  zs
//  Need <sysfunc.h> not <comms/sysfunc.h> for the QCD{SP,OC} build.
//
//  Revision 1.3.4.1  2003/09/24 12:49:57  mcneile
//  I have updated the include path.
//
//  Revision 1.3  2003/08/12 16:22:51  zs
//  Added Asqtad action parameters.
//
//  Revision 1.3  2001/07/03 17:01:02  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.2  2001/06/19 18:13:17  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: gjp.h,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/gjp.h,v $
//  $State: Exp $
//
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

#ifndef INCLUDED_GLOBAL_JOB_PARAMETER_H
#define INCLUDED_GLOBAL_JOB_PARAMETER_H     //!< Prevent multiple inclusion


CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/vector.h>
CPS_START_NAMESPACE
#ifdef PARALLEL
CPS_END_NAMESPACE
#include <sysfunc.h>
CPS_START_NAMESPACE
#endif

struct DoArg;

#ifdef PARALLEL
extern int gjp_local_axis[];
     // 1 if *_nodes = 1
     // 0 otherwise
     // Needed for fast access by communication routines
     // (some written in assembly).
     // It is set by GJP.Initialize.
     // {0,1,2,3,4} corresponds to {x,y,z,t,s}

extern SCUDir gjp_scu_dir[];
     // set to:  SCU_XP, SCU_XM, SCU_YP, SCU_YM,
     // SCU_ZP, SCU_ZM, SCU_TP, SCU_TM, s_p, s_m
     // where s_p, s_m is one of the SCU_*P, SCU_*M.
     // Needed by get_plus_data, get_minus_data and glb_sum.
     // This in combination with gjp_local_axis determines
     // the direction for communication.
     // It is set by GJP.Initialize.
     // {0,1,2,3,4} corresponds to {x,y,z,t,s}

extern int gjp_scu_wire_map[];
     // it gives the wire number for directions
     // 0-9 corresponding to
     // x+, x-, y+, y-, z+, z-, t+, t-, s+, s-
     // The local wires are set to 0 but it is
     // assumed that gjp_local_axis is used in conjunction
     // so that the local direction wire number is not
     // used. 

#endif

//! A container class for global parameters.
/*! An object of this class, called GJP, should be created at the highest
  scope (outside main). The header file declares GJP as external. All global
  values are then accessible via the methods of this object.  
*/


class GlobalJobParameter
{
 private:
  char *cname;    // Class name.

  // Be very careful when changing the order of the variables or
  // adding new variables in this class! Please refer to the comments 
  // of "int NodeSites(int dir)" for a more detailed explanation. 

  int x_node_sites;  // sites of a single node along the x direction
  int y_node_sites;  // sites of a single node along the y direction
  int z_node_sites;  // sites of a single node along the z direction
  int t_node_sites;  // sites of a single node along the t direction
  int s_node_sites;  // sites of a single node along the s direction
                     // (5th dir.), relevant to DWF only

  int x_nodes;      // Number of nodes along the x direction
  int y_nodes;      // Number of nodes along the y direction
  int z_nodes;      // Number of nodes along the z direction
  int t_nodes;      // Number of nodes along the t direction
  int s_nodes;      // Number of nodes along the s direction
                    // (5th dir.), relevant to DWF only

#ifdef PARALLEL
  SCUAxis s_axis;   // The machine axis on which the 5th direction
                    // axis is mapped. Relevant to DWF with s_nodes
                    // different than 1.
#endif

  int vol_node_sites;  // The number of sites (4-D) of a single node.
  int vol_sites;       // The number of sites (4-D) of the whole lattice
      
  int x_node_coor;  // "coordinate" of the node along the x direction
  int y_node_coor;  // "coordinate" of the node along the y direction
  int z_node_coor;  // "coordinate" of the node along the z direction
  int t_node_coor;  // "coordinate" of the node along the t direction
  int s_node_coor;  // "coordinate" of the node along the s direction
                    // (5th dir.), relevant to DWF only

  BndCndType x_bc;  // Boundary conditions of the whole lattice along x
  BndCndType y_bc;  // Boundary conditions of the whole lattice along y
  BndCndType z_bc;  // Boundary conditions of the whole lattice along z
  BndCndType t_bc;  // Boundary conditions of the whole lattice along t
  
  BndCndType x_node_bc;  // Boundary conditions on this node along x
  BndCndType y_node_bc;  // Boundary conditions on this node along y
  BndCndType z_node_bc;  // Boundary conditions on this node along z
  BndCndType t_node_bc;  // Boundary conditions on this node along t
  
  StartConfType start_conf_kind; // The kind of initial 
                                 // configuration

  Matrix *start_conf_load_addr;     // The address of
                       // the starting conf. to be used if 
                       // start_conf = START_CONF_MEM
		       // Otherwise it is 0.

  StartSeedType start_seed_kind;  // The kind of initial 
                                  // random generator seed

  int start_seed_value;   // The value of the starting seed
                          // if start_seed_kind
                          // is START_SEED_INPUT or 
                          // is START_SEED_INPUT_UNIFORM
                          // Derived from do_arg.start_seed_value.

  int colors;       // The number of colors.

  Float beta;       // The pure gauge action "beta"

  Float c_1 ;       // Related to the coefficient of the rectangle
                    // term in the pure gauge action.
                    // c_1 = 0 is the Wilson gauge action.
                    // c_1 = -0.05 is the tree level value.
                    // c_1 = -0.331 is the Iwasaki action.

  Float u_0;	    // The juvenile amphibian

  Float dwf_height; // The height of the domain wall

  Float dwf_a5_inv; // The inverse of the dwf 5th dir. lattice spacing 

  Float power_plaq_cutoff;
     // The cutoff parameter for the power
     // plaquete term in the pure gauge action
     // implemented by the GpowerPlaq class.

  int power_plaq_exponent; 
     // The exponent for the power
     // plaquete term in the pure gauge action
     // implemented by the GpowerPlaq class.

  int verbose_level;  
     // verbose_level = 0   : No output
     // verbose_level = 1   : ...
     // ???
     // verbose_level = 100 : Full output


  //------------------------------------------------------------------
  // Added by Ping
  //------------------------------------------------------------------
  // The following parameters are relevant to anisotropic lattices and
  // clover improvement.
  // xi as a prefix  indicates relevancy to anisotropic lattices.
  // xi as a postfix indicates relevancy to the special anisotropic direction.
  
  Float xi_bare;            // bare anisotropy
  int   xi_dir;             // the special anisotropic direction 
                            // [0,1,2,3] for [x,y,z,t]
  Float xi_v;               // bare velocity of light
  Float xi_v_xi;            // bare velocity of light in the special direction
  Float clover_coeff;       // The clover fermion coefficient, 1 at tree level
  Float clover_coeff_xi;    // clover term coefficient for plaquettes with
                            // links along the special anisotropic direction.
  
  // Added for Landau gauge fixing for anisotropic lattices
  Float xi_gfix;            // coefficient for Laudau gauge fixing

  //------------------------------------------------------------------
  // Added by Ping
  //------------------------------------------------------------------
  // The following two parameters are relevant to scu transfer frequency
  // and hardware global sum.
  int  gsum_fast_mode;      // 0[by default] for 25MHz. 50MHz otherwise.
  int  gsum_max_try;        // max num of tries of global sum, 2 by default

  Float power_rect_cutoff;
     // The cutoff parameter for the power
     // rectangle term in the pure gauge action
     // implemented by the GpowerRect class.

  int power_rect_exponent; 
     // The exponent for the power
     // rectangle term in the pure gauge action
     // implemented by the GpowerRect class.


  // Asqtad improved staggered action parameters
  
  Float asqtad_naik;
  Float asqtad_3staple;
  Float asqtad_5staple;
  Float asqtad_7staple;
  Float asqtad_lepage;
 
public:
  GlobalJobParameter();

  ~GlobalJobParameter();

  //! \todo These should be const methods, no?
  /*!\defgroup gjp_get_methods Methods that return the value of a global variable
    @{ */

  int NodeSites(int dir){ return (&x_node_sites)[dir]; }
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

  int XnodeSites(void)
      {return x_node_sites;}
  //!< Gets the dimension of the local lattice in the X direction.
  /*!<
    \return The size of the local lattice in the X direction.
  */

  int YnodeSites(void)
      {return y_node_sites;}
  //!< Gets the dimension of the local lattice in the Y direction.  
  /*!<
    \return The size of the local lattice in the Y direction.
  */

  int ZnodeSites(void)
      {return z_node_sites;}
  //!< Gets the dimension of the local lattice in the Z direction.
  /*!<
    \return The size of the local lattice in the Z direction.
  */

  int TnodeSites(void)
      {return t_node_sites;}
  //!< Gets the dimension of the local lattice in the T direction.
  /*!<
    \return The size of the local lattice in the T direction.
  */

  int SnodeSites(void)
      {return s_node_sites;}
  //!< Gets the dimension of the local lattice in the 5th direction.
  /*!<
    This is only relevant for Domain Wall Fermions.
    \return The size of the local lattice in the 5th direction.
  */

  int Nodes(int dir){ return (&x_nodes)[dir];}
  //!<Gets the dimension of the processor grid in a given direction.
  /*!<
    \param dir The direction in which to obtain the node grid
    size; 0, 1, 2, 3 or 4 corresponding to X, Y, Z, T or S (the latter
    is only relevant for Domain Wall Fermions).
    \return The size of the grid in direction \a dir.
  */

     // refer to the comments of NodeSites(int dir) above

  int Xnodes(void)
      {return x_nodes;}
  //!< Gets the dimension of the node grid in the X direction.  
  /*!<
    \return The size of the grid in the X direction.
  */

  int Ynodes(void)
      {return y_nodes;}
  //!< Gets the dimension of the node grid in the Y direction.  
  /*!<
    \return The size of the grid in the Y direction.
  */

  int Znodes(void)
      {return z_nodes;}
  //!< Gets the dimension of the node grid in the Z direction.  
  /*!<
    \return The size of the grid in the Z direction.
  */

  int Tnodes(void)
      {return t_nodes;}
  //!< Gets the dimension of the node grid in the T direction.  
  /*!<
    \return The size of the grid in the T direction.
  */

  int Snodes(void)
      {return s_nodes;}
  //!< Gets the dimension of the node grid in the 5th direction.  
  /*!<
    This is only relevant for Domain Wall Fermions.
    \return The size of the grid in the 5th direction.
  */

#ifdef PARALLEL
  SCUAxis Saxis(void)
  {return s_axis;}
     // Returns the value of s_axis =
     // the machine axis on which the 5th direction axis is  
     // mapped. Relevant to DWF with s_nodes different than 1.
#endif

  int VolNodeSites(void)
      {return vol_node_sites;}
  //!< Gets the local lattice volume.
  /*!<
    In a domain wall fermion context, where the lattice is 5-dimensional,
    this is the volume of the 4-dimensional slice perpendicular to the 5th
    direction.
    \return The number of lattice sites on the node.
   */

  int VolSites(void)
      {return vol_sites;}
  //!< Gets the global lattice volume.
  /*!<
    In a domain wall fermion context, where the lattice is 5-dimensional,
    this is the volume of the 4-dimensional slice perpendicular to the 5th
    direction.
    \return The number of lattice sites in the entire lattice..
   */

  int NodeCoor(int dir) {return (&x_node_coor)[dir];}
  //!< Gets the grid coordinate of this node in a given direction.
  /*!<
    \param dir The direction in which to obtain the node grid
    coordinate; 0, 1, 2, 3 or 4 corresponding to X, Y, Z, T or S (the latter
    is only relevant for Domain Wall Fermions).
    \return The grid coordinate of this node in direction \a dir.
  */    
  // refer to the comments of NodeSites(int dir) above

  int XnodeCoor(void)
      {return x_node_coor;}
  //!< Gets this nodes X direction grid coordinate.
  /*!<
    \return The grid coordinate of this node in the X direction.
  */

  int YnodeCoor(void)
      {return y_node_coor;}
  //!< Gets this nodes Y direction grid coordinate.
  /*!<
    \return The grid coordinate of this node in the Y direction.
  */

  int ZnodeCoor(void)
      {return z_node_coor;}
  //!< Gets this nodes Z direction grid coordinate.
  /*!<
    \return The grid coordinate of this node in the Z direction.
  */

  int TnodeCoor(void)
      {return t_node_coor;}
  //!< Gets this nodes T direction grid coordinate.
  /*!<
    \return The grid coordinate of this node in the T direction.
  */

  int SnodeCoor(void)
      {return s_node_coor;}
  //!< Gets this nodes 5th direction grid coordinate.
  /*!<
    This is only relevant for Domain Wall Fermions.    
    \return The grid coordinate of this node in the 5th direction.
  */

  BndCndType Bc(int dir){ return (&x_bc)[dir];}
  //!< Gets the global lattice boundary condition in a given direction.
  /*!< 
    \param dir The direction in which to obtain the boundary 
    condition; 0, 1, 2 or 3 corresponding to X, Y, Z or T.
    \return The type of boundary condition in direction \a dir.
  */
  // refer to the comments of NodeSites(int dir) above
  
  BndCndType Xbc(void)
      {return x_bc;}
  //!< Gets the global lattice boundary condition in the X direction.
  /*!<
    \return The type of global boundary condition along the X axis.
  */

  BndCndType Ybc(void)
      {return y_bc;}
  //!< Gets the global lattice boundary condition in the Y direction.
  /*!<
    \return The type of global boundary condition along the Y axis.
  */

  
  BndCndType Zbc(void)
      {return z_bc;}
  //!< Gets the global lattice boundary condition in the Z direction.
  /*!<
    \return The type of global boundary condition along the Z axis.
  */

  BndCndType Tbc(void)
      {return t_bc;}
  //!< Gets the global lattice boundary condition in the T direction.
  /*!<
    \return The type of global boundary condition along the T axis.
  */

  BndCndType NodeBc(int dir) { return (&x_node_bc)[dir];}
  //!< Gets the local lattice boundary condition in a given direction.
  /*!< 
    \param dir The direction in which to obtain the local boundary 
    condition; 0, 1, 2 or 3 corresponding to X, Y, Z or T.
    \return The type of boundary condition on this node in direction \a dir.
  */
  // refer to the comments of NodeSites(int dir) above

  BndCndType XnodeBc(void)
      {return x_node_bc;}
  //!< Gets the local lattice boundary condition in the X direction.
  /*!<
    \return The type of local boundary condition along the X axis.
  */
    
  BndCndType YnodeBc(void)
      {return y_node_bc;}
  //!< Gets the local lattice boundary condition in the Y direction.
  /*!<
    \return The type of local boundary condition along the Y axis.
  */
    
  BndCndType ZnodeBc(void)
      {return z_node_bc;}
  //!< Gets the local lattice boundary condition in the Z direction.
  /*!<
    \return The type of local boundary condition along the Z axis.
  */
    
  BndCndType TnodeBc(void)
      {return t_node_bc;}
  //!< Gets the local lattice boundary condition in the T direction.
  /*!<
    \return The type of local boundary condition along the T axis.
  */

  StartConfType StartConfKind(void)
      {return start_conf_kind;}
  //!< Gets the type of initial  gauge configuration.
  /*!<
    \return The type of initial gauge configuration.
  */    

  Matrix *StartConfLoadAddr(void)
      {return start_conf_load_addr;}
  //!< Gets the initial configuration.
  /*!<
    \return The address of the starting configuration
    if the gauge field starting type is \c START_CONF_MEM, 0 otherwise.
  */

  StartSeedType StartSeedKind(void)
      {return start_seed_kind;}
  //!< Gets the type of the initial RNG seed.
  /*!<
    \return The type of the initial RNG seed.
  */

  int StartSeedValue(void)
      {return start_seed_value;}
  //!< Gets the value of the starting seed.
  /*!<
    \return The value of the starting seed.
   */
  

  int Colors(void)
      {return colors;}
  //!< Gets the number of colours.
  /*!<
    \return The number of colours.
  */

  Float Beta(void)
      {return beta;}
  //!< Gets the "beta" parameter in the pure gauge action.
  
  /*!
    \return The coefficient of the plaquette term in the pure gauge action. .
  */
  Float C1(void)
      {return c_1;}
  //!< Gets c_1, the coefficient of the rectangle term in the pure gauge action.
  /*!<
- c_1 = 0 is the Wilson gauge action.
- c_1 = -0.05 is the tree level value.
- c_1 = -0.331 is the Iwasaki action.
    
    \return The coefficient of the rectangle term in the pure gauge action.
   */  

  Float u0(void)
      {return u_0;}
  //!< Gets the tadpole coefficient (the mean link).
  /*
    \return The tadpole coefficient.
  */


  Float DwfHeight(void)
      {return dwf_height;}
  //!< Gets the height of the domain wall.
  /*!<
    Obviously, only relevant for Domain Wall Fermions.
    \return The height of the domain wall.
  */

  Float DwfA5Inv(void)
      {return dwf_a5_inv;}
  //!< Gets the inverse of the 5th direction lattice spacing.
  /*!<
    Obviously, only relevant for Domain Wall Fermions.
    \return The inverse of the 5th direction lattice spacing.
  */
  

  //------------------------------------------------------------------
  // Added in by Ping for anisotropic lattices and clover improvement
  //------------------------------------------------------------------
  // The following parameters are relevant to anisotropic lattices and
  // clover improvement.
  // xi as a prefix  indicates relevancy to anisotropic lattices.
  // xi as a postfix indicates relevancy to the special anisotropic direction.

  Float XiBare()             {return xi_bare;}
  //!< Gets the bare lattice anisotropy.
  /*!<
    The anisotropy is 1 for an isotropic lattice.
    \return The bare anisotropy,
   */

  int   XiDir()              {return xi_dir;}
  //!< Gets the anisotropic direction.
  /*!<
    This will be one of 0, 1, 2 or 3 corresponding to the X, Y, Z and
    T directions. The default is 3.
    \return The anisotropic direction.
  */
 
  Float XiV()                {return xi_v;}
  //!< Gets the  bare speed of light.
  /*!<
    This is 1 for an isotropic lattice.
    \return The bare velocity of light.
  */
 
  Float XiVXi()              {return xi_v_xi;}
  //!<  Gets the bare speed of light in the anisotropic direction.
  /*
    This is 1 for an isotropic lattice.
    \return the bare speed of light  in the anisotropic direction.
  */
  
  Float CloverCoeff()        {return clover_coeff;}
  //!< Gets the clover coefficient.
  /*!<
    The coefficient of the clover term in the Sheikoleslami-Wohlert improved
    fermion action.
    This 1 for tree level improvement.
    \return The clover coefficient.
  */
 
  Float CloverCoeffXi()      {return clover_coeff_xi;}
  //!< Gets the anisotropic clover coefficient.
  /*!<
    The coefficient of the clover term with links in the anisotropic
    direction Sheikoleslami-Wohlert improved anisotropic fermion action.
    \return The anisotropic clover coefficient.
  */

  Float XiGfix()             {return xi_gfix;}
  //!< Gets the Landau gauge coefficient
  /*!<
    The coefficient for fixing to the Landau gauge on anisotropic
    lattices.  
    This is 1 for an isotropic lattice.
    \return The Landau gauge coefficient,
  */

  //------------------------------------------------------------------
  // Added in by Ping for global sum
  //------------------------------------------------------------------
  // The following two parameters are relevant to scu transfer frequency
  // and hardware global sum.

  int GsumFastMode()         {return gsum_fast_mode;}
  // Returns the selection of scu frequency. 
  // 0[by default] for 25MHz. 50MHz otherwise.

  int GsumMaxTry()           {return gsum_max_try;}
  // Returns max num of tries of global sum, 2 by default.

  Float PowerPlaqCutoff(void)
      {return power_plaq_cutoff;}
  //!< Gets the cut-off parameter of the power plaquette term in the pure gauge action.
  /*!< 
    \return The cut-off parameter.
  */

  int PowerPlaqExponent(void)
      {return power_plaq_exponent;}
  //!< Gets the exponent of the power plaquette term in the pure gauge action.
  /*!<
    \return The exponent.
  */

  int VerboseLevel(void)
      {return verbose_level;}
  //!< Gets the verbosity level.
  /*!<
    See the Verbose class for details.
    \return The verbosity level.
  */

  Float PowerRectCutoff(void)
      {return power_rect_cutoff;}
  //!< Gets the cut-off parameter of the power rectangle term in the pure gauge action.
  /*!<
    \return The cutoff parameter.
   */

  int PowerRectExponent(void)
      {return power_rect_exponent;}
  //!< Gets the exponent of the power rectangle term in the pure gauge action.
  /*!<
    \return The exponent.
  */

  // Asqtad improved staggered action parameters

  //! Gets the coefficient of the Naik term in the Asqtad improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float Naik_coeff(){ return asqtad_naik; }

  //! Gets the coefficient of the 3-staple term in the Asqtad improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float staple3_coeff(){ return asqtad_3staple; }

  //! Gets the coefficient of the 5-staple term in the Asqtad improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float staple5_coeff(){ return asqtad_5staple; }

  //! Gets the coefficient of the 7-staple term in the Asqtad improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float staple7_coeff(){ return asqtad_7staple; }
  
  //! Gets the coefficient of the Lepage term in the Asqtad improved staggered fermion action.
  /*!
    \return The coefficient.
  */
  Float Lepage_coeff(){ return asqtad_lepage; }

  
  /*! @} */
  
  
  /*!\defgroup gjp_set_methods Methods that set the value of a global variable
    @{ */

  //!  Initializes all global variables
  void Initialize(const DoArg& do_arg);

  void SnodeSites(int sites)
      {s_node_sites = sites;}
  //!< Sets the value of the dimension of the local lattice in the 5th direction.
  /*!<
    This is only relevant for Domain Wall Fermions.
    \param sites The dimension of the local lattice in the 5th direction.
  */

  //! Sets the global lattice boundary condition in the X direction.
  void Xbc(BndCndType bc);

  //! Sets the global lattice boundary condition in the Y direction.
  void Ybc(BndCndType bc);

  //! Sets the global lattice boundary condition in the Z direction.
  void Zbc(BndCndType bc);

  //! Sets the global lattice boundary condition in the T direction.
  void Tbc(BndCndType bc);

  void StartConfKind(StartConfType sc)
      {start_conf_kind = sc;}
  //!< Sets the type of initial  gauge configuration.
  /*!<
    \param sc The type of initial gauge configuration.
  */

  void StartSeedKind(StartSeedType ss)
      {start_seed_kind = ss;}
  //!< Sets the type of the initial RNG seed.
  /*!<
    \param sc The type of the initial RNG seed.
  */

  void DwfHeight(Float height)
      {dwf_height = height;}
  //!< Sets the height of the domain wall.
  /*!<
    Obviously, only relevant for Domain Wall Fermions.
    \param height The height of the domain wall.
  */

  void DwfA5Inv(Float a5_inv)
      {dwf_a5_inv = a5_inv;}
  //!< Sets the inverse of the 5th direction lattice spacing.
  /*!<
    Obviously, only relevant for Domain Wall Fermions.
    \param a5_inv The inverse of the 5th direction lattice spacing.
  */

     // Sets the inverse of the dwf 5th dir. lattice spacing.


  //------------------------------------------------------------------
  // Added in by Ping for global sum
  //------------------------------------------------------------------
  // The following two parameters are relevant to scu transfer frequency
  // and hardware global sum.

  void GsumFastMode(int i)         {gsum_fast_mode = i;}
  // Sets the selection of scu frequency. 
  // 0[by default] for 25MHz. 50MHz otherwise.

  void GsumMaxTry(int i)           {gsum_max_try = i;}
  // Sets max num of tries of global sum, 2 by default.

  /*! @} */

};


/*! An instance of the GlobalJobParameter class, named GJP, should be
  created at the highest scope (outside main). This external declaration
  allows access to all global variables.
*/
extern GlobalJobParameter GJP;





#endif

CPS_END_NAMESPACE
