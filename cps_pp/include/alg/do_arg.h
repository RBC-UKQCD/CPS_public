#include<config.h>
CPS_START_NAMESPACE
//----------------------------------------------------------------------
/*!\file
  \brief  Definition of the DoArg structure.
  
  $Id: do_arg.h,v 1.8 2004-06-02 09:36:38 zs Exp $
*/
//--------------------------------------------------------------------
/* CIM Sun Jul  6 23:30:27 GMT 1997 */
/* This file contains the  structure called DoArg where 
   the arguments on the #DO line are stored during program execution. */
//--------------------------------------------------------------------
#ifndef INCLUDED_DO_ARG_H
#define INCLUDED_DO_ARG_H



CPS_END_NAMESPACE
#include <util/gjp.h>
#include <util/vector.h>
CPS_START_NAMESPACE
#ifdef PARALLEL
CPS_END_NAMESPACE
#include <comms/sysfunc.h>
CPS_START_NAMESPACE
#endif

//! A structure to contain the initial values of run-time global variables.

struct DoArg {
  int x_node_sites;  //!< Local lattice dimension in the X direction
  int y_node_sites;  //!< Local lattice dimension in the Y direction
  int z_node_sites;  //!< Local lattice dimension in the Z direction
  int t_node_sites;  //!< Local lattice dimension in the T direction
  int s_node_sites;  //!< Local lattice dimension in the s (5th) direction
                     /*!< Relevant to domain wall fermions only. */

  int x_nodes;      //!< Number of nodes along the X direction
  int y_nodes;      //!< Number of nodes along the Y direction
  int z_nodes;      //!< Number of nodes along the Z direction
  int t_nodes;      //!< Number of nodes along the T direction
  int s_nodes;      //!< Number of nodes along the S (5th) direction
                    /*!<  Relevant to domain wall fermions only. */
#ifdef PARALLEL
  SCUAxis s_axis;   //!< The machine axis on which the 5th direction axis is mapped.
    /*!< Relevant only to domain wall fermions when the number
      of nodes along the 5th direction is greater than 1.
      It should take one of the values enumerated in ::SCUAxis.
    */
#endif

  BndCndType x_bc;  //!< Boundary condition in the x direction.
  BndCndType y_bc;  //!< Boundary condition in the y direction.
  BndCndType z_bc;  //!< Boundary condition in the z direction.
  BndCndType t_bc;  //!< Boundary condition in the t direction.
   
  StartConfType start_conf_kind; //!< The kind of initial configuration

  Matrix *start_conf_load_addr; //!< The address of the starting configuration
                                /*!< Used if ::start_conf_kind is
				  ::START_CONF_MEM or ::START_CONF_LOAD. */

  StartSeedType start_seed_kind;  //!< The kind of initial random number generator seed

  int start_seed_value;   //!< The value of the random number generator seed.
                          /*!< Used if ::start_seed_kind is ::START_SEED_INPUT
			    or ::START_SEED_INPUT_UNIFORM. */  

    
  Float beta;       //!< The pure gauge action "beta" parameter.

  Float c_1 ;       //!< The coefficient of the rectangle term in the pure gauge action.
/*!<
- c_1 = 0 is the Wilson gauge action.
- c_1 = -0.05 is the tree level value.
- c_1 = -0.331 is the Iwasaki action.
*/

  Float u0;         //!<The tadpole coefficient.

  Float dwf_height; //!< The height of the domain wall.

  Float dwf_a5_inv; //!< The inverse lattice spacing in the 5th direction for domain wall fermions. 

  Float power_plaq_cutoff;
     //!< The cutoff parameter in the power plaquette term in the pure gauge action.

  int power_plaq_exponent; 
     //!< The exponent in the power plaquette term in the pure gauge action.

  Float power_rect_cutoff;
     //!< The cutoff parameter for the power rectangle term in the pure gauge action.

  int power_rect_exponent; 
    //!< The exponent for the power rectangle term in the pure gauge action.
    

    
  int exec_task_list; //!< Not used


  //------------------------------------------------------------------
  // Added in by Ping for anisotropic lattices and clover improvement
  //------------------------------------------------------------------
  // The following parameters are relevant to anisotropic lattices and
  // clover improvement.
  // xi as a prefix  indicates relevancy to anisotropic lattices.
  // xi as a postfix indicates relevancy to the special anisotropic direction.
   
  Float xi_bare;            //!< The bare anisotropy.
    int   xi_dir;             //!< The special anisotropic direction
    /*!< One of 0, 1, 2, or 3 corresponding to  X, Y, Z or T respectively. */
  Float xi_v;               //!< The bare speed of light
  Float xi_v_xi;            //!< The bare speed of light in the direction of the anisotropy.
  Float clover_coeff;       //!< The coefficient of the clover term in the fermion action.
  Float clover_coeff_xi;    //!< The coefficient of the clover term for plaquettes with links along the anisotropic direction.
  // Added for Landau gauge fixing for anisotropic lattices
  Float xi_gfix;            //!< Coefficient for Landau gauge fixing.
 
  //------------------------------------------------------------------
  // Added in by Ping for global sum
  //------------------------------------------------------------------
  // The following two parameters are relevant to scu transfer frequency
  // and hardware global sum.
  int  gsum_fast_mode;      //!< 0 (by default) for 25MHz; 50MHz otherwise.
  int  gsum_max_try;        //!< Maximum number of tries of global sum; 2 by default.

    // Parameters for the asqtad improved staggered action.

    //! Coefficient of the Kogut-Susskind term in the Asqtad improved staggered fermion action.
    Float asqtad_KS;	
    //! Coefficient of the Naik term in the Asqtad improved staggered fermion action.
    Float asqtad_naik;	
    //! Coefficient of the 3-staple term in the Asqtad improved staggered fermion action.
    Float asqtad_3staple;
    //! Coefficient of the 5-staple term in the Asqtad improved staggered fermion action.
    Float asqtad_5staple;
    //! Coefficient of the 7-staple term in the Asqtad improved staggered fermion action.    
    Float asqtad_7staple;
    //! Coefficient of the Lepage term in the Asqtad improved staggered fermion action.
    Float asqtad_lepage; 
    

  //------------------------------------------------------------------
  // Added in by Ping for anisotropy, clover and global sum.
  //------------------------------------------------------------------
  // Default CTOR  [inline OK, since it will only be called once in practice.]

  DoArg()     :  
      s_node_sites(1),
	s_nodes(1),
    xi_bare(1.0),           
    xi_dir(3),              
    xi_v(1.0),              
    xi_v_xi(1.0),           
    clover_coeff(0.0),      
    clover_coeff_xi(0.0),   
    c_1(0.0),               
    u0(1.0),                
    xi_gfix(1.0),           
    gsum_fast_mode(0),      
    gsum_max_try(2),        
	dwf_a5_inv(1.0),
	asqtad_naik(0.0), asqtad_3staple(0.0),asqtad_5staple(0.0),
	asqtad_7staple(0.0),asqtad_lepage(0.0)
	{}
    
//!< Default values of some parameters.
/*!<
- Default local lattice dimension in the s (5th) direction: 1 (4-dim. lattice)
- Default processor grid dimension in the s (5th) direction: 1
- Default bare anisotropy: 1 (isotropic lattice)
- Default anisotropic direction: time	       
- Default bare speed of light: 1	       
- Default bare speed of light in the anisotropic direction:  1	       
- Default clover term coefficient: 0	       
- Default anisotropic clover term coefficient: 0  
- Default Iwazaki coefficient: 0		       
- Default tadpole coefficient: 1 (tree level value)
- Default coefficient for gauge fixing: 1       
- Default scu transfer freqency:   25MHz	       
- Default maximum number of global sum attempts:   2	       
- Default 5th direction inverse lattice spacing: 1
- Default Asqtad improved staggered fermion action parameters: 0.
*/
    
};



#endif /* !INCLUDED_DO_ARG_H */






CPS_END_NAMESPACE
