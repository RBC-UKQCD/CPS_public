#include<config.h>
CPS_START_NAMESPACE
/* CIM Sun Jul  6 23:30:27 GMT 1997 */
/* This file contains the  structure called DoArg where 
   the arguments on the #DO line are stored during program execution. */

#ifndef INCLUDED_DO_ARG_H
#define INCLUDED_DO_ARG_H

CPS_END_NAMESPACE
#include<util/gjp.h>
#include<util/vector.h>
#include<config.h>
CPS_START_NAMESPACE
#ifdef PARALLEL
CPS_END_NAMESPACE
#include <sysfunc.h>
CPS_START_NAMESPACE
#endif

struct DoArg {
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

  BndCndType x_bc;  // Boundary conditions along x
  BndCndType y_bc;  // Boundary conditions along y
  BndCndType z_bc;  // Boundary conditions along z
  BndCndType t_bc;  // Boundary conditions along t
   
  StartConfType start_conf_kind; // The kind of initial 
                                 // configuration

  Matrix *start_conf_load_addr; // The address of
                                // the starting configuration
                                // to be used if 
                                // start_conf = START_CONF_MEM

  StartSeedType start_seed_kind;  // The kind of initial 
                                  // random generator seed

  int start_seed_value;   // The initial value of seed
                          // if start_seed_kind
                          // is START_SEED_INPUT or
                          // is START_SEED_INPUT_UNIFORM

  int colors;       // The number of colors.

  Float beta;       // The pure gauge action "beta"

  Float c_1 ;       // Related to the coefficient of the rectangle
                    // term in the GimprRectFnone action.
                    // c_1 = 0 is the Wilson gauge action.
                    // c_1 = -0.05 is the tree level value.
                    // c_1 = -0.331 is the Iwasaki action.

  Float u0;         //the tadpole

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

  Float power_rect_cutoff;
     // The cutoff parameter for the power
     // rectangle term in the pure gauge action
     // implemented by the GpowerRect class.

  int power_rect_exponent; 
     // The exponent for the power
     // rectangle term in the pure gauge action
     // implemented by the GpowerRect class.

  int verbose_level;  
     // verbose_level = 0 : No output
     // verbose_level = 1 : ...

  int exec_task_list; // number of task list loops


  //------------------------------------------------------------------
  // Added in by Ping for anisotropic lattices and clover improvement
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
  Float xi_gfix;            // Coefficient for Landau gauge fixing
 
  //------------------------------------------------------------------
  // Added in by Ping for global sum
  //------------------------------------------------------------------
  // The following two parameters are relevant to scu transfer frequency
  // and hardware global sum.
  int  gsum_fast_mode;      // 0[by default] for 25MHz. 50MHz otherwise.
  int  gsum_max_try;        // max num of tries of global sum, 2 by default

  //------------------------------------------------------------------
  // Added in by Ping for anisotropy, clover and global sum.
  //------------------------------------------------------------------
  // Default CTOR  [inline OK, since it will only be called once in practice.]
  DoArg()     :  
    s_node_sites(1),        // default number of nodes in s:    1
    xi_bare(1.0),           // default bare anisotropy:         isotropic
    xi_dir(3),              // default anisotropic direction:   time
    xi_v(1.0),              // default bare velocity of light:  c
    xi_v_xi(1.0),           // default bare velocity of light:  c
    clover_coeff(0.0),      // default clover term coefficient: 0
    clover_coeff_xi(0.0),   // default clover term coefficient: 0
    c_1(0.0),               // default Iwazaki coefficient
    u0(1.0),                // default tadpole
    xi_gfix(1.0),           // default coefficient for gauge fixing: 1
    gsum_fast_mode(0),      // default scu transfer freqency:   25MHz
    gsum_max_try(2),        // default max num of gsum tries:   2
    dwf_a5_inv(1.0)         // default dwf inverse 5th dir. lattice spacing   
    {
    }
};

#endif /* !INCLUDED_DO_ARG_H */





CPS_END_NAMESPACE
