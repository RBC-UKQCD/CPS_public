class BGLAxisMap {
  int bgl_machine_dir_x; //!< Direction of the processor grid that X is "mapped" on. 
  int bgl_machine_dir_y; //!< Direction of the processor grid that Y is "mapped" on.
  int bgl_machine_dir_z; //!< Direction of the processor grid that Z is "mapped" on.  
  int bgl_machine_dir_t; //!< Direction of the processor grid that T is "mapped" on. 
};

class DoArg {
  int x_sites;  /*!< Global lattice dimension in the X direction*/
  int y_sites;  /*!< Global lattice dimension in the Y direction*/
  int z_sites;  /*!< Global lattice dimension in the Z direction*/
  int t_sites;  /*!< Global lattice dimension in the T direction*/
  int s_sites;  /*!< Global lattice dimension in the s (5th) direction*/

  int x_node_sites;  /*!< Local lattice dimension in the X direction*/
  int y_node_sites;  /*!< Local lattice dimension in the Y direction*/
  int z_node_sites;  /*!< Local lattice dimension in the Z direction*/
  int t_node_sites;  /*!< Local lattice dimension in the T direction*/
  int s_node_sites;  /*!< Local lattice dimension in the s (5th) direction*/
                     /*!< Relevant to domain wall fermions only. */

  int x_nodes;      /*!< Number of nodes along the X direction*/
  int y_nodes;      /*!< Number of nodes along the Y direction*/
  int z_nodes;      /*!< Number of nodes along the Z direction*/
  int t_nodes;      /*!< Number of nodes along the T direction*/
  int s_nodes;      /*!< Number of nodes along the S (5th) direction*/
                    /*!<  Relevant to domain wall fermions only. */
  int updates;
  int measurements;
  int measurefreq;
  int cg_reprod_freq; /* run reproducibility test in every nth time, never if 0 */
  BndCndType x_bc;  /*!< Boundary condition in the x direction.*/
  BndCndType y_bc;  /*!< Boundary condition in the y direction.*/
  BndCndType z_bc;  /*!< Boundary condition in the z direction.*/
  BndCndType t_bc;  /*!< Boundary condition in the t direction.*/
   
  StartConfType start_conf_kind; /*!< The kind of initial configuration*/
  unsigned long start_conf_load_addr;

  StartSeedType start_seed_kind;  /*!< The kind of initial random number generator seed*/
  string start_seed_filename<>;
  string start_conf_filename<>;
  int start_conf_alloc_flag;
  int wfm_alloc_flag;
  int wfm_send_alloc_flag;

  int start_seed_value;   /*!< The value of the random number generator seed.*/
                          /*!< Used if ::start_seed_kind is ::START_SEED_INPUT
			    or ::START_SEED_INPUT_UNIFORM. */  
    
  Float beta;       /*!< The pure gauge action "beta" parameter.*/

  Float c_1 ;       /*!< The coefficient of the rectangle term in the pure gauge action.*/
/*!<
- c_1 = 0 is the Wilson gauge action.
- c_1 = -0.05 is the tree level value.
- c_1 = -0.331 is the Iwasaki action.
*/

  Float u0;         /*!<The tadpole coefficient.*/

  Float dwf_height; /*!< The height of the domain wall.*/

  Float dwf_a5_inv; /*!< The inverse lattice spacing in the 5th direction for domain wall fermions. */

  Float power_plaq_cutoff;
     /*!< The cutoff parameter in the power plaquette term in the pure gauge action.*/

  int power_plaq_exponent; 
     /*!< The exponent in the power plaquette term in the pure gauge action.*/

  Float power_rect_cutoff;
     /*!< The cutoff parameter for the power rectangle term in the pure gauge action.*/

  int power_rect_exponent; 
    /*!< The exponent for the power rectangle term in the pure gauge action.*/
    

  int verbose_level;           /*!< Revived now!.*/
  int checksum_level;           /*!< Revived now!.*/
    
  int exec_task_list; /*!< Not used*/



  /* The following parameters are relevant to anisotropic lattices and*/
  /* clover improvement.*/
  /* xi as a prefix  indicates relevancy to anisotropic lattices.*/
  /* xi as a postfix indicates relevancy to the special anisotropic direction.*/
   
  Float xi_bare;            /*!< The bare anisotropy.*/
    int   xi_dir;             /*!< The special anisotropic direction*/
    /*!< One of 0, 1, 2, or 3 corresponding to  X, Y, Z or T respectively. */
  Float xi_v;               /*!< The bare speed of light*/
  Float xi_v_xi;            /*!< The bare speed of light in the direction of the anisotropy.*/
  Float clover_coeff;       /*!< The coefficient of the clover term in the fermion action.*/
  Float clover_coeff_xi;    /*!< The coefficient of the clover term for plaquettes with links along the anisotropic direction.*/
  /* Added for Landau gauge fixing for anisotropic lattices*/
  Float xi_gfix;            /*!< Coefficient for Landau gauge fixing.*/
  int gfix_chkb;	    /*!< Checkerboarding for the gauge fixing? */
 

    /* Parameters for the asqtad improved staggered action.*/

    /*! Coefficient of the Kogut-Susskind term in the Asqtad improved staggered fermion action.*/
    Float asqtad_KS;	
    /*! Coefficient of the Naik term in the Asqtad improved staggered fermion action.*/
    Float asqtad_naik;	
    /*! Coefficient of the 3-staple term in the Asqtad improved staggered fermion action.*/
    Float asqtad_3staple;
    /*! Coefficient of the 5-staple term in the Asqtad improved staggered fermion action.*/
    Float asqtad_5staple;
    /*! Coefficient of the 7-staple term in the Asqtad improved staggered fermion action.    */
    Float asqtad_7staple;
    /*! Coefficient of the Lepage term in the Asqtad improved staggered fermion action.*/
    Float asqtad_lepage; 

 // Parameters for the p4 improved staggered action.

    //! Coefficient of the Kogut-Susskind term in the P4 improved staggered fermion action.
    Float p4_KS;	
    //! Coefficient of the Naik term in the P4 improved staggered fermion action.
    Float p4_knight;	
    //! Coefficient of the 3-staple term in the P4 improved staggered fermion action.
    Float p4_3staple;
    //! Coefficient of the 5-staple term in the P4 improved staggered fermion action.
    Float p4_5staple;
    //! Coefficient of the 7-staple term in the P4 improved staggered fermion action.    
    Float p4_7staple;
    //! Coefficient of the Lepage term in the P4 improved staggered fermion action.
    Float p4_lepage; 
    
memfun  DoArg();
memfun  void SetupAsqTadU0(double u0);

/*!< Default values of some parameters.*/
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




