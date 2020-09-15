enum ContractionType {
  CONTRACTION_TYPE_LL_MESONS,
  CONTRACTION_TYPE_HL_MESONS,
  CONTRACTION_TYPE_O_VV_P_AA,
  CONTRACTION_TYPE_ALL_BILINEARS,
  CONTRACTION_TYPE_ALL_WALLSINK_BILINEARS_SPECIFIC_MOMENTUM,
  CONTRACTION_TYPE_FOURIER_PROP,
  CONTRACTION_TYPE_BILINEAR_VERTEX,
  CONTRACTION_TYPE_QUADRILINEAR_VERTEX,
  CONTRACTION_TYPE_TOPOLOGICAL_CHARGE,
  CONTRACTION_TYPE_MRES,
  CONTRACTION_TYPE_A2A_BILINEAR,
  CONTRACTION_TYPE_WILSON_FLOW,
  CONTRACTION_TYPE_K_TO_PIPI
};
struct ContractionTypeLLMesons{
 string prop_L<>;
 int sink_mom[3];
 string file<>;

 rpccommand GENERATE_PRINT_METHOD;
 rpccommand GENERATE_DEEPCOPY_METHOD;
};
struct ContractionTypeHLMesons{
 string prop_H<>;
 string prop_L<>;
 int sink_mom[3];
 string file<>;

 rpccommand GENERATE_PRINT_METHOD;
 rpccommand GENERATE_DEEPCOPY_METHOD;
};
struct ContractionTypeOVVpAA{
 string prop_H_t0<>;
 string prop_L_t0<>;
 string prop_H_t1<>;
 string prop_L_t1<>;
 string file<>;

 rpccommand GENERATE_PRINT_METHOD;
 rpccommand GENERATE_DEEPCOPY_METHOD;
};

/* Momentum in units of pi */
struct MomArg{
 Float p[3];

 rpccommand GENERATE_PRINT_METHOD;
 rpccommand GENERATE_DEEPCOPY_METHOD;
};
struct MomPairArg{
 Float p1[3];
 Float p2[3];

 rpccommand GENERATE_PRINT_METHOD;
 rpccommand GENERATE_DEEPCOPY_METHOD;
};

/*DaggerFlipMomentum is a dagger followed by a G-parity CC reln to swap the source momentum (only works for G-parity)*/
enum PropSuperscript {
  OpNone, 
  OpTranspose, 
  OpConj, 
  OpDagger, 
  OpFlipMomentum,
  OpTransposeFlipMomentum,
  OpConjFlipMomentum,
  OpDaggerFlipMomentum
};

/* tr(A op1(prop_1) B op2(prop_2)) for all A and B. 
   The total sink momentum is specified by the 'momenta' arg. The code loops over each given choice
   of sink momentum in the 'momenta' array.
   Operations applied to propagators chosen from the above.
*/
struct ContractionTypeAllBilinears{
  string prop_1<>;
  string prop_2<>;
  PropSuperscript op1;
  PropSuperscript op2;
  MomArg momenta<>;
  string file<>;

 rpccommand GENERATE_PRINT_METHOD;
 rpccommand GENERATE_DEEPCOPY_METHOD;
};

/* Same as above but using Fourier transformed gauge fixed propagators. User specified sink momentum applied to each of the two quarks. The total momentum is the sum of the 
   two specified.
   Note: As above, in practise the actual sink phase used when taking the Fourier transform of the first propagator is e^{-i p1 . x}. As the Hermitian
   conjugate of the resulting FT quark is taken, this ensures that the total momenta will be p1+p2.

   Option to allow cosine sinks, for which the phase is cos(p1_1 x_1)cos(p1_2 x_2)cos(p1_3 x_3) 

   Use this version, for example, when you want wallsink bilinears with zero total momentum but non-zero individual momenta.
*/

struct ContractionTypeAllWallSinkBilinearsSpecificMomentum{
  string prop_1<>;
  string prop_2<>;
  PropSuperscript op1;
  PropSuperscript op2;
  MomPairArg momenta<>;
  int cosine_sink;
  string file<>;
  
 rpccommand GENERATE_PRINT_METHOD;
 rpccommand GENERATE_DEEPCOPY_METHOD;
};


/* Fourier transformed propagator */
struct ContractionTypeFourierProp{
  string prop<>;
  int gauge_fix;
  MomArg momenta<>;
  string file<>;

 rpccommand GENERATE_PRINT_METHOD;
 rpccommand GENERATE_DEEPCOPY_METHOD;
};


/* prop_1^dag A prop_2 for all A */
struct ContractionTypeBilinearVertex{
  string prop_1<>;
  string prop_2<>;
  MomArg momenta<>;
  string file<>;

 rpccommand GENERATE_PRINT_METHOD;
 rpccommand GENERATE_DEEPCOPY_METHOD;
};

struct QuadrilinearSpinStructure{
 /*Note: Gamma are spin structures in the QDP index notation, Sigma are pauli matrix indices (0 is 2x2 unit matrix)*/
 /*Code will loop over each choice for each of the four matrices*/
 int Gamma1<>;
 int Gamma2<>;
 int Sigma1<>;
 int Sigma2<>;

 rpccommand GENERATE_PRINT_METHOD;
};

/* (prop_1^dag A prop_2) \otimes  (prop_3^dag B prop_4) for all A and B*/
struct ContractionTypeQuadrilinearVertex{
  string prop_1<>;
  string prop_2<>;
  string prop_3<>;
  string prop_4<>;

  MomArg momenta<>;
  string file<>;
  /*Specify which spin matrix choices to use - this is a good idea for G-parity*/
  /*If none are specified it will do all 16^2 possible combinations*/
  QuadrilinearSpinStructure spin_structs<>;
                                          
 rpccommand GENERATE_PRINT_METHOD;
 rpccommand GENERATE_DEEPCOPY_METHOD;
};

struct ContractionTypeTopologicalCharge{
  /*Each cycle performs 3 Ape smearing steps. Prior to each cycle the top charge is measured, and once more after all cycles have finished.*/
  int n_ape_smearing_cycles;
  /*SU(3) project the smeared links (I believe this should pretty much always be set to 1)*/ 
  int ape_smear_su3_project;

  /*Note, I would like to use ApeSmearArg directly but unions suck and do not allow members with constructors*/
  //! tolerance for the SU(3) projection
  Float ape_su3_proj_tolerance;
  //! smear in hyper-plane orthoganal to this direction
  int   ape_orthog;
  //! ape smearing coefficient 
  Float ape_coef;

  string file<>;
                                          
 rpccommand GENERATE_PRINT_METHOD;
 rpccommand GENERATE_DEEPCOPY_METHOD;
};

/*This is somewhat tricky for G-parity. We sink project to pi/L \vec{G}, which requires the propagator momentum to sum rather than cancel.  G^dag(-p) G(p)   so we need both G(-p) and G(p). 
  For cosine sources this is always true, but for other sources we must use the prop conj relation. This *requires* the source structure \eta (not including the phase factor) to obey    \sigma_2 C\gamma^5\eta^* = \eta\sigma_2 C\gamma^5
  which applies for all sources containing gauge links or the unit matrix.*/
       
struct ContractionTypeMres{
  string prop<>;  
  string file<>;
                                          
 rpccommand GENERATE_PRINT_METHOD;
 rpccommand GENERATE_DEEPCOPY_METHOD;
};

struct ContractionTypeWilsonFlow{
 /*Number of Wilson flow steps to perform*/
 int n_steps; 
 /*Wilson flow time increment*/
 Float time_step;

 string file<>;
                                          
 rpccommand GENERATE_PRINT_METHOD;
 rpccommand GENERATE_DEEPCOPY_METHOD;
};

enum A2ASmearingType {
  BOX_3D_SMEARING,
  EXPONENTIAL_3D_SMEARING
};

struct Box3dSmearing{
  int side_length;
                                          
 rpccommand GENERATE_PRINT_METHOD;
};
struct Exponential3dSmearing{
  Float radius;
                                          
 rpccommand GENERATE_PRINT_METHOD;
};

union A2ASmearing{
switch(A2ASmearingType type){
 case BOX_3D_SMEARING:
   Box3dSmearing box_3d_smearing;
 case EXPONENTIAL_3D_SMEARING:
   Exponential3dSmearing exponential_3d_smearing;
}
  rpccommand GENERATE_UNION_TYPEMAP;
  rpccommand GENERATE_DEEPCOPY_METHOD;
  rpccommand GENERATE_PRINT_METHOD;
};

struct MatIdxAndCoeff{
  int idx;
  Float coeff;
  rpccommand GENERATE_PRINT_METHOD;
};

struct ContractionTypeA2ABilinear{
  string prop_src_snk<>; /*Prop from source to sink*/
  string prop_snk_src<>; /*Prop from sink to source*/
  A2ASmearing source_smearing;
  A2ASmearing sink_smearing;

  /*Spin matrix indices run from 0..15 (QDP conventions). Resulting spin matrix is the linear combination of those matrices specified with the given coefficients*/
  /*Common examples are gamma^5 (15),  gamma^1 (1),  gamma^2 (2), gamma^3 (4), gamma^4 (8)*/
  MatIdxAndCoeff source_spin_matrix<>;
  MatIdxAndCoeff sink_spin_matrix<>;
  /*Flavor matrix indices run from 0..3 (unit matrix + 3 Pauli matrices)*/
  MatIdxAndCoeff source_flavor_matrix<>;
  MatIdxAndCoeff sink_flavor_matrix<>;

  string file<>;
                                          
 rpccommand GENERATE_PRINT_METHOD;
 rpccommand GENERATE_DEEPCOPY_METHOD;
};

struct ContractionTypeKtoPiPi{
  string prop_L<>; /*The light quark A2A propagator*/
  string prop_H<>; /*The heavy quark A2A propagator*/

  /*We wish to allow for multiple combinations of quark momentum for a given pion energy
  Thus for each pion, the user must specify the momenta for each quark in the form of a MomPairArg 'mom', 
  where mom.p1 is the momentum assigned to the 'v' field and mom.p2 to the 'w^\dagger' fields that form the A2A propagator
  The sum mom.p1 + mom.p2 is the total momentum of the pion*/

  MomPairArg p_qpi1;
  MomPairArg p_qpi2;	
  
  Float p_qK[3]; /*The momentum applied to the *strange* quark of the second pion. Minus this value is applied to the down quark such that their sum is 0*/
  /*NOTE: For G-parity the momentum components should all have the same sign and be odd-integer multiples of pi/2L*/

  int gparity_use_transconv_props; /*G-parity only: Transform the A2A propagators such that their Fourier transforms are translationally covariant*/

  A2ASmearing pion_source;
  A2ASmearing kaon_source;

  int t_sep_pi_k; /*Fixed separation between kaon and (closest) pion*/
  int t_sep_pion; /*Time separation between the two pions*/

  string file<>; /*Will have config idx appended*/

 rpccommand GENERATE_PRINT_METHOD;
 rpccommand GENERATE_DEEPCOPY_METHOD;
};


union GparityMeasurement{
switch(ContractionType type){
 case CONTRACTION_TYPE_LL_MESONS:
   ContractionTypeLLMesons contraction_type_ll_mesons;
 case CONTRACTION_TYPE_HL_MESONS:
   ContractionTypeHLMesons contraction_type_hl_mesons;
 case CONTRACTION_TYPE_O_VV_P_AA:
   ContractionTypeOVVpAA contraction_type_o_vv_p_aa;
 case CONTRACTION_TYPE_ALL_BILINEARS:
   ContractionTypeAllBilinears contraction_type_all_bilinears;
 case CONTRACTION_TYPE_ALL_WALLSINK_BILINEARS_SPECIFIC_MOMENTUM:
   ContractionTypeAllWallSinkBilinearsSpecificMomentum contraction_type_all_wallsink_bilinears_specific_momentum;
 case CONTRACTION_TYPE_FOURIER_PROP:
   ContractionTypeFourierProp contraction_type_fourier_prop;
 case CONTRACTION_TYPE_BILINEAR_VERTEX:
   ContractionTypeBilinearVertex contraction_type_bilinear_vertex;
 case CONTRACTION_TYPE_QUADRILINEAR_VERTEX:
   ContractionTypeQuadrilinearVertex contraction_type_quadrilinear_vertex;
 case CONTRACTION_TYPE_TOPOLOGICAL_CHARGE:
   ContractionTypeTopologicalCharge contraction_type_topological_charge;
 case CONTRACTION_TYPE_MRES:
   ContractionTypeMres contraction_type_mres;
 case CONTRACTION_TYPE_A2A_BILINEAR:
   ContractionTypeA2ABilinear contraction_type_a2a_bilinear;
 case CONTRACTION_TYPE_WILSON_FLOW:
   ContractionTypeWilsonFlow contraction_type_wilson_flow;
 case CONTRACTION_TYPE_K_TO_PIPI:
   ContractionTypeKtoPiPi contraction_type_k_to_pipi;
}
  rpccommand GENERATE_UNION_TYPEMAP;
  rpccommand GENERATE_DEEPCOPY_METHOD;
  rpccommand GENERATE_PRINT_METHOD;
};

class GparityContractArg{
  GparityMeasurement meas<>;
  string config_fmt<>; /* Should contain a %d which is replaced by a config index */
  int conf_start;
  int conf_incr;
  int conf_lessthan;
  FixGaugeArg fix_gauge; /* Gauge fixing - Defaults to FIX_GAUGE_NONE */

  memfun GparityContractArg();

 rpccommand GENERATE_DEEPCOPY_METHOD;
};


class GparityAMAarg{
  ContractionTypeAllBilinears bilinear_args<>;

  int exact_solve_timeslices<>; /* Specify the timeslices on which the 'rest' part is calculated */
  Float exact_precision;
  Float sloppy_precision;

  string config_fmt<>; /* Should contain a %d which is replaced by a config index */
  int conf_start;
  int conf_incr;
  int conf_lessthan;
  FixGaugeArg fix_gauge; /* Gauge fixing - Defaults to FIX_GAUGE_NONE */
  rpccommand GENERATE_DEEPCOPY_METHOD;
};

class GparityAMAbilBKarg{
  int exact_solve_timeslices<>; /* Specify the timeslices on which the 'rest' part is calculated */
  Float exact_precision;
  Float sloppy_precision;
  Float ml;
  Float mh;

  string config_fmt<>; /* Should contain a %d which is replaced by a config index */
  int conf_start;
  int conf_incr;
  int conf_lessthan;
  FixGaugeArg fix_gauge; /* Gauge fixing - Defaults to FIX_GAUGE_NONE */
  rpccommand GENERATE_DEEPCOPY_METHOD;
};

class GparityAMAarg2{
  int exact_solve_timeslices<>;
  int sloppy_solve_timeslices<>;
  int bk_tseps<>; /*K->K separations for BK. Note that you need to put sources on all t and t+tsep specified above*/
  Float exact_precision;
  Float sloppy_precision;
  Float ml;
  Float mh;

  string results_dir<>;
  string config_fmt<>; /* Should contain a %d which is replaced by a config index */
  string rng_fmt<>; /* Should contain a %d which is replaced by a config index */
  int conf_start;
  int conf_incr;
  int conf_lessthan;
  FixGaugeArg fix_gauge; /* Gauge fixing - Defaults to FIX_GAUGE_NONE */
  rpccommand GENERATE_DEEPCOPY_METHOD;
};
