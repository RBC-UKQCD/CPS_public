#include<config.h>
CPS_START_NAMESPACE
//
//  w_spect_arg.h:
//
//       structures relevant to the spectroscopy measurement of 
//       Wilson-type fermion.
//

#ifndef INCLUDED_W_SPECT_ARG_H
#define INCLUDED_W_SPECT_ARG_H

CPS_END_NAMESPACE
#include<alg/cg_arg.h>
CPS_START_NAMESPACE

//---------------------------------------------------------------------------
// The following specifies output format and filenames
//---------------------------------------------------------------------------
//
// BARYON_PAST:   same as in phys_v3.11.1 and older versions, which means 
//                information partially lost during time slice folding.
//                Print out only T/2 + 1 slices.
//                Fitting ansatz: cosh or sinh
//
// BARYON_FOLD:   average [up to a minus sign from bnd cnd for baryons]
//                between timeslice i of up Dirac projection and 
//                timeslice T - i of down projection.
//                Print out T slices.
//                Fitting ansatz: exp
//
// BARYON_RAW:    print out up and down projections seperately.
//                Print out T slices.
//                Fitting ansatz: exp
//
//
// Data Format in General
//        each row:
//             col 1 - measurement counter
//             col 2 - time slice index
//             col 3 - momentum index [0 for baryons]
//             col 4 - the real part of the hadron correlator
//             col 5 - the real part of the hadron correlator if BARYON_RAW
//        repeat above for all momenta from the same time slice
//             [only relevant to mesons] 
//        repeat above for T [BARYON_FOLD and BARYON_RAW] or 
//                         T/2+1 [mesons and BARYON_PAST] slices.
// 
//---------------------------------------------------------------------------

enum WbaryonFold {BARYON_FOLD, BARYON_RAW, BARYON_PAST};

struct WspectOutput 
{
  WbaryonFold      fold;    
  
  char *           cg;                  // datafile for cg iter & res
  char *           cg2;                 // need and extra one for second quark propagator?
  char *           pbp;                 // datafile for pbp and pbg5p
  char *           mid_point;           // datafile for mid-point sink
  char *           a0_p;             // datafile for <A_0 P>

  // the following is as in phys_v3.11.4.xiaodong, i.e. average over polarisation
  char *           a1;
  char *           b1;
  char *           pion;
  char *           pion_prime;
  char *           rho;
  char *           rho_prime;

// The a0 and a1's must be continguous because of part of the code
  char *           a0;                  // datafile for meson a0
  char *           a0_prime;
  char *           a1_x;                // x, y, z are to be understood
  char *           a1_y;                // as the three directions other
  char *           a1_z;                // than the propagation direction.
  char *           b1_x;
  char *           b1_y;
  char *           b1_z;
  char *           rho_x;
  char *           rho_y;
  char *           rho_z;
  char *           rho_x_prime;
  char *           rho_y_prime;
  char *           rho_z_prime;
  char *           nucleon;
  char *           nucleon_prime;
  char *           delta_x;
  char *           delta_y;
  char *           delta_z;
  char *           delta_t;

  WspectOutput() { pbp=0; mid_point=0; a0_p=0;}

};


//---------------------------------------------------------------------------
// The following specifies how the spectrum is measured.
//---------------------------------------------------------------------------
// One WspectArg needed for each SourceKind
//---------------------------------------------------------------------------
enum SourceKind { POINT_W = 0, 
		  WALL_W, 
		  BOX_W,  
                  JACOBI_W, //added by Thomas and Xiaodong
		  MAX_NUM_SINK  /*- N/A yet */
};


enum MomentumKind   {MOM_000 = 0,                 // average over permutaions
		     MOM_001, 
		     MOM_002, 
		     MOM_011, 
		     MOM_022, 
		     MOM_111, 
		     MOM_222, 
		     MAX_NUM_MOMENTA
};


//-----------------------------------------------------------------------------
// Related to extended mesons constructed from derivative operators
//-----------------------------------------------------------------------------
//Note: The actual derivative direction depends on propagation direction
//      See the end of w_quark.C for detail. 
//      DEV1DEV2 means DEV2 first DEV1 second !!
//The order is very important, changes will affect WspectExtendedMesons

//DEVIDEVI not used now

enum DEVOperatorKind{  
  //Group 0: Individual Operators
  UNIT=0,
  DEV1,
  DEV2,
  DEV3,
  DEV1DEV2,
  DEV2DEV1,
  DEV2DEV3,
  DEV3DEV2,
  DEV1DEV3,
  DEV3DEV1,
  DEV1DEV1,
  DEV2DEV2,
  DEV3DEV3,
  DEV_OP_NUM,
  
  //Naming convetion: F--first deriv   S--second deriv(followed by ANTISYM/SYM/DIAG)
  //Group 1: Sum over polarizations
  SUM_F,  //sum of first derivatives
  SUM_S_ANTISYM, //sum of second derivatives(antisymmtric combination: DIDJ-DJDI)
  SUM_S_SYM,     //sum of second derivatives(symmtric combination: DIDJ+DJDI)
  SUM_S_DIAG,    //sum of second derivatives(diagonal combination: DIDI)

  //Group 2: 
  SUM_F_S_ANTISYM, //SUM_F+SUM_S_ANTISYM
  //SUM_S_SYM  
  //SUM_S_DIAG

  //Group 3:
  //SUM_F_S_ANTISYM
  SUM_S_SYM_DIAG,
  //Group 4:
  SUM_UNIT_F_S_ANTISYM,
  //SUM_S_SYM_DIAG

  END_SUM_OP,

  //Related to operators with local B/E fields
  //group 1
  BEGIN_BE_OP,
  FB1_OP=0,FB2_OP,FB3_OP,
  FE1_OP,FE2_OP,FE3_OP,
  FUNIT_OP,  //for mixing
  
  //SUM operators
  //group 2
  SUM_MAGN_OP,
  SUM_ELEC_OP,
  
  //group 3
  SUM_MAGN_ELEC_OP,

  //
  END_BE_OP
};
  
  

enum WMesonOpKind{
  //each name starts with MO(Meson Operator)
  //Normal mesons
  /*  MO_a0,           MO_a0_prime,
  MO_a1_x,         MO_a1_y,         MO_a1_z,
  MO_b1_x,         MO_b1_y,         MO_b1_z,
  MO_rho_x,        MO_rho_y,        MO_rho_z,
  MO_rho_prime_x,  MO_rho_prime_y,  MO_rho_prime_z,
  MO_pion,         MO_pion_prime,
  */

  // extended mesons
  MO_a0xP_x, MO_a0xP_y, MO_a0xP_z,
  MO_pionxP_x, MO_pionxP_y, MO_pionxP_z,
  MO_a0_primexP_x, MO_a0_primexP_y, MO_a0_primexP_z,

  MO_rhoxP_A1, 
  MO_rhoxP_T1_x, MO_rhoxP_T1_y, MO_rhoxP_T1_z,
  MO_rhoxP_T2_x, MO_rhoxP_T2_y, MO_rhoxP_T2_z,

  MO_a1xP_A1,
  MO_a1xP_T2_x, MO_a1xP_T2_y, MO_a1xP_T2_z,
  MO_a1xP_E_1,MO_a1xP_E_2,

 
  MO_b1xP_T1_x, MO_b1xP_T1_y, MO_b1xP_T1_z,
 
 
  MO_b1xD_A2,
  MO_b1xD_T1_x, MO_b1xD_T1_y, MO_b1xD_T1_z,
  MO_b1xD_T2_x, MO_b1xD_T2_y, MO_b1xD_T2_z,
  MO_b1xD_E_1, MO_b1xD_E_2,

  MO_a0_primexD_x, MO_a0_primexD_y, MO_a0_primexD_z,

  MO_rhoxB_T1_x, MO_rhoxB_T1_y, MO_rhoxB_T1_z,
  MO_rhoxB_T2_x, MO_rhoxB_T2_y, MO_rhoxB_T2_z,

  MO_a1xB_A1,
  MO_a1xB_T1_x, MO_a1xB_T1_y, MO_a1xB_T1_z,
  MO_a1xB_T2_x, MO_a1xB_T2_y, MO_a1xB_T2_z, 
  MO_a1xD_A2,
  MO_a1xD_T1_x, MO_a1xD_T1_y, MO_a1xD_T1_z,
  MO_a1xD_T2_x, MO_a1xD_T2_y, MO_a1xD_T2_z, 
  MO_a1xD_E_1,  MO_a1xD_E_2,

  MO_rhoxD_A2,
  MO_rhoxD_T1_x,MO_rhoxD_T1_y,MO_rhoxD_T1_z,
  MO_rhoxD_T2_x,MO_rhoxD_T2_y,MO_rhoxD_T2_z,

  //added
  MO_pionxB_T1_x, MO_pionxB_T1_y, MO_pionxB_T1_z,
  MO_pionxD_T2_x, MO_pionxD_T2_y, MO_pionxD_T2_z,
  NUM_WMESON_OP_KIND
};

//--------------------------------
//enum WMesonState (with polarization)
//--------------------------------
enum WMesonState{
  //Normal Mesons
  /*  MS_a0,           MS_a0_prime,
  MS_a1_x,         MS_a1_y,         MS_a1_z,
  MS_b1_x,         MS_b1_y,         MS_b1_z,
  MS_rho_x,        MS_rho_y,        MS_rho_z,
  MS_rho_prime_x,  MS_rho_prime_y,  MS_rho_prime_z,
  MS_pion,         MS_pion_prime,
  */

  //ExtendedMesons
  
  MS_a0xP_x, MS_a0xP_y, MS_a0xP_z,
  MS_pionxP_x, MS_pionxP_y, MS_pionxP_z,
  MS_a0_primexP_x, MS_a0_primexP_y, MS_a0_primexP_z,

  MS_rhoxP_A1_1,
  MS_rhoxP_T1_x, MS_rhoxP_T1_y, MS_rhoxP_T1_z,
  MS_rhoxP_T2_x, MS_rhoxP_T2_y, MS_rhoxP_T2_z,
  MS_a1xP_A1_1,
  MS_a1xP_T2_x, MS_a1xP_T2_y, MS_a1xP_T2_z,
  MS_a1xP_E_1, MS_a1xP_E_2,

  MS_b1xP_T1_x, MS_b1xP_T1_y, MS_b1xP_T1_z,
  
  MS_b1xD_A2_1,
  MS_b1xD_T1_x, MS_b1xD_T1_y, MS_b1xD_T1_z,
  MS_b1xD_T2_x, MS_b1xD_T2_y, MS_b1xD_T2_z,
  MS_b1xD_E_1, MS_b1xD_E_2,

  MS_a0_primexD_x, MS_a0_primexD_y, MS_a0_primexD_z,

  MS_rhoxB_T1_x, MS_rhoxB_T1_y, MS_rhoxB_T1_z,
  MS_rhoxB_T2_x, MS_rhoxB_T2_y, MS_rhoxB_T2_z,

  MS_a1xB_A1_1,
  MS_a1xB_T1_x, MS_a1xB_T1_y, MS_a1xB_T1_z,
  MS_a1xB_T2_x, MS_a1xB_T2_y, MS_a1xB_T2_z,
  MS_a1xD_A2_1,
  MS_a1xD_T1_x,MS_a1xD_T1_y,MS_a1xD_T1_z,
  MS_a1xD_T2_x,MS_a1xD_T2_y,MS_a1xD_T2_z,
  MS_a1xD_E_1, MS_a1xD_E_2,

  MS_rhoxD_A2_1,
  MS_rhoxD_T1_x,MS_rhoxD_T1_y,MS_rhoxD_T1_z,
  MS_rhoxD_T2_x,MS_rhoxD_T2_y,MS_rhoxD_T2_z,
  
  MS_pionxB_T1_x,MS_pionxB_T1_y,MS_pionxB_T1_z,
  MS_pionxD_T2_x, MS_pionxD_T2_y, MS_pionxD_T2_z,
  // Mixing states, syntax:  source_mix_sink_polarisation
  // add desired combinations from mixing_terms.h here

  NUM_WMESON_STATE //should equal number of states in WGinfo.h?
};

//the average of all polarizations, written out to file
// the first EXTMESONS enum's have to be in the same order as in WspectOutput

#define MAX_FUZZING_C_NUM 10

enum WMesonOutputName{
//actual outputs(average over polarizarion)
  /*  a0,  a0_prime,
  a1,
  b1,
  rho,  rho_prime,
  pion,  pion_prime,
  */

  //Extended Meson
  a0xP,
  pionxP,
  a0_primexP,

  rhoxP_A1, 
  rhoxP_T1, 
  rhoxP_T2,

  a1xP_A1,
  a1xP_T2,
  a1xP_E,

  b1xP_T1, 

  b1xD_A2,
  b1xD_T1,
  b1xD_T2,
  b1xD_E,

  a0_primexD,

  rhoxB_T1,
  rhoxB_T2,


  a1xB_A1,
  a1xB_T1,
  a1xB_T2,

  a1xD_A2,
  a1xD_T1,
  a1xD_T2,
  a1xD_E,

  rhoxD_A2,
  rhoxD_T1,
  rhoxD_T2,

  pionxB_T1,
  pionxD_T2,
  // Mixing states
  // add desired combinations from mixing_terms.h here

  NUM_WMESON_OUTPUT // should be here when done with all 15 mesons
  //all polarizations
};

enum WMesonCategory{
  NORMALMESON,
  EXT_FIRSTDEV_MESON,
  EXT_SECONDDEV_SYM_MESON,
  EXT_SECONDDEV_ANTISYM_MESON,
  EXT_SECONDDEV_DIAG_MESON,
  MIXING
};


//------------------------------------------------------------------------
// Related to extended mesons constructed from local fields
//------------------------------------------------------------------------
//output index(average over all polarisations)
enum WExtMesonBEOutputName{
 BE_pionxB=0,
 BE_rhoxB_T1,
 NUM_WEXTMESON_BE_OUTPUT
};

//state IDs(with polarisation)
enum WExtMesonBEState{
  BE_MS_pionxB_x=0, BE_MS_pionxB_y, BE_MS_pionxB_z,
  BE_MS_rhoxB_T1_x, BE_MS_rhoxB_T1_y, BE_MS_rhoxB_T1_z,
  NUM_WEXTMESON_BE_STATES
};


//meson operator
enum WExtMesonBEOp{
  BE_MO_pionxB_x=0,BE_MO_pionxB_y, BE_MO_pionxB_z,
  BE_MO_rhoxB_T1_x, BE_MO_rhoxB_T1_y, BE_MO_rhoxB_T1_z,
  NUM_WEXTMESON_BE_OPS
};

//extmeson category

enum WExtMesonBECategory{
  ELEC_HYBRID_BE=0,
  MAG_HYBRID_BE,
  MIXING_BE
};

enum FieldTensorId{
  //used in state table
  //group 1
  FB1=0,FB2,FB3,
  FE1,FE2,FE3,
  NUM_FLDS,
  FUNIT,  //for mixing
  
  //SUM operators
  //group 2
  SUM_MAGN,
  SUM_ELEC,
  
  //group 3
  SUM_MAGN_ELEC,

  //
  NUM_FLD_OPS
  
};
  




//-------------------------------------------------------------------------
// struct WspectArg
//-------------------------------------------------------------------------

struct WspectArg {

  // The conjugate gradient argument for the quark propagator calculation
  CgArg cg;                     

  // Propagation direction. [0..3] as [x,y,z,t]
  int prop_dir;               

  // Apply multi-sinking from one calcuation of quark prop. - N/A yet.
  //  int num_sink;             
  
  // Apply multi-momenta projection.
  int num_mom;
				
  // Type of source for the calculation of the quark propagator
  SourceKind source_kind;       

  //tail of output file names 
  char *filetail;//eg. P, W, B4, JN30E2.50,... 

  //--------------------------
  //Paremeter of sources
  //--------------------------
  // Relevant only for BOX_W
  //
  // Note: 1. begin[prop_dir] and end[prop_dir] do not matter.
  //       2. In all directions including prop_dir, we check that
  //          0 <= begin[] <= end[] <= total_sites
  //       3. We check that at least in one direction other than prop_dir
  //          0 <= begin[] <  end[] <= total_sites
  int src_box_b[4];           
  int src_box_e[4];      


  // int snk_box_b[4];          - N/A yet. 
  // int snk_box_e[4];          - N/A yet.


  //added by Thomas and Xiaodong
  //relevant to Jacobi sources
  Float g_epsi;
  int g_n;
  int g_center[4]; //center of the Gaussian source
      

  // rescaling the source vector to get the full
  // range of single precision arithmetic
  // important for heavy quarks where the strong exponential
  // decay of propagators may be out of single precision range
  // default value = 1.0 (no rescaling)
  Float rescale_factor;

  // AOTS: time slices to put source on for quark propagators
  //       time slices are in arithmetic series [start, ...]
  //       num = 0 ==> no inversion will be done.
  int aots_num;                    
  int aots_start;                  
  int aots_step;                   


  //------------------------------
  //Mesurement control flags
  //------------------------------
  int baryons_on;                
  int normal_mesons_on;
  int extended_mesons_on;
  int extended_mesonsBE_on; //extended mesons from local B/E fields

  //related to extended_mesons using derivative operators
  //See enum DEVOpKind for explanation
  int extended_mesons_op_groupId; //use source operator sum only
  int extended_mesons_first_dev_on; 
  int extended_mesons_second_sym_dev_on; //do second derivative src_op
  int extended_mesons_second_antisym_dev_on;
  int extended_mesons_second_diag_dev_on;

  //Fuzzing parameters
  //-------------------------------
  int fuzzing_on; // 1/0(on/off)
  int sink_fuzzing_only; //1(yes): use fuzzed links at sink only. 0(no): use at both sink and source
  int fuzzing_level;
  int fuzzing_c_num;//number of coefficients to run
  Float fuzzing_c[MAX_FUZZING_C_NUM]; //multiplier
  int fuzzing_hits; //cabbobo hits


  //related to extended_mesonsBE
  int extended_mesonsBE_op_groupId; 
  int extended_mesonsBE_Elec_on;
  int extended_mesonsBE_Magn_on;

  //Fuzzing parameters for extended_mesonsBE
  int BEfuzzing_on; // 1/0(on/off)
  int BEfuzzing_level;
  int BEfuzzing_c_num;//number of coefficients to run
  Float BEfuzzing_c[MAX_FUZZING_C_NUM]; //multiplier
  int BEfuzzing_hits; //cabbobo hits

};

#endif /* !INCLUDED_W_SPECT_ARG_H */
CPS_END_NAMESPACE
