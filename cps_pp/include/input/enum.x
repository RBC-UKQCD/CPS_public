#ifndef VMLH
#include"precision.h"
#endif
/*!\file
  \brief  Magic numbers.

  $Id: enum.x,v 1.2 2004-12-11 20:57:33 chulwoo Exp $
*/
/*--------------------------------------------------------------------*/
/*  CVS keywords*/
/**/
/*  $Author: chulwoo $*/
/*  $Date: 2004-12-11 20:57:33 $*/
/*  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/input/enum.x,v 1.2 2004-12-11 20:57:33 chulwoo Exp $*/
/*  $Id: enum.x,v 1.2 2004-12-11 20:57:33 chulwoo Exp $*/
/*  $Name: not supported by cvs2svn $*/
/*  $Locker:  $*/
/*  $RCSfile: enum.x,v $*/
/*  $Revision: 1.2 $*/
/*  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/input/enum.x,v $*/
/*  $State: Exp $*/
/**/
/*--------------------------------------------------------------------*/




/*------------------------------------------------------------------*/
/*! The directions in the lattice*/
/*------------------------------------------------------------------*/
enum DirType {DIR_X,
	      DIR_Y,
	      DIR_Z,
	      DIR_T,
	      DIR_S};

/*------------------------------------------------------------------*/
/*! The types of fermion action*/
/*------------------------------------------------------------------*/
enum FclassType {
    F_CLASS_NONE,   
    F_CLASS_STAG,
    F_CLASS_WILSON,
    F_CLASS_CLOVER,
    F_CLASS_DWF,
    F_CLASS_ASQTAD,
    F_CLASS_P4
};


/*------------------------------------------------------------------*/
/*! The types of gauge action.*/
/*------------------------------------------------------------------*/
enum GclassType {
    G_CLASS_NONE,   
    G_CLASS_WILSON,
    G_CLASS_POWER_PLAQ,
    G_CLASS_IMPR_RECT,
    G_CLASS_POWER_RECT,
    G_CLASS_IMPR_OLSYM
};


/*------------------------------------------------------------------*/
/*! The lattice storage orders.*/
/*------------------------------------------------------------------*/
enum StrOrdType {
    CANONICAL = 0, /*!< Canonical storage order: %Lattice sites are ordered
		     so that the site with cartesian coordinates
		     (x, y, z, t) on a local lattice of dimensions
	     (N<sub>x</sub>, N<sub>y</sub>, N<sub>z</sub>, N<sub>t</sub>)
		     has the index
		     
		     n = x + N<sub>x</sub> y + N<sub>x</sub>N<sub>y</sub> z + N<sub>x</sub>N<sub>y</sub>N<sub>z</sub> t

		     <em>i.e</em> the x coordinate runs fastest, then y,
		     then z and t slowest.
		     The gauge link direction index is
		     0, 1, 2 or 3  for direction X, Y, Z and T respectively.

		     This order is expected by most Lattice functions.
		   */
    STAG      = 1,  /*!<
		      This order is used by actions with staggered fermions.
		      Gauge links are stored in the CANONICAL site order,
		      but they are multiplied by the staggered phases and
		      the hermitian conjugate of each link  is stored.
		      
		      Fermions are stored in an odd-even order, where lattice
		      sites are ordered by parity.
		      The parity of a site with cartesian coordinates
		      (x, y, z, t) is even(odd) if x+y+z+t is even(odd)
		      Even parity sites are numbered before odd parity sites.

			On a local lattice of dimensions
	     (N<sub>x</sub>, N<sub>y</sub>, N<sub>z</sub>, N<sub>t</sub>)
		     has the index

		     n = [t + N<sub>t</sub> x +
		     N<sub>t</sub> N<sub>x</sub> y +
		     N<sub>t</sub> N<sub>x</sub> N<sub>y</sub> z
   + N<sub>t</sub> N<sub>x</sub> N<sub>y</sub> N<sub>z</sub> ((x+y+z+t)%2)]/2

   		     <em>i.e</em> the t coordinate runs fastest, then x,
		     then y and z slowest.
		    */
    WILSON    = 2,  /*!<
		      This is used with actions with wilson-type fermions.
		      %Lattice sites for both gauge and fermion fields
		      are ordered by parity.
		      The parity of a site with cartesian coordinates
		      (x, y, z, t) is even(odd) if x+y+z+t is even(odd)

		      For gauge fields, 
		      even parity sites are numbered before odd parity sites.
      
		      On a local lattice of dimensions
	      (N<sub>x</sub>, N<sub>y</sub>, N<sub>z</sub>, N<sub>t</sub>)
		      each site has the index
		      n = [x-x%2 + N<sub>x</sub> y +
		      N<sub>x</sub> N<sub>y</sub> z +
		      N<sub>x</sub> N<sub>y</sub> N<sub>z</sub> t
     + N<sub>x</sub> N<sub>y</sub> N<sub>z</sub> N<sub>t</sub> ((x+y+z+t)%2)]/2

		      For fermion fields,
		      odd parity sites are numbered before even parity sites.

		      Each site has the index
		      n = [x-x%2 + N<sub>x</sub> y +
		      N<sub>x</sub> N<sub>y</sub> z +
		      N<sub>x</sub> N<sub>y</sub> N<sub>z</sub> t
     + N<sub>x</sub> N<sub>y</sub> N<sub>z</sub> N<sub>t</sub> ((x+y+z+t+1)%2)]/2
     
		     <em>i.e</em> in both cases, the x coordinate runs fastest,
		     then y, then z and t slowest.
     
		    */
    G_WILSON_HB  = 3 /*!< Storage order for the %Wilson gauge action heat bath:
		        Site ordering is CANONICAL but the link direction indices  
		       are 0, 1, 2 or 3  for direction T, X, Y and Z
		       respectively.
		       This order is expected by all QncWilsonHb functions.*/
};


/*------------------------------------------------------------------*/
/*! The  storage order conversion flags.*/
/*------------------------------------------------------------------*/
enum CnvFrmType {
    CNV_FRM_NO   = 0,   /*!< Do not convert fermion field */
    CNV_FRM_YES  = 1  /*!< Convert fermion field assuming that it is defined
			on the whole lattice \e i.e. one spinor per site and
			not just on the sites of a single checkerboard
			parity.*/
};

/*------------------------------------------------------------------*/
/*!  The dimensionality of the fermion field.*/
/*------------------------------------------------------------------*/
enum FermionFieldDimension { FOUR_D, FIVE_D };

/*------------------------------------------------------------------*/
/*! The kinds of preservation */
/*------------------------------------------------------------------*/
enum PreserveType {PRESERVE_NO  = 0,  /*!< Do not preserve. */
		   PRESERVE_YES = 1 /*!< Preserve. */
};

/*------------------------------------------------------------------*/
/*! The kinds of starting configurations.*/
/**/
/* START_CONF_ORD  Ordered start. After the configuration is set*/
/* the GJP.StartConfKind() is set to START_CONF_MEM.*/
/**/
/* START_CONF_DISORD  Disordered start. After the configuration*/
/* is set the GJP.StartConfKind() is set to START_CONF_MEM.*/
/**/
/* START_CONF_FILE Read from file. After the configuration is set*/
/* the GJP.StartConfKind() is set to START_CONF_MEM.*/
/**/
/* START_CONF_LOAD Memory is not allocated for the gauge field.*/
/* Instead the gauge_field pointer is set to GJP.StartConfLoadAddr().*/
/* After the configuration is set the GJP.StartConfKind() is set */
/* to START_CONF_MEM.*/
/**/
/* START_CONF_MEM Memory is allocated for the gauge field with */
/* pmalloc. When the program begins execution it checcks if the*/
/* pmalloc address is the same as the GJP.StartConfLoadAddr(). */
/* If not it exits with an error. Subsequent calls to the */
/* constructor with START_CONF_MEM do nothing except to set the */
/* gauge field to canonical order.*/
/**/
/* For implementation details see the Lattice constructor in */
/* util/lattice/lattice_base/lattice_base.C*/
/**/
/*------------------------------------------------------------------*/

enum StartConfType {
    START_CONF_ORD     = 0, /*!< Ordered start */
    START_CONF_DISORD  = 1, /*!< Disordered start */
    START_CONF_FILE    = 2, /*!< Read from file */
    START_CONF_LOAD    = 3, /*!< The gauge_field pointer is
			      set to the address specified in
			      the DoArg structure but no memory
			      is allocated for it.*/
    START_CONF_MEM     = 4  /*!< The gauge_field pointer is
			      set to the address specified in
			      the DoArg structure and memory
			      is allocated for it. */
};


/*------------------------------------------------------------------*/
/*! The possible kinds of initial RNG seeds.*/
/*------------------------------------------------------------------*/
enum StartSeedType {
    START_SEED_FIXED          = 0,
/*!< Seed has a fixed value but it is different on each node */
    START_SEED_FIXED_UNIFORM  = 1,
/*!< Seed has a fixed value and is the same on all nodes */
    START_SEED                = 2,
/*!< Seed has a time-dependent value and it is different on each node */
    START_SEED_UNIFORM        = 3,
/*!< Seed has a time-dependent value and is the same on all nodes */     
    START_SEED_INPUT          = 4,
/*!< Seed value comes from the DoArg structure and is different on all nodes*/
    START_SEED_INPUT_UNIFORM  = 5,
/*!< Seed value comes from the DoArg structure and is the same on all nodes*/ 
    START_SEED_INPUT_NODE     = 6
/*!< Seed structure comes from the DoArg structure and is different on all
  nodes. This should be used when in a farming mode. */     
};

/*------------------------------------------------------------------*/
/* The possible types of checker  boards (chequer boards?)*/
/*------------------------------------------------------------------*/
enum ChkbType {
    CHKB_EVEN =0,   /* Even checkerboard*/
    CHKB_ODD  =1};  /* Odd checkerboard*/


/*------------------------------------------------------------------*/
/*! Hermitian conjugate flag*/
/*------------------------------------------------------------------*/
enum DagType {
    DAG_NO  = 0,   /*!< No hermitian conjugate */
    DAG_YES = 1    /*!< Hermitian conjugate */
};


/*------------------------------------------------------------------*/
/*! The kinds of boundary conditions.*/
/*------------------------------------------------------------------*/
enum BndCndType {
    BND_CND_PRD,    /*!< Periodic */
    BND_CND_APRD    /*!< Antiperiodic */
};


/*------------------------------------------------------------------*/
/*! The possible kinds of gauge fixing*/
/*------------------------------------------------------------------*/
enum FixGaugeType{ 
    FIX_GAUGE_NONE = -2,   /*!< No gauge fixing */
    FIX_GAUGE_LANDAU = -1,  /*!< Fixing to Landau gauge */
    FIX_GAUGE_COULOMB_X = 0,  /*!< Fixing to Coulomb gauge */
    FIX_GAUGE_COULOMB_Y = 1,  /*!< Fixing to Coulomb gauge */
    FIX_GAUGE_COULOMB_Z = 2,  /*!< Fixing to Coulomb gauge */
    FIX_GAUGE_COULOMB_T = 3 /*!< Fixing to Coulomb gauge */
};		   

/*------------------------------------------------------------------*/
/*! The possible kinds of spin projections*/
/*------------------------------------------------------------------*/
enum SprojType {
    SPROJ_XM = 0,     /*!< (1 - gamma_0) projection */
    SPROJ_YM = 1,     /*!< (1 - gamma_1) projection */
    SPROJ_ZM = 2,     /*!< (1 - gamma_2) projection */
    SPROJ_TM = 3,     /*!< (1 - gamma_3) projection */
    SPROJ_XP = 4,     /*!< (1 + gamma_0) projection */
    SPROJ_YP = 5,     /*!< (1 + gamma_1) projection */
    SPROJ_ZP = 6,     /*!<(1 + gamma_2) projection */
    SPROJ_TP = 7   /*!< (1 + gamma_3)  projection */
};

/*------------------------------------------------------------------*/
/*! The possible kinds of Sigma projections*/
/*------------------------------------------------------------------*/
enum SigmaprojType {
    SIGMAPROJ_XY = 0,     /*!< Sigma_{0,1} projection */
    SIGMAPROJ_XZ = 1,     /*!< Sigma_{0,2} projection */
    SIGMAPROJ_XT = 2,     /*!< Sigma_{0,3} projection */
    SIGMAPROJ_YZ = 3,     /*!< Sigma_{1,2} projection */
    SIGMAPROJ_YT = 4,     /*!< Sigma_{1,3} projection */
    SIGMAPROJ_YX = 5,     /*!< Sigma_{1,0}     projection */
    SIGMAPROJ_ZT = 6,     /*!< Sigma_{2,3} projection */
    SIGMAPROJ_ZX = 7,     /*!< Sigma_{2,0} projection */
    SIGMAPROJ_ZY = 8,     /*!< Sigma_{2,1}   projection */
    SIGMAPROJ_TX = 9,     /*!< Sigma_{3,0} projection */
    SIGMAPROJ_TY =10,     /*!< Sigma_{3,1} projection */
    SIGMAPROJ_TZ =11     /*!< Sigma_{3,2}  projection */
};

/*------------------------------------------------------------------*/
/*! Operators, in terms of the fermion matrix \e M, of which the eigenvalues/vectors can be found.*/
/*------------------------------------------------------------------*/
enum RitzMatType {
    NONE,            /*!< No eigenvalues requested */
    MAT_HERM,        /*!< The hermitian matrix on the full lattice */
    MATPC_HERM,      /*!< The preconditioned hermitian matrix on a single parity. */
    MATPCDAG_MATPC,  /*!< The preconditioned \f$M^\dagger M\f$
		          on a single parity */
    NEG_MATPCDAG_MATPC,  /*!< The preconditioned \f$-M^\dagger M\f$
		          on a single parity */
    MATDAG_MAT,      /*!< \f$M^\dagger M\f$ on the full lattice */
    NEG_MATDAG_MAT,   /*!< \f$-M^\dagger M\f$ on the full lattice*/
    MATDAG_MAT_NORM,  /*!< \f$cM^\dagger M\f$ on the full lattice (normalised)*/
    NEG_MATDAG_MAT_NORM  /*!< \f$-cM^\dagger M\f$ on the full lattice (normalised)*/
};

/*------------------------------------------------------------------*/
/*! The types of rational approximation for RHMC*/
/*------------------------------------------------------------------*/
enum RatApproxType {
  CONSTANT,  /*!< CONSTANT - approximation is held constant*/
  DYNAMIC    /*!< DYNAMIC - approximation recalculated at the end of the trajectory*/
};
		  
/*------------------------------------------------------------------*/
/*! The types multishift solve to perform*/
/*------------------------------------------------------------------*/
enum MultiShiftSolveType {
  SINGLE,  /*!< SINGLE - solutions are summed to a common vector*/
  MULTI    /*!< MULTI - solutions are summed to separate vectors*/
};
enum WbaryonFold {BARYON_FOLD, BARYON_RAW, BARYON_PAST};
enum SourceKind { POINT_W = 0, 
		  WALL_W, 
		  BOX_W,  
                  JACOBI_W, /*added by Thomas and Xiaodong*/
		  MAX_NUM_SINK  /*- N/A yet */
};


enum MomentumKind   {MOM_000 = 0,                 /* average over permutaions*/
		     MOM_001, 
		     MOM_002, 
		     MOM_011, 
		     MOM_022, 
		     MOM_111, 
		     MOM_222, 
		     MAX_NUM_MOMENTA
};


/*-----------------------------------------------------------------------------*/
/* Related to extended mesons constructed from derivative operators*/
/*-----------------------------------------------------------------------------*/
/*Note: The actual derivative direction depends on propagation direction*/
/*      See the end of w_quark.C for detail. */
/*      DEV1DEV2 means DEV2 first DEV1 second !!*/
/*The order is very important, changes will affect WspectExtendedMesons*/

/*DEVIDEVI not used now*/

enum DEVOperatorKind{  
  /*Group 0: Individual Operators*/
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
  
  /*Naming convetion: F--first deriv   S--second deriv(followed by ANTISYM/SYM/DIAG)*/
  /*Group 1: Sum over polarizations*/
  SUM_F,  /*sum of first derivatives*/
  SUM_S_ANTISYM, /*sum of second derivatives(antisymmtric combination: DIDJ-DJDI)*/
  SUM_S_SYM,     /*sum of second derivatives(symmtric combination: DIDJ+DJDI)*/
  SUM_S_DIAG,    /*sum of second derivatives(diagonal combination: DIDI)*/

  /*Group 2: */
  SUM_F_S_ANTISYM, /*SUM_F+SUM_S_ANTISYM*/
  /*SUM_S_SYM  */
  /*SUM_S_DIAG*/

  /*Group 3:*/
  /*SUM_F_S_ANTISYM*/
  SUM_S_SYM_DIAG,
  /*Group 4:*/
  SUM_UNIT_F_S_ANTISYM,
  /*SUM_S_SYM_DIAG*/

  END_SUM_OP,

  /*Related to operators with local B/E fields*/
  /*group 1*/
  BEGIN_BE_OP,
  FB1_OP=0,FB2_OP,FB3_OP,
  FE1_OP,FE2_OP,FE3_OP,
  FUNIT_OP,  /*for mixing*/
  
  /*SUM operators*/
  /*group 2*/
  SUM_MAGN_OP,
  SUM_ELEC_OP,
  
  /*group 3*/
  SUM_MAGN_ELEC_OP,

  /**/
  END_BE_OP
};
  
  

enum WMesonOpKind{
  /*each name starts with MO(Meson Operator)*/
  /*Normal mesons*/
  /*  MO_a0,           MO_a0_prime,
  MO_a1_x,         MO_a1_y,         MO_a1_z,
  MO_b1_x,         MO_b1_y,         MO_b1_z,
  MO_rho_x,        MO_rho_y,        MO_rho_z,
  MO_rho_prime_x,  MO_rho_prime_y,  MO_rho_prime_z,
  MO_pion,         MO_pion_prime,
  */

  /* extended mesons*/
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

  /*added*/
  MO_pionxB_T1_x, MO_pionxB_T1_y, MO_pionxB_T1_z,
  MO_pionxD_T2_x, MO_pionxD_T2_y, MO_pionxD_T2_z,
  NUM_WMESON_OP_KIND
};

/*--------------------------------*/
/*enum WMesonState (with polarization)*/
/*--------------------------------*/
enum WMesonState{
  /*Normal Mesons*/
  /*  MS_a0,           MS_a0_prime,
  MS_a1_x,         MS_a1_y,         MS_a1_z,
  MS_b1_x,         MS_b1_y,         MS_b1_z,
  MS_rho_x,        MS_rho_y,        MS_rho_z,
  MS_rho_prime_x,  MS_rho_prime_y,  MS_rho_prime_z,
  MS_pion,         MS_pion_prime,
  */

  /*ExtendedMesons*/
  
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
  /* Mixing states, syntax:  source_mix_sink_polarisation*/
  /* add desired combinations from mixing_terms.h here*/

  NUM_WMESON_STATE /*should equal number of states in WGinfo.h?*/
};

/*the average of all polarizations, written out to file*/
/* the first EXTMESONS enum's have to be in the same order as in WspectOutput*/

#define MAX_FUZZING_C_NUM 10

enum WMesonOutputName{
/*actual outputs(average over polarizarion)*/
  /*  a0,  a0_prime,
  a1,
  b1,
  rho,  rho_prime,
  pion,  pion_prime,
  */

  /*Extended Meson*/
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
  /* Mixing states*/
  /* add desired combinations from mixing_terms.h here*/

  NUM_WMESON_OUTPUT /* should be here when done with all 15 mesons*/
  /*all polarizations*/
};

enum WMesonCategory{
  NORMALMESON,
  EXT_FIRSTDEV_MESON,
  EXT_SECONDDEV_SYM_MESON,
  EXT_SECONDDEV_ANTISYM_MESON,
  EXT_SECONDDEV_DIAG_MESON,
  MIXING
};


/*------------------------------------------------------------------------*/
/* Related to extended mesons constructed from local fields*/
/*------------------------------------------------------------------------*/
/*output index(average over all polarisations)*/
enum WExtMesonBEOutputName{
 BE_pionxB=0,
 BE_rhoxB_T1,
 NUM_WEXTMESON_BE_OUTPUT
};

/*state IDs(with polarisation)*/
enum WExtMesonBEState{
  BE_MS_pionxB_x=0, BE_MS_pionxB_y, BE_MS_pionxB_z,
  BE_MS_rhoxB_T1_x, BE_MS_rhoxB_T1_y, BE_MS_rhoxB_T1_z,
  NUM_WEXTMESON_BE_STATES
};


/*meson operator*/
enum WExtMesonBEOp{
  BE_MO_pionxB_x=0,BE_MO_pionxB_y, BE_MO_pionxB_z,
  BE_MO_rhoxB_T1_x, BE_MO_rhoxB_T1_y, BE_MO_rhoxB_T1_z,
  NUM_WEXTMESON_BE_OPS
};

/*extmeson category*/

enum WExtMesonBECategory{
  ELEC_HYBRID_BE=0,
  MAG_HYBRID_BE,
  MIXING_BE
};

enum FieldTensorId{
  /*used in state table*/
  /*group 1*/
  FB1=0,FB2,FB3,
  FE1,FE2,FE3,
  NUM_FLDS,
  FUNIT,  /*for mixing*/
  
  /*SUM operators*/
  /*group 2*/
  SUM_MAGN,
  SUM_ELEC,
  
  /*group 3*/
  SUM_MAGN_ELEC,

  /**/
  NUM_FLD_OPS
  
};
  

