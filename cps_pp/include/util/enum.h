#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Magic numbers.

  $Id: enum.h,v 1.10 2004-09-02 16:52:54 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-09-02 16:52:54 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/enum.h,v 1.10 2004-09-02 16:52:54 zs Exp $
//  $Id: enum.h,v 1.10 2004-09-02 16:52:54 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: enum.h,v $
//  $Revision: 1.10 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/enum.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

#ifndef INCLUDED_ENUM_H
#define INCLUDED_ENUM_H    //!< Prevent multiple inclusion


//------------------------------------------------------------------
//! The directions in the lattice
//------------------------------------------------------------------
enum DirType {DIR_X,
	      DIR_Y,
	      DIR_Z,
	      DIR_T,
	      DIR_S};

//------------------------------------------------------------------
//! The types of fermion action
//------------------------------------------------------------------
enum FclassType {
    F_CLASS_NONE,   
    F_CLASS_STAG,
    F_CLASS_WILSON,
    F_CLASS_CLOVER,
    F_CLASS_DWF,
    F_CLASS_ASQTAD
};


//------------------------------------------------------------------
//! The types of gauge action.
//------------------------------------------------------------------
enum GclassType {
    G_CLASS_NONE,   
    G_CLASS_WILSON,
    G_CLASS_POWER_PLAQ,
    G_CLASS_IMPR_RECT,
    G_CLASS_POWER_RECT,
    G_CLASS_IMPR_OLSYM
};


//------------------------------------------------------------------
//! The lattice storage orders.
//------------------------------------------------------------------
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


//------------------------------------------------------------------
//! The  storage order conversion flags.
//------------------------------------------------------------------
enum CnvFrmType {
    CNV_FRM_NO   = 0,   /*!< Do not convert fermion field */
    CNV_FRM_YES  = 1  /*!< Convert fermion field assuming that it is defined
			on the whole lattice \e i.e. one spinor per site and
			not just on the sites of a single checkerboard
			parity.*/
};

//------------------------------------------------------------------
//!  The dimensionality of the fermion field.
//------------------------------------------------------------------
enum FermionFieldDimension { FOUR_D, FIVE_D };

//------------------------------------------------------------------
//! The kinds of preservation 
//------------------------------------------------------------------
enum PreserveType {PRESERVE_NO  = 0,  /*!< Do not preserve. */
		   PRESERVE_YES = 1 /*!< Preserve. */
};

//------------------------------------------------------------------
//! The kinds of starting configurations.
//
// START_CONF_ORD  Ordered start. After the configuration is set
// the GJP.StartConfKind() is set to START_CONF_MEM.
//
// START_CONF_DISORD  Disordered start. After the configuration
// is set the GJP.StartConfKind() is set to START_CONF_MEM.
//
// START_CONF_FILE Read from file. After the configuration is set
// the GJP.StartConfKind() is set to START_CONF_MEM.
//
// START_CONF_LOAD Memory is not allocated for the gauge field.
// Instead the gauge_field pointer is set to GJP.StartConfLoadAddr().
// After the configuration is set the GJP.StartConfKind() is set 
// to START_CONF_MEM.
//
// START_CONF_MEM Memory is allocated for the gauge field with 
// pmalloc. When the program begins execution it checcks if the
// pmalloc address is the same as the GJP.StartConfLoadAddr(). 
// If not it exits with an error. Subsequent calls to the 
// constructor with START_CONF_MEM do nothing except to set the 
// gauge field to canonical order.
//
// For implementation details see the Lattice constructor in 
// util/lattice/lattice_base/lattice_base.C
//
//------------------------------------------------------------------

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


//------------------------------------------------------------------
//! The possible kinds of initial RNG seeds.
//------------------------------------------------------------------
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

//------------------------------------------------------------------
// The possible types of checker  boards (chequer boards?)
//------------------------------------------------------------------
enum ChkbType {
    CHKB_EVEN =0,   // Even checkerboard
    CHKB_ODD  =1};  // Odd checkerboard


//------------------------------------------------------------------
//! Hermitian conjugate flag
//------------------------------------------------------------------
enum DagType {
    DAG_NO  = 0,   /*!< No hermitian conjugate */
    DAG_YES = 1    /*!< Hermitian conjugate */
};


//------------------------------------------------------------------
//! The kinds of boundary conditions.
//------------------------------------------------------------------
enum BndCndType {
    BND_CND_PRD,    /*!< Periodic */
    BND_CND_APRD    /*!< Antiperiodic */
};


//------------------------------------------------------------------
//! The possible kinds of gauge fixing
//------------------------------------------------------------------
enum FixGaugeType{ 
    FIX_GAUGE_NONE = -2,   /*!< No gauge fixing */
    FIX_GAUGE_LANDAU = -1,  /*!< Fixing to Landau gauge */
    FIX_GAUGE_COULOMB_X = 0,  /*!< Fixing to Coulomb gauge */
    FIX_GAUGE_COULOMB_Y = 1,  /*!< Fixing to Coulomb gauge */
    FIX_GAUGE_COULOMB_Z = 2,  /*!< Fixing to Coulomb gauge */
    FIX_GAUGE_COULOMB_T = 3 /*!< Fixing to Coulomb gauge */
};		   

//------------------------------------------------------------------
//! The possible kinds of spin projections
//------------------------------------------------------------------
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

//------------------------------------------------------------------
//! The possible kinds of Sigma projections
//------------------------------------------------------------------
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

//------------------------------------------------------------------
//! Operators, in terms of the fermion matrix \e M, of which the eigenvalues/vectors can be found.
//------------------------------------------------------------------
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

//------------------------------------------------------------------
//! The types of rational approximation for RHMC
//------------------------------------------------------------------
enum RatApproxType {
  CONSTANT,  /*!< CONSTANT - approximation is held constant*/
  DYNAMIC    /*!< DYNAMIC - approximation recalculated at the end of the trajectory*/
};
		  
//------------------------------------------------------------------
//! The types multishift solve to perform
//------------------------------------------------------------------
enum MultiShiftSolveType {
  SINGLE,  /*!< SINGLE - solutions are summed to a common vector*/
  MULTI    /*!< MULTI - solutions are summed to separate vectors*/
};


#endif


CPS_END_NAMESPACE
