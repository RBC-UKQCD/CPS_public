#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/enum.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: enum.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.3  2002/12/04 17:16:27  zs
//  Merged the new 2^4 RNG into the code.
//  This new RNG is implemented in the LatRanGen class.
//  The following algorithm and utility classes are affected:
//
//  AlgEig                  Fdwf
//  AlgGheatBath            Fstag
//  AlgHmd                  GlobalJobParameter
//  AlgNoise                Lattice
//  AlgPbp                  Matrix
//  AlgThreept              RandomGenerator
//                          Vector
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
//  $RCSfile: enum.h,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/enum.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// enum.h
//
// Header file for enum types.
//
//------------------------------------------------------------------


#ifndef INCLUDED_ENUM_H
#define INCLUDED_ENUM_H


//------------------------------------------------------------------
// The possible kinds of directions
//------------------------------------------------------------------
enum DirType {DIR_X,
	      DIR_Y,
	      DIR_Z,
	      DIR_T,
	      DIR_S};

//------------------------------------------------------------------
// The possible kinds of fermion classes
//------------------------------------------------------------------
enum FclassType {F_CLASS_NONE,   
		 F_CLASS_STAG,
		 F_CLASS_WILSON,
		 F_CLASS_CLOVER,
		 F_CLASS_DWF};


//------------------------------------------------------------------
// The possible kinds of gauge classes
//------------------------------------------------------------------
enum GclassType {G_CLASS_NONE,   
                 G_CLASS_WILSON,
                 G_CLASS_POWER_PLAQ,
                 G_CLASS_IMPR_RECT,
                 G_CLASS_POWER_RECT,
                 G_CLASS_IMPR_OLSYM};


//------------------------------------------------------------------
// The possible kinds of lattice storage orders.
//------------------------------------------------------------------
enum StrOrdType {CANONICAL = 0, // Canonical storage order is
		                // expected by Lattice functions
		STAG      = 1,  // Staggered ferm. storage order is
                                // expected by DiracOpStag functions
	        WILSON    = 2,  // Wilson storage order is
                                // expected by DiracOpWilson 
                                // DiracOpClover and DiracOpDwf 
                                // functions
             G_WILSON_HB  = 3}; // Storage order for the 
                                // wilson gauge action heat bath.
                                // This order is expected by
                                // all QncWilsonHb functions.


//------------------------------------------------------------------
// The possible storage order conversion flags.
//------------------------------------------------------------------
enum CnvFrmType {CNV_FRM_NO   = 0,   // Do not convert fermion field
                 CNV_FRM_YES  = 1};  // Convert fermion field 
                                     // assuming that it is defined 
                                     // on the whole lattice i.e.
                                     // one spinor per site AND NOT
                                     // just on the sites of a 
                                     // single checkerboard.

//------------------------------------------------------------------
//  Divide FsiteSize by SnodeSites in RandGaussVector
//------------------------------------------------------------------
enum FermionFieldDimension { FOUR_D, FIVE_D };


//------------------------------------------------------------------
// The possible kinds of preservation
//------------------------------------------------------------------
enum PreserveType {PRESERVE_NO  = 0,  // Do not preserve.
		   PRESERVE_YES = 1}; // Preserve.


//------------------------------------------------------------------
// The possible kinds of starting configurations.
//
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
// START_CONF_LOAD Mamory is not allocated for the gauge filed.
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
enum StartConfType {START_CONF_ORD     = 0,
		    START_CONF_DISORD  = 1,
		    START_CONF_FILE    = 2,
		    START_CONF_LOAD    = 3,
		    START_CONF_MEM     = 4};


//------------------------------------------------------------------
// The possible kinds of starting seeds.
//
//
// START_SEED_FIXED It does not change between reboots but it is
// different on each node depending on the node coordinates. 
// It is derived from the qos function SeedST(). On node 0,0 it has
// the SeedST() value.
//
// START_SEED_FIXED_UNIFORM It does not change between reboots,
// is the same on all nodes and is equal to the qos SeedST() value.
//
// START_SEED It changes between reboots and it is 
// different on each node depending on the GJP node coordinates. 
// It is derived from the qos function SeedS(). On node 0,0 it has
// the SeedS() value.
//
// START_SEED_UNIFORM It changes between reboots, it is the 
// same on all nodes and is equal to the qos SeedS() value.
//
// START_SEED_INPUT It is different on each node depending on the 
// GJP node coordinates. It is derived from the DoArg 
// start_seed_value which is set by GJP.Initialize().
//
// START_SEED_INPUT_UNIFORM It is the same on each node. It is 
// derived from the DoArg start_seed_value which is set by 
// GJP.Initialize().
//
// START_SEED_INPUT_NODE It is different on each node depending 
// on the physical node coordinates (the ones given by the qos). 
// It is derived from the DoArg start_seed_value which is set by 
// GJP.Initialize(). If this option is used each processor has
// a different seed and therefore the machine can not be divided
// into partitions. This option is usefull only when each node is
// to be used as an individual computer with local communication.
//
// For implementation details see util/gjp/gjp.C
//
//------------------------------------------------------------------
enum StartSeedType {START_SEED_FIXED          = 0,   
		    START_SEED_FIXED_UNIFORM  = 1,
		    START_SEED                = 2,
		    START_SEED_UNIFORM        = 3,
		    START_SEED_INPUT          = 4,
		    START_SEED_INPUT_UNIFORM  = 5,
		    START_SEED_INPUT_NODE     = 6};


//------------------------------------------------------------------
// The possible types of checker boards
//------------------------------------------------------------------
enum ChkbType {CHKB_EVEN =0,   // Even checkerboard
               CHKB_ODD  =1};  // Odd checkerboard


//------------------------------------------------------------------
// The possible types of dagger
//------------------------------------------------------------------
enum DagType {DAG_NO  = 0,   // Do not take the dagger
              DAG_YES = 1};  // Take the dagger


//------------------------------------------------------------------
// The possible kinds of boundary conditions.
//------------------------------------------------------------------
enum BndCndType {BND_CND_PRD, BND_CND_APRD};


//------------------------------------------------------------------
// The possible kinds of gauges for FixGauge(..) routine
//------------------------------------------------------------------
enum FixGaugeType{ FIX_GAUGE_NONE = -2,
                   FIX_GAUGE_LANDAU = -1, 
		   FIX_GAUGE_COULOMB_X = 0, 
	           FIX_GAUGE_COULOMB_Y = 1, 
	           FIX_GAUGE_COULOMB_Z = 2, 
	           FIX_GAUGE_COULOMB_T = 3};


//------------------------------------------------------------------
// The possible kinds of spin projections
//------------------------------------------------------------------
enum SprojType {SPROJ_XM = 0,     // sproj with (1 - gamma_0)
		SPROJ_YM = 1,     // sproj with (1 - gamma_1)
		SPROJ_ZM = 2,     // sproj with (1 - gamma_2)
		SPROJ_TM = 3,     // sproj with (1 - gamma_3)
		SPROJ_XP = 4,     // sproj with (1 + gamma_0)
		SPROJ_YP = 5,     // sproj with (1 + gamma_1)
		SPROJ_ZP = 6,     // sproj with (1 + gamma_2)
		SPROJ_TP = 7 };   // sproj with (1 + gamma_3) 

//------------------------------------------------------------------
// The possible kinds of Sigma projections
//------------------------------------------------------------------
enum SigmaprojType {SIGMAPROJ_XY = 0,     // Sigmaproj with Sigma_{0,1}
		SIGMAPROJ_XZ = 1,     // Sigmaproj with Sigma_{0,2}
		SIGMAPROJ_XT = 2,     // Sigmaproj with Sigma_{0,3}
		SIGMAPROJ_YZ = 3,     // Sigmaproj with Sigma_{1,2}
		SIGMAPROJ_YT = 4,     // Sigmaproj with Sigma_{1,3}
		SIGMAPROJ_YX = 5,     // Sigmaproj with Sigma_{1,0}    
		SIGMAPROJ_ZT = 6,     // Sigmaproj with Sigma_{2,3}
		SIGMAPROJ_ZX = 7,     // Sigmaproj with Sigma_{2,0}
		SIGMAPROJ_ZY = 8,     // Sigmaproj with Sigma_{2,1}  
		SIGMAPROJ_TX = 9,     // Sigmaproj with Sigma_{3,0}
		SIGMAPROJ_TY =10,     // Sigmaproj with Sigma_{3,1}
		SIGMAPROJ_TZ =11 };   // Sigmaproj with Sigma_{3,2} 

//------------------------------------------------------------------
// The type of operators used in RitzEigMat and RitzMat
//------------------------------------------------------------------
enum RitzMatType {NONE,            // No eigenvalues requested
		  MAT_HERM,        // Use hermitian (full) matrix in RitzEig
		  MATPC_HERM,      // Use hermitian preconditioned mat in RitzEig
		  MATPCDAG_MATPC,  // Use preconditioned MatDag*Mat in Ritz
		  MATDAG_MAT,      // Use MatDag*Mat in Ritz
		  NEG_MATDAG_MAT}; // Use neg of MatDag*Mat in Ritz

#endif

CPS_END_NAMESPACE
