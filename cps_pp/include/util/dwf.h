#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/dwf.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: dwf.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:29  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:16  anj
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
//  $RCSfile: dwf.h,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/dwf.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// dwf.h
//
// C header file for the dwf fermion library
//
//------------------------------------------------------------------

#ifndef INCLUDED_DWF_H
#define INCLUDED_DWF_H

CPS_END_NAMESPACE
#include <util/wilson.h>
#include <util/vector.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// External declarations
//------------------------------------------------------------------

extern int dwfso_wire_map[];
// Set in dwf_int. For a given index 0-1 corresponding to
// S+, S-  it gives the corresponding wire.


//------------------------------------------------------------------
// Type definitions
//------------------------------------------------------------------
// The Dwf structure typedef
typedef struct{
  IFloat dwf_kappa;    // 1/[2*(5-dwf_height)]
  int vol_4d;         // The 4-dimensional volume   
  int ls;             // The extent of the 5th direction
  IFloat *frm_tmp1;    // Pointer to temorary fermion field 1
  IFloat *frm_tmp2;    // Pointer to temorary fermion field 2
  IFloat *comm_buf;    // Communications buffer (12 words for dwfso)
  Wilson *wilson_p;   // Pointer to the wilson structure
} Dwf;


//------------------------------------------------------------------
// Function declarations
//------------------------------------------------------------------


//------------------------------------------------------------------
// This routine performs all initializations needed before dwf 
// library funcs are called. It sets the addressing related arrays 
// and reserves memory for the needed temporary buffers. It only 
// needs to be called only once at the begining of the program 
// (or after a dwf_end call)before any number of calls to dwf 
// functions are made.
//------------------------------------------------------------------
void dwf_init(Dwf *dwf_p);             // pointer to Dwf struct.


//------------------------------------------------------------------
// This routine frees any memory reserved by dwf_init
//------------------------------------------------------------------
void dwf_end(Dwf *dwf_p);              // pointer to Dwf struct.


//------------------------------------------------------------------
// dwf_mdagm M^dag M where M is the fermion matrix.
// The in, out fields are defined on the checkerboard lattice.
// <out, in> = <dwf_mdagm*in, in> is returned in dot_prd.
//------------------------------------------------------------------
void  dwf_mdagm(Vector *out, 
		Matrix *gauge_field, 
		Vector *in, 
		Float *dot_prd,
		Float mass,
		Dwf *dwf_lib_arg);


//------------------------------------------------------------------
// dwf_dslash is the derivative part of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
// cb = 0/1 <--> even/odd checkerboard of in field.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//------------------------------------------------------------------
void dwf_dslash(Vector *out, 
		Matrix *gauge_field, 
		Vector *in,
		Float mass,
		int cb,
		int dag,
		Dwf *dwf_lib_arg);


//------------------------------------------------------------------
// dwf_m is the fermion matrix.  
// The in, out fields are defined on the checkerboard lattice.
//------------------------------------------------------------------
void  dwf_m(Vector *out, 
	    Matrix *gauge_field, 
	    Vector *in,
	    Float mass,
	    Dwf *dwf_lib_arg);


//------------------------------------------------------------------
// dwf_mdag is the dagger of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
//------------------------------------------------------------------
void dwf_mdag(Vector *out, 
	      Matrix *gauge_field, 
	      Vector *in,
	      Float mass,
	      Dwf *dwf_lib_arg);



//------------------------------------------------------------------
// dwf_dslash_4 is the derivative part of the space-time part of
// the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
// cb = 0/1 <--> even/odd checkerboard of in field.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//------------------------------------------------------------------
void dwf_dslash_4(Vector *out, 
		  Matrix *gauge_field, 
		  Vector *in,
		  int cb,
		  int dag,
		  Dwf *dwf_lib_arg);



//------------------------------------------------------------------
// dwf_dslash_5_plus is the derivative part of the 5th direction
// part of the fermion matrix. This routine accumulates the result
// on the out field 
// The in, out fields are defined on the checkerboard lattice.
// The action of this operator is the same for even/odd
// checkerboard fields because there is no gauge field along
// the 5th direction.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//------------------------------------------------------------------
void dwf_dslash_5_plus(Vector *out, 
		       Vector *in,
		       Float mass,
		       int dag,
		       Dwf *dwf_lib_arg);







#endif

CPS_END_NAMESPACE
