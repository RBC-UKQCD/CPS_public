#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:46:30 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/mobius.h,v 1.2 2013-04-05 17:46:30 chulwoo Exp $
//  $Id: mobius.h,v 1.2 2013-04-05 17:46:30 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: mobius.h,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/mobius.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// mobius.h
//
// C header file for the mobius fermion library
//
//------------------------------------------------------------------

#ifndef INCLUDED_MOBIUS_H
#define INCLUDED_MOBIUS_H

CPS_END_NAMESPACE
#include <util/wilson.h>
#include <util/vector.h>
#include <util/dwf.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// External declarations
//------------------------------------------------------------------

extern int mobiusso_wire_map[];
// Set in mobius_int. For a given index 0-1 corresponding to
// S+, S-  it gives the corresponding wire.


//------------------------------------------------------------------
// Type definitions
//------------------------------------------------------------------
// The Mobius structure typedef
#if 0 // use Dwf instead
typedef struct{
  IFloat mobius_kappa_b;    // 1/[2*(b*(4-dwf_height)+1)]
  IFloat mobius_kappa_c;    // 1/[2*(c*(4-dwf_height)-1)]
  int vol_4d;         // The 4-dimensional volume   
  int ls;             // The extent of the 5th direction
  IFloat *frm_tmp1;    // Pointer to temorary fermion field 1
  IFloat *frm_tmp2;    // Pointer to temorary fermion field 2

  // doesn't need any more
  //IFloat *frm_tmp3;    // Pointer to temorary fermion field 3
  
  IFloat *comm_buf;    // Communications buffer (12 words for mobiusso)
  Wilson *wilson_p;   // Pointer to the wilson structure
#if TARGET == QCDOC
  SCUDirArgIR *PlusArg[2];
  SCUDirArgIR *MinusArg[2];
  SCUDirArgMulti *Plus;
  SCUDirArgMulti *Minus;
#endif
} Mobius;
#endif

//------------------------------------------------------------------
// Function declarations
//------------------------------------------------------------------


//------------------------------------------------------------------
// This routine performs all initializations needed before mobius 
// library funcs are called. It sets the addressing related arrays 
// and reserves memory for the needed temporary buffers. It only 
// needs to be called only once at the begining of the program 
// (or after a mobius_end call)before any number of calls to mobius 
// functions are made.
//------------------------------------------------------------------
//void mobius_init(Dwf *mobius_p);             // pointer to Dwf struct.


//------------------------------------------------------------------
// This routine frees any memory reserved by mobius_init
//------------------------------------------------------------------
//void mobius_end(Dwf *mobius_p);              // pointer to Dwf struct.


//------------------------------------------------------------------
// mobius_mdagm M^dag M where M is the fermion matrix.
// The in, out fields are defined on the checkerboard lattice.
// <out, in> = <mobius_mdagm*in, in> is returned in dot_prd.
//------------------------------------------------------------------
void  mobius_mdagm(Vector *out, 
		   Matrix *gauge_field, 
		   Vector *in, 
		   Float *dot_prd,
		   Float mass,
		   Dwf *mobius_lib_arg);


// spectrum shifted version : (H-mu)(H-mu)
void    mobius_mdagm_shift(Vector *out, 
			   Matrix *gauge_field, 
			   Vector *in, 
			   Float *dot_prd,
			   Float mass,
			   Dwf *mobius_lib_arg,
			   Float shift);

//------------------------------------------------------------------
// mobius_dslash is the derivative part of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
// cb = 0/1 <--> even/odd checkerboard of in field.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//------------------------------------------------------------------
void mobius_dslash(Vector *out, 
		   Matrix *gauge_field, 
		   Vector *in,
		   Float mass,
		   int cb,
		   int dag,
		   Dwf *mobius_lib_arg);


//------------------------------------------------------------------
// mobius_m is the fermion matrix.  
// The in, out fields are defined on the checkerboard lattice.
//------------------------------------------------------------------
void  mobius_m(Vector *out, 
	       Matrix *gauge_field, 
	       Vector *in,
	       Float mass,
	       Dwf *mobius_lib_arg);
void  mobius_mdagm(Vector *out, 
		   Matrix *gauge_field, 
		   Vector *in,
		   Float mass,
		   Dwf *mobius_lib_arg);


//------------------------------------------------------------------
// mobius_mdag is the dagger of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
//------------------------------------------------------------------
void mobius_mdag(Vector *out, 
		 Matrix *gauge_field, 
		 Vector *in,
		 Float mass,
		 Dwf *mobius_lib_arg);


//------------------------------------------------------------------
// mobius_dslash_4 is the derivative part of the space-time part of
// the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
// cb = 0/1 <--> even/odd checkerboard of in field.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//------------------------------------------------------------------
void mobius_dslash_4(Vector *out, 
		     Matrix *gauge_field, 
		     Vector *in,
		     int cb,
		     int dag,
		     Dwf *mobius_lib_arg, Float mass);

void mobius_m5inv(Vector *out,
	       Vector *in, 
	       Float mass,
	       int dag,
	       Dwf *mobius_lib_arg);

// in-place version
void mobius_m5inv(Vector *inout, 
	       Float mass,
	       int dag,
	       Dwf *mobius_lib_arg);

void mobius_kappa_dslash_5_plus(Vector *out, 
				Vector *in, 
				Float mass,
				int dag, 
				Dwf *mobius_lib_arg,
				Float fact);
void mobius_dslash_5_plus(Vector *out, 
			  Vector *in, 
			  Float mass,
			  int dag, 
			  Dwf *mobius_lib_arg);
void mobius_kappa_dslash_5_plus_dag0(Vector *out,
                       Vector *in,
                       Float mass,
                       Dwf *mobius_lib_arg,
                       Float a_five_inv );
void mobius_kappa_dslash_5_plus_dag1(Vector *out,
                       Vector *in,
                       Float mass,
                       Dwf *mobius_lib_arg,
                       Float a_five_inv );

void ReflectAndMultGamma5( Vector *out, const Vector *in,  int nodevol, int ls);

void cPRL_plus(Vector *out, Vector *in, int dag, Float mass, Dwf *mobius_lib_arg);

void mobius_dminus(Vector *out, 
		   Matrix *gauge_field, 
		   Vector *in, 
		   int cb, 
		   int dag, 
		   Dwf *mobius_lib_arg);

#endif

CPS_END_NAMESPACE
