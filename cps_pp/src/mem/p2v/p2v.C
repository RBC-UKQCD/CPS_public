#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:05 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/mem/p2v/p2v.C,v 1.3 2004-06-04 21:14:05 chulwoo Exp $
//  $Id: p2v.C,v 1.3 2004-06-04 21:14:05 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/mem/p2v/p2v.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// p2v.C

CPS_END_NAMESPACE
#include <mem/p2v.h>
#include <string.h>		//memmove()
#include <util/verbose.h>
#include <util/error.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------//
// The following are defined in mem/include/link_p2v.lcf
//------------------------------------------------------------------//

// Vectror utilities library related
//------------------------------------------------------------------//
extern unsigned int vector_start;    // Physical starting address 
extern unsigned int vector_end;      // Physical ending address 
extern unsigned int vector_dest;     // Virtual starting address

// Wilson library related
//------------------------------------------------------------------//
extern unsigned int wfm0_start;    // Physical starting address 
extern unsigned int wfm0_end;      // Physical ending address 
extern unsigned int wfm0_dest;     // Virtual starting address
extern unsigned int wfm1_start;    // Physical starting address 
extern unsigned int wfm1_end;      // Physical ending address 
extern unsigned int wfm1_dest;     // Virtual starting address

// Staggered library related
//------------------------------------------------------------------//
extern unsigned int stag_ds_start;    // Physical starting address 
extern unsigned int stag_ds_end;      // Physical ending address 
extern unsigned int stag_ds_dest;     // Virtual starting address

// Gauge heat bath related
//---------------------------------------------------------------//
extern unsigned int ghb_start;    // Physical starting address 
extern unsigned int ghb_end;      // Physical ending address 
extern unsigned int ghb_dest;     // Virtual starting address


//------------------------------------------------------------------//
// Copies into CRAM the Vector utilities assembly code
//------------------------------------------------------------------//
void p2vVector(){
#if _TARTAN
  char *cname = " ";
  char *fname = "p2vVector()";
  VRB.Flow(cname, fname,
	   "Moving vector_util to CRAM, dest=%x, start= %x, end = %x, len = %x\n",
	   &vector_dest, &vector_start, &vector_end, &vector_end - &vector_start);
  memmove(&vector_dest, &vector_start, 
	  &vector_end - &vector_start);
#endif
}

//------------------------------------------------------------------//
// Copies into CRAM the Wilson library 
//------------------------------------------------------------------//
void p2vWilsonLib(){
#if _TARTAN
  char *cname = " ";
  char *fname = "p2vWilsoLib()";
  VRB.Flow(cname, fname,
	   "Moving wfm0 to CRAM0, dest=%x, start= %x, end = %x, len = %x\n",
	   &wfm0_dest, &wfm0_start, &wfm0_end, &wfm0_end - &wfm0_start);
  memmove(&wfm0_dest, &wfm0_start, 
	  &wfm0_end - &wfm0_start);
  
  VRB.Flow(cname, fname,
	   "Moving wfm1 to CRAM1, dest=%x, start= %x, end = %x, len = %x\n",
	   &wfm1_dest, &wfm1_start, &wfm1_end, &wfm1_end - &wfm1_start);
  memmove(&wfm1_dest, &wfm1_start, 
	  &wfm1_end - &wfm1_start);
#endif
}


//------------------------------------------------------------------//
// Copies into CRAM the Staggered dirac_serial assembly code
//------------------------------------------------------------------//
void p2vStagDs(){
#if _TARTAN
  char *cname = " ";
  char *fname = "p2vStagDs()";
  VRB.Flow(cname, fname,
	   "Moving stag_ds to CRAM, dest=%x, start= %x, end = %x, len = %x\n",
	   &stag_ds_dest, &stag_ds_start, &stag_ds_end, &stag_ds_end - &stag_ds_start);
  memmove(&stag_ds_dest, &stag_ds_start, 
	  &stag_ds_end - &stag_ds_start);
#endif
}


//---------------------------------------------------------------//
// Copies into CRAM the Gauge Heat Bath assembly code
//---------------------------------------------------------------//
void p2vGhb(){
#if _TARTAN
  char *cname = " ";
  char *fname = "p2vGhb()";

VRB.Debug(
"Moving ghb to CRAM, dest=%x, start= %x, end = %x, len = %x\n",
&ghb_dest, &ghb_start, &ghb_end, &ghb_end - &ghb_start);

  VRB.Flow(cname, fname,
   "Moving ghb to CRAM, dest=%x, start= %x, end = %x, len = %x\n",
   &ghb_dest, &ghb_start, &ghb_end, &ghb_end - &ghb_start);
  memmove(&ghb_dest, &ghb_start, 
	  &ghb_end - &ghb_start);

VRB.Debug("Move complete\n");
#endif
}






//------------------------------------------------------------------//
// Copies into CRAM the Clover Matrix Multiplication assembly code  
//------------------------------------------------------------------//
void p2vCloverLib() 
{
#ifdef _TARTAN
  char *cname = " ";
  char *fname = "p2vCloverLib()";

  extern unsigned int cfm_start;    // Physical starting address 
  extern unsigned int cfm_end;      // Physical ending address 
  extern unsigned int cfm_dest;     // Virtual starting address

  p2vWilsonLib();
  

  VRB.Debug("Moving cfm to CRAM, dest=%x, start= %x, end = %x, len = %x\n",
	    &cfm_dest, &cfm_start, &cfm_end, &cfm_end - &cfm_start);
  VRB.Flow(cname, fname,
	   "Moving cfm to CRAM, dest=%x, start= %x, end = %x, len = %x\n",
	   &cfm_dest, &cfm_start, &cfm_end, &cfm_end - &cfm_start);
  memmove(&cfm_dest, &cfm_start, 
	  &cfm_end - &cfm_start);
  VRB.Debug("Move complete\n");

#endif
} 



CPS_END_NAMESPACE
