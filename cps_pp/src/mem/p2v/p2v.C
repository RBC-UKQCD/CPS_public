#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:54 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/mem/p2v/p2v.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Id: p2v.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/07/03 17:00:48  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.3  2001/06/21 15:40:11  anj
//  Updated the _TARTAN ifdefs, using PARALLEL instead (where appropriate).Anj
//
//  Revision 1.2  2001/06/19 18:11:41  anj
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
//  Revision 1.2  2001/05/25 06:16:01  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: p2v.C,v $
//  $Revision: 1.2 $
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
