#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/nga_reg.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: nga_reg.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:12:13  anj
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
//  Revision 1.2  2001/05/25 06:16:03  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: nga_reg.h,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/nga_reg.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*
 *  nga_reg.h
 */

#ifndef INCLUDED_NGA_REG_H
#define INCLUDED_NGA_REG_H


/* for circular buffer */
#ifdef _TARTAN
enum{ BANK1_BASE = 0x900000, BANK2_BASE = 0xa00000,
      BANK3_BASE = 0xb00000, BANK4_BASE = 0xc00000,
      BANK_SIZE = 0x80000 };
#else
enum{ BANK1_BASE = 0, BANK2_BASE = 0,
      BANK3_BASE = 0, BANK4_BASE = 0,
      BANK_SIZE = 0 };
#endif

enum{ DSP_SCU_BASE = 0x813000 };

enum{ CBUF_CNTRL_REG_BASE = 0x815800 };

enum { NODE_ID_ADDR = 0x7ff00 };

// GRF: CRAM_SCRATCH_SIZE should be at least 72 for Staple to work
// and 90 for Rectangle to work.
enum{ CRAM_SCRATCH_SIZE = 96 };
#ifdef _TARTAN
enum{ CRAM_SCRATCH_ADDR = 0x809F95 - CRAM_SCRATCH_SIZE };
#else
// initialized in ../cbuf_sim/cram.C where memory is allocated
// in order to imitate the scratch cram.
extern unsigned int CRAM_SCRATCH_ADDR ; 
#endif

#endif


CPS_END_NAMESPACE
