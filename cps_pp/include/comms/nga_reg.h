#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008-09-18 14:42:27 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/nga_reg.h,v 1.7 2008-09-18 14:42:27 chulwoo Exp $
//  $Id: nga_reg.h,v 1.7 2008-09-18 14:42:27 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.7 $
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
extern unsigned long CRAM_SCRATCH_ADDR ; 
#endif

#endif


CPS_END_NAMESPACE
