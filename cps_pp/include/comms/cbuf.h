#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/cbuf.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: cbuf.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:12:12  anj
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
//  $RCSfile: cbuf.h,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/cbuf.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// cbuf.h
//
// Various routines that relate to the circular buffer.
//
//------------------------------------------------------------------

#ifndef INCLUDED_CBUF_H
#define INCLUDED_CBUF_H


//------------------------------------------------------------------
// saveCbufCntrlReg:
// Saves the contents of the circular buffer control registers
// in the cbuf_CntrlReg_reg static array
//------------------------------------------------------------------
void saveCbufCntrlReg(void);


//------------------------------------------------------------------
// restoreCbufCntrlReg:
// Restores the contents of the circular buffer control registers
// from the cbuf_CntrlReg_reg static array
//------------------------------------------------------------------
void restoreCbufCntrlReg(void);


//------------------------------------------------------------------
// setCbufCntrlReg(reg_no, value):
// Sets the contents of the circular buffer control register 
// reg_no to value.
//------------------------------------------------------------------
void setCbufCntrlReg(int reg_no, unsigned int value);


//------------------------------------------------------------------
// printCbufCntrlReg:
// Prints the contents of the circular buffer control registers
// under the control of VRB.Flow
//------------------------------------------------------------------
void printCbufCntrlReg(void);


#endif







CPS_END_NAMESPACE
