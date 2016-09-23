#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/08/18 11:57:36 $
//  $Header: /space/cvs/cps/cps++/include/comms/cbuf.h,v 1.4 2004/08/18 11:57:36 zs Exp $
//  $Id: cbuf.h,v 1.4 2004/08/18 11:57:36 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /space/cvs/cps/cps++/include/comms/cbuf.h,v $
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
