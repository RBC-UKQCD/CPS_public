#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2006-11-25 19:10:04 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qmp/cbuf/cbuf.C,v 1.2 2006-11-25 19:10:04 chulwoo Exp $
//  $Id: cbuf.C,v 1.2 2006-11-25 19:10:04 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qmp/cbuf/cbuf.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// 10/31/97
//
// cbuf:
//
// Routines that control the circular buffer.
//
//--------------------------------------------------------------------------

CPS_END_NAMESPACE
#include<comms/cbuf.h>
CPS_START_NAMESPACE


//--------------------------------------------------------------------------
// saveCbufCntrlReg:
// Saves the contents of the circular buffer control registers
// in the cbuf_cntrl_reg static array
//--------------------------------------------------------------------------
void saveCbufCntrlReg(void) {}


//--------------------------------------------------------------------------
// restoreCbufCntrlReg:
// Restores the contents of the circular buffer control registers
// from the cbuf_cntrl_reg static array
//--------------------------------------------------------------------------
void restoreCbufCntrlReg(void){}


//------------------------------------------------------------------
// setCbufCntrlReg(reg_no, value):
// Sets the contents of the circular buffer control register 
// reg_no to value.
//------------------------------------------------------------------
void setCbufCntrlReg(int reg_no, unsigned int value){}


//------------------------------------------------------------------
// printCbufCntrlReg:
// Prints the contents of the circular buffer control registers
// under the control of VRB.Flow
//------------------------------------------------------------------
void printCbufCntrlReg(void){}
CPS_END_NAMESPACE
