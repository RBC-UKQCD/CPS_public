#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/cbuf/cbuf.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: cbuf.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:11:42  anj
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
//  $RCSfile: cbuf.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/cbuf/cbuf.C,v $
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
#include<comms/nga_reg.h>
#include<util/verbose.h>
CPS_START_NAMESPACE

//--------------------------------------------------------------------------
// static array to save contents of circular buffer
//--------------------------------------------------------------------------
static unsigned cbuf_cntrl_reg[5];


//--------------------------------------------------------------------------
// saveCbufCntrlReg:
// Saves the contents of the circular buffer control registers
// in the cbuf_cntrl_reg static array
//--------------------------------------------------------------------------
void saveCbufCntrlReg(void) {
  int i;
  unsigned int *ptr;
  ptr = (unsigned *) CBUF_CNTRL_REG_BASE;

  for(i=0; i<5; i++){
    cbuf_cntrl_reg[i] = ptr[i];
  }

}


//--------------------------------------------------------------------------
// restoreCbufCntrlReg:
// Restores the contents of the circular buffer control registers
// from the cbuf_cntrl_reg static array
//--------------------------------------------------------------------------
void restoreCbufCntrlReg(void){
  int i;
  unsigned int *ptr;
  ptr = (unsigned *) CBUF_CNTRL_REG_BASE;

  for(i=0; i<5; i++){
    ptr[i] = cbuf_cntrl_reg[i];
  }

}


//------------------------------------------------------------------
// setCbufCntrlReg(reg_no, value):
// Sets the contents of the circular buffer control register 
// reg_no to value.
//------------------------------------------------------------------
void setCbufCntrlReg(int reg_no, unsigned int value){
  unsigned int *ptr;
  ptr = (unsigned *) CBUF_CNTRL_REG_BASE;
  ptr[reg_no] = value;
}


//------------------------------------------------------------------
// printCbufCntrlReg:
// Prints the contents of the circular buffer control registers
// under the control of VRB.Flow
//------------------------------------------------------------------
void printCbufCntrlReg(void){
  int i;
  unsigned int *ptr;
  ptr = (unsigned *) CBUF_CNTRL_REG_BASE;

  for(i=0; i<5; i++){
    VRB.Flow(" ","printCbufCntrlReg()",	"reg[%d], %x\n", i, ptr[i]);
  }
}
CPS_END_NAMESPACE
