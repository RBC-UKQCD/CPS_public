#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:05 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/scu_nos/bsm.C,v 1.3 2004-06-04 21:14:05 chulwoo Exp $
//  $Id: bsm.C,v 1.3 2004-06-04 21:14:05 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/scu_nos/bsm.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/****************************************************************
*  This routing allows for a simple c level control of the 	*
*  serial wires on the nga.  The function simply requires a	*
*  base address, a blocklength, stride, number of blocks, and a	*
*  wire number (0-7).  In addition the final argument is either	*
*  set to TRANSMIT or RECEIVE.  Once both sides are set up 	*
*  transmission automatically proceeds.				*
*****************************************************************/

CPS_END_NAMESPACE
#include<comms/scu.h>
#include<comms/nga_reg.h>
CPS_START_NAMESPACE

void bsm(IFloat* bad, int blklen, int std, int nblk, int wire, int xr)
{
  register int 	scu_base = DSP_SCU_BASE;

  while( *( (int*)(scu_base+0x10+wire) )  );

  *((int*)(scu_base+0x50+wire))=(nblk<<22)|(blklen<<12)|(std&0x0fff);
  *( (int*)(scu_base+0x48+wire) ) = (std>>12 );

  if(xr==TRANSMIT)
    *( (IFloat**)(scu_base+0x8+wire) ) = bad;
  else
    *( (IFloat**)(scu_base+0x0+wire) ) = bad;
}

CPS_END_NAMESPACE
