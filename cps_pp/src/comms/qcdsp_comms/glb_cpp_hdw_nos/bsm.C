#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/glb_cpp_hdw_nos/bsm.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: bsm.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:49:55  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:12:04  anj
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
//  Revision 1.2  2001/05/25 06:16:02  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: bsm.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/glb_cpp_hdw_nos/bsm.C,v $
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
