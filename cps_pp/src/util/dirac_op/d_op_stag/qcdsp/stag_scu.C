#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_stag/qcdsp/stag_scu.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: stag_scu.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:12:45  anj
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
//  Revision 1.2  2001/05/25 06:16:06  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: stag_scu.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_stag/qcdsp/stag_scu.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//-------------------------------------------------------------------
//    11/10/97  RDM
//
//  This is a copy of $top/util/dirac_op/d_op_wilson_opt/wilson_scu.C,
//  modified minimally for staggered fermions (names changed to avoid
//  conflict with Wilson).  The following comment is from the Wilson
//  version.
//   
//  A simple try at system calls for the wilson Dirac operator.
//  Likely slow, but extern "C" is easily callable from
//  assembly.
//
//  
//-------------------------------------------------------------------


CPS_END_NAMESPACE
#include <sysfunc.h>
#include <scu_dir_arg.h>
CPS_START_NAMESPACE

//-------------------------------------------------------------------
//  Permanent storage for arguments to C++ system calls.
//-------------------------------------------------------------------

SCUDirArg sarg[16];
SCUDirArg *sargpF[8];
SCUDirArg *sargpB[8];

//-------------------------------------------------------------------
//  This order is set to match the order in wfm_scu_init.asm
//-------------------------------------------------------------------

const SCUDir dir[] = {
  SCU_TP, SCU_TM, SCU_XP, SCU_XM, SCU_YP, SCU_YM, SCU_ZP, SCU_ZM,
  SCU_TP, SCU_TM, SCU_XP, SCU_XM, SCU_YP, SCU_YM, SCU_ZP, SCU_ZM,
};

const SCUXR xr[] = {
  SCU_SEND, SCU_REC, SCU_SEND, SCU_REC,// Used in forward communication
  SCU_SEND, SCU_REC, SCU_SEND, SCU_REC,
  SCU_REC, SCU_SEND, SCU_REC, SCU_SEND,// Used in backward communication
  SCU_REC, SCU_SEND, SCU_REC, SCU_SEND
};

//-------------------------------------------------------------------
//  We first set up forward and backward communication according to
//  Dong Chen's original conventions.  Forward means a send in the
//  positive direction, which is an even wire in Dong's convention.
//  Call this from dirac_init.  Only needs to be called once.
//-------------------------------------------------------------------

extern "C" void StagSCUSetup()
{
  
  int i;
  
  for ( i = 0; i < 16; i++ ) {

    sarg[i].Init( ( void * ) 0 , dir[i], xr[i], 0 );

    if ( i < 8 ) 
      sargpF[i] = &( sarg[i] );
    else 
      sargpB[i-8] = &( sarg[i] );

  }

}

//-------------------------------------------------------------------
//  This is really sloppy, but the staggered code already has
//  packed the DMA values.  We unpack and repack to use
//  system calls available for 5.1.1.  The pointer passed
//  in as an argument points to 16 values.  The first 8 are
//  the DMA values for wires 0 to 7.  The second 8 are the
//  addresses.
//-------------------------------------------------------------------

extern "C" void StagSCUCommForward( unsigned int *dma )
{
  
  int i;
  int j;

  for ( i = 0, j = 8; i < 8; i++, j++ )
    sarg[i].Reload( ( void * ) dma[j], ( dma[i] >> 12 ) & 0x3ff,
      ( dma[i] >> 22 ) & 0x3ff, dma[i] & 0xfff );

  SCUTrans ( sargpF, 8 );

}


//-------------------------------------------------------------------
//  This is really sloppy, but the staggered code already has
//  packed the DMA values.  We unpack and repack to use
//  system calls available for 5.1.1.  The pointer passed
//  in as an argument points to 16 values.  The first 8 are
//  the DMA values for wires 0 to 7.  The second 8 are the
//  addresses.
//-------------------------------------------------------------------

extern "C" void StagSCUCommBackward( unsigned int *dma )
{
  
  int i;
  int j;

  for ( i = 0, j = 8; i < 8; i++, j++ )
    sarg[j].Reload( ( void * ) dma[j], ( dma[i] >> 12 ) & 0x3ff,
      ( dma[i] >> 22 ) & 0x3ff, dma[i] & 0xfff );

  SCUTrans ( sargpB, 8 );
}

 
//-------------------------------------------------------------------
//-------------------------------------------------------------------

extern "C" void StagSCUComplete()
{
  SCUTransComplete();
}

CPS_END_NAMESPACE
