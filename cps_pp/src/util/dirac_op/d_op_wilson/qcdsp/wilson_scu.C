#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:08 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp/wilson_scu.C,v 1.3 2004-06-04 21:14:08 chulwoo Exp $
//  $Id: wilson_scu.C,v 1.3 2004-06-04 21:14:08 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: wilson_scu.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp/wilson_scu.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//-------------------------------------------------------------------
//    11/4/97  RDM
//
//
//  A simple try at system calls for the wilson Dirac operator.
//  Likely slow, but extern "C" is easily callable from
//  assembly.
//-------------------------------------------------------------------


CPS_END_NAMESPACE
#include <stdio.h>
#include <comms/sysfunc.h>
#include <scu_dir_arg.h>
#include<util/gjp.h>
CPS_START_NAMESPACE

//-------------------------------------------------------------------
//  Permanent storage for arguments to C++ system calls.
//-------------------------------------------------------------------

SCUDirArg arg[16];
SCUDirArg *argpF[8];
SCUDirArg *argpB[8];

//-------------------------------------------------------------------
//  This order is set to match the order in wfm_scu_init.asm
//-------------------------------------------------------------------

const SCUDir dir[] = {
  SCU_XP, SCU_XM, SCU_YP, SCU_YM, SCU_ZP, SCU_ZM, SCU_TP, SCU_TM,
  SCU_XP, SCU_XM, SCU_YP, SCU_YM, SCU_ZP, SCU_ZM, SCU_TP, SCU_TM
};

const SCUXR xr[] = {
  SCU_SEND, SCU_REC, SCU_SEND, SCU_REC,
  SCU_SEND, SCU_REC, SCU_SEND, SCU_REC,
  SCU_REC, SCU_SEND, SCU_REC, SCU_SEND,
  SCU_REC, SCU_SEND, SCU_REC, SCU_SEND
};

//-------------------------------------------------------------------
//  This is really sloppy, but the Wilson code already has
//  packed the DMA values.  We unpack and repack to use
//  system calls available for 5.1.1
//-------------------------------------------------------------------

extern "C" void WilsonSCUSetDMA( unsigned int X, unsigned int Y, 
				unsigned int Z, unsigned int T )
{
  int axis;
  
  unsigned int ui[16];

  ui[0] = X;
  ui[1] = X;
  ui[2] = Y;
  ui[3] = Y;
  ui[4] = Z;
  ui[5] = Z;
  ui[6] = T;
  ui[7] = T;
  ui[8] = X;
  ui[9] = X;
  ui[10] = Y;
  ui[11] = Y;
  ui[12] = Z;
  ui[13] = Z;
  ui[14] = T;
  ui[15] = T;

  int tot = 0;
  for(int k=0; k<4; k++){
    if(gjp_local_axis[k] != 1) tot = tot + 1;
  }
  tot = 2*tot;


  int j = 0;

  for (int i = 0; i < 16; i++ ) {
    
    axis = (i%8) / 2;

    if(gjp_local_axis[axis] != 1) {

      arg[j].Init( ( void * ) 0 , dir[i], xr[i],
		   ( ui[i] >> 12 ) & 0x3ff, ( ui[i] >> 22 ) & 0x3ff, ui[i] & 0xfff );
      
      if ( j < tot ) 
	argpF[j] = &( arg[j] );
      else 
	argpB[j-tot] = &( arg[j] );

      j++;
    }

  }

  SCUSetDMA ( argpF, j/2 );

}


//-------------------------------------------------------------------
//-------------------------------------------------------------------

extern "C" void WilsonSCUCommForward(
  void * XPaddr, void * YPaddr, void * ZPaddr, void * TPaddr, 
  void * XMaddr, void * YMaddr, void * ZMaddr, void * TMaddr  )
{
  int axis;
  void *addr[8];
  addr[0] = XPaddr;
  addr[1] = XMaddr;
  addr[2] = YPaddr;
  addr[3] = YMaddr;
  addr[4] = ZPaddr;
  addr[5] = ZMaddr;
  addr[6] = TPaddr;
  addr[7] = TMaddr;

  int j = 0;
  for(int i=0; i<8; i++){
    axis = i/2;

    if(gjp_local_axis[axis] != 1) {
      arg[j].Addr( addr[i] );
      j++;
    }

  }

  SCUTransAddr( argpF, j );
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------

extern "C" void WilsonSCUCommBackward(
  void * XMaddr, void * YMaddr, void * ZMaddr, void * TMaddr,
  void * XPaddr, void * YPaddr, void * ZPaddr, void * TPaddr )
{
  int axis;
  void *addr[8];
  addr[0] = XPaddr;
  addr[1] = XMaddr;
  addr[2] = YPaddr;
  addr[3] = YMaddr;
  addr[4] = ZPaddr;
  addr[5] = ZMaddr;
  addr[6] = TPaddr;
  addr[7] = TMaddr;

  int tot = 0;
  for(int k=0; k<4; k++){
    if(gjp_local_axis[k] != 1) tot = tot + 1;
  }
  tot = 2*tot;

  int j =0;
  for(int i=0; i<8; i++){
    axis = i/2;

    if(gjp_local_axis[axis] != 1) {
      arg[j+tot].Addr( addr[i] );
      j++;
    }

  }

  SCUTransAddr( argpB, j );
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------

extern "C" void WilsonSCUComplete()
{
  SCUTransComplete();
}

CPS_END_NAMESPACE
