#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/glb_hdw/global_sum_handler.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: global_sum_handler.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:12:09  anj
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
//  $RCSfile: global_sum_handler.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/glb_hdw/global_sum_handler.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include<global_sum_handler.h>
#include <stdio.h>                             // printf()
CPS_START_NAMESPACE




inline void GLOBAL_SUM_HANDLER_OVERFLOW(unsigned int err_info);
inline void GLOBAL_SUM_HANDLER_SCU_ERR(unsigned int err_info);
inline void PRINT_PATH();


//--------------------------------------------------------------------------
// union  ScuErrReg
//--------------------------------------------------------------------------

union  ScuErrReg{
  unsigned int _intval;
  struct
  {
    unsigned int recv_errs     :  4 ;
    unsigned int send_errs     :  4 ;
    unsigned int pass_thru_err :  1 ;
    unsigned int ack_err       :  1 ;
    unsigned int extra_recv    :  1 ;
  } _bitval;
};


//--------------------------------------------------------------------------
// void global_sum_handler(GLOBAL_SUM_ERR_SCOPE, unsigned int err_info)
//--------------------------------------------------------------------------
void global_sum_handler(GLOBAL_SUM_ERR_SCOPE err_type, unsigned int err_info)
{  

  const char *ERR_MSG = "This node has error";
  
  switch (err_type) {
  case FLT_OVERFLOW:
    GLOBAL_SUM_HANDLER_SCU_ERR(err_info);    
    break;    
  case SCU_ERR:    
    PRINT_PATH();  
    if (err_info) {
      ScuErrReg *scuErrReg = (ScuErrReg *)err_info;
      int local_err = 0;      
      for (int i = 0; i <= 8; ++i) {
	if (scuErrReg[i]._intval != 0)   {
	  // print the local error message for the first time
	  //---------------------------------------------------------------
	  if (local_err == 0) {
	    local_err = 1;
	    printf("\n%s\n", ERR_MSG);
	  }
	  // check scu err registors
	  //---------------------------------------------------------------
	  if (i < 8) {
	    int errs = scuErrReg[i]._bitval.send_errs;
	    if (errs != 0)       
	      printf("\nWire %d send PARITY ERROR (%d)\n", i, errs);
	    errs = scuErrReg[i]._bitval.recv_errs;
	    if (errs != 0)       
	      printf("\nWire %d recv PARITY ERROR (%d)\n", i, errs);
	    errs = scuErrReg[i]._bitval.extra_recv;      
	    if (errs != 0)   
	      printf("\nWire %d Extra Recv ERROR \n",i);
	    errs = scuErrReg[i]._bitval.ack_err;
	    if (errs != 0)       
	      printf("\nWire %d Acknowledge ERROR \n", i);
	    errs = scuErrReg[i]._bitval.pass_thru_err;
	    if (errs != 0)
	      printf("\nWire %d Passthrough ERROR \n", i);
	  }
	  // check scu status registors
	  //---------------------------------------------------------------
	  else {
	    unsigned int status = scuErrReg[i]._intval;
	    unsigned int bit = 1;
	    for (int wire = 0; wire < 8; ++wire) {
	      if ((status & bit) != 0)                
		printf("\nWire %d recv TIMEOUT.\n", wire);
	      bit = (bit << 1);
	      if ((status & bit) != 0)
		printf("\nWire %d send TIMEOUT.\n", wire);
	      bit = (bit << 1);
	    }
	  }
	}
      }
      if (local_err == 0) {
	printf("\n some other nodes screwed up.\n");
      }
    }
    GLOBAL_SUM_HANDLER_OVERFLOW(err_info);
    break;  
  }
}













 
//--------------------------------------------------------------------------
// GLOBAL_SUM_HANDLER_OVERFLOW  and   GLOBAL_SUM_HANDLER_SCU_ERR
//--------------------------------------------------------------------------
#ifdef MY_MINI_TEST_SYSTEM                     // simple test system 

CPS_END_NAMESPACE
#include <stdlib.h>                            // exit()
CPS_START_NAMESPACE
void GLOBAL_SUM_HANDLER_OVERFLOW(unsigned int err_info)
{
  printf("\nGSUM/BRDCST FAILED\n");    
  exit(-1);
}
inline void GLOBAL_SUM_HANDLER_SCU_ERR(unsigned int err_info)
{   
  printf("\nGLB_SUM OVERFLOW TRUE EXP %d\n", err_info);    
  exit(-1);
}


#else                                          // physics system


CPS_END_NAMESPACE
#include<util/error.h>
CPS_START_NAMESPACE
void GLOBAL_SUM_HANDLER_OVERFLOW(unsigned int err_info)
{
  ERR.Hardware("","", "\nGSUM/BRDCST OVERFLOW TRUE EXP %d\n", err_info);    
}
inline void GLOBAL_SUM_HANDLER_SCU_ERR(unsigned int err_info)
{
  ERR.General("", "", "GSUM OVERFLOW TRUE EXP %d\n", err_info);    
}

#endif        // #ifdef MY_MINI_TEST_SYSTEM




//--------------------------------------------------------------------------
// void PRINT_PATH();
//--------------------------------------------------------------------------
#ifndef GLOBAL_SUM_VERBOSE_ON               // if GLOBAL_SUM_VERBOSE_OFF
void PRINT_PATH()
{
}

#else                                       // if GLOBAL_SUM_VERBOSE_ON

CPS_END_NAMESPACE
#include<global_sum_info.h>
#include <sysfunc.h>
CPS_START_NAMESPACE


void PRINT_PATH()
{
  {
    int physSize[4] =  {SizeT(), SizeX(), SizeY(), SizeZ()};
    int thisPhysCoord[4] = { CoorT(), CoorX(), CoorY(), CoorZ()};
    int thisMachCoord[4];
      
    for (int physDim = 0; physDim < 4; ++physDim) {
      int machDim = SCURemap((SCUDir)(2 * physDim)) / 2; 
      thisMachCoord[machDim]   = thisPhysCoord[physDim];
    }
    printf("Coord[01, 23, 45, 67] = [%d, %d, %d, %d]\n",
	   thisMachCoord[0], thisMachCoord[1],
	   thisMachCoord[2], thisMachCoord[3]);
  }
  printf("max_mode: 0x%x\n", GLOBAL_SUM_INFO_BUF->max_mode);
  printf("brdcst_mode: 0x%x\n", GLOBAL_SUM_INFO_BUF->brdcst_mode);
  printf("send_wire: %d\n", GLOBAL_SUM_INFO_BUF->send_wire);
  printf("recv_wire: %d\n", GLOBAL_SUM_INFO_BUF->recv_wire);
  printf("scu_reset_mask: 0x%x\n", GLOBAL_SUM_INFO_BUF->scu_reset_mask);
  printf("max_num_try: %d\n", GLOBAL_SUM_INFO_BUF->max_num_try);
}

#endif                                      // if GLOBAL_SUM_VERBOSE_ON

CPS_END_NAMESPACE
