#include<config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*! stridetest. A simple test to ensure the stride is 
  correctly implemented in the MPI version by direct 
  comparison with the QCDSP output.

  Code originally from:
  George T. Fleming <gfleming@mps.ohio-state.edu>
  [recieved Wed, 11 Apr 2001 16:03:40 -0400]

  $Id: stridetest.C,v 1.2 2003-10-31 14:15:33 zs Exp $

  A.N.Jackson: ajackson@epcc.ed.ac.uk                       */
/*----------------------------------------------------------*/
CPS_END_NAMESPACE
#include <stdio.h>
#include <comms/sysfunc.h>
CPS_START_NAMESPACE

main() {
  int snd[4]  = { 1, 2, 3, 4 } ;
  int rcv1[4] = { 0, 0, 0, 0 } ;
  int rcv2[4] = { 0, 0, 0, 0 } ;
  
  //                                   blklen, numblk, stride
  SCUDirArg  x(snd,  SCU_XP, SCU_SEND, 4,           1,      1) ;
  SCUDirArg r1(rcv1, SCU_XM, SCU_REC,  4,           1,      1) ;
  SCUDirArg r2(rcv2, SCU_XM, SCU_REC,  1,           4,      1) ;  

  SCUTrans( &r1 ) ;
  SCUTrans( &x ) ;

  // HMM  SCUTransComplete();

  printf("[%i]snd  %i %i %i %i\n", UniqueID(), snd[0],  snd[1],  snd[2],  snd[3] ) ;
  printf("[%i]rcv1 %i %i %i %i\n", UniqueID(), rcv1[0], rcv1[1], rcv1[2], rcv1[3]) ;
  
  SCUTrans( &r2 ) ;
  SCUTrans( &x ) ;
  
  // HMM  SCUTransComplete();

  printf("[%i]rcv2 %i %i %i %i\n", UniqueID(), rcv1[0], rcv1[1], rcv1[2], rcv1[3]) ;
}

/* On QCDSP, this gives the following output:
--------------------------------------------------
snd  1 2 3 4
rcv1 1 2 3 4
rcv2 1 2 3 4
--------------------------------------------------
 */
CPS_END_NAMESPACE
