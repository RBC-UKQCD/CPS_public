#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Declarations of communications routines

  $Id: scu.h,v 1.2 2003-07-24 16:53:53 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/scu.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: scu.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:02  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:12:13  anj
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
//  $RCSfile: scu.h,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/scu.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*******************************************************************
*
* scu.h
*
* Various routines that contriol the scu. There may be several
* routines that perform similar functions.
*
*******************************************************************/

#ifndef INCLUDED_SCU_H
#define INCLUDED_SCU_H

CPS_END_NAMESPACE
#include <util/data_types.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// bsm
//------------------------------------------------------------------
#define TRANSMIT 'x'
#define RECEIVE 'r'
#define IDLE 'i'
extern void bsm(IFloat *, int , int , int , int , int  );

/*! \defgroup sendrecv Fairly high level send/receive routines
  \ingroup comms
  @{ */ 

//------------------------------------------------------------------
// getPlusData:
// mu = {0,1,2,3,4} corresponds to {x,y,z,t,s}
//------------------------------------------------------------------
//! Sends data in negative direction/receives data from positive direction.

extern void getPlusData(IFloat *rcv_buf, IFloat *send_buf, int len, int mu);
 

//------------------------------------------------------------------
// getMinusData:
// mu = {0,1,2,3,4} corresponds to {x,y,z,t,s}
//------------------------------------------------------------------
//! Sends data in positive direction/receives data from negative direction.
extern void getMinusData(IFloat* rcv_buf, IFloat* send_buf, int len, int mu);


//------------------------------------------------------------------
// getMinus2Data: get data from node "-mu, -nu"
// mu = {0,1,2,3} corresponds to {x,y,z,t}
//------------------------------------------------------------------
//! Sends data in positive direction/receives data from negative direction.
extern void 
getMinus2Data(IFloat* rcv_buf, IFloat* send_buf, int len, int mu, int nu);


//------------------------------------------------------------------
// getMinus3Data: get data from node "-mu, -nu, -rho" ( != dir )
// mu = {0,1,2,3} corresponds to {x,y,z,t}
//------------------------------------------------------------------
//! Sends data in positive direction/receives data from negative direction.
extern void 
getMinus3Data(IFloat* rcv_buf, IFloat* send_buf, int len, int dir);

/*! @} */

#endif







CPS_END_NAMESPACE
