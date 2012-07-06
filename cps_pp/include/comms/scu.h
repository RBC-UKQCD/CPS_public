#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Declarations of communications routines

  $Id: scu.h,v 1.7 2012-07-06 20:22:08 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-07-06 20:22:08 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/scu.h,v 1.7 2012-07-06 20:22:08 chulwoo Exp $
//  $Id: scu.h,v 1.7 2012-07-06 20:22:08 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.7 $
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
void getPlusData(IFloat *rcv_buf, IFloat *send_buf, int len, int mu);

//------------------------------------------------------------------
// getMinusData:
// mu = {0,1,2,3,4} corresponds to {x,y,z,t,s}
//------------------------------------------------------------------
//! Sends data in positive direction/receives data from negative direction.
void getMinusData(IFloat* rcv_buf, IFloat* send_buf, int len, int mu);


extern void getData(IFloat* rcv_buf, int rblklen, int rnumblk,int rstr,
		IFloat* send_buf, int sblklen, int snumblk,int sstr,
		int mu, int sign);

inline void getPlusData(IFloat* rcv_buf, int rblklen, int rnumblk,int rstr,
		IFloat* send_buf, int sblklen, int snumblk,int sstr,
		int mu, int sign){
	getData(rcv_buf,rblklen,rnumblk,rstr,send_buf,sblklen, snumblk,
sstr, mu, 1);
}
inline void getMiunsData(IFloat* rcv_buf, int rblklen, int rnumblk,int rstr,
		IFloat* send_buf, int sblklen, int snumblk,int sstr,
		int mu, int sign){
	getData(rcv_buf,rblklen,rnumblk,rstr,send_buf,sblklen, snumblk,
sstr, mu, -1);
}

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
