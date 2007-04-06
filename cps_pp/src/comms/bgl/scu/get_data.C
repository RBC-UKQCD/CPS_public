#include<config.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <util/gjp.h>
#include <comms/scu.h>
#include <comms/bgl_net.h>
#include <sys/bgl/bgl_sys_all.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// mu = {0,1,2,3,4} corresponds to {x,y,z,t,s}
//------------------------------------------------------------------

void getPlusData(IFloat *rcv_buf, IFloat *send_buf, int len, int mu)
{
  int dir_se;
  int dir_re;
  int i;
//  printf("Node %d:getPlusData(%p,%p,%d,%d)\n",UniqueID(),rcv_buf,send_buf,len,mu);

  // x+, y+, z+, t+
//  if( (mu == 0) || (mu == 1) || (mu == 2) || (mu == 3) ){
  if(GJP.Nodes(mu)>1){
    dir_se = (2 * mu) + 1;
    dir_re = (2 * mu);
    // Send-receive data

    IFloat *send_p=send_buf;
    IFloat *rcv_p= rcv_buf;
    while (len >24){
      BGLCPSTorus_send(dir_se, 24, send_p);
      BGLCPSTorus_recv(dir_re, 24, rcv_p);
      len -=24;
      send_p+=24;rcv_p+=24;   
    } 
    if (len>0){
    BGLCPSTorus_send(dir_se, len, send_p);
    BGLCPSTorus_recv(dir_re, len, rcv_p);
    }
  }

  // all other are local
  else {
    for(i = 0; i < len; ++i) {
      *rcv_buf++ = *send_buf++;
    }
  }

}


void getMinusData(IFloat* rcv_buf, IFloat* send_buf, int len, int mu)
{
  int dir_se;
  int dir_re;
  int i;

  // x-, y-, z-, t-
//  if( (mu == 0) || (mu == 1) || (mu == 2) || (mu == 3) ){
  if(GJP.Nodes(mu)>1){
    dir_se = (2 * mu);
    dir_re = (2 * mu) + 1;

    IFloat *send_p=send_buf;
    IFloat *rcv_p= rcv_buf;
    while (len >24){
      BGLCPSTorus_send(dir_se, 24, send_p);
      BGLCPSTorus_recv(dir_re, 24, rcv_p);
      len -=24;
      send_p+=24;rcv_p+=24;   
    } 
    if (len>0){
    BGLCPSTorus_send(dir_se, len, send_p);
    BGLCPSTorus_recv(dir_re, len, rcv_p);
    }
  }

  // all other are local
  else {
    for(i = 0; i < len; ++i) {
      *rcv_buf++ = *send_buf++;
    }
  }
}

void getMinus2Data(IFloat* rcv_buf, IFloat* send_buf, int len, int mu, int nu)
{
  for(int i = 0; i < len; ++i) {
    *rcv_buf++ = *send_buf++;
  }
}

void getMinus3Data(IFloat* rcv_buf, IFloat* send_buf, int len, int dir)
{
  for(int i = 0; i < len; ++i) {
    *rcv_buf++ = *send_buf++;
  }
}
CPS_END_NAMESPACE
