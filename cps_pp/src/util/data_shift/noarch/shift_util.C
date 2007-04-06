#include <config.h>
#include <util/gjp.h>
#include <util/error.h>
#include <util/data_shift.h>
#include <comms/scu.h>
CPS_START_NAMESPACE
void GlobalDataShift::Shift(int mu, int n_disp){
  if (n_disp==0) return;
 
  VRB.Flow(cname,"Shift(i,i)","mu=%d n_disp=%d\n",mu,n_disp);
//  SCUDir s_dir,r_dir;
  int a_disp;
  Float *send_p,*recv_p,*temp_p;
  if (n_disp>0){
   a_disp = n_disp;
//   s_dir = gjp_scu_dir[i*2];
//   r_dir = gjp_scu_dir[i*2+1];
  } else {
   a_disp = -n_disp;
//   s_dir = gjp_scu_dir[i*2+1];
//   r_dir = gjp_scu_dir[i*2];
  }

  send_p = (Float *)addr;
  recv_p = (Float *)temp_buf;

//  SCUDirArgIR Send(send_p,s_dir,SCU_SEND,data_len);
//  SCUDirArgIR Recv(recv_p,r_dir,SCU_REC,data_len);

//  sys_cacheflush(0);
  for(int i = 0;i<a_disp-1;i++){
//    Send.StartTrans();Recv.StartTrans();
//    Send.TransComplete();Recv.TransComplete();
    if (n_disp>0) getMinusData(recv_p,send_p,data_len/sizeof(Float),mu);
    else  getPlusData(recv_p,send_p,data_len/sizeof(Float),mu);
    temp_p = send_p;
    send_p = recv_p;
    recv_p = temp_p;
//    Send.Addr(send_p);
//    Recv.Addr(recv_p);
  }
    if (n_disp>0) getMinusData(recv_p,send_p,data_len/sizeof(Float),mu);
    else  getPlusData(recv_p,send_p,data_len/sizeof(Float),mu);
//  Send.StartTrans();Recv.StartTrans();
//  Send.TransComplete();Recv.TransComplete();

  if (recv_p != addr)
    memcpy(addr,recv_p,data_len);
}
CPS_END_NAMESPACE
