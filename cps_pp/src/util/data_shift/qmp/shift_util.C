#include <config.h>
#include <util/gjp.h>
#include <util/error.h>
#include <util/data_shift.h>
#include <qmp.h>
CPS_START_NAMESPACE
#ifndef USE_QMP
#define QMP
#endif

void GlobalDataShift::Shift(int direction, int n_disp){
  if (n_disp==0) return;
 
  SCUDir s_dir,r_dir;
  int a_disp;
  void *send_p,*recv_p,*temp_p;
#ifndef USE_QMP
  if (n_disp>0){
   a_disp = n_disp;
   s_dir = gjp_scu_dir[i*2];
   r_dir = gjp_scu_dir[i*2+1];
  } else {
   a_disp = -n_disp;
   s_dir = gjp_scu_dir[i*2+1];
   r_dir = gjp_scu_dir[i*2];
  }
#else
//  int direction = i;
  int sflag;
  if (n_disp > 0)
    sflag = +1;
  else
    sflag = -1;
#endif

  send_p = addr;
  recv_p = temp_buf;

#ifndef USE_QMP
  SCUDirArgIR Send(send_p,s_dir,SCU_SEND,data_len);
  SCUDirArgIR Recv(recv_p,r_dir,SCU_REC,data_len);
#else
  QMP_msgmem_t msgmem[2];
  QMP_msghandle_t msghandle[2];
  QMP_status_t status;
  QMP_msghandle_t multiple;
#endif

//  sys_cacheflush(0);
  for(int i = 0;i<a_disp-1;i++){

#ifndef USE_QMP
    Send.StartTrans();Recv.StartTrans();
    Send.TransComplete();Recv.TransComplete();
#else
    msgmem[0] = QMP_declare_msgmem((void *)send_p, data_len);
    msgmem[1] = QMP_declare_msgmem((void *)recv_p, data_len);
    msghandle[0] = QMP_declare_send_relative(msgmem[0], direction, sflag, 0);
    msghandle[1] = QMP_declare_receive_relative(msgmem[1], direction, -sflag, 0);
    multiple = QMP_declare_multiple(msghandle, 2);
    QMP_start(multiple);
    status = QMP_wait(multiple);
    if (status != QMP_SUCCESS)
      QMP_error("Error in GlobalDataShift::Shift:%s\n", QMP_error_string(status));
    QMP_free_msghandle(multiple);
    QMP_free_msgmem(msgmem[0]);
    QMP_free_msgmem(msgmem[1]);
#endif
  
    temp_p = send_p;
    send_p = recv_p;
    recv_p = temp_p;

#ifndef USE_QMP
    Send.Addr(send_p);
    Recv.Addr(recv_p);
#endif
  }
#ifndef USE_QMP
  Send.StartTrans();Recv.StartTrans();
  Send.TransComplete();Recv.TransComplete();
#else
  msgmem[0] = QMP_declare_msgmem((void *)send_p, data_len);
  msgmem[1] = QMP_declare_msgmem((void *)recv_p, data_len);
  msghandle[0] = QMP_declare_send_relative(msgmem[0], direction, sflag, 0);
  msghandle[1] = QMP_declare_receive_relative(msgmem[1], direction, -sflag, 0);
  multiple = QMP_declare_multiple(msghandle, 2);
  QMP_start(multiple);
  status = QMP_wait(multiple);
  if (status != QMP_SUCCESS)
    QMP_error("Error in GlobalDataShift::Shift:%s\n", QMP_error_string(status));
  QMP_free_msghandle(multiple);
  QMP_free_msgmem(msgmem[0]);
  QMP_free_msgmem(msgmem[1]);
#endif
  
  if (recv_p != addr)
    memcpy(addr,recv_p,data_len);
}
CPS_END_NAMESPACE
