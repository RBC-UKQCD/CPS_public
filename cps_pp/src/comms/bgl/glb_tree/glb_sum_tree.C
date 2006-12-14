#include<config.h>
CPS_START_NAMESPACE

//--------------------------------------------------------------------
/*
 *  glb_sum.C 
 *  Global sum using a torus nearest neighbor communications method.
 *  This is in double precision only.
 *  This is accurate for a torus only.
 *  This is a slow way to do global sums. Use it until a faster 
 *  "hardware" way is available.
 */
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <sys/bgl/bgl_sys_all.h>
#include <comms/bgl_net.h>
#include <comms/scu.h>
#include<comms/glb.h>
#include<util/gjp.h>

extern "C" void bgl_tree_double_sum(double *, double *);

CPS_START_NAMESPACE

static double transmit_buf[2] __attribute__((aligned(BGL_QUAD_ALIGNSIZE))); 
static double receive_buf[2]  __attribute__((aligned(BGL_QUAD_ALIGNSIZE)));
static double gsum_buf[2]     __attribute__((aligned(BGL_QUAD_ALIGNSIZE)));

void glb_sum(Float * gsum_p){
  double sum;
  double value;
  int pir;
  int dir_se;
  int dir_re;

  // Set internal direction
  // The direction is translated. You must use the values that 
  // will be translated to the internal directions (6,7)
  dir_se = bgl_cps_dir[6];
  dir_re = bgl_cps_dir[7];

  // Get the processor id
  pir = rts_get_processor_id();
  
  // Set the local value
  gsum_buf[0] = *gsum_p;

  // Get/send the data from/to the other cpu
  transmit_buf[0] = gsum_buf[0];
  BGLCPSTorus_send(dir_se, 2, transmit_buf);
  BGLCPSTorus_recv(dir_re, 2, receive_buf);

  // Add it up
  gsum_buf[0] += receive_buf[0];

  // If core 0 then do the global sum with the two node-added-value
  if(pir == 0){
    value = gsum_buf[0];
    bgl_tree_double_sum(&value, &sum);
    *gsum_p = sum;
  } else {
    sum = 0.0;
  }

  // Send sum to the other node
  transmit_buf[0] = sum;
  BGLCPSTorus_send(dir_se, 2, transmit_buf);
  BGLCPSTorus_recv(dir_re, 2, receive_buf);
  
  // If node 1 then set the result  *gsum_p = receive data
  if(pir == 1){
    *gsum_p = receive_buf[0];
  }


}


CPS_END_NAMESPACE
