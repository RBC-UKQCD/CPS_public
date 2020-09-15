#ifndef _MESONFIELD_CONTROLS_H
#define _MESONFIELD_CONTROLS_H

CPS_START_NAMESPACE

//Basic tunings performed on laptop (i7-4700mq 1 thread per core)
#define MF_COMPUTE_BI 32
#define MF_COMPUTE_BJ 12
#define MF_COMPUTE_BP 256

#undef USE_INNER_BLOCKING
#define MF_COMPUTE_BII -1
#define MF_COMPUTE_BJJ -1
#define MF_COMPUTE_BPP -1


#ifdef KNL_OPTIMIZATIONS
//Insert KNL tunings
#define MF_COMPUTE_BI 8
#define MF_COMPUTE_BJ 2
#define MF_COMPUTE_BP 1024

#define USE_INNER_BLOCKING
#define MF_COMPUTE_BII 8
#define MF_COMPUTE_BJJ 2
#define MF_COMPUTE_BPP 4
#endif

class BlockedMesonFieldArgs{ //Initialized to the above but can be modified by user
public:
  //Outer blocking
  static int bi;
  static int bj;
  static int bp;  

  //Inner blocking (if enabled)
  static int bii;
  static int bjj;
  static int bpp;
};


CPS_END_NAMESPACE


#endif
