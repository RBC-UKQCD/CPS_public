/*
 *
 *  Calls PAB's assembler routines to give an implementation
 *  of the wilson dslash. Aim is to give very high performance
 *  in a reasonably portable/retargettable library.
 */
#include <util/wfm.h>

#include "wfm_internal.h"
#include <sys/time.h>
#include <stdio.h>
CPS_START_NAMESPACE
 
#ifndef timersub
#define timersub(a, b, result)                                                \
  do {                                                                        \
    (result)->tv_sec = (a)->tv_sec - (b)->tv_sec;                             \
    (result)->tv_usec = (a)->tv_usec - (b)->tv_usec;                          \
    if ((result)->tv_usec < 0) {                                              \
      --(result)->tv_sec;                                                     \
      (result)->tv_usec += 1000000;                                           \
    }                                                                         \
  } while (0)
#endif
                                                                                             
static int calls = 0;

/*
 * This routine does the decompose
 */
unsigned wfm_tbl1,wfm_tbl2,wfm_tbl3,wfm_tbl4;


void wfm::decom(Float *psi, 
	       Float *u, 
	       int cb,
	       int dag)
{
  Float *gauge_par;
  /*Gauge args*/
  if ( cb == 1 ) { 
    gauge_par = (Float *)u + GAUGE_SIZE*vol;
  } else { 
    gauge_par = (Float *)u ;
  }
  if ( dag ) { 
    dec_hsu3_dag(psi,gauge_par,&vol,shift_table[cb]);
  } else { 
    dec_hsu3(psi,gauge_par,&vol,shift_table[cb]);
  }
  return;
}
extern "C" { 
  void rec_su3t(void *psi,void *gauge,void *chiin,void *len,unsigned *tims);
}
void wfm::recon(Float *chi, 
	       Float *u, 
	       int cb,
	       int dag)
{
  Float *gauge_notpar;
  /*Gauge args*/
  if ( cb == 0 ) { 
    gauge_notpar = (Float *)u + GAUGE_SIZE*vol;
  } else { 
    gauge_notpar = (Float *)u ;
  }
  unsigned tbls[2];
  if ( dag ) { 
    rec_su3_dag(chi,gauge_notpar,two_spinor,&vol);
  } else { 
    rec_su3t(chi,gauge_notpar,two_spinor,&vol,tbls);
    //    printf("REC: Time elapsed is %d cycles\n",tbls[1]-tbls[0]);
  }
  DiracOp::CGflops += 1320*vol;
  return;
}


void wfm::dslash(Float *chi, 
		 Float *u, 
		 Float *psi, 
		 int cb,
		 int dag)
{

  /*
   * To a first approximation, we simply
   * remap the arguments into a form acceptable
   * to the assembler, then call it
   */
  decom(psi,u,cb,dag);

  comm_start(cb);

  /*
   * Hackers: you could split here and do something else...
   * Such as DWF fifth dimension, or a clover term etc...
   */

  comm_complete(cb);

  recon(chi,u,cb,dag);

  return;
}

CPS_END_NAMESPACE






