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
extern unsigned long WfmFlops;

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
  if ( dag ) { 
    rec_su3_dag(chi,gauge_notpar,two_spinor,&vol);
  } else { 
    rec_su3(chi,gauge_notpar,two_spinor,&vol);
  }
  return;
}


void wfm::dslash(Float *chi, 
		 Float *u, 
		 Float *psi, 
		 int cb,
		 int dag)
{

  struct timeval t_start, t_stop;
  //  gettimeofday(&t_start,NULL);


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

  /*
  if ( calls ++ < 20 ) {
    gettimeofday(&t_stop,NULL);
    timersub(&t_stop,&t_start,&t_start);
    double flops= 1320 * vol ;
    double secs = t_start.tv_sec + 1.E-6 *t_start.tv_usec;
    printf("Wilson dslash: %f Mflops per node\n",flops/(secs*1000000) );
    printf("Psi %8.8x\n",psi);
    printf("Chi %8.8x\n",chi);
    printf("U   %8.8x\n",u);
    printf("2sp %8.8x\n",shift_table[0][0]);
  }
  */
  WfmFlops += 1320*vol;

  return;
}








