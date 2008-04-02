#include "wfm.h"
#include <stdio.h>
#include <stdlib.h>
#ifdef USE_QALLOC
#include <qcdoc_align.h>
#endif
#include <util/time.h>

/* qcdoc_align.h defines a cache touch macro which this code uses. However
   it doesnt get defined unless qcdoc_align.h is included. This causes
   non QALLOC targets to fail to build. I counter this by adding my own
   empty macro for now. Perhaps it would be better if another way to 
   get a cache-touch could be found 
*/
#ifndef cache_touch
#warning "Using empty cache_touch macro"
#define cache_touch(A) ({  })
#endif

/*=========================================================================*/
/* C binding for wilson library:                                           */
/*=========================================================================*/

/*
 * Two of these to allow for efficiently applying dslash to 
 * two 4d fermions simultaneously. On computes while the other 
 * communicates and vice versa.
 * 
 * At some point i should write an interface that has the Ls index
 * innermost. This would be very much more efficient on BG/L and PC's
 * due to the reuse of the loaded gauge field. For QCDOC, the bandwidth
 * is not problematic, so don't bother yet.
 */
static wfm StaticWilsonPAB[2];
static int WilsonLock[2];

void wfm_scope_check(int i);
void wfm_scope_assert(int i);
void wfm_init_internal(int num,WilsonArg *wilson_p);
void wfm_end_internal(int num);

#if 1
void wfm_scope_check(int i)
{
  if ( WilsonLock[i] ) { 
    printf("wilson_init_c: Oops - StaticWilsonPAB is in use");
    exit(-1);
  }
}
void wfm_scope_assert(int i)
{
  if ( !WilsonLock[i] ) { 
    printf("wfm_scope_assert: Oops - StaticWilsonPAB not initialised");
    exit(-1);
  }
}
#endif

void wfm_init_internal(int num,WilsonArg *wilson_p)
{
  wfm_scope_check(num);
  WilsonLock[num] = 1;
  StaticWilsonPAB[num].init(wilson_p);
}
void wfm_end_internal(int num)
{
  wfm_scope_assert(num);
  StaticWilsonPAB[num].end();
  WilsonLock[num] = 0;
}

extern "C" { 

  void wfm_init(WilsonArg *wilson_p)
  {
    wilson_p->instruction_reg_num = 1;
    wfm_init_internal(0,wilson_p);
  }
  
  void wfm_vec_init(WilsonArg *wilson_p)
  {
    wilson_p->instruction_reg_num = 1;
    wfm_init_internal(0,wilson_p);
    wilson_p->instruction_reg_num = 2;
    wfm_init_internal(1,wilson_p);
  }


  void wfm_end(struct WilsonArg *)
  {
    wfm_end_internal(0);
  }
  void wfm_vec_end(struct WilsonArg *wp)
  {
    wfm_end_internal(0);
    wfm_end_internal(1);
  }

  void wfm_mdagm(Float  *chi,        /* chi = MdagM(u) psi          */
		 Float  *u,          /* Gauge field                 */
		 Float  *psi,        /* chi = MdagM(u) psi          */
		 Float  *mp_sq_p,    /* pointer to Sum |M psi|^2    */
		 Float  Kappa)       /* Wilson's kappa parameter    */
  {
    wfm_scope_assert(0);
    StaticWilsonPAB[0].mdagm(chi,u,psi,mp_sq_p,Kappa);
  }

  void wfm_dslash(Float *chi, 
		  Float *u, 
		  Float *psi, 
		  int cb,
		  int dag)
  {
    wfm_scope_assert(0);
    StaticWilsonPAB[0].dslash(chi,u,psi,cb,dag);
  }

  void wfm_m(Float *chi, 
	     Float *u, 
	     Float *psi, 
	     Float kappa)
  {
    wfm_scope_assert(0);
    StaticWilsonPAB[0].m(chi,u,psi,kappa);
  }

  void wfm_mdag(Float *chi, 
		Float *u, 
		Float *psi, 
		Float kappa)
  {
    wfm_scope_assert(0);
    StaticWilsonPAB[0].mdag(chi,u,psi,kappa);
  }

} //extern "C"

void wfm_dslash_two( Float *chi0, Float *chi1, 
		     Float *u, 
		     Float *psi0, Float *psi1,
		     int cb0, int cb1, int dag)
{
   static unsigned long long called=0;
   static Float dslash_time=0.;
   static Float decom_time=0.;
   static Float recom_time=0.;
   static Float comm_time=0.;
   static Float comm_init_time=0.;

   called++;
   dslash_time -=CPS_NAMESPACE::dclock();
    wfm_scope_assert(0);
    wfm_scope_assert(1);

  decom_time -=CPS_NAMESPACE::dclock();
  cache_touch(psi0);
  cache_touch(psi0+4);
  cache_touch(psi0+8);
  cache_touch(psi0+12);
  cache_touch(psi0+16);
  cache_touch(psi0+20);
  StaticWilsonPAB[0].decom(psi0,u,cb0,dag);
  decom_time +=CPS_NAMESPACE::dclock();

  comm_init_time -=CPS_NAMESPACE::dclock();
  StaticWilsonPAB[0].comm_start(cb0);
  comm_init_time +=CPS_NAMESPACE::dclock();

  decom_time -=CPS_NAMESPACE::dclock();
  cache_touch(psi1);
  cache_touch(psi1+4);
  cache_touch(psi1+8);
  cache_touch(psi1+12);
  cache_touch(psi1+16);
  cache_touch(psi1+20);
  StaticWilsonPAB[1].decom(psi1,u,cb1,dag);
  decom_time +=CPS_NAMESPACE::dclock();

  comm_time -=CPS_NAMESPACE::dclock();
  StaticWilsonPAB[0].comm_complete(cb0);
  comm_time +=CPS_NAMESPACE::dclock();

  comm_init_time -=CPS_NAMESPACE::dclock();
  StaticWilsonPAB[1].comm_start(cb1);
  comm_init_time +=CPS_NAMESPACE::dclock();

  recom_time -=CPS_NAMESPACE::dclock();
  cache_touch(StaticWilsonPAB[0].two_spinor);
  cache_touch(StaticWilsonPAB[0].two_spinor+4);
  cache_touch(StaticWilsonPAB[0].two_spinor+8);
  StaticWilsonPAB[0].recon(chi0,u,cb0,dag);
  recom_time +=CPS_NAMESPACE::dclock();

  comm_time -=CPS_NAMESPACE::dclock();
  StaticWilsonPAB[1].comm_complete(cb1);
  comm_time +=CPS_NAMESPACE::dclock();

  recom_time -=CPS_NAMESPACE::dclock();
  cache_touch(StaticWilsonPAB[1].two_spinor);
  cache_touch(StaticWilsonPAB[1].two_spinor+4);
  cache_touch(StaticWilsonPAB[1].two_spinor+8);
  StaticWilsonPAB[1].recon(chi1,u,cb1,dag);
  recom_time +=CPS_NAMESPACE::dclock();

  dslash_time +=CPS_NAMESPACE::dclock();
  unsigned long long offset = 1000;
  unsigned long long interval = 10000;
  if (called%interval==offset){
    CPS_NAMESPACE::print_time("wfm_dslash_two","dslash_time",dslash_time/offset);
    CPS_NAMESPACE::print_time("wfm_dslash_two","decom_time",decom_time/offset);
    CPS_NAMESPACE::print_time("wfm_dslash_two","recom_time",recom_time/offset);
    CPS_NAMESPACE::print_time("wfm_dslash_two","comm_init_time",comm_init_time/offset);
    CPS_NAMESPACE::print_time("wfm_dslash_two","comm_time",comm_time/offset);
  }
  if (called%interval==0){
    called=0;
    decom_time=recom_time=comm_time=dslash_time=comm_init_time=0.;
  }

  return;
}

#if 1
void wfm_dslash_vec( int nvec,
		 Float *chis[],
		 Float *u, 
		 Float *psis[],
		 int cbs[],
		 int dag)
{
  int start;
  wfm_scope_assert(0);
  wfm_scope_assert(1);
  /*
   * Apply dslash to a set of spinors in
   * an efficient way.
   */
  if ( nvec &0x1 ) { 
    start = 1;
    StaticWilsonPAB[0].dslash(chis[0],u,psis[0],cbs[0],dag);
  } else { 
    start = 0;
  }
  for(int i=start;i<nvec;i+=2){ 
    wfm_dslash_two(chis[i],chis[i+1],
		      u,
		      psis[i],psis[i+1],
		      cbs[i] , cbs[i+1],
		      dag );
  }
}

#endif
