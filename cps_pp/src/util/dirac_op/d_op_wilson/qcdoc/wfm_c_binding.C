#include <util/wfm.h>
#include <stdio.h>
#include <stdlib.h>
CPS_START_NAMESPACE
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
static int pad[65536];
static wfm StaticWilsonPAB[2];
static int WilsonLock[2];

void wfm_scope_check(int i);
void wfm_scope_assert(int i);
void wfm_init_internal(int num,WilsonArg *wilson_p);
void wfm_end_internal(int num);

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

void wfm_init_internal(int num,WilsonArg *wilson_p)
{
//  printf("wfm_init_internal(%d,%p)\n",num,wilson_p);
  wfm_scope_check(num);
  WilsonLock[num] = 1;
  StaticWilsonPAB[num].init(wilson_p);
}
void wfm_end_internal(int num)
{
//  printf("wfm_end_internal(%d)\n",num);
  wfm_scope_assert(num);
  StaticWilsonPAB[num].end();
  WilsonLock[num] = 0;
}

extern "C" { 

  void wfm_init(WilsonArg *wilson_p)
  {
//    printf("wfm_init(%p)\n",wilson_p);
    wilson_p->instruction_reg_num = 10;
    wfm_init_internal(0,wilson_p);
  }
  
  void wfm_vec_init(WilsonArg *wilson_p)
  {
    wilson_p->instruction_reg_num = 10;
    wfm_init_internal(0,wilson_p);
    wilson_p->instruction_reg_num = 13;
    wfm_init_internal(1,wilson_p);
  }


  void wfm_end(struct WilsonArg *wp)
  {
//    printf("wfm_end(%p)\n",wp);
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

}


void wfm_dslash_begin( Float *chi0, 
		       Float *u, 
		       Float *psi0, 
		       int cb0, int dag)
{
  wfm_scope_assert(0);
  StaticWilsonPAB[0].decom(psi0,u,cb0,dag);
  StaticWilsonPAB[0].comm_start(cb0);
}

void wfm_dslash_end( Float *chi0, 
		     Float *u, 
		     Float *psi0, 
		     int cb0, int dag)
{
  wfm_scope_assert(0);
  StaticWilsonPAB[0].comm_complete(cb0);
  StaticWilsonPAB[0].recon(chi0,u,cb0,dag);
}


void wfm_dslash_two( Float *chi0, Float *chi1, 
		     Float *u, 
		     Float *psi0, Float *psi1,
		     int cb0, int cb1, int dag)
{
  wfm_scope_assert(0);
  wfm_scope_assert(1);

  StaticWilsonPAB[0].decom(psi0,u,cb0,dag);

  StaticWilsonPAB[0].comm_start(cb0);

  StaticWilsonPAB[1].decom(psi1,u,cb1,dag);

  StaticWilsonPAB[0].comm_complete(cb0);

  StaticWilsonPAB[1].comm_start(cb1);

  StaticWilsonPAB[0].recon(chi0,u,cb0,dag);

  StaticWilsonPAB[1].comm_complete(cb1);

  StaticWilsonPAB[1].recon(chi1,u,cb1,dag);

  return;
}
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

CPS_END_NAMESPACE
