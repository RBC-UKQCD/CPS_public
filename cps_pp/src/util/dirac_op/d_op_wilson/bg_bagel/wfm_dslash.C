/*
 *
 *  Calls PAB's assembler routines to give an implementation
 *  of the wilson dslash. Aim is to give very high performance
 *  in a reasonably portable/retargettable library.
 */
#include "wfm.h"
#include "wfm_internal.h"
#include <stdio.h>
#include <stdlib.h>
#ifdef USE_QALLOC
#include <qcdoc.h>
#include <ppc_lib.h>
#include <qcdoc_align.h>
#endif
#include <time.h>
#include <sys/time.h>
#include <util/dirac_op.h>
#ifndef cache_touch
#warning "Using empty cache_touch macro"
#define cache_touch(A)  ( { } )
#endif

#undef DEBUG_BENCHMARK_COMMS
#undef DEBUG_OUTPUT_VECTORS
/*
 * This routine does the decompose
 */
unsigned wfm_tbl1,wfm_tbl2,wfm_tbl3,wfm_tbl4;

static double complex_i[4] __attribute__ ((aligned(16))) = { 0.0, 1.0 } ;

#ifdef DEBUG_OUTPUT_VECTORS
static void site_to_global(int cb, int site, int gbl[4], int local_latt[4] )
{
  int base[4];
  int lsite[4];
  base[0] = CoorX()*local_latt[0];
  base[1] = CoorY()*local_latt[1];
  base[2] = CoorZ()*local_latt[2];
  base[3] = CoorT()*local_latt[3];

  int lx = local_latt[0]/2; 
  lsite[0] = site % (lx); site = (site-lsite[0])/lx;
  lsite[1] = site % (local_latt[1]); site = (site-lsite[1])/local_latt[1];
  lsite[2] = site % (local_latt[2]); site = (site-lsite[2])/local_latt[2];
  lsite[3] = site ;

  lsite[0] = lsite[0]*2 + cb;
  for (int i=0;i<4;i++) {
    gbl[i] = lsite[i]+base[i];
  }
}
#endif

void wfm::decom(Float *psi, 
		Float *u, 
		int cb,
		int dag)
{

  int lcb = (cb + base_parity) & 1;
  Float *gauge_par;
  /*Gauge args*/
  if ( cb == 1 ) { 
    gauge_par = (Float *)u + GAUGE_SIZE*vol;
  } else { 
    gauge_par = (Float *)u ;
  }

  nthread=1;
  /*Turn global checkerboard into the checkerboard within this node*/
  cb = lcb;

  int svol = vol/nthread;
  if ( (svol * nthread) != vol ) { 
    if ( isBoss() ) printf("Bagel threading model broke\n");
    exit(-1);
  }

  decom_internal_arg *args = new decom_internal_arg[nthread];

  args[0].cb = cb;
  args[0].dag = dag;
  args[0].vol = svol;
  args[0].bgl = WFM_BGL;
  args[0].sloppy = SloppyPrecision;

  for(int i=0;i<nthread;i++){
    args[i] = args[0];
    args[i].psi = psi         + svol*i*SPINOR_SIZE;
    args[i].u   = gauge_par   + svol*i*GAUGE_SIZE;
    args[i].shift_table = &shift_table[cb][i*svol*8];
  }

  for ( int i=1;i<nthread;i++) {
    thread_create(decom_internal,(void *)&args[i],i);
  }

  decom_internal(&args[0]);

  for ( int i=1;i<nthread;i++) {
    thread_join(i);
  }

//  printf("nthread=%d\n",nthread);
  delete [] args;
}

int is_stack(unsigned long);
int is_stack(unsigned long addr){
  if ( addr > 0x10000000  ) return 1;
  else return 0;
}
void *wfm::decom_internal(void *pooh)
{
#ifdef USE_THREADS_BGL
  rts_dcache_evict();
#endif
  decom_internal_arg *arg = (decom_internal_arg *)pooh;

  if ( arg->dag ) {

    if ( arg->bgl && arg->sloppy)  
      s_bgl_dec_hsu3_dag(arg->psi,arg->u,&arg->vol,arg->shift_table,complex_i);
    else if ( arg->bgl )
      bgl_dec_hsu3_dag(arg->psi,arg->u,&arg->vol,arg->shift_table,complex_i);
    else if ( arg->sloppy ) 
      s_dec_hsu3_dag(arg->psi,arg->u,&arg->vol,arg->shift_table);
    else
      dec_hsu3_dag(arg->psi,arg->u,&arg->vol,arg->shift_table);

  } else { 

    if ( arg->bgl && arg->sloppy ) 
      s_bgl_dec_hsu3(arg->psi,arg->u,&arg->vol,arg->shift_table,complex_i);
    else if ( arg->bgl )
      bgl_dec_hsu3(arg->psi,arg->u,&arg->vol,arg->shift_table,complex_i);
    else if ( arg->sloppy ) 
      s_dec_hsu3(arg->psi,arg->u,&arg->vol,arg->shift_table);
    else 
      dec_hsu3(arg->psi,arg->u,&arg->vol,arg->shift_table);

  }
#ifdef INCOHERENT_PILE_OF_CRAP_CACHE
  rts_dcache_evict();
#endif
  return ( void *) stdout;
}


void *wfm::recon_internal(void *pooh)
{  
  recon_internal_arg *arg = (recon_internal_arg *)pooh;
#ifdef INCOHERENT_PILE_OF_CRAP_CACHE
  rts_dcache_evict();
#endif

  if ( arg->dag ) { 

    if ( arg->bgl && arg->sloppy )
      s_bgl_rec_su3_dag(arg->chi,arg->u,arg->two_spinor,&arg->vol,complex_i);
    else if ( arg->bgl ){
//printf("recon internal: chi= %e u= %e two spinor= %e vol= %d \n", *(float*)arg->chi,*(float*)arg->u,*(float*)arg->two_spinor,arg->vol);fflush(stdout);
      bgl_rec_su3_dag(arg->chi,arg->u,arg->two_spinor,&arg->vol,complex_i);
    }
    else if( arg->sloppy)
      s_rec_su3_dag(arg->chi,arg->u,arg->two_spinor,&arg->vol);
    else 
      rec_su3_dag(arg->chi,arg->u,arg->two_spinor,&arg->vol);

  } else { 

    if ( arg->bgl && arg->sloppy )
      s_bgl_rec_su3(arg->chi,arg->u,arg->two_spinor,&arg->vol,complex_i);
    else if ( arg->bgl )
      bgl_rec_su3(arg->chi,arg->u,arg->two_spinor,&arg->vol,complex_i);
    else if ( arg->sloppy ) 
      s_rec_su3(arg->chi,arg->u,arg->two_spinor,&arg->vol);
    else
      rec_su3(arg->chi,arg->u,arg->two_spinor,&arg->vol);
  }
#ifdef INCOHERENT_PILE_OF_CRAP_CACHE
  rts_dcache_evict();
#endif

  return NULL;
}


void wfm::recon(Float *chi, 
		Float *u, 
		int cb,
		int dag)
{
  Float *gauge_notpar;
  /*Gauge args*/

  int lcb = (cb + base_parity) & 1;
  nthread=1;

  if ( cb == 0 ) { 
    gauge_notpar = (Float *)u + GAUGE_SIZE*vol;
  } else { 
    gauge_notpar = (Float *)u ;
  }
  cache_touch(gauge_notpar);

  cb = lcb;

  int svol = vol/nthread;
  if ( (svol * nthread) != vol ) { 
    if ( isBoss() ) printf("Bagel threading model broke\n");
    exit(-1);
  }

  recon_internal_arg *args = new recon_internal_arg[nthread];

  args[0].cb     = cb;
  args[0].dag     = dag;
  args[0].vol    = svol;
  args[0].bgl    = WFM_BGL;
  args[0].sloppy = SloppyPrecision;

  for(int i=0;i<nthread;i++){
    args[i]     = args[0];
    args[i].chi = chi               + svol*i*SPINOR_SIZE;
    args[i].u   = gauge_notpar      + svol*i*GAUGE_SIZE;
    args[i].two_spinor = two_spinor + (svol*i*8*PAD_HALF_SPINOR_SIZE*TwoSpinSize())/sizeof(Float);
  }

  for ( int i=1;i<nthread;i++) {
    thread_create(recon_internal,(void *)&args[i],i);
  }

  recon_internal(&args[0]);

  for ( int i=1;i<nthread;i++) {
    thread_join(i);
  }

  delete [] args;
  CPS_NAMESPACE::DiracOp::CGflops += 1320*vol;

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
  /*
   *Pull in the first Psi to cache early
   */
  cache_touch(psi);
  cache_touch(psi+4);
  cache_touch(psi+8);
  cache_touch(psi+12);
  cache_touch(psi+16);
  cache_touch(psi+20);
  decom(psi,u,cb,dag);

#ifdef DEBUG_BENCHMARK_COMMS
  double ndata = 2*2*allbound * 12 * TwoSpinSize() * 1.0E-6 * 100;
  struct timeval start,stop,delta;
  gettimeofday(&start,NULL);

  for(int i=0;i<100;i++) {
#endif
  comm_start(cb); 
  /*
   * Hackers: you could split here and do something else...
   * Such as DWF fith dimension, or a clover term etc...
   * Might as well pull in a few sites worth of pointer table 
   * while we're waiting for the comms
   */
  comm_complete(cb);
#ifdef DEBUG_BENCHMARK_COMMS
  }
  gettimeofday(&stop,NULL);
  timersub(&stop,&start,&delta);
  double seconds = delta.tv_usec * 1.0E-6 + delta.tv_sec;
  if ( isBoss() ) printf("Comms %le MB in %le seconds = %le MB/s\n",ndata,seconds,ndata/seconds);
  ndata = 2*2*allbound * 12 * TwoSpinSize() ;
  if ( isBoss() ) printf("ndata = %d \n",ndata);
#endif
#ifdef DEBUG_OUTPUT_VECTORS 
  static int file_cnt;
  {
    char buf[256];
  sprintf(buf,"2spin.%d.%d",UniqueID(),file_cnt);
  FILE *fp = fopen(buf,"w");
  for(int i=0;i<vol;i++) {
    for(int pm=0;pm<2;pm++){
      for(int mu=0;mu<4;mu++){
        int offset = interleave_site(pm,mu,i);
        for(int s=0;s<2;s++){
        for(int c=0;c<3;c++){
        for(int r=0;r<2;r++){
	  int scri;
          if ( WFM_BGL ) scri = r + s*6+c*2;        
	  else scri = r + s*2+c*4;        
          int gbl[4];
          site_to_global(cb, i, gbl, local_latt );
          if ( SloppyPrecision ) {
            float * pointer = (float *) two_spinor;
            fprintf(fp,"%d %d %d %d %d %d %d %d %d %e\n",gbl[0],gbl[1],gbl[2],gbl[3],
              pm,mu,s,c,r,pointer[PAD_HALF_SPINOR_SIZE*offset+scri]);
          } else { 
            Float * pointer = (Float *) two_spinor;
            fprintf(fp,"%d %d %d %d %d %d %d %d %d %e\n",gbl[0],gbl[1],gbl[2],gbl[3],
                    pm,mu,s,c,r,pointer[PAD_HALF_SPINOR_SIZE*offset+scri]);
          }
        }}}
      }
    }
  }
  fclose(fp);}
#endif 

  cache_touch(two_spinor);
  cache_touch(two_spinor+4);
  cache_touch(two_spinor+8); 
  recon(chi,u,cb,dag);
#ifdef DEBUG_OUTPUT_VECTORS
  {
  char buf[256];
  sprintf(buf,"recon.%d.%d",UniqueID(),file_cnt++);
  FILE *fp = fopen(buf,"w");
  for(int i=0;i<vol;i++) {
    for(int pm=0;pm<2;pm++){
      for(int mu=0;mu<4;mu++){
        int offset = interleave_site(pm,mu,i);
        for(int s=0;s<4;s++){
        for(int c=0;c<3;c++){
        for(int r=0;r<2;r++){
	  int scri;
          scri = r + s*6+c*2;        
	  Float * pointer = (Float *) chi;
          int gbl[4];
          site_to_global(cb, i, gbl, local_latt );
	  fprintf(fp,"%d %d %d %d %d %d %d %d %d %d %e\n",gbl[0],gbl[1],gbl[2],gbl[3],
                   i,pm,mu,s,c,r,pointer[SPINOR_SIZE*offset+scri]);
        }}}
      }
    }
  }
  fclose(fp);}
  exit(0);
#endif 
  return;
}


