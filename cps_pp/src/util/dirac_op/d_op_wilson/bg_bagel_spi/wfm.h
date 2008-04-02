#include <config.h>
/****************************************************************************/
/* Jan/2003                                                                 */
/* Peter Boyle.                                                             */
/*                                                                          */
/* wfm.h                                                                    */
/* C header file for the Wilson Fermion Matrix library                      */
/*                                                                          */
/****************************************************************************/

#ifndef INCLUDED_WFM_PAB_H
#define INCLUDED_WFM_PAB_H

#include <util/wilson.h>
#include "wfm_config.h"

#ifdef USE_COMMS_SPI 
#include <spi/bgp_SPI.h>
#endif


#ifdef USE_THREADS_BGL
#include <rts.h>
#endif
extern bool wfm_bgl_def ;
class WilsonArg {
public:
   int local_latt[4];
   int local_comm[4];
   int instruction_reg_num;
   bool SloppyPrecision;
   bool WFM_BGL;
   int CoreCount;

   WilsonArg () {
     CoreCount=1;
     SloppyPrecision=false;
#if 0
#ifdef BLUE_GENE_IS_DEFAULT 
     WFM_BGL = wfm_bgl_def;
#else 
     WFM_BGL = false;
#endif
# else
     WFM_BGL = true;
#endif
   }
};

struct decom_internal_arg {
  Float *psi;
  Float *u;
  int cb;
  int dag;
  int vol;
  unsigned long *shift_table;
  int sloppy;
  int bgl;
};
struct recon_internal_arg {
  Float *chi;
  Float *u;
  Float *two_spinor;
  int cb;
  int dag;
  int vol;
  int sloppy;
  int bgl;
};

#ifdef __cplusplus
extern "C" { 
#endif

/*
 * "C" linkage routines. These allocate internal 2spinor storage in a static
 * and scoped locked wilson struct. Interface constrains that only dirac
 * operation can be in progress at a time since the 2spinor storage would
 * be under conflict. 
 *
 * C++ interface will allows for starting multiple dirac
 * operations at a time, which enables us to overlap comms and compute
 * for two different 4d source vectors. Useful for overlap and domain wall
 * approaches.
 *
 */
void wfm_init(struct WilsonArg *); 
void wfm_end (struct WilsonArg *);
void wfm_vec_init(struct WilsonArg *wilson_p);
void wfm_vec_end(struct WilsonArg *);
/* 
 * Used to fill out the CPS wilson structure
 */
CPS_START_NAMESPACE
void wilson_compat_init(CPS_NAMESPACE::Wilson *wilson_p);
void wilson_compat_end(CPS_NAMESPACE::Wilson *wilson_p);
CPS_END_NAMESPACE

void wfm_mdagm(Float  *chi,        /* chi = MdagM(u) psi          */
	       Float  *u,          /* Gauge field                 */
	       Float  *psi,        /* chi = MdagM(u) psi          */
	       Float  *mp_sq_p,    /* pointer to Sum |M psi|^2    */
	       Float  Kappa);      /* Wilson's kappa parameter    */

void wfm_dslash(Float *chi, 
		Float *u, 
		Float *psi, 
		int cb,
		int dag);

void wfm_m(Float *chi, 
	   Float *u, 
	   Float *psi, 
	   Float kappa);

void wfm_mdag(Float *chi, 
	      Float *u, 
	      Float *psi, 
	      Float kappa);

void wfm_dslash_two( Float *chi0, Float *chi1,
		     Float *u,
		     Float *psi0, Float *psi1,
		     int cb0, int cb1, int dag);

void wfm_dslash_vec( int nvec,
		     Float *chis[],
		     Float *u,
		     Float *psis[],
		     int cbs[],
		     int dag);

#ifdef __cplusplus
}
#endif


#ifdef __cplusplus


#ifdef USE_COMMS_SCU
#include <qcdocos/scu_dir_arg.h>
#endif

/*--------------------------------------------------------------------------*/
/* C++ only type definitions                                                */
/*--------------------------------------------------------------------------*/

/* Wilson Class */
class wfm : public WilsonArg {

 public:

 virtual ~wfm(){}
  /* Constant definitions */
  enum {
// Hack to make CPS include/util/wilson.h and wfm.h
#ifndef INCLUDED_WILSON_H
    ND  = 4,
    SPINOR_SIZE = 24,
    HALF_SPINOR_SIZE = 12,
    BLOCK = 12,
    GAUGE_SIZE = 72,
#endif
    Nmu = 4,
    Ncb = 2,
    NMinusPlus = 2,
    Minus=0,
    Plus =1
  };


/*
 * Pointer to an array with addressing offsets into the (mu interleaved) 
 * two spinor. The index cb is the _source_ checkerboard, and destination
 * checkerboard is 1-cb.
 *
 * Entry n in shift_table[cb] where n = site *8 + pm*4 + mu
 *
 * gives the pointer to the spinor array of the opposite cb, with
 *
 *   new_site = (site - mu) if  pm ==0 == Minus
 *              (site + mu) if  pm ==1 == Plus
 *   
 *   The 2spinor fermion fields are interleaved, so that the index
 *   given by shift_table is then :
 *  
 *   newsite*8+pm*4+mu
 *
 * Kernels are:
 *
 *    dec_hsu3_8way():
 *                        (1+gamma_mu) PSI(site) -> shift_table[site*8 + mu]
 *          Umu^dag(site) (1-gamma_mu) PSI(site) -> shift_table[site*8+4+mu]
 *
 *    dec_dag_hsu3_8way():
 *                        (1-gamma_mu) PSI(site) -> shift_table[site*8 + mu]
 *          Umu^dag(site) (1+gamma_mu) PSI(site) -> shift_table[site*8+4+mu]
 *
 *    rec_su3_8way()
 *          PHI(site) = Umu(site) TO_4SPIN(+,mu,CHI(site*8 + mu))
 *                    +           TO_4SPIN(-,mu,CHI(site*8+4+mu))
 *
 *    rec_dag_su3_8way()
 *          PHI(site) = Umu(site) TO_4SPIN(+,mu,CHI(site*8 + mu))
 *                    +           TO_4SPIN(-,mu,CHI(site*8+4+mu))
 *
 *
 */

   unsigned long *shift_table[Ncb];
   unsigned long *face_table[Ncb][NMinusPlus][ND];

   int base_parity; // The global parity of local site 0 0 0 0 on this node

   int   vol;
   Float *spinor_tmp;/* temp spinor needed by mdagm                 */
   Float *two_spinor;/* point. array to 8 interleaved fwd proj half spinors  */
   Float *send_bufs[NMinusPlus][ND];
   Float *recv_bufs[NMinusPlus][ND];
   int  send_offset[NMinusPlus][ND];
   Float *XBaseAddr[Ncb][NMinusPlus];
   int nbound[ND];
   int allbound;

  //   SCUDirArgMulti *comm_f;
  //   SCUDirArgMulti *comm_b;

   /*
    * Public Methods
    */
  void init(WilsonArg *arg); 
  void end (void);

  /*
   * Behaviour control ...  padding and precision & which kernels
   */
  int PAD_HALF_SPINOR_SIZE ; // This is a function of Sloppy and Bgl

  bool SloppyPrecision ;
  bool WFM_BGL;
  int nthread ;

  void CoreCount(int);
  int TwoSpinSize() { 
    if ( SloppyPrecision ) return sizeof(float) ; 
    else return sizeof(Float);
  };

  void face_scatter(Float *TwoSpinor,
	   Float *RcvBuf, 
			   unsigned long*FaceTable,
			   unsigned long);

  /* Feb 2006
   * Multi-core support
   */
//  static void wfm::thread_create(void *(*func)(void *), void *arg,int id);
//  static void wfm::thread_join  (int id);
  static void thread_create(void *(*func)(void *), void *arg,int id);
  static void thread_join  (int id);


  void mdagm(Float  *chi,        /* chi = MdagM(u) psi          */
	     Float  *u,          /* Gauge field                 */
	     Float  *psi,        /* chi = MdagM(u) psi          */
	     Float  *mp_sq_p,    /* pointer to Sum |M psi|^2    */
	     Float  Kappa        /* Wilson's kappa parameter    */
	    );

  void dslash(Float *chi, 
	      Float *u, 
	      Float *psi, 
	      int cb,
	      int dag);

  void m(Float *chi, 
	 Float *u, 
	 Float *psi, 
	 Float kappa);
	 
  void mdag(Float *chi, 
	    Float *u, 
	    Float *psi, 
	    Float kappa);

  /*
   * Protected methods, can be overridden in derived classes.
   * Comms in this base class are virtual. 
   *
   * The base class just does a local copy.
   *
   * To implement comms for a new platform, simply
   * derive a platform specific class and insert your own calls.
   *
   */		

 //protected:

  /*
   * Wilson dslash may be split into four parts that allows 
   * the calculation of multiple operations (using multiple instances
   * of wfm) to be interleaved, with perfectly overlapping
   * communication and computation even on very small volumes.
   *
   * Generalisation from the implementation of wilson_dslash is obvious.
   */

  /*
   * C++ interface to the ASM kernels
   */
  void decom(Float *psi, 
	     Float *u, 
	     int cb,
	     int dag);


  void recon(Float *chi, 
	     Float *u, 
	     int cb,
	     int dag);

  /*Wrap simpler routines for BG/L co-routine pleasure*/
  static void *decom_internal(void *pooh);
  static void *recon_internal(void *pooh);



/*
 * COMMS support functions. These are virtual so they
 * can be overridden by a user without a recompile of the library.
 */
  virtual void comm_init(void);
  virtual void comm_end(void);
  virtual void comm_start(int cb);
  virtual void comm_complete(int cb);

  int isBoss(void);

  /*
   * BGL specific hackery.
   * A masterpiece of low level hacking in PAB's modest opinion
   */ 
  void mmu_optimise(void);
  void mmu_print(void);

#ifdef USE_COMMS_SPI 
  void spi_init();
  void spi_fwd_start(void);
  void spi_bwd_start(void);
  void spi_fwd_wait(void );
  void spi_bwd_wait(void );
#endif

 private:
/*
 * POINTER TABLE support functions
 */
  void pointers_init(void);
  void pointers_end(void);
  int  local_parity(int x,int y,int z,int t);
  int  local_psite(int addr[4],int latt[4]);
  int  interleave_site(int pm,int mu, int site);

  friend void wfm_dslash_two( Float *chi0, Float *chi1,
			      Float *u,
			      Float *psi0, Float *psi1,
			      int cb0, int cb1, int dag);

#ifdef USE_COMMS_SCU
/* FIXME. Should really derive a class wfmScu from wfm base class
 * and then
 * Need to arrange for the different instances to use different IR numbers
 */
  SCUDirArgIR *SendOps[Ncb];
  SCUDirArgIR *RecvOps[Ncb];
  SCUDirArgIR *DA_p[Ncb][16];
  SCUDirArgMulti *DA_multi;
  int LoadDirArgIRs[2];
#endif

  // Wonder what this does ? (BJ) perhaps it should be a 
  // SCU only beastie?
  int IR;

#ifdef USE_COMMS_QMP

  // Everything that is allocated for comms needs a QMP_mem_t
  QMP_mem_t* send_bufs_mem_t[NMinusPlus][ND];
  QMP_mem_t* recv_bufs_mem_t[NMinusPlus][ND];

  // We now receive the face directly into two spinor
  // so it has to be communicable, and so it needs a mem_t
  QMP_mem_t* two_spinor_mem_t;
  QMP_msgmem_t qmp_send_ops_msgmem_t[Ncb][8];
  QMP_msgmem_t qmp_recv_ops_msgmem_t[Ncb][8];
  QMP_msghandle_t qmp_multi_handles[Ncb];

  // This is a hack to catch if the code is being run in 
  // local comms mode. In that case allocating a multi of 
  // length 0 will barf...
  bool nonlocal_comms;

#endif

#ifdef USE_COMMS_SPI 
//DMA_CounterGroup_t *rec_counter_group;
//DMA_CounterGroup_t *inj_counter_group;
//DMA_InjFifoGroup_t *inj_fifo_group;
int comm_bytes_total;
#endif
  
};


#endif

#endif
