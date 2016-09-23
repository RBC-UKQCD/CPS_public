/****************************************************************************/
/*! \file
 \brief Definitions for the Wilson fermion matrix code

 $Id: wfm.h,v 1.7 2009/03/23 19:13:32 chulwoo Exp $
*/
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
#include <util/wfm_config.h>
CPS_START_NAMESPACE

struct WilsonArg {
   int local_latt[4];
   int local_comm[4];
   int instruction_reg_num;
//~~
//~~ twisted mass fermions: passes address of spinor_tmp from wfm_h 
//~~ for use by cps_compat.C
//~~
   Float *spinor_tmp;
//~~
};

#ifdef USE_COMMS_SCU
#include <qcdocos/scu_dir_arg.h>
#endif
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
  /*
   * Initialise for Wilson/Clover
   */
void wfm_init(struct WilsonArg *); 
void wfm_end (struct WilsonArg *);
  /*
   * Initialise for DWF
   */
void wfm_vec_init(WilsonArg *wilson_p);
void wfm_vec_end(struct WilsonArg *);
/*
 * Used to fill out the CPS wilson structure
 */

//~~
//~~ added for twisted mass fermions
//~~
void wilson_compat_init(Wilson *wilson_p, WilsonArg *wil);

void wilson_compat_end(Wilson *wilson_p);

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

void wfm_dslash_begin( Float *chi0, 
		       Float *u, 
		       Float *psi0, 
		       int cb0, int dag);
void wfm_dslash_end( Float *chi0, 
		     Float *u, 
		     Float *psi0, 
		     int cb0, int dag);

#ifdef __cplusplus
}
#endif


#ifdef __cplusplus

/*--------------------------------------------------------------------------*/
/* C++ only type definitions                                                */
/*--------------------------------------------------------------------------*/

/* Wilson Class */
class wfm : public WilsonArg {

 public:

  /* Constant definitions 
   * Some are defined in util/wilson.h
   */
  enum {
    //    ND  = 4,
    //    SPINOR_SIZE = 24,
    //    HALF_SPINOR_SIZE = 12,
    //    BLOCK = 12,
    //    PAD_HALF_SPINOR_SIZE = 16,
    //    GAUGE_SIZE = 72,
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

   int   vol;
   Float *spinor_tmp;/* temp spinor needed by mdagm; also redefined */
                     /* as Wilson af[0] and af[1]   s                */
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

  void mdagm(Float  *chi,        /* chi = MdagM(u) psi          */
	     Float  *u,          /* Gauge field                 */
	     Float  *psi,        /* chi = MdagM(u) psi          */
	     Float  *mp_sq_p,    /* pointer to Sum |M psi|^2    */
	     Float  Kappa        /* Wilson's kappa parameter    */
	    );
  ~wfm() {};
//  ~wfm() {printf("wfm %p has been destroyed\n",this);};
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

 protected:

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

/*
 * COMMS support functions. These are virtual so they
 * can be overridden by a user without a recompile of the library.
 */
  void comm_init(void);
  void comm_end(void);
  void comm_start(int cb);
  void comm_complete(int cb);

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

  friend void wfm_dslash_begin( Float *chi0, 
				Float *u, 
				Float *psi0, 
				int cb0, int dag);
  friend void wfm_dslash_end( Float *chi0, 
			      Float *u, 
			      Float *psi0, 
			      int cb0, int dag);
  int IR;
#ifdef USE_COMMS_SCU
/* FIXME. Should really derive a class wfmScu from wfm base class
 * and then 
 * Need to arrange for the different instances to use different IR numbers
 */
  SCUDirArgIR *SendOps[2];
  SCUDirArgIR *RecvOps[2];
  SCUDirArgIR *DA_p[2][16];
  SCUDirArgMulti *DA_multi;
  int LoadDirArgIRs[2];
#endif

};

#endif

CPS_END_NAMESPACE
#endif
