#include <config.h>
CPS_START_NAMESPACE
/*--------------------------------------------------------------------*/
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.
*/
/*--------------------------------------------------------------------*/

/***************************************************************************/
/*                                                                         */
/* wilson_dslash: It calculates chi = Dslash * psi, or                     */
/*                                chi = DslashDag * psi, where Dslassh is  */
/* the Wilson fermion Dslash matrix. Dslash is a function of the gauge     */
/* fields u.                                                               */
/* cb = 0/1 denotes input on even/odd checkerboards i.e.                   */
/* cb = 1 --> Dslash_EO, acts on an odd column vector and produces an      */
/* even column vector                                                      */
/* cb = 0 --> Dslash_OE, acts on an even column vector and produces an     */
/* odd column vector,                                                      */
/*                                                                         */
/* dag = 0/1  results into calculating Dslash or DslashDagger.             */
/* lx,ly,lz,lt is the lattice size.                                        */
/*                                                                         */
/* This routine is to be used with scalar machines.                        */
/*                                                                         */
/* WARNING:                                                                */
/*                                                                         */
/* This set of routines will work only if the node sublattices have        */
/* even number of sites in each direction.                                 */
/*                                                                         */
/***************************************************************************/

CPS_END_NAMESPACE
#include <util/data_types.h>
#include <util/vector.h>
#include <util/wilson.h>
#include <util/dirac_op.h>
#include <util/error.h>
#include <util/time_cps.h>
#include <comms/scu.h>
#include <sys/bgl/bgl_sys_all.h>
#if TARGET != BGL
inline double rts_get_timebase() {return 0;}
#endif
CPS_START_NAMESPACE

/***************************************************************************/
/* Function declarations                                                   */
/***************************************************************************/
void wfm_comm();

void wfm_cmat_two_spproj(IFloat *af0, IFloat *af1, IFloat *af2, IFloat *af3,
			 IFloat *ab0, IFloat *ab1, IFloat *ab2, IFloat *ab3,
			 IFloat *u, 
			 IFloat *psi, 
			 IFloat sign, 
			 Wilson *wilson_p, 
			 int cb);

void wfm_mat_trick(IFloat *chi, 
		   IFloat *u,
		   IFloat *af0, IFloat *af1, IFloat *af2, IFloat *af3,
		   IFloat sign, 
		   Wilson *wilson_p, 
		   int cb);

void wfm_trick(IFloat *chi, 
	       IFloat *ab0, IFloat *ab1, IFloat *ab2, IFloat *ab3,
	       IFloat sign, 
	       Wilson *wilson_p, 
	       int cb);



/***************************************************************************/
/* external variables                                                      */
/***************************************************************************/
double wfm_reg[64]   __attribute__((aligned(BGL_QUAD_ALIGNSIZE)));
static int count=0;


void wilson_dslash(IFloat *chi, 
		   IFloat *u, 
		   IFloat *psi, 
		   int cb,
		   int dag,
		   Wilson *wilson_p)
{
  char *fname = "wilson_dslash";
  IFloat sign;
  IFloat *u_eo[2];
  IFloat *af0;
  IFloat *af1;
  IFloat *af2;
  IFloat *af3;
  IFloat *ab0;
  IFloat *ab1;
  IFloat *ab2;
  IFloat *ab3;
  int cbn;
  unsigned long long start_time;
  unsigned long long stop_time;
  unsigned long long start_time_dslash;
  unsigned long long stop_time_dslash;
  unsigned long long cmat_two_spproj_time;
  unsigned long long comm_time;
  unsigned long long mat_trick_time;
  unsigned long long trick_time;
  unsigned long long dslash_time;
  int dslash_time_i;
  int cb_vol;
  IFloat dslash_perf;
  Float total_time=0.;
  
  start_time_dslash = rts_get_timebase();
  total_time -=dclock();

  /*------------------------------------------------------------------------*/
  /* Save double hummer floating point registers                            */
  /*------------------------------------------------------------------------*/
  QuadStore(wfm_reg  ,  0);
  QuadStore(wfm_reg+2,  1);
  QuadStore(wfm_reg+4,  2);
  QuadStore(wfm_reg+6,  3);
  QuadStore(wfm_reg+8,  4);
  QuadStore(wfm_reg+10, 5);
  QuadStore(wfm_reg+12, 6);
  QuadStore(wfm_reg+14, 7);
  QuadStore(wfm_reg+16, 8);
  QuadStore(wfm_reg+18, 9);  
  QuadStore(wfm_reg+20, 10);
  QuadStore(wfm_reg+22, 11);
  QuadStore(wfm_reg+24, 12);
  QuadStore(wfm_reg+26, 13);
  QuadStore(wfm_reg+28, 14);
  QuadStore(wfm_reg+30, 15);
  QuadStore(wfm_reg+32, 16);
  QuadStore(wfm_reg+34, 17);
  QuadStore(wfm_reg+36, 18);
  QuadStore(wfm_reg+38, 19);  
  QuadStore(wfm_reg+40, 20);
  QuadStore(wfm_reg+42, 21);
  QuadStore(wfm_reg+44, 22);
  QuadStore(wfm_reg+46, 23);
  QuadStore(wfm_reg+48, 24);
  QuadStore(wfm_reg+50, 25);
  QuadStore(wfm_reg+52, 26);
  QuadStore(wfm_reg+54, 27);
  QuadStore(wfm_reg+56, 28);
  QuadStore(wfm_reg+58, 29);
  QuadStore(wfm_reg+60, 30);
  QuadStore(wfm_reg+62, 31);

/*--------------------------------------------------------------------------*/
/* Initializations                                                          */
/*--------------------------------------------------------------------------*/
  /* For dag 0/1 set sign to +/- 1 */
  if(dag == 0)
    sign = 1.0;
  else if(dag == 1)
    sign = -1.0;
  else{
    ERR.General(" ",fname,"dag must be 0 or 1");
  }

  /* set cbn = !cb */
  if(cb == 0)
    cbn = 1;
  else if(cb == 1)
    cbn = 0;
  else{
    ERR.General(" ",fname,"cb must be 0 or 1");
  }


  /* Set the pointers for the gauge field on the even and odd checkerboards */
  u_eo[0] = u;
  u_eo[1] = u + GAUGE_SIZE * wilson_p->vol[0];

  /* Set the pointer to the pointer array pointing to the 4 forward half    */
  /* spinors af.                                                            */
  af0 = wilson_p->af[0];
  af1 = wilson_p->af[1];
  af2 = wilson_p->af[2];
  af3 = wilson_p->af[3];

  /* Set the pointer to the pointer array pointing to the 4 backward half   */
  /* spinors ab.                                                            */
  ab0 = wilson_p->ab[0];
  ab1 = wilson_p->ab[1];
  ab2 = wilson_p->ab[2];
  ab3 = wilson_p->ab[3];

  /*------------------------------------------------------------------------*/
  /* Spin project psi with (1-sign*gamma) into the 4 forward half           */
  /* spinors af[4]                                                          */
  /* Spin project psi with (1+sign*gamma) and do the color multiplication   */
  /* of the resulting half spinors with the complex conjugate gauge field.  */
  /* Put the result into the 4 backward half spinors ab[4].                 */
  /*------------------------------------------------------------------------*/
  start_time = rts_get_timebase();
  wfm_cmat_two_spproj(af0, af1, af2, af3, ab0, ab1, ab2, ab3, u_eo[cb], psi, sign, wilson_p, cb);
  stop_time = rts_get_timebase();
  cmat_two_spproj_time = stop_time - start_time;

  /*------------------------------------------------------------------------*/
  /* Communicate                                                            */
  /*------------------------------------------------------------------------*/
  start_time = rts_get_timebase();
  wfm_comm();
  stop_time = rts_get_timebase();
  comm_time = stop_time - start_time;

  /*------------------------------------------------------------------------*/
  /* Do the color multiplication for each forward half spinor af with the   */
  /* gauge field. Then expand the half spinors to full spinors (trick)      */
  /* add them and store their sum into the full spinor chi.                 */
  /*------------------------------------------------------------------------*/
  start_time = rts_get_timebase();
  wfm_mat_trick(chi, u_eo[cbn], af0, af1, af2, af3, sign, wilson_p, cb);
  stop_time = rts_get_timebase();
  mat_trick_time = stop_time - start_time;

  /*------------------------------------------------------------------------*/
  /* Expand the backward half spinors ab to full spinors (trick) add them   */
  /*and then add their sum into the full spinor chi.                        */
  /*------------------------------------------------------------------------*/
  start_time = rts_get_timebase();
  wfm_trick(chi, ab0, ab1, ab2, ab3, sign, wilson_p, cb);
  stop_time = rts_get_timebase();
  trick_time = stop_time - start_time;


  /*------------------------------------------------------------------------*/
  /* Restore double hummer floating point registers                         */
  /*------------------------------------------------------------------------*/
  QuadLoad(wfm_reg  ,  0);
  QuadLoad(wfm_reg+2,  1);
  QuadLoad(wfm_reg+4,  2);
  QuadLoad(wfm_reg+6,  3);
  QuadLoad(wfm_reg+8,  4);
  QuadLoad(wfm_reg+10, 5);
  QuadLoad(wfm_reg+12, 6);
  QuadLoad(wfm_reg+14, 7);
  QuadLoad(wfm_reg+16, 8);
  QuadLoad(wfm_reg+18, 9);  
  QuadLoad(wfm_reg+20, 10);
  QuadLoad(wfm_reg+22, 11);
  QuadLoad(wfm_reg+24, 12);
  QuadLoad(wfm_reg+26, 13);
  QuadLoad(wfm_reg+28, 14);
  QuadLoad(wfm_reg+30, 15);
  QuadLoad(wfm_reg+32, 16);
  QuadLoad(wfm_reg+34, 17);
  QuadLoad(wfm_reg+36, 18);
  QuadLoad(wfm_reg+38, 19);
  QuadLoad(wfm_reg+40, 20);
  QuadLoad(wfm_reg+42, 21);
  QuadLoad(wfm_reg+44, 22);
  QuadLoad(wfm_reg+46, 23);
  QuadLoad(wfm_reg+48, 24);
  QuadLoad(wfm_reg+50, 25);
  QuadLoad(wfm_reg+52, 26);
  QuadLoad(wfm_reg+54, 27);
  QuadLoad(wfm_reg+56, 28);
  QuadLoad(wfm_reg+58, 29);
  QuadLoad(wfm_reg+60, 30);
  QuadLoad(wfm_reg+62, 31);





  stop_time_dslash = rts_get_timebase();
  total_time +=dclock();


  /*------------------------------------------------------------------------*/
  /* ??? timing calc and prints                                             */
  /*------------------------------------------------------------------------*/

  dslash_time   = stop_time_dslash - start_time_dslash;
  dslash_time_i = dslash_time;
  cb_vol = wilson_p->vol[0];
  dslash_perf = dslash_time_i;
  dslash_perf = 1.0 / dslash_perf;
  dslash_perf = dslash_perf * cb_vol * 324 *100;

  DiracOp::CGflops += 1320*cb_vol;
#if 1
  if( (count%10000) == 20 && !UniqueID() ){
    printf("CMAT_TWO SPPROJ   TIME IN PCYCLES PER SITE = %llu\n", cmat_two_spproj_time/cb_vol);
    printf("COMM              TIME IN PCYCLES PER SITE = %llu\n", comm_time/cb_vol);
    printf("MAT_TRICK         TIME IN PCYCLES PER SITE = %llu\n", mat_trick_time/cb_vol);
    printf("TRICK             TIME IN PCYCLES PER SITE = %llu\n", trick_time/cb_vol);
    printf("DSLASH_eo         TIME IN PCYCLES PER SITE = %llu\n", dslash_time/cb_vol);
    printf("DSLASH_eo         PERFORMANCE     = %3.1f\%\n",   dslash_perf);
    print_flops("",fname,1320*cb_vol,total_time);
  }
#endif
  count++;

  /*------------------------------------------------------------------------*/

}
void wilson_dslash_two(Float *chi0, Float *chi1,
                   Float *u,
                   Float *psi0, Float *psi1,
                   int cb0, int cb1,
                   int dag,
                   Wilson *wp)
{
  wilson_dslash(chi0,u,psi0,cb0,dag,wp);
  wilson_dslash(chi1,u,psi1,cb1,dag,wp);
}

CPS_END_NAMESPACE
