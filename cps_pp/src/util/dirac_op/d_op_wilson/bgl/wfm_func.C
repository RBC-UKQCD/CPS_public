#include <config.h>
CPS_START_NAMESPACE
/***************************************************************************/
/*                                                                         */
/* Various optimized functions for the wilson operator                     */
/* The "in", out vectors are half-spinnors (12 components).                */
/* The "inf", outf are full-spinors (24 components).                       */
/*                                                                         */
/* The following registers are expected to be loaded as follows:           */
/* register 31 = 1                                                         */
/* register 30 = 0                                                         */
/* register 28 = sign for spproj, trick                                    */
/* register 29 = sign for cmat_spproj, mat_trick                           */
/*                                                                         */
/* These registers are noit altered by these routines.                     */
/***************************************************************************/

CPS_END_NAMESPACE
#include <util/data_types.h>
#include <util/vector.h>
#include <util/wilson.h>
#include <util/error.h>
#include <sys/bgl/bgl_sys_all.h>
#include <sys/bgl/bgl_sys.h>
#include <comms/scu.h>
CPS_START_NAMESPACE


#define U_P(r,row,col,d)  (u+(r+2*(row+3*(col+3*d))))

#define IN_P(r,c,s)       (in   +(r+2*(c+3*s)))
#define INF_P(r,c,s)      (inf  +(r+2*(c+3*s)))

#define OUT_P(r,c,s)      (out  +(r+2*(c+3*s)))
#define OUTF_P(r,c,s)     (outf +(r+2*(c+3*s)))

#define IN0_P(r,c,s)      (in0  +(r+2*(c+3*s)))
#define IN1_P(r,c,s)      (in1  +(r+2*(c+3*s)))
#define IN2_P(r,c,s)      (in2  +(r+2*(c+3*s)))
#define IN3_P(r,c,s)      (in3  +(r+2*(c+3*s)))

#define WFM_TMP0_P(r,c,s) (wfm_tmp0  +(r+2*(c+3*s)))
#define WFM_TMP1_P(r,c,s) (wfm_tmp1  +(r+2*(c+3*s)))
#define WFM_TMP2_P(r,c,s) (wfm_tmp2  +(r+2*(c+3*s)))
#define WFM_TMP3_P(r,c,s) (wfm_tmp3  +(r+2*(c+3*s)))

#define OUT0_P(r,c,s)     (out0 +(r+2*(c+3*s)))
#define OUT1_P(r,c,s)     (out1 +(r+2*(c+3*s)))
#define OUT2_P(r,c,s)     (out2 +(r+2*(c+3*s)))
#define OUT3_P(r,c,s)     (out3 +(r+2*(c+3*s)))



#define U(r,row,col,d)  *(u+(r+2*(row+3*(col+3*d))))

#define IN(r,c,s)       *(in   +(r+2*(c+3*s)))
#define INF(r,c,s)      *(inf  +(r+2*(c+3*s)))

#define OUT(r,c,s)      *(out  +(r+2*(c+3*s)))
#define OUTF(r,c,s)     *(outf +(r+2*(c+3*s)))

#define IN0(r,c,s)      *(in0  +(r+2*(c+3*s)))
#define IN1(r,c,s)      *(in1  +(r+2*(c+3*s)))
#define IN2(r,c,s)      *(in2  +(r+2*(c+3*s)))
#define IN3(r,c,s)      *(in3  +(r+2*(c+3*s)))

#define OUT0(r,c,s)     *(out0 +(r+2*(c+3*s)))
#define OUT1(r,c,s)     *(out1 +(r+2*(c+3*s)))
#define OUT2(r,c,s)     *(out2 +(r+2*(c+3*s)))
#define OUT3(r,c,s)     *(out3 +(r+2*(c+3*s)))


/***********************************************************************/
/* Color multiply the half spinors "in" with U_mu. Result in "out".    */
/* This routine prefetches the next mu of the gaufe field and the      */
/* half spinor inp using DCBT                                          */
/***********************************************************************/
void wilson_dslash_csmat(double *out, double *u, double *in, double *inp, int mu)
{
  register int o = 16;
  double *u_mu = U_P(0,0,0,mu);
  double *u_mut; // Prefetch pointer
  double *inpt;  // Prefetch pointer


  /* Interleaved loads,stores and calculations */
  QuadLoad (u_mu, 0);
  QuadLoadU(u_mu, 1, o);
  QuadLoadU(u_mu, 2, o);


  QuadLoad (in, 10);
  QuadLoadU(in, 11, o);
  QuadLoadU(in, 12, o);
  QuadLoadU(in, 13, o);

  asm volatile ( "fxcpmadd 22, 0, 10, 30");
  asm volatile ( "fxcpmadd 23, 1, 10, 30");
  asm volatile ( "fxcpmadd 24, 2, 10, 30");
  asm volatile ( "fxcpmadd 25, 0, 13, 30");
  QuadLoadU(in, 14, o);
  asm volatile ( "fxcpmadd 26, 1, 13, 30");
  asm volatile ( "fxcpmadd 27, 2, 13, 30");
  QuadLoadU(in, 15, o);

  asm volatile ( "fxcxnpma 16, 0, 10, 22");
  asm volatile ( "fxcxnpma 17, 1, 10, 23");
  QuadLoadU(u_mu, 3, o);
  asm volatile ( "fxcxnpma 18, 2, 10, 24");
  asm volatile ( "fxcxnpma 19, 0, 13, 25");
  QuadLoadU(u_mu, 4, o);
  asm volatile ( "fxcxnpma 20, 1, 13, 26");
  asm volatile ( "fxcxnpma 21, 2, 13, 27");
  QuadLoadU(u_mu, 5, o);
  
  asm volatile ( "fxcpmadd 22, 3, 11, 16");
  asm volatile ( "fxcpmadd 23, 4, 11, 17");
  QuadLoadU(u_mu, 6, o);
  asm volatile ( "fxcpmadd 24, 5, 11, 18");
  asm volatile ( "fxcpmadd 25, 3, 14, 19");
  QuadLoadU(u_mu, 7, o);
  asm volatile ( "fxcpmadd 26, 4, 14, 20");
  asm volatile ( "fxcpmadd 27, 5, 14, 21");
  QuadLoadU(u_mu, 8, o);

  u_mut = u_mu+18;
  asm volatile ( "fxcxnpma 16, 3, 11, 22");
  asm volatile ( "fxcxnpma 17, 4, 11, 23");
  DCBT(u_mut);
  u_mut = u_mu+22;
  asm volatile ( "fxcxnpma 18, 5, 11, 24");
  DCBT(u_mut);
  u_mut = u_mu+26;
  asm volatile ( "fxcxnpma 19, 3, 14, 25");
  DCBT(u_mut);
  u_mut = u_mu+30;
  asm volatile ( "fxcxnpma 20, 4, 14, 26");
  DCBT(u_mut);
  u_mut = u_mu+34;
  asm volatile ( "fxcxnpma 21, 5, 14, 27");
  DCBT(u_mut);

  asm volatile ( "fxcpmadd 22, 6, 12, 16");
  DCBT(inp);
  inpt = inp+4;
  asm volatile ( "fxcpmadd 23, 7, 12, 17");
  DCBT(inpt);
  inpt = inp+8;
  asm volatile ( "fxcpmadd 24, 8, 12, 18");
  DCBT(inpt);
  inpt = inp+12;
  asm volatile ( "fxcpmadd 25, 6, 15, 19");
  DCBT(inpt);
  inpt = inp+16;
  asm volatile ( "fxcpmadd 26, 7, 15, 20");
  DCBT(inpt);
  asm volatile ( "fxcpmadd 27, 8, 15, 21");
 
  asm volatile ( "fxcxnpma 16, 6, 12, 22");
  asm volatile ( "fxcxnpma 17, 7, 12, 23");
  asm volatile ( "fxcxnpma 18, 8, 12, 24");
  QuadStore(out, 16);
  asm volatile ( "fxcxnpma 19, 6, 15, 25");
  QuadStoreU(out, 17, o);
  asm volatile ( "fxcxnpma 20, 7, 15, 26");
  QuadStoreU(out, 18, o);
  asm volatile ( "fxcxnpma 21, 8, 15, 27");

  QuadStoreU(out, 19, o);
  QuadStoreU(out, 20, o);
  QuadStoreU(out, 21, o);

}



/***********************************************************************/
/* Color multiply the half spinors "in" with U_mu^dag. Result in "out".*/ 
/* This routine prefetches the next mu of the gaufe field using DCBT.  */
/***********************************************************************/
void wilson_dslash_csmatdag(double *out, double *u, double *in, int mu)
{
  double *u_mu = U_P(0,0,0,mu);
  double *u_mut; // Prefetch pointer

  /* Initial loads of U and the half spinor in */
  QuadLoad(U_P(0,0,0,mu), 0);
  QuadLoad(U_P(0,0,1,mu), 1);
  QuadLoad(U_P(0,0,2,mu), 2);

  QuadLoad(IN_P(0,0,0), 10);
  QuadLoad(IN_P(0,0,1), 13);
  QuadLoad(IN_P(0,1,0), 11);


  /* Interleaved loads,stores and calculations */
  asm volatile ( "fxcpmadd 22, 0, 10, 30");
  QuadLoad(IN_P(0,2,0), 12);
  asm volatile ( "fxcpmadd 23, 1, 10, 30");
  QuadLoad(IN_P(0,1,1), 14);
  asm volatile ( "fxcpmadd 24, 2, 10, 30");
  QuadLoad(IN_P(0,2,1), 15);
  asm volatile ( "fxcpmadd 25, 0, 13, 30");
  QuadLoad(U_P(0,1,0,mu), 3);
  asm volatile ( "fxcpmadd 26, 1, 13, 30");
  asm volatile ( "fxcpmadd 27, 2, 13, 30");
  QuadLoad(U_P(0,1,1,mu), 4);

  asm volatile ( "fxcxnsma 16, 0, 10, 22");
  asm volatile ( "fxcxnsma 17, 1, 10, 23");
  asm volatile ( "fxcxnsma 18, 2, 10, 24");
  QuadLoad(U_P(0,1,2,mu), 5);
  asm volatile ( "fxcxnsma 19, 0, 13, 25");
  asm volatile ( "fxcxnsma 20, 1, 13, 26");
  asm volatile ( "fxcxnsma 21, 2, 13, 27");
  QuadLoad(U_P(0,2,0,mu), 6);
  
  asm volatile ( "fxcpmadd 22, 3, 11, 16");
  asm volatile ( "fxcpmadd 23, 4, 11, 17");
  asm volatile ( "fxcpmadd 24, 5, 11, 18");
  QuadLoad(U_P(0,2,1,mu), 7);
  asm volatile ( "fxcpmadd 25, 3, 14, 19");
  asm volatile ( "fxcpmadd 26, 4, 14, 20");
  asm volatile ( "fxcpmadd 27, 5, 14, 21");
  QuadLoad(U_P(0,2,2,mu), 8);

  u_mut = u_mu+18;
  asm volatile ( "fxcxnsma 16, 3, 11, 22");
  asm volatile ( "fxcxnsma 17, 4, 11, 23");
  DCBT(u_mut);
  u_mut = u_mu+22;
  asm volatile ( "fxcxnsma 18, 5, 11, 24");
  DCBT(u_mut);
  u_mut = u_mu+26;
  asm volatile ( "fxcxnsma 19, 3, 14, 25");
  DCBT(u_mut);
  u_mut = u_mu+30;
  asm volatile ( "fxcxnsma 20, 4, 14, 26");
  DCBT(u_mut);
  u_mut = u_mu+34;
  asm volatile ( "fxcxnsma 21, 5, 14, 27");
  DCBT(u_mut);

  asm volatile ( "fxcpmadd 22, 6, 12, 16");
  asm volatile ( "fxcpmadd 23, 7, 12, 17");
  asm volatile ( "fxcpmadd 24, 8, 12, 18");
  asm volatile ( "fxcpmadd 25, 6, 15, 19");
  asm volatile ( "fxcpmadd 26, 7, 15, 20");
  asm volatile ( "fxcpmadd 27, 8, 15, 21");
 
  asm volatile ( "fxcxnsma 16, 6, 12, 22");
  asm volatile ( "fxcxnsma 17, 7, 12, 23");
  asm volatile ( "fxcxnsma 18, 8, 12, 24");
  QuadStore(OUT_P(0,0,0), 16);
  asm volatile ( "fxcxnsma 19, 6, 15, 25");
  QuadStore(OUT_P(0,1,0), 17);
  asm volatile ( "fxcxnsma 20, 7, 15, 26");
  QuadStore(OUT_P(0,2,0), 18);
  asm volatile ( "fxcxnsma 21, 8, 15, 27");

  /* final stores of the half spinor out */
  QuadStore(OUT_P(0,0,1), 19);
  QuadStore(OUT_P(0,1,1), 20);
  QuadStore(OUT_P(0,2,1), 21);

}



/***********************************************************************/
/* Spin project with [1 + sign * gamma] the full spinor "inf" onto the */
/* four half spinors "out0", "out1", "out2", "out3".                   */
/***********************************************************************/
/*
void wilson_dslash_spproj(double *out0, 
			  double *out1, 
			  double *out2, 
			  double *out3, 
			  double *inf)
{ 
  register int o = 16;

  QuadLoad( inf, 0);
  QuadLoadU(inf, 4, o);
  QuadLoadU(inf, 8, o);
  QuadLoadU(inf, 1, o);
  QuadLoadU(inf, 5, o);
  QuadLoadU(inf, 9, o);
  QuadLoadU(inf, 2, o);
  QuadLoadU(inf, 6, o);
  QuadLoadU(inf, 10, o);
  QuadLoadU(inf, 3, o);
  QuadLoadU(inf, 7, o);
  QuadLoadU(inf, 11, o);

  asm volatile ( "fxcxnpma 20, 28, 3, 0");
  asm volatile ( "fxcxnpma 21, 28, 2, 1");
  asm volatile ( "fpnmsub  22, 28, 3, 0");
  asm volatile ( "fpmadd   23, 28, 2, 1");
  asm volatile ( "fxcxnpma 24, 28, 2, 0");
  asm volatile ( "fxcxnsma 25, 28, 3, 1");
  QuadStore(OUT0_P(0,0,0), 20);
  asm volatile ( "fpmadd   26, 28, 2, 0");
  QuadStore(OUT0_P(0,0,1), 21);
  asm volatile ( "fpmadd   27, 28, 3, 1");
  QuadStore(OUT1_P(0,0,0), 22);
  QuadStore(OUT1_P(0,0,1), 23);
  QuadStore(OUT2_P(0,0,0), 24);
  QuadStore(OUT2_P(0,0,1), 25);
  QuadStore(OUT3_P(0,0,0), 26);
  QuadStore(OUT3_P(0,0,1), 27);

  asm volatile ( "fxcxnpma 20, 28, 7, 4");
  asm volatile ( "fxcxnpma 21, 28, 6, 5");
  asm volatile ( "fpnmsub  22, 28, 7, 4");
  asm volatile ( "fpmadd   23, 28, 6, 5");
  asm volatile ( "fxcxnpma 24, 28, 6, 4");
  asm volatile ( "fxcxnsma 25, 28, 7, 5");
  QuadStore(OUT0_P(0,1,0), 20);
  asm volatile ( "fpmadd   26, 28, 6, 4");
  QuadStore(OUT0_P(0,1,1), 21);
  asm volatile ( "fpmadd   27, 28, 7, 5");
  QuadStore(OUT1_P(0,1,0), 22);
  QuadStore(OUT1_P(0,1,1), 23);
  QuadStore(OUT2_P(0,1,0), 24);
  QuadStore(OUT2_P(0,1,1), 25);
  QuadStore(OUT3_P(0,1,0), 26);
  QuadStore(OUT3_P(0,1,1), 27);


  asm volatile ( "fxcxnpma 20, 28, 11, 8");
  asm volatile ( "fxcxnpma 21, 28, 10, 9");
  asm volatile ( "fpnmsub  22, 28, 11, 8");
  asm volatile ( "fpmadd   23, 28, 10, 9");
  asm volatile ( "fxcxnpma 24, 28, 10, 8");
  asm volatile ( "fxcxnsma 25, 28, 11, 9");
  QuadStore(OUT0_P(0,2,0), 20);
  asm volatile ( "fpmadd   26, 28, 10, 8");
  QuadStore(OUT0_P(0,2,1), 21);
  asm volatile ( "fpmadd   27, 28, 11, 9");
  QuadStore(OUT1_P(0,2,0), 22);
  QuadStore(OUT1_P(0,2,1), 23);
  QuadStore(OUT2_P(0,2,0), 24);
  QuadStore(OUT2_P(0,2,1), 25);
  QuadStore(OUT3_P(0,2,0), 26);
  QuadStore(OUT3_P(0,2,1), 27);

}
*/




/***********************************************************************/
/* Spin project with [1 + sign * gamma] the full spinor "inf" onto the */
/* four half spinors "out0", "out1", "out2", "out3".                   */
/***********************************************************************/
void wilson_dslash_spproj(double *out0, 
			  double *out1, 
			  double *out2, 
			  double *out3, 
			  double *inf)
{ 
  register int o = 16;
  
  /*-----------------------------------------------*/
  QuadLoad( inf, 0);
  QuadLoadU(inf, 4, o);
  QuadLoadU(inf, 8, o);
  QuadLoadU(inf, 1, o);
  QuadLoadU(inf, 5, o);
  QuadLoadU(inf, 9, o);
  QuadLoadU(inf, 2, o);
  QuadLoadU(inf, 6, o);
  QuadLoadU(inf, 10, o);
  QuadLoadU(inf, 3, o);
  QuadLoadU(inf, 7, o);
  QuadLoadU(inf, 11, o);

  /*-----------------------------------------------*/
  asm volatile ( "fxcxnpma 20, 28, 3, 0");
  asm volatile ( "fxcxnpma 22, 28, 7, 4");
  asm volatile ( "fxcxnpma 24, 28, 11, 8");
  asm volatile ( "fxcxnpma 21, 28, 2, 1");
  asm volatile ( "fxcxnpma 23, 28, 6, 5");
  asm volatile ( "fxcxnpma 25, 28, 10, 9");

  /*-----------------------------------------------*/
  QuadStore( out0, 20);
  asm volatile ( "fpnmsub  14, 28, 3, 0");
  QuadStoreU(out0, 22, o);
  asm volatile ( "fpnmsub  16, 28, 7, 4");
  QuadStoreU(out0, 24, o);
  asm volatile ( "fpnmsub  18, 28, 11, 8");
  QuadStoreU(out0, 21, o);
  asm volatile ( "fpmadd   15, 28, 2, 1");
  QuadStoreU(out0, 23, o);
  asm volatile ( "fpmadd   17, 28, 6, 5");
  QuadStoreU(out0, 25, o);
  asm volatile ( "fpmadd   19, 28, 10, 9");


  /*-----------------------------------------------*/
  QuadStore( out1, 14);
  asm volatile ( "fxcxnpma 20, 28, 2, 0");
  QuadStoreU(out1, 16, o);
  asm volatile ( "fxcxnpma 22, 28, 6, 4");
  QuadStoreU(out1, 18, o);
  asm volatile ( "fxcxnpma 24, 28, 10, 8");
  QuadStoreU(out1, 15, o);
  asm volatile ( "fxcxnsma 21, 28, 3, 1");
  QuadStoreU(out1, 17, o);
  asm volatile ( "fxcxnsma 23, 28, 7, 5");
  QuadStoreU(out1, 19, o);
  asm volatile ( "fxcxnsma 25, 28, 11, 9");


  /*-----------------------------------------------*/
  QuadStore( out2, 20);
  asm volatile ( "fpmadd   14, 28, 2, 0");
  QuadStoreU(out2, 22, o);
  asm volatile ( "fpmadd   16, 28, 6, 4");
  QuadStoreU(out2, 24, o);
  asm volatile ( "fpmadd   18, 28, 10, 8");
  QuadStoreU(out2, 21, o);
  asm volatile ( "fpmadd   15, 28, 3, 1");
  QuadStoreU(out2, 23, o);
  asm volatile ( "fpmadd   17, 28, 7, 5");
  QuadStoreU(out2, 25, o);
  asm volatile ( "fpmadd   19, 28, 11, 9");

  QuadStore( out3, 14);
  QuadStoreU(out3, 16, o);
  QuadStoreU(out3, 18, o);
  QuadStoreU(out3, 15, o);
  QuadStoreU(out3, 17, o);
  QuadStoreU(out3, 19, o);

}



/***********************************************************************/
/* Expand the half spinors "in" to a full spinor "outf" (trick), where */
/* the spin projection was done with [1+sign*gamma].                   */
/* If acc=0 store the result into the full spinor "outf".              */
/* If acc=1 accumulate the result into the full spinor "outf".         */
/* This routine prefetches the next in0 half spinor using DCBT.        */
/***********************************************************************/
void wilson_dslash_trick(double *outf, double *in0, double *in1, double *in2, double *in3, 
			 double sign, int accum)
{

  register int o0;
  register int o1 = 48;
  register int o2 = 96;
  register int o3 = 144;

 /* c=0 */
  QuadLoad(outf, 0);
  QuadLoadNU(outf, 1, o1);
  QuadLoadNU(outf, 2, o2);
  QuadLoadNU(outf, 3, o3);
  QuadLoad(in0, 10);
  QuadLoadNU(in0, 11, o1);
  QuadLoad(in1, 12);
  asm volatile ( "fpmadd   20,31,10,0");
  QuadLoadNU(in1, 13, o1);
  asm volatile ( "fpmadd   21,31,11,1");
  QuadLoad(in2, 14);
  asm volatile ( "fxcxnsma 22,29,11,2");
  QuadLoadNU(in2, 15, o1);
  asm volatile ( "fxcxnsma 23,29,10,3");
  asm volatile ( "fpmadd   0,31,12,20");
  asm volatile ( "fpmadd   1,31,13,21");
  asm volatile ( "fpmadd   2,29,13,22");
  asm volatile ( "fpnmsub  3,29,12,23");
  asm volatile ( "fpmadd   20,31,14,0");
  QuadLoad(in3, 16);
  asm volatile ( "fpmadd   21,31,15,1");
  QuadLoadNU(in3, 17, o1);
  asm volatile ( "fxcxnsma 22,29,14,2");
  DCBT(in0+12);
  asm volatile ( "fxcxnpma 23,29,15,3");
  DCBT(in0+16);
  asm volatile ( "fpmadd   0,31,16,20");
  asm volatile ( "fpmadd   1,31,17,21");
  QuadStore(outf, 0);
  DCBT(in0+20);
  asm volatile ( "fpmadd   2,29,16,22");
  QuadStoreNU(outf, 1, o1);
  asm volatile ( "fpmadd   3,29,17,23");
  QuadStoreNU(outf, 2, o2);
  QuadStoreNU(outf, 3, o3);
 

  o0 = 16;
  o1 = 64;
  o2 = 112;
  o3 = 160;

  /* c=1 */
  QuadLoadNU(outf, 0, o0);
  QuadLoadNU(outf, 1, o1);
  QuadLoadNU(outf, 2, o2);
  QuadLoadNU(outf, 3, o3);
  QuadLoadNU(in0, 10, o0);
  QuadLoadNU(in0, 11, o1);    
  QuadLoadNU(in1, 12, o0);
  QuadLoadNU(in1, 13, o1);
  QuadLoadNU(in2, 14, o0);
  asm volatile ( "fpmadd   20,31,10,0");
  asm volatile ( "fpmadd   21,31,11,1");
  asm volatile ( "fxcxnsma 22,29,11,2");
  QuadLoadNU(in2, 15, o1);
  asm volatile ( "fxcxnsma 23,29,10,3");
  asm volatile ( "fpmadd   0,31,12,20");
  asm volatile ( "fpmadd   1,31,13,21");
  QuadLoadNU(in3, 16, o0);
  asm volatile ( "fpmadd   2,29,13,22");
  asm volatile ( "fpnmsub  3,29,12,23");
  asm volatile ( "fpmadd   20,31,14,0");
  QuadLoadNU(in3, 17, o1);
  asm volatile ( "fpmadd   21,31,15,1");
  asm volatile ( "fxcxnsma 22,29,14,2");
  asm volatile ( "fxcxnpma 23,29,15,3");
  asm volatile ( "fpmadd   0,31,16,20");
  asm volatile ( "fpmadd   1,31,17,21");
  QuadStoreNU(outf, 0, o0);
  asm volatile ( "fpmadd   2,29,16,22");
  QuadStoreNU(outf, 1, o1);
  asm volatile ( "fpmadd   3,29,17,23");
  QuadStoreNU(outf, 2, o2);
  QuadStoreNU(outf, 3, o3);
   


  o0 = 32;
  o1 = 80;
  o2 = 128;
  o3 = 176;
  QuadLoadNU(outf, 0, o0);
  QuadLoadNU(outf, 1, o1);
  QuadLoadNU(outf, 2, o2);
  QuadLoadNU(outf, 3, o3);
  QuadLoadNU(in0, 10, o0);
  QuadLoadNU(in0, 11, o1);    
  QuadLoadNU(in1, 12, o0);
  QuadLoadNU(in1, 13, o1);
  QuadLoadNU(in2, 14, o0);
  asm volatile ( "fpmadd   20,31,10,0");
  asm volatile ( "fpmadd   21,31,11,1");
  asm volatile ( "fxcxnsma 22,29,11,2");
  QuadLoadNU(in2, 15, o1);
  asm volatile ( "fxcxnsma 23,29,10,3");
  asm volatile ( "fpmadd   0,31,12,20");
  asm volatile ( "fpmadd   1,31,13,21");
  QuadLoadNU(in3, 16, o0);
  asm volatile ( "fpmadd   2,29,13,22");
  asm volatile ( "fpnmsub  3,29,12,23");
  asm volatile ( "fpmadd   20,31,14,0");
  QuadLoadNU(in3, 17, o1);
  asm volatile ( "fpmadd   21,31,15,1");
  asm volatile ( "fxcxnsma 22,29,14,2");
  asm volatile ( "fxcxnpma 23,29,15,3");
  asm volatile ( "fpmadd   0,31,16,20");
  asm volatile ( "fpmadd   1,31,17,21");
  asm volatile ( "fpmadd   2,29,16,22");
  asm volatile ( "fpmadd   3,29,17,23");
  QuadStoreNU(outf, 0, o0);
  QuadStoreNU(outf, 1, o1);
  QuadStoreNU(outf, 2, o2);
  QuadStoreNU(outf, 3, o3);
  
}



/***********************************************************************/
/* Spin project with [1 + sign * gamma] the full spinor "inf" onto     */
/* four half spinors. Color multiply each of the half spinors mu with  */
/* u_mu. Store the result in out0, out1, ou2, out3.                    */
/***********************************************************************/
void wilson_dslash_cmat_spproj(double *out0, 
			       double *out1, 
			       double *out2, 
			       double *out3, 
			       double *u, 
			       double *inf)
{

  /*********************************************************************/
  /* mu=0                                                              */  
  /*********************************************************************/
  QuadLoad(INF_P(0,0,0), 9);
  QuadLoad(INF_P(0,0,1), 10);
  QuadLoad(INF_P(0,0,2), 11);
  QuadLoad(INF_P(0,0,3), 12);

  QuadLoad(U_P(0,0,0,0), 0);
  QuadLoad(U_P(0,0,1,0), 1);
  QuadLoad(U_P(0,0,2,0), 2);

  asm volatile ( "fxcxnpma 13, 29, 12, 9");
  asm volatile ( "fxcxnpma 14, 29, 11, 10");
  QuadLoad(INF_P(0,1,0), 9);
  QuadLoad(INF_P(0,1,1), 10);
  QuadLoad(INF_P(0,1,2), 11);
  QuadLoad(INF_P(0,1,3), 12);

  asm volatile ( "fxcpmadd 22, 0, 13, 30");
  asm volatile ( "fxcpmadd 23, 1, 13, 30");
  asm volatile ( "fxcpmadd 24, 2, 13, 30");
  asm volatile ( "fxcpmadd 25, 0, 14, 30");
  asm volatile ( "fxcpmadd 26, 1, 14, 30");
  asm volatile ( "fxcpmadd 27, 2, 14, 30");
  QuadLoad(U_P(0,1,0,0), 3);
  QuadLoad(U_P(0,1,1,0), 4);
  QuadLoad(U_P(0,1,2,0), 5);

  asm volatile ( "fxcxnsma 16, 0, 13, 22");
  asm volatile ( "fxcxnsma 17, 1, 13, 23");
  asm volatile ( "fxcxnsma 18, 2, 13, 24");
  asm volatile ( "fxcxnsma 19, 0, 14, 25");
  asm volatile ( "fxcxnsma 20, 1, 14, 26");
  asm volatile ( "fxcxnsma 21, 2, 14, 27");

  asm volatile ( "fxcxnpma 13, 29, 12, 9");
  asm volatile ( "fxcxnpma 14, 29, 11, 10");
  QuadLoad(INF_P(0,2,0), 9);
  QuadLoad(INF_P(0,2,1), 10);
  QuadLoad(INF_P(0,2,2), 11);
  QuadLoad(INF_P(0,2,3), 12);

  asm volatile ( "fxcpmadd 22, 3, 13, 16");
  asm volatile ( "fxcpmadd 23, 4, 13, 17");
  asm volatile ( "fxcpmadd 24, 5, 13, 18");
  asm volatile ( "fxcpmadd 25, 3, 14, 19");
  asm volatile ( "fxcpmadd 26, 4, 14, 20");
  asm volatile ( "fxcpmadd 27, 5, 14, 21");
  QuadLoad(U_P(0,2,0,0), 6);
  QuadLoad(U_P(0,2,1,0), 7);
  QuadLoad(U_P(0,2,2,0), 8);
  

  asm volatile ( "fxcxnsma 16, 3, 13, 22");
  asm volatile ( "fxcxnsma 17, 4, 13, 23");
  asm volatile ( "fxcxnsma 18, 5, 13, 24");
  asm volatile ( "fxcxnsma 19, 3, 14, 25");
  asm volatile ( "fxcxnsma 20, 4, 14, 26");
  asm volatile ( "fxcxnsma 21, 5, 14, 27");

  asm volatile ( "fxcxnpma 13, 29, 12, 9");
  asm volatile ( "fxcxnpma 14, 29, 11, 10");
  QuadLoad(INF_P(0,0,0), 9);
  QuadLoad(INF_P(0,0,1), 10);
  QuadLoad(INF_P(0,0,2), 11);
  QuadLoad(INF_P(0,0,3), 12);

  asm volatile ( "fxcpmadd 22, 6, 13, 16");
  asm volatile ( "fxcpmadd 23, 7, 13, 17");
  asm volatile ( "fxcpmadd 24, 8, 13, 18");
  asm volatile ( "fxcpmadd 25, 6, 14, 19");
  asm volatile ( "fxcpmadd 26, 7, 14, 20");
  asm volatile ( "fxcpmadd 27, 8, 14, 21");
  QuadLoad(U_P(0,0,0,1), 0);
  QuadLoad(U_P(0,0,1,1), 1);
  QuadLoad(U_P(0,0,2,1), 2);
 
  asm volatile ( "fxcxnsma 16, 6, 13, 22");
  asm volatile ( "fxcxnsma 17, 7, 13, 23");
  asm volatile ( "fxcxnsma 18, 8, 13, 24");
  asm volatile ( "fxcxnsma 19, 6, 14, 25");
  asm volatile ( "fxcxnsma 20, 7, 14, 26");
  asm volatile ( "fxcxnsma 21, 8, 14, 27");

  QuadStore(OUT0_P(0,0,0), 16);
  QuadStore(OUT0_P(0,1,0), 17);
  QuadStore(OUT0_P(0,2,0), 18);
  QuadStore(OUT0_P(0,0,1), 19);
  QuadStore(OUT0_P(0,1,1), 20);
  QuadStore(OUT0_P(0,2,1), 21);


  /*********************************************************************/
  /* mu=1                                                              */  
  /*********************************************************************/


  asm volatile ( "fpnmsub  13, 29, 12, 9");
  asm volatile ( "fpmadd   14, 29, 11, 10");
  QuadLoad(INF_P(0,1,0), 9);
  QuadLoad(INF_P(0,1,1), 10);
  QuadLoad(INF_P(0,1,2), 11);
  QuadLoad(INF_P(0,1,3), 12);

  asm volatile ( "fxcpmadd 22, 0, 13, 30");
  asm volatile ( "fxcpmadd 23, 1, 13, 30");
  asm volatile ( "fxcpmadd 24, 2, 13, 30");
  asm volatile ( "fxcpmadd 25, 0, 14, 30");
  asm volatile ( "fxcpmadd 26, 1, 14, 30");
  asm volatile ( "fxcpmadd 27, 2, 14, 30");
  QuadLoad(U_P(0,1,0,1), 3);
  QuadLoad(U_P(0,1,1,1), 4);
  QuadLoad(U_P(0,1,2,1), 5);

  asm volatile ( "fxcxnsma 16, 0, 13, 22");
  asm volatile ( "fxcxnsma 17, 1, 13, 23");
  asm volatile ( "fxcxnsma 18, 2, 13, 24");
  asm volatile ( "fxcxnsma 19, 0, 14, 25");
  asm volatile ( "fxcxnsma 20, 1, 14, 26");
  asm volatile ( "fxcxnsma 21, 2, 14, 27");
  

  asm volatile ( "fpnmsub  13, 29, 12, 9");
  asm volatile ( "fpmadd   14, 29, 11, 10");
  QuadLoad(INF_P(0,2,0), 9);
  QuadLoad(INF_P(0,2,1), 10);
  QuadLoad(INF_P(0,2,2), 11);
  QuadLoad(INF_P(0,2,3), 12);

  asm volatile ( "fxcpmadd 22, 3, 13, 16");
  asm volatile ( "fxcpmadd 23, 4, 13, 17");
  asm volatile ( "fxcpmadd 24, 5, 13, 18");
  asm volatile ( "fxcpmadd 25, 3, 14, 19");
  asm volatile ( "fxcpmadd 26, 4, 14, 20");
  asm volatile ( "fxcpmadd 27, 5, 14, 21");
  QuadLoad(U_P(0,2,0,1), 6);
  QuadLoad(U_P(0,2,1,1), 7);
  QuadLoad(U_P(0,2,2,1), 8);

  asm volatile ( "fxcxnsma 16, 3, 13, 22");
  asm volatile ( "fxcxnsma 17, 4, 13, 23");
  asm volatile ( "fxcxnsma 18, 5, 13, 24");
  asm volatile ( "fxcxnsma 19, 3, 14, 25");
  asm volatile ( "fxcxnsma 20, 4, 14, 26");
  asm volatile ( "fxcxnsma 21, 5, 14, 27");


  asm volatile ( "fpnmsub  13, 29, 12, 9");
  asm volatile ( "fpmadd   14, 29, 11, 10");
  QuadLoad(INF_P(0,0,0), 9);
  QuadLoad(INF_P(0,0,1), 10);
  QuadLoad(INF_P(0,0,2), 11);
  QuadLoad(INF_P(0,0,3), 12);

  asm volatile ( "fxcpmadd 22, 6, 13, 16");
  asm volatile ( "fxcpmadd 23, 7, 13, 17");
  asm volatile ( "fxcpmadd 24, 8, 13, 18");
  asm volatile ( "fxcpmadd 25, 6, 14, 19");
  asm volatile ( "fxcpmadd 26, 7, 14, 20");
  asm volatile ( "fxcpmadd 27, 8, 14, 21");
  QuadLoad(U_P(0,0,0,2), 0);
  QuadLoad(U_P(0,0,1,2), 1);
  QuadLoad(U_P(0,0,2,2), 2);
 
  asm volatile ( "fxcxnsma 16, 6, 13, 22");
  asm volatile ( "fxcxnsma 17, 7, 13, 23");
  asm volatile ( "fxcxnsma 18, 8, 13, 24");
  asm volatile ( "fxcxnsma 19, 6, 14, 25");
  asm volatile ( "fxcxnsma 20, 7, 14, 26");
  asm volatile ( "fxcxnsma 21, 8, 14, 27");


  QuadStore(OUT1_P(0,0,0), 16);
  QuadStore(OUT1_P(0,1,0), 17);
  QuadStore(OUT1_P(0,2,0), 18);
  QuadStore(OUT1_P(0,0,1), 19);
  QuadStore(OUT1_P(0,1,1), 20);
  QuadStore(OUT1_P(0,2,1), 21);



  /*********************************************************************/
  /* mu=2                                                              */  
  /*********************************************************************/
  asm volatile ( "fxcxnpma 13, 29, 11, 9");
  asm volatile ( "fxcxnsma 14, 29, 12, 10");
  QuadLoad(INF_P(0,1,0), 9);
  QuadLoad(INF_P(0,1,1), 10);
  QuadLoad(INF_P(0,1,2), 11);
  QuadLoad(INF_P(0,1,3), 12);

  asm volatile ( "fxcpmadd 22, 0, 13, 30");
  asm volatile ( "fxcpmadd 23, 1, 13, 30");
  asm volatile ( "fxcpmadd 24, 2, 13, 30");
  asm volatile ( "fxcpmadd 25, 0, 14, 30");
  asm volatile ( "fxcpmadd 26, 1, 14, 30");
  asm volatile ( "fxcpmadd 27, 2, 14, 30");
  QuadLoad(U_P(0,1,0,2), 3);
  QuadLoad(U_P(0,1,1,2), 4);
  QuadLoad(U_P(0,1,2,2), 5);

  asm volatile ( "fxcxnsma 16, 0, 13, 22");
  asm volatile ( "fxcxnsma 17, 1, 13, 23");
  asm volatile ( "fxcxnsma 18, 2, 13, 24");
  asm volatile ( "fxcxnsma 19, 0, 14, 25");
  asm volatile ( "fxcxnsma 20, 1, 14, 26");
  asm volatile ( "fxcxnsma 21, 2, 14, 27");
  

  asm volatile ( "fxcxnpma 13, 29, 11, 9");
  asm volatile ( "fxcxnsma 14, 29, 12, 10");
  QuadLoad(INF_P(0,2,0), 9);
  QuadLoad(INF_P(0,2,1), 10);
  QuadLoad(INF_P(0,2,2), 11);
  QuadLoad(INF_P(0,2,3), 12);

  asm volatile ( "fxcpmadd 22, 3, 13, 16");
  asm volatile ( "fxcpmadd 23, 4, 13, 17");
  asm volatile ( "fxcpmadd 24, 5, 13, 18");
  asm volatile ( "fxcpmadd 25, 3, 14, 19");
  asm volatile ( "fxcpmadd 26, 4, 14, 20");
  asm volatile ( "fxcpmadd 27, 5, 14, 21");
  QuadLoad(U_P(0,2,0,2), 6);
  QuadLoad(U_P(0,2,1,2), 7);
  QuadLoad(U_P(0,2,2,2), 8);


  asm volatile ( "fxcxnsma 16, 3, 13, 22");
  asm volatile ( "fxcxnsma 17, 4, 13, 23");
  asm volatile ( "fxcxnsma 18, 5, 13, 24");
  asm volatile ( "fxcxnsma 19, 3, 14, 25");
  asm volatile ( "fxcxnsma 20, 4, 14, 26");
  asm volatile ( "fxcxnsma 21, 5, 14, 27");


  asm volatile ( "fxcxnpma 13, 29, 11, 9");
  asm volatile ( "fxcxnsma 14, 29, 12, 10");
  QuadLoad(INF_P(0,0,0), 9);
  QuadLoad(INF_P(0,0,1), 10);
  QuadLoad(INF_P(0,0,2), 11);
  QuadLoad(INF_P(0,0,3), 12);

  asm volatile ( "fxcpmadd 22, 6, 13, 16");
  asm volatile ( "fxcpmadd 23, 7, 13, 17");
  asm volatile ( "fxcpmadd 24, 8, 13, 18");
  asm volatile ( "fxcpmadd 25, 6, 14, 19");
  asm volatile ( "fxcpmadd 26, 7, 14, 20");
  asm volatile ( "fxcpmadd 27, 8, 14, 21");
  QuadLoad(U_P(0,0,0,3), 0);
  QuadLoad(U_P(0,0,1,3), 1);
  QuadLoad(U_P(0,0,2,3), 2);
 
  asm volatile ( "fxcxnsma 16, 6, 13, 22");
  asm volatile ( "fxcxnsma 17, 7, 13, 23");
  asm volatile ( "fxcxnsma 18, 8, 13, 24");
  asm volatile ( "fxcxnsma 19, 6, 14, 25");
  asm volatile ( "fxcxnsma 20, 7, 14, 26");
  asm volatile ( "fxcxnsma 21, 8, 14, 27");

  QuadStore(OUT2_P(0,0,0), 16);
  QuadStore(OUT2_P(0,1,0), 17);
  QuadStore(OUT2_P(0,2,0), 18);
  QuadStore(OUT2_P(0,0,1), 19);
  QuadStore(OUT2_P(0,1,1), 20);
  QuadStore(OUT2_P(0,2,1), 21);


  /*********************************************************************/
  /* mu=3                                                              */  
  /*********************************************************************/
  asm volatile ( "fpmadd   13, 29, 11, 9");
  asm volatile ( "fpmadd   14, 29, 12, 10");
  QuadLoad(INF_P(0,1,0), 9);
  QuadLoad(INF_P(0,1,1), 10);
  QuadLoad(INF_P(0,1,2), 11);
  QuadLoad(INF_P(0,1,3), 12);

  asm volatile ( "fxcpmadd 22, 0, 13, 30");
  asm volatile ( "fxcpmadd 23, 1, 13, 30");
  asm volatile ( "fxcpmadd 24, 2, 13, 30");
  asm volatile ( "fxcpmadd 25, 0, 14, 30");
  asm volatile ( "fxcpmadd 26, 1, 14, 30");
  asm volatile ( "fxcpmadd 27, 2, 14, 30");
  QuadLoad(U_P(0,1,0,3), 3);
  QuadLoad(U_P(0,1,1,3), 4);
  QuadLoad(U_P(0,1,2,3), 5);

  asm volatile ( "fxcxnsma 16, 0, 13, 22");
  asm volatile ( "fxcxnsma 17, 1, 13, 23");
  asm volatile ( "fxcxnsma 18, 2, 13, 24");
  asm volatile ( "fxcxnsma 19, 0, 14, 25");
  asm volatile ( "fxcxnsma 20, 1, 14, 26");
  asm volatile ( "fxcxnsma 21, 2, 14, 27");
  

  asm volatile ( "fpmadd   13, 29, 11, 9");
  asm volatile ( "fpmadd   14, 29, 12, 10");
  QuadLoad(INF_P(0,2,0), 9);
  QuadLoad(INF_P(0,2,1), 10);
  QuadLoad(INF_P(0,2,2), 11);
  QuadLoad(INF_P(0,2,3), 12);

  asm volatile ( "fxcpmadd 22, 3, 13, 16");
  asm volatile ( "fxcpmadd 23, 4, 13, 17");
  asm volatile ( "fxcpmadd 24, 5, 13, 18");
  asm volatile ( "fxcpmadd 25, 3, 14, 19");
  asm volatile ( "fxcpmadd 26, 4, 14, 20");
  asm volatile ( "fxcpmadd 27, 5, 14, 21");
  QuadLoad(U_P(0,2,0,3), 6);
  QuadLoad(U_P(0,2,1,3), 7);
  QuadLoad(U_P(0,2,2,3), 8);

  asm volatile ( "fxcxnsma 16, 3, 13, 22");
  asm volatile ( "fxcxnsma 17, 4, 13, 23");
  asm volatile ( "fxcxnsma 18, 5, 13, 24");
  asm volatile ( "fxcxnsma 19, 3, 14, 25");
  asm volatile ( "fxcxnsma 20, 4, 14, 26");
  asm volatile ( "fxcxnsma 21, 5, 14, 27");


  asm volatile ( "fpmadd   13, 29, 11, 9");
  asm volatile ( "fpmadd   14, 29, 12, 10");

  asm volatile ( "fxcpmadd 22, 6, 13, 16");
  asm volatile ( "fxcpmadd 23, 7, 13, 17");
  asm volatile ( "fxcpmadd 24, 8, 13, 18");
  asm volatile ( "fxcpmadd 25, 6, 14, 19");
  asm volatile ( "fxcpmadd 26, 7, 14, 20");
  asm volatile ( "fxcpmadd 27, 8, 14, 21");
 
  asm volatile ( "fxcxnsma 16, 6, 13, 22");
  asm volatile ( "fxcxnsma 17, 7, 13, 23");
  asm volatile ( "fxcxnsma 18, 8, 13, 24");
  asm volatile ( "fxcxnsma 19, 6, 14, 25");
  asm volatile ( "fxcxnsma 20, 7, 14, 26");
  asm volatile ( "fxcxnsma 21, 8, 14, 27");


  QuadStore(OUT3_P(0,0,0), 16);
  QuadStore(OUT3_P(0,1,0), 17);
  QuadStore(OUT3_P(0,2,0), 18);
  QuadStore(OUT3_P(0,0,1), 19);
  QuadStore(OUT3_P(0,1,1), 20);
  QuadStore(OUT3_P(0,2,1), 21);


}




/***********************************************************************/
/* Color multiply each of the 4 half spinors in_mu with u_mu. Use the  */
/* result to reconstruct part of the final full spinor output outf     */
/* The pointer in0p is used to prefetch the next in0                   */
/***********************************************************************/
void wilson_dslash_mat_trick(double *outf, 
			     double *u, 
			     double *wfm_tmp0, 
			     double *wfm_tmp1, 
			     double *wfm_tmp2, 
			     double *wfm_tmp3, 
			     double *in0, 
			     double *in1, 
			     double *in2, 
			     double *in3,
			     double *in0p)
{
  register int o = 16;

  double *u_pp;  // gauge field refetch pointer
  double *in_pp; // input half spinor prefetch pointer


  /*********************************************************************/
  /* u_mu * half_spinor_mu, mu=0                                       */  
  /*********************************************************************/
  QuadLoad (u, 0);
  QuadLoadU(u, 1, o);
  QuadLoadU(u, 2, o);

  QuadLoad (in0, 10);
  QuadLoadU(in0, 11, o);
  QuadLoadU(in0, 12, o);
  QuadLoadU(in0, 13, o);

  asm volatile ( "fxcpmadd 22, 0, 10, 30");
  asm volatile ( "fxcpmadd 23, 1, 10, 30");
  asm volatile ( "fxcpmadd 24, 2, 10, 30");
  asm volatile ( "fxcpmadd 25, 0, 13, 30");
  QuadLoadU(in0, 14, o);
  asm volatile ( "fxcpmadd 26, 1, 13, 30");
  asm volatile ( "fxcpmadd 27, 2, 13, 30");
  QuadLoadU(in0, 15, o);

  asm volatile ( "fxcxnpma 16, 0, 10, 22");
  asm volatile ( "fxcxnpma 17, 1, 10, 23");
  QuadLoadU(u, 3, o);
  asm volatile ( "fxcxnpma 18, 2, 10, 24");
  asm volatile ( "fxcxnpma 19, 0, 13, 25");
  QuadLoadU(u, 4, o);
  asm volatile ( "fxcxnpma 20, 1, 13, 26");
  asm volatile ( "fxcxnpma 21, 2, 13, 27");
  QuadLoadU(u, 5, o);
  
  asm volatile ( "fxcpmadd 22, 3, 11, 16");
  asm volatile ( "fxcpmadd 23, 4, 11, 17");
  QuadLoadU(u, 6, o);
  asm volatile ( "fxcpmadd 24, 5, 11, 18");
  asm volatile ( "fxcpmadd 25, 3, 14, 19");
  QuadLoadU(u, 7, o);
  asm volatile ( "fxcpmadd 26, 4, 14, 20");
  asm volatile ( "fxcpmadd 27, 5, 14, 21");
  QuadLoadU(u, 8, o);

  u_pp = u+2;
  asm volatile ( "fxcxnpma 16, 3, 11, 22");
  asm volatile ( "fxcxnpma 17, 4, 11, 23");
  DCBT(u_pp);
  u_pp = u+6;
  asm volatile ( "fxcxnpma 18, 5, 11, 24");
  DCBT(u_pp);
  u_pp = u+10;
  asm volatile ( "fxcxnpma 19, 3, 14, 25");
  DCBT(u_pp);
  u_pp = u+14;
  asm volatile ( "fxcxnpma 20, 4, 14, 26");
  DCBT(u_pp);
  u_pp = u+18;
  asm volatile ( "fxcxnpma 21, 5, 14, 27");
  DCBT(u_pp);

  asm volatile ( "fxcpmadd 22, 6, 12, 16");
  DCBT(in1);
  in_pp = in1+4;
  asm volatile ( "fxcpmadd 23, 7, 12, 17");
  DCBT(in_pp);
  in_pp = in1+8;
  asm volatile ( "fxcpmadd 24, 8, 12, 18");
  DCBT(in_pp);
  in_pp = in1+12;
  asm volatile ( "fxcpmadd 25, 6, 15, 19");
  DCBT(in_pp);
  in_pp = in1+16;
  asm volatile ( "fxcpmadd 26, 7, 15, 20");
  DCBT(in_pp);
  asm volatile ( "fxcpmadd 27, 8, 15, 21");
 
  asm volatile ( "fxcxnpma 16, 6, 12, 22");
  asm volatile ( "fxcxnpma 17, 7, 12, 23");
  asm volatile ( "fxcxnpma 18, 8, 12, 24");
  asm volatile ( "fxcxnpma 19, 6, 15, 25");
  asm volatile ( "fxcxnpma 20, 7, 15, 26");
  asm volatile ( "fxcxnpma 21, 8, 15, 27");

  QuadLoadU(u, 0, o);
  QuadLoadU(u, 1, o);
  QuadLoadU(u, 2, o);
  QuadStore(WFM_TMP0_P(0,0,0), 16);
  QuadLoad (in1, 10);
  QuadStore(WFM_TMP0_P(0,1,0), 17);
  QuadLoadU(in1, 11, o);
  QuadStore(WFM_TMP0_P(0,2,0), 18);
  QuadLoadU(in1, 12, o);
  QuadStore(WFM_TMP0_P(0,0,1), 19);
  QuadLoadU(in1, 13, o);
  QuadStore(WFM_TMP0_P(0,1,1), 20);
  QuadStore(WFM_TMP0_P(0,2,1), 21);


  /*********************************************************************/
  /* u_mu * half_spinor_mu, mu=1                                       */  
  /*********************************************************************/


  asm volatile ( "fxcpmadd 22, 0, 10, 30");
  asm volatile ( "fxcpmadd 23, 1, 10, 30");
  asm volatile ( "fxcpmadd 24, 2, 10, 30");
  asm volatile ( "fxcpmadd 25, 0, 13, 30");
  QuadLoadU(in1, 14, o);
  asm volatile ( "fxcpmadd 26, 1, 13, 30");
  asm volatile ( "fxcpmadd 27, 2, 13, 30");
  QuadLoadU(in1, 15, o);

  asm volatile ( "fxcxnpma 16, 0, 10, 22");
  asm volatile ( "fxcxnpma 17, 1, 10, 23");
  QuadLoadU(u, 3, o);
  asm volatile ( "fxcxnpma 18, 2, 10, 24");
  asm volatile ( "fxcxnpma 19, 0, 13, 25");
  QuadLoadU(u, 4, o);
  asm volatile ( "fxcxnpma 20, 1, 13, 26");
  asm volatile ( "fxcxnpma 21, 2, 13, 27");
  QuadLoadU(u, 5, o);
  
  asm volatile ( "fxcpmadd 22, 3, 11, 16");
  asm volatile ( "fxcpmadd 23, 4, 11, 17");
  QuadLoadU(u, 6, o);
  asm volatile ( "fxcpmadd 24, 5, 11, 18");
  asm volatile ( "fxcpmadd 25, 3, 14, 19");
  QuadLoadU(u, 7, o);
  asm volatile ( "fxcpmadd 26, 4, 14, 20");
  asm volatile ( "fxcpmadd 27, 5, 14, 21");
  QuadLoadU(u, 8, o);

  u_pp = u+2;
  asm volatile ( "fxcxnpma 16, 3, 11, 22");
  asm volatile ( "fxcxnpma 17, 4, 11, 23");
  DCBT(u_pp);
  u_pp = u+6;
  asm volatile ( "fxcxnpma 18, 5, 11, 24");
  DCBT(u_pp);
  u_pp = u+10;
  asm volatile ( "fxcxnpma 19, 3, 14, 25");
  DCBT(u_pp);
  u_pp = u+14;
  asm volatile ( "fxcxnpma 20, 4, 14, 26");
  DCBT(u_pp);
  u_pp = u+18;
  asm volatile ( "fxcxnpma 21, 5, 14, 27");
  DCBT(u_pp);

  asm volatile ( "fxcpmadd 22, 6, 12, 16");
  DCBT(in2);
  in_pp = in2+4;
  asm volatile ( "fxcpmadd 23, 7, 12, 17");
  DCBT(in_pp);
  in_pp = in2+8;
  asm volatile ( "fxcpmadd 24, 8, 12, 18");
  DCBT(in_pp);
  in_pp = in2+12;
  asm volatile ( "fxcpmadd 25, 6, 15, 19");
  DCBT(in_pp);
  in_pp = in2+16;
  asm volatile ( "fxcpmadd 26, 7, 15, 20");
  DCBT(in_pp);
  asm volatile ( "fxcpmadd 27, 8, 15, 21");
 
  asm volatile ( "fxcxnpma 16, 6, 12, 22");
  asm volatile ( "fxcxnpma 17, 7, 12, 23");
  asm volatile ( "fxcxnpma 18, 8, 12, 24");
  asm volatile ( "fxcxnpma 19, 6, 15, 25");
  asm volatile ( "fxcxnpma 20, 7, 15, 26");
  asm volatile ( "fxcxnpma 21, 8, 15, 27");

  QuadLoadU(u, 0, o);
  QuadLoadU(u, 1, o);
  QuadLoadU(u, 2, o);
  QuadStore(WFM_TMP1_P(0,0,0), 16);
  QuadLoad (in2, 10);
  QuadStore(WFM_TMP1_P(0,1,0), 17);
  QuadLoadU(in2, 11, o);
  QuadStore(WFM_TMP1_P(0,2,0), 18);
  QuadLoadU(in2, 12, o);
  QuadStore(WFM_TMP1_P(0,0,1), 19);
  QuadLoadU(in2, 13, o);
  QuadStore(WFM_TMP1_P(0,1,1), 20);
  QuadStore(WFM_TMP1_P(0,2,1), 21);


  /*********************************************************************/
  /* u_mu * half_spinor_mu, mu=2                                       */  
  /*********************************************************************/


  asm volatile ( "fxcpmadd 22, 0, 10, 30");
  asm volatile ( "fxcpmadd 23, 1, 10, 30");
  asm volatile ( "fxcpmadd 24, 2, 10, 30");
  asm volatile ( "fxcpmadd 25, 0, 13, 30");
  QuadLoadU(in2, 14, o);
  asm volatile ( "fxcpmadd 26, 1, 13, 30");
  asm volatile ( "fxcpmadd 27, 2, 13, 30");
  QuadLoadU(in2, 15, o);

  asm volatile ( "fxcxnpma 16, 0, 10, 22");
  asm volatile ( "fxcxnpma 17, 1, 10, 23");
  QuadLoadU(u, 3, o);
  asm volatile ( "fxcxnpma 18, 2, 10, 24");
  asm volatile ( "fxcxnpma 19, 0, 13, 25");
  QuadLoadU(u, 4, o);
  asm volatile ( "fxcxnpma 20, 1, 13, 26");
  asm volatile ( "fxcxnpma 21, 2, 13, 27");
  QuadLoadU(u, 5, o);
  
  asm volatile ( "fxcpmadd 22, 3, 11, 16");
  asm volatile ( "fxcpmadd 23, 4, 11, 17");
  QuadLoadU(u, 6, o);
  asm volatile ( "fxcpmadd 24, 5, 11, 18");
  asm volatile ( "fxcpmadd 25, 3, 14, 19");
  QuadLoadU(u, 7, o);
  asm volatile ( "fxcpmadd 26, 4, 14, 20");
  asm volatile ( "fxcpmadd 27, 5, 14, 21");
  QuadLoadU(u, 8, o);

  u_pp = u+2;
  asm volatile ( "fxcxnpma 16, 3, 11, 22");
  asm volatile ( "fxcxnpma 17, 4, 11, 23");
  DCBT(u_pp);
  u_pp = u+6;
  asm volatile ( "fxcxnpma 18, 5, 11, 24");
  DCBT(u_pp);
  u_pp = u+10;
  asm volatile ( "fxcxnpma 19, 3, 14, 25");
  DCBT(u_pp);
  u_pp = u+14;
  asm volatile ( "fxcxnpma 20, 4, 14, 26");
  DCBT(u_pp);
  u_pp = u+18;
  asm volatile ( "fxcxnpma 21, 5, 14, 27");
  DCBT(u_pp);

  asm volatile ( "fxcpmadd 22, 6, 12, 16");
  DCBT(in3);
  in_pp = in3+4;
  asm volatile ( "fxcpmadd 23, 7, 12, 17");
  DCBT(in_pp);
  in_pp = in3+8;
  asm volatile ( "fxcpmadd 24, 8, 12, 18");
  DCBT(in_pp);
  in_pp = in3+12;
  asm volatile ( "fxcpmadd 25, 6, 15, 19");
  DCBT(in_pp);
  in_pp = in3+16;
  asm volatile ( "fxcpmadd 26, 7, 15, 20");
  DCBT(in_pp);
  asm volatile ( "fxcpmadd 27, 8, 15, 21");
 
  asm volatile ( "fxcxnpma 16, 6, 12, 22");
  asm volatile ( "fxcxnpma 17, 7, 12, 23");
  asm volatile ( "fxcxnpma 18, 8, 12, 24");
  asm volatile ( "fxcxnpma 19, 6, 15, 25");
  asm volatile ( "fxcxnpma 20, 7, 15, 26");
  asm volatile ( "fxcxnpma 21, 8, 15, 27");

  QuadLoadU(u, 0, o);
  QuadLoadU(u, 1, o);
  QuadLoadU(u, 2, o);

  QuadStore(WFM_TMP2_P(0,0,0), 16);
  QuadLoad (in3, 10);
  QuadStore(WFM_TMP2_P(0,1,0), 17);
  QuadLoadU(in3, 11, o);
  QuadStore(WFM_TMP2_P(0,2,0), 18);
  QuadLoadU(in3, 12, o);
  QuadStore(WFM_TMP2_P(0,0,1), 19);
  QuadLoadU(in3, 13, o);
  QuadStore(WFM_TMP2_P(0,1,1), 20);
  QuadStore(WFM_TMP2_P(0,2,1), 21);


  /*********************************************************************/
  /* u_mu * half_spinor_mu, mu=3                                       */  
  /*********************************************************************/


  asm volatile ( "fxcpmadd 22, 0, 10, 30");
  asm volatile ( "fxcpmadd 23, 1, 10, 30");
  asm volatile ( "fxcpmadd 24, 2, 10, 30");
  asm volatile ( "fxcpmadd 25, 0, 13, 30");
  QuadLoadU(in3, 14, o);
  asm volatile ( "fxcpmadd 26, 1, 13, 30");
  asm volatile ( "fxcpmadd 27, 2, 13, 30");
  QuadLoadU(in3, 15, o);

  asm volatile ( "fxcxnpma 16, 0, 10, 22");
  asm volatile ( "fxcxnpma 17, 1, 10, 23");
  QuadLoadU(u, 3, o);
  asm volatile ( "fxcxnpma 18, 2, 10, 24");
  asm volatile ( "fxcxnpma 19, 0, 13, 25");
  QuadLoadU(u, 4, o);
  asm volatile ( "fxcxnpma 20, 1, 13, 26");
  asm volatile ( "fxcxnpma 21, 2, 13, 27");
  QuadLoadU(u, 5, o);
  
  asm volatile ( "fxcpmadd 22, 3, 11, 16");
  asm volatile ( "fxcpmadd 23, 4, 11, 17");
  QuadLoadU(u, 6, o);
  asm volatile ( "fxcpmadd 24, 5, 11, 18");
  asm volatile ( "fxcpmadd 25, 3, 14, 19");
  QuadLoadU(u, 7, o);
  asm volatile ( "fxcpmadd 26, 4, 14, 20");
  asm volatile ( "fxcpmadd 27, 5, 14, 21");
  QuadLoadU(u, 8, o);

  //  u_pp = u+18;
  asm volatile ( "fxcxnpma 16, 3, 11, 22");
  asm volatile ( "fxcxnpma 17, 4, 11, 23");
  //  DCBT(u_pp);
  //  u_pp = u+22;
  asm volatile ( "fxcxnpma 18, 5, 11, 24");
  //  DCBT(u_pp);
  //  u_pp = u+26;
  asm volatile ( "fxcxnpma 19, 3, 14, 25");
  //  DCBT(u_pp);
  //  u_pp = u+30;
  asm volatile ( "fxcxnpma 20, 4, 14, 26");
  //  DCBT(u_pp);
  //  u_pp = u+34;
  asm volatile ( "fxcxnpma 21, 5, 14, 27");
  //  DCBT(u_pp);

  asm volatile ( "fxcpmadd 22, 6, 12, 16");
  DCBT(in0p);
  in_pp = in0p+4;
  asm volatile ( "fxcpmadd 23, 7, 12, 17");
  DCBT(in_pp);
  in_pp = in0p+8;
  asm volatile ( "fxcpmadd 24, 8, 12, 18");
  DCBT(in_pp);
  in_pp = in0p+12;
  asm volatile ( "fxcpmadd 25, 6, 15, 19");
  DCBT(in_pp);
  in_pp = in0p+16;
  asm volatile ( "fxcpmadd 26, 7, 15, 20");
  DCBT(in_pp);
  asm volatile ( "fxcpmadd 27, 8, 15, 21");
 
  asm volatile ( "fxcxnpma 16, 6, 12, 22");
  asm volatile ( "fxcxnpma 17, 7, 12, 23");
  asm volatile ( "fxcxnpma 18, 8, 12, 24");
  asm volatile ( "fxcxnpma 19, 6, 15, 25");
  asm volatile ( "fxcxnpma 20, 7, 15, 26");
  asm volatile ( "fxcxnpma 21, 8, 15, 27");

  QuadStore(WFM_TMP3_P(0,0,0), 16);
  QuadStore(WFM_TMP3_P(0,1,0), 17);
  QuadStore(WFM_TMP3_P(0,2,0), 18);
  QuadStore(WFM_TMP3_P(0,0,1), 19);
  QuadStore(WFM_TMP3_P(0,1,1), 20);
  QuadStore(WFM_TMP3_P(0,2,1), 21);


  /*********************************************************************/
  /* trick the half spinoros in the wfm_tmp0,1,2,3 arrays into the     */
  /* full outf spinor                                                  */
  /*********************************************************************/
  /* c=0 */
  asm volatile ( "fpmr 0, 30");
  asm volatile ( "fpmr 1, 30");
  asm volatile ( "fpmr 2, 30");
  asm volatile ( "fpmr 3, 30");

  QuadLoad(WFM_TMP0_P(0,0,0), 10);
  QuadLoad(WFM_TMP0_P(0,0,1), 11);    
  QuadLoad(WFM_TMP1_P(0,0,0), 12);
  QuadLoad(WFM_TMP1_P(0,0,1), 13);
  QuadLoad(WFM_TMP2_P(0,0,0), 14);
  asm volatile ( "fpmadd   20,31,10,0");
  asm volatile ( "fpmadd   21,31,11,1");
  asm volatile ( "fxcxnsma 22,29,11,2");
  QuadLoad(WFM_TMP2_P(0,0,1), 15);
  asm volatile ( "fxcxnsma 23,29,10,3");
  asm volatile ( "fpmadd   0,31,12,20");
  asm volatile ( "fpmadd   1,31,13,21");
  QuadLoad(WFM_TMP3_P(0,0,0), 16);
  asm volatile ( "fpmadd   2,29,13,22");
  asm volatile ( "fpnmsub  3,29,12,23");
  asm volatile ( "fpmadd   20,31,14,0");
  QuadLoad(WFM_TMP3_P(0,0,1), 17);
  asm volatile ( "fpmadd   21,31,15,1");
  DCBT(wfm_tmp0+12);
  asm volatile ( "fxcxnsma 22,29,14,2");
  asm volatile ( "fxcxnpma 23,29,15,3");
  DCBT(wfm_tmp0+16);
  asm volatile ( "fpmadd   0,31,16,20");
  asm volatile ( "fpmadd   1,31,17,21");
  DCBT(wfm_tmp0+20);
  asm volatile ( "fpmadd   2,29,16,22");
  asm volatile ( "fpmadd   3,29,17,23");
  QuadStore(OUTF_P(0,0,0), 0);
  QuadStore(OUTF_P(0,0,1), 1);
  QuadStore(OUTF_P(0,0,2), 2);
  QuadStore(OUTF_P(0,0,3), 3);
 
  
  /* c=1 */
  asm volatile ( "fpmr 0, 30");
  asm volatile ( "fpmr 1, 30");
  asm volatile ( "fpmr 2, 30");
  asm volatile ( "fpmr 3, 30");

  QuadLoad(WFM_TMP0_P(0,1,0), 10);
  QuadLoad(WFM_TMP0_P(0,1,1), 11);    
  QuadLoad(WFM_TMP1_P(0,1,0), 12);
  QuadLoad(WFM_TMP1_P(0,1,1), 13);
  QuadLoad(WFM_TMP2_P(0,1,0), 14);
  asm volatile ( "fpmadd   20,31,10,0");
  asm volatile ( "fpmadd   21,31,11,1");
  asm volatile ( "fxcxnsma 22,29,11,2");
  QuadLoad(WFM_TMP2_P(0,1,1), 15);
  asm volatile ( "fxcxnsma 23,29,10,3");
  asm volatile ( "fpmadd   0,31,12,20");
  asm volatile ( "fpmadd   1,31,13,21");
  QuadLoad(WFM_TMP3_P(0,1,0), 16);
  asm volatile ( "fpmadd   2,29,13,22");
  asm volatile ( "fpnmsub  3,29,12,23");
  asm volatile ( "fpmadd   20,31,14,0");
  QuadLoad(WFM_TMP3_P(0,1,1), 17);
  asm volatile ( "fpmadd   21,31,15,1");
  asm volatile ( "fxcxnsma 22,29,14,2");
  asm volatile ( "fxcxnpma 23,29,15,3");
  asm volatile ( "fpmadd   0,31,16,20");
  asm volatile ( "fpmadd   1,31,17,21");
  QuadStore(OUTF_P(0,1,0), 0);
  asm volatile ( "fpmadd   2,29,16,22");
  QuadStore(OUTF_P(0,1,1), 1);
  asm volatile ( "fpmadd   3,29,17,23");
  QuadStore(OUTF_P(0,1,2), 2);
  QuadStore(OUTF_P(0,1,3), 3);
   
  /* c=2 */
  asm volatile ( "fpmr 0, 30");
  asm volatile ( "fpmr 1, 30");
  asm volatile ( "fpmr 2, 30");
  asm volatile ( "fpmr 3, 30");

  QuadLoad(WFM_TMP0_P(0,2,0), 10);
  QuadLoad(WFM_TMP0_P(0,2,1), 11);    
  QuadLoad(WFM_TMP1_P(0,2,0), 12);
  QuadLoad(WFM_TMP1_P(0,2,1), 13);
  QuadLoad(WFM_TMP2_P(0,2,0), 14);
  asm volatile ( "fpmadd   20,31,10,0");
  asm volatile ( "fpmadd   21,31,11,1");
  asm volatile ( "fxcxnsma 22,29,11,2");
  QuadLoad(WFM_TMP2_P(0,2,1), 15);
  asm volatile ( "fxcxnsma 23,29,10,3");
  asm volatile ( "fpmadd   0,31,12,20");
  asm volatile ( "fpmadd   1,31,13,21");
  QuadLoad(WFM_TMP3_P(0,2,0), 16);
  asm volatile ( "fpmadd   2,29,13,22");
  asm volatile ( "fpnmsub  3,29,12,23");
  asm volatile ( "fpmadd   20,31,14,0");
  QuadLoad(WFM_TMP3_P(0,2,1), 17);
  asm volatile ( "fpmadd   21,31,15,1");
  asm volatile ( "fxcxnsma 22,29,14,2");
  asm volatile ( "fxcxnpma 23,29,15,3");
  asm volatile ( "fpmadd   0,31,16,20");
  asm volatile ( "fpmadd   1,31,17,21");
  asm volatile ( "fpmadd   2,29,16,22");
  asm volatile ( "fpmadd   3,29,17,23");
  QuadStore(OUTF_P(0,2,0), 0);
  QuadStore(OUTF_P(0,2,1), 1);
  QuadStore(OUTF_P(0,2,2), 2);
  QuadStore(OUTF_P(0,2,3), 3);
 
}



CPS_END_NAMESPACE
