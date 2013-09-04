#ifdef USE_SSE
#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-01-08 21:09:25 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/sse/sse-dwf_dslash_4.C,v 1.4 2013-01-08 21:09:25 chulwoo Exp $
//  $Id: sse-dwf_dslash_4.C,v 1.4 2013-01-08 21:09:25 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: sse-dwf_dslash_4.C,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/sse/sse-dwf_dslash_4.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// dwf_dslash_4.C
//
// dwf_dslash_4 is the derivative part of the space-time part of
// the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
// cb = 0/1 <--> even/odd checkerboard of in field.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
//#include <omp.h>
#include <util/omp_wrapper.h>
#include <pmmintrin.h>
#include <util/wilson.h>

#include<config.h>
#include<util/dwf.h>
#include<util/wilson.h>
#include<util/dirac_op.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/time_cps.h>

#include <comms/sysfunc_qmp.h>

// TIZB, restore later !
#pragma warning disable 592

#define __RESTRICT 
#define __RESTRICT


// #define AXG_WILSON


#define X_DIR_ON
#define Y_DIR_ON
#define Z_DIR_ON
#define T_DIR_ON

#define POS_DIR_ON
#define NEG_DIR_ON

//#define DEBUG_PRINT_DAG0
//#define DEBUG_PRINT_DAG1
//#define DEBUG_NOEDGE
#define DEBUG_MFLOPS (1*1);

// use macro version or inlines
#define USE_MACROS

//boundary communication on or not
#define BND_COMM

//#define PROFILE
#undef PROFILE

//#define PROF_printf(...)  printf(__VA_ARGS__);
//#define PROF_printf(...)  if (!UniqueID()) printf(__VA_ARGS__);
#define PROF_printf(...)  {}


// The sizes of the block are assumed to be even numbers
#define Y_BLOCK 16
#define Z_BLOCK 16


// USE register  t00-t11 to add into
//#define ADD2REG

// use HERN2 (doesn't use the temporal registers _a,_b,_c,_d)
//#define USE_HERN2


#ifndef AXG_WILSON
#define GAUGE_SIZE 72 /* re/im * colors*colors * dim */
#else
#define GAUGE_SIZE 54 /* re/im * colors*colors * dim */
#endif

CPS_START_NAMESPACE

void wilson_dslash_bnd_dag0(
                   IFloat *chi_p_f, 
		   IFloat *u_p_f, 
		   IFloat *psi_p_f, 
		   int cb,
		   Wilson *wilson_p);

void wilson_dslash_bnd_dag1(
                   IFloat *chi_p_f, 
		   IFloat *u_p_f, 
		   IFloat *psi_p_f, 
		   int cb,
		   Wilson *wilson_p);
  
void wilson_dslash_blk_dag0(
                   IFloat *chi_p_f, 
		   IFloat *u_p_f, 
		   IFloat *psi_p_f, 
		   int cb,
		   Wilson *wilson_p);

void wilson_dslash_blk_dag1(
                   IFloat *chi_p_f, 
		   IFloat *u_p_f, 
		   IFloat *psi_p_f, 
		   int cb,
		   Wilson *wilson_p);




void dwf_dslash_4(Vector *out, 
		  Matrix *gauge_field, 
		  Vector *in, 
		  int cb_5, 
		  int dag, 
		  Dwf *dwf_lib_arg)
{
  int i;
  int ls;
  IFloat *frm_in;
  IFloat *frm_out;
  IFloat *g_field;
  Wilson *wilson_p;
  //  int size_cb[2];
  int size_cb;
  int parity;

  
  //----------------------------------------------------------------
  // Initializations
  //----------------------------------------------------------------
  ls = dwf_lib_arg->ls;
  frm_in = (IFloat *) in;
  frm_out = (IFloat *) out;
  wilson_p = dwf_lib_arg->wilson_p;
  //  size_cb[0] = 24*wilson_p->vol[0]; // same
  //  size_cb[1] = 24*wilson_p->vol[1];  // same
   size_cb = 24*wilson_p->vol[1];  // same

   IFloat *u_p_f = (IFloat *) gauge_field;

#ifdef PROFILE
  //int NITR = wilson_p->NITR;
  const Float MultFlops=ls*wilson_p->MultFlops;
  const Float MultFlops_bnd=ls*wilson_p->MultFlops_bnd;
  const Float MultFlops_blk=ls*wilson_p->MultFlops_blk;
  const Float CPU_GHZ = wilson_p->CPU_GHZ;
  const int  num_threads= wilson_p->num_threads;
  Float time_bnd(0), time_blk(0), time_comm1(0), time_comm2(0);
#endif

   
  //----------------------------------------------------------------
  // Apply 4-dimensional Dslash
  //----------------------------------------------------------------
   for(int parity_5=0; parity_5<2;++parity_5){
     int cb = parity_5^cb_5;
     for(i=parity_5; i<ls; i+=2){
       // parity of 4-D checkerboard


       //------------------------------------------------------------
       //parity = (i + cb_5) % 2;
       //parity = (i&1)^cb_5;
       
       // Apply on 4-dim "parity" checkerboard part
       //------------------------------------------------------------
       //wilson_dslash(frm_out, g_field, frm_in, parity, dag, wilson_p);
       {
	 
	 IFloat *chi_p_f = frm_out + i*size_cb; 
	 IFloat *psi_p_f = frm_in  + i*size_cb;
	 
	 
	 const int   lx = wilson_p->ptr[0];
	 const int   ly = wilson_p->ptr[1];
	 const int   lz = wilson_p->ptr[2];
	 const int   lt = wilson_p->ptr[3];
	 const int   vol = wilson_p->vol[0];
	 
	 

	 IFloat **send_buf = wilson_p->send_buf;
	 
	 IFloat **recv_buf =wilson_p->recv_buf;
	 
	 
#ifdef PROFILE
	   Float time = -dclock();
#endif  
	   
	   /* Do the boundary parts */
	   
	   if(dag == 0 ) wilson_dslash_bnd_dag0(chi_p_f, u_p_f, psi_p_f, cb, wilson_p);
	   else          wilson_dslash_bnd_dag1(chi_p_f, u_p_f, psi_p_f, cb, wilson_p);

#ifdef PROFILE
	   time_bnd+= time+dclock();
	   time = -dclock();
#endif

#ifdef BND_COMM
	   
	   /*
	       Wait for finishing communication
	   */
	   for(int idir=0;idir<4;++idir)
	     if(GJP.Nodes(idir)!=1){
	       //printf("wait for communication for %d\n",i);
	       //fflush(stdin);
	       
	       QMP_status_t  status = QMP_wait(wilson_p->multiple[idir]);
	       if (status != QMP_SUCCESS)
		 QMP_error("Error in sse-dwf_dslash_4:%s\n",
			   QMP_error_string(status));
	       
	       status = QMP_wait(wilson_p->multiple[idir+4]);
	       if (status != QMP_SUCCESS)
		 QMP_error("Error in  sse-dwf_dslash_4:%s\n",
			   QMP_error_string(status));
	     }
#else
	   if(dag==1){//only do the communication for dag=1
	     	   /*
		     Wait for finishing communication
		   */
	     for(int idir=0;idir<4;++idir)
	       if(GJP.Nodes(idir)!=1){
		 //printf("wait for communication for %d\n",i);
		 //fflush(stdin);
		 
		 QMP_status_t  status = QMP_wait(wilson_p->multiple[idir]);
		 if (status != QMP_SUCCESS)
		   QMP_error("Error in  sse-dwf_dslash_4:%s\n",
			     QMP_error_string(status));
		 
		 status = QMP_wait(wilson_p->multiple[idir+4]);
		 if (status != QMP_SUCCESS)
		   QMP_error("Error in  sse-dwf_dslash_4:%s\n",
			     QMP_error_string(status));
	       }
	     
	   } else {
	     
	     int block[4];
	     block[0]=HALF_SPINOR_SIZE*ly*lz*lt/2;
	     block[1]=HALF_SPINOR_SIZE*lx*lz*lt/2;
	     block[2]=HALF_SPINOR_SIZE*lx*ly*lt/2;
	     block[3]=HALF_SPINOR_SIZE*lx*ly*lz/2;
	     
	     for(int idir=0;idir<4;++idir)
	       if(GJP.Nodes(idir)!=1)
		 for(int i=0;i<block[idir];++i) {
		   *(i+(IFloat*)(wilson_p->recv_buf[idir]))=0;
		   *(i+(IFloat*)(wilson_p->recv_buf[4+idir]))=0;
		 }
	   }
#endif

	   
	   
#ifdef PROFILE
	   time_comm2 += time+dclock();
	   time = -dclock();
#endif
	   
	   /* Do the bulk part */
	   
	   if(dag == 0) {
	     wilson_dslash_blk_dag0(chi_p_f, u_p_f, psi_p_f, cb, wilson_p);
	   }
	   else{
	     wilson_dslash_blk_dag1(chi_p_f, u_p_f, psi_p_f, cb, wilson_p);
	   }
	   
#ifdef PROFILE
	   time_blk+= time+dclock();
#endif
	 
       }
     }//loop(ls)
   }//lop(parity_5)

DiracOp::CGflops += 1320*wilson_p->vol[0]*ls;

#ifdef PROFILE
	 Float time_tot= time_bnd+time_blk+time_comm2;
	 const Float MultFlops_percent = MultFlops/time_tot
	   *1e-9/(CPU_GHZ * 4* num_threads)*100;
	 const Float MultFlops_blk_percent = MultFlops_blk/time_blk
	   *1e-9/(CPU_GHZ * 4* num_threads)*100;
	 const Float MultFlops_bnd_percent = MultFlops_bnd/time_bnd
	   *1e-9/(CPU_GHZ * 4* num_threads)*100;
	 
#if 0
	 PROF_printf("(cb5 %d dag %d) %e flops/ %e seconds = %e GFlops [ %4.2f %% peak] (%d threads)\n", 
		     cb_5,dag,(double)MultFlops,time_tot, MultFlops/time_tot*1e-9,
		     MultFlops_percent, num_threads);
#endif
	 
	 PROF_printf("  blk: %4.3e sec %4.3e GFlops (%2.1f %% tm, %2.1f %%cmp) [%4.2f %% peak]\n",
		     time_blk, MultFlops_blk/time_blk*1e-9, 
		     100*time_blk/time_tot, 100*(MultFlops_blk/MultFlops),
		     MultFlops_blk_percent);
	 PROF_printf("  bnd: %4.3e sec %4.3e GFlops (%2.1f %% tm, %2.1f %%cmp) [%4.2f %% peak]\n",
		     time_bnd, MultFlops_bnd/time_bnd*1e-9, 
		     100*time_bnd/time_tot, 100*(MultFlops_bnd/MultFlops),
		     MultFlops_bnd_percent);
	 PROF_printf("  comm2: %4.3e sec  (%2.1f %% tm)\n",
		     time_comm2, 100*time_comm2/time_tot);
	 
#endif

}



#include "sse-subs.h"

#include "sse-blk-dag0.h"
#include "sse-blk-dag1.h"


#define U(r,row,col,d)  *(u+(r+2*(row+3*(col+3*d))))
#define PSI(r,c,s)      *(psi +(r+2*(c+3*s)))
//! As above, but the vector is called chi
#define CHI(r,c,s)      *(chi +(r+2*(c+3*s)))
#define TMP(r,c,s)      *(tmp +(r+2*(c+3*s))) 
#define TMP1(r,c,s)     *(tmp1+(r+2*(c+3*s))) 
#define TMP2(r,c,s)     *(tmp2+(r+2*(c+3*s))) 
#define TMP3(r,c,s)     *(tmp3+(r+2*(c+3*s))) 
#define TMP4(r,c,s)     *(tmp4+(r+2*(c+3*s))) 
#define TMP5(r,c,s)     *(tmp5+(r+2*(c+3*s))) 
#define TMP6(r,c,s)     *(tmp6+(r+2*(c+3*s))) 
#define TMP7(r,c,s)     *(tmp7+(r+2*(c+3*s))) 
#define TMP8(r,c,s)     *(tmp8+(r+2*(c+3*s))) 
#define FBUF(r,c,s)     *(fbuf+(r+2*(c+3*s))) 

#include "sse-bnd-dag0.h"
#include "sse-bnd-dag1.h"


CPS_END_NAMESPACE
#endif
