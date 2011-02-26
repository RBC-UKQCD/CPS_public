#ifdef USE_SSE
#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/sse/dwf_dslash.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// dwf_dslash.C
//
// dwf_dslash is the derivative part of the fermion matrix. 
// The in, out fields are defined on the full lattice
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/dwf.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/wilson.h>
CPS_START_NAMESPACE

Float dclock(void);
//#define DEBUG_DWF_DSLASH(msg,a ...) do \
//     printf("[%05d] %s:%d:QMP/%s(): " msg, UniqueID(), \
//                             __FILE__,__LINE__,__FUNCTION__,##a); \
//  while(0);

#ifdef  DEBUG_DWF_DSLASH
#undef DEBUG_DWF_DSLASH
#define DEBUG_DWF_DSLASH(msg,a ...) do \
     printf("[%05d] " msg, UniqueID() ,##a); \
  while(0);
#else
#define DEBUG_DWF_DSLASH(msg,a ...) {}
#endif


void dwf_dslash(Vector *out, 
		Matrix *gauge_field, 
		Vector *in, 
		Float mass,
		int cb, 
		int dag, 
		Dwf *dwf_lib_arg)
{
//------------------------------------------------------------------
// Apply 4-dimensional Dslash
//------------------------------------------------------------------
//	if (!UniqueID()) DEBUG_DWF_DSLASH("dwf_dslash");
#define PROFILE

#ifdef PROFILE
  Wilson* wilson_p = dwf_lib_arg->wilson_p;
  const int ls = dwf_lib_arg->ls;
  const Float MultFlops=wilson_p->MultFlops * ls;
  //const Float MultFlops_bnd=wilson_p->MultFlops_bnd;
  //const Float MultFlops_blk=wilson_p->MultFlops_blk;
  const Float CPU_GHZ = wilson_p->CPU_GHZ;
  const int  num_threads= wilson_p->num_threads;

  Float time_4D=-dclock();
#endif
  
  dwf_dslash_4(out, gauge_field, in, cb, dag, dwf_lib_arg);



#ifdef PROFILE
  time_4D += dclock();
  Float MultFlops_percent = MultFlops/time_4D
    *1e-9/(CPU_GHZ * 4* num_threads)*100;

  if(!UniqueID())
  DEBUG_DWF_DSLASH("DWF4D (cb %d dag %d) %e flops/ %e seconds = %e GFlops [ %4.2f %% peak] (%d threads)\n", 
	 cb,dag,(double)MultFlops,time_4D, MultFlops/time_4D*1e-9,
	 MultFlops_percent, num_threads);

#endif

//------------------------------------------------------------------
// Apply 5th-direction Dslash
//------------------------------------------------------------------
#ifdef PROFILE
  const int ITER5= 1;//100;
  const Float MultFlops5= 12*ls* wilson_p->vol[0] * 2*2 * ITER5; //mul+add, re+im
  //
  // N.B.   we can kill mul by change of the varibale:
  //        psi'(s)  =  2^s  P_L psi  + 2^{-s}  P_R psi
  //  or something like that...
  //
  Float time_5D=-dclock();
#endif
  
  for(int i=0;i<ITER5;++i)
  dwf_dslash_5_plus(out, in, mass, dag, dwf_lib_arg);

#ifdef PROFILE
  time_5D += dclock();
  MultFlops_percent = MultFlops5/time_5D
    *1e-9/(CPU_GHZ * 4* num_threads)*100;

  if(!UniqueID())
  DEBUG_DWF_DSLASH("DWF5D (cb %d dag %d) %e flops/ %e seconds = %e GFlops [ %4.2f %% peak] (%d threads)\n", 
	 cb,dag,(double)MultFlops5,time_5D, MultFlops5/time_5D*1e-9,
	 MultFlops_percent, num_threads);

#endif



}

CPS_END_NAMESPACE
#endif
