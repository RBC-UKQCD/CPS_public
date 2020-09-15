#ifdef USE_SSE
#include <pmmintrin.h>
#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/noarch/mobius_dslash_5_plus-nonowait.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// mobius_dslash_5_plus.C
//
// mobius_dslash_5_plus is the derivative part of the 5th direction
// part of the fermion matrix. This routine accumulates the result
// on the out field 
// The in, out fields are defined on the checkerboard lattice.
// The action of this operator is the same for even/odd
// checkerboard fields because there is no gauge field along
// the 5th direction.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//
//
// Storage order for DWF fermions
//------------------------------------------------------------------
//  
//  |     |
//  | |r| |
//  | |i| |
//  |     |
//  | |r| | = |spin comp|
//  | |i| |
//  |     |
//  | |r| |
//  | |i| |
//  |     |
//  
//  
//  |             |
//  | |spin comp| |
//  |             |
//  |             |
//  | |spin comp| | 
//  |             | = |spinor|
//  |             |
//  | |spin comp| |
//  |             |
//  |             |
//  | |spin comp| |
//  |             |
//  
//  
//  |            |
//  |  |spinor|  |
//  |  |spinor|  |
//  |  |spinor|  |
//  |  |spinor|  |
//  |  |spinor|  |
//  |     .      | = |s-block|   The spinors are arranged in Wilson
//  |     .      |               order with odd - even 4d-checkerboard
//  |     .      |               storage.
//  |evn/odd vol |
//  |     .      |
//  |  |spinor|  |
//  |            |
//  
//  
//  |                |
//  | |s-block even| |  For even chckerboard
//  | |s-block odd|  |
//  | |s-block even| |
//  | |s-block odd|  |
//  |       .        |
//  |       .        |
//  |       .        |
//  |                |
//
//
//  |                |
//  | |s-block odd|  |  For odd chckerboard
//  | |s-block even| |
//  | |s-block odd|  |
//  | |s-block even| |
//  |       .        |
//  |       .        |
//  |       .        |
//  |                |
//
//------------------------------------------------------------------


CPS_END_NAMESPACE
#include<util/mobius.h>
#include<util/gjp.h>
#include<util/dirac_op.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/smalloc.h>



#ifdef USE_BLAS
#include<util/qblas_extend.h>
#endif
#include<comms/scu.h>
CPS_START_NAMESPACE

#include "../../d_op_dwf/sse/sse-subs-axpys.h"

#undef USE_BLAS

inline
void mobius_kappa_dslash_5_plus_dag0(Vector *out, 
			    Vector *in, 
			    Float mass,
			    Dwf *mobius_lib_arg,
			    Float a_five_inv
			    )
{
  //const IFloat two_over_a5 = 2.0 * a_five_inv;
  //const IFloat neg_mass_two_over_a5 = -2.0 * mass * a_five_inv;
  // for mobius, use 1/2*(1+-g5) for projectors
  // also, a_five_inv is a non-trivial factor depending on mobius kappa's or b,c coeff's,
  // depending on when it is called. Don't change name to avoid enbug.
  const IFloat two_over_a5 = 1.0 * a_five_inv;
  const IFloat neg_mass_two_over_a5 = -1.0 * mass * a_five_inv;
  const int local_ls    = GJP.SnodeSites(); 
  const int s_nodes     = GJP.Snodes();
  const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb   = mobius_lib_arg->vol_4d / 2;
  const int ls_stride   = 24 * vol_4d_cb;
  const int max_dex((local_ls-1)*vol_4d_cb);

  if (s_nodes != 1 || s_node_coor !=0) {
    // error!
  }
  
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
    int x;
    int idx;
    IFloat *f_in (reinterpret_cast<IFloat*>(in));
    IFloat *f_out(reinterpret_cast<IFloat*>(out));

    int dag_shift_pls;
    int dag_shift_min;
    int first12_target;
    int second12_target;

    register __m128d al;
    al = _mm_loaddup_pd(&two_over_a5);
    register __m128d nal;
    nal = _mm_loaddup_pd(&neg_mass_two_over_a5);

    DECLARE_AXPY_TMPVARS;


#if 0
    if(dag == 1){
      dag_shift_pls   =  12;
      dag_shift_min   =  0;
      first12_target  = -ls_stride;
      second12_target =  ls_stride+12;
    }
    else{
#endif
      dag_shift_pls   =  0;
      dag_shift_min   =  12;
      first12_target  =  ls_stride;
      second12_target = -ls_stride+12;
#if 0
    }
#endif
    // index of the start of the last slice in the fifth dimension
    const int the_end((local_ls-1)*ls_stride);



#ifdef _OPENMP
#pragma omp for 
#endif
    for(x=0; x<vol_4d_cb; x++)
      {
	const int shift(24*x);
#ifdef USE_BLAS
	cblas_daxpy(12,two_over_a5,
		    f_in  + shift + dag_shift_pls,
		    f_out + shift + dag_shift_pls + ls_stride );
	cblas_daxpy(12,neg_mass_two_over_a5,
		    f_in  + shift + dag_shift_min,
		    f_out + shift + dag_shift_min + the_end );
#else
	DAXPY12( al, 	    f_in  + shift + dag_shift_pls,
		    f_out + shift + dag_shift_pls + ls_stride );

	DAXPY12( nal, 	
		    f_in  + shift + dag_shift_min,
		    f_out + shift + dag_shift_min + the_end );
#endif

      }
    

#ifdef _OPENMP
#pragma omp for 
#endif
    for (idx=vol_4d_cb;idx<max_dex;idx+=4)
      {
	DAXPY12( al, f_in +      24*idx+0 ,  f_out+ 24*idx+0    +  first12_target );
	DAXPY12( al, f_in + 12 + 24*idx+0 ,  f_out+ 24*idx+0    + second12_target );
	DAXPY12( al, f_in +      24*idx+24 ,  f_out+ 24*idx+24  +  first12_target );
	DAXPY12( al, f_in + 12 + 24*idx+24 ,  f_out+ 24*idx+24  + second12_target );
	DAXPY12( al, f_in +      24*idx+48 ,  f_out+ 24*idx+48  +  first12_target );
	DAXPY12( al, f_in + 12 + 24*idx+48 ,  f_out+ 24*idx+48  + second12_target );
	DAXPY12( al, f_in +      24*idx+72 ,  f_out+ 24*idx+72  +  first12_target );
	DAXPY12( al, f_in + 12 + 24*idx+72 ,  f_out+ 24*idx+72  + second12_target );
      }

    
#ifdef _OPENMP
#pragma omp for 
#endif
    for(x=0; x<vol_4d_cb; x++)
      {
	const int shift(24*x);
#ifdef USE_BLAS
	cblas_daxpy(12,neg_mass_two_over_a5,
		    f_in  + shift + dag_shift_pls + the_end,
		    f_out + shift + dag_shift_pls );

	cblas_daxpy(12,two_over_a5,
		    f_in  + shift + dag_shift_min + the_end,
		    f_out + shift + dag_shift_min + the_end - ls_stride );
#else
	DAXPY12(nal,
		    f_in  + shift + dag_shift_pls + the_end,
		    f_out + shift + dag_shift_pls );

	DAXPY12(al,
		    f_in  + shift + dag_shift_min + the_end,
		    f_out + shift + dag_shift_min + the_end - ls_stride );
#endif
      }
    
  } // omp parallel

  DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
}

//-------------------------------------------------------------

inline
void mobius_kappa_dslash_5_plus_dag1(Vector *out, 
		       Vector *in, 
		       Float mass,
    		       Dwf *mobius_lib_arg,
		       Float a_five_inv )
{
  //const IFloat two_over_a5 = 2.0 * a_five_inv;
  //const IFloat neg_mass_two_over_a5 = -2.0 * mass * a_five_inv;
  // for mobius, use 1/2*(1+-g5) for projectors
  // also, a_five_inv is a non-trivial factor depending on mobius kappa's or b,c coeff's,
  // depending on when it is called. Don't change name to avoid enbug.
  const IFloat two_over_a5 = 1.0 * a_five_inv;
  const IFloat neg_mass_two_over_a5 = -1.0 * mass * a_five_inv;
  const int local_ls    = GJP.SnodeSites(); 
  const int s_nodes     = GJP.Snodes();
  const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb   = mobius_lib_arg->vol_4d / 2;
  const int ls_stride   = 24 * vol_4d_cb;
  const int max_dex((local_ls-1)*vol_4d_cb);



  if (s_nodes != 1 || s_node_coor !=0) {
    // error!
  }
  
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
    int x;
    int idx;
    IFloat *f_in (reinterpret_cast<IFloat*>(in));
    IFloat *f_out(reinterpret_cast<IFloat*>(out));

    int dag_shift_pls;
    int dag_shift_min;
    int first12_target;
    int second12_target;

    register __m128d al;
    al = _mm_loaddup_pd(&two_over_a5);
    register __m128d nal;
    nal = _mm_loaddup_pd(&neg_mass_two_over_a5);

    register __m128d x0, x1, x2, x3, x4, x5;		\
    register __m128d y0, y1, y2, y3, y4, y5;		\
    register __m128d z0, z1;				\


#if 0    
    if(dag == 1){
#endif
      dag_shift_pls   =  12;
      dag_shift_min   =  0;
      first12_target  = -ls_stride;
      second12_target =  ls_stride+12;
#if 0
    }
    else{
      dag_shift_pls   =  0;
      dag_shift_min   =  12;
      first12_target  =  ls_stride;
      second12_target = -ls_stride+12;
    }
#endif
    // index of the start of the last slice in the fifth dimension
    const int the_end((local_ls-1)*ls_stride);



#ifdef _OPENMP
#pragma omp for 
#endif
    for(x=0; x<vol_4d_cb; x++)
      {
	const int shift(24*x);
#ifdef USE_BLAS
	cblas_daxpy(12,two_over_a5,
		    f_in  + shift + dag_shift_pls,
		    f_out + shift + dag_shift_pls + ls_stride );
	cblas_daxpy(12,neg_mass_two_over_a5,
		    f_in  + shift + dag_shift_min,
		    f_out + shift + dag_shift_min + the_end );
#else
	DAXPY12( al, 	    f_in  + shift + dag_shift_pls,
		    f_out + shift + dag_shift_pls + ls_stride );

	DAXPY12( nal, 	
		    f_in  + shift + dag_shift_min,
		    f_out + shift + dag_shift_min + the_end );
#endif

      }
    
#ifdef _OPENMP   
#pragma omp for  
#endif
    for (idx=vol_4d_cb;idx<max_dex;idx++)
      {
	const int shift(24*idx);
#ifdef USE_BLAS
	cblas_daxpy(12,two_over_a5,
		    f_in  + shift,
		    f_out + shift + first12_target );

	cblas_daxpy(12,two_over_a5,
		    f_in  + shift + 12,
		    f_out + shift + second12_target );
#else
	DAXPY12(al,
		    f_in  + shift,
		    f_out + shift + first12_target );

	DAXPY12(al,
		    f_in  + shift + 12,
		    f_out + shift + second12_target );
#endif
      }
#ifdef _OPENMP    
#pragma omp for 
#endif
    for(x=0; x<vol_4d_cb; x++)
      {
	const int shift(24*x);
#ifdef USE_BLAS
	cblas_daxpy(12,neg_mass_two_over_a5,
		    f_in  + shift + dag_shift_pls + the_end,
		    f_out + shift + dag_shift_pls );

	cblas_daxpy(12,two_over_a5,
		    f_in  + shift + dag_shift_min + the_end,
		    f_out + shift + dag_shift_min + the_end - ls_stride );
#else
	DAXPY12(nal,
		    f_in  + shift + dag_shift_pls + the_end,
		    f_out + shift + dag_shift_pls );

	DAXPY12(al,
		    f_in  + shift + dag_shift_min + the_end,
		    f_out + shift + dag_shift_min + the_end - ls_stride );

#endif
      }
    
  } // omp parallel

  DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
}


void mobius_kappa_dslash_5_plus_a_five(Vector *out, 
			      Vector *in, 
			      Float mass,
			      int dag, 
			      Dwf *mobius_lib_arg,
			      Float a_five_inv )
{
  if (dag == 0) 
    mobius_kappa_dslash_5_plus_dag0(out, in, mass, mobius_lib_arg, a_five_inv);
  else
    mobius_kappa_dslash_5_plus_dag1(out, in, mass, mobius_lib_arg, a_five_inv);
}

void mobius_kappa_dslash_5_plus(Vector *out, 
		       Vector *in, 
		       Float mass,
		       int dag, 
		       Dwf *mobius_lib_arg, Float fact)
{
  if (dag == 0) 
    mobius_kappa_dslash_5_plus_dag0(out, in, mass, mobius_lib_arg, fact);
  else
    mobius_kappa_dslash_5_plus_dag1(out, in, mass, mobius_lib_arg, fact);
}

void mobius_dslash_5_plus(Vector *out, 
		       Vector *in, 
		       Float mass,
		       int dag, 
		       Dwf *mobius_lib_arg)
{
  if (dag == 0) 
    mobius_kappa_dslash_5_plus_dag0(out, in, mass, mobius_lib_arg, GJP.DwfA5Inv());
  else
    mobius_kappa_dslash_5_plus_dag1(out, in, mass, mobius_lib_arg, GJP.DwfA5Inv());
}


CPS_END_NAMESPACE
#endif
