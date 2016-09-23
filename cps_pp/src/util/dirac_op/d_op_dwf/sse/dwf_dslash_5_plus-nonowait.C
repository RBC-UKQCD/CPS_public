#ifdef USE_SSE
//#define USE_BLAS
#include <pmmintrin.h>
//#include <cblas.h>
#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_dwf/sse/dwf_dslash_5_plus-nonowait.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// dwf_dslash_5_plus.C
//
// dwf_dslash_5_plus is the derivative part of the 5th direction
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
#include<util/dwf.h>
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


#if 1
#define DAXPY( _A, _X, _Y, _x, _y )		\
    _x = _mm_load_pd( _X );			\
    _y = _mm_load_pd( _Y );			\
    _x = _mm_add_pd( _y, _mm_mul_pd(_A, _x) );	\
    _mm_store_pd( _Y, _x );			\

#else

#endif



#define DIST_XY 6
#define DAXPY12( _A, _X, _Y )			\
    DAXPY(_A, _X+0,  _Y+0, x0, y0);	\
    DAXPY(_A, _X+2,  _Y+2, x1, y1);		\
    DAXPY(_A, _X+4,  _Y+4, x2, y2);		\
    DAXPY(_A, _X+6,  _Y+6, x3, y3);		\
    DAXPY(_A, _X+8,  _Y+8, x4, y4);		\
    DAXPY(_A, _X+10, _Y+10, x5, y5);		\
    \
    _mm_prefetch((char *)(_X + DIST_XY*24), _MM_HINT_T0);	\
    _mm_prefetch((char *)(_Y + DIST_XY*24), _MM_HINT_T0);	\


#define DAXPYW( _A, _X, _Y, _W, _x, _y, _z )	\
    _x = _mm_load_pd( _X );			\
    _y = _mm_load_pd( _Y );			\
    _z= _mm_mul_pd(_A, _x);			\
    _z = _mm_add_pd( _y, _z);			\
    _mm_store_pd(_Y, _z);			\

#define DIST_XYW 6
#define DAXPYW12( _A, _X, _Y)				\
    DAXPYW(_A, _X+0,  _Y+0, w0, x0, y0, z0);		\
    DAXPYW(_A, _X+2,  _Y+2, w1, x1, y1, z1);		\
    DAXPYW(_A, _X+4,  _Y+4, w2, x2, y2, z0);		\
    DAXPYW(_A, _X+6,  _Y+6, w2, x3, y3, z1);		\
    DAXPYW(_A, _X+8,  _Y+8, w4, x4, y4, z0);		\
    DAXPYW(_A, _X+10, _Y+10, w5, x5, y5,z1);		\
									\
    _mm_prefetch((char *)(_X + DIST_XYW*24), _MM_HINT_T0);		\
    _mm_prefetch((char *)(_Y + DIST_XYW*24), _MM_HINT_T0);		\
    _mm_prefetch((char *)(_X + DIST_XYW*24+8), _MM_HINT_T0);		\
    _mm_prefetch((char *)(_Y + DIST_XYW*24+8), _MM_HINT_T0);		\


#ifdef USE_BLAS
#define DO24(_shift)		    \
    cblas_daxpy(12,two_over_a5,	    \
		f_in  + _shift,				\
		f_out + _shift + first12_target );	\
							\
    cblas_daxpy(12,two_over_a5,				\
		f_in  + _shift + 12,			\
		f_out + _shift + second12_target );	\
    
#else
#define DO24(_al, _x, _y)				\
    DAXPYW12(al, f_in  + _shift,			\
	     f_out + _shift + first12_target);		\
							\
    DAXPYW12(al, f_in  + _shift + 12,			\
	     f_out + _shift + second12_target);		\

#endif

//-----------------------------------------------------------

inline
void dwf_dslash_5_plus_dag0(Vector *out, 
		       Vector *in, 
		       Float mass,
		       Dwf *dwf_lib_arg)
{
  const IFloat two_over_a5 = 2.0 * GJP.DwfA5Inv();
  const IFloat neg_mass_two_over_a5 = -2.0 * mass * GJP.DwfA5Inv();
  const int local_ls    = GJP.SnodeSites(); 
  const int s_nodes     = GJP.Snodes();
  const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb   = dwf_lib_arg->vol_4d / 2;
  const int ls_stride   = 24 * vol_4d_cb;
  const int max_dex((local_ls-1)*vol_4d_cb);

  if (s_nodes != 1 || s_node_coor !=0) {
    // error!
  }
  
#pragma omp parallel default(shared)
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
    register __m128d z0, z1;			\




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



#pragma omp for 
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
    

#pragma omp for 
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

    
#pragma omp for 
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
void dwf_dslash_5_plus_dag1(Vector *out, 
		       Vector *in, 
		       Float mass,
		       Dwf *dwf_lib_arg)
{
  const IFloat two_over_a5 = 2.0 * GJP.DwfA5Inv();
  const IFloat neg_mass_two_over_a5 = -2.0 * mass * GJP.DwfA5Inv();
  const int local_ls    = GJP.SnodeSites(); 
  const int s_nodes     = GJP.Snodes();
  const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb   = dwf_lib_arg->vol_4d / 2;
  const int ls_stride   = 24 * vol_4d_cb;
  const int max_dex((local_ls-1)*vol_4d_cb);



  if (s_nodes != 1 || s_node_coor !=0) {
    // error!
  }
  
#pragma omp parallel default(shared)
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



#pragma omp for 
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
    
   
#pragma omp for  
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
    
#pragma omp for 
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


void dwf_dslash_5_plus(Vector *out, 
		       Vector *in, 
		       Float mass,
		       int dag, 
		       Dwf *dwf_lib_arg)
{
  if (dag == 0) 
    dwf_dslash_5_plus_dag0(out, in, mass, dwf_lib_arg);
  else
    dwf_dslash_5_plus_dag1(out, in, mass, dwf_lib_arg);
}


CPS_END_NAMESPACE
#endif
