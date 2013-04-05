#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//--------------------------------------------------------------------
//------------------------------------------------------------------
// This routine computes the inverse of the 5th dimension hopping term
// for DWF, plus the diagonal term for the whole 5 dim DWF Dirac operator.
// It is the inverse of 1 - kappa * dwf_dslash_5_plus
// where dwf_dslash_5_plus is the derivative part of the 5th direction
// part of the fermion matrix.
//
// It is part of the 4d odd-even preconditioned DWF operator
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



//  for 4d preconditioned operator,
//  |                |
//  | |s-block odd|  |  For odd chckerboard (same for even)
//  | |s-block odd|  |
//  | |s-block odd|  |
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
#include<comms/scu.h>

#include <pmmintrin.h>
#include "blas-subs.h"
#include "sse-subs-axpys.h"


CPS_START_NAMESPACE

/*
  To Chulwoo,

    There are several equivalent version for m5^{-1}.

    The before_loop_unrolling version, which has "two s-loop",
    is corresponding to the one in the note, and would be easy to read and modify.
    
    There are also unroll loop version which is trivially equivalent to the above, but lengthy.
    This is currently used. But you may want to switch to use the versions before unrolling
    as reference implementations. 

    At the end of this file, there are version with four s-loop, which corresponds to earlier version in the note.


    By changing the Float fact into array, we could easily implement the part of Mobius.
    
 */



#if 0
//----------------------------------------
//
//  dag 0 two s-loop version, before loop unrolling
//
inline 
void before_unrolling_dwf_m5inv_dag0( const Vector *inout,
				     const Float mass,
				     const Dwf *dwf_lib_arg)
{

  int x;
  int s;

// Initializations
//------------------------------------------------------------------
  const int ls = GJP.SnodeSites()*GJP.Snodes();
  // const int s_nodes = GJP.Snodes();
  // const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = dwf_lib_arg->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;
  const int f_size = 24 * vol_4d_cb * ls;

  const IFloat two_kappa = 2.0 * dwf_lib_arg->dwf_kappa;
  const IFloat inv_two_kappa = 1.0 / two_kappa;
  const Float  inv_d_last = 1.0 / ( 1.0 + pow(two_kappa, ls)*mass); // 1.0 / d_{ls-1}

  Float fact;
  
  //
  // constant pointers initialized to beginning of Vectors
  // FIXME : perhaps relative pointer computation is faster ?
  //
  //const Float* const f_in  = (IFloat *) in;
  Float* const f_out = (IFloat *) inout;

  time_elapse();

  //-----------------------------------------------------
  //  The first forward ls loop both for upper  and downer spinors
  //----------------------------------------------------

  // s = 0
  //-------
  for(x=0; x<vol_4d_cb; ++x) {
    // downer part  fout[x,ls-1] *= d_last;
    VEC_TIMESEQU_FLOAT(f_out+ 12+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
  }
  
  // FIXME  :   perhaps using array d[i] is faster ?
  fact = - two_kappa * mass * inv_d_last ; // - d0 / d_{ls-1}


  //register __m128d al;
  //al = _mm_loaddup_pd(&two_over_a5);
  //register __m128d nal;
  //nal = _mm_loaddup_pd(&neg_mass_two_over_a5);

  // s = 1 ... ls - 2
  //-------------------
  for(s=0; s<= ls-2 ;++s) {
    for(x=0; x<vol_4d_cb; ++x) {
      // upper part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, f_out + 24*(x+vol_4d_cb*s),f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, f_out + 12+24*(x+vol_4d_cb*s),f_out + 12+24*(x+vol_4d_cb*(ls-1)));
    }
    fact *= two_kappa;
  }

  DEBUG_DWF_DSLASH("first part %e\n", time_elapse());
  
  //-------------------------------------------------------------------
  //  The second backward ls loop both for upper  and downer spinors
  //-------------------------------------------------------------------
  // FIXME: perhaps x loop should also be the reverse order ?
  
  
  fact = - pow(two_kappa,ls-1) * mass * inv_d_last ; // - d_{ls-2} / d_{ls-1}

  // s = ls-2, ... ,  0
  //----------------------
  for(s=ls-2; s >=0 ; --s){
    for(x=0; x<vol_4d_cb; ++x) {
      // upper fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, f_out +24*(x+vol_4d_cb*(ls-1)), f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, f_out+12+24*(x+vol_4d_cb*(s+1)),f_out + 12+24*(x+vol_4d_cb*s));
    }
    //FIXME : which is faster ?
    //fact *= inv_two_kappa;
    fact /= two_kappa;
  }

  // s = ls - 1
  //-------------
  for(x=0; x<vol_4d_cb; ++x) {
    // upper   fou[x,ls-1] *= inv_d_last 
    VEC_TIMESEQU_FLOAT(f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
  }

  DEBUG_DWF_DSLASH("second part %e\n", time_elapse());

  DiracOp::CGflops+=vol_4d_cb*(ls*96-48);  
  //DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
}



//----------------------------------------
//
//  dag 1 two s-loop version, before loop unrolling
//
void before_unrolling_dwf_m5inv_dag1(Vector *inout, 
				     const Float mass,
				     const Dwf *dwf_lib_arg)
{

  int x;
  int s;

// Initializations
//------------------------------------------------------------------
  const int ls = GJP.SnodeSites()*GJP.Snodes();
  // const int s_nodes = GJP.Snodes();
  // const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = dwf_lib_arg->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;
  const int f_size = 24 * vol_4d_cb * ls;

  const IFloat two_kappa = 2.0 * dwf_lib_arg->dwf_kappa;
  const IFloat inv_two_kappa = 1.0 / two_kappa;
  const Float  inv_d_last = 1.0 / ( 1.0 + pow(two_kappa, ls)*mass); // 1.0 / d_{ls-1}

  Float fact;
  
  //
  // constant pointers initialized to beginning of Vectors
  // FIXME : perhaps relative pointer computation is faster ?
  //
  const Float* const f_in  = (IFloat *) inout;

  time_elapse();

  //-----------------------------------------------------
  //  The first forward ls loop both for upper  and downer spinors
  //----------------------------------------------------

  // s = 0
  //-------
  for(x=0; x<vol_4d_cb; ++x) {
    // upper part  fout[x,ls-1] *= d_last;
    VEC_TIMESEQU_FLOAT(f_out+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
  }
  
  // FIXME  :   perhaps using array d[i] is faster ?
  fact = - two_kappa * mass * inv_d_last ; // - d0 / d_{ls-1}


  // s = 1 ... ls - 2
  //-------------------
  for(s=0; s<= ls-2 ;++s) {
    for(x=0; x<vol_4d_cb; ++x) {
      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, f_out +24*(x+vol_4d_cb*s),f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, f_out+ 12 + 24*(x+vol_4d_cb*s),f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));
    }
    fact *= two_kappa;
  }

  //DEBUG_DWF_DSLASH("first part %e\n", time_elapse());
  
  //-------------------------------------------------------------------
  //  The second backward ls loop both for upper  and downer spinors
  //-------------------------------------------------------------------
  // FIXME: perhaps x loop should also be the reverse order ?
  
  
  fact = - pow(two_kappa,ls-1) * mass * inv_d_last ; // - d_{ls-2} / d_{ls-1}

  // s = ls-2, ... ,  0
  //----------------------
  for(s=ls-2; s >=0 ; --s){
    for(x=0; x<vol_4d_cb; ++x) {
      // upper fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, f_out+24*(x+vol_4d_cb*(s+1)),f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, f_out +12+24*(x+vol_4d_cb*(ls-1)), f_out + 12 + 24*(x+vol_4d_cb*s));
    }
    //FIXME : which is faster ?
    //fact *= inv_two_kappa;
    fact /= two_kappa;
  }

  // s = ls - 1
  //-------------
  for(x=0; x<vol_4d_cb; ++x) {
    // downer  fou[x,ls-1] *= inv_d_last 
    VEC_TIMESEQU_FLOAT(f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
  }

  //DEBUG_DWF_DSLASH("second part %e\n", time_elapse());

  DiracOp::CGflops+=vol_4d_cb*(ls*96-48);  
  //DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
}

#endif

//----------------------------------------
//
//  dag 0 two s-loop version, after loop unrolling
//
//  in-place version
//
void dwf_m5inv_dag0(Vector *inout,
		    const Float mass,
		    const Dwf *dwf_lib_arg)
{

  int x;
  int s;

  DECLARE_AXPY_TMPVARS;
  
  register __m128d _mm_fact, _mm_two_kappa;
  
  // Initializations
  //------------------------------------------------------------------
  const int ls = GJP.SnodeSites()*GJP.Snodes();
  // const int s_nodes = GJP.Snodes();
  // const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = dwf_lib_arg->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;
  const int f_size = 24 * vol_4d_cb * ls;

  const IFloat two_kappa = 2.0 * dwf_lib_arg->dwf_kappa;
  const IFloat inv_two_kappa = 1.0 / two_kappa;
  const Float  inv_d_last = 1.0 / ( 1.0 + pow(two_kappa, ls)*mass); // 1.0 / d_{ls-1}

  Float fact;
  
  //
  // constant pointers initialized to beginning of Vectors
  //
  Float* const f_out = (IFloat *) inout;

  //time_elapse();
  
  //-----------------------------------------------------
  //  The first forward ls loop both for upper  and downer spinors
  //----------------------------------------------------

  // s = 0
  //-------
  //#pragma ivdep
  for(x=0; x<vol_4d_cb; x+=4) {
    // downer part  fout[x,ls-1] *= d_last;
    VEC_TIMESEQU_FLOAT(f_out+ 12+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(24+f_out+ 12+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(48+f_out+ 12+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(72+f_out+ 12+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
  }
  
  //DEBUG_DWF_DSLASH("dag0 vec*=float %e\n", time_elapse());
  
  // FIXME  :   perhaps using array d[i] is faster ?
  fact = - two_kappa * mass * inv_d_last ; // - d0 / d_{ls-1}

  _mm_fact = _mm_loaddup_pd(&fact);
  _mm_two_kappa = _mm_loaddup_pd(&two_kappa);

  
#if 0
  // s = 1 ... ls - 2
  //-------------------
  for(s=0; s<= ls-2 ;++s) {
    //#pragma ivdep
    for(x=0; x<vol_4d_cb; x+=4) {
      // upper part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, f_out + 24*(x+vol_4d_cb*s),f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, f_out + 12+24*(x+vol_4d_cb*s),f_out + 12+24*(x+vol_4d_cb*(ls-1)));

      // upper part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, 24+f_out + 24*(x+vol_4d_cb*s),24+f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, 24+f_out + 12+24*(x+vol_4d_cb*s),24+f_out + 12+24*(x+vol_4d_cb*(ls-1)));

      // upper part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, 48+f_out + 24*(x+vol_4d_cb*s),48+f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, 48+f_out + 12+24*(x+vol_4d_cb*s),48+f_out + 12+24*(x+vol_4d_cb*(ls-1)));

      // upper part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, 72+f_out + 24*(x+vol_4d_cb*s),72+f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, 72+f_out + 12+24*(x+vol_4d_cb*s),72+f_out + 12+24*(x+vol_4d_cb*(ls-1)));
    }
    fact *= two_kappa;
  }
#else
  // s = 1 ... ls - 2
  //-------------------
  for(s=0; s<= ls-2 ;++s) {
    //#pragma ivdep
    for(x=0; x<vol_4d_cb; x+=4) {
      // upper part  fout[x,s+1] +=  two_kappa fout[x,s]
      DAXPY12( _mm_two_kappa, f_out + 24*(x+vol_4d_cb*s),f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  fact * fout[x,s]
      DAXPY12( _mm_fact, f_out + 12+24*(x+vol_4d_cb*s),f_out + 12+24*(x+vol_4d_cb*(ls-1)));

      // upper part  fout[x,s+1] +=  two_kappa fout[x,s]
      DAXPY12( _mm_two_kappa, 24+f_out + 24*(x+vol_4d_cb*s),24+f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  fact * fout[x,s]
      DAXPY12( _mm_fact, 24+f_out + 12+24*(x+vol_4d_cb*s),24+f_out + 12+24*(x+vol_4d_cb*(ls-1)));

      // upper part  fout[x,s+1] +=  two_kappa fout[x,s]
      DAXPY12( _mm_two_kappa, 48+f_out + 24*(x+vol_4d_cb*s),48+f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  fact * fout[x,s]
      DAXPY12( _mm_fact, 48+f_out + 12+24*(x+vol_4d_cb*s),48+f_out + 12+24*(x+vol_4d_cb*(ls-1)));

      // upper part  fout[x,s+1] +=  two_kappa fout[x,s]
      DAXPY12( _mm_two_kappa, 72+f_out + 24*(x+vol_4d_cb*s),72+f_out + 24*(x+vol_4d_cb*(s+1)));
      // downer part fout[x,ls-1] +=  fact * fout[x,s]
      DAXPY12( _mm_fact, 72+f_out + 12+24*(x+vol_4d_cb*s),72+f_out + 12+24*(x+vol_4d_cb*(ls-1)));
    }
    //fact *= two_kappa;
    _mm_fact = _mm_mul_pd(_mm_two_kappa, _mm_fact);
  }
#endif  

  //DEBUG_DWF_DSLASH("dag0 first part %e\n", time_elapse());
  
  //-------------------------------------------------------------------
  //  The second backward ls loop both for upper  and downer spinors
  //-------------------------------------------------------------------
  // FIXME: perhaps x loop should also be the reverse order ?
  
  
  fact = - pow(two_kappa,ls-1) * mass * inv_d_last ; // - d_{ls-2} / d_{ls-1}

  _mm_fact = _mm_loaddup_pd(&fact);
  
#if 0
  // s = ls-2, ... ,  0
  //----------------------
  for(s=ls-2; s >=0 ; --s){
    //#pragma ivdep
    for(x=0; x<vol_4d_cb; x+=4) {
      // upper fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, f_out +24*(x+vol_4d_cb*(ls-1)), f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, f_out+12+24*(x+vol_4d_cb*(s+1)),f_out + 12+24*(x+vol_4d_cb*s));
      
      // upper fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, 24+f_out +24*(x+vol_4d_cb*(ls-1)), 24+f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, 24+f_out+12+24*(x+vol_4d_cb*(s+1)),24+f_out + 12+24*(x+vol_4d_cb*s));
      
      // upper fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, 48+f_out +24*(x+vol_4d_cb*(ls-1)), 48+f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, 48+f_out+12+24*(x+vol_4d_cb*(s+1)),48+f_out + 12+24*(x+vol_4d_cb*s));

      // upper fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, 72+f_out +24*(x+vol_4d_cb*(ls-1)), 72+f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, 72+f_out+12+24*(x+vol_4d_cb*(s+1)),72+f_out + 12+24*(x+vol_4d_cb*s));

    }
    //FIXME : which is faster ?
    //fact *= inv_two_kappa;
    fact /= two_kappa;
  }
#else
    // s = ls-2, ... ,  0
  //----------------------
  for(s=ls-2; s >=0 ; --s){
    //#pragma ivdep
    for(x=0; x<vol_4d_cb; x+=4) {
      // upper fout[x,s] += fact* fout[x,ls-1]
      DAXPY12( _mm_fact, f_out +24*(x+vol_4d_cb*(ls-1)), f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += two_kappa fout[x,s+1]
      DAXPY12( _mm_two_kappa, f_out+12+24*(x+vol_4d_cb*(s+1)),f_out + 12+24*(x+vol_4d_cb*s));
      
      // upper fout[x,s] += fact* fout[x,ls-1]
      DAXPY12( _mm_fact, 24+f_out +24*(x+vol_4d_cb*(ls-1)), 24+f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += two_kappa fout[x,s+1]
      DAXPY12( _mm_two_kappa, 24+f_out+12+24*(x+vol_4d_cb*(s+1)),24+f_out + 12+24*(x+vol_4d_cb*s));
      
      // upper fout[x,s] += fact* fout[x,ls-1]
      DAXPY12( _mm_fact, 48+f_out +24*(x+vol_4d_cb*(ls-1)), 48+f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += two_kappa fout[x,s+1]
      DAXPY12( _mm_two_kappa, 48+f_out+12+24*(x+vol_4d_cb*(s+1)),48+f_out + 12+24*(x+vol_4d_cb*s));

      // upper fout[x,s] += fact* fout[x,ls-1]
      DAXPY12( _mm_fact, 72+f_out +24*(x+vol_4d_cb*(ls-1)), 72+f_out + 24*(x+vol_4d_cb*s));
      // downer fout[x,s] += two_kappa fout[x,s+1]
      DAXPY12( _mm_two_kappa, 72+f_out+12+24*(x+vol_4d_cb*(s+1)),72+f_out + 12+24*(x+vol_4d_cb*s));

    }
    //FIXME : which is faster ?
    //fact *= inv_two_kappa;
    //fact /= two_kappa;
    _mm_fact = _mm_div_pd( _mm_fact, _mm_two_kappa );
  }
#endif
  // s = ls - 1
  //-------------

  for(x=0; x<vol_4d_cb; x+=4) {
    // upper   fou[x,ls-1] *= inv_d_last 
    VEC_TIMESEQU_FLOAT(f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(24+f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);      
    VEC_TIMESEQU_FLOAT(48+f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);  
    VEC_TIMESEQU_FLOAT(72+f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);  
  }

  //DEBUG_DWF_DSLASH("dag0 second part %e\n", time_elapse());

  DiracOp::CGflops+=vol_4d_cb*(ls*96-48);
}

//----------------------------------------
//
//  dag 1 two s-loop version, after loop unrolling
//
//  in-place version
void dwf_m5inv_dag1(Vector *inout,
		    const Float mass,
		    const Dwf *dwf_lib_arg)
{

  int x;
  int s;

  DECLARE_AXPY_TMPVARS;
  
  register __m128d _mm_fact, _mm_two_kappa;
  
  // Initializations
//------------------------------------------------------------------
  const int ls = GJP.SnodeSites()*GJP.Snodes();
  // const int s_nodes = GJP.Snodes();
  // const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = dwf_lib_arg->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;
  const int f_size = 24 * vol_4d_cb * ls;

  const IFloat two_kappa = 2.0 * dwf_lib_arg->dwf_kappa;
  const IFloat inv_two_kappa = 1.0 / two_kappa;
  const Float  inv_d_last = 1.0 / ( 1.0 + pow(two_kappa, ls)*mass); // 1.0 / d_{ls-1}

  Float fact;
  
  // constant pointers initialized to beginning of Vectors

  Float* const f_out = (IFloat *) inout;
  
  time_elapse();

  //-----------------------------------------------------
  //  The first forward ls loop both for upper  and downer spinors
  //----------------------------------------------------

  // s = 0
  //-------
  for(x=0; x<vol_4d_cb; x+=4) {
    // upper part  fout[x,ls-1] *= d_last;
    VEC_TIMESEQU_FLOAT(f_out+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(24+f_out+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(48+f_out+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(72+f_out+24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
  }
  
  // FIXME  :   perhaps using array d[i] is faster ?
  fact = - two_kappa * mass * inv_d_last ; // - d0 / d_{ls-1}

  _mm_two_kappa = _mm_loaddup_pd(&two_kappa);
  _mm_fact = _mm_loaddup_pd(&fact);
  

#if 0
  // s = 1 ... ls - 2
  //-------------------
  for(s=0; s<= ls-2 ;++s) {
    for(x=0; x<vol_4d_cb; x+=4) {
      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, f_out +24*(x+vol_4d_cb*s),f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, f_out+ 12 + 24*(x+vol_4d_cb*s),f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));

      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, 24+f_out +24*(x+vol_4d_cb*s),24+f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, 24+f_out+ 12 + 24*(x+vol_4d_cb*s),24+f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));

      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, 48+f_out +24*(x+vol_4d_cb*s),48+f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, 48+f_out+ 12 + 24*(x+vol_4d_cb*s),48+f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));

      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      AXPY(12, fact, 72+f_out +24*(x+vol_4d_cb*s),72+f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  two_kappa fout[x,s]
      AXPY(12, two_kappa, 72+f_out+ 12 + 24*(x+vol_4d_cb*s),72+f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));
    }
    fact *= two_kappa;
  }
#else
    // s = 1 ... ls - 2
  //-------------------
  for(s=0; s<= ls-2 ;++s) {
    for(x=0; x<vol_4d_cb; x+=4) {
      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      DAXPY12( _mm_fact, f_out +24*(x+vol_4d_cb*s),f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  two_kappa fout[x,s]
      DAXPY12( _mm_two_kappa, f_out+ 12 + 24*(x+vol_4d_cb*s),f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));

      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      DAXPY12( _mm_fact, 24+f_out +24*(x+vol_4d_cb*s),24+f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  two_kappa fout[x,s]
      DAXPY12( _mm_two_kappa, 24+f_out+ 12 + 24*(x+vol_4d_cb*s),24+f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));

      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      DAXPY12( _mm_fact, 48+f_out +24*(x+vol_4d_cb*s),48+f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  two_kappa fout[x,s]
      DAXPY12( _mm_two_kappa, 48+f_out+ 12 + 24*(x+vol_4d_cb*s),48+f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));

      // upper part fout[x,ls-1] +=  fact * fout[x,s]
      DAXPY12( _mm_fact, 72+f_out +24*(x+vol_4d_cb*s),72+f_out +24*(x+vol_4d_cb*(ls-1)));
      // downer part  fout[x,s+1] +=  two_kappa fout[x,s]
      DAXPY12( _mm_two_kappa, 72+f_out+ 12 + 24*(x+vol_4d_cb*s),72+f_out+ 12 + 24*(x+vol_4d_cb*(s+1)));
    }
    //fact *= two_kappa;
    _mm_fact = _mm_mul_pd( _mm_two_kappa, _mm_fact);
  }
#endif
  //DEBUG_DWF_DSLASH("dag1 first part %e\n", time_elapse());
  
  //-------------------------------------------------------------------
  //  The second backward ls loop both for upper  and downer spinors
  //-------------------------------------------------------------------
  // FIXME: perhaps x loop should also be the reverse order ?
  
  
  fact = - pow(two_kappa,ls-1) * mass * inv_d_last ; // - d_{ls-2} / d_{ls-1}

  _mm_fact = _mm_loaddup_pd(&fact);
#if 0  
  // s = ls-2, ... ,  0
  //----------------------
  for(s=ls-2; s >=0 ; --s){
    for(x=0; x<vol_4d_cb; x+=4) {
      // upper fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, f_out+24*(x+vol_4d_cb*(s+1)),f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, f_out +12+24*(x+vol_4d_cb*(ls-1)), f_out + 12 + 24*(x+vol_4d_cb*s));

      // upper fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, 24+f_out+24*(x+vol_4d_cb*(s+1)),24+f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, 24+f_out +12+24*(x+vol_4d_cb*(ls-1)), 24+f_out + 12 + 24*(x+vol_4d_cb*s));

      // upper fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, 48+f_out+24*(x+vol_4d_cb*(s+1)),48+f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, 48+f_out +12+24*(x+vol_4d_cb*(ls-1)), 48+f_out + 12 + 24*(x+vol_4d_cb*s));

      // upper fout[x,s] += two_kappa fout[x,s+1]
      AXPY(12, two_kappa, 72+f_out+24*(x+vol_4d_cb*(s+1)),72+f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += fact* fout[x,ls-1]
      AXPY(12, fact, 72+f_out +12+24*(x+vol_4d_cb*(ls-1)), 72+f_out + 12 + 24*(x+vol_4d_cb*s));
}
    //FIXME : which is faster ?
    //fact *= inv_two_kappa;
    fact /= two_kappa;
  }
#else
    // s = ls-2, ... ,  0
  //----------------------
  for(s=ls-2; s >=0 ; --s){
    for(x=0; x<vol_4d_cb; x+=4) {
      // upper fout[x,s] += two_kappa fout[x,s+1]
      DAXPY12( _mm_two_kappa, f_out+24*(x+vol_4d_cb*(s+1)),f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += fact* fout[x,ls-1]
      DAXPY12( _mm_fact, f_out +12+24*(x+vol_4d_cb*(ls-1)), f_out + 12 + 24*(x+vol_4d_cb*s));

      // upper fout[x,s] += two_kappa fout[x,s+1]
      DAXPY12( _mm_two_kappa, 24+f_out+24*(x+vol_4d_cb*(s+1)),24+f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += fact* fout[x,ls-1]
      DAXPY12( _mm_fact, 24+f_out +12+24*(x+vol_4d_cb*(ls-1)), 24+f_out + 12 + 24*(x+vol_4d_cb*s));

      // upper fout[x,s] += two_kappa fout[x,s+1]
      DAXPY12( _mm_two_kappa, 48+f_out+24*(x+vol_4d_cb*(s+1)),48+f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += fact* fout[x,ls-1]
      DAXPY12( _mm_fact, 48+f_out +12+24*(x+vol_4d_cb*(ls-1)), 48+f_out + 12 + 24*(x+vol_4d_cb*s));

      // upper fout[x,s] += two_kappa fout[x,s+1]
      DAXPY12( _mm_two_kappa, 72+f_out+24*(x+vol_4d_cb*(s+1)),72+f_out +24*(x+vol_4d_cb*s));
      // downer fout[x,s] += fact* fout[x,ls-1]
      DAXPY12( _mm_fact, 72+f_out +12+24*(x+vol_4d_cb*(ls-1)), 72+f_out + 12 + 24*(x+vol_4d_cb*s));
}
    //FIXME : which is faster ?
    //fact *= inv_two_kappa;
    //fact /= two_kappa;
    _mm_fact = _mm_div_pd(_mm_fact, _mm_two_kappa);
  }
#endif
  
  // s = ls - 1
  //-------------
  for(x=0; x<vol_4d_cb; x+=4) {
    // downer  fou[x,ls-1] *= inv_d_last 
    VEC_TIMESEQU_FLOAT(f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(24+f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(48+f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(72+f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
  }

  //DEBUG_DWF_DSLASH("dag1 second part %e\n", time_elapse());
  //DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
  DiracOp::CGflops+=vol_4d_cb*(ls*96-48);
}


//----------------------------------------
//
//  dag 0 Modified  four s-loop version
//
void four_s_loop_dwf_m5inv_dag0(Vector *out, 
		    const Vector *in,
		    const Float mass,
		    const Dwf *dwf_lib_arg)

{

  const int dag=0;
  
  int x;
  int s;

// Initializations
//------------------------------------------------------------------
  const int local_ls = GJP.SnodeSites(); 
  const int s_nodes = GJP.Snodes();
  const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = dwf_lib_arg->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;
  const int f_size = 24 * vol_4d_cb * local_ls;
  IFloat *f_in;
  IFloat *f_out;
  IFloat *f_temp;
  const IFloat two_kappa = 2.0 * dwf_lib_arg->dwf_kappa;
  Float fact1;
  const Float fact2  = 1.0/(1.0+mass*pow(two_kappa,local_ls));
  Float fact;
  


  time_elapse();
  // [1 + gamma_5] term (if dag=1 [1 - gamma_5] term)
  //
  //------------------------------------------------------------------
  // pointers initialized to beginning of Vectors
  f_in  = (IFloat *) in;
  f_out = (IFloat *) out;

#define RELIDX
  // copy upper 2 spin of v (in) into u (out) 
  for(s=0; s < local_ls; s++){
    for(x=0; x<vol_4d_cb; x++){


      MOVE_FLOAT(f_out + 24*(x+vol_4d_cb*s), f_in + 24*(x+vol_4d_cb*s), 12);

      
      //moveFloat(f_out , f_in , 12);
      //f_in  =  f_in + 24;
      //f_out = f_out + 24;
    }
  }

  // finish s=0 term
  //----------------
  // (re)set pointers
  //f_out = (IFloat *) out;
  f_temp = 0; //dwf_lib_arg->frm_tmp3;//not initialized!

  ERR.NotImplemented("","m5inv","dwf_lib_arg->frm_tmp3 is not allocated");
  
  //initialize f_temp (s=0)
  for(x=0; x<vol_4d_cb; x++){
#ifndef USE_BLAS
    moveFloat(f_temp+24*x, f_out+24*x, 12);
#else
    cblas_dcopy(12, f_out+24*x, f_temp+24*x);
#endif
    //f_temp  =  f_temp + 24;
    //f_out = f_out + 24;
  }
  //f_temp = f_temp-ls_stride;//back to the beginning
  //f_in += ls_stride;
  //f_out is at s=1, so don't do anything
  printf("copy %e\n", time_elapse());
  
  for(s=1; s < local_ls; s++){
    fact1 = -mass*pow(two_kappa,local_ls-s);

    for(x=0; x<vol_4d_cb; x+=4){

#ifndef USE_BLAS      
      fTimesV1PlusV2(f_temp, fact1, f_out+24*(x+vol_4d_cb*s), f_temp+24*x, 12);
      fTimesV1PlusV2(f_temp, fact1, f_out+24*(1+x+vol_4d_cb*s), f_temp+24*(x+1), 12);
      fTimesV1PlusV2(f_temp, fact1, f_out+24*(2+x+vol_4d_cb*s), f_temp+24*(x+2), 12);
      fTimesV1PlusV2(f_temp, fact1, f_out+24*(3+x+vol_4d_cb*s), f_temp+24*(x+3), 12);
#else
      cblas_daxpy(12, fact1, f_out+24*(x+vol_4d_cb*s), 1, f_temp+24*x, 1);
      cblas_daxpy(12, fact1, f_out+24*(1+x+vol_4d_cb*s), 1, f_temp+24*(x+1), 1);
      cblas_daxpy(12, fact1, f_out+24*(2+x+vol_4d_cb*s), 1, f_temp+24*(x+2), 1);
      cblas_daxpy(12, fact1, f_out+24*(3+x+vol_4d_cb*s), 1, f_temp+24*(x+3), 1);
#endif
      
      //f_out = f_out + 24;
      //f_temp = f_temp + 24;
      //f_in = f_in + 24;
    }
    //f_temp = f_temp-ls_stride;//back to the beginning
  }
  // reset f_out & f_in
  //f_out = f_out - ls_stride*local_ls;
  //f_in = f_in - ls_stride*local_ls;
    
  printf("first part %e\n",time_elapse());
  
  //finish s=0 f_temp, end of back sub
  //fact2 = 1.0/(1.0+mass*pow(two_kappa,local_ls));
  for(x=0; x<vol_4d_cb; x+=4){
#ifndef USE_BLAS
    vecTimesEquFloat(f_temp, fact2, 12);
    f_temp = f_temp + 24;
    vecTimesEquFloat(f_temp, fact2, 12);
    f_temp = f_temp + 24;
    vecTimesEquFloat(f_temp, fact2, 12);
    f_temp = f_temp + 24;
    vecTimesEquFloat(f_temp, fact2, 12);
    f_temp = f_temp + 24;
#else
    cblas_dscal(12, fact2, f_temp,1);
    cblas_dscal(12, fact2, f_temp+24,1);
    cblas_dscal(12, fact2, f_temp+48,1);
    cblas_dscal(12, fact2, f_temp+72,1);
    f_temp = f_temp + 96;
#endif
  }
  // reset f_temp
  f_temp=f_temp-ls_stride;
  // copy temp s=0 into out 
  for(x=0; x<vol_4d_cb; x++){
    moveFloat(f_out, f_temp, 12);
    f_out = f_out + 24;
    f_temp = f_temp + 24;
  }
  // Do forward substitution (s=0 already done)
  for(s=1; s < local_ls; s++){
    f_temp=f_out-ls_stride;

    for(x=0; x<vol_4d_cb; x+=4){
#ifndef USE_BLAS
      fTimesV1PlusV2(f_out, two_kappa, f_temp, f_out, 12);
      f_out = f_out + 24;f_temp = f_temp + 24;
      fTimesV1PlusV2(f_out, two_kappa, f_temp, f_out, 12);
      f_out = f_out + 24;f_temp = f_temp + 24;
      fTimesV1PlusV2(f_out, two_kappa, f_temp, f_out, 12);
      f_out = f_out + 24;f_temp = f_temp + 24;
      fTimesV1PlusV2(f_out, two_kappa, f_temp, f_out, 12);
      f_out = f_out + 24;f_temp = f_temp + 24;
#else
      cblas_daxpy(12, two_kappa, f_temp,1, f_out,1);
      cblas_daxpy(12, two_kappa, f_temp+24,1, f_out+24,1);
      cblas_daxpy(12, two_kappa, f_temp+48,1, f_out+48,1);
      cblas_daxpy(12, two_kappa, f_temp+72,1, f_out+72,1);
      f_out = f_out + 96;f_temp = f_temp + 96;
#endif
    }
  }

  printf("second part %e\n",time_elapse());
  
  // end PR * M5,R^-1

  //    return;
  // [1 - gamma_5] term (if dag=1 [1 - gamma_5] term)
  //-------------------------------------------------

  f_in  = (IFloat *) in;
  f_out = (IFloat *) out;
  // lower two spin components
  if(dag == 0){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }
  // start backward substitution
  // copy v (in) into u (out) 
  for(s=0; s < local_ls; s++){
    for(x=0; x<vol_4d_cb; x++){
      moveFloat(f_out, f_in, 12);
      f_in  =  f_in + 24;
      f_out = f_out + 24;
    }
  }

  // finish intermediate vector, store in f_temp=f_out
  //----------------
  f_out = (IFloat *) out;
  f_temp = f_out;
  //  if(dag == 0){
    f_out = f_out + 12;
    f_temp = f_temp + 12;
    //  }

  f_out = f_out + (local_ls)*ls_stride;
  f_temp = f_temp + (local_ls-1)*ls_stride;
  for(s=0; s<local_ls-1; s++){//just counting, pointers decremented below
    fact = two_kappa;
    for(x=0; x<vol_4d_cb; x++){
      f_out = f_out - 24;
      f_temp = f_temp - 24;
      fTimesV1PlusV2(f_temp, fact, f_out, f_temp, 12);
    }
  }

  //Now forward substitute
  f_out = (IFloat *) out;
  //  if(dag == 0){
    f_out = f_out + 12;
    //  }
  //f_temp = f_temp - ls_stride;

  fact = 1.0/(1.0+pow(two_kappa,local_ls)*mass);
  for(x=0; x<vol_4d_cb; x++){
    vecTimesEquFloat(f_temp, fact, 12);
    f_temp = f_temp + 24;
  }
  f_temp = f_temp - ls_stride;

  // Do forward substitution
  f_out = f_out+ls_stride;
  for(s=1; s < local_ls; s++){
    fact = -mass*pow(two_kappa,local_ls-s);
    for(x=0; x<vol_4d_cb; x++){
      fTimesV1PlusV2(f_out, fact, f_temp, f_out, 12);
      f_out = f_out + 24;
      f_temp = f_temp + 24;
    }
    // back to beginning
    f_temp = f_temp - ls_stride;
  }

  // end PL * M5,L^-1

  //DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
}

//----------------------------------------
//
//  dag 1  Modified four s-loop version
//
void four_s_loop_dwf_m5inv_dag1(Vector *out, 
	       Vector *in,
	       Float mass,
	       Dwf *dwf_lib_arg)
{

  const int dag=1;
  
  int x;
  int s;

// Initializations
//------------------------------------------------------------------
  const int local_ls = GJP.SnodeSites(); 
  const int s_nodes = GJP.Snodes();
  const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = dwf_lib_arg->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;
  const int f_size = 24 * vol_4d_cb * local_ls;
  IFloat *f_in;
  IFloat *f_out;
  IFloat *f_temp;
  const IFloat two_kappa = 2.0 * dwf_lib_arg->dwf_kappa;
  Float fact1;
  const Float fact2  = 1.0/(1.0+mass*pow(two_kappa,local_ls));
  Float fact;
  


  // [1 + gamma_5] term (if dag=1 [1 - gamma_5] term)
  //
  //------------------------------------------------------------------
  // pointers initialized to beginning of Vectors
  f_in  = (IFloat *) in;
  f_out = (IFloat *) out;
  if(dag == 1){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }
  // copy upper 2 spin of v (in) into u (out) 
  for(s=0; s < local_ls; s++){
    for(x=0; x<vol_4d_cb; x++){
      moveFloat(f_out, f_in, 12);
      f_in  =  f_in + 24;
      f_out = f_out + 24;
    }
  }

  // finish s=0 term
  //----------------
  // (re)set pointers
  f_out = (IFloat *) out;
  f_temp = 0;//dwf_lib_arg->frm_tmp3;//not initialized!
  ERR.NotImplemented("","m5inv","dwf_lib_arg->frm_tmp3 is not allocated");
  
  if(dag == 1){
    f_out  =  f_out + 12;
    f_temp = f_temp + 12;
  }

  //initialize f_temp (s=0)
  for(x=0; x<vol_4d_cb; x++){
#ifndef USE_BLAS
    moveFloat(f_temp, f_out, 12);
#else
    cblas_dcopy(12, f_out, f_temp);
#endif
    f_temp  =  f_temp + 24;
    f_out = f_out + 24;
  }
  f_temp = f_temp-ls_stride;//back to the beginning
  f_in += ls_stride;
  //f_out is at s=1, so don't do anything
  for(s=1; s < local_ls; s++){
    fact1 = -mass*pow(two_kappa,local_ls-s);
    for(x=0; x<vol_4d_cb; x++){
#ifndef USE_BLAS      
      fTimesV1PlusV2(f_temp, fact1, f_out, f_temp, 12);
#else
      cblas_daxpy(12, fact1, f_out, 1, f_temp, 1);
#endif
      f_out = f_out + 24;
      f_temp = f_temp + 24;
      f_in = f_in + 24;
    }
    f_temp = f_temp-ls_stride;//back to the beginning
  }
  // reset f_out & f_in
  f_out = f_out - ls_stride*local_ls;
  f_in = f_in - ls_stride*local_ls;
    
  //finish s=0 f_temp, end of back sub
  //fact2 = 1.0/(1.0+mass*pow(two_kappa,local_ls));
  for(x=0; x<vol_4d_cb; x++){
#ifndef USE_BLAS
    vecTimesEquFloat(f_temp, fact2, 12);
#else
    cblas_dscal(12, fact2, f_temp,1);
#endif
    f_temp = f_temp + 24;
  }
  // reset f_temp
  f_temp=f_temp-ls_stride;
  // copy temp s=0 into out 
  for(x=0; x<vol_4d_cb; x++){
    moveFloat(f_out, f_temp, 12);
    f_out = f_out + 24;
    f_temp = f_temp + 24;
  }
  // Do forward substitution (s=0 already done)
  for(s=1; s < local_ls; s++){
    f_temp=f_out-ls_stride;

    for(x=0; x<vol_4d_cb; x++){
#ifndef USE_BLAS
      fTimesV1PlusV2(f_out, two_kappa, f_temp, f_out, 12);
#else
      cblas_daxpy(12, two_kappa, f_temp,1, f_out,1);
#endif
      f_out = f_out + 24;
      f_temp = f_temp + 24;
    }
  }

  // end PR * M5,R^-1

  //    return;
  // [1 - gamma_5] term (if dag=1 [1 - gamma_5] term)
  //-------------------------------------------------

  f_in  = (IFloat *) in;
  f_out = (IFloat *) out;
  // lower two spin components
  if(dag == 0){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }
  // start backward substitution
  // copy v (in) into u (out) 
  for(s=0; s < local_ls; s++){
    for(x=0; x<vol_4d_cb; x++){
      moveFloat(f_out, f_in, 12);
      f_in  =  f_in + 24;
      f_out = f_out + 24;
    }
  }

  // finish intermediate vector, store in f_temp=f_out
  //----------------
  f_out = (IFloat *) out;
  f_temp = f_out;
  if(dag == 0){
    f_out = f_out + 12;
    f_temp = f_temp + 12;
  }

  f_out = f_out + (local_ls)*ls_stride;
  f_temp = f_temp + (local_ls-1)*ls_stride;
  for(s=0; s<local_ls-1; s++){//just counting, pointers decremented below
    fact = two_kappa;
    for(x=0; x<vol_4d_cb; x++){
      f_out = f_out - 24;
      f_temp = f_temp - 24;
      fTimesV1PlusV2(f_temp, fact, f_out, f_temp, 12);
    }
  }

  //Now forward substitute
  f_out = (IFloat *) out;
  if(dag == 0){
    f_out = f_out + 12;
  }
  //f_temp = f_temp - ls_stride;

  fact = 1.0/(1.0+pow(two_kappa,local_ls)*mass);
  for(x=0; x<vol_4d_cb; x++){
    vecTimesEquFloat(f_temp, fact, 12);
    f_temp = f_temp + 24;
  }
  f_temp = f_temp - ls_stride;

  // Do forward substitution
  f_out = f_out+ls_stride;
  for(s=1; s < local_ls; s++){
    fact = -mass*pow(two_kappa,local_ls-s);
    for(x=0; x<vol_4d_cb; x++){
      fTimesV1PlusV2(f_out, fact, f_temp, f_out, 12);
      f_out = f_out + 24;
      f_temp = f_temp + 24;
    }
    // back to beginning
    f_temp = f_temp - ls_stride;
  }

  // end PL * M5,L^-1

  //DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
}



//----------------------------------------

void dwf_m5inv(Vector *inout,
	       Float mass,
	       int dag,
	       Dwf *dwf_lib_arg)
{
  if(dag==0)
    dwf_m5inv_dag0(inout, mass, dwf_lib_arg);
  else 
    dwf_m5inv_dag1(inout, mass, dwf_lib_arg);  
}

//----------------------------------------

void dwf_m5inv(Vector *out, Vector *in,
	       Float mass,
	       int dag,
	       Dwf *dwf_lib_arg)
{
  const int vol_4d_cb = dwf_lib_arg->vol_4d / 2;
  const int f_size = 24 * dwf_lib_arg->ls * vol_4d_cb;

  moveFloat( (Float*)out, (Float*)in, f_size );

  if(dag==0)
    dwf_m5inv_dag0(out, mass, dwf_lib_arg);
  else 
    dwf_m5inv_dag1(out, mass, dwf_lib_arg);  
}


//----------------------------------------
//
//  The original four s-loop version by Tom
//
void four_s_loop_dwf_m5inv(Vector *out, 
			   Vector *in,
			   Float mass,
			   int dag,
			   Dwf *dwf_lib_arg)
{

  if( GJP.Snodes()!=1) ERR.General("DiracOpDwf", "dwf_m5inv", "Does not work for spread-out Ls\n");

  int x;
  int s;

// Initializations
//------------------------------------------------------------------
  int local_ls = GJP.SnodeSites(); 
  int s_nodes = GJP.Snodes();
  int s_node_coor = GJP.SnodeCoor();
  int vol_4d_cb = dwf_lib_arg->vol_4d / 2;
  int ls_stride = 24 * vol_4d_cb;
  int f_size = 24 * vol_4d_cb * local_ls;
  IFloat *f_in;
  IFloat *f_out;
  IFloat *f_temp;
  IFloat two_kappa = 2.0 * dwf_lib_arg->dwf_kappa;
  Float fact;
// [1 + gamma_5] term (if dag=1 [1 - gamma_5] term)
//
//------------------------------------------------------------------
  // pointers initialized to beginning of Vectors
  f_in  = (IFloat *) in;
  f_out = (IFloat *) out;
  if(dag == 1){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }
  // copy upper 2 spin of v (in) into u (out) 
  for(s=0; s < local_ls; s++){
    for(x=0; x<vol_4d_cb; x++){
      moveFloat(f_out, f_in, 12);
      f_in  =  f_in + 24;
      f_out = f_out + 24;
    }
  }

  // finish s=0 term
  //----------------
  // (re)set pointers
  f_out = (IFloat *) out;
  f_temp = 0;//dwf_lib_arg->frm_tmp3;//not initialized!
  ERR.NotImplemented("","m5inv","dwf_lib_arg->frm_tmp3 is not allocated");
  if(dag == 1){
    f_out  =  f_out + 12;
    f_temp = f_temp + 12;
  }

  //initialize f_temp (s=0)
  for(x=0; x<vol_4d_cb; x++){
      moveFloat(f_temp, f_out, 12);
      f_temp  =  f_temp + 24;
      f_out = f_out + 24;
  }
  f_temp = f_temp-ls_stride;//back to the beginning
  //f_out is at s=1, so don't do anything
  for(s=1; s < local_ls; s++){
    fact = -mass*pow(two_kappa,local_ls-s);
    for(x=0; x<vol_4d_cb; x++){
      fTimesV1PlusV2(f_temp, fact, f_out, f_temp, 12);
      f_out = f_out + 24;
      f_temp = f_temp + 24;
    }
    f_temp = f_temp-ls_stride;//back to the beginning
  }
  // reset f_out
  f_out = f_out - ls_stride*local_ls;

  //finish s=0 f_temp, end of back sub
  fact = 1.0/(1.0+mass*pow(two_kappa,local_ls));
  for(x=0; x<vol_4d_cb; x++){
    vecTimesEquFloat(f_temp, fact, 12);
    f_temp = f_temp + 24;
  }
  // reset f_temp
  f_temp=f_temp-ls_stride;
  // copy temp s=0 into out 
  for(x=0; x<vol_4d_cb; x++){
    moveFloat(f_out, f_temp, 12);
    f_out = f_out + 24;
    f_temp = f_temp + 24;
  }
  // Do forward substitution (s=0 already done)
  for(s=1; s < local_ls; s++){
    f_temp=f_out-ls_stride;
    for(x=0; x<vol_4d_cb; x++){
      fTimesV1PlusV2(f_out, two_kappa, f_temp, f_out, 12);
      f_out = f_out + 24;
      f_temp = f_temp + 24;
    }
  }

  // end PR * M5,R^-1

  // [1 - gamma_5] term (if dag=1 [1 - gamma_5] term)
  //-------------------------------------------------

  f_in  = (IFloat *) in;
  f_out = (IFloat *) out;
  // lower two spin components
  if(dag == 0){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }
  // start backward substitution
  // copy v (in) into u (out) 
  for(s=0; s < local_ls; s++){
    for(x=0; x<vol_4d_cb; x++){
      moveFloat(f_out, f_in, 12);
      f_in  =  f_in + 24;
      f_out = f_out + 24;
    }
  }

  // finish intermediate vector, store in f_temp=f_out
  //----------------
  f_out = (IFloat *) out;
  f_temp = f_out;
  if(dag == 0){
    f_out = f_out + 12;
    f_temp = f_temp + 12;
  }

  f_out = f_out + (local_ls)*ls_stride;
  f_temp = f_temp + (local_ls-1)*ls_stride;
  for(s=0; s<local_ls-1; s++){//just counting, pointers decremented below
    fact = two_kappa;
    for(x=0; x<vol_4d_cb; x++){
      f_out = f_out - 24;
      f_temp = f_temp - 24;
      fTimesV1PlusV2(f_temp, fact, f_out, f_temp, 12);
    }
  }

  //Now forward substitute
  f_out = (IFloat *) out;
  if(dag == 0){
    f_out = f_out + 12;
  }
  //f_temp = f_temp - ls_stride;

  fact = 1.0/(1.0+pow(two_kappa,local_ls)*mass);
  for(x=0; x<vol_4d_cb; x++){
    vecTimesEquFloat(f_temp, fact, 12);
    f_temp = f_temp + 24;
  }
  f_temp = f_temp - ls_stride;

  // Do forward substitution
  f_out = f_out+ls_stride;
  for(s=1; s < local_ls; s++){
    fact = -mass*pow(two_kappa,local_ls-s);
    for(x=0; x<vol_4d_cb; x++){
      fTimesV1PlusV2(f_out, fact, f_temp, f_out, 12);
      f_out = f_out + 24;
      f_temp = f_temp + 24;
    }
    // back to beginning
    f_temp = f_temp - ls_stride;
  }

  // end PL * M5,L^-1

  //DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
}



CPS_END_NAMESPACE


