#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//--------------------------------------------------------------------
//------------------------------------------------------------------
// This routine computes the inverse of the 5th dimension hopping term
// for MOBIUS, plus the diagonal term for the whole 5 dim MOBIUS Dirac operator.
// It is the inverse of 1 - kappa * mobius_dslash_5_plus
// where mobius_dslash_5_plus is the derivative part of the 5th direction
// part of the fermion matrix.
//
// It is part of the 4d odd-even preconditioned MOBIUS operator
//
// Storage order for MOBIUS fermions
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
#include<util/mobius.h>
#include<util/gjp.h>
#include<util/dirac_op.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/smalloc.h>
#include<comms/scu.h>

#include "blas-subs.h"

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





//----------------------------------------
//
//  dag 0 two s-loop version, after loop unrolling
//
void mobius_m5inv_dag0(Vector *inout, 
		       const Float mass,
		       const Dwf *mobius_lib_arg)
{

  int x;
  int s;

  // Initializations
  //------------------------------------------------------------------
  const int ls = GJP.SnodeSites()*GJP.Snodes();
  // const int s_nodes = GJP.Snodes();
  // const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = mobius_lib_arg->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;
  const int f_size = 24 * vol_4d_cb * ls;

  const IFloat two_kappa = -mobius_lib_arg->mobius_kappa_b/mobius_lib_arg->mobius_kappa_c;
  const IFloat inv_two_kappa = 1.0 / two_kappa;
  const Float  inv_d_last = 1.0 / ( 1.0 + pow(two_kappa, ls)*mass); // 1.0 / d_{ls-1}

  Float fact;

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
  
  // FIXME  :   perhaps using array d[i] is faster ?
  fact = - two_kappa * mass * inv_d_last ; // - d0 / d_{ls-1}

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

  //DEBUG_MOBIUS_DSLASH("dag0 first part %e\n", time_elapse());
  
  //-------------------------------------------------------------------
  //  The second backward ls loop both for upper  and downer spinors
  //-------------------------------------------------------------------
  // FIXME: perhaps x loop should also be the reverse order ?
  
  
  fact = - pow(two_kappa,ls-1) * mass * inv_d_last ; // - d_{ls-2} / d_{ls-1}

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

  // s = ls - 1
  //-------------

  for(x=0; x<vol_4d_cb; x+=4) {
    // upper   fou[x,ls-1] *= inv_d_last 
    VEC_TIMESEQU_FLOAT(f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(24+f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);      
    VEC_TIMESEQU_FLOAT(48+f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);  
    VEC_TIMESEQU_FLOAT(72+f_out+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);  
  }

  //DEBUG_MOBIUS_DSLASH("dag0 second part %e\n", time_elapse());
  DiracOp::CGflops+=vol_4d_cb*(ls*96-48);  
  //DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
}

//----------------------------------------
//
//  dag 1 two s-loop version, after loop unrolling
//
void mobius_m5inv_dag1(Vector *inout, 
		    const Float mass,
		    const Dwf *mobius_lib_arg)
{

  int x;
  int s;

// Initializations
//------------------------------------------------------------------
  const int ls = GJP.SnodeSites()*GJP.Snodes();
  // const int s_nodes = GJP.Snodes();
  // const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = mobius_lib_arg->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;
  const int f_size = 24 * vol_4d_cb * ls;

  const IFloat two_kappa = -mobius_lib_arg->mobius_kappa_b/mobius_lib_arg->mobius_kappa_c;
  const IFloat inv_two_kappa = 1.0 / two_kappa;
  const Float  inv_d_last = 1.0 / ( 1.0 + pow(two_kappa, ls)*mass); // 1.0 / d_{ls-1}

  Float fact;
  
  Float* const f_out = (IFloat *) inout;

  //time_elapse();
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

  //DEBUG_MOBIUS_DSLASH("dag1 first part %e\n", time_elapse());
  
  //-------------------------------------------------------------------
  //  The second backward ls loop both for upper  and downer spinors
  //-------------------------------------------------------------------
  // FIXME: perhaps x loop should also be the reverse order ?
  
  
  fact = - pow(two_kappa,ls-1) * mass * inv_d_last ; // - d_{ls-2} / d_{ls-1}

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

  // s = ls - 1
  //-------------
  for(x=0; x<vol_4d_cb; x+=4) {
    // downer  fou[x,ls-1] *= inv_d_last 
    VEC_TIMESEQU_FLOAT(f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(24+f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(48+f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
    VEC_TIMESEQU_FLOAT(72+f_out+12+ 24*(x+vol_4d_cb*(ls-1)), inv_d_last, 12);
  }

  //DEBUG_MOBIUS_DSLASH("dag1 second part %e\n", time_elapse());
  
  //DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
  DiracOp::CGflops+=vol_4d_cb*(ls*96-48);
}






//----------------------------------------

void mobius_m5inv(Vector *inout,
	       Float mass,
	       int dag,
	       Dwf *mobius_lib_arg)
{
  if(dag==0)
    mobius_m5inv_dag0(inout, mass, mobius_lib_arg);
  else 
    mobius_m5inv_dag1(inout, mass, mobius_lib_arg);  
}

//----------------------------------------

void mobius_m5inv(Vector *out, Vector *in,
		  Float mass,
		  int dag,
		  Dwf *mobius_lib_arg)
{
  const int vol_4d_cb = mobius_lib_arg->vol_4d / 2;
  const int f_size = 24 * mobius_lib_arg->ls * vol_4d_cb;

  moveFloat( (Float*)out, (Float*)in, f_size );

  if(dag==0)
    mobius_m5inv_dag0(out, mass, mobius_lib_arg);
  else 
    mobius_m5inv_dag1(out, mass, mobius_lib_arg);  
}




CPS_END_NAMESPACE
