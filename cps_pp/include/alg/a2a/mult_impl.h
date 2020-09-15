#ifndef _MULT_IMPL
#define _MULT_IMPL

#include<alg/a2a/gsl_wrapper.h>

//Options for _mult_impl:
//MULT_IMPL_BASIC        :not using any external libraries
//MULT_IMPL_BLOCK_BASIC  :blocked matrix implementation without external libraries
//MULT_IMPL_GSL          :blocked matrix implementation using GSL BLAS (and I think other BLAS can be slotted in by linking to the appropriate libraries)
//MULT_IMPL_GRID         :using Grid library SIMD intrinsics with a hand-crafted wrapper

#define MULT_IMPL_GSL
//#define MULT_IMPL_GRID

#if defined(MULT_IMPL_BASIC)
#  include<alg/a2a/mult_impl_basic.h>
#elif defined(MULT_IMPL_BLOCK_BASIC)
#  include<alg/a2a/mult_impl_block_basic.h>
#elif defined(MULT_IMPL_GSL)
#  include<alg/a2a/mult_impl_gsl.h>
#elif defined(MULT_IMPL_GRID)
#  include<alg/a2a/mult_impl_grid.h>
#else
#  error Must specify a MULT_IMPL_* in alg/a2a/mult_impl.h
#endif

//Matrix product of meson field pairs
//out(t1,t4) = l(t1,t2) * r(t3,t4)     (The stored timeslices are only used to unpack TimePackedIndex so it doesn't matter if t2 and t3 are thrown away; their indices are contracted over hence the times are not needed)
template<typename mf_Policies, 
	 template <typename> class lA2AfieldL,  template <typename> class lA2AfieldR,
	 template <typename> class rA2AfieldL,  template <typename> class rA2AfieldR
	 >
void mult(A2AmesonField<mf_Policies,lA2AfieldL,rA2AfieldR> &out, const A2AmesonField<mf_Policies,lA2AfieldL,lA2AfieldR> &l, const A2AmesonField<mf_Policies,rA2AfieldL,rA2AfieldR> &r, const bool node_local){
  _mult_impl<mf_Policies,lA2AfieldL,lA2AfieldR,rA2AfieldL,rA2AfieldR>::mult(out,l,r, node_local);
}


#endif
