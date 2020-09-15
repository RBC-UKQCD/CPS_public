#ifndef _MULT_IMPL_GSL_H
#define _MULT_IMPL_GSL_H

CPS_END_NAMESPACE
#include<alg/a2a/gsl_wrapper.h>
CPS_START_NAMESPACE

//Implementations for meson field contractions
template<typename mf_Policies, 
	 template <typename> class lA2AfieldL,  template <typename> class lA2AfieldR,
	 template <typename> class rA2AfieldL,  template <typename> class rA2AfieldR
	 >
class _mult_impl{ //necessary to avoid an annoying ambigous overload when mesonfield friends mult
public:
  typedef gsl_wrapper<typename mf_Policies::ScalarComplexType::value_type> gw;

  //Matrix product of meson field pairs
  //out(t1,t4) = l(t1,t2) * r(t3,t4)     (The stored timeslices are only used to unpack TimePackedIndex so it doesn't matter if t2 and t3 are thrown away; their indices are contracted over hence the times are not needed)

  inline static int nearest_divisor(const int of, const int base_divisor){
    //printf("nearest_divisor of %d, base_divisor %d\n", of, base_divisor); fflush(stdout);
    assert(base_divisor > 0);
    if(of % base_divisor == 0) return base_divisor;

    int nearest_below = base_divisor;    
    bool no_nearest_below = false;
    while(of % nearest_below != 0){ 
      --nearest_below;
      if(nearest_below == 0){ no_nearest_below = true; break; }
    }
    int nearest_above = base_divisor;
    bool no_nearest_above = false;
    while(of % nearest_above !=0){
      ++nearest_above;
      if(nearest_above == of){ no_nearest_above = true; break; }
    }
    if(no_nearest_above && no_nearest_below) return of;
    if(no_nearest_below) return nearest_above;
    if(no_nearest_above) return nearest_below;
    
    int sep_above = nearest_above - base_divisor;
    int sep_below = base_divisor - nearest_below;
    return sep_above < sep_below ? nearest_above : nearest_below;
  }
  
  static void mult(A2AmesonField<mf_Policies,lA2AfieldL,rA2AfieldR> &out, const A2AmesonField<mf_Policies,lA2AfieldL,lA2AfieldR> &l, const A2AmesonField<mf_Policies,rA2AfieldL,rA2AfieldR> &r, const bool node_local){
    typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
    typedef typename ScalarComplexType::value_type mf_Float;
    
    assert( (void*)&out != (void*)&l || (void*)&out != (void*)&r );

    if(! l.getColParams().paramsEqual( r.getRowParams() ) ){
      if(!UniqueID()){
	printf("mult():  Illegal matrix product: underlying vector parameters must match\n"); fflush(stdout);
	std::cout << "left-column: " << l.getColParams().print() << "\n";
	std::cout << "right-row: " << r.getRowParams().print() << "\n";
	std::cout.flush();
      }
      exit(-1);
    }

    out.setup(l.getRowParams(),r.getColParams(), l.tl, r.tr ); //zeroes output, so safe to re-use
  
    int ni = l.getNrows();
    int nk = r.getNcols();

    typedef typename A2AmesonField<mf_Policies,lA2AfieldL,lA2AfieldR>::RightDilutionType LeftDilutionType;
    typedef typename A2AmesonField<mf_Policies,rA2AfieldL,rA2AfieldR>::LeftDilutionType RightDilutionType;

    ModeContractionIndices<LeftDilutionType,RightDilutionType> j_ind2(l.getColParams()); //these maps could be cached somewhere
    
    modeIndexSet lmodeparams; lmodeparams.time = l.tr;
    modeIndexSet rmodeparams; rmodeparams.time = r.tl;
    
    int nj = j_ind2.getNindices(lmodeparams,rmodeparams);

    int jlmap[nj], jrmap[nj];
    for(int j = 0; j < nj; j++)
      j_ind2.getBothIndices(jlmap[j],jrmap[j],j,lmodeparams,rmodeparams);

    //Try a blocked matrix multiply
    //Because ni, nj are different and not necessarily multiples of a common blocking we need to dynamically choose the block size
    int nodes = 1;
    for(int i=0;i<5;i++) nodes *= GJP.Nodes(i);

    int compute_elements = omp_get_max_threads() * ( node_local ? 1 : nodes );

    //Want the total number of blocks to be close to the number of compute elements = (number of nodes)*(number of threads)
    //We shouldn't just take the cubed-root though because quite often the number of indices differs substantially
    //We want   ni0 * nj0 * nk0 = nodes
    //and the ratios to be approximately the same between the number of blocks and the number of indices
    //Take ratios wrt smallest so these are always >=1

    int smallest = ni;  
    if(nj < smallest) smallest = nj;
    if(nk < smallest) smallest = nk;

    int ratios[3] = {ni/smallest, nj/smallest, nk/smallest};
    int base = (int)pow( compute_elements/ratios[0]/ratios[1]/ratios[2], 1/3.);  //compute_element
    if(!base) ++base;
    
    int ni0 = nearest_divisor(ni, ratios[0]*base);
    int nj0 = nearest_divisor(nj, ratios[1]*base);
    int nk0 = nearest_divisor(nk, ratios[2]*base);

    assert(ni % ni0 == 0);
    assert(nj % nj0 == 0);
    assert(nk % nk0 == 0);
    
    int bi = ni/ni0;
    int bj = nj/nj0; 
    int bk = nk/nk0;
    
    //parallelize ijk
    int work = ni0 * nj0 * nk0;
    int node_work, node_off; bool do_work;
    getNodeWork(work,node_work,node_off,do_work,node_local);

    //if(!UniqueID()) printf("mult sizes %d %d %d block sizes %d %d %d, num blocks %d %d %d. Work %d, node_work %d\n",ni,nj,nk,bi,bj,bk,ni0,nj0,nk0,work,node_work);

    if(do_work){    
      Float t1 = dclock();

      //complex mult  re = re*re - im*im, im = re*im + im*re   //6 flops
      //complex add   2 flops

      Float flops_total = Float(ni)*Float(nk)*Float(nj)*8.;

      A2AmesonField<mf_Policies,lA2AfieldL,lA2AfieldR> lreord;
      A2AmesonField<mf_Policies,rA2AfieldL,rA2AfieldR> rreord;
#ifndef MEMTEST_MODE
      r.rowReorder(rreord,jrmap,nj);
      l.colReorder(lreord,jlmap,nj);
#endif
      
      typename gw::matrix_complex *lreord_gsl = gw::matrix_complex_alloc(ni,nj);
      typename gw::matrix_complex *rreord_gsl = gw::matrix_complex_alloc(nj,nk);
      
#ifndef MEMTEST_MODE
      
#pragma omp parallel for
      for(int i=0;i<ni;i++)
	for(int j=0;j<nj;j++){
	  const ScalarComplexType & el = lreord(i, j);
	  mf_Float *el_gsl = (mf_Float*)gw::matrix_complex_ptr(lreord_gsl,i,j);
	  *(el_gsl++) = std::real(el);
	  *(el_gsl) = std::imag(el);
	}

#pragma omp parallel for
      for(int j=0;j<nj;j++)
	for(int k=0;k<nk;k++){
	  const ScalarComplexType & el = rreord(j, k);
	  mf_Float *el_gsl = (mf_Float*)gw::matrix_complex_ptr(rreord_gsl,j,k);
	  *(el_gsl++) = std::real(el);
	  *(el_gsl) = std::imag(el);
	}
      
#endif
      
      static const int lcol_stride = 1;      
      int rrow_stride = rreord.getNcols();

#pragma omp parallel for
      for(int i0j0k0 = node_off; i0j0k0 < node_off + node_work; ++i0j0k0){
	int rem = i0j0k0;
	int k0 = rem % nk0; rem /= nk0;
	int j0 = rem % nj0; rem /= nj0;
	int i0 = rem;
	i0 *= bi; j0 *= bj; k0 *= bk;

	typename gw::complex tmp;
	typename gw::matrix_complex *tmp_out = gw::matrix_complex_alloc(bi,bk);

	typename gw::matrix_complex_const_view ijblock_view = gw::matrix_complex_const_submatrix(lreord_gsl,i0,j0,bi,bj);
	typename gw::matrix_complex_const_view jkblock_view = gw::matrix_complex_const_submatrix(rreord_gsl,j0,k0,bj,bk);

	const typename gw::matrix_complex *const ijblock = &ijblock_view.matrix; //gw::matrix_complex_alloc(bi,bj);
	const typename gw::matrix_complex *const jkblock = &jkblock_view.matrix;  //gw::matrix_complex_alloc(bj,bk);
	
	typename gw::complex one; GSL_SET_COMPLEX(&one,1.0,0.0);
	typename gw::complex zero; GSL_SET_COMPLEX(&zero,0.0,0.0);

#ifndef MEMTEST_MODE
	gw::matrix_complex_set_zero(tmp_out);
	gw::blas_gemm(CblasNoTrans, CblasNoTrans, one, ijblock, jkblock, zero, tmp_out);
	
	for(int i=0;i<bi;i++) 
	  for(int k=0;k<bk;k++){
	    mf_Float const* el = (mf_Float const*)gw::matrix_complex_ptr(tmp_out,i,k);
	    mf_Float(&out_el)[2] = reinterpret_cast<mf_Float(&)[2]>(out(i0+i,k0+k));
#pragma omp atomic
	    out_el[0] += *(el++);
#pragma omp atomic
	    out_el[1] += *(el);
	  }
#endif
	
	gw::matrix_complex_free(tmp_out);
      }


      Float t2 = dclock();

      Float flops_per_sec = flops_total/(t2-t1);
      //if(!UniqueID()) printf("node mult flops/s %g  (time %f total flops %g)\n",flops_per_sec,t2-t1,flops_total);

      gw::matrix_complex_free(lreord_gsl);
      gw::matrix_complex_free(rreord_gsl);

    }
    Float time = -dclock();
    if(!node_local) out.nodeSum();
    time += dclock();
    //if(!UniqueID()) printf("mult comms time %g s\n",time);
  }
  
};

#endif
