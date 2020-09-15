#ifndef _MULT_IMPL_BLOCK_BASIC_H
#define _MULT_IMPL_BLOCK_BASIC_H

//Implementations for meson field contractions
template<typename mf_Policies, 
	 template <typename> class lA2AfieldL,  template <typename> class lA2AfieldR,
	 template <typename> class rA2AfieldL,  template <typename> class rA2AfieldR
	 >
class _mult_impl{ //necessary to avoid an annoying ambigous overload when mesonfield friends mult
public:
  //Matrix product of meson field pairs
  //out(t1,t4) = l(t1,t2) * r(t3,t4)     (The stored timeslices are only used to unpack TimePackedIndex so it doesn't matter if t2 and t3 are thrown away; their indices are contracted over hence the times are not needed)
  static void mult(A2AmesonField<mf_Policies,lA2AfieldL,rA2AfieldR> &out, const A2AmesonField<mf_Policies,lA2AfieldL,lA2AfieldR> &l, const A2AmesonField<mf_Policies,rA2AfieldL,rA2AfieldR> &r, const bool node_local){
    typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
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
    int bmax = 128; //base block size; actual blocks this size or smaller
    int bi = bmax, bj = bmax, bk = bmax;
    while( ni % bi != 0 ) --bi;
    while( nj % bj != 0 ) --bj;
    while( nk % bk != 0 ) --bk;

    //TEST
    //bi = ni/16; bj = nj/16; bk = nk/16;

    int ni0 = ni/bi, nj0 = nj/bj, nk0 = nk/bk;
    if(!UniqueID()) printf("mult sizes %d %d %d block sizes %d %d %d, num blocks %d %d %d\n",ni,nj,nk,bi,bj,bk,ni0,nj0,nk0);
    assert(ni0 * bi == ni);
    assert(nj0 * bj == nj);
    assert(nk0 * bk == nk);
    
    //parallelize ijk
    int work = ni0 * nj0 * nk0;
    int node_work, node_off; bool do_work;
    getNodeWork(work,node_work,node_off,do_work,node_local);

    if(do_work){    
      Float t1 = dclock();

      //complex mult  re = re*re - im*im, im = re*im + im*re   //6 flops
      //complex add   2 flops

      Float flops_total = Float(ni)*Float(nk)*Float(nj)*8.;

      A2AmesonField<mf_Policies,lA2AfieldL,lA2AfieldR> lreord;
      l.colReorder(lreord,jlmap,nj);
      A2AmesonField<mf_Policies,rA2AfieldL,rA2AfieldR> rreord;
      r.rowReorder(rreord,jrmap,nj);
      
      static const int lcol_stride = 1;      
      int rrow_stride = rreord.getNcols();

#pragma omp parallel for
      for(int i0j0k0 = node_off; i0j0k0 < node_off + node_work; ++i0j0k0){
	int rem = i0j0k0;
	int k0 = rem % nk0; rem /= nk0;
	int j0 = rem % nj0; rem /= nj0;
	int i0 = rem;
	i0 *= bi; j0 *= bj; k0 *= bk;

	ScalarComplexType ijblock[bi][bj];
	for(int i=0;i<bi;i++) for(int j=0;j<bj;j++) ijblock[i][j] = lreord(i0+i, j0+j);
	
	//std::complex<mf_Float> jkblock[bj][bk];
	//for(int j=0;j<bj;j++) for(int k=0;k<bk;k++) jkblock[j][k] = rreord(j0+j, k0+k);

	ScalarComplexType kjblock[bk][bj];
	for(int j=0;j<bj;j++) for(int k=0;k<bk;k++) kjblock[k][j] = rreord(j0+j, k0+k); //inplace transpose to speed things up

	for(int i=i0; i<i0+bi; ++i){
	  for(int k=k0; k<k0+bk; ++k){
	    for(int jc = 0; jc < bj; ++jc){
	      //out(i,k) += ijblock[i-i0][jc] * jkblock[jc][k-k0];
	      out(i,k) += ijblock[i-i0][jc] * kjblock[k-k0][jc];
	    }
	  }
	}
      }

      Float t2 = dclock();

      Float flops_per_sec = flops_total/(t2-t1);
      if(!UniqueID()) printf("node mult flops/s %g  (time %f total flops %g)\n",flops_per_sec,t2-t1,flops_total);

    }
    if(!node_local) out.nodeSum();
  }
};

#endif
