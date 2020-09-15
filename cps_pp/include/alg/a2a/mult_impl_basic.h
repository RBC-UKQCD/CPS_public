#ifndef _MULT_IMPL_BASIC_H
#define _MULT_IMPL_BASIC_H

#define DO_REORDER  //Testing on my laptop indicated its better to reorder the matrices to improve cache usage

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

    int work = ni*nk;
    int node_work, node_off; bool do_work;
    getNodeWork(work,node_work,node_off,do_work,node_local);

    typedef typename A2AmesonField<mf_Policies,lA2AfieldL,lA2AfieldR>::RightDilutionType LeftDilutionType;
    typedef typename A2AmesonField<mf_Policies,rA2AfieldL,rA2AfieldR>::LeftDilutionType RightDilutionType;

    ModeContractionIndices<LeftDilutionType,RightDilutionType> j_ind2(l.getColParams());
    
    if(do_work){
      Float time = -dclock();

      modeIndexSet lmodeparams; lmodeparams.time = l.tr;
      modeIndexSet rmodeparams; rmodeparams.time = r.tl;
	
      int nj = j_ind2.getNindices(lmodeparams,rmodeparams);

      //complex mult  re = re*re - im*im, im = re*im + im*re   //6 flops
      //complex add   2 flops

      Float flops_total = Float(ni)*Float(nk)*Float(nj)*8.;

      int jlmap[nj], jrmap[nj];
      for(int j = 0; j < nj; j++)
	j_ind2.getBothIndices(jlmap[j],jrmap[j],j,lmodeparams,rmodeparams);

#  ifndef DO_REORDER

#pragma omp parallel for
      for(int ik = node_off; ik < node_off + node_work; ++ik){
	int i = ik % ni;
	int k = ik / ni;
	for(int j = 0; j < nj; j++)
	  out(i,k) += l(i,jlmap[j]) * r(jrmap[j],k);
      }
    

#  else

      A2AmesonField<mf_Policies,lA2AfieldL,lA2AfieldR> lreord;
      l.colReorder(lreord,jlmap,nj);
      A2AmesonField<mf_Policies,rA2AfieldL,rA2AfieldR> rreord;
      r.rowReorder(rreord,jrmap,nj);

      //A2AmesonField<mf_Policies,rA2AfieldR,rA2AfieldL> rreord_T;
      //rreord.transpose(rreord_T); //more efficient memory access

      static const int lcol_stride = 1;      
      int rrow_stride = rreord.getNcols();

#pragma omp parallel for
      for(int ik = node_off; ik < node_off + node_work; ++ik){
	int i = ik % ni;
	int k = ik / ni;
	ScalarComplexType const* lbase = &lreord(i,0);
	ScalarComplexType const* rbase = &rreord(0,k);
	//std::complex<mf_Complex> const* rbase = &rreord_T(k,0);

	for(int j = 0; j < nj; ++j){
	  out(i,k) += (*lbase)*(*rbase);
	  lbase += lcol_stride;
	  rbase += rrow_stride;
	  //++lbase;
	  //++rbase;
	}

      }
#   endif
      
      time += dclock();
      Float flops_per_sec = flops_total/time;
      if(!UniqueID()) printf("node mult flops/s %g  (time %f total flops %g)\n",flops_per_sec,time,flops_total);
    }
  

    if(!node_local) out.nodeSum();
  }
};


#endif
