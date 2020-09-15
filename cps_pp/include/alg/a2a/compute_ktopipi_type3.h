#ifndef _COMPUTE_KTOPIPI_TYPE3_H
#define _COMPUTE_KTOPIPI_TYPE3_H

CPS_START_NAMESPACE

//TYPE 3 and MIX 3
//Each contraction of this type is made up of different trace combinations of two objects (below [but not in the code!] for simplicity we ignore the fact that the two vectors in the 
//meson fields are allowed to vary in position relative to each other):
//1) \prop^L(x_op,x_pi1) S_2 \prop^L(x_pi1,x_pi2) S_2 \prop^L(x_pi2,x_K) \gamma^5 \prop^H(x_K,x_op)
//2) \prop^L(x_op,x_op)   OR   \prop^H(x_op,x_op)

//We use g5-hermiticity on the strange prop in part 1 (but not part 2 where it appears)
//1) \prop^L(x_op,x_pi1) S_2 \prop^L(x_pi1,x_pi2) S_2 \prop^L(x_pi2,x_K)  [ \prop^H(x_op,x_K) ]^dag \gamma^5
// = vL(x_op) [[ wL^dag(x_pi1) S_2 vL(x_pi1) ]] [[ wL^dag(x_pi2) S_2 vL(x_pi2) ]] [[ wL^dag(x_K) wH(x_K) ]] vH^dag(x_op) \gamma^5
//  where [[ ]] indicate meson fields
  
//2) vL(x_op) wL^dag(x_op)   or  vH(x_op) wH^dag(x_op)


//Run inside threaded environment

template<typename mf_Policies>
void ComputeKtoPiPiGparity<mf_Policies>::type3_contract(ResultsContainerType &result, const int t_K, const int t_dis, const int thread_id, 
						     const SCFmat part1[2], const SCFmat &part2_L, const SCFmat &part2_H){
#ifndef MEMTEST_MODE
  static const int con_off = 13; //index of first contraction in set

  for(int mu=0;mu<4;mu++){ //sum over mu here
    for(int gcombidx=0;gcombidx<8;gcombidx++){
      const SCFmat &G1 = Gamma1<ComplexType>(gcombidx,mu);
      const SCFmat &G2 = Gamma2<ComplexType>(gcombidx,mu);

      for(int pt1_pion=0; pt1_pion<2; pt1_pion++){  //which pion comes first in part 1?
	SCFmat G1_pt1 = part1[pt1_pion]; //= G1*part1[pt1_pion];
	multGammaLeft(G1_pt1,1,gcombidx,mu);

	CPScolorMatrix<ComplexType> tr_sf_G1_pt1 = G1_pt1.SpinFlavorTrace();

	SCFmat G2_pt2_L = part2_L; //= G2*part2_L;
	multGammaLeft(G2_pt2_L,2,gcombidx,mu);

	CPScolorMatrix<ComplexType> tr_sf_G2_pt2_L = G2_pt2_L.SpinFlavorTrace();

	SCFmat G2_pt2_H = part2_H; // = G2*part2_H;
	multGammaLeft(G2_pt2_H,2,gcombidx,mu);

	CPScolorMatrix<ComplexType> tr_sf_G2_pt2_H = G2_pt2_H.SpinFlavorTrace();

	SCFmat ctrans_G2_pt2_L(G2_pt2_L); //speedup by transposing part 1
	ctrans_G2_pt2_L.TransposeColor();
		
	CPSspinMatrix<CPSflavorMatrix<ComplexType> > tr_c_G1_pt1 = G1_pt1.ColorTrace();

#define C(IDX) result(t_K,t_dis,IDX-con_off,gcombidx,thread_id)	      

	//First 6 have a light-quark loop
	C(13) += G1_pt1.Trace() * G2_pt2_L.Trace();
	C(14) += Trace( tr_sf_G1_pt1 , Transpose(tr_sf_G2_pt2_L) );
	C(15) += Trace( tr_sf_G1_pt1 , tr_sf_G2_pt2_L );
	C(16) += Trace( G1_pt1 , G2_pt2_L );
	C(17) += Trace( G1_pt1 , ctrans_G2_pt2_L );
	C(18) += Trace( tr_c_G1_pt1 , G2_pt2_L.ColorTrace() );
	      
	//Second 4 have strange loop
	C(19) += G1_pt1.Trace() * G2_pt2_H.Trace();	      
	C(20) += Trace( tr_sf_G1_pt1 , tr_sf_G2_pt2_H );
	C(21) += Trace( G1_pt1 , G2_pt2_H );
	C(22) += Trace( tr_c_G1_pt1 , G2_pt2_H.ColorTrace() );

#undef C 
      }
    }
  }
#endif
}

template<typename mf_Policies>
void ComputeKtoPiPiGparity<mf_Policies>::type3_compute_mfproducts(std::vector<std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > > &con_pi1_pi2_k,
							       std::vector<std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > > &con_pi2_pi1_k,							     
							       const std::vector<int> &tsep_k_pi, const int tsep_pion, const int tstep, const std::vector<ThreeMomentum> &p_pi_1_all,
							       const std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > &mf_kaon, MesonFieldMomentumContainer<mf_Policies> &mf_pions,
							       const int Lt, const int tpi_sampled, const int ntsep_k_pi){
  //Form the product of the three meson fields
  //con_*_*_k = [[ wL^dag(y) S_2 vL(y) ]] [[ wL^dag(z) S_2 vL(z) ]] [[ wL^dag(x_K) wH(x_K) ]]
  if(!UniqueID()){ printf("Computing con_*_*_k\n"); fflush(stdout); }

  resize_2d(con_pi1_pi2_k,tpi_sampled,ntsep_k_pi);
  resize_2d(con_pi2_pi1_k,tpi_sampled,ntsep_k_pi);

  A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> tmp_WV;
  A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> tmp_pi1pi2, tmp_pi2pi1;

  int nmom = p_pi_1_all.size();
  for(int pidx=0;pidx<nmom;pidx++){ //Average over momentum orientations
    const ThreeMomentum &p_pi_1 = p_pi_1_all[pidx];
    ThreeMomentum p_pi_2 = -p_pi_1;

    std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> > &mf_pi1 = mf_pions.get(p_pi_1); //*mf_pi1_ptr;
    std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> > &mf_pi2 = mf_pions.get(p_pi_2); //*mf_pi2_ptr;
#ifdef NODE_DISTRIBUTE_MESONFIELDS
    nodeGetMany(2,&mf_pi1,&mf_pi2);
#endif
    //nodeGetPionMf(mf_pi1,mf_pi2);

    //for(int tpi1=0;tpi1<Lt;tpi1+=tstep){ //my sensible ordering
    for(int t_pi1_lin = 1; t_pi1_lin <= Lt; t_pi1_lin += tstep){ //Daiqian's weird ordering
      int tpi1 = modLt(t_pi1_lin,Lt);

      int tpi1_idx = tpi1/tstep;
      int tpi2 = modLt(tpi1  + tsep_pion, Lt);

      for(int tkp = 0; tkp < ntsep_k_pi; tkp++){
	int tk = modLt(tpi1 - tsep_k_pi[tkp], Lt);
      
	A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> *into_pi1_pi2 = pidx == 0 ? &con_pi1_pi2_k[tpi1_idx][tkp] : &tmp_pi1pi2;
	A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> *into_pi2_pi1 = pidx == 0 ? &con_pi2_pi1_k[tpi1_idx][tkp] : &tmp_pi2pi1;
	  
	mult(tmp_WV, mf_pi1[tpi1], mf_pi2[tpi2]); //we can re-use this from type 2 I think
	mult(*into_pi1_pi2, tmp_WV, mf_kaon[tk]);
	  
	mult(tmp_WV, mf_pi2[tpi2], mf_pi1[tpi1]);
	mult(*into_pi2_pi1, tmp_WV, mf_kaon[tk]);
	  
	if(pidx > 0){
	  con_pi1_pi2_k[tpi1_idx][tkp].plus_equals(tmp_pi1pi2,true);
	  con_pi2_pi1_k[tpi1_idx][tkp].plus_equals(tmp_pi2pi1,true);
	}

	//NB time coordinate of con_*_* is the time coordinate of pi1 (that closest to the kaon)
      }
    }
#ifdef NODE_DISTRIBUTE_MESONFIELDS
    nodeDistributeMany(2,&mf_pi1,&mf_pi2);
#endif

    //nodeDistributePionMf(mf_pi1,mf_pi2);
  }
  if(nmom > 1)
    for(int t=0;t<tpi_sampled;t++)
      for(int tkp = 0; tkp < ntsep_k_pi; tkp++){
	con_pi1_pi2_k[t][tkp].times_equals(1./nmom);
	con_pi2_pi1_k[t][tkp].times_equals(1./nmom);
      }
  if(!UniqueID()){ printf("Finished computing con_pi_pi_k\n"); fflush(stdout); }
}

template<typename mf_Policies>
void ComputeKtoPiPiGparity<mf_Policies>::type3_mult_vMv_setup(mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorWfftw,A2AvectorV> &mult_vMv_split_part1_pi1_pi2,
							   mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorWfftw,A2AvectorV> &mult_vMv_split_part1_pi2_pi1,
							   const A2AvectorV<mf_Policies> & vL, const A2AvectorV<mf_Policies> & vH,
							   const std::vector<std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > > &con_pi1_pi2_k,
							   const std::vector<std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > > &con_pi2_pi1_k,
							   const int top_loc, const int t_pi1_idx, const int tkp){

  //Split the vector-mesonfield outer product into two stages where in the first we reorder the mesonfield to optimize cache hits
  int top_glb = top_loc  + GJP.TnodeCoor()*GJP.TnodeSites();
  mult_vMv_split_part1_pi1_pi2.setup(vL,con_pi1_pi2_k[t_pi1_idx][tkp],vH, top_glb);
  mult_vMv_split_part1_pi2_pi1.setup(vL,con_pi2_pi1_k[t_pi1_idx][tkp],vH, top_glb);
}



template<typename mf_Policies>
void ComputeKtoPiPiGparity<mf_Policies>::type3_precompute_part1(SCFmatVector &mult_vMv_contracted_part1_pi1_pi2,
							     SCFmatVector &mult_vMv_contracted_part1_pi2_pi1,
							     mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorWfftw,A2AvectorV> &mult_vMv_split_part1_pi1_pi2,
							     mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorWfftw,A2AvectorV> &mult_vMv_split_part1_pi2_pi1,
							     const int top_loc, const int t_pi1_idx, const int tkp){

  //Contract on all 3d sites on this node with fixed operator time coord top_glb into a canonically ordered output vector
  int top_glb = top_loc  + GJP.TnodeCoor()*GJP.TnodeSites();

  mult_vMv_split_part1_pi1_pi2.contract(mult_vMv_contracted_part1_pi1_pi2, false, true);
  mult_vMv_split_part1_pi1_pi2.free_mem(); //we don't need these any more
  
  mult_vMv_split_part1_pi2_pi1.contract(mult_vMv_contracted_part1_pi2_pi1, false, true);
  mult_vMv_split_part1_pi2_pi1.free_mem();
}


//This version averages over multiple pion momentum configurations. Use to project onto A1 representation at run-time. Saves a lot of time!
//This version also overlaps computation for multiple K->pi separations. Result should be an array of ResultsContainerType the same size as the vector 'tsep_k_pi'
template<typename mf_Policies>
void ComputeKtoPiPiGparity<mf_Policies>::type3(ResultsContainerType result[], MixDiagResultsContainerType mix3[],
					    const std::vector<int> &tsep_k_pi, const int &tsep_pion, const int &tstep, const std::vector<ThreeMomentum> &p_pi_1_all, 
					    const std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > &mf_kaon, MesonFieldMomentumContainer<mf_Policies> &mf_pions,
					    const A2AvectorV<mf_Policies> & vL, const A2AvectorV<mf_Policies> & vH, 
					    const A2AvectorW<mf_Policies> & wL, const A2AvectorW<mf_Policies> & wH){

  SCFmat mix3_Gamma[2];
  mix3_Gamma[0].unit().pr(F0).gr(-5);
  mix3_Gamma[1].unit().pr(F1).gr(-5).timesMinusOne();
  
  const int Lt = GJP.Tnodes()*GJP.TnodeSites();
  const int tpi_sampled = Lt/tstep;
  const int ntsep_k_pi = tsep_k_pi.size();
  
  static const int n_contract = 10; //ten type3 diagrams
  static const int con_off = 13; //index of first contraction in set
  const int nthread = omp_get_max_threads();
  
  for(int tkp=0;tkp<tsep_k_pi.size();tkp++){
    result[tkp].resize(n_contract,nthread); //it will be thread-reduced before this method ends
    mix3[tkp].resize(nthread);
  }
  
  const int size_3d = vL.getMode(0).nodeSites(0)*vL.getMode(0).nodeSites(1)*vL.getMode(0).nodeSites(2);
  
  //Form the product of the three meson fields
  //con_*_*_k = [[ wL^dag(y) S_2 vL(y) ]] [[ wL^dag(z) S_2 vL(z) ]] [[ wL^dag(x_K) wH(x_K) ]]
  std::vector<std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > > con_pi1_pi2_k; // [tpi1_idx][tkp]
  std::vector<std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > > con_pi2_pi1_k;  
  type3_compute_mfproducts(con_pi1_pi2_k,con_pi2_pi1_k,tsep_k_pi,tsep_pion,tstep,p_pi_1_all,mf_kaon,mf_pions,Lt,tpi_sampled,ntsep_k_pi);
  
  //Determine which local operator timeslices are actually used
  std::vector< std::vector<bool> > node_top_used(tpi_sampled);
  for(int t_pi1_lin = 1; t_pi1_lin <= Lt; t_pi1_lin += tstep){ //Daiqian's weird ordering
    int t_pi1 = modLt(t_pi1_lin,Lt);   int t_pi1_idx = t_pi1 / tstep;
    getUsedTimeslices(node_top_used[t_pi1_idx],tsep_k_pi,t_pi1);
  }
  
  for(int top_loc = 0; top_loc < GJP.TnodeSites(); top_loc++){
    const int top_glb = top_loc  + GJP.TnodeCoor()*GJP.TnodeSites();
    
    //Construct part 2 (independent of kaon position):
    //vL(x_op) wL^dag(x_op)   or  vH(x_op) wH^dag(x_op)
    //Loop over Q_i insertion location. Each node naturally has its own sublattice to work on. Thread over sites in usual way
    SCFmatVector part2_L(size_3d); //[x3d]
    SCFmatVector part2_H(size_3d); //[x3d]

#pragma omp parallel for
    for(int xop3d_loc = 0; xop3d_loc < size_3d; xop3d_loc++){
      int thread_id = omp_get_thread_num();
      
      //Construct part 2 (independent of kaon position):
      //vL(x_op) wL^dag(x_op)   or  vH(x_op) wH^dag(x_op)
      mult(part2_L[xop3d_loc], vL, wL, xop3d_loc, top_loc, false, true);
      mult(part2_H[xop3d_loc], vH, wH, xop3d_loc, top_loc, false, true);
    }

    //for(int t_pi1 = 0; t_pi1 < Lt; t_pi1+=tstep){ //my sensible ordering
    for(int t_pi1_lin = 1; t_pi1_lin <= Lt; t_pi1_lin += tstep){ //Daiqian's weird ordering
      int t_pi1 = modLt(t_pi1_lin,Lt);
      
      int t_pi1_idx = t_pi1/tstep;
      int t_pi2 = modLt(t_pi1 + tsep_pion, Lt);
      
      if(!node_top_used[t_pi1_idx][top_loc]) continue; //if this timeslice is not used for any K->pi sep, skip it
      
      for(int tkp = 0; tkp < ntsep_k_pi; tkp++){	  
	int t_K = modLt(t_pi1 - tsep_k_pi[tkp], Lt);
	int t_dis = modLt(top_glb - t_K, Lt); //distance between kaon and operator is the output time coordinate
	
	if(t_dis >= tsep_k_pi[tkp] || t_dis == 0) continue; //don't bother computing operator insertion locations outside of the region between the kaon and first pion or on top of either operator

	//Split the vector-mesonfield outer product into two stages where in the first we reorder the mesonfield to optimize cache hits
#ifndef DISABLE_TYPE3_SPLIT_VMV
	mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorWfftw,A2AvectorV> mult_vMv_split_part1_pi1_pi2;
	mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorWfftw,A2AvectorV> mult_vMv_split_part1_pi2_pi1;
	type3_mult_vMv_setup(mult_vMv_split_part1_pi1_pi2,mult_vMv_split_part1_pi2_pi1,vL,vH,con_pi1_pi2_k,con_pi2_pi1_k,top_loc,t_pi1_idx,tkp);
	
# ifndef DISABLE_TYPE3_PRECOMPUTE
	//Contract on all 3d sites on this node with fixed operator time coord top_glb into a canonically ordered output vector
	SCFmatVector mult_vMv_contracted_part1_pi1_pi2;  //[x3d];
	SCFmatVector mult_vMv_contracted_part1_pi2_pi1;  //[x3d];
	type3_precompute_part1(mult_vMv_contracted_part1_pi1_pi2, mult_vMv_contracted_part1_pi2_pi1, mult_vMv_split_part1_pi1_pi2, mult_vMv_split_part1_pi2_pi1, top_loc,t_pi1_idx,tkp);
# endif
#endif	
	//Now loop over Q_i insertion location. Each node naturally has its own sublattice to work on. Thread over sites in usual way
#pragma omp parallel for
	for(int xop3d_loc = 0; xop3d_loc < size_3d; xop3d_loc++){
	  int thread_id = omp_get_thread_num();
      
	  //Construct part 1:
	  // = vL(x_op) [[ wL^dag(x_pi1) S_2 vL(x_pi1) ]] [[ wL^dag(x_pi2) S_2 vL(x_pi2) ]] [[ wL^dag(x_K) wH(x_K) ]] vH^dag(x_op) \gamma^5
	  SCFmat part1[2]; 
	  
#if defined(DISABLE_TYPE3_SPLIT_VMV)
	  mult(part1[0], vL, con_pi1_pi2_k[t_pi1_idx][tkp], vH, xop3d_loc, top_loc, false, true);
	  mult(part1[1], vL, con_pi2_pi1_k[t_pi1_idx][tkp], vH, xop3d_loc, top_loc, false, true);
#elif defined(DISABLE_TYPE3_PRECOMPUTE)
	  mult_vMv_split_part1_pi1_pi2.contract(part1[0],xop3d_loc, false, true);
	  mult_vMv_split_part1_pi2_pi1.contract(part1[1],xop3d_loc, false, true);
#else
	  part1[0] = mult_vMv_contracted_part1_pi1_pi2[xop3d_loc];
	  part1[1] = mult_vMv_contracted_part1_pi2_pi1[xop3d_loc];
#endif

	  part1[0].gr(-5);
	  part1[1].gr(-5);
	  
	  type3_contract(result[tkp], t_K, t_dis, thread_id, part1, part2_L[xop3d_loc], part2_H[xop3d_loc]);
	  
	  //Compute mix3 diagram
	  //These are identical to the type3 diagrams but without the internal quark loop, and with the vertex replaced with a pseudoscalar vertex
	  for(int mix3_gidx=0; mix3_gidx<2; mix3_gidx++){
	    for(int pt1_pion=0; pt1_pion<2; pt1_pion++){
#ifndef MEMTEST_MODE
#define M mix3[tkp](t_K,t_dis,mix3_gidx,thread_id)
	      M += Trace( part1[pt1_pion], mix3_Gamma[mix3_gidx] );
#undef M
#endif
	    }
	  }
	}//xop3d loop
      }//tkp loop    
    }//t_pi1 loop
  }//top_loc loop

  for(int tkp = 0; tkp < ntsep_k_pi; tkp++){
    result[tkp].threadSum();
    result[tkp].nodeSum();
      
    mix3[tkp].threadSum();
    mix3[tkp].nodeSum();
      
#ifndef DAIQIAN_COMPATIBILITY_MODE
    result[tkp] *= Float(0.5); //coefficient of 0.5 associated with average over pt2 pion ordering
    mix3[tkp] *= Float(0.5);
#endif
  }

}

CPS_END_NAMESPACE
#endif
