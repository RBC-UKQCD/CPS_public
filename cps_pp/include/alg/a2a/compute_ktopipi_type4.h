#ifndef _COMPUTE_KTOPIPI_TYPE4_H
#define _COMPUTE_KTOPIPI_TYPE4_H

CPS_START_NAMESPACE

//Each contraction of this type is made up of different trace combinations of two objects (below for simplicity we ignore the fact that the two vectors in the meson fields are allowed to vary in position relative to each other):
//1) \prop^L(x_op,x_K) \gamma^5 \prop^H(x_K,x_op)
//   we use g5-hermiticity on the strange prop
//  \prop^L(x_op,x_K)  [ \prop^H(x_op,x_K) ]^dag \gamma^5
//= vL(x_op) [[ wL^dag(x_K) wH(x_K) ]] vH^dag(x_op) \gamma_5
    
//2) \prop^L(x_op,x_op)   OR   \prop^H(x_op,x_op)
//   = vL(x_op) wL^dag(x_op)   or  vH(x_op) wH^dag(x_op)

//This disconnected 'bubble' is computed during the pipi calculation and therefore there is no need to recompute
template<typename mf_Policies>
void ComputeKtoPiPiGparity<mf_Policies>::type4_contract(ResultsContainerType &result, const int t_K, const int t_dis, const int thread_id, 
						     const SCFmat &part1, const SCFmat &part2_L, const SCFmat &part2_H){
#ifndef MEMTEST_MODE
  static const int con_off = 23; //index of first contraction in set

  for(int mu=0;mu<4;mu++){ //sum over mu here
    for(int gcombidx=0;gcombidx<8;gcombidx++){
      const SCFmat &G1 = Gamma1<ComplexType>(gcombidx,mu);
      const SCFmat &G2 = Gamma2<ComplexType>(gcombidx,mu);

      //SCFmat G1_pt1 = G1*part1;
      SCFmat G1_pt1 = part1;
      multGammaLeft(G1_pt1,1,gcombidx,mu);

      CPScolorMatrix<ComplexType> tr_sf_G1_pt1 = G1_pt1.SpinFlavorTrace();

      //SCFmat G2_pt2_L = G2*part2_L;
      SCFmat G2_pt2_L = part2_L;
      multGammaLeft(G2_pt2_L,2,gcombidx,mu);

      CPScolorMatrix<ComplexType> tr_sf_G2_pt2_L = G2_pt2_L.SpinFlavorTrace();

      //SCFmat G2_pt2_H = G2*part2_H;
      SCFmat G2_pt2_H = part2_H;
      multGammaLeft(G2_pt2_H,2,gcombidx,mu);

      CPScolorMatrix<ComplexType> tr_sf_G2_pt2_H = G2_pt2_H.SpinFlavorTrace();
	    
      SCFmat ctrans_G2_pt2_L(G2_pt2_L); //speedup by transposing part 1
      ctrans_G2_pt2_L.TransposeColor();
	
      CPSspinMatrix<CPSflavorMatrix<ComplexType> > tr_c_G1_pt1 = G1_pt1.ColorTrace();
	
#define C(IDX) result(t_K,t_dis,IDX-con_off,gcombidx,thread_id)	      

      //First 6 have a light-quark loop
      C(23) += G1_pt1.Trace() * G2_pt2_L.Trace();

      C(24) += Trace( tr_sf_G1_pt1 , Transpose(tr_sf_G2_pt2_L) );
      C(25) += Trace( tr_sf_G1_pt1 , tr_sf_G2_pt2_L );
      C(26) += Trace( G1_pt1 , G2_pt2_L );
      C(27) += Trace( G1_pt1 , ctrans_G2_pt2_L );
      C(28) += Trace( tr_c_G1_pt1 , G2_pt2_L.ColorTrace() );
	      
      //Second 4 have strange loop
      C(29) += G1_pt1.Trace() * G2_pt2_H.Trace();	      
      C(30) += Trace( tr_sf_G1_pt1 , tr_sf_G2_pt2_H );
      C(31) += Trace( G1_pt1 , G2_pt2_H );
      C(32) += Trace( tr_c_G1_pt1 , G2_pt2_H.ColorTrace() );

#undef C	     	    
    }
  }
#endif
}

template<typename mf_Policies>
void ComputeKtoPiPiGparity<mf_Policies>::type4_mult_vMv_setup(std::vector<mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorWfftw,A2AvectorV> > &mult_vMv_split_part1,
							   const std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > &mf_kaon,
							   const A2AvectorV<mf_Policies> & vL, const A2AvectorV<mf_Policies> & vH,
							   const int top_loc, const int tstep, const int Lt){
  mult_vMv_split_part1.resize(Lt/tstep); //[tKidx]
  int top_glb = top_loc  + GJP.TnodeCoor()*GJP.TnodeSites();

  for(int tkidx=0; tkidx < Lt/tstep; tkidx++)
    mult_vMv_split_part1[tkidx].setup(vL,mf_kaon[tkidx*tstep],vH,top_glb);
}

template<typename mf_Policies>
void ComputeKtoPiPiGparity<mf_Policies>::type4_precompute_part1(std::vector<SCFmatVector> &mult_vMv_contracted_part1,
							     std::vector<mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorWfftw,A2AvectorV> > &mult_vMv_split_part1,
							     const int top_loc, const int tstep, const int Lt){

  mult_vMv_contracted_part1.resize(Lt/tstep); //[tKidx]

  for(int tkidx=0; tkidx < Lt/tstep; tkidx++){
    mult_vMv_split_part1[tkidx].contract(mult_vMv_contracted_part1[tkidx],false,true);
    mult_vMv_split_part1[tkidx].free_mem();
  }
}





template<typename mf_Policies>
void ComputeKtoPiPiGparity<mf_Policies>::type4(ResultsContainerType &result, MixDiagResultsContainerType &mix4,
					    const int &tstep,
					    const std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > &mf_kaon,
					    const A2AvectorV<mf_Policies> & vL, const A2AvectorV<mf_Policies> & vH, 
					    const A2AvectorW<mf_Policies> & wL, const A2AvectorW<mf_Policies> & wH){
  
  
  SCFmat mix4_Gamma[2];
  mix4_Gamma[0].unit().pr(F0).gr(-5);
  mix4_Gamma[1].unit().pr(F1).gr(-5).timesMinusOne();
  
  //CK: the loop term could be re-used from type3
  int Lt = GJP.Tnodes()*GJP.TnodeSites();
  assert(Lt % tstep == 0);
  
  static const int n_contract = 10; //ten type4 diagrams
  static const int con_off = 23; //index of first contraction in set
  int nthread = omp_get_max_threads();

  result.resize(n_contract,nthread); //it will be thread-reduced before this method ends
  mix4.resize(nthread);

  int size_3d = vL.getMode(0).nodeSites(0)*vL.getMode(0).nodeSites(1)*vL.getMode(0).nodeSites(2);

  for(int top_loc = 0; top_loc < GJP.TnodeSites(); top_loc++){
    const int top_glb = top_loc  + GJP.TnodeCoor()*GJP.TnodeSites();

#ifndef DISABLE_TYPE4_SPLIT_VMV
    std::vector<mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorWfftw,A2AvectorV> > mult_vMv_split_part1; //[tkidx in Lt/tstep]
    type4_mult_vMv_setup(mult_vMv_split_part1,mf_kaon,vL,vH,top_loc,tstep,Lt);

# ifndef DISABLE_TYPE4_PRECOMPUTE
    std::vector<SCFmatVector > mult_vMv_contracted_part1; //[tkidx in Lt/tstep][x3d]
    type4_precompute_part1(mult_vMv_contracted_part1,mult_vMv_split_part1,top_loc,tstep,Lt);
# endif
#endif

    //Now loop over Q_i insertion location. Each node naturally has its own sublattice to work on. Thread over sites in usual way
#pragma omp parallel for
    for(int xop3d_loc = 0; xop3d_loc < size_3d; xop3d_loc++){
      int thread_id = omp_get_thread_num();

      for(int t_K = 0; t_K < Lt; t_K += tstep){ //global times
	int t_dis = modLt(top_glb - t_K, Lt); //distance between kaon and operator is the output time coordinate

	//Construct part 1:
	// = vL(x_op) [[ wL^dag(x_K) wH(x_K) ]] vH^dag(x_op) \gamma_5
	SCFmat part1;

#if defined(DISABLE_TYPE4_SPLIT_VMV)
	mult(part1, vL, mf_kaon[t_K], vH, xop3d_loc, top_loc, false, true);
#elif defined(DISABLE_TYPE4_PRECOMPUTE)
	mult_vMv_split_part1[t_K/tstep].contract(part1,xop3d_loc,false,true);
#else
	part1 = mult_vMv_contracted_part1[t_K/tstep][xop3d_loc];
#endif

	part1.gr(-5);
		
	//Construct part 2:
	//vL(x_op) wL^dag(x_op)   or  vH(x_op) wH^dag(x_op)  (CK: should re-use these from type-3)
	SCFmat part2_L, part2_H;
	mult(part2_L, vL, wL, xop3d_loc, top_loc, false, true);
	mult(part2_H, vH, wH, xop3d_loc, top_loc, false, true);

	type4_contract(result,t_K,t_dis,thread_id,part1,part2_L,part2_H);

	//Compute mix4 diagram
	//These are identical to the type4 diagrams but without the quark loop, and with the vertex replaced with a pseudoscalar vertex
	for(int mix4_gidx=0; mix4_gidx<2; mix4_gidx++){
#ifndef MEMTEST_MODE
#define M mix4(t_K,t_dis,mix4_gidx,thread_id)
	  //M += ( part1 * mix4_Gamma[mix4_gidx] ).Trace();
	  M += Trace( part1 , mix4_Gamma[mix4_gidx] );
#undef M
#endif
	}

      }//t_K loop

    }//xop3d loop
  }//top_loc loop

  result.threadSum();
  result.nodeSum();

  mix4.threadSum();
  mix4.nodeSum();
}

CPS_END_NAMESPACE

#endif
