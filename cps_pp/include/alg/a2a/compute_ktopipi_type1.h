#ifndef _COMPUTE_KTOPIPI_TYPE1_H
#define _COMPUTE_KTOPIPI_TYPE1_H

CPS_START_NAMESPACE

//TYPE 1
//Each contraction of this type is made up of different trace combinations of two objects:
//1) \sum_{ \vec x }  \Gamma_1 \prop^L(x_op,x) S_2 \prop^L(x,x_op)    
//2) \sum_{ \vec y, \vec x_K }  \Gamma_2 \prop^L(x_op,y) S_2 \prop^L(y,x_K) \gamma^5 \prop^H(x_K,x_op)
    
//We have to average over the case where object 1 is associated with pi1 (i.e. terminates at the timeslice closer to the kaon)
//and the case where object 1 is associated with pi1 (i.e. terminates at the timeslice further from the kaon)

//Part 2 in terms of v and w is
//\sum_{ \vec y, \vec x_K }  \Gamma_2 vL_i(x_op;y_4) [[ wL_i^dag(y) S_2 vL_j(y;t_K) ]] [[ wL_j^dag(x_K)\gamma^5 vH_k(x_K;top) ]] wH_k^dag(x_op)
//We can use g5-hermiticity to reduce the number of W (stochastic) fields that lie at the operator location:
//\sum_{ \vec y, \vec x_K }  \Gamma_2 vL_i(x_op;y_4) [[ wL_i^dag(y) S_2 vL_j(y;t_K) ]] [[ wL_j^dag(x_K)\gamma^5 \gamma^5 wH_k(x_K) ) ]] vH_k^\dagger(x_op;t_K)\gamma^5

//where [[ ]] indicate meson fields

//Part 1 is
//\sum_{ \vec x }  \Gamma_1 vL_i(x_op; x_4) [[wL_i^dag(x) S_2 vL_j(x;top)]] wL_j^dag(x_op)


//Run inside threaded environment
template<typename mf_Policies>
void ComputeKtoPiPiGparity<mf_Policies>::type1_contract(ResultsContainerType &result, const int t_K, const int t_dis, const int thread_id, const SCFmat part1[2], const SCFmat part2[2]){
#ifndef MEMTEST_MODE
  static const int n_contract = 6; //six type1 diagrams
  static const int con_off = 1; //index of first contraction in set
  for(int pt1_pion=0; pt1_pion<2; pt1_pion++){ //which pion is associated with part 1?
    int pt2_pion = (pt1_pion + 1) % 2;

    for(int mu=0;mu<4;mu++){ //sum over mu here
      for(int gcombidx=0;gcombidx<8;gcombidx++){
	SCFmat G1_pt1 = part1[pt1_pion];
	multGammaLeft(G1_pt1,1,gcombidx,mu);

	CPScolorMatrix<ComplexType> tr_sf_G1_pt1 = G1_pt1.SpinFlavorTrace();

	SCFmat G2_pt2 = part2[pt2_pion];
	multGammaLeft(G2_pt2,2,gcombidx,mu);

	CPScolorMatrix<ComplexType> tr_sf_G2_pt2 = G2_pt2.SpinFlavorTrace();

	SCFmat ctrans_G2_pt2(G2_pt2);
	ctrans_G2_pt2.TransposeColor();

#define C(IDX) result(t_K,t_dis,IDX-con_off,gcombidx,thread_id)	      
	C(1) += G1_pt1.Trace() * G2_pt2.Trace();
	C(2) += Trace( tr_sf_G1_pt1, Transpose(tr_sf_G2_pt2) );
	C(3) += Trace( tr_sf_G1_pt1 , tr_sf_G2_pt2 );
	C(4) += Trace( G1_pt1, G2_pt2 );
	C(5) += Trace( G1_pt1, ctrans_G2_pt2 );
	C(6) += Trace( G1_pt1.ColorTrace() , G2_pt2.ColorTrace() );
#undef C	     
      }
    }
  }
#endif
}

template<typename mf_Policies>
void ComputeKtoPiPiGparity<mf_Policies>::generateRandomOffsets(std::vector<OneFlavorIntegerField*> &random_fields, const std::vector<int> &tsep_k_pi, const int tstep, const int xyzStep){
  int Lt = GJP.Tnodes()*GJP.TnodeSites();
  int size_3d = GJP.VolNodeSites()/GJP.TnodeSites();

  static const NullObject n;
  random_fields.resize(Lt, NULL); //indexed by t_pi_1

  if(xyzStep == 1) return; //no need for random field

  LRG.SetInterval(1.0,0.0);

#ifdef DAIQIAN_EVIL_RANDOM_SITE_OFFSET
  int maxDeltaT = tsep_k_pi[0];
  for(int i=0;i<tsep_k_pi.size();i++) if(tsep_k_pi[i] > maxDeltaT) maxDeltaT = tsep_k_pi[i];

  int t_pi_1_start = GJP.TnodeSites() * GJP.TnodeCoor() + 1;
  int t_pi_1_shift = GJP.TnodeSites() + maxDeltaT - 2; 
  IFloat hi, lo; UniformRandomGenerator::GetInterval(hi,lo);

  for(int shift = 0; shift <= t_pi_1_shift; shift += tstep) {
    int tpi1 = modLt(t_pi_1_start + shift,Lt);
    if(!UniqueID()){ printf("Assigning random field for tpi1=%d\n",tpi1); fflush(stdout); }
    //if(random_fields[tpi1].get() != NULL) continue; //YES I KNOW THIS WILL OVERWRITE SOME ALREADY GENERATED RANDOM NUMBERS

    if(random_fields[tpi1] == NULL){
      random_fields[tpi1] = new OneFlavorIntegerField(n);
      random_fields[tpi1]->zero();
    }

    for(int t_lcl = 0; t_lcl < GJP.TnodeSites(); t_lcl++) {
      int t_op = t_lcl + GJP.TnodeSites() * GJP.TnodeCoor();
      if(t_op >= t_pi_1_start + shift || t_op <= t_pi_1_start + shift - maxDeltaT) continue; //you see why it gives me a headache? These coordinates are not modulo the lattice size

      for(int i = 0; i < size_3d / xyzStep; i++){
	int x = i*xyzStep + GJP.VolNodeSites()/GJP.TnodeSites()*t_lcl; //Generate in same order as Daiqian. He doesn't set the hypercube RNG so it uses whichever one was used last! I know this sucks.
	*(random_fields[tpi1]->site_ptr(x)) = ((int)(LRG.Urand(FOUR_D) * xyzStep)) % xyzStep;
      }
    }
  }
#else
  //for(int t_pi1 = 0; t_pi1 < Lt; t_pi1 += tstep){ //my sensible ordering
  for(int t_pi1_lin = 1; t_pi1_lin <= Lt; t_pi1_lin += tstep){ //Daiqian's weird ordering
    int t_pi1 = modLt(t_pi1_lin,Lt);

    random_fields[t_pi1] = new OneFlavorIntegerField(n);
    random_fields[t_pi1]->zero();
    
    for(int t=0;t<GJP.TnodeSites();t++){
      for(int i = 0; i < size_3d / xyzStep; i++){
	int x = i*xyzStep + GJP.VolNodeSites()/GJP.TnodeSites()*t;
	LRG.AssignGenerator(x);
	*(random_fields[t_pi1].site_ptr(x)) = ((int)(LRG.Urand(FOUR_D) * xyzStep)) % xyzStep;
      }
    }
  }
#endif
}



template<typename mf_Policies>
void ComputeKtoPiPiGparity<mf_Policies>::type1_compute_mfproducts(std::vector<std::vector< A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > > &con_pi1_K,
							       std::vector<std::vector< A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > > &con_pi2_K,
							       const std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> > &mf_pi1,
							       const std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> > &mf_pi2,
							       const std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > &mf_kaon, const MesonFieldMomentumContainer<mf_Policies> &mf_pions,
							       const std::vector<int> &tsep_k_pi, const int tsep_pion, const int Lt, const int ntsep_k_pi){
  //The two meson field are independent of x_op so we can pregenerate them for each y_4, top    
  //Form contraction  con_pi_K(y_4) =   [[ wL_i^dag(y_4) S_2 vL_j(y_4;t_K) ]] [[ wL_j^dag(t_K) wH_k(t_K) ) ]]
  //Compute contraction for each K->pi separation. Try to reuse as there will be some overlap.
  con_pi1_K.resize(Lt); //[tpi][tsep_k_pi]
  con_pi2_K.resize(Lt);
    
  if(!UniqueID()){ printf("Computing con_pi_K\n"); fflush(stdout); }
    
  for(int tpi=0;tpi<Lt;tpi++){
    if(!UniqueID()){ printf("tpi = %d\n",tpi); fflush(stdout); }
    con_pi1_K[tpi].resize(ntsep_k_pi);
    con_pi2_K[tpi].resize(ntsep_k_pi);

    for(int tkpi_idx=0;tkpi_idx<ntsep_k_pi;tkpi_idx++){
      int tk_pi1 = modLt(tpi - tsep_k_pi[tkpi_idx], Lt);
      int tk_pi2 = modLt(tpi - tsep_pion - tsep_k_pi[tkpi_idx], Lt); //on the further timeslice

      mult(con_pi1_K[tpi][tkpi_idx], mf_pi1[tpi], mf_kaon[tk_pi1]); //node and thread distributed
      mult(con_pi2_K[tpi][tkpi_idx], mf_pi2[tpi], mf_kaon[tk_pi2]);
    }
    //NB time coordinate of con_pi1_K, con_pi2_K is the time of the respective pion
  }
  if(!UniqueID()){ printf("Finished computing con_pi_K\n"); fflush(stdout); }
}


template<typename mf_Policies>
void ComputeKtoPiPiGparity<mf_Policies>::type1_mult_vMv_setup(mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorVfftw,A2AvectorW> &mult_vMv_split_part1_pi1,
							   mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorVfftw,A2AvectorW> &mult_vMv_split_part1_pi2,
							   std::vector<mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorWfftw,A2AvectorV> > &mult_vMv_split_part2_pi1,
							   std::vector<mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorWfftw,A2AvectorV> > &mult_vMv_split_part2_pi2,
							   const std::vector<std::vector< A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > > &con_pi1_K,
							   const std::vector<std::vector< A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > > &con_pi2_K,
							   const std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> > &mf_pi1,
							   const std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> > &mf_pi2,							   
							   const A2AvectorV<mf_Policies> & vL, const A2AvectorV<mf_Policies> & vH, 
							   const A2AvectorW<mf_Policies> & wL,
							   const ModeContractionIndices<StandardIndexDilution,TimePackedIndexDilution> &i_ind_vw,
							   const ModeContractionIndices<StandardIndexDilution,FullyPackedIndexDilution> &j_ind_vw,
							   const ModeContractionIndices<TimePackedIndexDilution,StandardIndexDilution> &j_ind_wv,
							   const int top_loc, const int t_pi1, const int t_pi2, 
							   const int Lt, const std::vector<int> &tsep_k_pi, const int ntsep_k_pi, const int t_K_all[], const std::vector<bool> &node_top_used){
  assert(node_top_used[top_loc]);

  //Split the vector-mesonfield outer product into two stages where in the first we reorder the mesonfield to optimize cache hits
  mult_vMv_split_part2_pi1.resize(ntsep_k_pi);
  mult_vMv_split_part2_pi2.resize(ntsep_k_pi);
    
  int top_glb = top_loc  + GJP.TnodeCoor()*GJP.TnodeSites();
  
  mult_vMv_split_part1_pi1.setup( vL, mf_pi1[t_pi1], wL, top_glb, i_ind_vw, j_ind_vw);
  mult_vMv_split_part1_pi2.setup( vL, mf_pi2[t_pi2], wL, top_glb, i_ind_vw, j_ind_vw);
  
#pragma omp parallel for
  for(int tkpi_idx=0;tkpi_idx<ntsep_k_pi;tkpi_idx++){
    int t_dis = modLt(top_glb - t_K_all[tkpi_idx], Lt);
    if(t_dis >= tsep_k_pi[tkpi_idx] || t_dis == 0) continue; //save time by skipping setup on those we are not going to use

    mult_vMv_split_part2_pi1[tkpi_idx].setup( vL, con_pi1_K[t_pi1][tkpi_idx], vH, top_glb, i_ind_vw, j_ind_wv);
    mult_vMv_split_part2_pi2[tkpi_idx].setup( vL, con_pi2_K[t_pi2][tkpi_idx], vH, top_glb, i_ind_vw, j_ind_wv);
  }
}





template<typename mf_Policies>
void ComputeKtoPiPiGparity<mf_Policies>::type1_precompute_part1_part2(SCFmatVector &mult_vMv_contracted_part1_pi1,
								   SCFmatVector &mult_vMv_contracted_part1_pi2,
								   std::vector<SCFmatVector > &mult_vMv_contracted_part2_pi1,
								   std::vector<SCFmatVector > &mult_vMv_contracted_part2_pi2,
								   mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorVfftw,A2AvectorW> &mult_vMv_split_part1_pi1,
								   mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorVfftw,A2AvectorW> &mult_vMv_split_part1_pi2,
								   std::vector<mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorWfftw,A2AvectorV> > &mult_vMv_split_part2_pi1,
								   std::vector<mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorWfftw,A2AvectorV> > &mult_vMv_split_part2_pi2,
								   const int top_loc, const int Lt, const std::vector<int> &tsep_k_pi, const int ntsep_k_pi, const int t_K_all[], const std::vector<bool> &node_top_used){


  assert(node_top_used[top_loc]);

  //Contract on all 3d sites on this node with fixed operator time coord top_glb into a canonically ordered output vector
  mult_vMv_contracted_part2_pi1.resize(ntsep_k_pi);
  mult_vMv_contracted_part2_pi2.resize(ntsep_k_pi);

  mult_vMv_split_part1_pi1.contract(mult_vMv_contracted_part1_pi1,false,true);
  mult_vMv_split_part1_pi1.free_mem();
  
  mult_vMv_split_part1_pi2.contract(mult_vMv_contracted_part1_pi2,false,true);
  mult_vMv_split_part1_pi2.free_mem();

  int top_glb = top_loc  + GJP.TnodeCoor()*GJP.TnodeSites();

  for(int tkpi_idx=0;tkpi_idx<ntsep_k_pi;tkpi_idx++){
    int t_dis = modLt(top_glb - t_K_all[tkpi_idx], Lt);
    if(t_dis >= tsep_k_pi[tkpi_idx] || t_dis == 0) continue; //save time by skipping setup on those we are not going to use

    mult_vMv_split_part2_pi1[tkpi_idx].contract(mult_vMv_contracted_part2_pi1[tkpi_idx],false,true);
    mult_vMv_split_part2_pi1[tkpi_idx].free_mem();
    
    mult_vMv_split_part2_pi2[tkpi_idx].contract(mult_vMv_contracted_part2_pi2[tkpi_idx],false,true);
    mult_vMv_split_part2_pi2[tkpi_idx].free_mem();
  }
}



//xyzStep is the 3d spatial sampling frequency. If set to anything other than one it will compute the diagram
//for a random site within the block (site, site+xyzStep) in canonical ordering. Daiqian's original implementation is machine-size dependent, but for repro I had to add an option to do it his way 

//This version overlaps computation for multiple K->pi separations. Result should be an array of ResultsContainerType the same size as the vector 'tsep_k_pi'
template<typename mf_Policies>
void ComputeKtoPiPiGparity<mf_Policies>::type1(ResultsContainerType result[],
					    const std::vector<int> &tsep_k_pi, const int tsep_pion, const int tstep, const int xyzStep, const ThreeMomentum &p_pi_1, 
					    const std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > &mf_kaon, MesonFieldMomentumContainer<mf_Policies> &mf_pions,
					    const A2AvectorV<mf_Policies> & vL, const A2AvectorV<mf_Policies> & vH, 
					    const A2AvectorW<mf_Policies> & wL, const A2AvectorW<mf_Policies> & wH){

  //Precompute mode mappings
  ModeContractionIndices<StandardIndexDilution,TimePackedIndexDilution> i_ind_vw(vL);
  ModeContractionIndices<StandardIndexDilution,FullyPackedIndexDilution> j_ind_vw(wL);
  ModeContractionIndices<TimePackedIndexDilution,StandardIndexDilution> j_ind_wv(vH);

  const int Lt = GJP.Tnodes()*GJP.TnodeSites();
  const int ntsep_k_pi = tsep_k_pi.size();
    
  const ThreeMomentum p_pi_2 = -p_pi_1;

  static const int n_contract = 6; //six type1 diagrams
  static const int con_off = 1; //index of first contraction in set
  int nthread = omp_get_max_threads();

  for(int i=0;i<ntsep_k_pi;i++)
    result[i].resize(n_contract,nthread); //it will be thread-reduced before this method ends. Resize also zeroes 'result'
    
  const int size_3d = vL.getMode(0).nodeSites(0)*vL.getMode(0).nodeSites(1)*vL.getMode(0).nodeSites(2);

  std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> > &mf_pi1 = mf_pions.get(p_pi_1); //*mf_pi1_ptr;
  std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> > &mf_pi2 = mf_pions.get(p_pi_2); //*mf_pi2_ptr;
#ifdef NODE_DISTRIBUTE_MESONFIELDS
  if(!UniqueID()) printf("Memory prior to fetching meson fields type1 K->pipi:\n");
  printMem();
  nodeGetMany(2,&mf_pi1,&mf_pi2);
  if(!UniqueID()) printf("Memory after fetching meson fields type1 K->pipi:\n");
  printMem();
#endif

  //nodeGetPionMf(mf_pi1,mf_pi2);
  
  std::vector<std::vector< A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > > con_pi1_K(Lt); //[tpi][tsep_k_pi]
  std::vector<std::vector< A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorWfftw> > > con_pi2_K(Lt);
    
  type1_compute_mfproducts(con_pi1_K,con_pi2_K,mf_pi1,mf_pi2,mf_kaon,mf_pions,tsep_k_pi,tsep_pion,Lt,ntsep_k_pi);

  if(!UniqueID()) printf("Memory after computing mfproducts type1 K->pipi:\n");
  printMem();

  //Generate the random offsets in the same order as Daiqian did (an order which gives me a headache to think about)
  std::vector<OneFlavorIntegerField*> random_fields; //indexed by t_pi_1
  generateRandomOffsets(random_fields,tsep_k_pi,tstep,xyzStep);

  if(!UniqueID()) printf("Memory after generating random offsets type1 K->pipi:\n");
  printMem();

  //for(int t_pi1=0;t_pi1<GJP.TnodeSites();t_pi1++) if(!UniqueID()){ printf("Random field ptr t_pi1=%d -> %p\n",t_pi1, random_fields[t_pi1]); fflush(stdout); } 


  //for(int t_pi1 = 0; t_pi1 < Lt; t_pi1 += tstep){ //my sensible ordering
  for(int t_pi1_lin = 1; t_pi1_lin <= Lt; t_pi1_lin += tstep){ //Daiqian's weird ordering
    int t_pi1 = modLt(t_pi1_lin,Lt);

    int t_pi2 = modLt(t_pi1 + tsep_pion, Lt);

    //Using the pion timeslices, get tK for each separation
    int t_K_all[ntsep_k_pi]; 
    for(int tkpi_idx=0;tkpi_idx<ntsep_k_pi;tkpi_idx++) 
      t_K_all[tkpi_idx] = modLt(t_pi1 - tsep_k_pi[tkpi_idx], Lt);

    //Determine what node timeslices are actually needed
    std::vector<bool> node_top_used;
    getUsedTimeslices(node_top_used,tsep_k_pi,t_pi1);

    for(int top_loc = 0; top_loc < GJP.TnodeSites(); top_loc++){
      if(!node_top_used[top_loc]) continue;
      const int top_glb = top_loc  + GJP.TnodeCoor()*GJP.TnodeSites();

      //Split the vector-mesonfield outer product into two stages where in the first we reorder the mesonfield to optimize cache hits
#ifndef DISABLE_TYPE1_SPLIT_VMV
      mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorVfftw,A2AvectorW> mult_vMv_split_part1_pi1;
      mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorVfftw,A2AvectorW> mult_vMv_split_part1_pi2;
      std::vector<mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorWfftw,A2AvectorV> > mult_vMv_split_part2_pi1; //[ntsep_k_pi];
      std::vector<mult_vMv_split<mf_Policies,A2AvectorV,A2AvectorWfftw,A2AvectorWfftw,A2AvectorV> > mult_vMv_split_part2_pi2; //[ntsep_k_pi];
      
      type1_mult_vMv_setup(mult_vMv_split_part1_pi1,mult_vMv_split_part1_pi2,mult_vMv_split_part2_pi1,mult_vMv_split_part2_pi2,
			   con_pi1_K,con_pi2_K,mf_pi1,mf_pi2,vL,vH,wL,i_ind_vw,j_ind_vw,j_ind_wv,top_loc, t_pi1,t_pi2, Lt, tsep_k_pi, ntsep_k_pi,t_K_all,node_top_used);

# ifndef DISABLE_TYPE1_PRECOMPUTE
      //Contract on all 3d sites on this node with fixed operator time coord top_glb into a canonically ordered output vector
      SCFmatVector mult_vMv_contracted_part1_pi1; //[x3d];
      SCFmatVector mult_vMv_contracted_part1_pi2; //[x3d];
      
      std::vector<SCFmatVector> mult_vMv_contracted_part2_pi1; //[ntsep_k_pi][x3d];
      std::vector<SCFmatVector> mult_vMv_contracted_part2_pi2; //[ntsep_k_pi][x3d];

    
      type1_precompute_part1_part2(mult_vMv_contracted_part1_pi1,mult_vMv_contracted_part1_pi2,mult_vMv_contracted_part2_pi1,mult_vMv_contracted_part2_pi2,mult_vMv_split_part1_pi1,mult_vMv_split_part1_pi2,
				   mult_vMv_split_part2_pi1,mult_vMv_split_part2_pi2, top_loc, Lt,tsep_k_pi,ntsep_k_pi, t_K_all,node_top_used);
# endif
#endif

      if(is_grid_vector_complex<typename mf_Policies::ComplexType>::value && xyzStep > 1) ERR.General("ComputeKtoPiPiGparity<mf_Policies>","type1","Random offset not implemented for Grid-vectorized fields\n");
      if(xyzStep > 1 && random_fields[t_pi1] == NULL) ERR.General("ComputeKtoPiPiGparity<mf_Policies>","type1","Random field not initialized for t_pi1=%d (got %p) on node %d\n",t_pi1,random_fields[t_pi1],UniqueID());
      OneFlavorIntegerField const* randoff = random_fields[t_pi1];

      //Now loop over Q_i insertion location. Each node naturally has its own sublattice to work on. Thread over sites in usual way
#pragma omp parallel for
      for(int xop3d_loc_base = 0; xop3d_loc_base < size_3d; xop3d_loc_base+=xyzStep){
	int thread_id = omp_get_thread_num();

	int dx = (xyzStep > 1 ? *(randoff->site_ptr(xop3d_loc_base+size_3d*top_loc)) : 0);
	int xop3d_loc = xop3d_loc_base + dx;
	assert(xop3d_loc < GJP.VolNodeSites()/GJP.TnodeSites());

	//if(!UniqueID()){ printf("top_loc=%d, xop3d_loc_base=%d, dx=%d, xop3d_loc=%d\n",top_loc,xop3d_loc_base,dx,xop3d_loc); fflush(stdout); }

	//Construct part 1:
	//\Gamma_1 vL_i(x_op; x_4) [[\sum_{\vec x} wL_i^dag(x) S_2 vL_j(x;top)]] wL_j^dag(x_op)
	SCFmat part1[2]; //part1 goes from insertion to pi1, pi2 (x_4 = t_pi1, t_pi2)
	  
#if defined(DISABLE_TYPE1_SPLIT_VMV)
	mult(part1[0], vL, mf_pi1[t_pi1], wL, xop3d_loc, top_loc, false, true);
	mult(part1[1], vL, mf_pi2[t_pi2], wL, xop3d_loc, top_loc, false, true);
#elif defined(DISABLE_TYPE1_PRECOMPUTE)
	mult_vMv_split_part1_pi1.contract(part1[0],xop3d_loc,false, true);
	mult_vMv_split_part1_pi2.contract(part1[1],xop3d_loc,false, true);
#else
	part1[0] = mult_vMv_contracted_part1_pi1[xop3d_loc];
	part1[1] = mult_vMv_contracted_part1_pi2[xop3d_loc];
#endif

	for(int tkidx =0; tkidx< ntsep_k_pi; tkidx++){
	  int t_K = t_K_all[tkidx];
	  int t_dis = modLt(top_glb - t_K, Lt); //distance between kaon and operator is the output time coordinate
	  if(t_dis >= tsep_k_pi[tkidx] || t_dis == 0) continue; //don't bother computing operator insertion locations outside of the region between the kaon and first pion or on top of either operator
	  	
	  //Construct part 2:
	  //\Gamma_2 vL_i(x_op;y_4) [[ wL_i^dag(y) S_2 vL_j(y;t_K) ]] [[ wL_j^dag(x_K)\gamma^5 \gamma^5 wH_k(x_K) ) ]] vH_k^\dagger(x_op;t_K)\gamma^5
	  SCFmat part2[2]; //part2 has pi1, pi2 at y_4 = t_pi1, t_pi2

#if defined(DISABLE_TYPE1_SPLIT_VMV)
	  mult(part2[0], vL, con_pi1_K[t_pi1][tkidx], vH, xop3d_loc, top_loc, false, true);
	  mult(part2[1], vL, con_pi2_K[t_pi2][tkidx], vH, xop3d_loc, top_loc, false, true);
#elif defined(DISABLE_TYPE1_PRECOMPUTE)
	  mult_vMv_split_part2_pi1[tkidx].contract(part2[0],xop3d_loc,false, true);
	  mult_vMv_split_part2_pi2[tkidx].contract(part2[1],xop3d_loc,false, true);
#else
	  part2[0] = mult_vMv_contracted_part2_pi1[tkidx][xop3d_loc];
	  part2[1] = mult_vMv_contracted_part2_pi2[tkidx][xop3d_loc];
#endif

	  part2[0].gr(-5); //right multiply by g5
	  part2[1].gr(-5);

	  type1_contract(result[tkidx],t_K,t_dis,thread_id,part1,part2);
	}//end of loop over k->pi seps

      }//xop3d loop
    }//top_loc loop
  }//tpi loop

  if(!UniqueID()) printf("Memory before finishing up type1 K->pipi:\n");
  printMem();


  if(!UniqueID()){ printf("Type 1 finishing up results\n"); fflush(stdout); }
  for(int tkpi_idx =0; tkpi_idx< ntsep_k_pi; tkpi_idx++){
    result[tkpi_idx].threadSum();
    result[tkpi_idx].nodeSum();
    result[tkpi_idx] *= Float(0.5); //coefficient of 0.5 associated with average over pt1_pion loop
    result[tkpi_idx] *= Float(xyzStep); //correct for spatial blocking
  }
  if(!UniqueID()) printf("Memory after finishing up type1 K->pipi:\n");
  printMem();

  for(int i=0;i<random_fields.size();i++) if(random_fields[i] != NULL) delete random_fields[i];

  if(!UniqueID()) printf("Memory after deleting random fields type1 K->pipi:\n");
  printMem();
  

#ifdef NODE_DISTRIBUTE_MESONFIELDS
  nodeDistributeMany(2,&mf_pi1,&mf_pi2);
  if(!UniqueID()) printf("Memory after redistributing meson fields type1 K->pipi:\n");
  printMem();
#endif

  //nodeDistributePionMf(mf_pi1,mf_pi2);
}



CPS_END_NAMESPACE

#endif
