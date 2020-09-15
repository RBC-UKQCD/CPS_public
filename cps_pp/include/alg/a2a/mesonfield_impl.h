#ifndef MESON_FIELD_IMPL
#define MESON_FIELD_IMPL

#include<alg/a2a/mult_impl.h>
#include<alg/a2a/mult_vMv_split.h>
#include<alg/a2a/mult_vMv_split_grid.h>
#include<alg/a2a/mult_vMv_impl.h>
#include<alg/a2a/mult_vv_impl.h>
#include<alg/a2a/mesonfield_io.h>
#include<alg/a2a/mesonfield_compute_impl.h>

template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::plus_equals(const A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR> &with, const bool parallel){
  if(nmodes_l != with.nmodes_l || nmodes_r != with.nmodes_r || 
     !lindexdilution.paramsEqual(with.lindexdilution) || !rindexdilution.paramsEqual(with.rindexdilution) ){
    ERR.General("A2AmesonField","plus_equals(..)","Second meson field must have the same underlying parameters\n");
  }
  if(parallel){
#pragma omp_parallel for
    for(int i=0;i<fsize;i++) mf[i] += with.mf[i];
  }else{
    for(int i=0;i<fsize;i++) mf[i] += with.mf[i];
  }
}

template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::times_equals(const ScalarComplexType f,const bool parallel){
  if(parallel){
#pragma omp_parallel for
    for(int i=0;i<fsize;i++) mf[i] *= f;			       
  }else{
    for(int i=0;i<fsize;i++) mf[i] *= f;
  }
}


template<typename T>
inline std::complex<T> complexAvg(const std::complex<T>&a, const std::complex<T> &b){
  return (a+b)/T(2.0);
}

//Replace this meson field with the average of this and a second field, 'with'
template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::average(const A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR> &with, const bool parallel){
  if(nmodes_l != with.nmodes_l || nmodes_r != with.nmodes_r || 
     !lindexdilution.paramsEqual(with.lindexdilution) || !rindexdilution.paramsEqual(with.rindexdilution) ){
    ERR.General("A2AmesonField","average(..)","Second meson field must have the same underlying parameters\n");
  }
  if(parallel){
#pragma omp_parallel for
    for(int i=0;i<fsize;i++) mf[i] = complexAvg(mf[i],with.mf[i]);//(mf[i] + with.mf[i])/2.0;
  }else{
    for(int i=0;i<fsize;i++) mf[i] = complexAvg(mf[i],with.mf[i]);//(mf[i] + with.mf[i])/2.0;
  }
}

//Reorder the rows so that all the elements in idx_map are sequential. Indices not in map are ignored. Use at your own risk
template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::rowReorder(A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR> &into, const int idx_map[], int map_size, bool parallel) const{
  into.setup(lindexdilution, rindexdilution, tl, tr);

#define DOIT \
    int irow = idx_map[i]; \
    for(int j=0;j<nmodes_r;j++) \
      into(i,j) = (*this)(irow,j);

  if(parallel){
#pragma omp parallel for
    for(int i=0;i<map_size;i++){
      DOIT;
    }
  }else{
    for(int i=0;i<map_size;i++){
      DOIT;
    }
  }
#undef DOIT

}
template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::colReorder(A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR> &into, const int idx_map[], int map_size, bool parallel) const{
  into.setup(lindexdilution, rindexdilution, tl, tr);

#define DOIT \
    for(int j=0;j<map_size;j++){ \
      int jcol = idx_map[j]; \
      into(i,j) = (*this)(i,jcol); \
    }

  if(parallel){
#pragma omp parallel for
    for(int i=0;i<nmodes_l;i++){
      DOIT;
    }
  }else{
    for(int i=0;i<nmodes_l;i++){
      DOIT;
    } 
  }
}



//Do a column reorder but where we pack the row indices to exclude those not used (as indicated by input bool array)
//Output as a GSL matrix
template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
typename gsl_wrapper<typename mf_Policies::ScalarComplexType::value_type>::matrix_complex * A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::GSLpackedColReorder(const int idx_map[], int map_size, bool rowidx_used[], typename gsl_wrapper<typename ScalarComplexType::value_type>::matrix_complex *reuse ) const{
  typedef gsl_wrapper<typename ScalarComplexType::value_type> gw;
  assert(sizeof(typename gw::complex) == sizeof(ScalarComplexType));
  int rows = nmodes_l;
  int cols = nmodes_r;

  int nrows_used = 0;
  for(int i_full=0;i_full<rows;i_full++) if(rowidx_used[i_full]) nrows_used++;

  typename gw::matrix_complex *M_packed;
  if(reuse!=NULL){
    M_packed = reuse;
    M_packed->size1 = nrows_used;
    M_packed->size2 = M_packed->tda = map_size;
  }else M_packed = gw::matrix_complex_alloc(nrows_used,map_size);

  //Look for contiguous blocks in the idx_map we can take advantage of
  std::vector<std::pair<int,int> > blocks;
  find_contiguous_blocks(blocks,idx_map,map_size);

  int i_packed = 0;

  for(int i_full=0;i_full<rows;i_full++){
    if(rowidx_used[i_full]){
      ScalarComplexType const* mf_row_base = mf + nmodes_r*i_full; //meson field are row major so columns are contiguous
      typename gw::complex* row_base = gw::matrix_complex_ptr(M_packed,i_packed,0); //GSL matrix are also row major
      for(int b=0;b<blocks.size();b++){
	ScalarComplexType const* block_ptr = mf_row_base + idx_map[blocks[b].first];
	memcpy((void*)row_base,(void*)block_ptr,blocks[b].second*sizeof(ScalarComplexType));
	row_base += blocks[b].second;
      }
      i_packed++;
    }
  }

  return M_packed;
}

#ifdef USE_GRID
  //Do a column reorder but where we pack the row indices to exclude those not used (as indicated by input bool array)
  //Output to a linearized matrix of Grid SIMD vectors where we have splatted the scalar onto all SIMD lanes
  //Does not set the size of the output vector, allowing reuse of a previously allocated vector providing it's large enough
template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::splatPackedColReorder(Grid::Vector<typename mf_Policies::ComplexType> &into, const int idx_map[], int map_size, bool rowidx_used[], bool do_resize) const{
  typedef typename mf_Policies::ComplexType VectorComplexType;
  int full_rows = nmodes_l;
  int full_cols = nmodes_r;

  int nrows_used = 0;
  for(int i_full=0;i_full<full_rows;i_full++) if(rowidx_used[i_full]) nrows_used++;

  if(do_resize) into.resize(nrows_used*map_size);

  //Look for contiguous blocks in the idx_map we can take advantage of
  std::vector<std::pair<int,int> > blocks;
  find_contiguous_blocks(blocks,idx_map,map_size);

  int i_packed = 0;

  for(int i_full=0;i_full<full_rows;i_full++){
    if(rowidx_used[i_full]){
      ScalarComplexType const* mf_row_base = mf + nmodes_r*i_full;
      VectorComplexType* row_base = &into[map_size*i_packed];

      for(int b=0;b<blocks.size();b++){
	ScalarComplexType const* block_ptr = mf_row_base + idx_map[blocks[b].first];
	for(int bb=0;bb<blocks[b].second;bb++)
	  vsplat(row_base[bb],block_ptr[bb]);
	row_base += blocks[b].second;
      }
      i_packed++;
    }
  }
}
template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::scalarPackedColReorder(Grid::Vector<typename mf_Policies::ScalarComplexType> &into, const int idx_map[], int map_size, bool rowidx_used[], bool do_resize) const{
  int full_rows = nmodes_l;
  int full_cols = nmodes_r;

  int nrows_used = 0;
  for(int i_full=0;i_full<full_rows;i_full++) if(rowidx_used[i_full]) nrows_used++;

  if(do_resize) into.resize(nrows_used*map_size);

  //Look for contiguous blocks in the idx_map we can take advantage of
  std::vector<std::pair<int,int> > blocks;
  find_contiguous_blocks(blocks,idx_map,map_size);

  int i_packed = 0;

  for(int i_full=0;i_full<full_rows;i_full++){
    if(rowidx_used[i_full]){
      ScalarComplexType const* mf_row_base = mf + nmodes_r*i_full;
      ScalarComplexType* row_base = &into[map_size*i_packed];

      for(int b=0;b<blocks.size();b++){
	ScalarComplexType const* block_ptr = mf_row_base + idx_map[blocks[b].first];
	for(int bb=0;bb<blocks[b].second;bb++)
	  row_base[bb] = block_ptr[bb];
	row_base += blocks[b].second;
      }
      i_packed++;
    }
  }
}
#endif




template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::transpose(A2AmesonField<mf_Policies,A2AfieldR,A2AfieldL> &into) const{
  assert( (void*)this != (void*)&into );
  into.setup(rindexdilution, lindexdilution, tr, tl);
#pragma omp parallel for
  for(int i=0;i<nmodes_l;i++)
    for(int j=0;j<nmodes_r;j++)
      into(j,i) = (*this)(i,j);
}

//Compute   l^ij(t1,t2) r^ji(t3,t4)
//Threaded but *node local*
template<typename mf_Policies, 
	 template <typename> class lA2AfieldL,  template <typename> class lA2AfieldR,
	 template <typename> class rA2AfieldL,  template <typename> class rA2AfieldR
	 >
typename mf_Policies::ScalarComplexType trace(const A2AmesonField<mf_Policies,lA2AfieldL,lA2AfieldR> &l, const A2AmesonField<mf_Policies,rA2AfieldL,rA2AfieldR> &r){
  //Check the indices match
  if(! l.getRowParams().paramsEqual( r.getColParams() ) || ! l.getColParams().paramsEqual( r.getRowParams() ) )
    ERR.General("","trace(const A2AmesonField<mf_Policies,lA2AfieldL,lA2AfieldR> &, const A2AmesonField<mf_Policies,rA2AfieldL,rA2AfieldR> &)","Illegal matrix product: underlying vector parameters must match\n");
  typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
  ScalarComplexType into(0,0);

  typedef typename A2AmesonField<mf_Policies,lA2AfieldL,lA2AfieldR>::LeftDilutionType DilType0;
  typedef typename A2AmesonField<mf_Policies,lA2AfieldL,lA2AfieldR>::RightDilutionType DilType1;
  typedef typename A2AmesonField<mf_Policies,rA2AfieldL,rA2AfieldR>::LeftDilutionType DilType2;
  typedef typename A2AmesonField<mf_Policies,rA2AfieldL,rA2AfieldR>::RightDilutionType DilType3;

  ModeContractionIndices<DilType0,DilType3> i_ind(l.getRowParams());
  ModeContractionIndices<DilType1,DilType2> j_ind(r.getRowParams());

  int times[4] = { l.getRowTimeslice(), l.getColTimeslice(), r.getRowTimeslice(), r.getColTimeslice() };

  //W * W is only non-zero when the timeslice upon which we evaluate them are equal
  const int n_threads = omp_get_max_threads();
  std::vector<ScalarComplexType> ret_vec(n_threads,(0.,0.));
    
  modeIndexSet lip; lip.time = times[0];
  modeIndexSet rip; rip.time = times[3];

  modeIndexSet ljp; ljp.time = times[1];
  modeIndexSet rjp; rjp.time = times[2];

  const int &ni = i_ind.getNindices(lip,rip); //how many indices to loop over
  const int &nj = j_ind.getNindices(ljp,rjp);

#ifndef MEMTEST_MODE
  
#pragma omp parallel for
  for(int i = 0; i < ni; i++){
    int id = omp_get_thread_num();
    int li = i_ind.getLeftIndex(i,lip,rip);
    int ri = i_ind.getRightIndex(i,lip,rip);

    for(int j = 0; j < nj; j++){
      int lj = j_ind.getLeftIndex(j,ljp,rjp);
      int rj = j_ind.getRightIndex(j,ljp,rjp);
      
      ret_vec[id] += l(li,lj) *  r(rj,ri);
    }
  }
  for(int i=0;i<n_threads;i++) into += ret_vec[i];
  
#endif
  
  return into;
}




//Compute   l^ij(t1,t2) r^ji(t3,t4) for all t1, t4  and place into matrix element t1,t4
//This is both threaded and distributed over nodes
template<typename mf_Policies, 
	 template <typename> class lA2AfieldL,  template <typename> class lA2AfieldR,
	 template <typename> class rA2AfieldL,  template <typename> class rA2AfieldR
	 >
void trace(fMatrix<typename mf_Policies::ScalarComplexType> &into, const std::vector<A2AmesonField<mf_Policies,lA2AfieldL,lA2AfieldR> > &l, const std::vector<A2AmesonField<mf_Policies,rA2AfieldL,rA2AfieldR> > &r){
  //Distribute load over all nodes
  int lsize = l.size();
  int rsize = r.size();

  into.resize(lsize,rsize);

  int nodes = 1; for(int i=0;i<5;i++) nodes *= GJP.Nodes(i);
  int work = lsize*rsize;

  bool do_work = true;
  if(nodes > work){
    nodes = work; if(UniqueID() >= work) do_work = false; //too many nodes, at least for this parallelization. Might want to consider parallelizing in a different way!
  }

  int node_work = work/nodes;
  if(node_work * nodes < work && !UniqueID()) node_work += work - node_work * nodes; //node 0 mops up remainder

  int node_off = UniqueID()*node_work;

#ifndef MEMTEST_MODE
  if(do_work){
    for(int tt=node_off; tt<node_off + node_work; tt++){
      int rem = tt;
      int tsnk = rem % lsize; rem /= lsize; //sink time
      int tsrc = rem; //source time

      into(tsnk,tsrc) = trace(l[tsnk],r[tsrc]);
    }
  }
  into.nodeSum(); //give all nodes a copy
#endif
}


struct nodeDistributeCounter{
  //Get a rank idx that cycles between 0... nodes-1
  static int getNext(){
    static int cur = 0;
    static int nodes = -1;
    if(nodes == -1){
      nodes = 1; for(int i=0;i<5;i++) nodes *= GJP.Nodes(i);
    }
    int out = cur;
    cur = (cur + 1) % nodes;
    return out;
  }
};


//Delete all the data associated with this meson field apart from on node with UniqueID 'node'. The node index is saved so that the data can be later retrieved.
template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::nodeDistribute(int node_uniqueid){
  if(node_uniqueid == -1) node_uniqueid = nodeDistributeCounter::getNext(); //draw the next node index from the pool

  int nodes = 1; for(int i=0;i<5;i++) nodes *= GJP.Nodes(i);
  if(node_uniqueid < 0 || node_uniqueid >= nodes) ERR.General("A2AmesonField","nodeDistribute","Invalid node rank %d\n", node_uniqueid);

#ifndef USE_MPI
  if(nodes > 1) ERR.General("A2AmesonField","nodeDistribute","Implementation requires MPI\n");
  //For one node, don't do anything
#else
  //if(!UniqueID()){ printf("nodeDistribute start with destination node UniqueID %d\n",node_uniqueid); fflush(stdout); }
  MPI_Barrier(MPI_COMM_WORLD);

  //MPI rank and CPS UniqueID may not be the same. Translate
  int my_rank;
  int ret = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if(ret != MPI_SUCCESS) ERR.General("A2AmesonField","nodeDistribute","Comm_rank failed\n");

  if(mf == NULL){ printf("Error nodeDistribute: Rank %d has null pointer!\n",my_rank); fflush(stdout); exit(-1); }

  //printf("UniqueID %d -> MPI rank %d\n",UniqueID(),my_rank); fflush(stdout);

  int rank_tmp = (UniqueID() == node_uniqueid ? my_rank : 0);
  ret = MPI_Allreduce(&rank_tmp,&node_mpi_rank, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); //node is now the MPI rank corresponding to UniqueID == _node
  if(ret != MPI_SUCCESS) ERR.General("A2AmesonField","nodeDistribute","Reduce failed\n");
  MPI_Barrier(MPI_COMM_WORLD);

  //printf("UniqueID %d (MPI rank %d) got storage rank %d\n",UniqueID(),my_rank,node_mpi_rank); fflush(stdout);

  if(UniqueID() != node_uniqueid){
    //printf("UniqueID %d (MPI rank %d) free'd memory\n",UniqueID(),my_rank); fflush(stdout);
    free(mf); mf = NULL;
  }//else{ printf("UniqueID %d (MPI rank %d) is storage node, not freeing\n",UniqueID(),my_rank); fflush(stdout); }

  //if(!UniqueID()) printf("A2AmesonField::nodeDistribute %f MB stored on node %d (MPI rank %d)\n",(double)byte_size()/(1024.*1024.), node_uniqueid, node_mpi_rank);
  //if(my_rank == node_mpi_rank) printf("A2AmesonField::nodeDistribute I am node with MPI rank %d and I have UniqueID %d, my first elem remains %f\n",my_rank,UniqueID(),mf[0]);
#endif
}

#ifdef USE_MPI

template<typename mf_Float>
struct getMPIdataType{
};
template<>
struct getMPIdataType<double>{
  static MPI_Datatype doit(){ return MPI_DOUBLE; }
};
template<>
struct getMPIdataType<float>{
  static MPI_Datatype doit(){ return MPI_FLOAT; }
};

#endif


//Get back the data. After the call, all nodes will have a complete copy
template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::nodeGet(){
  typedef typename ScalarComplexType::value_type mf_Float;
  if(node_mpi_rank == -1) return; //already on all nodes
#ifndef USE_MPI
  int nodes = 1; for(int i=0;i<5;i++) nodes *= GJP.Nodes(i);
  if(nodes > 1) ERR.General("A2AmesonField","nodeGet","Implementation requires MPI\n");
#else
  int mpi_rank; //get this node's MPI rank
  int ret = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if(ret != MPI_SUCCESS) ERR.General("A2AmesonField","nodeGet","Comm_rank failed\n");

  if(mpi_rank != node_mpi_rank){
    //if(mf != NULL) printf("rank %d pointer should be NULL but it isn't!\n",mpi_rank); fflush(stdout);
    mf = (ScalarComplexType*)malloc(byte_size());  
    if(mf == NULL){ printf("rank %d failed to allocate memory!\n",mpi_rank); fflush(stdout); exit(-1); }
    //printf("rank %d allocated memory\n",mpi_rank); fflush(stdout);
  }//else{ printf("rank %d is root, first element of data %f\n",mpi_rank,mf[0]); fflush(stdout); }
  
  MPI_Datatype dtype = getMPIdataType<mf_Float>::doit();
  int dsize;
  MPI_Type_size(dtype,&dsize);
  if(sizeof(mf_Float) != dsize){
    if(!UniqueID()){ printf("MPI size mismatch, want %d, got %d\n",sizeof(mf_Float),dsize); fflush(stdout); }
    MPI_Barrier(MPI_COMM_WORLD);
    exit(-1);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  ret = MPI_Bcast((void*)mf, 2*fsize, dtype, node_mpi_rank, MPI_COMM_WORLD);
  if(ret != MPI_SUCCESS){
    printf("rank %d Bcast fail\n",mpi_rank,mf[0]); fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    exit(-1);
  }
  //printf("rank %d completed Bcast, first element %f\n",mpi_rank,mf[0]); fflush(stdout);

  //if(mpi_rank == node_mpi_rank) printf("A2AmesonField::nodeGet %f MB gathered from node %d (MPI rank %d)\n",(double)byte_size()/(1024.*1024.), UniqueID(), node_mpi_rank);

  node_mpi_rank = -1; //now on all nodes
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}





#endif








