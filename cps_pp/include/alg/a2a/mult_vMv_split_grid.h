#ifndef _MULT_VMV_SPLIT_GRID_H
#define _MULT_VMV_SPLIT_GRID_H
#ifdef USE_GRID

//Try to save memory at the cost of some performance
#define VMV_SPLIT_GRID_MEM_SAVE
//Don't splat-vectorize packed mesonfield at beginning, instead do it at the point of use
//#define VMV_SPLIT_GRID_STREAMING_SPLAT
//Do blocked matrix multiplication
#define VMV_BLOCKED_MATRIX_MULT

template<typename mf_Policies>
class multiply_M_r_op_grid: public SCFoperation<Grid::Vector<
#ifdef VMV_SPLIT_GRID_STREAMING_SPLAT
						  typename mf_Policies::ScalarComplexType
#else
						  typename mf_Policies::ComplexType
#endif
						  > >{
  typedef typename mf_Policies::ComplexType VectorComplexType;
  
#ifdef VMV_SPLIT_GRID_STREAMING_SPLAT
  typedef typename mf_Policies::ScalarComplexType MComplexType;
#else
  typedef typename mf_Policies::ComplexType MComplexType;
#endif
  
  std::vector<Grid::Vector<VectorComplexType> >& Mr;
  const std::vector<Grid::Vector<VectorComplexType> >& rreord;
  std::vector<int> const* i_packed_unmap_all;
  
public:
  multiply_M_r_op_grid(std::vector<Grid::Vector<VectorComplexType> >& _Mr, const std::vector<Grid::Vector<VectorComplexType> >& _rreord,
		  std::vector<int> const* _i_packed_unmap_all): Mr(_Mr), rreord(_rreord),i_packed_unmap_all(_i_packed_unmap_all){
  }
  ~multiply_M_r_op_grid(){
  }
  
  void operator()(const Grid::Vector<MComplexType>& M_packed, const int scf, const int rows, const int cols){
#ifndef MEMTEST_MODE
    const std::vector<int> &i_packed_unmap = i_packed_unmap_all[scf];
# ifdef VMV_SPLIT_GRID_STREAMING_SPLAT
    VectorComplexType tmp;
# endif

# ifdef VMV_BLOCKED_MATRIX_MULT
    int block_width_max = cols;
    int block_height_max = 32; //4;

    //Blocked matrix multiply
    for(int i0=0; i0<rows; i0+=block_height_max){
      int iblock_size = std::min(rows - i0, block_height_max);
	
      for(int j0=0; j0<cols; j0+=block_width_max){
    	int jblock_size = std::min(cols - j0, block_width_max);
	
    	for(int ii=0;ii<iblock_size;ii++){
    	  VectorComplexType &into = Mr[scf][i_packed_unmap[i0+ii]];
    	  for(int jj=0;jj<jblock_size;jj++){
#  ifdef VMV_SPLIT_GRID_STREAMING_SPLAT
	    vsplat(tmp,M_packed[cols*(i0+ii) +j0+jj]);
	    into = into + tmp*rreord[scf][j0+jj];
#  else
    	    into = into + M_packed[cols*(i0+ii) +j0+jj]*rreord[scf][j0+jj];
#  endif
	  }
    	}
      }
    }

# else
    
    for(int i=0;i<rows;i++){
      VectorComplexType &into = Mr[scf][i_packed_unmap[i]];
      zeroit(into);
      for(int j=0;j<cols;j++){
#  ifdef VMV_SPLIT_GRID_STREAMING_SPLAT
	vsplat(tmp,M_packed[cols*i + j]);
	into = into + tmp*rreord[scf][j];
#  else
    	into = into + M_packed[cols*i + j]*rreord[scf][j];
#  endif
      }
    }
    
# endif
    
#endif
    
  }
};


template<typename mf_Policies, 
	 template <typename> class lA2AfieldL,  template <typename> class lA2AfieldR,
	 template <typename> class rA2AfieldL,  template <typename> class rA2AfieldR>
class multiply_M_r_singlescf_op_grid: public SCFoperation<Grid::Vector<
#ifdef VMV_SPLIT_GRID_STREAMING_SPLAT
						  typename mf_Policies::ScalarComplexType
#else
						  typename mf_Policies::ComplexType
#endif
							    > >{
  typedef typename mf_Policies::ComplexType VectorComplexType;

#ifdef VMV_SPLIT_GRID_STREAMING_SPLAT
  typedef typename mf_Policies::ScalarComplexType MComplexType;
#else
  typedef typename mf_Policies::ComplexType MComplexType;
#endif
  
  const int* work;
  const int* off;
  std::vector< std::vector<Grid::Vector<VectorComplexType> > > &Mr;
  std::vector< std::vector<Grid::Vector<VectorComplexType> > > &rreord;
  
  mult_vMv_split_v<mf_Policies,lA2AfieldL,lA2AfieldR,rA2AfieldL,rA2AfieldR, grid_vector_complex_mark> const* split_obj;
public:
  multiply_M_r_singlescf_op_grid(const int* _work, const int* _off, std::vector<  std::vector<Grid::Vector<VectorComplexType> > > &_Mr, std::vector< std::vector<Grid::Vector<VectorComplexType> > > &_rreord,mult_vMv_split_v<mf_Policies,lA2AfieldL,lA2AfieldR,rA2AfieldL,rA2AfieldR, grid_vector_complex_mark> const* _split_obj): work(_work),off(_off),Mr(_Mr),rreord(_rreord),split_obj(_split_obj){}
  
  void operator()(const Grid::Vector<MComplexType>& M_packed, const int scf, const int rows, const int cols){
#ifndef MEMTEST_MODE
#pragma omp parallel
    {
      int me = omp_get_thread_num();
      split_obj->multiply_M_r_singlescf(Mr,rreord,M_packed,off[me], work[me],scf);
    }
#endif
  }
};



//For local outer contraction of meson field by two vectors we can save a lot of time by column reordering the meson field to improve cache use. 
//Save even more time by doing this outside the site loop (it makes no reference to the 3d position, only the time at which the vectors
//are evaluated)
template<typename mf_Policies, 
	 template <typename> class lA2AfieldL,  template <typename> class lA2AfieldR,
	 template <typename> class rA2AfieldL,  template <typename> class rA2AfieldR
	 >
class mult_vMv_split_v<mf_Policies, lA2AfieldL, lA2AfieldR, rA2AfieldL, rA2AfieldR, grid_vector_complex_mark>:
public mult_vMv_split_base<mf_Policies,lA2AfieldL,lA2AfieldR,rA2AfieldL,rA2AfieldR>{
  typedef mult_vMv_split_base<mf_Policies,lA2AfieldL,lA2AfieldR,rA2AfieldL,rA2AfieldR> Base;
  
  //Note:
  //il is the index of l, 
  //ir is the row index of M, 
  //jl is the column index of M and 
  //jr is the index of r
  typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
  typedef typename mf_Policies::ComplexType VectorComplexType;

#ifdef VMV_SPLIT_GRID_STREAMING_SPLAT
  typedef typename mf_Policies::ScalarComplexType MComplexType;
#else
  typedef typename mf_Policies::ComplexType MComplexType;
#endif
  
  typedef typename Base::iLeftDilutionType iLeftDilutionType;
  typedef typename Base::iRightDilutionType iRightDilutionType;
  
  typedef typename Base::jLeftDilutionType jLeftDilutionType;    
  typedef typename Base::jRightDilutionType jRightDilutionType;

  const static int nscf = 2*3*4;

  //Packed matrices
#ifdef VMV_SPLIT_GRID_MEM_SAVE
  Grid::Vector<MComplexType> mf_reord_lo_lo; //shared nl*nl submatrix
  Grid::Vector<MComplexType> mf_reord_lo_hi[nscf]; //the nl * nh[scf] submatrix
  Grid::Vector<MComplexType> mf_reord_hi_lo[nscf]; //the nh[scf] * nl submatrix
  Grid::Vector<MComplexType> mf_reord_hi_hi[nscf]; //the nh[scf] * nh[scf] submatrix
#else
  std::vector<Grid::Vector<MComplexType> > mf_reord; //vector of linearized matrices in packed format where only the rows used are stored. One matrix for each spin/color/flavor combination of the vector r
#endif

  std::vector<int> i_packed_unmap_all[nscf];
  int logical_sites_3d;
 
  bool setup_called;

  template<typename , 
	   template <typename> class,  template <typename> class,
	   template <typename> class,  template <typename> class>
  friend class multiply_M_r_singlescf_op_grid;

  void constructPackedMloopSCF(SCFoperation<Grid::Vector<MComplexType> > &op) const{
#ifdef VMV_SPLIT_GRID_MEM_SAVE
    int nl_row = this->Mptr->getRowParams().getNl();
    int nl_col = this->Mptr->getColParams().getNl();
    int nj_max = 0;
    bool nj_all_same = true;
    for(int scf=0;scf<nscf;scf++){
      if(this->nj[scf] > nj_max) nj_max = this->nj[scf];
      if(this->nj[scf] != this->nj[0]) nj_all_same = false;
    }
    Grid::Vector<MComplexType> M_packed(this->nrows_used * nj_max); //linearized matrix  j + nj_max*i. Acts as buffer. Contains enough space for largest nj

    if(nj_all_same) pokeSubmatrix<MComplexType>( M_packed.data(),mf_reord_lo_lo.data(), this->nrows_used, nj_max, 0, 0, nl_row, nl_col);
#endif
    
    //M * r
    for(int scf=0;scf<nscf;scf++){
      int nj_this = this->nj[scf]; //vector size
      
#ifdef VMV_SPLIT_GRID_MEM_SAVE
      int nh_row = this->nrows_used - nl_row;
      int nh_col = nj_this - nl_col;
      
      if(!nj_all_same) pokeSubmatrix<MComplexType>(M_packed.data(), mf_reord_lo_lo.data(), this->nrows_used, nj_this, 0, 0, nl_row, nl_col);
      pokeSubmatrix<MComplexType>(M_packed.data(), mf_reord_lo_hi[scf].data(), this->nrows_used, nj_this, 0, nl_col, nl_row, nh_col);
      pokeSubmatrix<MComplexType>(M_packed.data(), mf_reord_hi_lo[scf].data(), this->nrows_used, nj_this, nl_row, 0, nh_row, nl_col);
      pokeSubmatrix<MComplexType>(M_packed.data(), mf_reord_hi_hi[scf].data(), this->nrows_used, nj_this, nl_row, nl_col, nh_row, nh_col);
#else
      const Grid::Vector<MComplexType>& M_packed = mf_reord[scf]; //scope for reuse here
#endif

      op(M_packed, scf, this->nrows_used, nj_this);
    }  
  }

  
  //off is the 3d site offset for the start of the internal site loop, and work is the number of sites to iterate over
  //Mr is part of an array of length nsites. For each site there are nscf=24 vectors of the appropriate number of modes
  void multiply_M_r(std::vector<Grid::Vector<VectorComplexType> > &Mr, const std::vector<Grid::Vector<VectorComplexType> >& rreord) const{
    multiply_M_r_op_grid<mf_Policies> op(Mr, rreord, this->i_packed_unmap_all);
    constructPackedMloopSCF(op);
  }

  
  //off is the 3d site offset for the start of the internal site loop, and work is the number of sites to iterate over 
  //M_packed is the Mesonfield in packed format.
  void multiply_M_r_singlescf(std::vector<std::vector<Grid::Vector<VectorComplexType> > >& Mr, const std::vector<std::vector<Grid::Vector<VectorComplexType> > >& rreord, 
			      const Grid::Vector<MComplexType>& M_packed,
			      const int off, const int work, const int scf) const{
    //M * r
    int nj_this = this->nj[scf]; //vector size
    const std::vector<int> &i_packed_unmap = i_packed_unmap_all[scf];

    int rows = this->nrows_used;
    int cols = nj_this;

#ifdef VMV_SPLIT_GRID_STREAMING_SPLAT
    VectorComplexType tmp;
#endif

#ifdef VMV_BLOCKED_MATRIX_MULT
    int block_width_max =  cols;
    int block_height_max = 8; //4;
          
    for(int i0=0; i0<rows; i0+=block_height_max){
      int iblock_size = std::min(rows - i0, block_height_max);
      
      for(int j0=0; j0<cols; j0+=block_width_max){
    	int jblock_size = std::min(cols - j0, block_width_max);

    	for(int s=off;s<off+work;s++){
    	  VectorComplexType const* base = &rreord[s][scf][j0];
# ifndef MEMTEST_MODE
    	  for(int i_packed=0;i_packed < iblock_size; i_packed++){
    	    VectorComplexType &into = Mr[s][scf][ i_packed_unmap[i0+i_packed] ];
    	    zeroit(into);
	    
    	    for(int j_packed=0;j_packed<jblock_size;j_packed++){
#  ifdef VMV_SPLIT_GRID_STREAMING_SPLAT
	      vsplat(tmp, M_packed[cols*(i0+i_packed)+j0+j_packed]);
	      into = into + tmp*base[j_packed];
#  else
    	      into = into + M_packed[cols*(i0+i_packed)+j0+j_packed]*base[j_packed];
#  endif
	    }
    	  }
# endif
    	}
      }
    }
#else
    
    for(int i=0;i<rows;i++){
      int i_full = i_packed_unmap[i];

      for(int s=off;s<off+work;s++){
# ifndef MEMTEST_MODE
	VectorComplexType &into = Mr[s][scf][i_full];
	zeroit(into);
	
	for(int j=0;j<cols;j++){
#  ifdef VMV_SPLIT_GRID_STREAMING_SPLAT
	  vsplat(tmp,M_packed[cols*i+j]);
	  into = into + tmp * rreord[s][scf][j];
#  else
	  into = into + M_packed[cols*i+j] * rreord[s][scf][j];
#  endif
	}
# endif
      }
    }
	  
#endif

  }
  
  void site_multiply_l_Mr(CPSspinColorFlavorMatrix<VectorComplexType> &out, 
			  const std::vector<Grid::Vector<VectorComplexType> > &lreord,
			  const std::vector<Grid::Vector<VectorComplexType> > &Mr) const{
    //Vector vector multiplication l*(M*r)
    for(int sl=0;sl<4;sl++){
      for(int sr=0;sr<4;sr++){
	for(int cl=0;cl<3;cl++){
	  for(int cr=0;cr<3;cr++){
	    for(int fl=0;fl<2;fl++){
	      int scfl = fl + 2*(cl + 3*sl);
	      int ni_this = this->ni[scfl];

	      Grid::Vector<VectorComplexType> const& lbase = lreord[scfl];
	      const std::vector<std::pair<int,int> > &blocks = this->blocks_scf[scfl];
	      
	      for(int fr=0;fr<2;fr++){
		int scfr = fr + 2*(cr + 3*sr);

		VectorComplexType &into = out(sl,sr)(cl,cr)(fl,fr);
		zeroit(into);

		VectorComplexType const* Mr_base = &Mr[scfr][0];

		int loff = 0;
		for(int b=0;b<blocks.size();b++){
		  VectorComplexType const* Mr_block_ptr = Mr_base + this->irmap[scfl][blocks[b].first]; //Mr is not packed, lreord is. Cycle over blocks of consecutive elements
#ifndef MEMTEST_MODE
		  for(int i=0;i<blocks[b].second;i++)
		    into = into + lbase[loff++] * Mr_block_ptr[i];		  
#endif
		}
	      }
	    }
	  }
	}
      }
    }
  }


public:

  mult_vMv_split_v(): setup_called(false){

  }

  void free_mem(){
    this->free_mem_base();
  }


  //This should be called outside the site loop (and outside any parallel region)
  //top_glb is the time in global lattice coordinates.
  void setup(const lA2AfieldL<mf_Policies> &l,  const A2AmesonField<mf_Policies,lA2AfieldR,rA2AfieldL> &M, const rA2AfieldR<mf_Policies> &r, const int &_top_glb){
    //Precompute index mappings
    ModeContractionIndices<iLeftDilutionType,iRightDilutionType> i_ind(l);
    ModeContractionIndices<jLeftDilutionType,jRightDilutionType> j_ind(r);
    setup(l,M,r,_top_glb,i_ind,j_ind);
  }
  void setup(const lA2AfieldL<mf_Policies> &l,  const A2AmesonField<mf_Policies,lA2AfieldR,rA2AfieldL> &M, const rA2AfieldR<mf_Policies> &r, const int &_top_glb, 
	     const ModeContractionIndices<iLeftDilutionType,iRightDilutionType> &i_ind, const ModeContractionIndices<jLeftDilutionType,jRightDilutionType>& j_ind){
    this->setup_base(l,M,r,_top_glb,i_ind,j_ind);

    logical_sites_3d = l.getMode(0).nodeSites(0)*l.getMode(0).nodeSites(1)*l.getMode(0).nodeSites(2);
#ifndef VMV_SPLIT_GRID_MEM_SAVE
    mf_reord.resize(nscf);
#endif    

    //Not all rows or columns of M are used, so lets use a packed matrix

#ifdef VMV_SPLIT_GRID_MEM_SAVE
    Grid::Vector<MComplexType> mf_reord_scf;
#endif
    
    for(int scf=0;scf<nscf;scf++){
      int nj_this = this->nj[scf];
      std::vector<int> &jlmap_this = this->jlmap[scf];

#ifndef VMV_SPLIT_GRID_MEM_SAVE
      Grid::Vector<MComplexType> &mf_reord_scf = mf_reord[scf];
#endif

#ifdef VMV_SPLIT_GRID_STREAMING_SPLAT
      M.scalarPackedColReorder(mf_reord_scf, jlmap_this.data(), nj_this, this->rowidx_used);
#else
      M.splatPackedColReorder(mf_reord_scf, jlmap_this.data(), nj_this, this->rowidx_used);
#endif

#ifdef VMV_SPLIT_GRID_MEM_SAVE
      int nl_row = M.getRowParams().getNl();
      int nl_col = M.getColParams().getNl();
      int nh_row = this->nrows_used - nl_row;
      int nh_col = nj_this - nl_col;

      if(scf == 0){
	mf_reord_lo_lo.resize(nl_row*nl_col);
	getSubmatrix<MComplexType >(mf_reord_lo_lo.data(), mf_reord_scf.data(), this->nrows_used, nj_this, 0, 0, nl_row, nl_col);
      }
      mf_reord_lo_hi[scf].resize(nl_row*nh_col);
      getSubmatrix<MComplexType >(mf_reord_lo_hi[scf].data(), mf_reord_scf.data(), this->nrows_used, nj_this, 0, nl_col, nl_row, nh_col);

      mf_reord_hi_lo[scf].resize(nh_row*nl_col);
      getSubmatrix<MComplexType >(mf_reord_hi_lo[scf].data(), mf_reord_scf.data(), this->nrows_used, nj_this, nl_row, 0, nh_row, nl_col);

      mf_reord_hi_hi[scf].resize(nh_row*nh_col);
      getSubmatrix<MComplexType >(mf_reord_hi_hi[scf].data(), mf_reord_scf.data(), this->nrows_used, nj_this, nl_row, nl_col, nh_row, nh_col);
#endif
             
      //Store the map between packed and full indices
      int i_packed = 0;
      i_packed_unmap_all[scf].resize(this->nrows_used);
      for(int i_full=0;i_full<this->Mrows;i_full++) 
	if(this->rowidx_used[i_full]) i_packed_unmap_all[scf][i_packed++] = i_full;
    }
    
    setup_called = true;
  }

public:
  //Contract on all 3d sites on this node with fixed operator time coord top_glb into a canonically ordered output vector
  void contract(Grid::Vector<CPSspinColorFlavorMatrix<VectorComplexType> > &out, const bool conj_l, const bool conj_r) const{
    int top = this->top_glb - GJP.TnodeSites()*GJP.TnodeCoor();
    assert(top >= 0 && top < GJP.TnodeSites()); //make sure you use this method on the appropriate node!

    int sites_3d = logical_sites_3d;

    out.resize(sites_3d);

    std::vector< std::vector<Grid::Vector<VectorComplexType> > > lreord(sites_3d); //[3d site][scf][reordered mode]
    std::vector< std::vector<Grid::Vector<VectorComplexType> > > rreord(sites_3d);

    std::vector<  std::vector<Grid::Vector<VectorComplexType> > > Mr(sites_3d); //[3d site][scf][M row]

    //Run everything in parallel environment to avoid thread creation overheads
    int work[omp_get_max_threads()], off[omp_get_max_threads()];

#pragma omp parallel
    {
      int me = omp_get_thread_num();
      int team = omp_get_num_threads();
   
      //Reorder rows and columns of left and right fields such that they can be accessed sequentially
      thread_work(work[me], off[me], sites_3d, me, team);

      for(int s=off[me];s<off[me]+work[me];s++){
	int site4dop = s + sites_3d*top;
	site_reorder_lr(lreord[s],rreord[s],conj_l,conj_r,site4dop);
      }
      
      for(int s=off[me];s<off[me]+work[me];s++){
	Mr[s].resize(nscf); 
	for(int scf=0;scf<nscf;scf++){
	  Mr[s][scf].resize(this->Mrows);
	  for(int i=0;i<this->Mrows;i++)
	    zeroit(Mr[s][scf][i]);
	}
      }
    }

    multiply_M_r_singlescf_op_grid<mf_Policies,lA2AfieldL,lA2AfieldR,rA2AfieldL,rA2AfieldR> op(work,off,Mr,rreord,this);
    constructPackedMloopSCF(op);

#pragma omp parallel
    {
      int me = omp_get_thread_num();
      
      //Vector vector multiplication l*(M*r)
      for(int x3d=off[me];x3d<off[me]+work[me];x3d++)
	site_multiply_l_Mr(out[x3d], lreord[x3d], Mr[x3d]);
      
    } //end of parallel region



  }//end of method


  //Run inside a threaded/parallelized loop over 3d sites. xop is a 3d coordinate!
  void contract(CPSspinColorFlavorMatrix<VectorComplexType> &out, const int &xop, const bool &conj_l, const bool &conj_r) const{
    int top = this->top_glb - GJP.TnodeSites()*GJP.TnodeCoor();
    assert(top >= 0 && top < GJP.TnodeSites()); //make sure you use this method on the appropriate node!

    std::vector<Grid::Vector<VectorComplexType > > lreord; //[scf][reordered mode]
    std::vector<Grid::Vector<VectorComplexType > > rreord;

    std::vector<Grid::Vector<VectorComplexType > > Mr(nscf); //[scf][M row]
    for(int scf=0;scf<nscf;scf++){
      Mr[scf].resize(this->Mrows);
      for(int i=0;i<this->Mrows;i++)
	zeroit(Mr[scf][i]);
    }
    int sites_3d = logical_sites_3d;
    int site4dop = xop + sites_3d*top;
    this->site_reorder_lr(lreord,rreord,conj_l,conj_r,site4dop);

    //M * r
    multiply_M_r(Mr,rreord);

    //Vector vector multiplication l*(M*r)
    site_multiply_l_Mr(out, lreord, Mr);
  }





};




#endif

#endif
