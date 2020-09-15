#ifndef _MULT_VMV_SPLIT_H
#define _MULT_VMV_SPLIT_H

//Try to save memory at the cost of some performance
#define VMV_SPLIT_MEM_SAVE

template<typename ComplexMatrixType>
class SCFoperation{
public:
  virtual void operator()(const ComplexMatrixType& M, const int scf, const int rows, const int cols) = 0;
};
template<typename ScalarComplexType>
class multiply_M_r_op: public SCFoperation<typename gsl_wrapper<typename ScalarComplexType::value_type>::matrix_complex>{
  typedef typename ScalarComplexType::value_type mf_Float;
  typedef gsl_wrapper<mf_Float> gw;
  
  std::vector<std::vector<ScalarComplexType> >& Mr;
  const std::vector<std::vector<ScalarComplexType> >& rreord;
  std::vector<int> const* i_packed_unmap_all; //array of vectors, one for each scf
  
  //Internal
  typename gw::vector_complex* Mr_packed;
  typename gw::complex one;
  typename gw::complex zero;

public:
  multiply_M_r_op(std::vector<std::vector<ScalarComplexType> >& _Mr, const std::vector<std::vector<ScalarComplexType> >& _rreord,
		  std::vector<int> const* _i_packed_unmap_all, const int _nrows_used): Mr(_Mr), rreord(_rreord),i_packed_unmap_all(_i_packed_unmap_all){
    Mr_packed = gw::vector_complex_alloc(_nrows_used);
    GSL_SET_COMPLEX(&one,1.0,0.0);
    GSL_SET_COMPLEX(&zero,0.0,0.0);
  }
  ~multiply_M_r_op(){
    gw::vector_complex_free(Mr_packed);
  }
  
  void operator()(const typename gw::matrix_complex& M_packed, const int scf, const int rows, const int cols){
    const std::vector<int> &i_packed_unmap = i_packed_unmap_all[scf];
    
    int block_width_max =  cols;
    int block_height_max = 8; //4;
          
    for(int i0=0; i0<rows; i0+=block_height_max){
      int iblock_size = std::min(rows - i0, block_height_max);
      for(int j0=0; j0<cols; j0+=block_width_max){
	int jblock_size = std::min(cols - j0, block_width_max);
	
	typename gw::matrix_complex_const_view submatrix = gw::matrix_complex_const_submatrix(&M_packed, i0, j0, iblock_size, jblock_size);
	
	mf_Float const* base = (mf_Float const*)&rreord[scf][j0];
	typename gw::block_complex_struct block;
	block.data = base;
	block.size = jblock_size;
	  
	typename gw::vector_complex rgsl;
	rgsl.block = &block;
	rgsl.data = base;
	rgsl.stride = 1;
	rgsl.owner = 1;
	rgsl.size = jblock_size;
	  
	Mr_packed->size = iblock_size;
	Mr_packed->block->size = iblock_size;
	  
	gw::blas_gemv(CblasNoTrans, one, &submatrix.matrix, &rgsl, zero, Mr_packed);
	  
	typename gw::complex tmp;
	  
	for(int i_packed=0;i_packed < iblock_size; i_packed++){
	  mf_Float(&tmp)[2] = reinterpret_cast<mf_Float(&)[2]>(Mr[scf][ i_packed_unmap[i0+i_packed] ]);
	  mf_Float *t = Mr_packed->data + 2*i_packed*Mr_packed->stride;
	  tmp[0] += *t++; tmp[1] += *t;	       
	}		    
      }
    }
  }

};

template<typename mf_Policies, 
	 template <typename> class lA2AfieldL,  template <typename> class lA2AfieldR,
	 template <typename> class rA2AfieldL,  template <typename> class rA2AfieldR,
	 typename ComplexClass
	 >
class mult_vMv_split_v{};

template<typename mf_Policies, 
	 template <typename> class lA2AfieldL,  template <typename> class lA2AfieldR,
	 template <typename> class rA2AfieldL,  template <typename> class rA2AfieldR	
	 >
class mult_vMv_split: public mult_vMv_split_v<mf_Policies,lA2AfieldL,lA2AfieldR,rA2AfieldL,rA2AfieldR,typename ComplexClassify<typename mf_Policies::ComplexType>::type>{};


template<typename mf_Policies, 
	 template <typename> class lA2AfieldL,  template <typename> class lA2AfieldR,
	 template <typename> class rA2AfieldL,  template <typename> class rA2AfieldR>
class multiply_M_r_singlescf_op: public SCFoperation<typename gsl_wrapper<typename mf_Policies::ScalarComplexType::value_type>::matrix_complex>{
  typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
  typedef gsl_wrapper<typename ScalarComplexType::value_type> gw;
  const int* work; //one for each thread
  const int* off; //one for each thread
  std::vector<  std::vector<std::vector<ScalarComplexType> > > &Mr;
  std::vector< std::vector<std::vector<ScalarComplexType> > > &rreord;
  
  mult_vMv_split_v<mf_Policies,lA2AfieldL,lA2AfieldR,rA2AfieldL,rA2AfieldR,complex_double_or_float_mark> const* split_obj;
public:
  multiply_M_r_singlescf_op(const int* _work, const int* _off, std::vector<  std::vector<std::vector<ScalarComplexType> > > &_Mr, std::vector< std::vector<std::vector<ScalarComplexType> > > &_rreord,mult_vMv_split_v<mf_Policies,lA2AfieldL,lA2AfieldR,rA2AfieldL,rA2AfieldR,complex_double_or_float_mark> const* _split_obj): work(_work),off(_off),Mr(_Mr),rreord(_rreord),split_obj(_split_obj){}
  
  void operator()(const typename gw::matrix_complex& M_packed, const int scf, const int rows, const int cols){
#pragma omp parallel
    {
      int me = omp_get_thread_num();
      split_obj->multiply_M_r_singlescf(Mr,rreord,M_packed,off[me], work[me],scf);
    }
  }
};


template<typename mf_Policies, 
	 template <typename> class lA2AfieldL,  template <typename> class lA2AfieldR,
	 template <typename> class rA2AfieldL,  template <typename> class rA2AfieldR
	 >
class mult_vMv_split_base{
protected:
  typedef typename lA2AfieldL<mf_Policies>::DilutionType iLeftDilutionType;
  typedef typename A2AmesonField<mf_Policies,lA2AfieldR,rA2AfieldL>::LeftDilutionType iRightDilutionType;
  
  typedef typename A2AmesonField<mf_Policies,lA2AfieldR,rA2AfieldL>::RightDilutionType jLeftDilutionType;    
  typedef typename rA2AfieldR<mf_Policies>::DilutionType jRightDilutionType;

  typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
  typedef typename mf_Policies::ComplexType ComplexType;
  
  const static int nscf = 2*3*4;

  //Mapping information
  int ni[nscf], nj[nscf]; //mapping f+2*(c+3*s)
  std::vector<int> ilmap[nscf], irmap[nscf];
  std::vector<int> jlmap[nscf], jrmap[nscf];
    
  std::vector< std::vector<std::pair<int,int> > > blocks_scf; //[scf]  contiguous blocks of 'i' indices

  //Info of the row packing
  bool* rowidx_used;
  int nrows_used; //number of packed rows in the output

  const lA2AfieldL<mf_Policies> *lptr;
  const A2AmesonField<mf_Policies,lA2AfieldR,rA2AfieldL> *Mptr;
  const rA2AfieldR<mf_Policies> *rptr;
  int top_glb;

  int Mrows, Mcols;  

  mult_vMv_split_base():rowidx_used(NULL){}
  
  void setup_base(const lA2AfieldL<mf_Policies> &l,  const A2AmesonField<mf_Policies,lA2AfieldR,rA2AfieldL> &M, const rA2AfieldR<mf_Policies> &r, const int &_top_glb, 
	     const ModeContractionIndices<iLeftDilutionType,iRightDilutionType> &i_ind, const ModeContractionIndices<jLeftDilutionType,jRightDilutionType>& j_ind){
    lptr = &l; rptr = &r; Mptr = &M; top_glb = _top_glb;
  
    modeIndexSet ilp, irp, jlp, jrp;
    ilp.time = top_glb;
    irp.time = M.getRowTimeslice();
    
    jlp.time = M.getColTimeslice();
    jrp.time = top_glb;

    Mrows = M.getNrows();
    Mcols = M.getNcols();

    if(rowidx_used != NULL) free(rowidx_used);
    rowidx_used = (bool*)malloc(Mrows*sizeof(bool)); //Is a particular row of M actually used?
    for(int i=0;i<Mrows;i++) rowidx_used[i] = false;

    //Store maps
    for(int s=0;s<4;s++){
      for(int c=0;c<3;c++){
	int sc = c + 3*s;
	ilp.spin_color = jrp.spin_color = sc;
    	for(int f=0;f<2;f++){
	  ilp.flavor = jrp.flavor = f;

	  int scf = f + 2*ilp.spin_color;

	  //i index
	  int ni_this = i_ind.getNindices(ilp,irp);
	  ni[scf] = ni_this;
	  
	  std::vector<int> &ilmap_this = ilmap[scf]; ilmap_this.resize(ni_this);
	  std::vector<int> &irmap_this = irmap[scf]; irmap_this.resize(ni_this);

	  for(int i = 0; i < ni_this; i++){
	    i_ind.getBothIndices(ilmap_this[i],irmap_this[i],i,ilp,irp);
	    rowidx_used[ irmap_this[i] ] = true; //this row index is used
	  }

	  //j index
	  int nj_this = j_ind.getNindices(jlp,jrp);
	  nj[scf] = nj_this;
	  
	  std::vector<int> &jlmap_this = jlmap[scf]; jlmap_this.resize(nj_this);
	  std::vector<int> &jrmap_this = jrmap[scf]; jrmap_this.resize(nj_this);

	  for(int j = 0; j < nj_this; j++)
	    j_ind.getBothIndices(jlmap_this[j],jrmap_this[j],j,jlp,jrp);
	}
      }
    }
    
    nrows_used = 0;
    for(int i=0;i<Mrows;i++)
      if(rowidx_used[i]) ++nrows_used;
    
    //Get contiguous blocks of i indices
    blocks_scf.resize(nscf);
    for(int scfl=0;scfl<nscf;scfl++){
      int ni_this = ni[scfl];
      find_contiguous_blocks(blocks_scf[scfl],&irmap[scfl][0],ni_this);
    }
  }

  template<typename MatrixVectorComplex>
  void site_reorder_lr(MatrixVectorComplex &lreord,   //[scf][reordered mode]
		       MatrixVectorComplex &rreord,
		       const bool conj_l, const bool conj_r, const int site4dop) const{    
    lreord.resize(nscf); rreord.resize(nscf);

    for(int sc=0;sc<12;sc++){
      for(int f=0;f<2;f++){
	int scf = f + 2*sc;

	//i index
	int ni_this = this->ni[scf];
	const std::vector<int> &ilmap_this = this->ilmap[scf];
	lreord[scf].resize(ni_this);

	for(int i = 0; i < ni_this; i++){
	  const ComplexType &lval_tmp = this->lptr->nativeElem(ilmap_this[i], site4dop, sc, f);
#ifndef MEMTEST_MODE
	  lreord[scf][i] = conj_l ? cconj(lval_tmp) : lval_tmp;
#endif
	}

	//j index
	int nj_this = this->nj[scf];
	const std::vector<int> &jrmap_this = this->jrmap[scf]; //jrmap_this.resize(nj_this);

	rreord[scf].resize(nj_this);
	for(int j = 0; j < nj_this; j++){
	  const ComplexType &rval_tmp = this->rptr->nativeElem(jrmap_this[j], site4dop, sc, f);
#ifndef MEMTEST_MODE
	  rreord[scf][j] = conj_r ? cconj(rval_tmp) : rval_tmp;
#endif
	}

      }
    }
  }

#define FREEIT(A) if(A != NULL){ free(A); A=NULL; }
  
  void free_mem_base(){
    FREEIT(rowidx_used);
  }
  ~mult_vMv_split_base(){
    FREEIT(rowidx_used);
  }
  
};


//For local outer contraction of meson field by two vectors we can save a lot of time by column reordering the meson field to improve cache use. 
//Save even more time by doing this outside the site loop (it makes no reference to the 3d position, only the time at which the vectors
//are evaluated)
template<typename mf_Policies, 
	 template <typename> class lA2AfieldL,  template <typename> class lA2AfieldR,
	 template <typename> class rA2AfieldL,  template <typename> class rA2AfieldR
	 >
class mult_vMv_split_v<mf_Policies,lA2AfieldL,lA2AfieldR,rA2AfieldL,rA2AfieldR, complex_double_or_float_mark>: public mult_vMv_split_base<mf_Policies,lA2AfieldL,lA2AfieldR,rA2AfieldL,rA2AfieldR>{
  typedef mult_vMv_split_base<mf_Policies,lA2AfieldL,lA2AfieldR,rA2AfieldL,rA2AfieldR> Base;
  
  //Note:
  //il is the index of l, 
  //ir is the row index of M, 
  //jl is the column index of M and 
  //jr is the index of r
  typedef typename Base::ScalarComplexType ScalarComplexType;
  
  typedef typename Base::iLeftDilutionType iLeftDilutionType;
  typedef typename Base::iRightDilutionType iRightDilutionType;
  
  typedef typename Base::jLeftDilutionType jLeftDilutionType;    
  typedef typename Base::jRightDilutionType jRightDilutionType;

  typedef typename ScalarComplexType::value_type mf_Float;
  
  //Packed matrices
  typedef gsl_wrapper<mf_Float> gw;

  const static int nscf = 2*3*4;
  
#ifdef VMV_SPLIT_MEM_SAVE
  ScalarComplexType *mf_reord_lo_lo; //shared nl*nl submatrix
  ScalarComplexType *mf_reord_lo_hi[nscf]; //the nl * nh[scf] submatrix
  ScalarComplexType *mf_reord_hi_lo[nscf]; //the nh[scf] * nl submatrix
  ScalarComplexType *mf_reord_hi_hi[nscf]; //the nh[scf] * nh[scf] submatrix
#else
  std::vector<typename gw::matrix_complex*> mf_reord; //vector of gsl matrices in packed format where only the rows used are stored. One matrix for each spin/color/flavor combination of the vector r
#endif

  std::vector<int> i_packed_unmap_all[nscf];
  
  bool setup_called;

  template<typename , 
	   template <typename> class,  template <typename> class,
	   template <typename> class,  template <typename> class>
  friend class multiply_M_r_singlescf_op;

  void constructPackedMloopSCF(SCFoperation<typename gw::matrix_complex> &op){
#ifdef VMV_SPLIT_MEM_SAVE
    int nl_row = this->Mptr->getRowParams().getNl();
    int nl_col = this->Mptr->getColParams().getNl();
    int nj_max = 0;
    bool nj_all_same = true;
    for(int scf=0;scf<nscf;scf++){
      if(this->nj[scf] > nj_max) nj_max = this->nj[scf];
      if(this->nj[scf] != this->nj[0]) nj_all_same = false;
    }
    
    typename gw::matrix_complex* M_packed = gw::matrix_complex_alloc(this->nrows_used,nj_max); //use as buffer
    if(nj_all_same) pokeSubmatrix<ScalarComplexType>( (ScalarComplexType*)M_packed->data, (const ScalarComplexType*)mf_reord_lo_lo, this->nrows_used, nj_max, 0, 0, nl_row, nl_col);
#endif
    
    //M * r
    for(int scf=0;scf<nscf;scf++){
      int nj_this = this->nj[scf]; //vector size
      
#ifdef VMV_SPLIT_MEM_SAVE
      int nh_row = this->nrows_used - nl_row;
      int nh_col = nj_this - nl_col;
      M_packed->size2 = nj_this;
      
      if(!nj_all_same) pokeSubmatrix<ScalarComplexType>( (ScalarComplexType*)M_packed->data, (const ScalarComplexType*)mf_reord_lo_lo, this->nrows_used, nj_this, 0, 0, nl_row, nl_col);
      pokeSubmatrix<ScalarComplexType>( (ScalarComplexType*)M_packed->data, (const ScalarComplexType*)mf_reord_lo_hi[scf], this->nrows_used, nj_this, 0, nl_col, nl_row, nh_col);
      pokeSubmatrix<ScalarComplexType>( (ScalarComplexType*)M_packed->data, (const ScalarComplexType*)mf_reord_hi_lo[scf], this->nrows_used, nj_this, nl_row, 0, nh_row, nl_col);
      pokeSubmatrix<ScalarComplexType>( (ScalarComplexType*)M_packed->data, (const ScalarComplexType*)mf_reord_hi_hi[scf], this->nrows_used, nj_this, nl_row, nl_col, nh_row, nh_col);
#else
      typename gw::matrix_complex* M_packed = mf_reord[scf]; //scope for reuse here
#endif

      op(*M_packed, scf, this->nrows_used, nj_this);
    }
      
#ifdef VMV_SPLIT_MEM_SAVE
    gw::matrix_complex_free(M_packed);
#endif      
  }

  
  void multiply_M_r(std::vector<std::vector<ScalarComplexType> >& Mr, const std::vector<std::vector<ScalarComplexType> >& rreord) const{
    multiply_M_r_op<ScalarComplexType> op(Mr, rreord, this->i_packed_unmap_all, this->nrows_used);
    constructPackedMloopSCF(op);
  }

  //off is the 3d site offset for the start of the internal site loop, and work is the number of sites to iterate over 
  //M_packed is the Mesonfield in packed format.
  void multiply_M_r_singlescf(std::vector<std::vector<std::vector<ScalarComplexType> > >& Mr, const std::vector<std::vector<std::vector<ScalarComplexType> > >& rreord, 
			      const typename gw::matrix_complex & M_packed,
			      const int off, const int work, const int scf) const{
    typename gw::vector_complex* Mr_packed = gw::vector_complex_alloc(this->nrows_used);
    typename gw::complex one; GSL_SET_COMPLEX(&one,1.0,0.0);
    typename gw::complex zero; GSL_SET_COMPLEX(&zero,0.0,0.0);

    //M * r
    int nj_this = this->nj[scf]; //vector size
    const std::vector<int> &i_packed_unmap = this->i_packed_unmap_all[scf];
    
    size_t block_width_max =  M_packed.size2;
    size_t block_height_max = 8; //4;
    
    for(int j0=0; j0<M_packed.size2; j0+=block_width_max){ //columns on outer loop as GSL matrices are row major
      int jblock_size = std::min(M_packed.size2 - j0, block_width_max);
      
      for(int i0=0; i0<M_packed.size1; i0+=block_height_max){
	int iblock_size = std::min(M_packed.size1 - i0, block_height_max);
	
	//if(!me) printf("i0=%d j0=%d  iblock_size=%d jblock_size=%d total rows=%d cols=%d\n",i0,j0,iblock_size,jblock_size,M_packed->size1,M_packed->size2);
	typename gw::matrix_complex_const_view submatrix = gw::matrix_complex_const_submatrix(&M_packed, i0, j0, iblock_size, jblock_size);
	
	for(int s=off;s<off+work;s++){
	  mf_Float* base = (mf_Float*)&rreord[s][scf][j0];
	  typename gw::block_complex_struct block;
	  block.data = base;
	  block.size = jblock_size;
	  
	  typename gw::vector_complex rgsl;
	  rgsl.block = &block;
	  rgsl.data = base;
	  rgsl.stride = 1;
	  rgsl.owner = 1;
	  rgsl.size = jblock_size;
	  
	  Mr_packed->size = iblock_size;
	  Mr_packed->block->size = iblock_size;
	  
	  gw::blas_gemv(CblasNoTrans, one, &submatrix.matrix, &rgsl, zero, Mr_packed);
	  
	  typename gw::complex tmp;
	  
	  for(int i_packed=0;i_packed < iblock_size; i_packed++){
	    mf_Float(&tmp)[2] = reinterpret_cast<mf_Float(&)[2]>(Mr[s][scf][ i_packed_unmap[i0+i_packed] ]);
	    mf_Float *t = Mr_packed->data + 2*i_packed*Mr_packed->stride;
	    tmp[0] += *t++; tmp[1] += *t;	       
	  }
	}	    
      }
    }

    gw::vector_complex_free(Mr_packed);
  }
  
  void site_multiply_l_Mr(CPSspinColorFlavorMatrix<ScalarComplexType> &out, 
			  const std::vector<std::vector<ScalarComplexType> > &lreord,
			  const std::vector<std::vector<ScalarComplexType> > &Mr,
			  typename gw::vector_complex* Mr_gsl_buffer) const{
    //Vector vector multiplication l*(M*r)
    for(int sl=0;sl<4;sl++){
      for(int cl=0;cl<3;cl++){
	for(int fl=0;fl<2;fl++){
	  int scfl = fl + 2*(cl + 3*sl);
	  int ni_this = this->ni[scfl];

	  mf_Float* base = (mf_Float*)&lreord[scfl][0];

	  typename gw::block_complex_struct block;
	  block.data = base;
	  block.size = ni_this;
	
	  typename gw::vector_complex lreord_gsl;
	  lreord_gsl.block = &block;
	  lreord_gsl.data = base;
	  lreord_gsl.stride = 1;
	  lreord_gsl.owner = 1;
	  lreord_gsl.size = ni_this;

	  const std::vector<std::pair<int,int> > &blocks = this->blocks_scf[scfl];

	  for(int sr=0;sr<4;sr++){
	    for(int cr=0;cr<3;cr++){
	      for(int fr=0;fr<2;fr++){
		ScalarComplexType &into = out(sl,sr)(cl,cr)(fl,fr);

		int scfr = fr + 2*(cr + 3*sr);

		typename gw::vector_complex* Mr_gsl = Mr_gsl_buffer;
		Mr_gsl->size = ni_this;
		Mr_gsl->block->size = ni_this;

		ScalarComplexType const* Mr_base = &Mr[scfr][0];

		for(int b=0;b<blocks.size();b++){
		  ScalarComplexType const* block_ptr = Mr_base + this->irmap[scfl][blocks[b].first];
		  mf_Float *t = Mr_gsl->data + 2*blocks[b].first*Mr_gsl->stride;
		  memcpy((void*)t, (void*)block_ptr, 2*blocks[b].second*sizeof(mf_Float));
		}
		typename gw::complex dot;
		gw::blas_dotu(&lreord_gsl, Mr_gsl, &dot);
		    
		reinterpret_cast<typename ScalarComplexType::value_type(&)[2]>(into)[0] = GSL_REAL(dot);
		reinterpret_cast<typename ScalarComplexType::value_type(&)[2]>(into)[1] = GSL_IMAG(dot);
	      }
	    }
	  }
	}
      }
    }
  }


public:

  mult_vMv_split_v(): setup_called(false){
#ifdef VMV_SPLIT_MEM_SAVE
    mf_reord_lo_lo = NULL;
    for(int scf=0;scf<nscf;scf++){
      mf_reord_lo_hi[scf] = NULL;
      mf_reord_hi_lo[scf] = NULL;
      mf_reord_hi_hi[scf] = NULL;
    }
#endif
  }

  void free_mem(){
    this->free_mem_base();
#ifdef VMV_SPLIT_MEM_SAVE
    FREEIT(mf_reord_lo_lo);

    for(int scf=0;scf<nscf;scf++){
    FREEIT(mf_reord_lo_hi[scf]);
    FREEIT(mf_reord_hi_lo[scf]);
    FREEIT(mf_reord_hi_hi[scf]);
    }
#else
    for(int i=0;i<mf_reord.size();i++) gw::matrix_complex_free(mf_reord[i]);
#endif

  }

  ~mult_vMv_split_v(){
#ifdef VMV_SPLIT_MEM_SAVE
    FREEIT(mf_reord_lo_lo);

    for(int scf=0;scf<nscf;scf++){
    FREEIT(mf_reord_lo_hi[scf]);
    FREEIT(mf_reord_hi_lo[scf]);
    FREEIT(mf_reord_hi_hi[scf]);
    }
#else
    for(int i=0;i<mf_reord.size();i++) gw::matrix_complex_free(mf_reord[i]);
#endif
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
    setup_base(l,M,r,_top_glb,i_ind,j_ind);

    //Use GSL BLAS
#ifndef VMV_SPLIT_MEM_SAVE
    mf_reord.resize(nscf);
#endif    
    assert(sizeof(typename gw::complex) == sizeof(std::complex<mf_Float>) ); 

    //Not all rows or columns of M are used, so lets use a packed matrix
    for(int scf=0;scf<nscf;scf++){
      int nj_this = this->nj[scf];
      std::vector<int> &jlmap_this = this->jlmap[scf];
      
      typename gw::matrix_complex* mf_scf_reord = M.GSLpackedColReorder(&jlmap_this.front(), nj_this, this->rowidx_used); //packs the GSL matrix
#ifdef VMV_SPLIT_MEM_SAVE
      int nl_row = M.getRowParams().getNl();
      int nl_col = M.getColParams().getNl();
      int nh_row = this->nrows_used - nl_row;
      int nh_col = nj_this - nl_col;

      if(scf == 0){
	mf_reord_lo_lo = (ScalarComplexType*)malloc(nl_row*nl_col*sizeof(ScalarComplexType));
	getSubmatrix<ScalarComplexType >(mf_reord_lo_lo, (const ScalarComplexType*)mf_scf_reord->data, this->nrows_used, nj_this, 0, 0, nl_row, nl_col);
      }
      mf_reord_lo_hi[scf] = (ScalarComplexType*)malloc(nl_row*nh_col*sizeof(ScalarComplexType));
      getSubmatrix<ScalarComplexType >(mf_reord_lo_hi[scf], (const ScalarComplexType*)mf_scf_reord->data, this->nrows_used, nj_this, 0, nl_col, nl_row, nh_col);

      mf_reord_hi_lo[scf] = (ScalarComplexType*)malloc(nh_row*nl_col*sizeof(ScalarComplexType));
      getSubmatrix<ScalarComplexType >(mf_reord_hi_lo[scf], (const ScalarComplexType*)mf_scf_reord->data, this->nrows_used, nj_this, nl_row, 0, nh_row, nl_col);

      mf_reord_hi_hi[scf] = (ScalarComplexType*)malloc(nh_row*nh_col*sizeof(ScalarComplexType));
      getSubmatrix<ScalarComplexType >(mf_reord_hi_hi[scf], (const ScalarComplexType*)mf_scf_reord->data, this->nrows_used, nj_this, nl_row, nl_col, nh_row, nh_col);

      gw::matrix_complex_free(mf_scf_reord);
#else
      mf_reord[scf] = mf_scf_reord;
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
  void contract(std::vector<CPSspinColorFlavorMatrix<ScalarComplexType>> &out, const bool conj_l, const bool conj_r) const{
    int top = this->top_glb - GJP.TnodeSites()*GJP.TnodeCoor();
    assert(top >= 0 && top < GJP.TnodeSites()); //make sure you use this method on the appropriate node!

    int sites_3d = GJP.VolNodeSites()/GJP.TnodeSites();

    out.resize(sites_3d);

    std::vector< std::vector<std::vector<ScalarComplexType> > > lreord(sites_3d); //[3d site][scf][reordered mode]
    std::vector< std::vector<std::vector<ScalarComplexType> > > rreord(sites_3d);

    assert(sizeof(typename gw::complex) == sizeof(std::complex<mf_Float>) ); 
    typedef gsl_wrapper<mf_Float> gw;
    
    std::vector<  std::vector<std::vector<ScalarComplexType> > > Mr(sites_3d); //[3d site][scf][M row]

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
	    Mr[s][scf][i] = 0.0;
	}
      }
    }

    multiply_M_r_singlescf_op<mf_Policies,lA2AfieldL,lA2AfieldR,rA2AfieldL,rA2AfieldR> op(work,off,Mr,rreord,this);
    constructPackedMloopSCF(op);

    #pragma omp parallel
    {
      int me = omp_get_thread_num();
      //M * r
      //multiply_M_r(&Mr[0],&rreord[0],off,work);
      
      //Vector vector multiplication l*(M*r)
      typename gw::vector_complex* Mr_gsl_buffer = gw::vector_complex_alloc(this->Mrows);
      for(int x3d=off[me];x3d<off[me]+work[me];x3d++)
	site_multiply_l_Mr(out[x3d], lreord[x3d], Mr[x3d], Mr_gsl_buffer);

      gw::vector_complex_free(Mr_gsl_buffer);
      
    } //end of parallel region



  }//end of method


  //Run inside a threaded/parallelized loop over 3d sites. xop is a 3d coordinate!
  void contract(CPSspinColorFlavorMatrix<ScalarComplexType> &out, const int &xop, const bool &conj_l, const bool &conj_r) const{
    int top = this->top_glb - GJP.TnodeSites()*GJP.TnodeCoor();
    assert(top >= 0 && top < GJP.TnodeSites()); //make sure you use this method on the appropriate node!

    std::vector<std::vector<ScalarComplexType > > lreord; //[scf][reordered mode]
    std::vector<std::vector<ScalarComplexType > > rreord;

    assert(sizeof(typename gw::complex) == sizeof(ScalarComplexType) ); 
    typedef gsl_wrapper<mf_Float> gw;
    
    std::vector<std::vector<ScalarComplexType > > Mr(nscf); //[scf][M row]
    for(int scf=0;scf<nscf;scf++){
      Mr[scf].resize(this->Mrows);
      for(int i=0;i<this->Mrows;i++)
	Mr[scf][i] = 0.0;
    }
    int sites_3d = GJP.VolNodeSites()/GJP.TnodeSites();
    int site4dop = xop + sites_3d*top;
    site_reorder_lr(lreord,rreord,conj_l,conj_r,site4dop);

    //M * r
    multiply_M_r(Mr,rreord);

    //Vector vector multiplication l*(M*r)
    typename gw::vector_complex* Mr_gsl_buffer = gw::vector_complex_alloc(this->Mrows);
    site_multiply_l_Mr(out, lreord, Mr, Mr_gsl_buffer);

    gw::vector_complex_free(Mr_gsl_buffer);
  }




};






#endif
