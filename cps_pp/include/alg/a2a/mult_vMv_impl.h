#ifndef _MULT_VMV_IMPL
#define _MULT_VMV_IMPL

//Vector mesonfield outer product implementation

template<typename mf_Policies, 
	 template <typename> class lA2AfieldL,  template <typename> class lA2AfieldR,
	 template <typename> class rA2AfieldL,  template <typename> class rA2AfieldR,
	 typename ComplexClass>
class _mult_vMv_impl_v{};

template<typename mf_Policies, 
	 template <typename> class lA2AfieldL,  template <typename> class lA2AfieldR,
	 template <typename> class rA2AfieldL,  template <typename> class rA2AfieldR
	 >
class _mult_vMv_impl_v<mf_Policies,lA2AfieldL,lA2AfieldR,rA2AfieldL,rA2AfieldR,complex_double_or_float_mark>{ //necessary to avoid an annoying ambigous overload when mesonfield friends mult
public:
  typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
  //Form SpinColorFlavorMatrix prod1 = vL_i(\vec xop, top ; tpi2) [\sum_{\vec xpi2} wL_i^dag(\vec xpi2, tpi2) S2 vL_j(\vec xpi2, tpi2; top)] wL_j^dag(\vec xop,top)

  // l^i(xop,top) M^ij(tl,tr) r^j(xop,top)
  //argument xop is the *local* 3d site index in canonical ordering, top is the *local* time coordinate
  // Node local and unthreaded
  static void mult(CPSspinColorFlavorMatrix<ScalarComplexType> &out, const lA2AfieldL<mf_Policies> &l,  const A2AmesonField<mf_Policies,lA2AfieldR,rA2AfieldL> &M, const rA2AfieldR<mf_Policies> &r, const int &xop, const int &top, const bool &conj_l, const bool &conj_r){
    typedef typename lA2AfieldL<mf_Policies>::DilutionType iLeftDilutionType;
    typedef typename A2AmesonField<mf_Policies,lA2AfieldR,rA2AfieldL>::LeftDilutionType iRightDilutionType;

    typedef typename A2AmesonField<mf_Policies,lA2AfieldR,rA2AfieldL>::RightDilutionType jLeftDilutionType;    
    typedef typename rA2AfieldR<mf_Policies>::DilutionType jRightDilutionType;

    out.zero();

    int top_glb = top+GJP.TnodeSites()*GJP.TnodeCoor();

    //Precompute index mappings
    ModeContractionIndices<iLeftDilutionType,iRightDilutionType> i_ind(l);
    ModeContractionIndices<jLeftDilutionType,jRightDilutionType> j_ind(r);

    modeIndexSet ilp, irp, jlp, jrp;
    ilp.time = top_glb;
    irp.time = M.getRowTimeslice();
    
    jlp.time = M.getColTimeslice();
    jrp.time = top_glb;

    int site4dop = xop + GJP.VolNodeSites()/GJP.TnodeSites()*top;

    const static int nscf = 2*3*4;
    const static int complex_scf_vect_size = nscf*2;

    int ni[nscf], nj[nscf]; //mapping f+2*(c+3*s)
    std::vector<int> ilmap[nscf], irmap[nscf], jlmap[nscf], jrmap[nscf];

    //Reorder rows and columns such that they can be accessed sequentially
    std::vector<ScalarComplexType> lreord[nscf];
    std::vector<ScalarComplexType> rreord[nscf];

    int Mrows = M.getNrows();
    int Mcols = M.getNcols();

    //Are particular row and column indices of M actually used?
    bool rowidx_used[Mrows]; for(int i=0;i<Mrows;i++) rowidx_used[i] = false;
    bool colidx_used[Mcols]; for(int i=0;i<Mcols;i++) colidx_used[i] = false;
    
    //Note:
    //il is the index of l, 
    //ir is the row index of M, 
    //jl is the column index of M and 
    //jr is the index of r

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

	  //M.rowReorder(rowreord[scf], &irmap_this.front(), ni_this);
	  lreord[scf].resize(ni_this);
	  for(int i = 0; i < ni_this; i++){
	    const ScalarComplexType &lval_tmp = l.nativeElem(ilmap_this[i], site4dop, sc, f);
	    lreord[scf][i] = conj_l ? std::conj(lval_tmp) : lval_tmp;
	  }

	  //j index
	  int nj_this = j_ind.getNindices(jlp,jrp);
	  nj[scf] = nj_this;
	  
	  std::vector<int> &jlmap_this = jlmap[scf]; jlmap_this.resize(nj_this);
	  std::vector<int> &jrmap_this = jrmap[scf]; jrmap_this.resize(nj_this);

	  for(int j = 0; j < nj_this; j++){
	    j_ind.getBothIndices(jlmap_this[j],jrmap_this[j],j,jlp,jrp);
	    colidx_used[ jlmap_this[j] ] = true;
	  }

	  //M.colReorder(colreord[scf], &jlmap_this.front(), nj_this);
	  rreord[scf].resize(nj_this);
	  for(int j = 0; j < nj_this; j++){
	    const ScalarComplexType &rval_tmp = r.nativeElem(jrmap_this[j], site4dop, sc, f);
	    rreord[scf][j] = conj_r ? std::conj(rval_tmp) : rval_tmp;
	  }

	}
      }
    }	  


    //Matrix vector multiplication  M*r
    ScalarComplexType Mr[Mrows][nscf];

    //Use GSL BLAS
    typedef gsl_wrapper<typename ScalarComplexType::value_type> gw;

    assert(sizeof(typename gw::complex) == sizeof(ScalarComplexType) );

    typename gw::complex tmp;

    //Not all rows or columns of M are used, so lets use a packed matrix
    int nrows_used = 0;
    for(int i_full=0;i_full<Mrows;i_full++) if(rowidx_used[i_full]) nrows_used++;

    typename gw::complex one; GSL_SET_COMPLEX(&one,1.0,0.0);
    typename gw::complex zero; GSL_SET_COMPLEX(&zero,0.0,0.0);

    typename gw::vector_complex *Mr_packed = gw::vector_complex_alloc(nrows_used);
    typename gw::matrix_complex *M_packed_buffer = gw::matrix_complex_alloc(nrows_used,Mcols); //matrix cannot be larger than this so we can reuse the memory

    for(int scf=0;scf<nscf;scf++){
      int nj_this = nj[scf];
      
      std::vector<int> &jlmap_this = jlmap[scf];
      
      typename gw::matrix_complex * M_packed = M.GSLpackedColReorder(&jlmap_this.front(), nj_this, rowidx_used, M_packed_buffer); //packs the GSL matrix
             
      int i_packed = 0;
      int i_packed_unmap[nrows_used];
      for(int i_full=0;i_full<Mrows;i_full++)
    	if(rowidx_used[i_full]) i_packed_unmap[i_packed++] = i_full;
      
      typename ScalarComplexType::value_type* base = (typename ScalarComplexType::value_type*)&rreord[scf][0];
      typename gw::block_complex_struct block;
      block.data = base;
      block.size = nj_this;

      typename gw::vector_complex rgsl;
      rgsl.block = &block;
      rgsl.data = base;
      rgsl.stride = 1;
      rgsl.owner = 1;
      rgsl.size = nj_this;

      gw::blas_gemv(CblasNoTrans, one, M_packed, &rgsl, zero, Mr_packed);

      for(int i_packed=0;i_packed < nrows_used; i_packed++){
    	tmp = gw::vector_complex_get(Mr_packed,i_packed);
    	Mr[ i_packed_unmap[i_packed] ][scf] = ScalarComplexType( GSL_REAL(tmp), GSL_IMAG(tmp) );
      }
    }
    gw::vector_complex_free(Mr_packed);
    gw::matrix_complex_free(M_packed_buffer);

    //Vector vector multiplication l*(M*r)
    for(int sl=0;sl<4;sl++){
      for(int sr=0;sr<4;sr++){
	for(int cl=0;cl<3;cl++){
	  for(int cr=0;cr<3;cr++){
	    for(int fl=0;fl<2;fl++){
	      int scfl = fl + 2*(cl + 3*sl);
	      int ni_this = ni[scfl];
	      for(int fr=0;fr<2;fr++){
		int scfr = fr + 2*(cr + 3*sr);
		ScalarComplexType &into = out(sl,sr)(cl,cr)(fl,fr);
		for(int i=0;i<ni_this;i++){
		  int i_full = irmap[scfl][i];
		  into += lreord[scfl][i] * Mr[i_full][scfr];
		}
	      }
	    }
	  }
	}
      }
    }	    

  }

  static void mult_slow(CPSspinColorFlavorMatrix<ScalarComplexType> &out, const lA2AfieldL<mf_Policies> &l,  const A2AmesonField<mf_Policies,lA2AfieldR,rA2AfieldL> &M, const rA2AfieldR<mf_Policies> &r, const int &xop, const int &top, const bool &conj_l, const bool &conj_r){

    int site4dop = xop + GJP.VolNodeSites()/GJP.TnodeSites()*top;

    A2Aparams i_params(l), j_params(r);
    StandardIndexDilution idil(i_params), jdil(j_params);
    
    int ni = idil.getNmodes();
    int nj = jdil.getNmodes();

    out.zero();

    for(int sl=0;sl<4;sl++){
      for(int sr=0;sr<4;sr++){
	for(int cl=0;cl<3;cl++){
	  for(int cr=0;cr<3;cr++){
	    for(int fl=0;fl<2;fl++){	  
	      for(int fr=0;fr<2;fr++){

		for(int i=0;i<ni;i++){
		  const ScalarComplexType &lval_tmp = l.elem(i,xop,top,cl+3*sl,fl);
		  ScalarComplexType lval = conj_l ? std::conj(lval_tmp) : lval_tmp;
		  
  		  for(int j=0;j<nj;j++){
  		    const ScalarComplexType &rval_tmp = r.elem(j,xop,top,cr+3*sr,fr);
  		    ScalarComplexType rval = conj_r ? std::conj(rval_tmp) : rval_tmp;

		    const ScalarComplexType &Mval = M.elem(i,j);
		    ScalarComplexType delta = lval * Mval * rval;
  		    out(sl,sr)(cl,cr)(fl,fr) += delta;
  		  }
  		}
  	      }
  	    }
  	  }
  	}
      }
    }

  }
};




#ifdef USE_GRID

template<typename mf_Policies, 
	 template <typename> class lA2AfieldL,  template <typename> class lA2AfieldR,
	 template <typename> class rA2AfieldL,  template <typename> class rA2AfieldR>
class _mult_vMv_impl_v<mf_Policies,lA2AfieldL,lA2AfieldR,rA2AfieldL,rA2AfieldR,grid_vector_complex_mark>{ //for SIMD vectorized W and V vectors
public:
  typedef typename mf_Policies::ComplexType VectorComplexType;
  typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
  //Form SpinColorFlavorMatrix prod1 = vL_i(\vec xop, top ; tpi2) [\sum_{\vec xpi2} wL_i^dag(\vec xpi2, tpi2) S2 vL_j(\vec xpi2, tpi2; top)] wL_j^dag(\vec xop,top)

  // l^i(xop,top) M^ij(tl,tr) r^j(xop,top)
  //argument xop is the *local* 3d site index of the reduced (logical) lattice in canonical ordering, top is the *local* time coordinate
  //it is assumed that the lattice is not SIMD vectorized in the time direction
  // Node local and unthreaded
  static void mult(CPSspinColorFlavorMatrix<VectorComplexType> &out, const lA2AfieldL<mf_Policies> &l,  const A2AmesonField<mf_Policies,lA2AfieldR,rA2AfieldL> &M, const rA2AfieldR<mf_Policies> &r, const int &xop, const int &top, const bool &conj_l, const bool &conj_r){
    typedef typename lA2AfieldL<mf_Policies>::DilutionType iLeftDilutionType;
    typedef typename A2AmesonField<mf_Policies,lA2AfieldR,rA2AfieldL>::LeftDilutionType iRightDilutionType;

    typedef typename A2AmesonField<mf_Policies,lA2AfieldR,rA2AfieldL>::RightDilutionType jLeftDilutionType;    
    typedef typename rA2AfieldR<mf_Policies>::DilutionType jRightDilutionType;

    assert(l.getMode(0).SIMDlogicalNodes(3) == 1);
    
    out.zero();
    
    int top_glb = top+GJP.TnodeSites()*GJP.TnodeCoor();

    //Precompute index mappings
    ModeContractionIndices<iLeftDilutionType,iRightDilutionType> i_ind(l);
    ModeContractionIndices<jLeftDilutionType,jRightDilutionType> j_ind(r);

    modeIndexSet ilp, irp, jlp, jrp;
    ilp.time = top_glb;
    irp.time = M.getRowTimeslice();
    
    jlp.time = M.getColTimeslice();
    jrp.time = top_glb;

    int site4dop = l.getMode(0).threeToFour(xop, top);

    const static int nscf = 2*3*4;
    int ni[nscf], nj[nscf]; //mapping f+2*(c+3*s)
    std::vector<int> ilmap[nscf], irmap[nscf], jlmap[nscf], jrmap[nscf];

    const int Mrows = M.getNrows();

    //Are particular row indices of M actually used?
    bool rowidx_used[Mrows]; for(int i=0;i<Mrows;i++) rowidx_used[i] = false;

    //Reorder rows and columns such that they can be accessed sequentially
    Grid::Vector<VectorComplexType> lreord[nscf];
    Grid::Vector<VectorComplexType> rreord[nscf];
    
    //Note:
    //il is the index of l, 
    //ir is the row index of M, 
    //jl is the column index of M and 
    //jr is the index of r

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

	  lreord[scf].resize(ni_this);
    	  for(int i = 0; i < ni_this; i++){
    	    i_ind.getBothIndices(ilmap_this[i],irmap_this[i],i,ilp,irp);
	    rowidx_used[ irmap_this[i] ] = true; //this row index of Mr is used
	    
	    const VectorComplexType &lval_tmp = l.nativeElem(ilmap_this[i], site4dop, sc, f);
	    lreord[scf][i] = conj_l ? Grid::conjugate(lval_tmp) : lval_tmp;
    	  }

    	  //j index
    	  int nj_this = j_ind.getNindices(jlp,jrp);
    	  nj[scf] = nj_this;
	  
    	  std::vector<int> &jlmap_this = jlmap[scf]; jlmap_this.resize(nj_this);
    	  std::vector<int> &jrmap_this = jrmap[scf]; jrmap_this.resize(nj_this);

	  rreord[scf].resize(nj_this);
    	  for(int j = 0; j < nj_this; j++){
    	    j_ind.getBothIndices(jlmap_this[j],jrmap_this[j],j,jlp,jrp);

	    const VectorComplexType &rval_tmp = r.nativeElem(jrmap_this[j], site4dop, sc, f);
	    rreord[scf][j] = conj_r ? Grid::conjugate(rval_tmp) : rval_tmp;
    	  }

     	}
      }
    }	  


    //Matrix vector multiplication  M*r contracted on mode index j. Only do it for rows that are actually used
    VectorComplexType Mr[Mrows][nscf];
    VectorComplexType tmp_v;
    for(int i=0;i<Mrows;i++){
      if(!rowidx_used[i]) continue;
      
      for(int scf=0;scf<nscf;scf++){
	int sc = scf/2; int f = scf % 2;
	int nj_this = nj[scf];
	zeroit(Mr[i][scf]);
	
	for(int j=0;j<nj_this;j++){
	  Grid::vsplat(tmp_v, M(i, jlmap[scf][j]) );
	  //const VectorComplexType &relem = r.nativeElem(jrmap[scf][j], site4dop, sc, f);
	  //Mr[i][scf] = Mr[i][scf] + tmp_v * (conj_r ? Grid::conjugate(relem) : relem);
	  
	  Mr[i][scf] = Mr[i][scf] + tmp_v * rreord[scf][j];
	}
      }
    }

    //Vector vector multiplication l*(M*r)
    for(int sl=0;sl<4;sl++){
      for(int sr=0;sr<4;sr++){
	for(int cl=0;cl<3;cl++){
	  int scl = cl + 3*sl;
	  for(int cr=0;cr<3;cr++){
	    for(int fl=0;fl<2;fl++){
	      int scfl = fl + 2*scl;
	      int ni_this = ni[scfl];
	      for(int fr=0;fr<2;fr++){
		int scfr = fr + 2*(cr + 3*sr);		

		VectorComplexType &into = out(sl,sr)(cl,cr)(fl,fr);

		for(int i=0;i<ni_this;i++){
		  //const VectorComplexType lelem = l.nativeElem(ilmap[scfl][i], site4dop, scl, fl);
		  //into = into + (conj_l ? Grid::conjugate(lelem) : lelem) * Mr[irmap[scfl][i]][scfr];
		  
		  into = into + lreord[scfl][i]*Mr[irmap[scfl][i]][scfr];
		}
	      }
	    }
	  }
	}
      }
    }	    

  }




  static void mult_slow(CPSspinColorFlavorMatrix<VectorComplexType> &out, const lA2AfieldL<mf_Policies> &l,  const A2AmesonField<mf_Policies,lA2AfieldR,rA2AfieldL> &M, const rA2AfieldR<mf_Policies> &r, const int &xop, const int &top, const bool &conj_l, const bool &conj_r){
    assert(l.getMode(0).SIMDlogicalNodes(3) == 1);
    
    int site4dop = l.getMode(0).threeToFour(xop,top);

    A2Aparams i_params(l), j_params(r);
    StandardIndexDilution idil(i_params), jdil(j_params);
    
    int ni = idil.getNmodes();
    int nj = jdil.getNmodes();

    out.zero();

    VectorComplexType tmp;
    
    for(int sl=0;sl<4;sl++){
      for(int cl=0;cl<3;cl++){
    	for(int fl=0;fl<2;fl++){	  

	  for(int sr=0;sr<4;sr++){
	    for(int cr=0;cr<3;cr++){
	      for(int fr=0;fr<2;fr++){

		for(int i=0;i<ni;i++){

		  const VectorComplexType &lval_tmp = l.elem(i,xop,top,cl+3*sl,fl);
		  VectorComplexType lval = conj_l ? Grid::conjugate(lval_tmp) : lval_tmp;
		  
  		  for(int j=0;j<nj;j++){
  		    const VectorComplexType &rval_tmp = r.elem(j,xop,top,cr+3*sr,fr);
  		    VectorComplexType rval = conj_r ? Grid::conjugate(rval_tmp) : rval_tmp;

		    const ScalarComplexType &Mval = M.elem(i,j);
		    Grid::vsplat(tmp, Mval);
		    
		    VectorComplexType delta = lval * tmp * rval;
  		    out(sl,sr)(cl,cr)(fl,fr) = out(sl,sr)(cl,cr)(fl,fr) + delta;
  		  }
  		}
  	      }
  	    }
  	  }
  	}
      }
    }

  }
};

// l^i(xop,top) M^ij r^j(xop,top)
template<typename mf_Policies, 
	 template <typename> class lA2Afield,  
	 template <typename> class MA2AfieldL,  template <typename> class MA2AfieldR,
	 template <typename> class rA2Afield  
	 >
void mult(CPSspinColorFlavorMatrix<typename mf_Policies::ComplexType> &out, const lA2Afield<mf_Policies> &l,  const A2AmesonField<mf_Policies,MA2AfieldL,MA2AfieldR> &M, const rA2Afield<mf_Policies> &r, const int &xop, const int &top, const bool &conj_l, const bool &conj_r){
  _mult_vMv_impl_v<mf_Policies,lA2Afield,MA2AfieldL,MA2AfieldR,rA2Afield, typename ComplexClassify<typename mf_Policies::ComplexType>::type >::mult(out,l,M,r,xop,top,conj_l,conj_r);
}
template<typename mf_Policies, 
	 template <typename> class lA2Afield,  
	 template <typename> class MA2AfieldL,  template <typename> class MA2AfieldR,
	 template <typename> class rA2Afield  
	 >
void mult_slow(CPSspinColorFlavorMatrix<typename mf_Policies::ComplexType> &out, const lA2Afield<mf_Policies> &l,  const A2AmesonField<mf_Policies,MA2AfieldL,MA2AfieldR> &M, const rA2Afield<mf_Policies> &r, const int &xop, const int &top, const bool &conj_l, const bool &conj_r){
  _mult_vMv_impl_v<mf_Policies,lA2Afield,MA2AfieldL,MA2AfieldR,rA2Afield, typename ComplexClassify<typename mf_Policies::ComplexType>::type >::mult_slow(out,l,M,r,xop,top,conj_l,conj_r);
}


#endif



#endif
