#ifndef _INNER_PRODUCT_GRID
#define _INNER_PRODUCT_GRID

//Grid implementations of performance critical components of inner_mult
#ifdef USE_GRID
#include<alg/a2a/grid_simdlinalg.h>
CPS_START_NAMESPACE

template<typename mf_Complex, bool conj_left, bool conj_right>
struct grid_mul{};

template<typename mf_Complex>
struct grid_mul<mf_Complex,true,false>{
  typedef typename Grid_A2A::SIMDvtype< mf_Complex >::Vtype Vtype;
  inline static void doit(Vtype &into, const Vtype &l, const Vtype &r){
    into += Grid::conjugate(l) * r; 
  }
};
template<typename mf_Complex>
struct grid_mul<mf_Complex,false,true>{
  typedef typename Grid_A2A::SIMDvtype< mf_Complex >::Vtype Vtype;
  inline static void doit(Vtype &into, const Vtype &l, const Vtype &r){
    into += l * Grid::conjugate(r); 
  }
};
template<typename mf_Complex>
struct grid_mul<mf_Complex,true,true>{
  typedef typename Grid_A2A::SIMDvtype< mf_Complex >::Vtype Vtype;
  inline static void doit(Vtype &into, const Vtype &l, const Vtype &r){
    into += Grid::conjugate(l * r); 
  }
};
template<typename mf_Complex>
struct grid_mul<mf_Complex,false,false>{
  typedef typename Grid_A2A::SIMDvtype< mf_Complex >::Vtype Vtype;
  inline static void doit(Vtype &into, const Vtype &l, const Vtype &r){
    into += l * r; 
  }
};

template<typename mf_Complex, bool conj_left, bool conj_right>
struct grid_g5contract{

  inline static void doit(std::complex<double> &v3, const mf_Complex *const l, const mf_Complex *const r){
    const static int sc_size =12;
    const static int half_sc = 6;
    typedef typename Grid_A2A::SIMDvtype<mf_Complex>::Vtype Vtype;
    const int nsimd = Vtype::Nsimd();
    const int overspill_amnt = sc_size % nsimd; //how many elements don't fit in nice SIMD blocks
    const int mid_block = sc_size / 2 / nsimd; //index of block corresponding to middle of the 12 elements
    const int mid_subidx = sc_size / 2 - mid_block * nsimd;
    const int nblocks = sc_size / nsimd;

    Vtype vl, vr, sol;
    zeroit(sol);
    mf_Complex *il = const_cast<mf_Complex *>(l);
    mf_Complex *ir = const_cast<mf_Complex *>(r);

    //First set of blocks
    for(int i=0;i<mid_block;i++){
      vset(vl,il);
      vset(vr,ir);
      il += nsimd;
      ir += nsimd;
      grid_mul<mf_Complex,conj_left,conj_right>::doit(sol,vl,vr);
    }

    //Middle block
    mf_Complex block[nsimd];
    {
      for(int i=0;i<mid_subidx;i++)
	block[i] = *(il++);
      for(int i=mid_subidx;i<nsimd;i++)
	block[i] = -*(il++);
      vset(vl,block);
      vset(vr,ir);
      ir += nsimd;
      grid_mul<mf_Complex,conj_left,conj_right>::doit(sol,vl,vr);
    }

    //Remaining packed blocks
    for(int i=mid_block+1;i<nblocks;i++){
      vset(vl,il);
      vl = -vl;
      vset(vr,ir);
      il += nsimd;
      ir += nsimd;
      grid_mul<mf_Complex,conj_left,conj_right>::doit(sol,vl,vr);
    }

    //Overspill
    if(overspill_amnt > 0){
      for(int i=overspill_amnt;i<nsimd;i++) block[i] = 0.;

      for(int i=0;i<overspill_amnt;i++) 
	block[i] = il[i];
      vset(vl,block);
      vl = -vl;
      for(int i=0;i<overspill_amnt;i++) 
	block[i] = ir[i];      
      vset(vr,block);
      grid_mul<mf_Complex,conj_left,conj_right>::doit(sol,vl,vr);
    }
    v3 = Reduce(sol);
  }
};




//Version of the above for G-parity that produces a flavor matrix reusing vectors where possible
template<typename mf_Complex, bool conj_left, bool conj_right>
struct grid_scf_contract{
  typedef typename Grid_A2A::SIMDvtype<mf_Complex>::Vtype Vtype;
  
  inline static void grid_g5_forml(Vtype l_f[], const int f, const SCFvectorPtr<mf_Complex> &l, 
				   mf_Complex block[], const int nsimd, const int overspill_amnt, const int mid_block, const int mid_subidx, const int nblocks ){
    mf_Complex *il = const_cast<mf_Complex*>(l.getPtr(f));
    if(nsimd == 4){ //avx256 with complex<float> or avx512 with complex<double>
      vset(l_f[0],il); //first 4 complexes
      block[0] = il[4]; //second 4
      block[1] = il[5];
      block[2] = -il[6];
      block[3] = -il[7];
      vset(l_f[1],block);
      vset(l_f[2],il+8); //last 4
      l_f[2] = -l_f[2];
    }else{
      //First block
      for(int i=0;i<mid_block;i++){
	vset(l_f[i],il);
	il += nsimd;
      }
      //Mid block
      for(int i=0;i<mid_subidx;i++)
	block[i] = *(il++);
      for(int i=mid_subidx;i<nsimd;i++)
	block[i] = -*(il++);
      vset(l_f[mid_block],block);
      
      //Remaining packed blocks
      for(int i=mid_block+1;i<nblocks;i++){
	vset(l_f[i],il);
	l_f[i] = -l_f[i];
	il += nsimd;
      }
      
      //Overspill
      if(overspill_amnt > 0){
	for(int i=0;i<overspill_amnt;i++) 
	  block[i] = -il[i];
	for(int i=overspill_amnt;i<nsimd;i++) block[i] = 0.;
	vset(l_f[nblocks],block);
      }
    }
  }
  inline static void grid_g5_formr(Vtype r_f[], const int f, const SCFvectorPtr<mf_Complex> &r, 
				   mf_Complex block[], const int nsimd, const int overspill_amnt, const int nblocks ){
    mf_Complex *ir = const_cast<mf_Complex*>(r.getPtr(f));
      
    //Packed blocks
    for(int i=0;i<nblocks;i++){
      vset(r_f[i],ir);
      ir += nsimd;
    }
    //Overspill
    if(overspill_amnt > 0){
      for(int i=0;i<overspill_amnt;i++) 
	block[i] = ir[i];
      for(int i=overspill_amnt;i<nsimd;i++) block[i] = 0.;
      vset(r_f[nblocks],block);
    }
  }

  static void grid_g5con(FlavorMatrix &lMr, const SCFvectorPtr<mf_Complex> &l, const SCFvectorPtr<mf_Complex> &r){
    const int sc_size = 12;
    const int nsimd = Vtype::Nsimd();
    const int overspill_amnt = sc_size % nsimd; //how many elements don't fit in nice SIMD blocks
    const int mid_block = sc_size / 2 / nsimd; //index of block corresponding to middle of the 12 elements
    const int mid_subidx = sc_size / 2 - mid_block * nsimd;
    const int nblocks = sc_size / nsimd;
    const int nblocks_inc_overspill = nblocks + (overspill_amnt > 0 ? 1:0);
    
    Vtype l_0[nblocks_inc_overspill];
    Vtype l_1[nblocks_inc_overspill];
    Vtype r_0[nblocks_inc_overspill];
    Vtype r_1[nblocks_inc_overspill];

    mf_Complex block[nsimd];

    if(!l.isZero(0))
      grid_g5_forml(l_0,0,l, block, nsimd, overspill_amnt, mid_block, mid_subidx, nblocks );
    if(!l.isZero(1))
      grid_g5_forml(l_1,1,l, block, nsimd, overspill_amnt, mid_block, mid_subidx, nblocks );
    if(!r.isZero(0))
      grid_g5_formr(r_0,0,r, block, nsimd, overspill_amnt, nblocks );
    if(!r.isZero(1))
      grid_g5_formr(r_1,1,r, block, nsimd, overspill_amnt, nblocks );    
      
    Vtype sol;
    //0,0
    zeroit(sol);
    if(!l.isZero(0) && !r.isZero(0)){
      for(int i=0;i<nblocks_inc_overspill;i++)
	grid_mul<mf_Complex,conj_left,conj_right>::doit(sol,l_0[i],r_0[i]);
      lMr(0,0) = Reduce(sol);
    }
    //0,1
    zeroit(sol);
    if(!l.isZero(0) && !r.isZero(1)){
      for(int i=0;i<nblocks_inc_overspill;i++)
	grid_mul<mf_Complex,conj_left,conj_right>::doit(sol,l_0[i],r_1[i]);
      lMr(0,1) = Reduce(sol);
    }
    //1,0
    zeroit(sol);
    if(!l.isZero(1) && !r.isZero(0)){
      for(int i=0;i<nblocks_inc_overspill;i++)
	grid_mul<mf_Complex,conj_left,conj_right>::doit(sol,l_1[i],r_0[i]);
      lMr(1,0) = Reduce(sol);
    }
    //1,1
    zeroit(sol);
    if(!l.isZero(1) && !r.isZero(1)){
      for(int i=0;i<nblocks_inc_overspill;i++)
	grid_mul<mf_Complex,conj_left,conj_right>::doit(sol,l_1[i],r_1[i]);
      lMr(1,1) = Reduce(sol);
    }
  }


  static void grid_unitcon(FlavorMatrix &lMr, const SCFvectorPtr<mf_Complex> &l, const SCFvectorPtr<mf_Complex> &r){
    const int sc_size = 12;
    const int nsimd = Vtype::Nsimd();
    const int overspill_amnt = sc_size % nsimd; //how many elements don't fit in nice SIMD blocks
    const int mid_block = sc_size / 2 / nsimd; //index of block corresponding to middle of the 12 elements
    const int mid_subidx = sc_size / 2 - mid_block * nsimd;
    const int nblocks = sc_size / nsimd;
    const int nblocks_inc_overspill = nblocks + (overspill_amnt > 0 ? 1:0);
    
    Vtype l_0[nblocks_inc_overspill];
    Vtype l_1[nblocks_inc_overspill];
    Vtype r_0[nblocks_inc_overspill];
    Vtype r_1[nblocks_inc_overspill];

    mf_Complex block[nsimd];

    if(!l.isZero(0))
      grid_g5_formr(l_0,0,l, block, nsimd, overspill_amnt, nblocks );
    if(!l.isZero(1))
      grid_g5_formr(l_1,1,l, block, nsimd, overspill_amnt, nblocks );
    if(!r.isZero(0))
      grid_g5_formr(r_0,0,r, block, nsimd, overspill_amnt, nblocks );
    if(!r.isZero(1))
      grid_g5_formr(r_1,1,r, block, nsimd, overspill_amnt, nblocks );    
      
    Vtype sol;
    //0,0
    zeroit(sol);
    if(!l.isZero(0) && !r.isZero(0)){
      for(int i=0;i<nblocks_inc_overspill;i++)
	grid_mul<mf_Complex,conj_left,conj_right>::doit(sol,l_0[i],r_0[i]);
      lMr(0,0) = Reduce(sol);
    }
    //0,1
    zeroit(sol);
    if(!l.isZero(0) && !r.isZero(1)){
      for(int i=0;i<nblocks_inc_overspill;i++)
	grid_mul<mf_Complex,conj_left,conj_right>::doit(sol,l_0[i],r_1[i]);
      lMr(0,1) = Reduce(sol);
    }
    //1,0
    zeroit(sol);
    if(!l.isZero(1) && !r.isZero(0)){
      for(int i=0;i<nblocks_inc_overspill;i++)
	grid_mul<mf_Complex,conj_left,conj_right>::doit(sol,l_1[i],r_0[i]);
      lMr(1,0) = Reduce(sol);
    }
    //1,1
    zeroit(sol);
    if(!l.isZero(1) && !r.isZero(1)){
      for(int i=0;i<nblocks_inc_overspill;i++)
	grid_mul<mf_Complex,conj_left,conj_right>::doit(sol,l_1[i],r_1[i]);
      lMr(1,1) = Reduce(sol);
    }
  }
};



template<int smatidx,typename mf_Complex, bool conj_left, bool conj_right>
struct grid_scf_contract_select{};

template<typename mf_Complex, bool conj_left, bool conj_right>
struct grid_scf_contract_select<15,mf_Complex,conj_left,conj_right>{
  inline static void doit(FlavorMatrix &lMr, const SCFvectorPtr<mf_Complex> &l, const SCFvectorPtr<mf_Complex> &r){
    grid_scf_contract<mf_Complex,conj_left,conj_right>::grid_g5con(lMr,l,r);
  }
};
template<typename mf_Complex, bool conj_left, bool conj_right>
struct grid_scf_contract_select<0,mf_Complex,conj_left,conj_right>{
  inline static void doit(FlavorMatrix &lMr, const SCFvectorPtr<mf_Complex> &l, const SCFvectorPtr<mf_Complex> &r){
    grid_scf_contract<mf_Complex,conj_left,conj_right>::grid_unitcon(lMr,l,r);
  }
};





CPS_END_NAMESPACE
#endif
#endif

