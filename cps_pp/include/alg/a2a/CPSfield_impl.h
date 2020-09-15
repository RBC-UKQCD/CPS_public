#ifndef CPS_FIELD_IMPL
#define CPS_FIELD_IMPL

//Implementations of CPSfield.h


//Real-reduce for norm2
template<typename T>
struct normdefs{};

template<>
struct normdefs<double>{
  inline static double real_reduce(const double in){ return in; }
  inline static double conjugate(const double in){ return in; }
};
template<>
struct normdefs<float>{
  inline static double real_reduce(const float in){ return in; }
  inline static float conjugate(const float in){ return in; }
};
template<typename T>
struct normdefs<std::complex<T> >{
  inline static double real_reduce(const std::complex<T> in){ return in.real(); }
  inline static std::complex<T> conjugate(const std::complex<T> in){ return std::conj(in); }
};
#ifdef USE_GRID
template<>
struct normdefs<Grid::vRealD>{
  inline static double real_reduce(const Grid::vRealD in){ return Reduce(in); }
  inline static Grid::vRealD conjugate(const Grid::vRealD in){ return in; }
};
template<>
struct normdefs<Grid::vRealF>{
  inline static double real_reduce(const Grid::vRealF in){ return Reduce(in); }
  inline static Grid::vRealF conjugate(const Grid::vRealF in){ return in; }
};
template<>
struct normdefs<Grid::vComplexD>{
  inline static double real_reduce(const Grid::vComplexD in){ return std::real(Reduce(in)); }
  inline static Grid::vComplexD conjugate(const Grid::vComplexD in){ return Grid::conjugate(in); }
};
template<>
struct normdefs<Grid::vComplexF>{
  inline static double real_reduce(const Grid::vComplexF in){ return std::real(Reduce(in)); }
  inline static Grid::vComplexF conjugate(const Grid::vComplexF in){ return Grid::conjugate(in); }
};
#endif

template< typename SiteType, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
double CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy>::norm2() const{
  SiteType accum[omp_get_max_threads()];
  memset(accum, 0, omp_get_max_threads()*sizeof(SiteType));
#pragma omp parallel for  
  for(int i=0;i<this->nfsites();i++){
    SiteType const *site = this->fsite_ptr(i);
    for(int s=0;s<SiteSize;s++)
      accum[omp_get_thread_num()] = accum[omp_get_thread_num()] + normdefs<SiteType>::conjugate(site[s])*site[s];
  }
  SiteType total;
  memset(&total, 0, sizeof(SiteType));

  for(int i=0;i<omp_get_max_threads();i++)
    total = total + accum[i];

  double final = normdefs<SiteType>::real_reduce(total);
  glb_sum_five(&final);
  return final;
}




#ifdef USE_GRID

template<typename T,typename CPScomplex>
struct GridTensorConvert{};

template<typename complex_scalar, typename CPScomplex>
struct GridTensorConvert<Grid::QCD::iSpinColourVector<complex_scalar>, CPScomplex>{
  static_assert(!Grid::isSIMDvectorized<complex_scalar>::value && Grid::isComplex<complex_scalar>::value, "Only applies to scalar complex types");

  //12-component complex spin-color vector
  //We have assured the input is not SIMD vectorized so the output type is the same
  inline static void doit(CPScomplex* cps, const Grid::QCD::iSpinColourVector<complex_scalar> &grid, const int f){
    for(int s=0;s<Grid::QCD::Ns;s++)
      for(int c=0;c<Grid::QCD::Nc;c++)
	*cps++ = grid()(s)(c);
  }
  inline static void doit(Grid::QCD::iSpinColourVector<complex_scalar> &grid, CPScomplex const* cps, const int f){
    for(int s=0;s<Grid::QCD::Ns;s++)
      for(int c=0;c<Grid::QCD::Nc;c++)
	grid()(s)(c) = *cps++;
  }
};
template<typename complex_scalar, typename CPScomplex>
struct GridTensorConvert<Grid::QCD::iGparitySpinColourVector<complex_scalar>, CPScomplex>{
  static_assert(!Grid::isSIMDvectorized<complex_scalar>::value && Grid::isComplex<complex_scalar>::value, "Only applies to scalar complex types");

  //12-component complex spin-color vector
  //We have assured the input is not SIMD vectorized so the output type is the same
  inline static void doit(CPScomplex* cps, const Grid::QCD::iGparitySpinColourVector<complex_scalar> &grid, const int f){
    for(int s=0;s<Grid::QCD::Ns;s++)
      for(int c=0;c<Grid::QCD::Nc;c++)
	*cps++ = grid(f)(s)(c);
  }
  inline static void doit(Grid::QCD::iGparitySpinColourVector<complex_scalar> &grid, CPScomplex const* cps, const int f){
    for(int s=0;s<Grid::QCD::Ns;s++)
      for(int c=0;c<Grid::QCD::Nc;c++)
  	grid(f)(s)(c) = *cps++;
  }
};
template<typename complex_scalar, typename CPScomplex>
struct GridTensorConvert<Grid::QCD::iLorentzColourMatrix<complex_scalar>, CPScomplex>{
  static_assert(!Grid::isSIMDvectorized<complex_scalar>::value && Grid::isComplex<complex_scalar>::value, "Only applies to scalar complex types");

  //Gauge field  mu=0..3  3*3 complex
  //We have assured the input is not SIMD vectorized so the output type is the same
  inline static void doit(CPScomplex* cps, const Grid::QCD::iLorentzColourMatrix<complex_scalar> &grid, const int f){
    for(int mu=0;mu<4;mu++)
      for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
	  *cps++ = grid(mu)()(i,j);
  }
  inline static void doit(Grid::QCD::iLorentzColourMatrix<complex_scalar> &grid, CPScomplex const* cps, const int f){
    for(int mu=0;mu<4;mu++)
      for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
	  grid(mu)()(i,j) = *cps++;
  }
};


template<int Ndim>
struct dimensionMap{};

template<>
struct dimensionMap<5>{
  const int cps_to_grid[5] = {1,2,3,4,0};
  const int grid_to_cps[5] = {4,0,1,2,3};
};
template<>
struct dimensionMap<4>{
  const int cps_to_grid[4] = {0,1,2,3};
  const int grid_to_cps[4] = {0,1,2,3};
};

template<typename Type, int SiteSize, typename DimPol, typename FlavPol, typename AllocPol,
	 typename GridField, typename ComplexClass>
class CPSfieldGridImpex{};

template<typename Type, int SiteSize, typename DimPol, typename FlavPol, typename AllocPol,
	 typename GridField>
class CPSfieldGridImpex<Type,SiteSize,DimPol,FlavPol,AllocPol,GridField,complex_double_or_float_mark>{
  typedef CPSfield<Type,SiteSize,DimPol,FlavPol,AllocPol> CPSfieldType;

public:

  static void import(CPSfieldType &into, const GridField &from){
    const int Nd = DimPol::EuclideanDimension;
    assert(Nd == from._grid->Nd());
    dimensionMap<CPSfieldType::EuclideanDimension> dim_map;

    typedef typename Grid::GridTypeMapper<typename GridField::vector_object>::scalar_object sobj;
#pragma omp parallel for
    for(int site=0;site<into.nsites();site++){
      std::vector<int> x(Nd);
      into.siteUnmap(site, &x[0]);

      std::vector<int> grid_x(Nd);
      for(int i=0;i<Nd;i++)
	grid_x[ dim_map.cps_to_grid[i] ] = x[i];

      sobj siteGrid; //contains both flavors if Gparity
      peekLocalSite(siteGrid,from,grid_x);

      for(int f=0;f<into.nflavors();f++){
	typename CPSfieldType::FieldSiteType *cps = into.site_ptr(site,f);
	GridTensorConvert<sobj, typename CPSfieldType::FieldSiteType>::doit(cps, siteGrid, f);
      }      
    }
  }
  
  static void exportit(GridField &into, const CPSfieldType &from){
    const int Nd = DimPol::EuclideanDimension;
    assert(Nd == into._grid->Nd());
    dimensionMap<CPSfieldType::EuclideanDimension> dim_map;
  
    typedef typename Grid::GridTypeMapper<typename GridField::vector_object>::scalar_object sobj;
    int nthread = omp_get_max_threads();
    
    int nsimd = into._grid->Nsimd();
    std::vector<std::vector<sobj> > tstore(nthread,std::vector<sobj>(nsimd)); //thread-individual temp storage for Grid-converted tensors
    std::vector<std::vector<sobj*> > tstore_ptrs(nthread,std::vector<sobj*>(nsimd));
    for(int i=0;i<nthread;i++)
      for(int j=0;j<nsimd;j++)
	tstore_ptrs[i][j] = &tstore[i][j];

    std::vector<std::vector<int> > out_icoor(nsimd); //store inner coordinate offsets
    for(int i=0;i<nsimd;i++){
      out_icoor[i].resize(Nd);
      into._grid->iCoorFromIindex(out_icoor[i], i);
    }
#pragma omp parallel for
    for(int out_oidx=0;out_oidx<into._grid->oSites();out_oidx++){
      int me = omp_get_thread_num();
      std::vector<int> out_ocoor(Nd);
      into._grid->oCoorFromOindex(out_ocoor, out_oidx);

      std::vector<int> lcoor(Nd);
      std::vector<int> lcoor_cps(Nd);
      
      for(int lane=0; lane < nsimd; lane++){
	for(int mu=0;mu<Nd;mu++){
	  lcoor[mu] = out_ocoor[mu] + into._grid->_rdimensions[mu]*out_icoor[lane][mu];
	  lcoor_cps[ dim_map.grid_to_cps[mu] ] = lcoor[mu];
	}
	int cps_site = from.siteMap(&lcoor_cps[0]);
	for(int f=0;f<from.nflavors();f++){
	  typename CPSfieldType::FieldSiteType const* cps = from.site_ptr(cps_site,f);
	  GridTensorConvert<sobj, typename CPSfieldType::FieldSiteType>::doit(tstore[me][lane], cps, f);
	}
      }
      merge(into._odata[out_oidx], tstore_ptrs[me], 0);
    }
  }
  
};

#ifdef USE_GRID
template<typename Type, int SiteSize, typename DimPol, typename FlavPol, typename AllocPol,
	 typename GridField>
class CPSfieldGridImpex<Type,SiteSize,DimPol,FlavPol,AllocPol,GridField,grid_vector_complex_mark>{
  typedef CPSfield<Type,SiteSize,DimPol,FlavPol,AllocPol> CPSfieldType;

public:

  static void import(CPSfieldType &into, const GridField &from){
    const int Nd = DimPol::EuclideanDimension;
    assert(Nd == from._grid->Nd());
    typedef typename Grid::GridTypeMapper<Type>::scalar_type CPSscalarType;
    typedef typename ComplexClassify<CPSscalarType>::type CPSscalarTypeClass;
    
    //Create temp CPS unvectorized field
    typedef typename StandardDimensionPolicy<DimPol::EuclideanDimension>::type CPSscalarDimPol;
    NullObject n;
    CPSfield<CPSscalarType,SiteSize,CPSscalarDimPol,FlavPol,StandardAllocPolicy> cps_unpacked(n);

    CPSfieldGridImpex<CPSscalarType,SiteSize,CPSscalarDimPol,FlavPol,StandardAllocPolicy,GridField, CPSscalarTypeClass>::import(cps_unpacked,from);
    into.importField(cps_unpacked);
  }
  
  static void exportit(GridField &into, const CPSfieldType &from){
    const int Nd = DimPol::EuclideanDimension;
    assert(Nd == into._grid->Nd());
    typedef typename Grid::GridTypeMapper<Type>::scalar_type CPSscalarType;
    typedef typename ComplexClassify<CPSscalarType>::type CPSscalarTypeClass;

    //Create temp CPS unvectorized field
    typedef typename StandardDimensionPolicy<DimPol::EuclideanDimension>::type CPSscalarDimPol;
    NullObject n;
    CPSfield<CPSscalarType,SiteSize,CPSscalarDimPol,FlavPol,StandardAllocPolicy> cps_unpacked(n);
    cps_unpacked.importField(from);
    CPSfieldGridImpex<CPSscalarType,SiteSize,CPSscalarDimPol,FlavPol,StandardAllocPolicy,GridField, CPSscalarTypeClass>::exportit(into, cps_unpacked);
  }
};
#endif





template< typename SiteType, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
template<typename GridField>
void  CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy>::importGridField(const GridField &grid){
  typedef typename ComplexClassify<SiteType>::type ComplexClass;
  CPSfieldGridImpex<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy,GridField,ComplexClass>::import(*this, grid);
}
template< typename SiteType, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
template<typename GridField>
void  CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy>::exportGridField(GridField &grid) const{
  typedef typename ComplexClassify<SiteType>::type ComplexClass;
  CPSfieldGridImpex<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy,GridField,ComplexClass>::exportit(grid,*this);
}
#endif

template<int SiteSize,
	 typename TypeA, typename DimPolA, typename FlavPolA, typename AllocPolA,
	 typename TypeB, typename DimPolB, typename FlavPolB, typename AllocPolB,
	 typename Enable = void>
class CPSfieldCopy;

//Generic copy. SiteSize and number of Euclidean dimensions must be the same
#ifdef USE_GRID
#define CONDITION sameDim<DimPolA,DimPolB>::val && !Grid::is_simd<TypeA>::value && !Grid::is_simd<TypeB>::value
#else
#define CONDITION sameDim<DimPolA,DimPolB>::val
#endif

template<int SiteSize,
	 typename TypeA, typename DimPolA, typename FlavPolA, typename AllocPolA,
	 typename TypeB, typename DimPolB, typename FlavPolB, typename AllocPolB>
class CPSfieldCopy<SiteSize,TypeA,DimPolA,FlavPolA,AllocPolA, TypeB,DimPolB,FlavPolB,AllocPolB, typename my_enable_if<CONDITION,void>::type>{
public: 
  static void copy(CPSfield<TypeA,SiteSize,DimPolA,FlavPolA,AllocPolA> &into,
		   const CPSfield<TypeB,SiteSize,DimPolB,FlavPolB,AllocPolB> &from){
    assert(into.nfsites() == from.nfsites()); //should be true in # Euclidean dimensions the same, but not guaranteed
    
    #pragma omp parallel for
    for(int fs=0;fs<into.nfsites();fs++){
      int x[5], f; into.fsiteUnmap(fs,x,f); //doesn't matter if the linearization differs between the two
      TypeA* toptr = into.fsite_ptr(fs);
      TypeB const* fromptr = from.site_ptr(x,f);
      for(int i=0;i<SiteSize;i++) toptr[i] = fromptr[i];
    }
  }
};
#undef CONDITION

#ifdef USE_GRID

std::string vtostring(const int* v, const int ndim){
  std::ostringstream os;
  os << '(';
  for(int i=0;i<ndim-1;i++) os << v[i] << ", ";
  os << v[ndim-1] << ')';
  return os.str();
}

//TypeA is Grid_simd type
#define CONDITION sameDim<DimPolA,DimPolB>::val && Grid::is_simd<GridSIMDTypeA>::value && !Grid::is_simd<TypeB>::value

template<int SiteSize,
	 typename GridSIMDTypeA, typename DimPolA, typename FlavPolA, typename AllocPolA,
	 typename TypeB, typename DimPolB, typename FlavPolB, typename AllocPolB>
class CPSfieldCopy<SiteSize,
		   GridSIMDTypeA, DimPolA, FlavPolA, AllocPolA,
		   TypeB, DimPolB, FlavPolB, AllocPolB, typename my_enable_if<CONDITION,void>::type>
{
public:
  static void copy(CPSfield<GridSIMDTypeA,SiteSize,DimPolA,FlavPolA,AllocPolA> &into,
		   const CPSfield<TypeB,SiteSize,DimPolB,FlavPolB,AllocPolB> &from){
    const int nsimd = GridSIMDTypeA::Nsimd();
    const int ndim = DimPolA::EuclideanDimension;
    if(from.nfsites()/nsimd != into.nfsites()) ERR.General("CPSfieldCopy","copy(<SIMD field> &into, const <non-SIMD field> &from)","Expected from.nfsites/nsimd = into.nfsites, got %d/%d (=%d) != %d\n",from.nfsites(),nsimd, from.nfsites()/nsimd, into.nfsites());
    
    std::vector<std::vector<int> > packed_offsets(nsimd,std::vector<int>(ndim));
    for(int i=0;i<nsimd;i++){
      into.SIMDunmap(i,&packed_offsets[i][0]);
    }
    
#pragma omp parallel for
    for(int fs=0;fs<into.nfsites();fs++){
      int x[ndim], f; into.fsiteUnmap(fs,x,f);
      GridSIMDTypeA* toptr = into.fsite_ptr(fs);

      //x is the root coordinate corresponding to SIMD packed index 0      
      std::vector<TypeB const*> ptrs(nsimd);
      ptrs[0] = from.site_ptr(x,f);
      
      int xx[ndim];
      for(int i=1;i<nsimd;i++){
	for(int d=0;d<ndim;d++)
	  xx[d] = x[d] + packed_offsets[i][d];  //xx = x + offset
	ptrs[i] = from.site_ptr(xx,f);
      }
      into.SIMDpack(toptr, ptrs, SiteSize);
    }
  }
};
#undef CONDITION

//TypeB is Grid_simd type
#define CONDITION sameDim<DimPolA,DimPolB>::val && !Grid::is_simd<TypeA>::value && Grid::is_simd<GridSIMDTypeB>::value

template<int SiteSize,
	 typename TypeA, typename DimPolA, typename FlavPolA, typename AllocPolA,
	 typename GridSIMDTypeB, typename DimPolB, typename FlavPolB, typename AllocPolB>
class CPSfieldCopy<SiteSize,
		   TypeA, DimPolA, FlavPolA, AllocPolA,
		   GridSIMDTypeB, DimPolB, FlavPolB, AllocPolB, typename my_enable_if<CONDITION,void>::type>
{
public:
  static void copy(CPSfield<TypeA,SiteSize,DimPolA,FlavPolA,AllocPolA> &into,
		   const CPSfield<GridSIMDTypeB,SiteSize,DimPolB,FlavPolB,AllocPolB> &from){
    const int nsimd = GridSIMDTypeB::Nsimd();
    const int ndim = DimPolA::EuclideanDimension;
    if(into.nfsites()/nsimd != from.nfsites()) ERR.General("CPSfieldCopy","copy(<non-SIMD field> &into, const <SIMD-field> &from)","Expected into.nfsites/nsimd = from.nfsites, got %d/%d (=%d) != %d\n",into.nfsites(),nsimd, into.nfsites()/nsimd, from.nfsites());

    std::vector<std::vector<int> > packed_offsets(nsimd,std::vector<int>(ndim));
    for(int i=0;i<nsimd;i++) from.SIMDunmap(i,&packed_offsets[i][0]);

#pragma omp parallel for
    for(int fs=0;fs<from.nfsites();fs++){
      int x[ndim], f; from.fsiteUnmap(fs,x,f);
      GridSIMDTypeB const* fromptr = from.fsite_ptr(fs);

      //x is the root coordinate corresponding to SIMD packed index 0
      std::vector<TypeA*> ptrs(nsimd);
      ptrs[0] = into.site_ptr(x,f);
      
      int xx[ndim];
      for(int i=1;i<nsimd;i++){
	for(int d=0;d<ndim;d++)
	  xx[d] = x[d] + packed_offsets[i][d];  //xx = x + offset

	ptrs[i] = into.site_ptr(xx,f);
      }
      from.SIMDunpack(ptrs, fromptr, SiteSize);
    }
  }
};
#undef CONDITION

#endif

template< typename SiteType, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
template< typename extSiteType, typename extDimPol, typename extFlavPol, typename extAllocPol>
void CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy>::importField(const CPSfield<extSiteType,SiteSize,extDimPol,extFlavPol,extAllocPol> &r){
  CPSfieldCopy<SiteSize,
	       SiteType,DimensionPolicy,FlavorPolicy,AllocPolicy,
	       extSiteType, extDimPol, extFlavPol, extAllocPol>::copy(*this,r);
}
template< typename SiteType, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
template< typename extSiteType, typename extDimPol, typename extFlavPol, typename extAllocPol>
void CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy>::exportField(CPSfield<extSiteType,SiteSize,extDimPol,extFlavPol,extAllocPol> &r) const{
  CPSfieldCopy<SiteSize,
	       extSiteType, extDimPol, extFlavPol, extAllocPol,
	       SiteType,DimensionPolicy,FlavorPolicy,AllocPolicy>::copy(r,*this);
}




template<typename SiteType>
class _testRandom{
public:
  static void rand(SiteType* f, int fsize, const Float hi, const Float lo){
    for(int i=0;i<fsize;i++) f[i] = LRG.Urand(hi,lo,FOUR_D);
  }
};
template<typename T>
class _testRandom<std::complex<T> >{
public:
  static void rand(std::complex<T>* f, int fsize, const Float hi, const Float lo){
    assert(sizeof(std::complex<T>) == 2*sizeof(T));
    T* ff = (T*)f;
    for(int i=0;i<2*fsize;i++) ff[i] = LRG.Urand(hi,lo,FOUR_D);
  }
};


//Set each float to a uniform random number in the specified range.
//WARNING: Uses only the current RNG in LRG, and does not change this based on site. This is therefore only useful for testing*
template< typename SiteType, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy>::testRandom(const Float hi, const Float lo){
  _testRandom<SiteType>::rand(this->f,this->fsize,hi,lo);
}
template< typename SiteType, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
//Set this field to the average of this and a second field, r
void CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy>::average(const CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &r, const bool &parallel){
  //The beauty of having the ordering baked into the policy class is that we implicitly *know* the ordering of the second field, so we can just loop over the floats in a dumb way
  if(parallel){
#pragma omp parallel for
    for(int i=0;i<fsize;i++) f[i] = (f[i] + r.f[i])/2.0;
  }else{
    for(int i=0;i<fsize;i++) f[i] = (f[i] + r.f[i])/2.0;
  }
}



struct _gauge_fix_site_op_impl{
  
  template< typename mf_Complex, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy, typename my_enable_if<_equal<typename ComplexClassify<mf_Complex>::type,complex_double_or_float_mark>::value,int>::type = 0>
  inline static void gauge_fix_site_op(CPSfermion<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy> &field, const int x4d[], const int &f, Lattice &lat, const bool dagger){
    typedef typename mf_Complex::value_type mf_Float;
    int i = x4d[0] + GJP.XnodeSites()*( x4d[1] + GJP.YnodeSites()* ( x4d[2] + GJP.ZnodeSites()*x4d[3] ) );
    mf_Complex tmp[3];
    const Matrix* gfmat = lat.FixGaugeMatrix(i,f);
    mf_Complex* sc_base = (mf_Complex*)field.site_ptr(x4d,f); //if Dimension < 4 the site_ptr method will ignore the remaining indices. Make sure this is what you want
    for(int s=0;s<4;s++){
      memcpy(tmp, sc_base + 3 * s, 3 * sizeof(mf_Complex));
      if(!dagger)
	colorMatrixMultiplyVector<mf_Float,Float>( (mf_Float*)(sc_base + 3*s), (Float*)gfmat, (mf_Float*)tmp);
      else
	colorMatrixDaggerMultiplyVector<mf_Float,Float>( (mf_Float*)(sc_base + 3*s), (Float*)gfmat, (mf_Float*)tmp);      
    }
  }
#ifdef USE_GRID
  template< typename mf_Complex, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy, typename my_enable_if<_equal<typename ComplexClassify<mf_Complex>::type,grid_vector_complex_mark>::value,int>::type = 0>
  inline static void gauge_fix_site_op(CPSfermion<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy> &field, const int x4d[], const int &f, Lattice &lat, const bool dagger){
    //x4d is an outer site index
    int nsimd = field.Nsimd();
    int ndim = DimensionPolicy::EuclideanDimension;
    assert(ndim == 4);

    //Assemble pointers to the GF matrices for each lane
    std::vector<cps::Complex*> gf_base_ptrs(nsimd);
    int x4d_lane[4];
    int lane_off[4];
    
    for(int lane=0;lane<nsimd;lane++){
      field.SIMDunmap(lane, lane_off);		      
      for(int xx=0;xx<4;xx++) x4d_lane[xx] = x4d[xx] + lane_off[xx];
      int gf_off = x4d_lane[0] + GJP.XnodeSites()*( x4d_lane[1] + GJP.YnodeSites()* ( x4d_lane[2] + GJP.ZnodeSites()*x4d_lane[3] ) );
      gf_base_ptrs[lane] = (cps::Complex*)lat.FixGaugeMatrix(gf_off,f);
    }


    //Poke the GFmatrix elements into SIMD vector objects
    typedef typename mf_Complex::scalar_type stype;
    stype* buf = (stype*)memalign(128, nsimd*sizeof(stype));

    mf_Complex gfmat[3][3];
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){

	for(int lane=0;lane<nsimd;lane++)
	  buf[lane] = *(gf_base_ptrs[lane] + j + 3*i);
	vset(gfmat[i][j], buf);
      }
    }

    free(buf);

    //Do the matrix multiplication
    mf_Complex* tmp = (mf_Complex*)memalign(128, 3*sizeof(mf_Complex));
    mf_Complex* sc_base = field.site_ptr(x4d,f); 
    for(int s=0;s<4;s++){
      mf_Complex* s_base = sc_base + 3 * s;
      memcpy(tmp, s_base, 3 * sizeof(mf_Complex));
      if(!dagger)
	for(int i=0;i<3;i++)
	  s_base[i] = gfmat[i][0]*tmp[0] + gfmat[i][1]*tmp[1] + gfmat[i][2]*tmp[2];
      else
	for(int i=0;i<3;i++)
	  s_base[i] = conjugate(gfmat[0][i])*tmp[0] + conjugate(gfmat[1][i])*tmp[1] + conjugate(gfmat[2][i])*tmp[2];
    }
    free(tmp);
  }
#endif
  
};



//Apply gauge fixing matrices to the field
template< typename mf_Complex, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void CPSfermion<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>::gauge_fix_site_op(const int x4d[], const int &f, Lattice &lat, const bool dagger){
  _gauge_fix_site_op_impl::gauge_fix_site_op(*this, x4d, f, lat,dagger);
}

template< typename mf_Complex, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void CPSfermion<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>::getMomentumUnits(double punits[3]){
  for(int i=0;i<3;i++){
    int fac;
    if(GJP.Bc(i) == BND_CND_PRD) fac = 1;
    else if(GJP.Bc(i) == BND_CND_APRD) fac = 2;
    else if(GJP.Bc(i) == BND_CND_GPARITY) fac = 4;
    else{ ERR.General("CPSfermion","getMomentumUnits","Unknown boundary condition"); }

    punits[i] = 6.283185308/(GJP.NodeSites(i)*GJP.Nodes(i)*fac); // 2pi/(fac*L)
  }
}

//Apply the phase exp(-ip.x) to each site of this vector, where p is a *three momentum*
//The units of the momentum are 2pi/L for periodic BCs, pi/L for antiperiodic BCs and pi/2L for G-parity BCs
//x_lcl is the site in node lattice coords. 3 or more dimensions (those after 3 are ignored)

struct _apply_phase_site_op_impl{
  template< typename mf_Complex, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy, typename my_enable_if<_equal<typename ComplexClassify<mf_Complex>::type,complex_double_or_float_mark>::value,int>::type = 0>
  inline static void apply_phase_site_op(CPSfermion<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy> &field, const int x_lcl[], const int &flav, const int p[], const double punits[]){
    assert(DimensionPolicy::EuclideanDimension >= 3);
    
    int x_glb[DimensionPolicy::EuclideanDimension]; for(int i=0;i<DimensionPolicy::EuclideanDimension;i++) x_glb[i] = x_lcl[i] + GJP.NodeCoor(i)*GJP.NodeSites(i);
    
    double phi = 0;
    for(int i=0;i<3;i++) phi += p[i]*punits[i]*x_glb[i];
    std::complex<double> phase( cos(phi), -sin(phi) );
    mf_Complex phase_prec(phase);

    mf_Complex *base = field.site_ptr(x_lcl,flav);
    for(int sc=0;sc<12;sc++){
      mf_Complex* v = base + sc;
      (*v) *= phase_prec;
    }
  }

#ifdef USE_GRID
  
  template< typename mf_Complex, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy, typename my_enable_if<_equal<typename ComplexClassify<mf_Complex>::type,grid_vector_complex_mark>::value,int>::type = 0>
  inline static void apply_phase_site_op(CPSfermion<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy> &field, const int x_lcl[], const int &flav, const int p[], const double punits[]){
    assert(DimensionPolicy::EuclideanDimension >= 3);

    int nsimd = field.Nsimd();

    typedef typename mf_Complex::scalar_type stype;
    stype* buf = (stype*)memalign(128, nsimd*sizeof(stype));

    int lane_off[DimensionPolicy::EuclideanDimension];
    int x_gbl_lane[DimensionPolicy::EuclideanDimension];
    
    for(int lane = 0; lane < nsimd; lane++){
      field.SIMDunmap(lane, lane_off);
      for(int xx=0;xx<DimensionPolicy::EuclideanDimension;xx++) x_gbl_lane[xx] = x_lcl[xx] + lane_off[xx] + GJP.NodeCoor(xx)*GJP.NodeSites(xx);
      
      double phi = 0;
      for(int i=0;i<3;i++) phi += p[i]*punits[i]*x_gbl_lane[i];

      buf[lane] = stype( cos(phi), -sin(phi) );
    }

    mf_Complex vphase;
    vset(vphase, buf);
    free(buf);

    mf_Complex* base = field.site_ptr(x_lcl,flav);
    for(int sc=0;sc<12;sc++){
      mf_Complex* v = base + sc;
      *v = vphase * (*v);
    }
  }

#endif
};




template< typename mf_Complex, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void CPSfermion<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>::apply_phase_site_op(const int x_lcl[], const int &flav, const int p[], const double punits[]){
  _apply_phase_site_op_impl::apply_phase_site_op(*this, x_lcl, flav, p, punits);
}  


//Apply gauge fixing matrices to the field
template< typename mf_Complex, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void CPSfermion4D<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>::gauge_fix_site_op(int fi, Lattice &lat,const bool dagger){
  int x4d[4]; int f; this->fsiteUnmap(fi,x4d,f);
  CPSfermion<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>::gauge_fix_site_op(x4d,f,lat,dagger);
}
template< typename mf_Complex, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void CPSfermion4D<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>::gaugeFix(Lattice &lat, const bool parallel, const bool dagger){
  if(parallel){
#pragma omp parallel for
    for(int fi=0;fi<this->nfsites();fi++)
      gauge_fix_site_op(fi,lat,dagger);
  }else{
    for(int fi=0;fi<this->nfsites();fi++)
      gauge_fix_site_op(fi,lat,dagger);
  }
}


//Apply the phase exp(-ip.x) to each site of this vector, where p is a *three momentum*
//The units of the momentum are 2pi/L for periodic BCs, pi/L for antiperiodic BCs and pi/2L for G-parity BCs
template< typename mf_Complex, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void CPSfermion4D<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>::apply_phase_site_op(int sf,const int p[],double punits[]){
  int x[this->EuclideanDimension]; int f; this->fsiteUnmap(sf,x,f);
  CPSfermion<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>::apply_phase_site_op(x,f,p,punits);
}

template< typename mf_Complex, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void CPSfermion4D<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>::applyPhase(const int p[], const bool &parallel){
  const char *fname = "apply_phase(int p[])";

  double punits[3];
  CPSfermion<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>::getMomentumUnits(punits);
  
  if(parallel){
#pragma omp parallel for
    for(int sf=0;sf<this->nfsites();sf++)
      apply_phase_site_op(sf,p,punits);
  }else{
    for(int sf=0;sf<this->nfsites();sf++)
      apply_phase_site_op(sf,p,punits);
  }
}


//Set the real and imaginary parts to uniform random numbers drawn from the appropriate local RNGs
template< typename mf_Complex, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void CPSfermion4D<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>::setUniformRandom(const Float &hi, const Float &lo){
  typedef typename mf_Complex::value_type mf_Float;
  LRG.SetInterval(hi,lo);
  for(int i = 0; i < this->sites*this->flavors; ++i) {
    int flav = i / this->sites;
    int st = i % this->sites;

    LRG.AssignGenerator(st,flav);
    mf_Float *p = (mf_Float*)this->site_ptr(st,flav);

    for(int site_lcl_off=0;site_lcl_off<2*FieldSiteSize;site_lcl_off++)
      *(p++) = LRG.Urand(FOUR_D);
  }
}

template< typename mf_Complex, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void CPSfermion4D<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>::setGaussianRandom(){
  typedef typename mf_Complex::value_type mf_Float;
  for(int i = 0; i < this->sites*this->flavors; ++i) {
    int flav = i / this->sites;
    int st = i % this->sites;

    LRG.AssignGenerator(st,flav);
    mf_Float *p = (mf_Float*)this->site_ptr(st,flav);

    for(int site_lcl_off=0;site_lcl_off<2*FieldSiteSize;site_lcl_off++)
      *(p++) = LRG.Grand(FOUR_D);
  }
}

template< typename mf_Complex, typename FlavorPolicy, typename AllocPolicy>
void CPSfermion5D<mf_Complex,FlavorPolicy,AllocPolicy>::setGaussianRandom(){
  typedef typename mf_Complex::value_type mf_Float;
  for(int i = 0; i < this->sites*this->flavors; ++i) {
    int flav = i / this->sites;
    int st = i % this->sites;

    LRG.AssignGenerator(st,flav);
    mf_Float *p = (mf_Float*)this->site_ptr(st,flav);

    for(int site_lcl_off=0;site_lcl_off<2*FieldSiteSize;site_lcl_off++)
      *(p++) = LRG.Grand(FIVE_D);
  }
}




//Gauge fix 3D fermion field with dynamic info type
template< typename mf_Complex, typename FlavorPolicy, typename AllocPolicy>
struct _ferm3d_gfix_impl{

  static void gaugeFix(CPSfermion3D<mf_Complex,FlavorPolicy,AllocPolicy> &field, Lattice &lat, const typename GaugeFix3DInfo<FlavorPolicy>::InfoType &t, const bool &parallel){
    if(GJP.Gparity() && field.nflavors() == 1) ERR.General("CPSfermion3D","gaugeFix(Lattice &, const int &, const bool &)","For one flavor fields with G-parity enabled, to gauge fix we need to know the flavor of this field\n");

#define LOOP								\
    for(int fi=0;fi<field.nfsites();fi++){				\
      int x4d[4]; int f; field.fsiteUnmap(fi,x4d,f);			\
      x4d[3] = t;							\
      field.CPSfermion<mf_Complex,SpatialPolicy,FlavorPolicy>::gauge_fix_site_op(x4d,f,lat); \
    }

    if(parallel){
#pragma omp parallel for
      LOOP;
    }else{
      LOOP;
    }
#undef LOOP
  }

};
//Partial specialization for one flavor. We must provide the flavor index for the gauge fixing matrix, i.e. the flavor that this field represents
template< typename mf_Complex, typename AllocPolicy>
struct _ferm3d_gfix_impl<mf_Complex,FixedFlavorPolicy<1>,AllocPolicy>{
  static void gaugeFix(CPSfermion3D<mf_Complex,FixedFlavorPolicy<1>,AllocPolicy> &field, Lattice &lat, const typename GaugeFix3DInfo<FixedFlavorPolicy<1> >::InfoType &time_flav, const bool &parallel){
    printf("_ferm3d_gfix_impl::gauge_fix with time=%d, flav=%d\n",time_flav.first,time_flav.second);
    typedef typename mf_Complex::value_type mf_Float;

#define SITE_OP								\
    int x4d[4]; field.siteUnmap(i,x4d);		\
    x4d[3] = time_flav.first;						\
    int gfmat_site = x4d[0] + GJP.XnodeSites()*( x4d[1] + GJP.YnodeSites()* ( x4d[2] + GJP.ZnodeSites()*x4d[3] )); \
    mf_Complex tmp[3];							\
    const Matrix* gfmat = lat.FixGaugeMatrix(gfmat_site,time_flav.second);	\
    mf_Complex* sc_base = field.site_ptr(x4d);			\
    for(int s=0;s<4;s++){						\
      memcpy(tmp, sc_base + 3 * s, 3 * sizeof(mf_Complex));		\
      colorMatrixMultiplyVector<mf_Float,Float>( (mf_Float*)(sc_base + 3*s), (Float*)gfmat, (mf_Float*)tmp); \
    }									

    if(parallel){
#pragma omp parallel for
      for(int i=0;i<field.nsites();i++){
	SITE_OP;
      }
    }else{
      for(int i=0;i<field.nsites();i++){
	SITE_OP;
      }
    }
#undef SITE_OP

  }


};


template< typename mf_Complex, typename FlavorPolicy, typename AllocPolicy>
void CPSfermion3D<mf_Complex,FlavorPolicy,AllocPolicy>::gaugeFix(Lattice &lat, const typename GaugeFix3DInfo<FlavorPolicy>::InfoType &t, const bool &parallel){
    _ferm3d_gfix_impl<mf_Complex,FlavorPolicy,AllocPolicy>::gaugeFix(*this,lat,t,parallel);
}


//Apply the phase exp(-ip.x) to each site of this vector, where p is a *three momentum*
//The units of the momentum are 2pi/L for periodic BCs, pi/L for antiperiodic BCs and pi/2L for G-parity BCs
template< typename mf_Complex, typename FlavorPolicy, typename AllocPolicy>
void CPSfermion3D<mf_Complex,FlavorPolicy,AllocPolicy>::apply_phase_site_op(const int &sf,const int p[],double punits[]){
  int x[this->Dimension]; int f; this->fsiteUnmap(sf,x,f);
  CPSfermion<mf_Complex,SpatialPolicy,FlavorPolicy,AllocPolicy>::apply_phase_site_op(x,f,p,punits);
}

template< typename mf_Complex, typename FlavorPolicy, typename AllocPolicy>
void CPSfermion3D<mf_Complex,FlavorPolicy,AllocPolicy>::applyPhase(const int p[], const bool &parallel){
  const char *fname = "apply_phase(int p[])";

  double punits[3];
  CPSfermion<mf_Complex,SpatialPolicy,FlavorPolicy>::getMomentumUnits(punits);
  
  if(parallel){
#pragma omp parallel for
    for(int sf=0;sf<this->nfsites();sf++)
      apply_phase_site_op(sf,p,punits);
  }else{
    for(int sf=0;sf<this->nfsites();sf++)
      apply_phase_site_op(sf,p,punits);
  }
}



















//Make a random complex scalar field of type
template< typename mf_Complex, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void CPScomplex4D<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>::setRandom(const RandomType &type){
  LRG.SetInterval(1, 0);
  for(int i = 0; i < this->sites*this->flavors; ++i) {
    int flav = i / this->sites;
    int st = i % this->sites;

    LRG.AssignGenerator(st,flav);
    mf_Complex *p = this->site_ptr(st,flav);
    RandomComplex<mf_Complex>::rand(p,type,FOUR_D);
  }
}

//Set the real and imaginary parts to uniform random numbers drawn from the appropriate local RNGs
template< typename mf_Complex, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void CPScomplex4D<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>::setUniformRandom(const Float &hi, const Float &lo){
  typedef typename mf_Complex::value_type mf_Float;
  LRG.SetInterval(hi,lo);
  for(int i = 0; i < this->sites*this->flavors; ++i) {
    int flav = i / this->sites;
    int st = i % this->sites;

    LRG.AssignGenerator(st,flav);
    mf_Float *p = (mf_Float*)this->site_ptr(st,flav);

    for(int i=0;i<2;i++)
      *(p++) = LRG.Urand(FOUR_D);
  }
}

 
//Perform the FFT
template< typename mf_Complex, typename FlavorPolicy, typename AllocPolicy>
void CPSglobalComplexSpatial<mf_Complex,FlavorPolicy,AllocPolicy>::fft(){
  typedef typename mf_Complex::value_type mf_Float;
  const int fft_dim[3] = {this->glb_size[2], this->glb_size[1], this->glb_size[0]};
  const int size_3d_glb = fft_dim[0] * fft_dim[1] * fft_dim[2];

  size_t this_floatsize = this->size() * 2;
  
  typename FFTWwrapper<mf_Float>::complexType* fft_mem = FFTWwrapper<mf_Float>::alloc_complex(this_floatsize);
  
  memcpy((void *)fft_mem, this->ptr(), this_floatsize*sizeof(mf_Float));

  //Plan creation is expensive, so make it static
  static typename FFTWwrapper<mf_Float>::planType plan_src;
  static bool init = false;
  if(!init){
    plan_src = FFTWwrapper<mf_Float>::plan_many_dft(3, fft_dim, 1,
						    fft_mem, NULL, 1, size_3d_glb,
						    fft_mem, NULL, 1, size_3d_glb,
						    FFTW_FORWARD, FFTW_ESTIMATE);
    init = true;
  }

  for(int f = 0; f < this->nflavors(); f++) {
    int off = f * size_3d_glb;
    FFTWwrapper<mf_Float>::execute_dft(plan_src, fft_mem + off, fft_mem + off);
  }

  memcpy((void *)this->ptr(), (void*)fft_mem, this_floatsize*sizeof(mf_Float));
  FFTWwrapper<mf_Float>::free(fft_mem);

  //FFTWwrapper<mf_Float>::cleanup(); //Don't need to cleanup, it doesn't have the function I initially thought
}
  
  

//Scatter to a local field

template< typename mf_Complex, typename FlavorPolicy, typename AllocPolicy,
	  typename extComplex, typename extDimPolicy, typename extAllocPolicy,
	  typename complex_class, int extEuclDim>
struct _CPSglobalComplexSpatial_scatter_impl{};

//Standard implementation for std::complex
template< typename mf_Complex, typename FlavorPolicy, typename AllocPolicy,
	  typename extComplex, typename extDimPolicy, typename extAllocPolicy>
struct _CPSglobalComplexSpatial_scatter_impl<mf_Complex,FlavorPolicy,AllocPolicy,  extComplex, extDimPolicy, extAllocPolicy, complex_double_or_float_mark, 3>{
  static void doit(CPSfield<extComplex,1,extDimPolicy,FlavorPolicy,extAllocPolicy> &to, const CPSglobalComplexSpatial<mf_Complex,FlavorPolicy,AllocPolicy> &from){
    const char *fname = "scatter(...)";
    int orig[3]; for(int i=0;i<3;i++) orig[i] = GJP.NodeSites(i)*GJP.NodeCoor(i);

#pragma omp parallel for
    for(int i=0;i<to.nfsites();i++){
      int x[3]; int flavor;  to.fsiteUnmap(i,x,flavor); //unmap the target coordinate
      for(int j=0;j<3;j++) x[j] += orig[j]; //global coord

      extComplex* tosite = to.fsite_ptr(i);
      mf_Complex const* fromsite = from.site_ptr(x,flavor);

      *tosite = *fromsite;
    }	
  }
};

#ifdef USE_GRID

//Implementation for Grid vector complex types
template< typename mf_Complex, typename FlavorPolicy, typename AllocPolicy,
	  typename extComplex, typename extDimPolicy, typename extAllocPolicy>
struct _CPSglobalComplexSpatial_scatter_impl<mf_Complex,FlavorPolicy,AllocPolicy,  extComplex, extDimPolicy, extAllocPolicy, grid_vector_complex_mark, 3>{
  static void doit(CPSfield<extComplex,1,extDimPolicy,FlavorPolicy,extAllocPolicy> &to, const CPSglobalComplexSpatial<mf_Complex,FlavorPolicy,AllocPolicy> &from){
    const char *fname = "scatter(...)";
    int orig[3]; for(int i=0;i<3;i++) orig[i] = GJP.NodeSites(i)*GJP.NodeCoor(i);

    const int ndim = 3;
    int nsimd = extComplex::Nsimd();
    std::vector<std::vector<int> > packed_offsets(nsimd,std::vector<int>(ndim)); //get the vector offsets for the different SIMD packed sites
    for(int i=0;i<nsimd;i++) to.SIMDunmap(i,&packed_offsets[i][0]);
    
#pragma omp parallel for
    for(int i=0;i<to.nfsites();i++){
      int x[3]; int flavor;  to.fsiteUnmap(i,x,flavor); //unmap the target coordinate. This is a root coordinate, we need to construct the other offsets
      for(int j=0;j<3;j++) x[j] += orig[j]; //global coord

      extComplex* toptr = to.fsite_ptr(i); 

      //x is the root coordinate corresponding to SIMD packed index 0      
      std::vector<mf_Complex const*> ptrs(nsimd);
      ptrs[0] = from.site_ptr(x,flavor);
      
      int xx[ndim];
      for(int i=1;i<nsimd;i++){
	for(int d=0;d<ndim;d++)
	  xx[d] = x[d] + packed_offsets[i][d];  //xx = x + offset
	ptrs[i] = from.site_ptr(xx,flavor);
      }
      to.SIMDpack(toptr, ptrs, 1);
    }	
  }
};

#endif


template< typename mf_Complex, typename FlavorPolicy, typename AllocPolicy>
template<typename extComplex, typename extDimPolicy, typename extAllocPolicy>
void CPSglobalComplexSpatial<mf_Complex,FlavorPolicy,AllocPolicy>::scatter(CPSfield<extComplex,1,extDimPolicy,FlavorPolicy,extAllocPolicy> &to) const{
  _CPSglobalComplexSpatial_scatter_impl<mf_Complex,FlavorPolicy,AllocPolicy,
					extComplex,extDimPolicy,extAllocPolicy,
					typename ComplexClassify<extComplex>::type,
					extDimPolicy::EuclideanDimension>::doit(to,*this);
}



template< typename SiteType, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy,
	  typename extSiteType, typename extDimPol, typename extAllocPol,
	  typename my_enable_if<intEq<DimensionPolicy::EuclideanDimension,extDimPol::EuclideanDimension>::val, int>::type = 0>
struct _gather_scatter_impl{
  typedef typename DimensionPolicy::EquivalentLocalPolicy EquivalentLocalPolicy;
  
  static void gather(CPSfieldGlobalInOneDir<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &into, const CPSfield<extSiteType,SiteSize,extDimPol,FlavorPolicy,extAllocPol> &from){
    NullObject n;
    CPSfield<SiteType,SiteSize,EquivalentLocalPolicy,FlavorPolicy,AllocPolicy> tmp(n);
    tmp.importField(from);
    _gather_scatter_impl<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy,
			 SiteType, EquivalentLocalPolicy, AllocPolicy>::gather(into, tmp);    
  }
  static void scatter(CPSfield<extSiteType,SiteSize,extDimPol,FlavorPolicy,extAllocPol> &to, const CPSfieldGlobalInOneDir<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &from){
    NullObject n;
    CPSfield<SiteType,SiteSize,EquivalentLocalPolicy,FlavorPolicy,AllocPolicy> tmp(n);
    _gather_scatter_impl<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy,
			 SiteType, EquivalentLocalPolicy, AllocPolicy>::scatter(tmp, from);
    to.importField(tmp);
  }
    
};
template< typename SiteType, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy,
	  typename my_enable_if<intEq<DimensionPolicy::EuclideanDimension,DimensionPolicy::EquivalentLocalPolicy::EuclideanDimension>::val, int>::type test>
struct _gather_scatter_impl<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy,
		    SiteType, typename DimensionPolicy::EquivalentLocalPolicy, AllocPolicy, test>{
  typedef typename DimensionPolicy::EquivalentLocalPolicy LocalDimensionPolicy;

  static void gather(CPSfieldGlobalInOneDir<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &into, const CPSfield<SiteType,SiteSize,LocalDimensionPolicy,FlavorPolicy,AllocPolicy> &from){
    assert(LocalDimensionPolicy::EuclideanDimension == DimensionPolicy::EuclideanDimension);
    const int &dir = into.getDir();

    const char *fname = "gather(...)";
    NullObject nullobj;
    CPSfield<SiteType,SiteSize,LocalDimensionPolicy,FlavorPolicy,AllocPolicy> tmp1(nullobj);
    CPSfield<SiteType,SiteSize,LocalDimensionPolicy,FlavorPolicy,AllocPolicy> tmp2(nullobj);
    CPSfield<SiteType,SiteSize,LocalDimensionPolicy,FlavorPolicy,AllocPolicy>* send = const_cast<CPSfield<SiteType,SiteSize,LocalDimensionPolicy,FlavorPolicy,AllocPolicy>* >(&from);
    CPSfield<SiteType,SiteSize,LocalDimensionPolicy,FlavorPolicy,AllocPolicy>* recv = &tmp2;

    int cur_dir_origin = GJP.NodeSites(dir)*GJP.NodeCoor(dir);    
    int size_in_Float = from.size() * sizeof(SiteType) / sizeof(IFloat); //getPlusData measures the send/recv size in units of sizeof(IFloat)

    int nshift = GJP.Nodes(dir);

    for(int shift = 0; shift < nshift; shift++){
#pragma omp parallel for
      for(int i=0;i<send->nfsites();i++){
	int x[DimensionPolicy::EuclideanDimension]; int flavor;  send->fsiteUnmap(i,x,flavor); //unmap the buffer coordinate
	x[dir] += cur_dir_origin; //now a global coordinate in the dir direction

	SiteType* tosite = into.site_ptr(x,flavor);
	SiteType* fromsite = send->fsite_ptr(i);

	memcpy((void*)tosite, (void*)fromsite, into.siteSize()*sizeof(SiteType));
      }	

      if(shift != nshift-1){
	getPlusData((IFloat*)recv->ptr(), (IFloat*)send->ptr(), size_in_Float, dir);
	cur_dir_origin += GJP.NodeSites(dir);
	cur_dir_origin %= (GJP.NodeSites(dir)*GJP.Nodes(dir));

	if(shift == 0){
	  recv = &tmp1;
	  send = &tmp2;
	}else std::swap(send,recv);
      }
    }    
  }

  static void scatter(CPSfield<SiteType,SiteSize,LocalDimensionPolicy,FlavorPolicy,AllocPolicy> &to, const CPSfieldGlobalInOneDir<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &from){
    assert(LocalDimensionPolicy::EuclideanDimension == DimensionPolicy::EuclideanDimension);
    
    const int &dir = from.getDir();
    
    const char *fname = "scatter(...)";
    int cur_dir_origin = GJP.NodeSites(dir)*GJP.NodeCoor(dir);

#pragma omp parallel for
    for(int i=0;i<to.nfsites();i++){
      int x[DimensionPolicy::EuclideanDimension]; int flavor;  to.fsiteUnmap(i,x, flavor); //unmap the target coordinate
      x[dir] += cur_dir_origin; //now a global coordinate in the dir direction
      
      SiteType* tosite = to.fsite_ptr(i);
      SiteType const* fromsite = from.site_ptr(x,flavor);
      
      memcpy((void*)tosite, (void*)fromsite, from.siteSize()*sizeof(SiteType));
    }
  }
  
};




//Gather up the row. Involves internode communication
template< typename SiteType, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
template<typename extSiteType, typename extDimPol, typename extAllocPol>
void CPSfieldGlobalInOneDir<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy>::gather(const CPSfield<extSiteType,SiteSize,extDimPol,FlavorPolicy,extAllocPol> &from){
  _gather_scatter_impl<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy,
	       extSiteType, extDimPol, extAllocPol>::gather(*this, from);  
}


//Scatter back out. Involves no communication
template< typename SiteType, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
template<typename extSiteType, typename extDimPol, typename extAllocPol>
void CPSfieldGlobalInOneDir<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy>::scatter(CPSfield<extSiteType,SiteSize,extDimPol,FlavorPolicy,extAllocPol> &to) const{
  _gather_scatter_impl<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy,
		       extSiteType, extDimPol, extAllocPol>::scatter(to, *this);
}

#define FFT_MULTI

#ifndef FFT_MULTI

//Perform a fast Fourier transform along the principal direction
//NOTE: This won't work correctly if the DimensionPolicy does not use canonical ordering: FIXME
//Assumes SiteType is a std::complex type
template< typename SiteType, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void CPSfieldGlobalInOneDir<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy>::fft(const bool inverse_transform){
  const int dir = this->getDir();
  const char* fname = "fft()";
  
  //We do a large number of simple linear FFTs. This field has its principal direction as the fastest changing index so this is nice and easy
  int sc_size = this->siteSize(); //we have to assume the sites comprise complex numbers
  int size_1d_glb = GJP.NodeSites(dir) * GJP.Nodes(dir);
  const int n_fft = this->nsites() / GJP.NodeSites(dir) * sc_size * this->nflavors();

  //Plan creation is expensive, so make it static and only re-create if the field size changes
  //Create a plan for each direction because we can have non-cubic spatial volumes
  static FFTplanContainer<typename SiteType::value_type> plan_f[4];
  static bool plan_init = false;
  static int plan_sc_size;
  static bool plan_inv_trans;
  if(!plan_init || sc_size != plan_sc_size || inverse_transform != plan_inv_trans){ //recreate/create
    typename FFTWwrapper<typename SiteType::value_type>::complexType *tmp_f; //I don't think it actually does anything with this

    for(int i=0;i<4;i++){
      int size_i = GJP.NodeSites(i) * GJP.Nodes(i);

      plan_f[i].setPlan(1, &size_i, 1, 
			tmp_f, NULL, sc_size, size_i * sc_size,
			tmp_f, NULL, sc_size, size_i * sc_size,
			inverse_transform ? FFTW_BACKWARD : FFTW_FORWARD, FFTW_ESTIMATE);  
    }
    plan_sc_size = sc_size;
    plan_inv_trans = inverse_transform;
    plan_init = true;
  }
    
  typename FFTWwrapper<typename SiteType::value_type>::complexType *fftw_mem = FFTWwrapper<typename SiteType::value_type>::alloc_complex(size_1d_glb * n_fft);
    
  memcpy((void *)fftw_mem, this->ptr(), this->size()*sizeof(SiteType));
#pragma omp parallel for
  for(int n = 0; n < n_fft; n++) {
    int sc_id = n % sc_size;
    int chunk_id = n / sc_size; //3d block index
    int off = size_1d_glb * sc_size * chunk_id + sc_id;
    FFTWwrapper<typename SiteType::value_type>::execute_dft(plan_f[dir].getPlan(), fftw_mem + off, fftw_mem + off); 
  }

  //FFTWwrapper<SiteType>::cleanup(); //I think this actually destroys existing plans!

  if(!inverse_transform) memcpy(this->ptr(), (void *)fftw_mem, this->size()*sizeof(SiteType));
  else for(int i=0;i<this->size();i++) this->ptr()[i] = *( (SiteType*)fftw_mem+i )/double(size_1d_glb);
  
  FFTWwrapper<typename SiteType::value_type>::free(fftw_mem);
}

#else

template< typename SiteType, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void CPSfieldGlobalInOneDir<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy>::fft(const bool inverse_transform){
  const int dir = this->getDir();
  const char* fname = "fft()";
  
  //We do a large number of simple linear FFTs. This field has its principal direction as the fastest changing index so this is nice and easy
  int sc_size = this->siteSize(); //we have to assume the sites comprise complex numbers
  int size_1d_glb = GJP.NodeSites(dir) * GJP.Nodes(dir);
  const int n_fft = this->nsites() / GJP.NodeSites(dir) * this->nflavors();

  //Plan creation is expensive, so make it static and only re-create if the field size changes
  //Create a plan for each direction because we can have non-cubic spatial volumes
  static FFTplanContainer<typename SiteType::value_type> plan_f[4];
  static bool plan_init = false;
  static int plan_sc_size;
  static bool plan_inv_trans;
  if(!plan_init || sc_size != plan_sc_size || inverse_transform != plan_inv_trans){ //recreate/create
    typename FFTWwrapper<typename SiteType::value_type>::complexType *tmp_f; //I don't think it actually does anything with this

    for(int i=0;i<4;i++){
      int size_i = GJP.NodeSites(i) * GJP.Nodes(i);

      plan_f[i].setPlan(1, &size_i, sc_size, 
			tmp_f, NULL, sc_size, 1,
			tmp_f, NULL, sc_size, 1,
			inverse_transform ? FFTW_BACKWARD : FFTW_FORWARD, FFTW_ESTIMATE);  
    }
    plan_sc_size = sc_size;
    plan_inv_trans = inverse_transform;
    plan_init = true;
  }

  typename FFTWwrapper<typename SiteType::value_type>::complexType *fftw_mem = FFTWwrapper<typename SiteType::value_type>::alloc_complex(size_1d_glb * n_fft * sc_size);
    
  memcpy((void *)fftw_mem, this->ptr(), this->size()*sizeof(SiteType));
#pragma omp parallel for
  for(int n = 0; n < n_fft; n++) {
    int chunk_id = n; //3d block index
    int off = size_1d_glb * sc_size * chunk_id;
    FFTWwrapper<typename SiteType::value_type>::execute_dft(plan_f[dir].getPlan(), fftw_mem + off, fftw_mem + off); 
  }

  if(!inverse_transform) memcpy(this->ptr(), (void *)fftw_mem, this->size()*sizeof(SiteType));
  else for(int i=0;i<this->size();i++) this->ptr()[i] = *( (SiteType*)fftw_mem+i )/double(size_1d_glb);
  
  FFTWwrapper<typename SiteType::value_type>::free(fftw_mem);
}

#endif





#endif
