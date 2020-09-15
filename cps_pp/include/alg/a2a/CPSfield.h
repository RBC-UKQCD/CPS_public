#ifndef CPS_FIELD_H
#define CPS_FIELD_H

#include<algorithm>
#include<comms/scu.h>
#include<alg/a2a/fftw_wrapper.h>
#include<alg/a2a/utils.h>
#include<alg/a2a/CPSfield_policies.h>
CPS_START_NAMESPACE 

typedef std::complex<float> ComplexF;
typedef std::complex<double> ComplexD;

//A wrapper for a CPS-style field. Most functionality is generic so it can do quite a lot of cool things
template< typename SiteType, int SiteSize, typename DimensionPolicy, typename FlavorPolicy = DynamicFlavorPolicy, typename AllocPolicy = StandardAllocPolicy>
class CPSfield: public DimensionPolicy, public FlavorPolicy, public AllocPolicy{
  SiteType* f;
protected:
  int sites; //number of Euclidean sites
  int flavors; //number of flavors
  int fsites; //number of generalized sites (including flavor)

  int fsize; //number of SiteType in the array = SiteSize * fsites
  
  void alloc(){
    this->_alloc((void**)&f, fsize*sizeof(SiteType));
  }
  void freemem(){
    this->_free((void*)f);
  }

public:
  enum { FieldSiteSize = SiteSize };
  typedef SiteType FieldSiteType;
  typedef DimensionPolicy FieldDimensionPolicy;
  typedef FlavorPolicy FieldFlavorPolicy;
  typedef AllocPolicy FieldAllocPolicy;
  
  typedef typename DimensionPolicy::ParamType InputParamType;

  CPSfield(const InputParamType &params): DimensionPolicy(params){
    this->setFlavors(flavors); //from FlavorPolicy
    this->setSites(sites,fsites,flavors); //from DimensionPolicy
    fsize = fsites * SiteSize;
    alloc();
  }
  CPSfield(const CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &r): fsize(r.fsize), flavors(r.flavors),sites(r.sites),fsites(r.fsites), DimensionPolicy(r){
    alloc();
    memcpy(f,r.f,sizeof(SiteType) * fsize);
  }

  //Copy from external pointer. Make sure you set the params and policies correctly because it has no way of bounds checking
  CPSfield(SiteType const* copyme, const InputParamType &params): DimensionPolicy(params){
    this->setFlavors(flavors); //from FlavorPolicy
    this->setSites(sites,fsites,flavors); //from DimensionPolicy
    fsize = fsites * SiteSize;
    alloc();
    memcpy(f,copyme,sizeof(SiteType) * fsize);
  }
  
  //Set the field to zero
  void zero(){
    memset(f, 0, sizeof(SiteType) * fsize);      
  }

  CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &operator=(const CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &r){
    static_cast<DimensionPolicy&>(*this) = r; //copy policy info
    
    sites = r.sites;
    fsites = r.fsites;
    flavors = r.flavors;

    int old_fsize = fsize;
    fsize = r.fsize;

    if(fsize != old_fsize){
     freemem();
      alloc();
    }
    memcpy(f,r.f,sizeof(SiteType) * fsize);
    return *this;
  }

  static std::size_t byte_size(const InputParamType &params){
    CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,NullAllocPolicy> tmp(params); //doesn't allocate
    std::size_t out = SiteSize * sizeof(SiteType);
    return tmp.nfsites() * out;
  }
  std::size_t byte_size() const{
    return this->nfsites() * SiteSize * sizeof(SiteType);
  }
  
  //Set each element to a uniform random number in the specified range.
  //WARNING: Uses only the current RNG in LRG, and does not change this based on site. This is therefore only useful for testing*
  void testRandom(const Float hi = 0.5, const Float lo = -0.5);

  inline const int nsites() const{ return sites; }
  inline const int nflavors() const{ return flavors; }
  inline const int nfsites() const{ return fsites; } //number of generalized sites including flavor

  //Number of SiteType per site
  inline const int siteSize() const{ return SiteSize; }

  //Number of SiteType in field
  inline const int size() const{ return fsize; }

  //Accessors
  inline SiteType* ptr(){ return f; }
  inline SiteType const* ptr() const{ return f; }

  //Accessors *do not check bounds*
  //int fsite is the linearized N-dimensional site/flavorcoordinate with the mapping specified by the policy class
  inline int fsite_offset(const int fsite) const{ return SiteSize*fsite; }
  
  inline SiteType* fsite_ptr(const int fsite){  //fsite is in the internal flavor/Euclidean mapping of the DimensionPolicy. Use only if you know what you are doing
    return f + SiteSize*fsite;
  }
  inline SiteType const* fsite_ptr(const int fsite) const{  //fsite is in the internal flavor/Euclidean mapping of the DimensionPolicy. Use only if you know what you are doing
    return f + SiteSize*fsite;
  }

  //int site is the linearized N-dimension Euclidean coordinate with mapping specified by the policy class
  inline int site_offset(const int site, const int flav = 0) const{ return SiteSize*this->siteFsiteConvert(site,flav); }
  inline int site_offset(const int x[], const int flav = 0) const{ return SiteSize*this->fsiteMap(x,flav); }

  inline SiteType* site_ptr(const int site, const int flav = 0){  //site is in the internal Euclidean mapping of the DimensionPolicy
    return f + SiteSize*this->siteFsiteConvert(site,flav);
  }
  inline SiteType* site_ptr(const int x[], const int flav = 0){ 
    return f + SiteSize*this->fsiteMap(x,flav);
  }    

  inline SiteType const* site_ptr(const int site, const int flav = 0) const{  //site is in the internal Euclidean mapping of the DimensionPolicy
    return f + SiteSize*this->siteFsiteConvert(site,flav);
  }
  inline SiteType const* site_ptr(const int x[], const int flav = 0) const{ 
    return f + SiteSize*this->fsiteMap(x,flav);
  }    
 
  inline int flav_offset() const{ return SiteSize*this->fsiteFlavorOffset(); } //pointer offset between flavors

  //Set this field to the average of this and a second field, r
  void average(const CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &r, const bool &parallel = true);

  //Self destruct initialized (no more sfree!!)
  virtual ~CPSfield(){
    freemem();
  }

  //Import an export field with arbitrary DimensionPolicy (must have same Euclidean dimension!) and precision. Must have same SiteSize and FlavorPolicy
  template< typename extSiteType, typename extDimPol, typename extFlavPol, typename extAllocPol>
  void importField(const CPSfield<extSiteType,SiteSize,extDimPol,extFlavPol,extAllocPol> &r);

  template< typename extSiteType, typename extDimPol, typename extFlavPol, typename extAllocPol>
  void exportField(CPSfield<extSiteType,SiteSize,extDimPol,extFlavPol,extAllocPol> &r) const;

  bool equals(const CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &r) const{
    for(int i=0;i<fsize;i++)
      if(f[i] != r.f[i]) return false;
    return true;
  }

#define CONDITION is_double_or_float<typename extField::FieldSiteType>::value \
  && is_double_or_float<SiteType>::value \
  && _equal<DimensionPolicy, typename extField::FieldDimensionPolicy>::value \
  && _equal<FlavorPolicy, typename extField::FieldFlavorPolicy>::value
  
  template<typename extField>
  bool equals(const extField &r, typename my_enable_if<CONDITION,const double>::type tolerance) const{
    for(int i=0;i<fsize;i++){
      if( fabs(f[i] - r.f[i]) > tolerance) return false;
    }
    return true;
  }
#undef CONDITION
  
#define CONDITION is_complex_double_or_float<typename extField::FieldSiteType>::value \
  && is_complex_double_or_float<SiteType>::value \
  && _equal<DimensionPolicy, typename extField::FieldDimensionPolicy>::value \
  && _equal<FlavorPolicy, typename extField::FieldFlavorPolicy>::value
  
  template<typename extField>
  bool equals(const extField &r, typename my_enable_if<CONDITION,const double>::type tolerance, bool verbose = false) const{
    for(int i=0;i<fsize;i++){
      if( fabs(f[i].real() - r.f[i].real()) > tolerance || fabs(f[i].imag() - r.f[i].imag()) > tolerance ){
	if(verbose && !UniqueID()){
	  int rem = i;
	  int s = rem % SiteSize; rem /= SiteSize;
	  int x = rem % sites; rem /= sites;
	  int flav = rem;
	  int coor[DimensionPolicy::EuclideanDimension]; this->siteUnmap(x,coor);
	  std::ostringstream os; for(int a=0;a<DimensionPolicy::EuclideanDimension;a++) os << coor[a] << " ";
	  std::string coor_str = os.str();
	  
	  printf("Err: off %d  [s=%d coor=(%s) f=%d] this[%g,%g] vs that[%g,%g] : diff [%g,%g]\n",i, s,coor_str.c_str(),flav,
		 f[i].real(),f[i].imag(),r.f[i].real(),r.f[i].imag(),fabs(f[i].real()-r.f[i].real()), fabs(f[i].imag()-r.f[i].imag()) );
	}
	return false;
      }
    }
    return true;
  }
#undef CONDITION

#ifdef USE_GRID
  
#define CONDITION _equal<  typename ComplexClassify<typename extField::FieldSiteType>::type  ,  grid_vector_complex_mark>::value \
  && _equal<typename ComplexClassify<SiteType>::type,grid_vector_complex_mark>::value \
  && _equal<DimensionPolicy, typename extField::FieldDimensionPolicy>::value \
  && _equal<FlavorPolicy, typename extField::FieldFlavorPolicy>::value

  template<typename extField>
  bool equals(const extField &r, typename my_enable_if<CONDITION,const double>::type tolerance, bool verbose = false) const{
    typedef typename SiteType::scalar_type ThisScalarType;
    typedef typename extField::FieldSiteType::scalar_type ThatScalarType;
    typedef typename DimensionPolicy::EquivalentScalarPolicy ScalarDimPol;
    NullObject null_obj;
    CPSfield<ThisScalarType,SiteSize,ScalarDimPol, FlavorPolicy, StandardAllocPolicy> tmp_this(null_obj);
    CPSfield<ThatScalarType,SiteSize,ScalarDimPol, FlavorPolicy, StandardAllocPolicy> tmp_that(null_obj);
    tmp_this.importField(*this);
    tmp_that.importField(r);
    return tmp_this.equals(tmp_that,tolerance,verbose);
  }
  
#undef CONDITION
#endif
  
  double norm2() const;
  
#ifdef USE_GRID
  //Import for Grid Lattice<blah> types
  template<typename GridField>
  void importGridField(const GridField &grid);
  
  template<typename GridField>
  void exportGridField(GridField &grid) const;
#endif

  CPSfield & operator+=(const CPSfield &r){
#pragma omp parallel for
    for(int i=0;i<fsize;i++) f[i] += r.f[i];
    return *this;
  }
  CPSfield & operator-=(const CPSfield &r){
#pragma omp parallel for
    for(int i=0;i<fsize;i++) f[i] -= r.f[i];
    return *this;
  }

  CPSfield operator+(const CPSfield &r) const{
    CPSfield out(*this); out += r;
    return out;
  }
  CPSfield operator-(const CPSfield &r){
    CPSfield out(*this); out -= r;
    return out;
  }
};

#define INHERIT_TYPEDEFS(...) \
  typedef typename __VA_ARGS__::FieldSiteType FieldSiteType; \
  typedef typename __VA_ARGS__::FieldDimensionPolicy FieldDimensionPolicy; \
  typedef typename __VA_ARGS__::FieldFlavorPolicy FieldFlavorPolicy; \
  typedef typename __VA_ARGS__::FieldAllocPolicy FieldAllocPolicy; \
  typedef typename __VA_ARGS__::InputParamType InputParamType; \
  enum { FieldSiteSize = __VA_ARGS__::FieldSiteSize }


#define DEFINE_ADDSUB_DERIVED(DerivedType) \
  DerivedType & operator+=(const DerivedType &r){ \
    this->CPSfield<FieldSiteType,FieldSiteSize,FieldDimensionPolicy,FieldFlavorPolicy,FieldAllocPolicy>::operator+=(r); return *this; \
  } \
  DerivedType & operator-=(const DerivedType &r){ \
    this->CPSfield<FieldSiteType,FieldSiteSize,FieldDimensionPolicy,FieldFlavorPolicy,FieldAllocPolicy>::operator-=(r); return *this; \
  } \
  DerivedType operator+(const DerivedType &r) const{ \
    DerivedType out(*this); out += r; \
    return out; \
  } \
  DerivedType operator-(const DerivedType &r){ \
    DerivedType out(*this); out -= r; \
    return out; \
  }




template< typename mf_Complex, typename DimensionPolicy, typename FlavorPolicy = DynamicFlavorPolicy, typename AllocPolicy = StandardAllocPolicy>
class CPSfermion: public CPSfield<mf_Complex,12,DimensionPolicy,FlavorPolicy,AllocPolicy>{
protected:
  void gauge_fix_site_op(const int x4d[], const int &f, Lattice &lat,const bool dagger = false);

  static void getMomentumUnits(double punits[3]);

  //Apply the phase exp(-ip.x) to each site of this vector, where p is a *three momentum*
  //The units of the momentum are 2pi/L for periodic BCs, pi/L for antiperiodic BCs and pi/2L for G-parity BCs
  //x_lcl is the site in node lattice coords
  void apply_phase_site_op(const int x_lcl[], const int &flav, const int p[], const double punits[]);

public:
  INHERIT_TYPEDEFS(CPSfield<mf_Complex,12,DimensionPolicy,FlavorPolicy,AllocPolicy>);
  
  CPSfermion(): CPSfield<mf_Complex,12,DimensionPolicy,FlavorPolicy,AllocPolicy>(NullObject()){} //default constructor won't compile if policies need arguments
  CPSfermion(const InputParamType &params): CPSfield<mf_Complex,12,DimensionPolicy,FlavorPolicy,AllocPolicy>(params){}
  CPSfermion(const CPSfermion<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy> &r): CPSfield<mf_Complex,12,DimensionPolicy,FlavorPolicy,AllocPolicy>(r){}
};

template<typename FlavorPolicy>
struct GaugeFix3DInfo{};

template<>
struct GaugeFix3DInfo<DynamicFlavorPolicy>{
  typedef int InfoType;
};
template<>
struct GaugeFix3DInfo<FixedFlavorPolicy<2> >{
  typedef int InfoType;
};
template<>
struct GaugeFix3DInfo<FixedFlavorPolicy<1> >{
  typedef std::pair<int,int> InfoType; //time, flavor (latter ignored if no GPBC)
};

template< typename mf_Complex, typename FlavorPolicy = DynamicFlavorPolicy, typename AllocPolicy = StandardAllocPolicy>
class CPSfermion3D: public CPSfermion<mf_Complex,SpatialPolicy,FlavorPolicy,AllocPolicy>{
  void apply_phase_site_op(const int &sf,const int p[],double punits[]);

  template< typename mf_Complex2, typename FlavorPolicy2>
  friend struct _ferm3d_gfix_impl;
public:
  INHERIT_TYPEDEFS(CPSfermion<mf_Complex,SpatialPolicy,FlavorPolicy,AllocPolicy>);
  
  CPSfermion3D(): CPSfermion<mf_Complex,SpatialPolicy,FlavorPolicy>(){}
  CPSfermion3D(const CPSfermion3D<mf_Complex> &r): CPSfermion<mf_Complex,SpatialPolicy,FlavorPolicy,AllocPolicy>(r){}

  //Apply gauge fixing matrices to the field
  //Because this is a 3d field we must also provide a time coordinate.
  //If the field is one flavor we must also provide the flavor
  //We make the field_info type dynamic based on the FlavorPolicy for this reason (pretty cool!)
  void gaugeFix(Lattice &lat, const typename GaugeFix3DInfo<FlavorPolicy>::InfoType &field_info, const bool &parallel);

  //Apply the phase exp(-ip.x) to each site of this vector, where p is a *three momentum*
  //The units of the momentum are 2pi/L for periodic BCs, pi/L for antiperiodic BCs and pi/2L for G-parity BCs
  void applyPhase(const int p[], const bool &parallel);

  DEFINE_ADDSUB_DERIVED(CPSfermion3D);
};

template< typename mf_Complex, typename DimensionPolicy = FourDpolicy, typename FlavorPolicy = DynamicFlavorPolicy, typename AllocPolicy = StandardAllocPolicy>
class CPSfermion4D: public CPSfermion<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>{
  typename my_enable_if<DimensionPolicy::EuclideanDimension == 4, int>::type dummy;
  void gauge_fix_site_op(int fi, Lattice &lat, const bool dagger = false);
  void apply_phase_site_op(int sf,const int p[],double punits[]);
public:  
  INHERIT_TYPEDEFS(CPSfermion<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>);
  
  CPSfermion4D(): CPSfermion<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>(){}
  CPSfermion4D(const InputParamType &params): CPSfermion<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>(params){}
  CPSfermion4D(const CPSfermion4D<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy> &r): CPSfermion<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy>(r){}

  //Apply gauge fixing matrices to the field. 
  //NOTE: This does not work correctly for GPBC and FlavorPolicy==FixedFlavorPolicy<1> because we need to provide the flavor 
  //that this field represents to obtain the gauge-fixing matrix. I fixed this for CPSfermion3D and a similar implementation will work here
  //dagger = true  applied V^\dagger to the vector to invert a previous gauge fix
  void gaugeFix(Lattice &lat, const bool parallel, const bool dagger = false);

  //Apply the phase exp(-ip.x) to each site of this vector, where p is a *three momentum*
  //The units of the momentum are 2pi/L for periodic BCs, pi/L for antiperiodic BCs and pi/2L for G-parity BCs
  void applyPhase(const int p[], const bool &parallel);

  //Set the real and imaginary parts to uniform random numbers drawn from the appropriate local RNGs
  void setUniformRandom(const Float &hi = 0.5, const Float &lo = -0.5);

  void setGaussianRandom();

  DEFINE_ADDSUB_DERIVED(CPSfermion4D);
};

template< typename mf_Complex, typename FlavorPolicy = DynamicFlavorPolicy, typename AllocPolicy = StandardAllocPolicy>
class CPSfermion5D: public CPSfield<mf_Complex,12,FiveDpolicy,FlavorPolicy,AllocPolicy>{
public:
  INHERIT_TYPEDEFS(CPSfield<mf_Complex,12,FiveDpolicy,FlavorPolicy,AllocPolicy>);
  
  CPSfermion5D(): CPSfield<mf_Complex,12,FiveDpolicy,FlavorPolicy,AllocPolicy>(NullObject()){}
  CPSfermion5D(const CPSfermion5D<mf_Complex,FlavorPolicy,AllocPolicy> &r): CPSfield<mf_Complex,12,FiveDpolicy,FlavorPolicy,AllocPolicy>(r){}
  
#ifdef USE_BFM
private:
  template<typename FloatExt>
  void impexFermion(Fermion_t bfm_field, const int cb, const int do_import, bfm_qdp<FloatExt> &dwf){
    if(this->flavors == 2) assert(dwf.gparity);

    const int sc_incr = dwf.nsimd() * 2; //stride between spin-color indices
    FloatExt * bb = (FloatExt*)bfm_field;

    typedef typename mf_Complex::value_type mf_Float;
    
#pragma omp parallel for
    for(int fs=0;fs<this->fsites;fs++){
      int x[5], f; this->fsiteUnmap(fs);
      if( (x[0]+x[1]+x[2]+x[3] + (dwf.precon_5d ? x[4] : 0)) % 2 == cb){
	mf_Float* cps_base = (mf_Float*)this->fsite_ptr(fs);

	int bidx_off = dwf.gparity ? 
	  dwf.bagel_gparity_idx5d(x, x[4], 0, 0, 12, 1, f) :
	  dwf.bagel_idx5d(x, x[4], 0, 0, 12, 1);

	FloatExt * bfm_base = bb + bidx_off;

	for(int i=0;i<12;i++)
	  for(int reim=0;reim<2;reim++)
	    if(do_import)
	      *(cps_base + 2*i + reim) = *(bfm_base + 2*sc_incr*i + reim);
	    else
	      *(bfm_base + 2*sc_incr*i + reim) = *(cps_base + 2*i + reim);
      }
    }
  }
public:
  enum { FieldSiteSize = 12 };
	 
  template<typename FloatExt>
  void importFermion(const Fermion_t bfm_field, const int cb, bfm_qdp<FloatExt> &dwf){
    impexFermion<FloatExt>(const_cast<Fermion_t>(bfm_field), cb, 1, dwf);
  }
  template<typename FloatExt>
  void exportFermion(const Fermion_t bfm_field, const int cb, bfm_qdp<FloatExt> &dwf) const{
    const_cast<CPSfermion5D<mf_Complex,FlavorPolicy,AllocPolicy>*>(this)->impexFermion<FloatExt>(bfm_field, cb, 0, dwf);
  }
#endif

  void setGaussianRandom();
  
  DEFINE_ADDSUB_DERIVED(CPSfermion5D);
};


template< typename mf_Complex, typename DimensionPolicy = FourDpolicy, typename FlavorPolicy = DynamicFlavorPolicy, typename AllocPolicy = StandardAllocPolicy>
class CPScomplex4D: public CPSfield<mf_Complex,1,DimensionPolicy,FlavorPolicy,AllocPolicy>{
  typename my_enable_if<DimensionPolicy::EuclideanDimension == 4, int>::type dummy;
public:
  INHERIT_TYPEDEFS(CPSfield<mf_Complex,1,DimensionPolicy,FlavorPolicy,AllocPolicy>);
  
  CPScomplex4D(): CPSfield<mf_Complex,1,DimensionPolicy,FlavorPolicy,AllocPolicy>(NullObject()){}
  CPScomplex4D(const InputParamType &params): CPSfield<mf_Complex,1,DimensionPolicy,FlavorPolicy,AllocPolicy>(params){}
  CPScomplex4D(const CPScomplex4D<mf_Complex,DimensionPolicy,FlavorPolicy,AllocPolicy> &r): CPSfield<mf_Complex,1,DimensionPolicy,FlavorPolicy,AllocPolicy>(r){}

  //Make a random complex scalar field of type
  void setRandom(const RandomType &type);

  //Set the real and imaginary parts to uniform random numbers drawn from the appropriate local RNGs
  void setUniformRandom(const Float &hi = 0.5, const Float &lo = -0.5);

  DEFINE_ADDSUB_DERIVED(CPScomplex4D);
};

//3d complex number field
template< typename mf_Complex, typename FlavorPolicy = DynamicFlavorPolicy, typename AllocPolicy = StandardAllocPolicy>
class CPScomplexSpatial: public CPSfield<mf_Complex,1,SpatialPolicy,FlavorPolicy,AllocPolicy>{
public:
  INHERIT_TYPEDEFS(CPSfield<mf_Complex,1,SpatialPolicy,FlavorPolicy,AllocPolicy>);
  
  CPScomplexSpatial(): CPSfield<mf_Complex,1,SpatialPolicy,FlavorPolicy,AllocPolicy>(NullObject()){}
  CPScomplexSpatial(const CPScomplexSpatial<mf_Complex,FlavorPolicy,AllocPolicy> &r): CPSfield<mf_Complex,1,SpatialPolicy,FlavorPolicy,AllocPolicy>(r){}

  DEFINE_ADDSUB_DERIVED(CPScomplexSpatial);
};

//Lattice-spanning 'global' 3d complex field
template< typename mf_Complex, typename FlavorPolicy = DynamicFlavorPolicy, typename AllocPolicy = StandardAllocPolicy>
class CPSglobalComplexSpatial: public CPSfield<mf_Complex,1,GlobalSpatialPolicy,FlavorPolicy,AllocPolicy>{
  
  template< typename _mf_Complex, typename _FlavorPolicy, typename _AllocPolicy,
	    typename extComplex, typename extDimPolicy, typename extAllocPolicy,
	    typename complex_class, int extEuclDim>
  friend struct _CPSglobalComplexSpatial_scatter_impl;
public:
  INHERIT_TYPEDEFS(CPSfield<mf_Complex,1,GlobalSpatialPolicy,FlavorPolicy,AllocPolicy>);
  
  CPSglobalComplexSpatial(): CPSfield<mf_Complex,1,GlobalSpatialPolicy,FlavorPolicy,AllocPolicy>(NullObject()){}
  CPSglobalComplexSpatial(const CPSglobalComplexSpatial<mf_Complex,FlavorPolicy,AllocPolicy> &r): CPSfield<mf_Complex,1,GlobalSpatialPolicy,FlavorPolicy,AllocPolicy>(r){}
  
  //Perform the FFT
  void fft();

  //Scatter to a local field
  template<typename extComplex, typename extDimPolicy, typename extAllocPolicy>
  void scatter(CPSfield<extComplex,1,extDimPolicy,FlavorPolicy,extAllocPolicy> &to) const;

  DEFINE_ADDSUB_DERIVED(CPSglobalComplexSpatial);
};


//This field contains an entire row of sub-lattices along a particular dimension. Every node along that row contains an identical copy
template< typename SiteType, int SiteSize, typename DimensionPolicy, typename FlavorPolicy = DynamicFlavorPolicy, typename AllocPolicy = StandardAllocPolicy>
class CPSfieldGlobalInOneDir: public CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy>{
public:
  INHERIT_TYPEDEFS(CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy>);
  
  CPSfieldGlobalInOneDir(const int &dir): CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy>(dir){}
  CPSfieldGlobalInOneDir(const CPSfieldGlobalInOneDir<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &r): CPSfield<SiteType,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy>(r){}

  //Gather up the row. Involves internode communication
  template<typename extSiteType, typename extDimPol, typename extAllocPol>
  void gather(const CPSfield<extSiteType,SiteSize,extDimPol,FlavorPolicy,extAllocPol> &from);

  //Scatter back out. Involves no communication
  template<typename extSiteType, typename extDimPol, typename extAllocPol>
  void scatter(CPSfield<extSiteType,SiteSize,extDimPol,FlavorPolicy,extAllocPol> &to) const;

  //Perform a fast Fourier transform along the principal direction. It currently assumes the DimensionPolicy has the sites mapped in canonical ordering
  void fft(const bool inverse_transform = false);

  DEFINE_ADDSUB_DERIVED(CPSfieldGlobalInOneDir);
};

template< typename mf_Complex, typename FlavorPolicy = DynamicFlavorPolicy, typename AllocPolicy = StandardAllocPolicy>
class CPSfermion4DglobalInOneDir: public CPSfieldGlobalInOneDir<mf_Complex,12,FourDglobalInOneDir,FlavorPolicy,AllocPolicy>{
public:
  INHERIT_TYPEDEFS(CPSfieldGlobalInOneDir<mf_Complex,12,FourDglobalInOneDir,FlavorPolicy,AllocPolicy>);
  
  CPSfermion4DglobalInOneDir(const int &dir): CPSfieldGlobalInOneDir<mf_Complex,12,FourDglobalInOneDir,FlavorPolicy,AllocPolicy>(dir){}
  CPSfermion4DglobalInOneDir(const CPSfermion4DglobalInOneDir<mf_Complex,FlavorPolicy,AllocPolicy> &r): CPSfieldGlobalInOneDir<mf_Complex,12,FourDglobalInOneDir,FlavorPolicy,AllocPolicy>(r){}

  DEFINE_ADDSUB_DERIVED(CPSfermion4DglobalInOneDir);
};
template< typename mf_Complex, typename FlavorPolicy = DynamicFlavorPolicy, typename AllocPolicy = StandardAllocPolicy>
class CPSfermion3DglobalInOneDir: public CPSfieldGlobalInOneDir<mf_Complex,12,ThreeDglobalInOneDir,FlavorPolicy,AllocPolicy>{
public:
  INHERIT_TYPEDEFS(CPSfieldGlobalInOneDir<mf_Complex,12,ThreeDglobalInOneDir,FlavorPolicy,AllocPolicy>);
  
  CPSfermion3DglobalInOneDir(const int &dir): CPSfieldGlobalInOneDir<mf_Complex,12,ThreeDglobalInOneDir,FlavorPolicy,AllocPolicy>(dir){}
  CPSfermion3DglobalInOneDir(const CPSfermion3DglobalInOneDir<mf_Complex,FlavorPolicy,AllocPolicy> &r): CPSfieldGlobalInOneDir<mf_Complex,12,ThreeDglobalInOneDir,FlavorPolicy,AllocPolicy>(r){}

  DEFINE_ADDSUB_DERIVED(CPSfermion3DglobalInOneDir);
};


////////Checkerboarded types/////////////
template< typename mf_Complex, typename CBpolicy, typename FlavorPolicy = DynamicFlavorPolicy, typename AllocPolicy = StandardAllocPolicy>
class CPSfermion5Dprec: public CPSfield<mf_Complex,12,FiveDevenOddpolicy<CBpolicy>,FlavorPolicy,AllocPolicy>{
public:
  INHERIT_TYPEDEFS(CPSfield<mf_Complex,12,FiveDevenOddpolicy<CBpolicy>,FlavorPolicy,AllocPolicy>);
  
  CPSfermion5Dprec(): CPSfield<mf_Complex,12,FiveDevenOddpolicy<CBpolicy>,FlavorPolicy,AllocPolicy>(NullObject()){}
  CPSfermion5Dprec(const CPSfermion5Dprec<mf_Complex,CBpolicy,FlavorPolicy,AllocPolicy> &r): CPSfield<mf_Complex,12,FiveDevenOddpolicy<CBpolicy>,FlavorPolicy,AllocPolicy>(r){}

  DEFINE_ADDSUB_DERIVED(CPSfermion5Dprec);
};


template< typename mf_Complex, typename FlavorPolicy = DynamicFlavorPolicy, typename AllocPolicy = StandardAllocPolicy>
class CPSfermion5Dcb4Deven: public CPSfermion5Dprec<mf_Complex,CheckerBoard<4,0>,FlavorPolicy,AllocPolicy>{
public:
  INHERIT_TYPEDEFS(CPSfermion5Dprec<mf_Complex,CheckerBoard<4,0>,FlavorPolicy,AllocPolicy>);

  CPSfermion5Dcb4Deven(): CPSfermion5Dprec<mf_Complex,CheckerBoard<4,0>,FlavorPolicy,AllocPolicy>(){}
  CPSfermion5Dcb4Deven(const CPSfermion5Dcb4Deven<mf_Complex,FlavorPolicy,AllocPolicy> &r): CPSfermion5Dprec<mf_Complex,CheckerBoard<4,0>,FlavorPolicy,AllocPolicy>(r){}

  DEFINE_ADDSUB_DERIVED(CPSfermion5Dcb4Deven);
};
template< typename mf_Complex, typename FlavorPolicy = DynamicFlavorPolicy, typename AllocPolicy = StandardAllocPolicy>
class CPSfermion5Dcb4Dodd: public CPSfermion5Dprec<mf_Complex,CheckerBoard<4,1>,FlavorPolicy,AllocPolicy>{
public:
  INHERIT_TYPEDEFS(CPSfermion5Dprec<mf_Complex,CheckerBoard<4,1>,FlavorPolicy,AllocPolicy>);
  
  CPSfermion5Dcb4Dodd(): CPSfermion5Dprec<mf_Complex,CheckerBoard<4,1>,FlavorPolicy,AllocPolicy>(){}
  CPSfermion5Dcb4Dodd(const CPSfermion5Dcb4Dodd<mf_Complex,FlavorPolicy,AllocPolicy> &r): CPSfermion5Dprec<mf_Complex,CheckerBoard<4,1>,FlavorPolicy,AllocPolicy>(r){}

  DEFINE_ADDSUB_DERIVED(CPSfermion5Dcb4Dodd);
};





#include<alg/a2a/CPSfield_impl.h>

CPS_END_NAMESPACE
#endif
